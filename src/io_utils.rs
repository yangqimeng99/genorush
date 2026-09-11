use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, IsTerminal, Read, Write};
use std::num::NonZero;
use std::path::Path;

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use noodles_bgzf as bgzf;
use rayon::prelude::*;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Buffer size for the reader that sits directly on the file or on stdin.
/// Large enough that a `Stdin` handle's internal per-`read` lock is amortised
/// over a big chunk rather than paid per small read.
const IO_BUF: usize = 256 * 1024;

/// The conventional stand-in for stdin/stdout on a command line. A path
/// equal to this is never opened as a file.
pub const STDIO: &str = "-";

/// Whether `path` refers to stdin/stdout rather than a file on disk.
pub fn is_stdio(path: &Path) -> bool {
    path.as_os_str() == STDIO
}

/// Rejects a set of paths that would all have to be the same stream.
///
/// Two inputs reading `-` would race for the same stdin, each getting an
/// arbitrary half of it; two outputs writing `-` would interleave two
/// different record streams into one pipe. Both produce plausible-looking
/// garbage rather than an error, so they are refused up front. `role`
/// names what is being checked, for the error message.
pub fn ensure_one_stdio_at_most(paths: &[&Path], role: &str) -> Result<()> {
    let n = paths.iter().filter(|p| is_stdio(p)).count();
    if n > 1 {
        bail!(
            "{n} {role}s were given as `-`, but there is only one stdin/stdout to go around; \
             at most one {role} can use `-`, the rest must be files"
        );
    }
    Ok(())
}

/// How an output should be compressed, for outputs whose path can't say.
#[derive(Debug, Clone, Copy, Default)]
pub struct OutputOpts {
    /// Force gzip compression regardless of the output path's extension.
    /// Set by the global `--gzip` flag, which is the only way to ask for
    /// compressed output on stdout (`-o -` has no extension to inspect).
    pub gzip: bool,

    /// Write BGZF instead of plain gzip, and compress whether or not the
    /// path asked for it. Set by the global `--bgzf` flag.
    ///
    /// BGZF is gzip -- any reader that takes `.gz` takes this -- but it is
    /// gzip constrained so that it can be indexed: every member holds at
    /// most 64 KiB of uncompressed data and records its own compressed size,
    /// and the stream ends with a distinguished empty block that proves it
    /// was not truncated. That is what `tabix`, `csi`, and every random-access
    /// reader in the ecosystem require. A `.vcf.gz` written as ordinary gzip
    /// looks correct, decompresses correctly, and cannot be indexed -- which
    /// is exactly the kind of quiet dead end worth a flag of its own.
    pub bgzf: bool,
}

/// Opens `path` for reading, transparently decompressing gzip/bgzip input,
/// or reads stdin when `path` is `-`.
///
/// Detection is by magic bytes rather than file extension, so gzip data
/// piped through a renamed file (or through a pipe, which has no name at
/// all) still works. The two magic bytes are *peeked* via `BufRead::fill_buf`
/// rather than read-and-rewound: a pipe cannot seek back, and peeking works
/// the same on both kinds of input, so there is one code path instead of
/// two. `MultiGzDecoder` is required (not `GzDecoder`) because
/// bgzip-compressed references, and this tool's own `BlockWriter` output,
/// are valid concatenated multi-member gzip streams.
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let raw: Box<dyn Read + Send> = if is_stdio(path) {
        // `Stdin`, not `StdinLock`: the lock guard is not `Send`, and these
        // readers get moved onto their own thread by `spawn_reader`.
        Box::new(io::stdin())
    } else {
        Box::new(
            File::open(path)
                .with_context(|| format!("failed to open input file: {}", path.display()))?,
        )
    };

    let mut buf = BufReader::with_capacity(IO_BUF, raw);
    let is_gzip = {
        let head = buf
            .fill_buf()
            .with_context(|| format!("failed to read from: {}", display(path)))?;
        head.len() >= 2 && head[..2] == GZIP_MAGIC
    };

    if is_gzip {
        Ok(Box::new(BufReader::with_capacity(
            IO_BUF,
            MultiGzDecoder::new(buf),
        )))
    } else {
        Ok(Box::new(buf))
    }
}

/// How a path should be described in messages: `-` is not a filename.
pub fn display(path: &Path) -> String {
    if is_stdio(path) {
        "<stdin/stdout>".to_string()
    } else {
        path.display().to_string()
    }
}

/// Reads at most `max_lines` lines from `reader` into `out`, trimming
/// leading/trailing whitespace from each line (mirrors Python's `str.strip()`
/// semantics used by the original script). Returns the number of lines read;
/// 0 means EOF.
pub fn read_line_chunk(
    reader: &mut dyn BufRead,
    out: &mut Vec<String>,
    max_lines: usize,
) -> Result<usize> {
    out.clear();
    let mut buf = String::new();
    let mut n_read = 0;
    while n_read < max_lines {
        buf.clear();
        let n = reader.read_line(&mut buf).context("failed reading line")?;
        if n == 0 {
            break;
        }
        out.push(buf.trim().to_string());
        n_read += 1;
    }
    Ok(n_read)
}

/// A writer for commands that already buffer output in batches and want
/// gzip compression to scale with `-j` instead of being a single-threaded
/// bottleneck.
///
/// A standard gzip stream can't be *decompressed* in parallel (DEFLATE's
/// back-references make it inherently sequential), but nothing stops
/// *compressing* independent chunks of input in parallel and concatenating
/// the results: RFC 1952 defines a gzip file as a sequence of one or more
/// independently-decodable "members", and any conforming reader (including
/// this project's own `open_reader`/`MultiGzDecoder`, plus `gzip`, `zcat`,
/// and every bioinformatics tool that accepts `.gz` input) reads a
/// concatenation of members exactly as if it were one. This is the same
/// technique `pigz` and `bgzip` use for their own compression speedups.
///
/// `write_blocks` takes each caller-provided block, compresses it into its
/// own gzip member (in parallel across blocks, via rayon), and writes the
/// members to the sink in the same order the blocks were given -- so
/// output is deterministic and byte-order-preserving despite the
/// compression happening out of order across threads.
///
/// The sink is a trait object so the same batching and parallel compression
/// applies whether the destination is a file or stdout. Blocks are large
/// (hundreds of KB up), so the one dynamic call per block costs nothing
/// measurable next to the compression it wraps.
pub enum BlockWriter {
    Gzip {
        sink: BufWriter<Box<dyn Write + Send>>,
        level: Compression,
    },
    /// Held in an `Option` so that finishing -- which is what writes the EOF
    /// block -- happens exactly once however many times `flush` is called.
    Bgzf(Option<bgzf::io::MultithreadedWriter<BufWriter<Box<dyn Write + Send>>>>),
    Plain(BufWriter<Box<dyn Write + Send>>),
}

/// Whether this write should be refused: gzip bytes aimed at a terminal.
///
/// Outputs default to stdout, so forgetting to pipe or redirect is easy.
/// Plain text landing in the terminal is the usual Unix outcome and stays
/// allowed -- that is what looking at a few records means. A gzip stream is
/// different: it is binary, it can garble the terminal, and no one asks for
/// it on purpose. Split out from the writer so the rule can be tested
/// without a terminal attached.
fn is_binary_to_terminal(stdout_bound: bool, gzip: bool, stdout_is_tty: bool) -> bool {
    stdout_bound && gzip && stdout_is_tty
}

/// Opens `path` for batched writing, or stdout when `path` is `-`.
///
/// Output is gzip-compressed when the path ends in `.gz` or when
/// `opts.gzip` is set. Stdout has no extension to inspect, so `--gzip` is
/// the only way to ask for compressed output there. See `BlockWriter` for
/// why this is not a plain `Write`.
pub fn open_block_writer(path: &Path, opts: OutputOpts) -> Result<BlockWriter> {
    let named_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    let is_gz = opts.gzip || opts.bgzf || named_gz;

    if is_binary_to_terminal(is_stdio(path), is_gz, io::stdout().is_terminal()) {
        bail!(
            "refusing to write gzip-compressed output to the terminal. Redirect it to a file \
             or pipe it into another command, or drop -z/--gzip to get plain text here"
        );
    }

    let sink: Box<dyn Write + Send> = if is_stdio(path) {
        // `Stdout` wraps a `LineWriter`, but `write_blocks` hands it whole
        // blocks at a time, so that costs one newline scan per block rather
        // than a flush per line.
        Box::new(io::stdout())
    } else {
        Box::new(
            File::create(path)
                .with_context(|| format!("failed to create output file: {}", path.display()))?,
        )
    };
    let sink = BufWriter::with_capacity(IO_BUF, sink);

    if opts.bgzf {
        // BGZF does its own framing -- 64 KiB members, each recording its
        // compressed size -- so the blocks this writer is handed are just
        // bytes to it. Its workers come from the same `-j` budget as
        // everything else.
        let workers = NonZero::new(rayon::current_num_threads()).unwrap_or(NonZero::<usize>::MIN);
        Ok(BlockWriter::Bgzf(Some(
            bgzf::io::MultithreadedWriter::with_worker_count(workers, sink),
        )))
    } else if is_gz {
        Ok(BlockWriter::Gzip {
            sink,
            level: Compression::default(),
        })
    } else {
        Ok(BlockWriter::Plain(sink))
    }
}

impl BlockWriter {
    /// Writes `blocks` in order. Empty blocks are skipped (no point paying
    /// a gzip member's ~20-byte fixed overhead for zero content, e.g. a
    /// low-sampling-rate chunk that happened to select nothing in one
    /// sub-split).
    pub fn write_blocks(&mut self, blocks: Vec<Vec<u8>>) -> Result<()> {
        match self {
            BlockWriter::Gzip { sink, level } => {
                let level = *level;
                let compressed: Vec<Vec<u8>> = blocks
                    .into_par_iter()
                    .filter(|b| !b.is_empty())
                    .map(|block| -> Result<Vec<u8>> {
                        let mut enc = GzEncoder::new(Vec::new(), level);
                        enc.write_all(&block)?;
                        Ok(enc.finish()?)
                    })
                    .collect::<Result<Vec<_>>>()?;
                for member in compressed {
                    sink.write_all(&member)?;
                }
                Ok(())
            }
            BlockWriter::Bgzf(writer) => {
                let writer = writer
                    .as_mut()
                    .context("BGZF output was already finished")?;
                for block in blocks {
                    if !block.is_empty() {
                        writer.write_all(&block)?;
                    }
                }
                Ok(())
            }
            BlockWriter::Plain(w) => {
                for block in blocks {
                    if !block.is_empty() {
                        w.write_all(&block)?;
                    }
                }
                Ok(())
            }
        }
    }

    /// Finishes the output. For BGZF this is what appends the EOF block that
    /// marks the stream complete, so it must happen before the process ends
    /// -- a BGZF file without it reads fine until something checks, and then
    /// reads as truncated.
    pub fn flush(&mut self) -> io::Result<()> {
        match self {
            BlockWriter::Gzip { sink, .. } => sink.flush(),
            BlockWriter::Bgzf(writer) => match writer.take() {
                Some(mut w) => {
                    let mut sink = w.finish()?;
                    sink.flush()
                }
                None => Ok(()),
            },
            BlockWriter::Plain(w) => w.flush(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::testutil::scratch_path as scratch;
    use std::path::PathBuf;

    #[test]
    fn recognizes_the_stdio_path() {
        assert!(is_stdio(Path::new("-")));
        assert!(!is_stdio(Path::new("./-")));
        assert!(!is_stdio(Path::new("-.fastq")));
        assert!(!is_stdio(Path::new("reads.fq")));
    }

    #[test]
    fn allows_at_most_one_stdio_per_role() {
        let dash = PathBuf::from("-");
        let file = PathBuf::from("reads.fq");
        assert!(ensure_one_stdio_at_most(&[&dash, &file], "input").is_ok());
        assert!(ensure_one_stdio_at_most(&[&file, &file], "input").is_ok());
        let err = ensure_one_stdio_at_most(&[&dash, &dash], "input").unwrap_err();
        assert!(
            err.to_string().contains("only one stdin/stdout"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn describes_stdio_without_pretending_it_is_a_file() {
        assert_eq!(display(Path::new("-")), "<stdin/stdout>");
        assert_eq!(display(Path::new("reads.fq")), "reads.fq");
    }

    fn round_trip(path: &Path, opts: OutputOpts, payload: &[u8]) -> (bool, Vec<u8>) {
        let mut w = open_block_writer(path, opts).unwrap();
        w.write_blocks(vec![payload.to_vec()]).unwrap();
        w.flush().unwrap();

        let raw = std::fs::read(path).unwrap();
        let compressed = raw.len() >= 2 && raw[..2] == GZIP_MAGIC;

        // Reading back goes through the same magic-byte peek every command
        // uses, so this also covers detection without a seekable rewind.
        let mut reader = open_reader(path).unwrap();
        let mut back = Vec::new();
        reader.read_to_end(&mut back).unwrap();
        (compressed, back)
    }

    #[test]
    fn only_refuses_compressed_output_to_a_terminal() {
        // gzip + stdout + tty is the one combination worth refusing.
        assert!(is_binary_to_terminal(true, true, true));
        // Plain text to a terminal is how you look at a few records.
        assert!(!is_binary_to_terminal(true, false, true));
        // Piped or redirected: the reader is another program, not a screen.
        assert!(!is_binary_to_terminal(true, true, false));
        // A file output is never the terminal, whatever stdout is doing.
        assert!(!is_binary_to_terminal(false, true, true));
    }

    #[test]
    fn plain_output_stays_plain_and_reads_back() {
        let p = scratch("plain");
        let (compressed, back) = round_trip(&p, OutputOpts::default(), b"@r\nACGT\n+\nIIII\n");
        assert!(
            !compressed,
            "no .gz extension and no --gzip: must stay plain"
        );
        assert_eq!(back, b"@r\nACGT\n+\nIIII\n");
    }

    #[test]
    fn gz_extension_compresses_and_reads_back() {
        let p = scratch("ext").with_extension("gz");
        let (compressed, back) = round_trip(&p, OutputOpts::default(), b"ACGTACGTACGT\n");
        assert!(compressed, ".gz extension must compress");
        assert_eq!(back, b"ACGTACGTACGT\n");
    }

    #[test]
    fn gzip_flag_compresses_regardless_of_extension() {
        let p = scratch("forced");
        let (compressed, back) = round_trip(
            &p,
            OutputOpts {
                gzip: true,
                bgzf: false,
            },
            b"ACGTACGTACGT\n",
        );
        assert!(compressed, "--gzip must compress even without a .gz name");
        assert_eq!(back, b"ACGTACGTACGT\n");
    }

    /// The 28 bytes every BGZF stream ends with: an empty member that says
    /// "this file is complete". A reader that does not find it knows the
    /// file was truncated -- which is the whole reason it exists.
    const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    fn looks_gzip(bytes: &[u8]) -> bool {
        bytes.len() >= 2 && bytes[..2] == GZIP_MAGIC
    }

    fn write_with(path: &Path, opts: OutputOpts, payload: &[u8]) -> Vec<u8> {
        let mut w = open_block_writer(path, opts).unwrap();
        w.write_blocks(vec![payload.to_vec()]).unwrap();
        w.flush().unwrap();
        std::fs::read(path).unwrap()
    }

    #[test]
    fn bgzf_output_is_gzip_that_declares_its_block_sizes() {
        let p = scratch("bgzf_shape");
        let raw = write_with(
            &p,
            OutputOpts {
                gzip: false,
                bgzf: true,
            },
            b"ACGTACGTACGT\n",
        );

        assert!(looks_gzip(&raw), "BGZF is gzip");
        // FEXTRA must be set, and the extra field must carry the `BC`
        // subfield holding the member's own compressed size. Without it a
        // reader cannot jump from member to member, which is the entire
        // difference between this and ordinary gzip.
        assert_eq!(raw[3] & 0x04, 0x04, "FEXTRA flag should be set");
        assert_eq!(
            &raw[12..14],
            b"BC",
            "the BC subfield should follow the header"
        );
    }

    #[test]
    fn bgzf_output_ends_with_the_marker_that_proves_it_is_complete() {
        let p = scratch("bgzf_eof");
        let raw = write_with(
            &p,
            OutputOpts {
                gzip: false,
                bgzf: true,
            },
            b"ACGT\n",
        );
        assert!(raw.len() > BGZF_EOF.len());
        assert_eq!(
            &raw[raw.len() - BGZF_EOF.len()..],
            &BGZF_EOF,
            "a BGZF stream without its EOF block reads as truncated"
        );
    }

    #[test]
    fn plain_gzip_output_has_no_bgzf_framing() {
        // The contrast is the point: this file decompresses identically and
        // cannot be indexed.
        let p = scratch("plain_gz_shape").with_extension("gz");
        let raw = write_with(&p, OutputOpts::default(), b"ACGT\n");
        assert!(looks_gzip(&raw));
        assert_ne!(raw[3] & 0x04, 0x04, "plain gzip sets no extra field here");
        assert!(!raw.ends_with(&BGZF_EOF));
    }

    #[test]
    fn bgzf_output_reads_back_through_the_ordinary_reader() {
        let p = scratch("bgzf_roundtrip");
        let payload: Vec<u8> = (0..200_000u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        write_with(
            &p,
            OutputOpts {
                gzip: false,
                bgzf: true,
            },
            &payload,
        );

        // Spanning many 64 KiB members: nothing downstream should need to
        // know this was BGZF rather than gzip.
        let mut back = Vec::new();
        open_reader(&p).unwrap().read_to_end(&mut back).unwrap();
        assert_eq!(back, payload);
    }

    #[test]
    fn bgzf_compresses_even_when_the_path_says_nothing() {
        let p = scratch("bgzf_unnamed");
        let raw = write_with(
            &p,
            OutputOpts {
                gzip: false,
                bgzf: true,
            },
            b"ACGT\n",
        );
        assert!(looks_gzip(&raw), "--bgzf implies compression");
    }

    #[test]
    fn finishing_bgzf_twice_is_harmless() {
        // Commands call flush at the end of their own paths; finishing twice
        // must not append a second EOF block or fail.
        let p = scratch("bgzf_double_finish");
        let mut w = open_block_writer(
            &p,
            OutputOpts {
                gzip: false,
                bgzf: true,
            },
        )
        .unwrap();
        w.write_blocks(vec![b"ACGT\n".to_vec()]).unwrap();
        w.flush().unwrap();
        w.flush().unwrap();
        let raw = std::fs::read(&p).unwrap();
        assert_eq!(&raw[raw.len() - BGZF_EOF.len()..], &BGZF_EOF);
        assert!(
            !raw[..raw.len() - BGZF_EOF.len()].ends_with(&BGZF_EOF),
            "the EOF block should appear once"
        );
    }

    #[test]
    fn multi_member_output_reads_back_as_one_stream() {
        // Every block becomes its own gzip member; a reader must see them
        // as a single continuous stream.
        let p = scratch("members").with_extension("gz");
        let mut w = open_block_writer(&p, OutputOpts::default()).unwrap();
        w.write_blocks(vec![b"first\n".to_vec(), Vec::new(), b"second\n".to_vec()])
            .unwrap();
        w.write_blocks(vec![b"third\n".to_vec()]).unwrap();
        w.flush().unwrap();

        let mut back = String::new();
        open_reader(&p).unwrap().read_to_string(&mut back).unwrap();
        assert_eq!(back, "first\nsecond\nthird\n");
    }
}
