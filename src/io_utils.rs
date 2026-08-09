use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Opens `path` for reading, transparently decompressing gzip/bgzip input.
/// Detection is by magic bytes rather than file extension, so gzip data
/// piped through a renamed file still works. `MultiGzDecoder` is required
/// (not `GzDecoder`) because bgzip-compressed genome references are valid
/// concatenated multi-member gzip streams.
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let mut file = File::open(path)
        .with_context(|| format!("failed to open input file: {}", path.display()))?;

    let mut magic = [0u8; 2];
    let read_n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;

    if read_n == 2 && magic == GZIP_MAGIC {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Opens `path` for writing. If the path ends in `.gz`, output is gzip
/// compressed on the fly; otherwise it is written as plain text.
pub fn open_writer(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path)
        .with_context(|| format!("failed to create output file: {}", path.display()))?;
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    if is_gz {
        Ok(Box::new(GzEncoder::new(
            BufWriter::new(file),
            Compression::default(),
        )))
    } else {
        Ok(Box::new(BufWriter::new(file)))
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

pub fn write_lines(writer: &mut dyn Write, lines: &[String]) -> io::Result<()> {
    for line in lines {
        writer.write_all(line.as_bytes())?;
        writer.write_all(b"\n")?;
    }
    Ok(())
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
/// members to the file in the same order the blocks were given -- so
/// output is deterministic and byte-order-preserving despite the
/// compression happening out of order across threads.
pub enum BlockWriter {
    Gzip {
        file: BufWriter<File>,
        level: Compression,
    },
    Plain(BufWriter<File>),
}

/// Opens `path` for batched writing, gzip-compressing in parallel members
/// if the path ends in `.gz`. See `BlockWriter` for why this differs from
/// `open_writer`.
pub fn open_block_writer(path: &Path) -> Result<BlockWriter> {
    let file = File::create(path)
        .with_context(|| format!("failed to create output file: {}", path.display()))?;
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    if is_gz {
        Ok(BlockWriter::Gzip {
            file: BufWriter::new(file),
            level: Compression::default(),
        })
    } else {
        Ok(BlockWriter::Plain(BufWriter::new(file)))
    }
}

impl BlockWriter {
    /// Writes `blocks` in order. Empty blocks are skipped (no point paying
    /// a gzip member's ~20-byte fixed overhead for zero content, e.g. a
    /// low-sampling-rate chunk that happened to select nothing in one
    /// sub-split).
    pub fn write_blocks(&mut self, blocks: Vec<Vec<u8>>) -> Result<()> {
        match self {
            BlockWriter::Gzip { file, level } => {
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
                    file.write_all(&member)?;
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

    pub fn flush(&mut self) -> io::Result<()> {
        match self {
            BlockWriter::Gzip { file, .. } => file.flush(),
            BlockWriter::Plain(w) => w.flush(),
        }
    }
}
