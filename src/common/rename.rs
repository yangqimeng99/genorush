//! Shared engine behind every `<format> rename` leaf command: load an
//! old-name -> new-name mapping table, then stream the input through a
//! format-specific per-line transform in parallel chunks. Format-specific
//! code (fastx::rename, gff::rename) only supplies the transform closure.

use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{ensure, Context, Result};
use clap::Args;
use rayon::prelude::*;

use crate::io_utils::{
    display, ensure_one_stdio_at_most, open_block_writer, open_reader, read_line_chunk,
    BlockWriter, OutputOpts, STDIO,
};

#[derive(Args, Debug)]
pub struct RenameCommonArgs {
    /// Input file, or `-` to read stdin. Gzip/bgzip is auto-detected from
    /// the data itself, so piped-in compressed input works too.
    #[arg(value_name = "FILE")]
    pub input: PathBuf,

    /// Mapping file: two whitespace-separated columns per line, `new_name old_name`.
    /// May itself be gzip compressed, and may be `-` (but not when the
    /// input is also `-`).
    #[arg(short = 'n', long)]
    pub name: PathBuf,

    /// Output file. Defaults to stdout, so the result can be piped straight
    /// into the next command; `-` means the same thing written out.
    /// Gzip-compressed if the path ends in `.gz` or `-z/--gzip` is passed.
    #[arg(short = 'o', long, default_value = STDIO)]
    pub output: PathBuf,

    /// Number of lines processed per parallel batch. Bounds peak memory
    /// independently of input size; raise it for fewer, larger batches on
    /// machines with headroom, lower it if memory is tight.
    #[arg(long, default_value_t = 200_000)]
    pub chunk_lines: usize,
}

pub fn load_name_dict(path: &std::path::Path) -> Result<HashMap<String, String>> {
    use std::io::BufRead;

    let reader = open_reader(path)?;
    let mut map = HashMap::new();
    for (i, raw) in reader.lines().enumerate() {
        let raw = raw.with_context(|| format!("failed reading name file line {}", i + 1))?;
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split_whitespace();
        let new_name = fields.next().with_context(|| {
            format!(
                "name file line {}: expected 2 columns, got: {:?}",
                i + 1,
                line
            )
        })?;
        let old_name = fields.next().with_context(|| {
            format!(
                "name file line {}: expected 2 columns, got: {:?}",
                i + 1,
                line
            )
        })?;
        map.insert(old_name.to_string(), new_name.to_string());
    }
    Ok(map)
}

/// Formats one batch of lines straight into the byte blocks that
/// `BlockWriter::write_blocks` compresses in parallel.
///
/// The obvious shape for this -- transform each line into a `String`,
/// collect them, then join -- allocates once per line and copies every line
/// twice. Neither is needed: each rayon task owns one output buffer and
/// appends into it, so an unchanged line is a `memcpy` and a rewritten one
/// is a few appends, with no per-line allocation either way. The buffers are
/// sized up front from the lines they will hold, so they never grow.
fn format_into_blocks(
    lines: &[String],
    dict: &HashMap<String, String>,
    transform: &(impl Fn(&str, &HashMap<String, String>, &mut Vec<u8>) + Sync),
) -> Vec<Vec<u8>> {
    if lines.is_empty() {
        return Vec::new();
    }
    let n = rayon::current_num_threads().max(1).min(lines.len());
    let group_size = lines.len().div_ceil(n);
    lines
        .par_chunks(group_size)
        .map(|group| {
            // A rewritten line can differ in length from its input, but not
            // by much -- a name swap -- so this is the right order of
            // magnitude and usually exact.
            let mut buf = Vec::with_capacity(group.iter().map(|l| l.len() + 1).sum());
            for line in group {
                transform(line, dict, &mut buf);
                buf.push(b'\n');
            }
            buf
        })
        .collect()
}

/// Streams `args.input` to `args.output`, applying `transform` to every line
/// in parallel batches of `args.chunk_lines`. `transform` receives the
/// trimmed line and the loaded name dictionary.
///
/// `transform` appends the transformed line to a buffer rather than
/// returning a new one. Returning `String` allocated once per line, and
/// returning `Cow` only helped the lines that pass through unchanged -- on a
/// GFF whose seqids are all in the mapping, every line is rewritten and the
/// `Cow` was pure overhead. Writing into a caller-owned buffer costs nothing
/// in either direction.
pub fn run(
    args: &RenameCommonArgs,
    transform: impl Fn(&str, &HashMap<String, String>, &mut Vec<u8>) + Sync,
    opts: OutputOpts,
) -> Result<()> {
    ensure!(args.chunk_lines > 0, "--chunk-lines must be > 0");
    ensure_one_stdio_at_most(&[&args.input, &args.name], "input")?;

    let start = Instant::now();
    log::info!("loading name mapping from {}", display(&args.name));
    let dict = load_name_dict(&args.name)?;
    log::info!("loaded {} name mappings", dict.len());

    log::info!(
        "processing {} -> {}",
        display(&args.input),
        display(&args.output)
    );

    let mut reader = open_reader(&args.input)?;
    let mut writer: BlockWriter = open_block_writer(&args.output, opts)?;

    let mut chunk = Vec::with_capacity(args.chunk_lines);
    let mut total_lines: u64 = 0;
    loop {
        let n = read_line_chunk(reader.as_mut(), &mut chunk, args.chunk_lines)?;
        if n == 0 {
            break;
        }
        writer.write_blocks(format_into_blocks(&chunk, &dict, &transform))?;
        total_lines += n as u64;
    }
    writer.flush().context("failed to flush output")?;

    log::info!(
        "done: {} lines processed in {:.2?}",
        total_lines,
        start.elapsed()
    );
    Ok(())
}
