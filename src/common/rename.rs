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

use crate::io_utils::{open_reader, open_writer, read_line_chunk, write_lines};

#[derive(Args, Debug)]
pub struct RenameCommonArgs {
    /// Input file. Gzip/bgzip is auto-detected regardless of extension.
    #[arg(value_name = "FILE")]
    pub input: PathBuf,

    /// Mapping file: two whitespace-separated columns per line, `new_name old_name`.
    /// May itself be gzip compressed.
    #[arg(short = 'n', long)]
    pub name: PathBuf,

    /// Output file. Written gzip-compressed if the path ends in `.gz`.
    #[arg(short = 'o', long)]
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
        let new_name = fields
            .next()
            .with_context(|| format!("name file line {}: expected 2 columns, got: {:?}", i + 1, line))?;
        let old_name = fields
            .next()
            .with_context(|| format!("name file line {}: expected 2 columns, got: {:?}", i + 1, line))?;
        map.insert(old_name.to_string(), new_name.to_string());
    }
    Ok(map)
}

/// Streams `args.input` to `args.output`, applying `transform` to every line
/// in parallel batches of `args.chunk_lines`. `transform` receives the
/// trimmed line and the loaded name dictionary.
pub fn run(args: &RenameCommonArgs, transform: impl Fn(&str, &HashMap<String, String>) -> String + Sync) -> Result<()> {
    ensure!(args.chunk_lines > 0, "--chunk-lines must be > 0");

    let start = Instant::now();
    log::info!("loading name mapping from {}", args.name.display());
    let dict = load_name_dict(&args.name)?;
    log::info!("loaded {} name mappings", dict.len());

    log::info!("processing {} -> {}", args.input.display(), args.output.display());

    let mut reader = open_reader(&args.input)?;
    let mut writer = open_writer(&args.output)?;

    let mut chunk = Vec::with_capacity(args.chunk_lines);
    let mut total_lines: u64 = 0;
    loop {
        let n = read_line_chunk(reader.as_mut(), &mut chunk, args.chunk_lines)?;
        if n == 0 {
            break;
        }
        let out: Vec<String> = chunk.par_iter().map(|l| transform(l, &dict)).collect();
        write_lines(writer.as_mut(), &out)?;
        total_lines += n as u64;
    }
    writer.flush().context("failed to flush output")?;

    log::info!("done: {} lines processed in {:.2?}", total_lines, start.elapsed());
    Ok(())
}
