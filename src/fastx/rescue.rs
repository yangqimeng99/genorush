//! `fastx rescue`: recover the leading run of clean, well-formed FASTQ
//! records from a truncated or corrupted file (or pair of mate files),
//! discarding only the wreckage from the point of failure onward.
//!
//! Interrupted downloads are the common case: a `.fq.gz` that stops midway
//! through leaves a gzip stream with a truncated final deflate block (or a
//! bad trailing CRC). `flate2` surfaces that as a read error partway through
//! decompression rather than a clean EOF -- everything decoded *before* that
//! error is still perfectly good data, and for paired-end reads, any pair
//! where both mates were fully decoded and still line up is still usable
//! for alignment. This command reuses the exact concurrent-reader and
//! pairing infrastructure from `fastx sample`
//! (`common::fastq::{spawn_reader, recv_pair_step}`), but where `sample`
//! treats a read error, a count mismatch, or an ID mismatch as fatal, this
//! command treats each of them as simply "the point where usable data
//! ends" -- it stops there, keeps everything decoded up to that point, and
//! reports what it recovered instead of aborting.
//!
//! Confirmed-good records are buffered in `--chunk-records`-sized batches
//! and handed to the parallel-compressing `BlockWriter` (`io_utils`,
//! `common::fastq::format_into_blocks`) exactly like `fastx sample`/`fastx
//! cat` -- rescue just also flushes whatever partial batch is left over
//! the moment corruption is hit or the input ends, since that data needs
//! to make it to disk regardless of where it fell relative to a chunk
//! boundary.
//!
//! This only catches *structural* corruption (truncation, bad gzip, a
//! header/plus-line/seq-qual-length violation). Bit-level corruption that
//! leaves a record structurally well-formed but with wrong sequence/quality
//! content is undetectable here -- FASTQ carries no per-record checksum.

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{ensure, Result};
use clap::Args;

use crate::common::fastq::{
    format_into_blocks, recv_pair_step, spawn_reader, FastqRecord, PairStep,
};
use crate::io_utils::{open_block_writer, BlockWriter};

#[derive(Args, Debug)]
pub struct RescueArgs {
    /// Input FASTQ (.fastq/.fq, gzip/bgzip auto-detected), possibly truncated/corrupted.
    /// Read 1 (R1) in paired-end mode.
    #[arg(short = 'i', long = "in1", value_name = "FILE")]
    in1: PathBuf,

    /// Read 2 (R2) mate file. Presence of this flag switches to paired-end mode.
    #[arg(short = 'I', long = "in2", value_name = "FILE")]
    in2: Option<PathBuf>,

    /// Output for the rescued read 1 / single-end reads. Gzip-compressed if
    /// the path ends in `.gz`.
    #[arg(short = 'o', long = "out1", value_name = "FILE")]
    out1: PathBuf,

    /// Output for the rescued read 2. Required when -I/--in2 is given.
    #[arg(short = 'O', long = "out2", value_name = "FILE", requires = "in2")]
    out2: Option<PathBuf>,

    /// Don't stop at the first R1/R2 ID mismatch; only stop on a read error
    /// or a read-count mismatch. Only use this if your headers don't follow
    /// standard `/1`+`/2` or Illumina `1:...`+`2:...` mate-suffix conventions
    /// (otherwise an ID mismatch is usually itself the corruption point).
    #[arg(long)]
    no_pair_check: bool,

    /// Records processed per parallel compression batch.
    #[arg(long, default_value_t = 50_000)]
    chunk_records: usize,
}

pub fn run(args: RescueArgs) -> Result<()> {
    if args.in2.is_some() {
        ensure!(
            args.out2.is_some(),
            "-O/--out2 is required when -I/--in2 is given"
        );
    }

    match (&args.in2, &args.out2) {
        (Some(in2), Some(out2)) => run_pe(&args, in2, out2),
        _ => run_se(&args),
    }
}

/// Exits the process with a code reflecting the outcome, since a partial
/// rescue is a distinct, script-detectable outcome from a clean read and
/// from an unusable input: 0 = read cleanly end to end, 3 = corruption
/// found but something was rescued, 1 = corruption found and nothing was
/// salvageable.
fn finish(corrupted: bool, rescued: u64) -> Result<()> {
    if corrupted && rescued == 0 {
        std::process::exit(1);
    } else if corrupted {
        std::process::exit(3);
    }
    Ok(())
}

fn run_se(args: &RescueArgs) -> Result<()> {
    let start = Instant::now();
    log::info!("rescuing {} -> {}", args.in1.display(), args.out1.display());
    let rx = spawn_reader(args.in1.clone())?;
    let mut writer = open_block_writer(&args.out1)?;

    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(args.chunk_records);
    let mut rescued: u64 = 0;
    let mut corrupted = false;
    loop {
        match rx.recv() {
            Ok(Ok(rec)) => {
                chunk.push(rec);
                rescued += 1;
                if chunk.len() >= args.chunk_records {
                    let refs: Vec<&FastqRecord> = chunk.iter().collect();
                    writer.write_blocks(format_into_blocks(&refs)?)?;
                    chunk.clear();
                }
            }
            Ok(Err(e)) => {
                log::warn!("stopping at read #{}: {e}", rescued + 1);
                corrupted = true;
                break;
            }
            Err(_) => break, // clean, synchronized EOF
        }
    }
    if !chunk.is_empty() {
        let refs: Vec<&FastqRecord> = chunk.iter().collect();
        writer.write_blocks(format_into_blocks(&refs)?)?;
    }
    writer.flush()?;

    if corrupted {
        log::warn!("input looks truncated/corrupted: rescued {rescued} clean read(s) before the failure point");
    } else {
        log::info!("input read cleanly: {rescued} read(s), no corruption detected");
    }
    log::info!("done in {:.2?}", start.elapsed());

    finish(corrupted, rescued)
}

#[allow(clippy::too_many_arguments)]
fn flush_pe_chunk(
    chunk: &[(FastqRecord, FastqRecord)],
    w1: &mut BlockWriter,
    w2: &mut BlockWriter,
) -> Result<()> {
    let r1s: Vec<&FastqRecord> = chunk.iter().map(|(a, _)| a).collect();
    let r2s: Vec<&FastqRecord> = chunk.iter().map(|(_, b)| b).collect();
    w1.write_blocks(format_into_blocks(&r1s)?)?;
    w2.write_blocks(format_into_blocks(&r2s)?)?;
    Ok(())
}

fn run_pe(args: &RescueArgs, in2: &Path, out2: &Path) -> Result<()> {
    let start = Instant::now();
    log::info!(
        "rescuing paired-end {} + {} -> {} + {}",
        args.in1.display(),
        in2.display(),
        args.out1.display(),
        out2.display()
    );
    let rx1 = spawn_reader(args.in1.clone())?;
    let rx2 = spawn_reader(in2.to_path_buf())?;
    let mut w1 = open_block_writer(&args.out1)?;
    let mut w2 = open_block_writer(out2)?;
    let check_ids = !args.no_pair_check;

    let mut chunk: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(args.chunk_records);
    let mut rescued: u64 = 0;
    let mut corrupted = false;
    loop {
        match recv_pair_step(&rx1, &rx2) {
            PairStep::Pair { r1, r2, ids_match } => {
                if check_ids && !ids_match {
                    log::warn!(
                        "stopping at pair #{}: R1/R2 IDs diverge ({:?} vs {:?}) -- likely the desync point \
                         (pass --no-pair-check to keep going past ID mismatches)",
                        rescued + 1,
                        r1.base_id(),
                        r2.base_id()
                    );
                    corrupted = true;
                    break;
                }
                chunk.push((r1, r2));
                rescued += 1;
                if chunk.len() >= args.chunk_records {
                    flush_pe_chunk(&chunk, &mut w1, &mut w2)?;
                    chunk.clear();
                }
            }
            PairStep::Eof => break, // clean, synchronized EOF on both mates
            PairStep::CountMismatch => {
                log::warn!(
                    "stopping at pair #{}: R1 and R2 have a different number of reads from here on",
                    rescued + 1
                );
                corrupted = true;
                break;
            }
            PairStep::ReadError(e) => {
                log::warn!("stopping at pair #{}: {e}", rescued + 1);
                corrupted = true;
                break;
            }
        }
    }
    if !chunk.is_empty() {
        flush_pe_chunk(&chunk, &mut w1, &mut w2)?;
    }
    w1.flush()?;
    w2.flush()?;

    if corrupted {
        log::warn!("input looks truncated/corrupted: rescued {rescued} clean read pair(s) before the failure point");
    } else {
        log::info!("input read cleanly: {rescued} read pair(s), no corruption detected");
    }
    log::info!("done in {:.2?}", start.elapsed());

    finish(corrupted, rescued)
}
