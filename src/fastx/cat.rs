//! `fastx cat`: concatenate FASTQ files from repeated sequencing runs of
//! the same biological sample, checking for duplicate read IDs along the way.
//!
//! Concatenating raw FASTQ from multiple lanes/flowcells of one sample is
//! standard practice and normally safe: real Illumina read IDs encode
//! flowcell/lane/tile/coordinates, so genuine cross-run ID collisions are
//! not expected. The realistic failure mode isn't the ID scheme -- it's
//! operator error: the same file accidentally included twice in the file
//! list (a typo'd path, a copy-pasted glob that matched more than
//! intended). That mistake is completely silent if nothing checks for it,
//! and it quietly inflates coverage/duplicates data going into downstream
//! alignment or variant calling. This command hashes every read ID
//! (`common::hash::fnv1a`) as it streams through, and aborts with the
//! specific source files and record positions involved the moment a
//! duplicate shows up, rather than concatenating first and leaving the
//! problem for something else to notice later (or never).
//!
//! Paired-end mode additionally re-checks R1/R2 pairing *within* each
//! source file pair as it streams through (reusing `recv_pair_step`,
//! the same mechanism `fastx sample`/`fastx rescue` use) -- catching a
//! corrupt or mismatched individual run, not just cross-run duplicates.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, ensure, Result};
use clap::Args;

use crate::common::fastq::{
    format_into_blocks, recv_pair_step, spawn_reader, FastqRecord, PairStep,
};
use crate::common::hash::fnv1a;
use crate::io_utils::{open_block_writer, BlockWriter};

#[derive(Args, Debug)]
pub struct CatArgs {
    /// Read 1 (or single-end) input file. Repeat once per source run, in
    /// the order they should be concatenated.
    #[arg(long = "r1", value_name = "FILE", required = true)]
    r1: Vec<PathBuf>,

    /// Read 2 mate file, one per --r1 in the same order. Presence switches
    /// to paired-end mode.
    #[arg(long = "r2", value_name = "FILE")]
    r2: Vec<PathBuf>,

    /// Output for read 1 / single-end reads. Gzip-compressed if the path ends in `.gz`.
    #[arg(short = 'o', long = "out1", value_name = "FILE")]
    out1: PathBuf,

    /// Output for read 2. Required when --r2 is given.
    #[arg(short = 'O', long = "out2", value_name = "FILE")]
    out2: Option<PathBuf>,

    /// Don't check for duplicate read IDs across the concatenated files.
    /// Only disable this if you're confident the inputs are genuinely
    /// distinct (e.g. a platform that doesn't guarantee globally unique IDs).
    #[arg(long)]
    allow_duplicate_ids: bool,

    /// Records processed per parallel compression batch.
    #[arg(long, default_value_t = 50_000)]
    chunk_records: usize,
}

pub fn run(args: CatArgs) -> Result<()> {
    ensure!(!args.r1.is_empty(), "at least one --r1 input is required");
    if args.r2.is_empty() {
        ensure!(
            args.out2.is_none(),
            "-O/--out2 was given but no --r2 inputs were provided"
        );
        run_se(&args)
    } else {
        ensure!(
            args.r2.len() == args.r1.len(),
            "--r1 and --r2 must be given the same number of times ({} vs {})",
            args.r1.len(),
            args.r2.len()
        );
        ensure!(
            args.out2.is_some(),
            "-O/--out2 is required when --r2 is given"
        );
        run_pe(&args)
    }
}

/// Records the source file + local record index where a hash was first
/// seen, so a duplicate hit can report exactly what to go check.
type SeenIds = HashMap<u64, (PathBuf, u64)>;

fn check_duplicate(seen: &mut SeenIds, id: &str, path: &Path, local_idx: u64) -> Result<()> {
    let h = fnv1a(id.as_bytes());
    if let Some((prev_path, prev_idx)) = seen.insert(h, (path.to_path_buf(), local_idx)) {
        bail!(
            "duplicate read ID {id:?}: first seen in {} (record #{prev_idx}), again in {} (record #{local_idx}) -- \
             did you accidentally include the same file twice? pass --allow-duplicate-ids to skip this check",
            prev_path.display(),
            path.display()
        );
    }
    Ok(())
}

fn run_se(args: &CatArgs) -> Result<()> {
    let start = Instant::now();
    let mut writer = open_block_writer(&args.out1)?;
    let mut seen = SeenIds::new();
    let check_ids = !args.allow_duplicate_ids;
    let mut total: u64 = 0;

    for (file_idx, path) in args.r1.iter().enumerate() {
        log::info!(
            "cat-ing source {}/{}: {}",
            file_idx + 1,
            args.r1.len(),
            path.display()
        );
        let rx = spawn_reader(path.clone())?;
        let mut local_idx: u64 = 0;
        let mut chunk: Vec<FastqRecord> = Vec::with_capacity(args.chunk_records);
        loop {
            chunk.clear();
            for r in rx.iter().take(args.chunk_records) {
                let rec = r?;
                if check_ids {
                    check_duplicate(&mut seen, rec.base_id(), path, local_idx)?;
                }
                chunk.push(rec);
                local_idx += 1;
            }
            if chunk.is_empty() {
                break;
            }
            let refs: Vec<&FastqRecord> = chunk.iter().collect();
            writer.write_blocks(format_into_blocks(&refs)?)?;
        }
        total += local_idx;
    }
    writer.flush()?;

    log::info!(
        "concatenated {} source file(s), {total} read(s) total, in {:.2?}",
        args.r1.len(),
        start.elapsed()
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn cat_one_pe_source(
    r1_path: &Path,
    r2_path: &Path,
    w1: &mut BlockWriter,
    w2: &mut BlockWriter,
    seen: &mut SeenIds,
    check_ids: bool,
    chunk_records: usize,
) -> Result<u64> {
    let rx1 = spawn_reader(r1_path.to_path_buf())?;
    let rx2 = spawn_reader(r2_path.to_path_buf())?;
    let mut local_idx: u64 = 0;
    let mut chunk: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(chunk_records);
    loop {
        chunk.clear();
        for _ in 0..chunk_records {
            match recv_pair_step(&rx1, &rx2) {
                PairStep::Pair { r1, r2, ids_match } => {
                    if !ids_match {
                        bail!(
                            "within source pair {} + {}, read 1/2 desync at local pair #{local_idx}: \
                             IDs {:?} vs {:?} do not match",
                            r1_path.display(),
                            r2_path.display(),
                            r1.base_id(),
                            r2.base_id()
                        );
                    }
                    if check_ids {
                        check_duplicate(seen, r1.base_id(), r1_path, local_idx)?;
                    }
                    chunk.push((r1, r2));
                    local_idx += 1;
                }
                PairStep::Eof => break,
                PairStep::CountMismatch => bail!(
                    "source pair {} + {} have different numbers of reads",
                    r1_path.display(),
                    r2_path.display()
                ),
                PairStep::ReadError(e) => return Err(e),
            }
        }
        if chunk.is_empty() {
            break;
        }
        let r1s: Vec<&FastqRecord> = chunk.iter().map(|(a, _)| a).collect();
        let r2s: Vec<&FastqRecord> = chunk.iter().map(|(_, b)| b).collect();
        w1.write_blocks(format_into_blocks(&r1s)?)?;
        w2.write_blocks(format_into_blocks(&r2s)?)?;
    }
    Ok(local_idx)
}

fn run_pe(args: &CatArgs) -> Result<()> {
    let start = Instant::now();
    let mut w1 = open_block_writer(&args.out1)?;
    let mut w2 = open_block_writer(args.out2.as_ref().expect("checked by run()"))?;
    let mut seen = SeenIds::new();
    let check_ids = !args.allow_duplicate_ids;
    let mut total: u64 = 0;

    for (file_idx, (r1_path, r2_path)) in args.r1.iter().zip(args.r2.iter()).enumerate() {
        log::info!(
            "cat-ing source {}/{}: {} + {}",
            file_idx + 1,
            args.r1.len(),
            r1_path.display(),
            r2_path.display()
        );
        total += cat_one_pe_source(
            r1_path,
            r2_path,
            &mut w1,
            &mut w2,
            &mut seen,
            check_ids,
            args.chunk_records,
        )?;
    }
    w1.flush()?;
    w2.flush()?;

    log::info!(
        "concatenated {} source pair(s), {total} read pair(s) total, in {:.2?}",
        args.r1.len(),
        start.elapsed()
    );
    Ok(())
}
