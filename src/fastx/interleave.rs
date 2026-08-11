//! `fastx interleave`: merge R1/R2 into a single standard interleaved
//! FASTQ (`R1,R2,R1,R2,...`).
//!
//! Unlike `fastx deinterleave`, this direction has no ambiguity to detect:
//! the output layout is ours to choose, so it's always proper interleaving
//! -- the convention tools like `bwa mem -p` expect. (A naive `cat R1 R2`
//! concatenation is not "interleaving" in any useful sense and is already
//! one shell command away if that's what's wanted; see `fastx cat` if the
//! goal is merging multiple *same-mate* files from repeated sequencing
//! runs, which is a different problem with its own duplicate-ID check.)
//!
//! Reuses `spawn_reader`/`recv_pair_step` for concurrent, pairing-checked
//! reading (same infrastructure and `--no-pair-check` escape hatch as
//! `fastx sample`/`fastx rescue`), and `format_into_blocks` for parallel
//! gzip output: each buffered chunk of pairs is flattened into
//! `[r1, r2, r1, r2, ...]` before splitting into compression blocks, so
//! block boundaries never need to respect pair boundaries -- decompression
//! just concatenates all blocks back into one continuous, correctly
//! ordered stream regardless of where they were cut.

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{bail, Result};
use clap::Args;

use crate::common::fastq::{
    format_into_blocks, recv_pair_step, spawn_reader, FastqRecord, PairStep,
};
use crate::io_utils::open_block_writer;

#[derive(Args, Debug)]
pub struct InterleaveArgs {
    /// Read 1 (R1) FASTQ (.fastq/.fq, gzip/bgzip auto-detected).
    #[arg(short = 'i', long = "in1", value_name = "FILE")]
    in1: PathBuf,

    /// Read 2 (R2) mate file.
    #[arg(short = 'I', long = "in2", value_name = "FILE")]
    in2: PathBuf,

    /// Output interleaved FASTQ. Gzip-compressed if the path ends in `.gz`.
    #[arg(short = 'o', long = "out", value_name = "FILE")]
    output: PathBuf,

    /// Skip verifying that R1/R2 read IDs correspond to the same pair at
    /// each position. Only disable this if your headers don't follow
    /// standard `/1`+`/2` or Illumina `1:...`+`2:...` mate-suffix conventions.
    #[arg(long)]
    no_pair_check: bool,

    /// Pairs processed per parallel compression batch.
    #[arg(long, default_value_t = 50_000)]
    chunk_records: usize,
}

pub fn run(args: InterleaveArgs) -> Result<()> {
    let start = Instant::now();
    log::info!(
        "interleaving {} + {} -> {}",
        args.in1.display(),
        args.in2.display(),
        args.output.display()
    );

    let rx1 = spawn_reader(args.in1.clone())?;
    let rx2 = spawn_reader(args.in2.clone())?;
    let mut writer = open_block_writer(&args.output)?;
    let check_ids = !args.no_pair_check;

    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(args.chunk_records * 2);
    let mut total: u64 = 0;
    loop {
        chunk.clear();
        for _ in 0..args.chunk_records {
            match recv_pair_step(&rx1, &rx2) {
                PairStep::Pair { r1, r2, ids_match } => {
                    if check_ids && !ids_match {
                        bail!(
                            "read 1/2 desync at pair #{total}: IDs {:?} vs {:?} do not match; \
                             files are not properly paired (pass --no-pair-check to override)",
                            r1.base_id(),
                            r2.base_id()
                        );
                    }
                    chunk.push(r1);
                    chunk.push(r2);
                }
                PairStep::Eof => break,
                PairStep::CountMismatch => bail!(
                    "read 1/2 have different numbers of reads (mismatch detected at pair #{total}); \
                     files are not properly paired"
                ),
                PairStep::ReadError(e) => return Err(e),
            }
        }
        if chunk.is_empty() {
            break;
        }
        let refs: Vec<&FastqRecord> = chunk.iter().collect();
        writer.write_blocks(format_into_blocks(&refs)?)?;
        total += (chunk.len() / 2) as u64;
    }
    writer.flush()?;

    log::info!("interleaved {total} read pairs in {:.2?}", start.elapsed());
    Ok(())
}
