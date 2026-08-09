//! `fastx sample`: downsample FASTQ reads by proportion or exact count.
//!
//! Unlike `seqkit sample`, which has no notion of paired-end input (you run
//! it twice, once per mate file, and rely on passing the identical `-s` seed
//! both times so the two independent runs *happen* to make the same
//! keep/discard decision at every read index), this command reads both
//! mates in one process and makes a single decision per read *pair*. Pairing
//! is therefore correct by construction rather than by convention, and a
//! desynced R1/R2 (different read counts, or mismatched IDs) is caught and
//! reported instead of silently producing broken pairs.
//!
//! `-p` (proportion) sampling computes each record's keep/discard decision
//! as a pure function of (seed, global index) via `common::rng`, so batches
//! can be evaluated in parallel with rayon and give the same result
//! regardless of batch size or thread count. `-n` (exact count) uses
//! single-pass reservoir sampling (Algorithm R): unlike `seqkit sample -n`,
//! which either loads the entire input into memory or re-reads the file
//! twice, memory use here is O(n) in the sample size, not the input size,
//! in a single pass.

use std::io::Write as _;
use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver};
use std::thread;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use anyhow::{bail, ensure, Result};
use clap::Args;
use rayon::prelude::*;

use crate::common::fastq::{read_fastq_record, FastqRecord};
use crate::common::rng::{deterministic_f64, SplitMix64};
use crate::io_utils::{open_reader, open_writer};

#[derive(Args, Debug)]
pub struct SampleArgs {
    /// Input FASTQ (.fastq/.fq, gzip/bgzip auto-detected). Read 1 (R1) in paired-end mode.
    #[arg(short = 'i', long = "in1", value_name = "FILE")]
    in1: PathBuf,

    /// Read 2 (R2) mate file. Presence of this flag switches to paired-end mode.
    #[arg(short = 'I', long = "in2", value_name = "FILE")]
    in2: Option<PathBuf>,

    /// Output for read 1 / single-end reads. Gzip-compressed if the path ends in `.gz`.
    #[arg(short = 'o', long = "out1", value_name = "FILE")]
    out1: PathBuf,

    /// Output for read 2. Required when -I/--in2 is given.
    #[arg(short = 'O', long = "out2", value_name = "FILE", requires = "in2")]
    out2: Option<PathBuf>,

    /// Sample by proportion in (0, 1], e.g. 0.1 for 10%. Single pass, O(1) memory,
    /// parallelized across --chunk-records batches. Mutually exclusive with -n.
    #[arg(short = 'p', long, conflicts_with = "number")]
    proportion: Option<f64>,

    /// Sample an exact number of reads (read pairs, in paired-end mode) via
    /// single-pass reservoir sampling. Memory is O(number), not O(total reads).
    #[arg(short = 'n', long, conflicts_with = "proportion")]
    number: Option<u64>,

    /// Random seed for reproducible sampling. In paired-end mode both mates are
    /// drawn from one decision stream keyed by pair index, so pairing is always
    /// consistent regardless of seed — the seed only controls which pairs are picked.
    #[arg(short = 's', long = "seed", default_value_t = 42)]
    seed: u64,

    /// Use a time-based seed instead of -s/--seed (non-reproducible). Takes priority over -s.
    #[arg(short = 'r', long)]
    non_deterministic: bool,

    /// Skip verifying that R1/R2 read IDs correspond to the same pair at each
    /// position. Only disable this if your headers don't follow standard
    /// `/1`+`/2` or Illumina `1:...`+`2:...` mate-suffix conventions.
    #[arg(long)]
    no_pair_check: bool,

    /// Records (or pairs) processed per parallel batch in proportion mode.
    #[arg(long, default_value_t = 50_000)]
    chunk_records: usize,
}

fn effective_seed(args: &SampleArgs) -> u64 {
    if args.non_deterministic {
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before UNIX epoch")
            .as_nanos() as u64
    } else {
        args.seed
    }
}

/// Spawns a background thread that decompresses and parses `path`, streaming
/// parsed records out through a bounded channel. Running this on its own
/// thread is what lets R1 and R2 decompress concurrently in paired-end mode
/// instead of sequentially, as two separate `seqkit sample` invocations would.
fn spawn_reader(path: PathBuf) -> Result<Receiver<Result<FastqRecord>>> {
    let mut reader = open_reader(&path)?;
    let (tx, rx) = mpsc::sync_channel::<Result<FastqRecord>>(4096);
    thread::spawn(move || {
        let mut line_no: u64 = 1;
        loop {
            match read_fastq_record(reader.as_mut(), line_no) {
                Ok(Some(rec)) => {
                    line_no += 4;
                    if tx.send(Ok(rec)).is_err() {
                        break;
                    }
                }
                Ok(None) => break,
                Err(e) => {
                    let _ = tx.send(Err(e));
                    break;
                }
            }
        }
    });
    Ok(rx)
}

pub fn run(args: SampleArgs) -> Result<()> {
    ensure!(
        args.proportion.is_some() || args.number.is_some(),
        "must specify one of -p/--proportion or -n/--number"
    );
    if let Some(p) = args.proportion {
        ensure!(p > 0.0 && p <= 1.0, "-p/--proportion must be in (0, 1], got {p}");
    }
    if let Some(n) = args.number {
        ensure!(n > 0, "-n/--number must be > 0");
    }
    ensure!(args.chunk_records > 0, "--chunk-records must be > 0");
    if args.in2.is_some() {
        ensure!(args.out2.is_some(), "-O/--out2 is required when -I/--in2 is given");
    }

    let seed = effective_seed(&args);
    log::info!(
        "random seed: {seed}{}",
        if args.non_deterministic { " (time-based)" } else { "" }
    );

    match (&args.in2, &args.out2) {
        (Some(in2), Some(out2)) => run_pe(&args, seed, in2, out2),
        _ => run_se(&args, seed),
    }
}

fn run_se(args: &SampleArgs, seed: u64) -> Result<()> {
    let start = Instant::now();
    log::info!("sampling {} -> {}", args.in1.display(), args.out1.display());
    let rx = spawn_reader(args.in1.clone())?;
    let mut writer = open_writer(&args.out1)?;

    let (total, kept) = if let Some(p) = args.proportion {
        run_proportion_se(&rx, writer.as_mut(), seed, p, args.chunk_records)?
    } else {
        run_reservoir_se(&rx, writer.as_mut(), seed, args.number.unwrap())?
    };
    writer.flush()?;

    log::info!(
        "sampled {kept}/{total} reads ({:.2}%) in {:.2?}",
        100.0 * kept as f64 / total.max(1) as f64,
        start.elapsed()
    );
    Ok(())
}

fn run_proportion_se(
    rx: &Receiver<Result<FastqRecord>>,
    writer: &mut dyn std::io::Write,
    seed: u64,
    p: f64,
    chunk_size: usize,
) -> Result<(u64, u64)> {
    let mut total: u64 = 0;
    let mut kept: u64 = 0;
    let mut buf: Vec<FastqRecord> = Vec::with_capacity(chunk_size);
    loop {
        buf.clear();
        for r in rx.iter().take(chunk_size) {
            buf.push(r?);
        }
        if buf.is_empty() {
            break;
        }
        let base_idx = total;
        let selected: Vec<&FastqRecord> = buf
            .par_iter()
            .enumerate()
            .filter(|(i, _)| deterministic_f64(seed, base_idx + *i as u64) <= p)
            .map(|(_, r)| r)
            .collect();
        for r in &selected {
            r.write_to(writer)?;
        }
        kept += selected.len() as u64;
        total += buf.len() as u64;
    }
    Ok((total, kept))
}

fn run_reservoir_se(rx: &Receiver<Result<FastqRecord>>, writer: &mut dyn std::io::Write, seed: u64, n: u64) -> Result<(u64, u64)> {
    let n = n as usize;
    let mut rng = SplitMix64::new(seed);
    let mut reservoir: Vec<Option<(u64, FastqRecord)>> = vec![None; n];
    let mut total: u64 = 0;

    for r in rx.iter() {
        let rec = r?;
        let i = total;
        if (i as usize) < n {
            reservoir[i as usize] = Some((i, rec));
        } else {
            let j = rng.next_below(i + 1);
            if (j as usize) < n {
                reservoir[j as usize] = Some((i, rec));
            }
        }
        total += 1;
    }

    let mut chosen: Vec<(u64, FastqRecord)> = reservoir.into_iter().flatten().collect();
    chosen.sort_by_key(|(idx, _)| *idx);
    for (_, rec) in &chosen {
        rec.write_to(writer)?;
    }
    let kept = chosen.len() as u64;
    Ok((total, kept))
}

/// Receives one record from each of the two mate channels and checks they
/// stay in lockstep: both ending at the same time (clean EOF), and — unless
/// disabled — both records referring to the same underlying read pair.
fn recv_pair(
    rx1: &Receiver<Result<FastqRecord>>,
    rx2: &Receiver<Result<FastqRecord>>,
    idx: u64,
    check_ids: bool,
) -> Result<Option<(FastqRecord, FastqRecord)>> {
    match (rx1.recv(), rx2.recv()) {
        (Err(_), Err(_)) => Ok(None),
        (Ok(r1), Ok(r2)) => {
            let r1 = r1?;
            let r2 = r2?;
            if check_ids && r1.base_id() != r2.base_id() {
                bail!(
                    "read 1/2 desync at pair #{idx}: IDs {:?} vs {:?} do not match; \
                     files are not properly paired (pass --no-pair-check to override)",
                    r1.base_id(),
                    r2.base_id()
                );
            }
            Ok(Some((r1, r2)))
        }
        _ => bail!(
            "read 1/2 have different numbers of reads (mismatch detected at pair #{idx}); \
             files are not properly paired"
        ),
    }
}

fn run_pe(args: &SampleArgs, seed: u64, in2: &std::path::Path, out2: &std::path::Path) -> Result<()> {
    let start = Instant::now();
    log::info!(
        "sampling paired-end {} + {} -> {} + {}",
        args.in1.display(),
        in2.display(),
        args.out1.display(),
        out2.display()
    );
    let rx1 = spawn_reader(args.in1.clone())?;
    let rx2 = spawn_reader(in2.to_path_buf())?;
    let mut w1 = open_writer(&args.out1)?;
    let mut w2 = open_writer(out2)?;
    let check_ids = !args.no_pair_check;

    let (total, kept) = if let Some(p) = args.proportion {
        run_proportion_pe(&rx1, &rx2, w1.as_mut(), w2.as_mut(), seed, p, args.chunk_records, check_ids)?
    } else {
        run_reservoir_pe(&rx1, &rx2, w1.as_mut(), w2.as_mut(), seed, args.number.unwrap(), check_ids)?
    };
    w1.flush()?;
    w2.flush()?;

    log::info!(
        "sampled {kept}/{total} read pairs ({:.2}%) in {:.2?}",
        100.0 * kept as f64 / total.max(1) as f64,
        start.elapsed()
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_proportion_pe(
    rx1: &Receiver<Result<FastqRecord>>,
    rx2: &Receiver<Result<FastqRecord>>,
    w1: &mut dyn std::io::Write,
    w2: &mut dyn std::io::Write,
    seed: u64,
    p: f64,
    chunk_size: usize,
    check_ids: bool,
) -> Result<(u64, u64)> {
    let mut total: u64 = 0;
    let mut kept: u64 = 0;
    let mut buf: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(chunk_size);
    loop {
        buf.clear();
        for _ in 0..chunk_size {
            match recv_pair(rx1, rx2, total + buf.len() as u64, check_ids)? {
                Some(pair) => buf.push(pair),
                None => break,
            }
        }
        if buf.is_empty() {
            break;
        }
        let base_idx = total;
        let selected: Vec<&(FastqRecord, FastqRecord)> = buf
            .par_iter()
            .enumerate()
            .filter(|(i, _)| deterministic_f64(seed, base_idx + *i as u64) <= p)
            .map(|(_, pair)| pair)
            .collect();
        for (r1, r2) in &selected {
            r1.write_to(w1)?;
            r2.write_to(w2)?;
        }
        kept += selected.len() as u64;
        total += buf.len() as u64;
    }
    Ok((total, kept))
}

fn run_reservoir_pe(
    rx1: &Receiver<Result<FastqRecord>>,
    rx2: &Receiver<Result<FastqRecord>>,
    w1: &mut dyn std::io::Write,
    w2: &mut dyn std::io::Write,
    seed: u64,
    n: u64,
    check_ids: bool,
) -> Result<(u64, u64)> {
    let n = n as usize;
    let mut rng = SplitMix64::new(seed);
    let mut reservoir: Vec<Option<(u64, FastqRecord, FastqRecord)>> = vec![None; n];
    let mut total: u64 = 0;

    while let Some((r1, r2)) = recv_pair(rx1, rx2, total, check_ids)? {
        let i = total;
        if (i as usize) < n {
            reservoir[i as usize] = Some((i, r1, r2));
        } else {
            let j = rng.next_below(i + 1);
            if (j as usize) < n {
                reservoir[j as usize] = Some((i, r1, r2));
            }
        }
        total += 1;
    }

    let mut chosen: Vec<(u64, FastqRecord, FastqRecord)> = reservoir.into_iter().flatten().collect();
    chosen.sort_by_key(|(idx, _, _)| *idx);
    for (_, r1, r2) in &chosen {
        r1.write_to(w1)?;
        r2.write_to(w2)?;
    }
    let kept = chosen.len() as u64;
    Ok((total, kept))
}
