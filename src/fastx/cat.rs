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

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, ensure, Result};
use clap::Args;

use crate::common::fastq::{
    format_into_blocks, recv_pair_step, spawn_reader, FastqRecord, PairStep,
};
use crate::common::hash::{fnv1a, BuildIdHasher};
use crate::io_utils::{
    display, ensure_one_stdio_at_most, is_stdio, open_block_writer, BlockWriter, OutputOpts, STDIO,
};

#[derive(Args, Debug)]
pub struct CatArgs {
    /// Read 1 (or single-end) input file. Repeat once per source run, in
    /// the order they should be concatenated. At most one input across
    /// --r1/--r2 may be `-`, since there is only one stdin.
    #[arg(long = "r1", value_name = "FILE", required = true)]
    r1: Vec<PathBuf>,

    /// Read 2 mate file, one per --r1 in the same order. Presence switches
    /// to paired-end mode.
    #[arg(long = "r2", value_name = "FILE")]
    r2: Vec<PathBuf>,

    /// Output for read 1 / single-end reads. Defaults to stdout.
    /// Gzip-compressed if the path ends in `.gz` or `-z/--gzip` is passed.
    #[arg(short = 'o', long = "out1", value_name = "FILE", default_value = STDIO)]
    out1: PathBuf,

    /// Output for read 2. Required when --r2 is given. Only one of the two
    /// outputs can be `-`: both would land in the same pipe.
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

pub fn run(args: CatArgs, opts: OutputOpts) -> Result<()> {
    ensure!(!args.r1.is_empty(), "at least one --r1 input is required");
    let inputs: Vec<&Path> = args
        .r1
        .iter()
        .chain(args.r2.iter())
        .map(|p| p.as_path())
        .collect();
    ensure_one_stdio_at_most(&inputs, "input")?;
    let mut outputs: Vec<&Path> = vec![args.out1.as_path()];
    outputs.extend(args.out2.as_deref());
    ensure_one_stdio_at_most(&outputs, "output")?;
    if args.r2.is_empty() {
        ensure!(
            args.out2.is_none(),
            "-O/--out2 was given but no --r2 inputs were provided"
        );
        run_se(&args, opts)
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
        run_pe(&args, opts)
    }
}

/// Every read ID seen so far, as `fnv1a` hashes and nothing else.
///
/// This set is the command's whole memory footprint, and it grows with the
/// input: one entry per read, held until the run ends. What it holds
/// therefore matters at a scale nothing else here does. An earlier version
/// mapped each hash to the source path and record index where it was first
/// seen, so that a duplicate could be reported precisely -- 154 bytes per
/// read once the per-record `PathBuf` allocation and the table's growth
/// spikes were counted, or 42 GB for a 30x bovine WGS sample. A bare
/// `HashSet` of hashes costs about 10 bytes per read, and the first-seen
/// position is recovered by re-reading the inputs at the moment a repeat is
/// found -- see `locate_first`.
type SeenIds = HashSet<u64, BuildIdHasher>;

/// Where `id` first appears across `sources`: which source, by position in
/// the list, and which record within it (0-based).
///
/// Re-reading the inputs is only ever done once a hash repeat has been found,
/// which on real data means the command is about to stop. Paying for it there
/// buys two things: the error can still name the first occurrence even though
/// only hashes are kept, and a 64-bit collision between two genuinely
/// different IDs becomes distinguishable from a real duplicate rather than
/// aborting a correct run.
///
/// `None` means the answer is unavailable -- one of the inputs is stdin and
/// cannot be read a second time -- or that the ID was not found at all.
fn locate_first(sources: &[PathBuf], id: &str) -> Result<Option<(usize, u64)>> {
    if sources.iter().any(|p| is_stdio(p)) {
        return Ok(None);
    }
    for (source_idx, path) in sources.iter().enumerate() {
        let rx = spawn_reader(path.clone())?;
        for (local_idx, r) in (0_u64..).zip(rx.iter()) {
            if r?.base_id() == id {
                return Ok(Some((source_idx, local_idx)));
            }
        }
    }
    Ok(None)
}

fn check_duplicate(
    seen: &mut SeenIds,
    sources: &[PathBuf],
    source_idx: usize,
    id: &str,
    local_idx: u64,
) -> Result<()> {
    if seen.insert(fnv1a(id.as_bytes())) {
        return Ok(());
    }
    match locate_first(sources, id)? {
        // The only record carrying this ID is the one in hand, so nothing was
        // actually repeated: two different IDs hashed to the same u64.
        // Astronomically rare, and not a reason to fail a good run.
        Some(first) if first == (source_idx, local_idx) => {
            log::debug!(
                "hash collision on read ID {id:?} at record #{}; not a duplicate",
                local_idx + 1
            );
            Ok(())
        }
        Some((first_src, first_idx)) => bail!(
            "duplicate read ID {id:?}: first seen in {} (record #{}), again in {} (record #{}) -- \
             did you accidentally include the same file twice? pass --allow-duplicate-ids to skip this check",
            display(&sources[first_src]),
            first_idx + 1,
            display(&sources[source_idx]),
            local_idx + 1
        ),
        None => bail!(
            "duplicate read ID {id:?}: seen again in {} (record #{}), having already appeared \
             earlier in the inputs -- the earlier position can't be reported because stdin \
             cannot be re-read. Did you accidentally include the same file twice? \
             pass --allow-duplicate-ids to skip this check",
            display(&sources[source_idx]),
            local_idx + 1
        ),
    }
}

fn run_se(args: &CatArgs, opts: OutputOpts) -> Result<()> {
    let start = Instant::now();
    let mut writer = open_block_writer(&args.out1, opts)?;
    let mut seen = SeenIds::default();
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
                    check_duplicate(&mut seen, &args.r1, file_idx, rec.base_id(), local_idx)?;
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
    sources: &[PathBuf],
    source_idx: usize,
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
                            "within source pair {} + {}, read 1/2 desync at local pair #{}: \
                             IDs {:?} vs {:?} do not match",
                            r1_path.display(),
                            r2_path.display(),
                            local_idx + 1,
                            r1.base_id(),
                            r2.base_id()
                        );
                    }
                    if check_ids {
                        check_duplicate(seen, sources, source_idx, r1.base_id(), local_idx)?;
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

fn run_pe(args: &CatArgs, opts: OutputOpts) -> Result<()> {
    let start = Instant::now();
    let mut w1 = open_block_writer(&args.out1, opts)?;
    let mut w2 = open_block_writer(args.out2.as_ref().expect("checked by run()"), opts)?;
    let mut seen = SeenIds::default();
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
            &args.r1,
            file_idx,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::testutil::{headers_of, scratch_path as scratch, write_fastq};

    /// Builds the arguments a command line would have produced. `chunk_records`
    /// is deliberately tiny so multi-chunk paths are exercised by fixtures of
    /// a handful of records.
    fn args(
        r1: Vec<PathBuf>,
        r2: Vec<PathBuf>,
        out1: PathBuf,
        out2: Option<PathBuf>,
        allow_duplicate_ids: bool,
    ) -> CatArgs {
        CatArgs {
            r1,
            r2,
            out1,
            out2,
            allow_duplicate_ids,
            chunk_records: 3,
        }
    }

    fn source(name: &str, headers: &[&str]) -> PathBuf {
        let p = scratch(name);
        write_fastq(&p, headers);
        p
    }

    // ---- single-end ----

    #[test]
    fn distinct_sources_concatenate_in_order() {
        let a = source("run1", &["@a1/1", "@a2/1", "@a3/1", "@a4/1"]);
        let b = source("run2", &["@b1/1", "@b2/1"]);
        let out = scratch("out");

        run(
            args(vec![a, b], vec![], out.clone(), None, false),
            OutputOpts::default(),
        )
        .expect("distinct sources should concatenate");

        assert_eq!(
            headers_of(&out),
            ["@a1/1", "@a2/1", "@a3/1", "@a4/1", "@b1/1", "@b2/1"]
        );
    }

    #[test]
    fn a_duplicate_id_across_sources_is_rejected() {
        let a = source("dup_a", &["@x1/1", "@x2/1", "@x3/1"]);
        let b = source("dup_b", &["@y1/1", "@y2/1", "@x2/1"]);
        let out = scratch("out");

        let err = run(
            args(vec![a.clone(), b.clone()], vec![], out, None, false),
            OutputOpts::default(),
        )
        .expect_err("a repeated read ID must stop the run");
        let msg = err.to_string();

        assert!(msg.contains("duplicate read ID"), "unexpected error: {msg}");
        assert!(
            msg.contains("x2"),
            "the offending ID should be named: {msg}"
        );
        // Both ends of the collision have to be reportable, or the user has
        // no way to tell which input to go and look at.
        assert!(
            msg.contains(a.to_str().unwrap()) && msg.contains(b.to_str().unwrap()),
            "both source files should be named: {msg}"
        );
        // 1-based, and the two ends really are at different positions:
        // second record of the first file, third of the second.
        assert!(
            msg.contains("record #2") && msg.contains("record #3"),
            "positions should be reported 1-based and distinctly: {msg}"
        );
        assert!(
            msg.contains("--allow-duplicate-ids"),
            "the error should name the escape hatch: {msg}"
        );
    }

    #[test]
    fn the_same_file_listed_twice_is_rejected() {
        // The realistic operator error this check exists for: a typo'd path
        // or a glob that matched more than intended, silently doubling
        // coverage.
        let a = source("twice", &["@r1/1", "@r2/1"]);
        let out = scratch("out");

        let err = run(
            args(vec![a.clone(), a], vec![], out, None, false),
            OutputOpts::default(),
        )
        .expect_err("the same file twice must stop the run");
        assert!(
            err.to_string()
                .contains("did you accidentally include the same file twice"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn a_duplicate_within_one_source_is_rejected() {
        let a = source("self_dup", &["@r1/1", "@r2/1", "@r1/1"]);
        let out = scratch("out");

        let err = run(
            args(vec![a], vec![], out, None, false),
            OutputOpts::default(),
        )
        .expect_err("a file that repeats an ID internally must stop the run");
        assert!(
            err.to_string().contains("duplicate read ID"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn allow_duplicate_ids_skips_the_check_entirely() {
        let a = source("allowed", &["@r1/1", "@r2/1"]);
        let out = scratch("out");

        run(
            args(vec![a.clone(), a], vec![], out.clone(), None, true),
            OutputOpts::default(),
        )
        .expect("--allow-duplicate-ids should let it through");

        assert_eq!(headers_of(&out), ["@r1/1", "@r2/1", "@r1/1", "@r2/1"]);
    }

    #[test]
    fn mate_markers_do_not_make_two_reads_distinct() {
        // `/1` and `/2` are stripped by base_id(), so a single-end
        // concatenation of a file with its own mate file is still a duplicate
        // of every pair id -- which is what catches "I listed R1 and R2 as
        // two runs of the same mate".
        let r1 = source("mm_r1", &["@p1/1", "@p2/1"]);
        let r2 = source("mm_r2", &["@p1/2", "@p2/2"]);
        let out = scratch("out");

        let err = run(
            args(vec![r1, r2], vec![], out, None, false),
            OutputOpts::default(),
        )
        .expect_err("R1 and R2 of the same pairs share their IDs");
        assert!(
            err.to_string().contains("duplicate read ID"),
            "unexpected error: {err}"
        );
    }

    // ---- paired-end ----

    #[test]
    fn paired_sources_concatenate_both_mates_in_order() {
        let a1 = source("pe_a1", &["@a1/1", "@a2/1"]);
        let a2 = source("pe_a2", &["@a1/2", "@a2/2"]);
        let b1 = source("pe_b1", &["@b1/1", "@b2/1"]);
        let b2 = source("pe_b2", &["@b1/2", "@b2/2"]);
        let out1 = scratch("out1");
        let out2 = scratch("out2");

        run(
            args(
                vec![a1, b1],
                vec![a2, b2],
                out1.clone(),
                Some(out2.clone()),
                false,
            ),
            OutputOpts::default(),
        )
        .expect("paired sources should concatenate");

        assert_eq!(headers_of(&out1), ["@a1/1", "@a2/1", "@b1/1", "@b2/1"]);
        assert_eq!(headers_of(&out2), ["@a1/2", "@a2/2", "@b1/2", "@b2/2"]);
    }

    #[test]
    fn a_duplicate_pair_id_across_paired_sources_is_rejected() {
        let a1 = source("pedup_a1", &["@p1/1", "@p2/1"]);
        let a2 = source("pedup_a2", &["@p1/2", "@p2/2"]);
        let b1 = source("pedup_b1", &["@p2/1"]);
        let b2 = source("pedup_b2", &["@p2/2"]);
        let out1 = scratch("out1");
        let out2 = scratch("out2");

        let err = run(
            args(vec![a1, b1], vec![a2, b2], out1, Some(out2), false),
            OutputOpts::default(),
        )
        .expect_err("a repeated pair ID must stop the run");
        assert!(
            err.to_string().contains("duplicate read ID"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn desynced_mates_within_one_source_are_rejected() {
        // Checked per source pair, not just across runs: a single corrupt
        // run should be caught too.
        let r1 = source("desync_r1", &["@p1/1", "@p2/1", "@p3/1"]);
        let r2 = source("desync_r2", &["@p1/2", "@WRONG/2", "@p3/2"]);
        let out1 = scratch("out1");
        let out2 = scratch("out2");

        let err = run(
            args(vec![r1], vec![r2], out1, Some(out2), false),
            OutputOpts::default(),
        )
        .expect_err("mismatched mate IDs must stop the run");
        assert!(
            err.to_string().contains("read 1/2 desync"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn a_mate_count_mismatch_within_one_source_is_rejected() {
        let r1 = source("count_r1", &["@p1/1", "@p2/1", "@p3/1"]);
        let r2 = source("count_r2", &["@p1/2", "@p2/2"]);
        let out1 = scratch("out1");
        let out2 = scratch("out2");

        let err = run(
            args(vec![r1], vec![r2], out1, Some(out2), false),
            OutputOpts::default(),
        )
        .expect_err("unequal mate counts must stop the run");
        assert!(
            err.to_string().contains("different numbers of reads"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn mismatched_numbers_of_r1_and_r2_inputs_are_rejected() {
        let a1 = source("uneven_a1", &["@a/1"]);
        let b1 = source("uneven_b1", &["@b/1"]);
        let a2 = source("uneven_a2", &["@a/2"]);
        let out1 = scratch("out1");
        let out2 = scratch("out2");

        let err = run(
            args(vec![a1, b1], vec![a2], out1, Some(out2), false),
            OutputOpts::default(),
        )
        .expect_err("two --r1 and one --r2 is not a pairing");
        assert!(
            err.to_string()
                .contains("must be given the same number of times"),
            "unexpected error: {err}"
        );
    }
}
