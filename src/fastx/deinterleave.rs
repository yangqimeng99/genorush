//! `fastx deinterleave`: split a merged paired-end FASTQ back into R1/R2.
//!
//! "Merged" is ambiguous in the wild, in two independent ways.
//!
//! First, the *order*: proper interleaved files alternate `R1,R2,R1,R2,...`,
//! but plenty of files in circulation are just
//! `cat R1.fastq R2.fastq > merged.fastq` -- every R1 record, then every R2
//! record, back to back. These two layouts are not the same format wearing
//! different clothes; a splitter that assumes interleaved and gets a
//! concatenated file silently produces garbage (every "pair" is two
//! unrelated reads).
//!
//! Second, and more fundamentally: *whether position says anything at all*.
//! Both layouts above are positional -- which mate a record belongs to is a
//! function of where it sits. Real files exist where that is simply false.
//! An SRA-derived file observed in practice (SRR17458599, 399,065,686
//! records) is laid out as one 24,174,089-record run of R2, then all
//! 199,532,843 R1 records, then the remaining 175,358,754 R2 records --
//! three runs, with a globally sequential read numbering that gives the two
//! mates of a pair *different* ids (`@SRR17458599.1 1/2` and
//! `@SRR17458599.24174090 24174090/1`). No positional rule splits that file,
//! and no id comparison pairs it up. What does work is the marker each
//! record carries in its own header: `/1` or `/2`.
//!
//! So this command supports both families:
//!
//! - **Positional** (`--layout interleaved` / `--layout concat`): split by
//!   index. Fast, single pass for `interleaved`; `concat` needs the record
//!   count up front, so it pays a cheap count-only pre-pass.
//! - **By marker** (`--layout by-suffix`): route each record by the `/1` or
//!   `/2` marker in its header, ignoring position entirely. Single pass,
//!   O(1) memory, and immune to any ordering the file happens to use.
//!
//! `--layout auto` (the default) picks between them by probing the first few
//! thousand records: adjacent records sharing a read id means interleaved;
//! otherwise, every record carrying a mate marker means `by-suffix`. Only
//! when neither holds -- no markers *and* no shared ids -- does it fall back
//! to reading the whole input once to hash every `FastqRecord::base_id()`
//! (`common::hash::fnv1a`) and test the two positional hypotheses against
//! those hashes, which is the only case that needs 8 bytes of memory per
//! record. Earlier versions always took that path; on the 14.7 GB file above
//! that cost 5 minutes and 3.05 GiB of RSS only to conclude that neither
//! hypothesis held.
//!
//! Probing the head decides the layout but does not *prove* it, so the
//! chosen splitter re-checks its own assumption on every record as it goes
//! (interleaved: the two records of each pair must share a base id;
//! by-suffix: every record must carry a marker) and stops with a specific
//! record number if the assumption ever breaks. The trade is deliberate: a
//! violation is caught after some output has already been written rather
//! than before anything is written, in exchange for not reading multi-GB
//! inputs twice. Such a run exits non-zero and its outputs must be
//! discarded; `--layout` skips the probe entirely when the layout is
//! already known.
//!
//! Every read here goes through `common::fastq::spawn_reader` rather than
//! `io_utils::open_reader` directly, even though there's only ever one
//! input file (no second mate to decompress concurrently with). The reason
//! is the splitters, not detection: while they dispatch a chunk's worth of
//! records to the parallel-compressing `BlockWriter`, the reader thread
//! keeps decompressing/parsing the *next* chunk into its channel buffer in
//! the background instead of the main thread sitting idle waiting for
//! compression to finish -- the same overlap `fastx sample`/`fastx cat`
//! already get from reading two mates concurrently, available here even for
//! a single file.

use std::path::{Path, PathBuf};
use std::sync::mpsc::Receiver;
use std::time::Instant;

use anyhow::{bail, ensure, Result};
use clap::{Args, ValueEnum};

use crate::common::fastq::{format_into_blocks, spawn_reader, FastqRecord, Mate};
use crate::common::hash::fnv1a;
use crate::io_utils::{
    display, ensure_one_stdio_at_most, is_stdio, open_block_writer, BlockWriter, OutputOpts,
};

/// How many records the `auto` probe reads before deciding. Large enough
/// that an alternating-mate pattern can't hold by coincidence, small enough
/// that the probe is imperceptible next to the split itself.
const PROBE_RECORDS: usize = 4096;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum Layout {
    Auto,
    Interleaved,
    Concat,
    BySuffix,
}

#[derive(Args, Debug)]
pub struct DeinterleaveArgs {
    /// Merged input FASTQ (.fastq/.fq, gzip/bgzip auto-detected), or `-`
    /// to read stdin. A pipe can only be read once, so stdin rules out the
    /// layouts that need a second pass: `--layout concat` always, and
    /// `--layout auto` only if the head of the stream turns out to be
    /// inconclusive.
    #[arg(short = 'i', long = "in", value_name = "FILE")]
    input: PathBuf,

    /// Output for read 1, or `-` for stdout. Gzip-compressed if the path
    /// ends in `.gz` or `-z/--gzip` is passed.
    #[arg(short = 'o', long = "out1", value_name = "FILE")]
    out1: PathBuf,

    /// Output for read 2, or `-` for stdout -- though only one of the two
    /// outputs can be `-`, since both would land in the same pipe.
    /// Gzip-compressed if the path ends in `.gz` or `-z/--gzip` is passed.
    #[arg(short = 'O', long = "out2", value_name = "FILE")]
    out2: PathBuf,

    /// How R1/R2 were merged.
    ///
    /// `auto` (default) probes the first few thousand records and picks
    /// between `interleaved` and `by-suffix`, falling back to a full-input
    /// scan only when the head is inconclusive. `interleaved` and `concat`
    /// split by position (`R1,R2,R1,R2,...` and all-R1-then-all-R2
    /// respectively). `by-suffix` ignores position entirely and routes each
    /// record by the `/1` or `/2` marker in its own header -- the only mode
    /// that handles files whose mates are neither alternating, nor in two
    /// clean halves, nor sharing a read id.
    #[arg(long, value_enum, default_value_t = Layout::Auto)]
    layout: Layout,

    /// Skip the layout self-check the splitters run while writing (paired
    /// base ids for `interleaved`, per-half mate consistency for `concat`).
    /// Only disable this if your headers don't follow standard `/1`+`/2` or
    /// Illumina `1:...`+`2:...` mate conventions and the check produces
    /// false alarms. Has no effect on `by-suffix`, whose routing *is* the
    /// header check.
    #[arg(long)]
    no_pair_check: bool,

    /// Records processed per parallel compression batch.
    #[arg(long, default_value_t = 50_000)]
    chunk_records: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DetectedLayout {
    Interleaved,
    Concat,
}

/// A record stream that can have records handed back to it.
///
/// The probe has to read records to form an opinion, and those records are
/// still part of the input. On a file the splitter could simply re-open and
/// start over, which is what an earlier version did; on stdin there is no
/// starting over. Buffering the probed records here and replaying them
/// ahead of the channel makes one pass serve both, and removes the second
/// open in the file case as a side effect.
struct Records {
    preamble: std::vec::IntoIter<FastqRecord>,
    rx: Receiver<Result<FastqRecord>>,
}

impl Records {
    fn open(path: &Path) -> Result<Self> {
        Ok(Self {
            preamble: Vec::new().into_iter(),
            rx: spawn_reader(path.to_path_buf())?,
        })
    }

    fn next_record(&mut self) -> Option<Result<FastqRecord>> {
        if let Some(rec) = self.preamble.next() {
            return Some(Ok(rec));
        }
        self.rx.recv().ok()
    }

    /// Puts already-read records back at the front of the stream. Only
    /// valid while nothing is pending, which is the case right after a
    /// probe has drained everything it buffered.
    fn replay(&mut self, records: Vec<FastqRecord>) {
        debug_assert_eq!(self.preamble.len(), 0);
        self.preamble = records.into_iter();
    }

    /// Refills `chunk` with up to `n` records, leaving it empty at EOF.
    fn fill(&mut self, chunk: &mut Vec<FastqRecord>, n: usize) -> Result<()> {
        chunk.clear();
        while chunk.len() < n {
            match self.next_record() {
                Some(r) => chunk.push(r?),
                None => break,
            }
        }
        Ok(())
    }
}

/// What the first `PROBE_RECORDS` records say about the file's layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct HeadProbe {
    records: usize,
    /// Every probed record carried a `/1` or `/2` marker.
    all_have_marker: bool,
    /// Every complete adjacent pair `(2i, 2i+1)` shared a base id, i.e. the
    /// head looks interleaved. False when fewer than two records were read.
    pairs_share_id: bool,
}

/// Reads at most `k` records off the front of `src`, reports what they imply
/// about the layout, and replays them so the split still sees them. Cheap:
/// it decompresses only as far as it needs to.
fn probe_head(src: &mut Records, k: usize) -> Result<HeadProbe> {
    let mut buf: Vec<FastqRecord> = Vec::with_capacity(k);
    while buf.len() < k {
        match src.next_record() {
            Some(r) => buf.push(r?),
            None => break,
        }
    }

    let records = buf.len();
    let all_have_marker = buf.iter().all(|r| r.mate_suffix().is_some());
    let pairs_share_id = buf
        .as_chunks::<2>()
        .0
        .iter()
        .all(|pair| pair[0].base_id() == pair[1].base_id());

    src.replay(buf);
    Ok(HeadProbe {
        records,
        all_have_marker: all_have_marker && records > 0,
        pairs_share_id: pairs_share_id && records >= 2,
    })
}

/// Reads the whole input once, hashing every record's base ID, and tests
/// both positional layout hypotheses against those hashes. Returns the layout
/// and the total record count (needed by the concat split, computed here for
/// free). Only reached when the head probe was inconclusive, since this is
/// the one code path whose memory grows with the input (8 bytes/record).
fn detect_layout(path: &Path) -> Result<(DetectedLayout, usize)> {
    let rx = spawn_reader(path.to_path_buf())?;
    let mut hashes: Vec<u64> = Vec::new();
    for r in rx.iter() {
        hashes.push(fnv1a(r?.base_id().as_bytes()));
    }

    let n = hashes.len();
    ensure!(n > 0, "input is empty, nothing to deinterleave");
    ensure!(
        n.is_multiple_of(2),
        "input has an odd number of records ({n}); a merged paired-end file must have an even count"
    );
    let half = n / 2;

    let interleaved_ok = (0..half).all(|i| hashes[2 * i] == hashes[2 * i + 1]);
    let concat_ok = (0..half).all(|i| hashes[i] == hashes[i + half]);

    match (interleaved_ok, concat_ok) {
        (true, false) => Ok((DetectedLayout::Interleaved, n)),
        (false, true) => Ok((DetectedLayout::Concat, n)),
        (true, true) => bail!(
            "input matches both the interleaved and R1-then-R2-concatenation layouts \
             (likely a tiny or degenerate file); pass --layout explicitly to disambiguate"
        ),
        (false, false) => {
            let first_interleaved_mismatch =
                (0..half).find(|&i| hashes[2 * i] != hashes[2 * i + 1]);
            let first_concat_mismatch = (0..half).find(|&i| hashes[i] != hashes[i + half]);
            bail!(
                "could not determine how R1/R2 were merged: not a clean interleaved file \
                 (first mismatched pair at index {first_interleaved_mismatch:?}), not a clean \
                 R1-then-R2 concatenation (first mismatch at index {first_concat_mismatch:?}), \
                 and the records carry no /1 or /2 mate markers to route by instead; \
                 the file may be corrupt or use a merge convention this tool doesn't recognize -- \
                 pass --layout explicitly if you already know which one it is"
            );
        }
    }
}

/// Counts records without hashing -- used for `--layout concat` given
/// explicitly, which still needs the midpoint but not the (unnecessary)
/// hypothesis test detection would otherwise do.
fn count_records(path: &Path) -> Result<usize> {
    let rx = spawn_reader(path.to_path_buf())?;
    let mut n = 0usize;
    for r in rx.iter() {
        r?;
        n += 1;
    }
    Ok(n)
}

/// Single streaming pass: record `i` (0-indexed, global) goes to R1 if
/// even, R2 if odd. No pre-pass needed -- unlike concat, interleaved
/// doesn't need to know the total count up front.
///
/// `chunk_records` is rounded up to an even number so a pair never straddles
/// two chunks, which keeps the `check_ids` verification a within-chunk
/// comparison of two borrowed records instead of having to carry an owned
/// id across the boundary.
fn split_interleaved(
    src: &mut Records,
    w1: &mut BlockWriter,
    w2: &mut BlockWriter,
    chunk_records: usize,
    check_ids: bool,
) -> Result<(u64, u64)> {
    let chunk_records = chunk_records + (chunk_records % 2);
    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(chunk_records);
    let mut total: u64 = 0;
    loop {
        src.fill(&mut chunk, chunk_records)?;
        if chunk.is_empty() {
            break;
        }
        if check_ids {
            for (i, pair) in chunk.as_chunks::<2>().0.iter().enumerate() {
                if pair[0].base_id() != pair[1].base_id() {
                    bail!(
                        "not an interleaved file: records #{} and #{} do not share a read id \
                         ({:?} vs {:?}). If the mates are laid out some other way, try \
                         --layout by-suffix (routes by the /1 and /2 header markers) or \
                         --layout concat; pass --no-pair-check to split by position anyway. \
                         The output files are incomplete and must be discarded.",
                        total + (2 * i) as u64 + 1,
                        total + (2 * i) as u64 + 2,
                        pair[0].base_id(),
                        pair[1].base_id()
                    );
                }
            }
        }
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        for (i, rec) in chunk.iter().enumerate() {
            if (total + i as u64).is_multiple_of(2) {
                r1.push(rec);
            } else {
                r2.push(rec);
            }
        }
        w1.write_blocks(format_into_blocks(&r1)?)?;
        w2.write_blocks(format_into_blocks(&r2)?)?;
        total += chunk.len() as u64;
    }
    ensure!(
        total.is_multiple_of(2),
        "input has an odd number of records ({total}); not a valid interleaved paired file"
    );
    Ok((total / 2, total / 2))
}

/// Requires the midpoint (`n`) up front: record `i` goes to R1 if `i < n/2`,
/// R2 otherwise. A chunk straddling the midpoint is split within itself, so
/// chunk boundaries don't need to align with it.
///
/// With `check_mates` on, every record that carries a `/1` or `/2` marker
/// must agree with the other records in its half, and the two halves must
/// disagree with each other. That is exactly the check that catches a file
/// whose mate runs don't start and end at the midpoint -- the failure mode
/// that otherwise produces a full-looking, silently wrong pair of outputs.
fn split_concat(
    src: &mut Records,
    w1: &mut BlockWriter,
    w2: &mut BlockWriter,
    chunk_records: usize,
    n: usize,
    check_mates: bool,
) -> Result<(u64, u64)> {
    ensure!(
        n.is_multiple_of(2),
        "input has an odd number of records ({n}); a merged paired-end file must have an even count"
    );
    let half = n / 2;
    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(chunk_records);
    let mut total: usize = 0;
    // The mate each half has settled on, learned from the first marked
    // record seen in that half.
    let mut half_mate: [Option<Mate>; 2] = [None, None];
    loop {
        src.fill(&mut chunk, chunk_records)?;
        if chunk.is_empty() {
            break;
        }
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        for (i, rec) in chunk.iter().enumerate() {
            let idx = total + i;
            let side = usize::from(idx >= half);
            if check_mates {
                if let Some(mate) = rec.mate_suffix() {
                    match half_mate[side] {
                        None => half_mate[side] = Some(mate),
                        Some(settled) if settled != mate => bail!(
                            "not an R1-then-R2 concatenation: record #{} is {:?} but the \
                             {} half started out as {:?}, so the mate runs do not line up with \
                             the midpoint at record #{}. Try --layout by-suffix, which routes by \
                             the /1 and /2 header markers and does not care where the runs \
                             start; pass --no-pair-check to split by position anyway. The \
                             output files are incomplete and must be discarded.",
                            idx + 1,
                            mate,
                            if side == 0 { "first" } else { "second" },
                            settled,
                            half + 1
                        ),
                        Some(_) => {}
                    }
                }
            }
            if side == 0 {
                r1.push(rec);
            } else {
                r2.push(rec);
            }
        }
        w1.write_blocks(format_into_blocks(&r1)?)?;
        w2.write_blocks(format_into_blocks(&r2)?)?;
        total += chunk.len();
    }
    ensure!(
        total == n,
        "record count changed between passes ({n} then {total}); input may have changed while running"
    );
    if let [Some(first), Some(second)] = half_mate {
        // Within-half consistency alone would accept two R1 files stapled
        // together, which is not a merged paired-end file at all.
        ensure!(
            first != second,
            "not an R1-then-R2 concatenation: both halves are {first:?}, so this file does not \
             contain two different mates. The output files are incomplete and must be discarded."
        );
        if first == Mate::R2 {
            log::warn!(
                "the first half of this file is R2 and the second half is R1, so --out1 now \
                 holds the R2 reads and --out2 the R1 reads. --layout by-suffix would label \
                 them by their own headers instead."
            );
        }
    }
    Ok((half as u64, (n - half) as u64))
}

/// Routes every record by the `/1` or `/2` marker in its own header, with no
/// reference to position at all. Single pass, O(1) memory beyond one chunk.
///
/// A record with no marker is fatal: this mode exists precisely because
/// position carries no information in such files, so there is nothing to
/// fall back on for that record, and silently dropping it (as the shell
/// `awk` one-liners that inspired this mode do) loses data.
fn split_by_suffix(
    src: &mut Records,
    w1: &mut BlockWriter,
    w2: &mut BlockWriter,
    chunk_records: usize,
) -> Result<(u64, u64)> {
    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(chunk_records);
    let mut total: u64 = 0;
    let mut n1: u64 = 0;
    let mut n2: u64 = 0;
    loop {
        src.fill(&mut chunk, chunk_records)?;
        if chunk.is_empty() {
            break;
        }
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        for (i, rec) in chunk.iter().enumerate() {
            match rec.mate_suffix() {
                Some(Mate::R1) => r1.push(rec),
                Some(Mate::R2) => r2.push(rec),
                None => bail!(
                    "record #{} carries no /1 or /2 mate marker in its header: {:?}. \
                     --layout by-suffix has nothing else to route it by. If this file is \
                     positional after all, use --layout interleaved or --layout concat. \
                     The output files are incomplete and must be discarded.",
                    total + i as u64 + 1,
                    rec.header
                ),
            }
        }
        n1 += r1.len() as u64;
        n2 += r2.len() as u64;
        w1.write_blocks(format_into_blocks(&r1)?)?;
        w2.write_blocks(format_into_blocks(&r2)?)?;
        total += chunk.len() as u64;
    }
    Ok((n1, n2))
}

/// A layout that needs the record count up front, or a second look at the
/// input, cannot be served by a pipe: there is no second look. Reject that
/// combination with a message that says what to do instead, rather than
/// half-reading stdin and failing somewhere less obvious.
fn ensure_rereadable(input: &Path, what: &str) -> Result<()> {
    ensure!(
        !is_stdio(input),
        "{what}, which is impossible when the input is stdin -- a pipe can only be read once. \
         Either give the input as a file, or pass a single-pass --layout \
         (interleaved or by-suffix) if you already know how the file is laid out."
    );
    Ok(())
}

pub fn run(args: DeinterleaveArgs, opts: OutputOpts) -> Result<()> {
    ensure!(args.chunk_records > 0, "--chunk-records must be > 0");
    ensure_one_stdio_at_most(&[&args.out1, &args.out2], "output")?;

    let start = Instant::now();
    let mut w1 = open_block_writer(&args.out1, opts)?;
    let mut w2 = open_block_writer(&args.out2, opts)?;
    let check = !args.no_pair_check;
    log::info!(
        "splitting {} -> {} + {}",
        display(&args.input),
        display(&args.out1),
        display(&args.out2)
    );

    let (n1, n2) = match args.layout {
        Layout::Interleaved => {
            log::info!("layout: interleaved (explicit)");
            let mut src = Records::open(&args.input)?;
            split_interleaved(&mut src, &mut w1, &mut w2, args.chunk_records, check)?
        }
        Layout::Concat => {
            ensure_rereadable(
                &args.input,
                "--layout concat needs the total record count before it can find the midpoint, \
                 so it reads the input twice",
            )?;
            log::info!("layout: concat (explicit); counting records first");
            let n = count_records(&args.input)?;
            let mut src = Records::open(&args.input)?;
            split_concat(&mut src, &mut w1, &mut w2, args.chunk_records, n, check)?
        }
        Layout::BySuffix => {
            log::info!("layout: by-suffix (explicit)");
            let mut src = Records::open(&args.input)?;
            split_by_suffix(&mut src, &mut w1, &mut w2, args.chunk_records)?
        }
        Layout::Auto => {
            let mut src = Records::open(&args.input)?;
            let probe = probe_head(&mut src, PROBE_RECORDS)?;
            ensure!(probe.records > 0, "input is empty, nothing to deinterleave");
            if probe.pairs_share_id {
                log::info!(
                    "probed {} records: adjacent records share a read id -- splitting as \
                     interleaved, verifying each pair while writing",
                    probe.records
                );
                split_interleaved(&mut src, &mut w1, &mut w2, args.chunk_records, check)?
            } else if probe.all_have_marker {
                log::info!(
                    "probed {} records: no shared ids, but every record carries a /1 or /2 \
                     marker -- routing by marker (--layout by-suffix)",
                    probe.records
                );
                split_by_suffix(&mut src, &mut w1, &mut w2, args.chunk_records)?
            } else {
                drop(src);
                ensure_rereadable(
                    &args.input,
                    "the first records carry no mate markers and no shared ids, so the layout \
                     can only be settled by scanning the whole input and then splitting it",
                )?;
                log::info!(
                    "probed {} records: inconclusive (no mate markers, no shared ids) -- \
                     falling back to scanning the whole input to detect the layout",
                    probe.records
                );
                let (layout, n) = detect_layout(&args.input)?;
                log::info!("detected layout: {layout:?} ({n} total records)");
                let mut src = Records::open(&args.input)?;
                match layout {
                    DetectedLayout::Interleaved => {
                        split_interleaved(&mut src, &mut w1, &mut w2, args.chunk_records, check)?
                    }
                    DetectedLayout::Concat => {
                        split_concat(&mut src, &mut w1, &mut w2, args.chunk_records, n, check)?
                    }
                }
            }
        }
    };
    w1.flush()?;
    w2.flush()?;

    log::info!(
        "split {} records into {n1}/{n2} (R1/R2) in {:.2?}",
        n1 + n2,
        start.elapsed()
    );
    ensure!(n1 + n2 > 0, "input is empty, nothing to deinterleave");
    ensure!(
        n1 == n2,
        "R1 and R2 ended up with different numbers of records ({n1} vs {n2}); both output \
         files were written in full, but they cannot be treated as positionally paired \
         downstream -- check whether the input is really a merged paired-end file"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::fastq::read_fastq_record;
    use std::fs;
    use std::io::Cursor;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// Test scratch space lives under the build directory rather than the
    /// system temp dir: it is already git-ignored, it is guaranteed writable
    /// wherever `cargo` itself can write, and it keeps test artifacts next to
    /// the build they came from instead of scattered in /tmp.
    fn scratch(name: &str) -> PathBuf {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("target")
            .join("test-scratch");
        fs::create_dir_all(&dir).expect("failed to create test scratch dir");
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        dir.join(format!("{name}.{n}.{}.fastq", std::process::id()))
    }

    fn write_fastq(path: &Path, headers: &[&str]) {
        let mut body = String::new();
        for h in headers {
            body.push_str(&format!("{h}\nACGT\n+\nIIII\n"));
        }
        fs::write(path, body).expect("failed to write test fastq");
    }

    fn headers_of(path: &Path) -> Vec<String> {
        let body = fs::read_to_string(path).expect("failed to read output");
        let mut reader = Cursor::new(body.into_bytes());
        let mut out = Vec::new();
        let mut line_no = 1;
        while let Some(rec) = read_fastq_record(&mut reader, line_no).expect("malformed output") {
            out.push(rec.header);
            line_no += 4;
        }
        out
    }

    /// Runs a split over a temporary input and returns the two outputs'
    /// headers, so tests can assert on routing without caring about bodies.
    fn split(headers: &[&str], layout: Layout, check: bool) -> Result<(Vec<String>, Vec<String>)> {
        let input = scratch("in");
        let out1 = scratch("out1");
        let out2 = scratch("out2");
        write_fastq(&input, headers);
        let mut w1 = open_block_writer(&out1, OutputOpts::default())?;
        let mut w2 = open_block_writer(&out2, OutputOpts::default())?;
        let result = (|| -> Result<()> {
            let mut src = Records::open(&input)?;
            match layout {
                Layout::Interleaved => {
                    split_interleaved(&mut src, &mut w1, &mut w2, 4, check)?;
                }
                Layout::Concat => {
                    let n = count_records(&input)?;
                    split_concat(&mut src, &mut w1, &mut w2, 4, n, check)?;
                }
                Layout::BySuffix => {
                    split_by_suffix(&mut src, &mut w1, &mut w2, 4)?;
                }
                Layout::Auto => unreachable!("tests dispatch layouts explicitly"),
            }
            Ok(())
        })();
        w1.flush()?;
        w2.flush()?;
        result?;
        Ok((headers_of(&out1), headers_of(&out2)))
    }

    // ---- existing positional behavior, previously untested ----

    #[test]
    fn interleaved_splits_alternating_records() {
        let (r1, r2) = split(&["@a/1", "@a/2", "@b/1", "@b/2"], Layout::Interleaved, true).unwrap();
        assert_eq!(r1, ["@a/1", "@b/1"]);
        assert_eq!(r2, ["@a/2", "@b/2"]);
    }

    #[test]
    fn interleaved_rejects_pairs_that_do_not_share_an_id() {
        let err = split(&["@a/1", "@b/2", "@c/1", "@d/2"], Layout::Interleaved, true)
            .expect_err("mismatched pair must be rejected");
        assert!(
            err.to_string().contains("not an interleaved file"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn interleaved_no_pair_check_splits_by_position_regardless() {
        let (r1, r2) = split(
            &["@a/1", "@b/2", "@c/1", "@d/2"],
            Layout::Interleaved,
            false,
        )
        .unwrap();
        assert_eq!(r1, ["@a/1", "@c/1"]);
        assert_eq!(r2, ["@b/2", "@d/2"]);
    }

    #[test]
    fn concat_splits_at_the_midpoint() {
        let (r1, r2) = split(&["@a/1", "@b/1", "@a/2", "@b/2"], Layout::Concat, true).unwrap();
        assert_eq!(r1, ["@a/1", "@b/1"]);
        assert_eq!(r2, ["@a/2", "@b/2"]);
    }

    #[test]
    fn concat_splits_across_a_chunk_boundary() {
        // 10 records with chunk size 4: the midpoint (5) falls inside the
        // second chunk, so the chunk has to be split within itself.
        let headers: Vec<String> = (0..5)
            .map(|i| format!("@r{i}/1"))
            .chain((0..5).map(|i| format!("@r{i}/2")))
            .collect();
        let refs: Vec<&str> = headers.iter().map(|s| s.as_str()).collect();
        let (r1, r2) = split(&refs, Layout::Concat, true).unwrap();
        assert_eq!(r1, ["@r0/1", "@r1/1", "@r2/1", "@r3/1", "@r4/1"]);
        assert_eq!(r2, ["@r0/2", "@r1/2", "@r2/2", "@r3/2", "@r4/2"]);
    }

    #[test]
    fn concat_rejects_mate_runs_that_miss_the_midpoint() {
        // One R2, then two R1, then one R2: mate runs exist but don't line up
        // with the midpoint. This is the shape that used to produce a
        // silently wrong split.
        let err = split(&["@a/2", "@b/1", "@c/1", "@d/2"], Layout::Concat, true)
            .expect_err("misaligned mate runs must be rejected");
        assert!(
            err.to_string().contains("not an R1-then-R2 concatenation"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn concat_rejects_two_halves_of_the_same_mate() {
        // Two R1 files stapled together: each half is internally consistent,
        // so only the across-halves check catches it.
        let err = split(&["@a/1", "@b/1", "@c/1", "@d/1"], Layout::Concat, true)
            .expect_err("two same-mate halves must be rejected");
        assert!(
            err.to_string().contains("both halves are"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn concat_without_markers_still_splits_by_position() {
        let (r1, r2) = split(&["@a", "@b", "@a", "@b"], Layout::Concat, true).unwrap();
        assert_eq!(r1, ["@a", "@b"]);
        assert_eq!(r2, ["@a", "@b"]);
    }

    // ---- by-suffix routing ----

    #[test]
    fn by_suffix_routes_regardless_of_order() {
        // The SRR17458599 shape in miniature: an R2 run, then an R1 run, then
        // another R2 run, with every record numbered differently.
        let (r1, r2) = split(
            &["@SRR.1 1/2", "@SRR.2 2/1", "@SRR.3 3/1", "@SRR.4 4/2"],
            Layout::BySuffix,
            true,
        )
        .unwrap();
        assert_eq!(r1, ["@SRR.2 2/1", "@SRR.3 3/1"]);
        assert_eq!(r2, ["@SRR.1 1/2", "@SRR.4 4/2"]);
    }

    #[test]
    fn by_suffix_preserves_input_order_within_each_mate() {
        let headers: Vec<String> = (0..12)
            .map(|i| format!("@r{i} {i}/{}", if i % 3 == 0 { 1 } else { 2 }))
            .collect();
        let refs: Vec<&str> = headers.iter().map(|s| s.as_str()).collect();
        let (r1, r2) = split(&refs, Layout::BySuffix, true).unwrap();
        assert_eq!(r1, ["@r0 0/1", "@r3 3/1", "@r6 6/1", "@r9 9/1"]);
        assert_eq!(r2.len(), 8);
        assert_eq!(r2[0], "@r1 1/2");
        assert_eq!(r2[7], "@r11 11/2");
    }

    #[test]
    fn by_suffix_rejects_a_record_without_a_marker() {
        let err = split(&["@a/1", "@b 1:N:0:AT", "@c/2"], Layout::BySuffix, true)
            .expect_err("a record with no marker must be rejected");
        assert!(
            err.to_string().contains("carries no /1 or /2 mate marker"),
            "unexpected error: {err}"
        );
    }

    // ---- head probe ----

    fn probe(headers: &[&str]) -> HeadProbe {
        let input = scratch("probe");
        write_fastq(&input, headers);
        let mut src = Records::open(&input).unwrap();
        probe_head(&mut src, PROBE_RECORDS).unwrap()
    }

    /// A probe must hand back everything it consumed, or the split that
    /// follows it silently loses the head of the input -- the whole reason
    /// the stream is replayable.
    #[test]
    fn probe_replays_every_record_it_consumed() {
        let input = scratch("replay");
        let headers = ["@a/1", "@a/2", "@b/1", "@b/2"];
        write_fastq(&input, &headers);
        let mut src = Records::open(&input).unwrap();
        probe_head(&mut src, PROBE_RECORDS).unwrap();

        let mut seen = Vec::new();
        while let Some(r) = src.next_record() {
            seen.push(r.unwrap().header);
        }
        assert_eq!(seen, headers);
    }

    #[test]
    fn probe_recognizes_an_interleaved_head() {
        let p = probe(&["@a/1", "@a/2", "@b/1", "@b/2"]);
        assert!(p.pairs_share_id);
        assert!(p.all_have_marker);
        assert_eq!(p.records, 4);
    }

    #[test]
    fn probe_recognizes_a_marker_only_head() {
        // SRA shape: no shared ids anywhere, but every record self-identifies.
        let p = probe(&["@SRR.1 1/2", "@SRR.2 2/2", "@SRR.3 3/2", "@SRR.4 4/2"]);
        assert!(!p.pairs_share_id);
        assert!(p.all_have_marker);
    }

    #[test]
    fn probe_reports_inconclusive_head() {
        let p = probe(&["@a", "@b", "@c", "@d"]);
        assert!(!p.pairs_share_id);
        assert!(!p.all_have_marker);
    }

    #[test]
    fn probe_of_a_single_record_cannot_claim_interleaving() {
        let p = probe(&["@a/1"]);
        assert!(!p.pairs_share_id);
        assert!(p.all_have_marker);
        assert_eq!(p.records, 1);
    }
}
