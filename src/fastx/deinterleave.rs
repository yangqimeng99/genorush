//! `fastx deinterleave`: split a merged paired-end FASTQ back into R1/R2.
//!
//! "Merged" is ambiguous in the wild: proper interleaved files alternate
//! `R1,R2,R1,R2,...`, but plenty of files in circulation are just
//! `cat R1.fastq R2.fastq > merged.fastq` -- every R1 record, then every R2
//! record, back to back. These two layouts are not the same format wearing
//! different clothes; a splitter that assumes interleaved and gets a
//! concatenated file silently produces garbage (every "pair" is two
//! unrelated reads). This command detects which layout it's actually
//! looking at instead of assuming.
//!
//! Detection needs the *total* record count before it can even test the
//! concatenation hypothesis (record `i` must match record `i + n/2`), so
//! `--layout auto` (the default) reads the whole input once to hash every
//! `FastqRecord::base_id()` (`common::hash::fnv1a`, 8 bytes/record, not the
//! full string) and test both hypotheses, then reads it a second time to
//! actually perform the split. If neither hypothesis holds cleanly, this
//! refuses to guess and reports exactly where each one first broke down.
//! `--layout interleaved`/`--layout concat` skip detection for a known
//! layout -- `interleaved` needs only a single streaming pass with no
//! extra memory; `concat` still needs the record count up front (a
//! count-only pass, cheaper than full detection since it skips hashing).
//!
//! All of the above is *positional*: which mate a record belongs to is a
//! function of where it sits. Real files exist where that premise is simply
//! false. An SRA-derived file observed in practice (SRR17458599,
//! 399,065,686 records) is laid out as one 24,174,089-record run of R2,
//! then all 199,532,843 R1 records, then the remaining 175,358,754 R2
//! records -- three runs, with a globally sequential read numbering that
//! gives the two mates of a pair *different* ids (`@SRR17458599.1 1/2` and
//! `@SRR17458599.24174090 24174090/1`). No positional rule splits that
//! file, and no id comparison pairs it up. What does work is the marker
//! each record carries in its own header, which is what `--layout
//! by-suffix` routes on: a single pass, O(1) memory, immune to whatever
//! ordering the file happens to use.
//!
//! Every read here goes through `common::fastq::spawn_reader` rather than
//! `io_utils::open_reader` directly, even though there's only ever one
//! input file (no second mate to decompress concurrently with). The reason
//! is `split_interleaved`/`split_concat`, not detection: while those
//! functions dispatch a chunk's worth of records to the parallel-compressing
//! `BlockWriter`, the reader thread keeps decompressing/parsing the *next*
//! chunk into its channel buffer in the background instead of the main
//! thread sitting idle waiting for compression to finish -- the same
//! overlap `fastx sample`/`fastx cat` already get from reading two mates
//! concurrently, available here even for a single file. (`detect_layout`/
//! `count_records` also use it, for consistency and because it lets go of
//! manual `line_no` bookkeeping, but since those two functions do almost no
//! work besides reading, the wall-clock benefit there is close to zero --
//! decompression is already the sole bottleneck either way.)

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, ensure, Result};
use clap::{Args, ValueEnum};

use crate::common::fastq::{format_into_blocks, spawn_reader, FastqRecord, Mate};
use crate::common::hash::fnv1a;
use crate::io_utils::{open_block_writer, BlockWriter};

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum Layout {
    Auto,
    Interleaved,
    Concat,
    BySuffix,
}

#[derive(Args, Debug)]
pub struct DeinterleaveArgs {
    /// Merged input FASTQ (.fastq/.fq, gzip/bgzip auto-detected).
    #[arg(short = 'i', long = "in", value_name = "FILE")]
    input: PathBuf,

    /// Output for read 1. Gzip-compressed if the path ends in `.gz`.
    #[arg(short = 'o', long = "out1", value_name = "FILE")]
    out1: PathBuf,

    /// Output for read 2. Gzip-compressed if the path ends in `.gz`.
    #[arg(short = 'O', long = "out2", value_name = "FILE")]
    out2: PathBuf,

    /// How R1/R2 were merged.
    ///
    /// `auto` (default) detects it by reading the file once.
    /// `interleaved` and `concat` skip detection and split by position
    /// (`R1,R2,R1,R2,...` and all-R1-then-all-R2 respectively).
    /// `by-suffix` ignores position entirely and routes each record by the
    /// `/1` or `/2` marker in its own header -- the only mode that handles
    /// files whose mates are neither alternating, nor in two clean halves,
    /// nor sharing a read id.
    #[arg(long, value_enum, default_value_t = Layout::Auto)]
    layout: Layout,

    /// Records processed per parallel compression batch.
    #[arg(long, default_value_t = 50_000)]
    chunk_records: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DetectedLayout {
    Interleaved,
    Concat,
}

/// Reads the whole input once, hashing every record's base ID, and tests
/// both layout hypotheses against those hashes. Returns the layout and the
/// total record count (needed by the concat split, computed here for free).
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
                 (first mismatched pair at index {first_interleaved_mismatch:?}) and not a clean \
                 R1-then-R2 concatenation (first mismatch at index {first_concat_mismatch:?}); \
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
fn split_interleaved(
    input: &Path,
    w1: &mut BlockWriter,
    w2: &mut BlockWriter,
    chunk_records: usize,
) -> Result<u64> {
    let rx = spawn_reader(input.to_path_buf())?;
    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(chunk_records);
    let mut total: u64 = 0;
    loop {
        chunk.clear();
        for r in rx.iter().take(chunk_records) {
            chunk.push(r?);
        }
        if chunk.is_empty() {
            break;
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
    Ok(total)
}

/// Requires the midpoint (`n`) up front: record `i` goes to R1 if `i < n/2`,
/// R2 otherwise. A chunk straddling the midpoint is split within itself, so
/// chunk boundaries don't need to align with it.
fn split_concat(
    input: &Path,
    w1: &mut BlockWriter,
    w2: &mut BlockWriter,
    chunk_records: usize,
    n: usize,
) -> Result<u64> {
    ensure!(
        n.is_multiple_of(2),
        "input has an odd number of records ({n}); a merged paired-end file must have an even count"
    );
    let half = n / 2;
    let rx = spawn_reader(input.to_path_buf())?;
    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(chunk_records);
    let mut total: usize = 0;
    loop {
        chunk.clear();
        for r in rx.iter().take(chunk_records) {
            chunk.push(r?);
        }
        if chunk.is_empty() {
            break;
        }
        let mut r1 = Vec::new();
        let mut r2 = Vec::new();
        for (i, rec) in chunk.iter().enumerate() {
            if total + i < half {
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
    Ok(total as u64)
}

/// Routes every record by the `/1` or `/2` marker in its own header, with no
/// reference to position at all. Single pass, O(1) memory beyond one chunk.
/// Returns the per-mate record counts, which -- unlike a positional split --
/// are not equal by construction.
///
/// A record with no marker is fatal: this mode exists precisely because
/// position carries no information in such files, so there is nothing to
/// fall back on for that record, and silently dropping it (as the shell
/// `awk` one-liners that inspired this mode do) loses data.
fn split_by_suffix(
    input: &Path,
    w1: &mut BlockWriter,
    w2: &mut BlockWriter,
    chunk_records: usize,
) -> Result<(u64, u64)> {
    let rx = spawn_reader(input.to_path_buf())?;
    let mut chunk: Vec<FastqRecord> = Vec::with_capacity(chunk_records);
    let mut total: u64 = 0;
    let mut n1: u64 = 0;
    let mut n2: u64 = 0;
    loop {
        chunk.clear();
        for r in rx.iter().take(chunk_records) {
            chunk.push(r?);
        }
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

pub fn run(args: DeinterleaveArgs) -> Result<()> {
    let start = Instant::now();
    let mut w1 = open_block_writer(&args.out1)?;
    let mut w2 = open_block_writer(&args.out2)?;

    let (n1, n2) = match args.layout {
        Layout::Interleaved => {
            log::info!("layout: interleaved (explicit)");
            let total = split_interleaved(&args.input, &mut w1, &mut w2, args.chunk_records)?;
            (total / 2, total / 2)
        }
        Layout::Concat => {
            log::info!("layout: concat (explicit); counting records first");
            let n = count_records(&args.input)?;
            let total = split_concat(&args.input, &mut w1, &mut w2, args.chunk_records, n)?;
            (total / 2, total / 2)
        }
        Layout::BySuffix => {
            log::info!("layout: by-suffix (explicit)");
            split_by_suffix(&args.input, &mut w1, &mut w2, args.chunk_records)?
        }
        Layout::Auto => {
            log::info!("detecting merge layout (reads the input once)");
            let (layout, n) = detect_layout(&args.input)?;
            log::info!("detected layout: {layout:?} ({n} total records)");
            let total = match layout {
                DetectedLayout::Interleaved => {
                    split_interleaved(&args.input, &mut w1, &mut w2, args.chunk_records)?
                }
                DetectedLayout::Concat => {
                    split_concat(&args.input, &mut w1, &mut w2, args.chunk_records, n)?
                }
            };
            (total / 2, total / 2)
        }
    };
    w1.flush()?;
    w2.flush()?;

    log::info!(
        "split {} records into {n1}/{n2} (R1/R2) in {:.2?}",
        n1 + n2,
        start.elapsed()
    );
    // Positional splits get this for free; by-suffix cannot, and unequal
    // counts mean the outputs -- complete as they are -- can't be treated as
    // positionally paired downstream.
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

    /// Runs a by-suffix split over a temporary input and returns the two
    /// outputs' headers, so tests can assert on routing without caring about
    /// record bodies.
    fn split_suffix(headers: &[&str]) -> Result<(Vec<String>, Vec<String>)> {
        let input = scratch("in");
        let out1 = scratch("out1");
        let out2 = scratch("out2");
        write_fastq(&input, headers);
        let mut w1 = open_block_writer(&out1)?;
        let mut w2 = open_block_writer(&out2)?;
        let result = split_by_suffix(&input, &mut w1, &mut w2, 4);
        w1.flush()?;
        w2.flush()?;
        result?;
        Ok((headers_of(&out1), headers_of(&out2)))
    }

    #[test]
    fn by_suffix_routes_regardless_of_order() {
        // The SRR17458599 shape in miniature: an R2 run, then an R1 run, then
        // another R2 run, with every record numbered differently.
        let (r1, r2) =
            split_suffix(&["@SRR.1 1/2", "@SRR.2 2/1", "@SRR.3 3/1", "@SRR.4 4/2"]).unwrap();
        assert_eq!(r1, ["@SRR.2 2/1", "@SRR.3 3/1"]);
        assert_eq!(r2, ["@SRR.1 1/2", "@SRR.4 4/2"]);
    }

    #[test]
    fn by_suffix_preserves_input_order_within_each_mate() {
        let headers: Vec<String> = (0..12)
            .map(|i| format!("@r{i} {i}/{}", if i % 3 == 0 { 1 } else { 2 }))
            .collect();
        let refs: Vec<&str> = headers.iter().map(|s| s.as_str()).collect();
        let (r1, r2) = split_suffix(&refs).unwrap();
        assert_eq!(r1, ["@r0 0/1", "@r3 3/1", "@r6 6/1", "@r9 9/1"]);
        assert_eq!(r2.len(), 8);
        assert_eq!(r2[0], "@r1 1/2");
        assert_eq!(r2[7], "@r11 11/2");
    }

    #[test]
    fn by_suffix_rejects_a_record_without_a_marker() {
        let err = split_suffix(&["@a/1", "@b 1:N:0:AT", "@c/2"])
            .expect_err("a record with no marker must be rejected");
        assert!(
            err.to_string().contains("carries no /1 or /2 mate marker"),
            "unexpected error: {err}"
        );
    }
}
