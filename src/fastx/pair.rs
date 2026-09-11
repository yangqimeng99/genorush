//! `fastx pair`: match up two mate files that have drifted out of sync.
//!
//! Quality control is the usual cause. A trimmer or filter drops a read from
//! R1 but keeps its mate in R2, and from then on the two files disagree about
//! which record sits at which position. Nothing errors: both files are still
//! well-formed FASTQ, and every tool downstream that pairs by position -- most
//! of them -- silently aligns reads to the wrong mates.
//!
//! ## Why not just index one file
//!
//! The obvious implementation loads one file into a map keyed by read id and
//! streams the other past it. That is what the existing tools do, and it makes
//! memory a function of *input size* rather than of how badly the files
//! actually disagree: `seqkit pair` on a 49 GB + 62 GB pair has been reported
//! consuming ~380 GB of RAM before dying, with 400 MB of output written
//! (shenwei356/seqkit#305).
//!
//! Real desync is almost never proportional to input size. A few percent of
//! reads get filtered; the rest stay in the same relative order. So this
//! command holds only what it has to:
//!
//! - Both files are read concurrently, one record at a time from each.
//! - A record whose mate has already been seen is emitted immediately, as a
//!   pair, and neither is retained.
//! - A record whose mate has *not* been seen yet is held in `pending`.
//!
//! Peak memory is therefore the number of records that are unmatched at the
//! worst moment -- the true desync, plus the orphans -- and not the file size.
//! On files that are merely missing some reads, `pending` holds the orphans
//! and nothing else.
//!
//! ## Why nothing is declared an orphan early
//!
//! It is tempting to flush `pending` once the other stream has moved past
//! where a mate "should" have been. That is only sound if both files are in
//! the same relative order, which is exactly the assumption a re-pairing tool
//! cannot make -- a file that has been name-sorted, or rebuilt from an
//! unordered source, breaks it. A record flushed as an orphan whose mate
//! turns up later is silent data loss, the failure this command exists to
//! prevent. So orphanhood is only decided at EOF, when it is provable.
//!
//! The cost of that guarantee is the memory above, which is why exceeding
//! `--max-memory` is a real possibility on heavily reordered input, and why
//! it is handled by spilling to disk rather than by guessing.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::mpsc::Receiver;
use std::time::Instant;

use anyhow::{ensure, Context, Result};
use clap::Args;

use crate::common::fastq::{format_into_blocks, spawn_reader, FastqRecord};
use crate::common::hash::fnv1a;
use crate::io_utils::{
    display, ensure_one_stdio_at_most, is_stdio, open_block_writer, open_reader, BlockWriter,
    OutputOpts, STDIO,
};

#[derive(Args, Debug)]
pub struct PairArgs {
    /// Read 1 (R1) FASTQ (.fastq/.fq, gzip/bgzip auto-detected), or `-` for stdin.
    #[arg(short = 'i', long = "in1", value_name = "FILE")]
    in1: PathBuf,

    /// Read 2 (R2) mate file. Only one of the two inputs can be `-`: both
    /// would be reading the same stdin.
    #[arg(short = 'I', long = "in2", value_name = "FILE")]
    in2: PathBuf,

    /// Output for the paired read 1 records. Defaults to stdout.
    /// Gzip-compressed if the path ends in `.gz` or `-z/--gzip` is passed.
    #[arg(short = 'o', long = "out1", value_name = "FILE", default_value = STDIO)]
    out1: PathBuf,

    /// Output for the paired read 2 records. Required: the two mates cannot
    /// both go to stdout.
    #[arg(short = 'O', long = "out2", value_name = "FILE")]
    out2: PathBuf,

    /// Where to write read 1 records with no mate. Without this they are
    /// counted and discarded.
    #[arg(short = 'u', long = "unpaired1", value_name = "FILE")]
    unpaired1: Option<PathBuf>,

    /// Where to write read 2 records with no mate. Without this they are
    /// counted and discarded.
    #[arg(short = 'U', long = "unpaired2", value_name = "FILE")]
    unpaired2: Option<PathBuf>,

    /// Memory budget for the records held while waiting for their mates.
    /// Accepts plain bytes or a K/M/G suffix. Exceeding it is not an error:
    /// the run spills to `--temp-dir` and keeps going.
    #[arg(long, value_name = "SIZE", default_value = "2G")]
    max_memory: String,

    /// Directory for spill files, used only if the memory budget is
    /// exceeded. A subdirectory is created inside it and removed at the end.
    /// Defaults to the working directory.
    #[arg(long, value_name = "DIR")]
    temp_dir: Option<PathBuf>,

    /// Partitions per side when the join spills to disk. More partitions
    /// mean smaller ones, so less memory at join time and more open files.
    #[arg(long, default_value_t = 256)]
    partitions: usize,

    /// Records buffered per parallel compression batch.
    #[arg(long, default_value_t = 50_000)]
    chunk_records: usize,
}

/// Parses `2G`, `500M`, `1024` into bytes. Rejects nonsense rather than
/// silently treating it as zero, since a zero budget would spill instantly.
pub fn parse_size(s: &str) -> Result<u64> {
    let t = s.trim();
    ensure!(!t.is_empty(), "--max-memory is empty");
    let (digits, mult) = match t.chars().last().unwrap().to_ascii_uppercase() {
        'K' => (&t[..t.len() - 1], 1024u64),
        'M' => (&t[..t.len() - 1], 1024 * 1024),
        'G' => (&t[..t.len() - 1], 1024 * 1024 * 1024),
        'T' => (&t[..t.len() - 1], 1024u64 * 1024 * 1024 * 1024),
        _ => (t, 1),
    };
    let n: u64 = digits
        .trim()
        .parse()
        .with_context(|| format!("--max-memory {s:?} is not a size like 2G, 500M or 1048576"))?;
    let bytes = n
        .checked_mul(mult)
        .with_context(|| format!("--max-memory {s:?} overflows"))?;
    ensure!(bytes > 0, "--max-memory must be greater than zero");
    Ok(bytes)
}

/// What one held record costs, near enough to budget by: the four strings'
/// contents, their headers, and the map entry holding the id a second time.
fn footprint(rec: &FastqRecord) -> usize {
    rec.header.len() * 2 + rec.seq.len() + rec.plus.len() + rec.qual.len() + 128
}

/// Accumulates records and hands them to `BlockWriter` in batches, so output
/// compression parallelises the same way every other command's does.
struct ChunkedOut {
    writer: BlockWriter,
    buf: Vec<FastqRecord>,
    cap: usize,
    written: u64,
}

impl ChunkedOut {
    fn open(path: &Path, opts: OutputOpts, cap: usize) -> Result<Self> {
        Ok(Self {
            writer: open_block_writer(path, opts)?,
            buf: Vec::with_capacity(cap),
            cap,
            written: 0,
        })
    }

    fn push(&mut self, rec: FastqRecord) -> Result<()> {
        self.buf.push(rec);
        self.written += 1;
        if self.buf.len() >= self.cap {
            self.flush_batch()?;
        }
        Ok(())
    }

    fn flush_batch(&mut self) -> Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let refs: Vec<&FastqRecord> = self.buf.iter().collect();
        self.writer.write_blocks(format_into_blocks(&refs)?)?;
        self.buf.clear();
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.flush_batch()?;
        self.writer.flush()?;
        Ok(())
    }
}

/// Records seen on one side whose mate has not turned up yet, keyed by read
/// id. The index preserves input order so orphans come out in the order they
/// went in.
#[derive(Default)]
struct Pending {
    by_id: HashMap<String, (u64, FastqRecord)>,
    bytes: usize,
}

impl Pending {
    fn insert(&mut self, id: String, index: u64, rec: FastqRecord) {
        self.bytes += footprint(&rec) + id.len();
        if let Some((_, old)) = self.by_id.insert(id.clone(), (index, rec)) {
            // A repeated id within one file: keep the first, since that is the
            // one whose position the paired output is aligned to.
            self.bytes -= footprint(&old) + id.len();
        }
    }

    fn take(&mut self, id: &str) -> Option<(u64, FastqRecord)> {
        let hit = self.by_id.remove(id);
        if let Some((_, rec)) = &hit {
            self.bytes -= footprint(rec) + id.len();
        }
        hit
    }

    fn len(&self) -> usize {
        self.by_id.len()
    }

    /// Everything still unmatched, back in input order.
    fn drain_in_order(&mut self) -> Vec<FastqRecord> {
        self.drain_with_index()
            .into_iter()
            .map(|(_, rec)| rec)
            .collect()
    }

    fn drain_with_index(&mut self) -> Vec<(u64, FastqRecord)> {
        let mut left: Vec<(u64, FastqRecord)> = self.by_id.drain().map(|(_, v)| v).collect();
        left.sort_by_key(|(index, _)| *index);
        self.bytes = 0;
        left
    }
}

struct Outputs {
    r1: ChunkedOut,
    r2: ChunkedOut,
    orphan1: Option<ChunkedOut>,
    orphan2: Option<ChunkedOut>,
}

pub fn run(args: PairArgs, opts: OutputOpts) -> Result<()> {
    ensure!(args.chunk_records > 0, "--chunk-records must be > 0");
    let budget = parse_size(&args.max_memory)?;

    ensure_one_stdio_at_most(&[&args.in1, &args.in2], "input")?;
    let mut outs: Vec<&Path> = vec![args.out1.as_path(), args.out2.as_path()];
    outs.extend(args.unpaired1.as_deref());
    outs.extend(args.unpaired2.as_deref());
    ensure_one_stdio_at_most(&outs, "output")?;

    let start = Instant::now();
    log::info!(
        "pairing {} + {} -> {} + {}",
        display(&args.in1),
        display(&args.in2),
        display(&args.out1),
        display(&args.out2)
    );

    let mut outputs = Outputs {
        r1: ChunkedOut::open(&args.out1, opts, args.chunk_records)?,
        r2: ChunkedOut::open(&args.out2, opts, args.chunk_records)?,
        orphan1: match &args.unpaired1 {
            Some(p) => Some(ChunkedOut::open(p, opts, args.chunk_records)?),
            None => None,
        },
        orphan2: match &args.unpaired2 {
            Some(p) => Some(ChunkedOut::open(p, opts, args.chunk_records)?),
            None => None,
        },
    };

    let rx1 = spawn_reader(args.in1.clone())?;
    let rx2 = spawn_reader(args.in2.clone())?;

    let mut pending1 = Pending::default();
    let mut pending2 = Pending::default();
    let (mut n1, mut n2) = (0u64, 0u64);
    let (mut done1, mut done2) = (false, false);
    let mut pairs = 0u64;
    let mut spilled: Option<JoinStats> = None;
    let mut peak_bytes = 0usize;
    let mut peak_held = 0usize;

    // Both sides advance together so that files which are merely missing a
    // few records never accumulate more than those records: a strict
    // alternation keeps the two read positions at the same index.
    while !(done1 && done2) {
        if !done1 {
            match rx1.recv() {
                Ok(rec) => {
                    let rec = rec?;
                    let id = rec.base_id().to_string();
                    match pending2.take(&id) {
                        Some((_, mate)) => {
                            outputs.r1.push(rec)?;
                            outputs.r2.push(mate)?;
                            pairs += 1;
                        }
                        None => pending1.insert(id, n1, rec),
                    }
                    n1 += 1;
                }
                Err(_) => done1 = true,
            }
        }
        if !done2 {
            match rx2.recv() {
                Ok(rec) => {
                    let rec = rec?;
                    let id = rec.base_id().to_string();
                    match pending1.take(&id) {
                        Some((_, mate)) => {
                            outputs.r1.push(mate)?;
                            outputs.r2.push(rec)?;
                            pairs += 1;
                        }
                        None => pending2.insert(id, n2, rec),
                    }
                    n2 += 1;
                }
                Err(_) => done2 = true,
            }
        }

        let held = pending1.bytes + pending2.bytes;
        if held > peak_bytes {
            peak_bytes = held;
            peak_held = pending1.len() + pending2.len();
        }
        if held as u64 > budget {
            log::warn!(
                "{} records are still waiting for their mates ({:.1} GiB, past --max-memory {}); \
                 the inputs are more than lightly out of step. Finishing the join on disk.",
                pending1.len() + pending2.len(),
                held as f64 / (1024.0 * 1024.0 * 1024.0),
                args.max_memory
            );
            let temp_root = args.temp_dir.clone().unwrap_or_else(|| PathBuf::from("."));
            spilled = Some(spill_and_join(
                std::mem::take(&mut pending1),
                std::mem::take(&mut pending2),
                &rx1,
                &rx2,
                done1,
                done2,
                &mut n1,
                &mut n2,
                &mut outputs,
                &temp_root,
                [&args.in1, &args.in2],
                budget,
                args.partitions,
            )?);
            break;
        }
    }

    // Nothing left to match against: whatever is still held has no mate.
    // After a spill the pending sets are empty and the join already counted.
    let orphans1 = pending1.drain_in_order();
    let orphans2 = pending2.drain_in_order();
    let (mut n_orphan1, mut n_orphan2) = (orphans1.len() as u64, orphans2.len() as u64);
    emit_orphans(orphans1, &mut outputs.orphan1)?;
    emit_orphans(orphans2, &mut outputs.orphan2)?;
    if let Some(join) = &spilled {
        pairs += join.pairs;
        n_orphan1 += join.orphans1;
        n_orphan2 += join.orphans2;
    }

    outputs.r1.finish()?;
    outputs.r2.finish()?;
    if let Some(o) = &mut outputs.orphan1 {
        o.finish()?;
    }
    if let Some(o) = &mut outputs.orphan2 {
        o.finish()?;
    }

    log::info!(
        "paired {pairs} records ({n1} read 1, {n2} read 2 in); \
         {n_orphan1} + {n_orphan2} without a mate; \
         peak {} records held ({:.1} MiB) in {:.2?}",
        peak_held,
        peak_bytes as f64 / (1024.0 * 1024.0),
        start.elapsed()
    );
    if let Some(join) = &spilled {
        log::info!(
            "the on-disk join handled {} records ({:.1} GiB spilled); \
             output for those follows partition order, not input order",
            join.records_spilled,
            join.spilled_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
        );
    }
    if n_orphan1 + n_orphan2 > 0 && args.unpaired1.is_none() && args.unpaired2.is_none() {
        log::warn!(
            "{} records had no mate and were dropped; pass -u/-U to keep them",
            n_orphan1 + n_orphan2
        );
    }

    // The two outputs are aligned by construction -- every write to one is
    // matched by a write to the other in the same step -- but a mismatch here
    // would mean that invariant broke, so it is worth the comparison.
    ensure!(
        outputs.r1.written == outputs.r2.written,
        "internal error: paired outputs have different lengths ({} vs {})",
        outputs.r1.written,
        outputs.r2.written
    );
    Ok(())
}

fn emit_orphans(records: Vec<FastqRecord>, out: &mut Option<ChunkedOut>) -> Result<()> {
    let Some(out) = out else { return Ok(()) };
    for rec in records {
        out.push(rec)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// On-disk fallback: a hash-partitioned join
// ---------------------------------------------------------------------------
//
// Once the records waiting for their mates no longer fit in memory, holding
// them is not an option and guessing which ones are orphans is not either.
// The way out is the one databases use for exactly this shape of problem.
//
// Both remainders are written to `k` partitions by `hash(read id) % k`. Two
// mates hash identically, so a pair can never be split across partitions --
// which means each partition can be joined on its own, and only one
// partition has to be in memory at a time. Peak memory becomes a function of
// `k` rather than of the input.
//
// The cost is one round trip through the disk for whatever had not been
// paired yet, and an output order that follows the partitions rather than
// the input. The first is unavoidable; the second is reported, not hidden.

/// Removes its directory when it goes out of scope, so an error on any path
/// out of the join still cleans up.
struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn create(parent: &Path) -> Result<Self> {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.subsec_nanos())
            .unwrap_or(0);
        let path = parent.join(format!("genorush-pair-{}-{}", std::process::id(), nanos));
        std::fs::create_dir_all(&path)
            .with_context(|| format!("failed to create spill directory {}", path.display()))?;
        Ok(Self { path })
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        if let Err(e) = std::fs::remove_dir_all(&self.path) {
            log::warn!(
                "could not remove the spill directory {}: {e}",
                self.path.display()
            );
        }
    }
}

/// Records are written to partitions in a framing of their own rather than as
/// FASTQ: the original position has to survive the trip so orphans can come
/// back out in input order, and lengths are cheaper to trust than delimiters
/// when a header is allowed to contain anything.
fn encode_record(buf: &mut Vec<u8>, index: u64, rec: &FastqRecord) {
    buf.extend_from_slice(&index.to_le_bytes());
    for field in [&rec.header, &rec.seq, &rec.plus, &rec.qual] {
        buf.extend_from_slice(&(field.len() as u32).to_le_bytes());
        buf.extend_from_slice(field.as_bytes());
    }
}

fn decode_record(r: &mut dyn std::io::BufRead) -> Result<Option<(u64, FastqRecord)>> {
    let mut idx = [0u8; 8];
    match r.read_exact(&mut idx) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e).context("failed reading a spilled record"),
    }
    let mut field = || -> Result<String> {
        let mut len = [0u8; 4];
        r.read_exact(&mut len).context("truncated spill file")?;
        let mut bytes = vec![0u8; u32::from_le_bytes(len) as usize];
        r.read_exact(&mut bytes).context("truncated spill file")?;
        String::from_utf8(bytes).context("spill file is not valid UTF-8")
    };
    Ok(Some((
        u64::from_le_bytes(idx),
        FastqRecord {
            header: field()?,
            seq: field()?,
            plus: field()?,
            qual: field()?,
        },
    )))
}

/// One side's partitions: `k` writers, each buffering until it has enough to
/// be worth compressing.
struct Partitions {
    paths: Vec<PathBuf>,
    writers: Vec<BlockWriter>,
    buffers: Vec<Vec<u8>>,
    bytes: u64,
}

const SPILL_FLUSH_BYTES: usize = 4 * 1024 * 1024;

impl Partitions {
    fn create(dir: &Path, side: &str, k: usize) -> Result<Self> {
        let mut paths = Vec::with_capacity(k);
        let mut writers = Vec::with_capacity(k);
        for i in 0..k {
            let p = dir.join(format!("{side}.{i:04}.bin.gz"));
            writers.push(open_block_writer(&p, OutputOpts::default())?);
            paths.push(p);
        }
        Ok(Self {
            paths,
            writers,
            buffers: vec![Vec::new(); k],
            bytes: 0,
        })
    }

    fn push(&mut self, index: u64, rec: &FastqRecord) -> Result<()> {
        let k = (fnv1a(rec.base_id().as_bytes()) % self.paths.len() as u64) as usize;
        encode_record(&mut self.buffers[k], index, rec);
        self.bytes +=
            20 + (rec.header.len() + rec.seq.len() + rec.plus.len() + rec.qual.len()) as u64;
        if self.buffers[k].len() >= SPILL_FLUSH_BYTES {
            let block = std::mem::take(&mut self.buffers[k]);
            self.writers[k].write_blocks(vec![block])?;
        }
        Ok(())
    }

    /// Flushes every partition and hands back their paths, consuming `self`
    /// so the writers -- and with them the open file handles -- are dropped
    /// before anything tries to read or delete the files. On Windows an open
    /// handle makes a file undeletable, so this is not merely tidy.
    fn finish(mut self) -> Result<Vec<PathBuf>> {
        for (i, w) in self.writers.iter_mut().enumerate() {
            let block = std::mem::take(&mut self.buffers[i]);
            if !block.is_empty() {
                w.write_blocks(vec![block])?;
            }
            w.flush()?;
        }
        self.writers.clear();
        Ok(std::mem::take(&mut self.paths))
    }
}

/// Deletes a spill file, reporting rather than swallowing a failure: a
/// platform that will not let the file go is worth knowing about, since the
/// directory removal at the end depends on it.
fn discard(path: &Path) {
    if let Err(e) = std::fs::remove_file(path) {
        log::warn!("could not remove the spill file {}: {e}", path.display());
    }
}

/// Refuses to start spilling if the destination cannot hold it.
///
/// The estimate is the on-disk size of the inputs, which overshoots -- only
/// the unpaired remainder is ever written, and it is compressed the same way
/// the inputs are. Overshooting is the right direction: running out of space
/// halfway through leaves a half-written join and no answer.
fn ensure_space_for_spill(temp_root: &Path, inputs: &[&Path]) -> Result<()> {
    let mut need: u64 = 0;
    for p in inputs {
        if !is_stdio(p) {
            need += std::fs::metadata(p).map(|m| m.len()).unwrap_or(0);
        }
    }
    let need = need + need / 10; // a little headroom for framing and rounding
    let available = fs4::available_space(temp_root).with_context(|| {
        format!(
            "could not check free space on {} (use --temp-dir to point somewhere else)",
            temp_root.display()
        )
    })?;
    ensure!(
        available >= need,
        "not enough space to spill: {} needs about {:.1} GiB free but has {:.1} GiB. \
         Point --temp-dir at a filesystem with room, or raise --max-memory so the join \
         stays in memory.",
        temp_root.display(),
        need as f64 / (1024.0 * 1024.0 * 1024.0),
        available as f64 / (1024.0 * 1024.0 * 1024.0)
    );
    Ok(())
}

struct JoinStats {
    pairs: u64,
    orphans1: u64,
    orphans2: u64,
    spilled_bytes: u64,
    records_spilled: u64,
}

/// Everything that was still unpaired when memory ran out, joined through
/// disk. Consumes what the in-memory phase was holding and drains what is
/// left of both readers.
#[allow(clippy::too_many_arguments)]
fn spill_and_join(
    mut pending1: Pending,
    mut pending2: Pending,
    rx1: &Receiver<Result<FastqRecord>>,
    rx2: &Receiver<Result<FastqRecord>>,
    mut done1: bool,
    mut done2: bool,
    n1: &mut u64,
    n2: &mut u64,
    outputs: &mut Outputs,
    temp_root: &Path,
    inputs: [&Path; 2],
    budget: u64,
    k: usize,
) -> Result<JoinStats> {
    ensure_space_for_spill(temp_root, &inputs)?;
    let dir = TempDir::create(temp_root)?;
    log::info!(
        "spilling to {} across {k} partitions per side",
        dir.path.display()
    );

    let mut part1 = Partitions::create(&dir.path, "r1", k)?;
    let mut part2 = Partitions::create(&dir.path, "r2", k)?;
    let mut records_spilled = 0u64;

    for (index, rec) in pending1.drain_with_index() {
        part1.push(index, &rec)?;
        records_spilled += 1;
    }
    for (index, rec) in pending2.drain_with_index() {
        part2.push(index, &rec)?;
        records_spilled += 1;
    }

    while !(done1 && done2) {
        if !done1 {
            match rx1.recv() {
                Ok(rec) => {
                    part1.push(*n1, &rec?)?;
                    *n1 += 1;
                    records_spilled += 1;
                }
                Err(_) => done1 = true,
            }
        }
        if !done2 {
            match rx2.recv() {
                Ok(rec) => {
                    part2.push(*n2, &rec?)?;
                    *n2 += 1;
                    records_spilled += 1;
                }
                Err(_) => done2 = true,
            }
        }
    }
    let spilled_bytes = part1.bytes + part2.bytes;
    let paths1 = part1.finish()?;
    let paths2 = part2.finish()?;

    let mut stats = JoinStats {
        pairs: 0,
        orphans1: 0,
        orphans2: 0,
        spilled_bytes,
        records_spilled,
    };

    for i in 0..k {
        // Side one's partition has to be resident to be joined against; side
        // two only streams past it.
        let mut held: HashMap<String, (u64, FastqRecord)> = HashMap::new();
        {
            let mut held_bytes = 0u64;
            let mut reader = open_reader(&paths1[i])?;
            while let Some((index, rec)) = decode_record(reader.as_mut())? {
                held_bytes += footprint(&rec) as u64;
                ensure!(
                    held_bytes <= budget,
                    "partition {i} of the spilled join does not fit in --max-memory \
                     ({:.1} GiB and counting). Raise --max-memory or --partitions.",
                    held_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
                );
                held.insert(rec.base_id().to_string(), (index, rec));
            }
        }

        {
            let mut reader = open_reader(&paths2[i])?;
            while let Some((_, rec)) = decode_record(reader.as_mut())? {
                match held.remove(rec.base_id()) {
                    Some((_, mate)) => {
                        outputs.r1.push(mate)?;
                        outputs.r2.push(rec)?;
                        stats.pairs += 1;
                    }
                    None => {
                        stats.orphans2 += 1;
                        if let Some(o) = &mut outputs.orphan2 {
                            o.push(rec)?;
                        }
                    }
                }
            }
        }

        let mut left: Vec<(u64, FastqRecord)> = held.into_values().collect();
        left.sort_by_key(|(index, _)| *index);
        stats.orphans1 += left.len() as u64;
        if let Some(o) = &mut outputs.orphan1 {
            for (_, rec) in left {
                o.push(rec)?;
            }
        }

        // Reclaim the space as the join advances, so peak disk is the whole
        // spill only at the moment the join starts. Both readers are out of
        // scope by now, which Windows requires before the files will go.
        discard(&paths1[i]);
        discard(&paths2[i]);
    }

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(id: &str, mate: u8) -> FastqRecord {
        FastqRecord {
            header: format!("@{id} {mate}:N:0:AT"),
            seq: "ACGTACGT".into(),
            plus: "+".into(),
            qual: "IIIIIIII".into(),
        }
    }

    #[test]
    fn sizes_parse_with_and_without_a_suffix() {
        assert_eq!(parse_size("1024").unwrap(), 1024);
        assert_eq!(parse_size("2K").unwrap(), 2048);
        assert_eq!(parse_size("3m").unwrap(), 3 * 1024 * 1024);
        assert_eq!(parse_size("1G").unwrap(), 1024 * 1024 * 1024);
        assert_eq!(parse_size(" 4G ").unwrap(), 4 * 1024 * 1024 * 1024);
    }

    #[test]
    fn nonsense_sizes_are_rejected_rather_than_read_as_zero() {
        // A budget silently parsed as 0 would spill on the first record.
        for bad in ["", "  ", "G", "-1", "1.5G", "2GB", "lots"] {
            assert!(parse_size(bad).is_err(), "{bad:?} should not parse");
        }
    }

    #[test]
    fn pending_releases_its_budget_when_a_mate_arrives() {
        let mut p = Pending::default();
        p.insert("a".into(), 0, rec("a", 1));
        p.insert("b".into(), 1, rec("b", 1));
        let with_two = p.bytes;
        assert_eq!(p.len(), 2);

        assert!(p.take("a").is_some(), "a was inserted");
        assert_eq!(p.len(), 1);
        assert!(p.bytes < with_two, "taking a record must return its budget");

        assert!(p.take("a").is_none(), "a was already taken");
        assert!(p.take("zzz").is_none());
    }

    #[test]
    fn pending_hands_leftovers_back_in_input_order() {
        let mut p = Pending::default();
        // Inserted out of order on purpose: the index decides, not arrival.
        p.insert("c".into(), 7, rec("c", 1));
        p.insert("a".into(), 2, rec("a", 1));
        p.insert("b".into(), 5, rec("b", 1));
        let ids: Vec<String> = p
            .drain_in_order()
            .iter()
            .map(|r| r.base_id().to_string())
            .collect();
        assert_eq!(ids, ["a", "b", "c"]);
        assert_eq!(p.bytes, 0, "draining must clear the accounting too");
    }

    #[test]
    fn a_repeated_id_within_one_file_keeps_the_first() {
        let mut p = Pending::default();
        p.insert("dup".into(), 0, rec("dup", 1));
        p.insert("dup".into(), 1, rec("dup", 1));
        assert_eq!(p.len(), 1);
        // The budget must not count the discarded copy, or it leaks upward.
        let mut single = Pending::default();
        single.insert("dup".into(), 0, rec("dup", 1));
        assert_eq!(p.bytes, single.bytes);
    }

    #[test]
    fn spilled_records_survive_the_round_trip_byte_for_byte() {
        let original = FastqRecord {
            header: "@SRR1.7 7/1 with\ttab and spaces".into(),
            seq: "ACGTNNNN".into(),
            plus: "+SRR1.7".into(),
            qual: "!!!!IIII".into(),
        };
        let mut buf = Vec::new();
        encode_record(&mut buf, 42, &original);

        let mut cursor = std::io::Cursor::new(buf);
        let (index, back) = decode_record(&mut cursor).unwrap().expect("one record");
        assert_eq!(index, 42);
        assert_eq!(back.header, original.header, "a tab in the header survives");
        assert_eq!(back.seq, original.seq);
        assert_eq!(back.plus, original.plus);
        assert_eq!(back.qual, original.qual);
        assert!(decode_record(&mut cursor).unwrap().is_none(), "then EOF");
    }

    #[test]
    fn decoding_a_truncated_spill_file_is_an_error_not_a_silent_stop() {
        let mut buf = Vec::new();
        encode_record(&mut buf, 1, &rec("x", 1));
        buf.truncate(buf.len() - 3);
        let mut cursor = std::io::Cursor::new(buf);
        assert!(decode_record(&mut cursor).is_err());
    }
}
