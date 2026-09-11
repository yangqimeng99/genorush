# GenoRush

[中文说明](README.zh.md)

A fast, natively multi-threaded, cross-platform command-line toolkit for
bioinformatics data wrangling, written in Rust. Built in the spirit of
[seqkit](https://github.com/shenwei356/seqkit): a single static binary,
one command per task, no runtime dependencies.

Status: early stage, actively growing. A few commands exist today; more are
planned under new categories (`vcf`, `sv`, ...) as needs come up.

## Install

Prebuilt binaries for Linux (static musl, works on any distro/glibc version), macOS
(Apple Silicon), and Windows are attached to
[GitHub Releases](https://github.com/yangqimeng99/genorush/releases).
Download the archive for your platform, extract, and put `genorush` (or
`genorush.exe`) on your `PATH`.

Or build from source (requires the [Rust toolchain](https://rustup.rs)):

```bash
git clone https://github.com/yangqimeng99/genorush.git
cd genorush
cargo build --release
./target/release/genorush --help
```

## Commands

```
genorush <category> <action> [options]
```

| Category | Action   | Does |
|----------|----------|------|
| `fastx`  | `rename` | Rename sequence names in a FASTA file via a mapping table |
| `gff`    | `rename` | Rename the seqid column in a GFF/GTF file via a mapping table |
| `fastx`  | `sample` | Downsample FASTQ reads by proportion or exact count, single- or paired-end |
| `fastx`  | `rescue` | Recover the leading run of clean reads from a truncated/corrupted FASTQ, single- or paired-end |
| `fastx`  | `interleave` | Merge R1/R2 into a single standard interleaved FASTQ |
| `fastx`  | `deinterleave` | Split a merged FASTQ back into R1/R2, by position (interleaved or R1-then-R2 concatenated) or by the `/1`+`/2` header markers |
| `fastx`  | `cat` | Concatenate FASTQ from repeated sequencing runs, checking for duplicate read IDs |
| `fastx`  | `pair` | Match up two mate files that have drifted out of sync, in memory proportional to the drift |

Every subcommand accepts a global `-j/--threads` flag (default: `1`; pass
`0` to use all logical cores).

### Pipes

**Output goes to stdout unless you name a file**, and every input accepts
`-` for stdin, so commands compose with the rest of a pipeline:

```bash
# no -o: sampled reads go straight to an aligner, no temporary file
genorush fastx sample -i reads.fq.gz -p 0.1 -s 42 | bwa mem ref.fa -

# interleave and align in one go
genorush fastx interleave -i R1.fq.gz -I R2.fq.gz | bwa mem -p ref.fa -

# read a merged FASTQ off a pipe, write both mates to files
zcat merged.fq.gz | genorush fastx deinterleave -i - -o R1.fq.gz -O R2.fq.gz

# both ends piped, output gzip-compressed
cat genome.fa | genorush fastx rename - -n map.tsv -z > renamed.fa.gz
```

`-o -` is still accepted and means exactly what omitting it does. Commands
with two outputs (`fastx deinterleave`, and the paired-end modes of
`sample`/`rescue`/`cat`) default only their first output to stdout; the
second has to be a file, since two record streams cannot share one pipe.

### BGZF

`-b/--bgzf` writes BGZF instead of ordinary gzip:

```bash
genorush fastx rename genome.fa -n map.tsv -o genome.renamed.fa.gz --bgzf
samtools faidx genome.renamed.fa.gz        # works
```

BGZF *is* gzip — anything that reads `.gz` reads this — but written so it can
be indexed: bounded members that record their own compressed size, and a
final empty block that proves the stream is whole. Without it, `samtools
faidx` on the same content says *"Cannot index files compressed with gzip,
please use bgzip"*, and `tabix` will say the same about a `.vcf.gz`. A file
that decompresses perfectly and cannot be indexed is the kind of dead end
worth its own flag.

It costs a little: members are capped at 64 KiB, so compression is marginally
worse than the large blocks `-z` uses. Take it when something downstream will
want to seek.

Compressed *input* needs no flag: gzip is recognised from the data itself,
so a pipe carrying `.gz` bytes just works. Compressed *output* can't be
inferred the same way — a pipe has no `.gz` extension to inspect — so
`-z/--gzip` is how a pipeline asks for it. For file outputs nothing
changes: a path ending in `.gz` is still compressed on its own.

Two consequences of there being exactly one stdin and one stdout:

- At most one input per command may be `-`, and at most one output. Two
  inputs reading the same stdin would each get an arbitrary half of it; two
  outputs sharing one pipe would interleave two record streams. Both are
  rejected rather than silently producing plausible garbage.
- `fastx deinterleave` needs to read its input twice for `--layout concat`
  (to find the midpoint) and for `--layout auto` when the head of the file
  is inconclusive. A pipe can only be read once, so those combinations are
  refused with a message pointing at the single-pass layouts
  (`--layout interleaved`, `--layout by-suffix`). Everything else, this
  command included, is single-pass and pipes fine.

A downstream that stops reading (`... | head`) ends the run cleanly rather
than reporting a broken-pipe failure. Gzip output aimed at a terminal is
refused, since defaulting to stdout makes a forgotten redirect easy and
binary on a terminal is never what was wanted; plain text to a terminal is
left alone, because that is how you look at a few records.

### `fastx rename` / `gff rename`

```bash
genorush fastx rename genome.fa  -n name_map.tsv -o renamed.fa
genorush gff   rename genes.gff  -n name_map.tsv -o renamed.gff.gz
```

Gzip/bgzip input is auto-detected by content, not by file extension.
Output is gzip-compressed automatically when the output path ends in
`.gz`. See [`docs/en/rename.md`](docs/en/rename.md) for the full design
writeup, including exactly how this compares to the Python script it
replaces.

### `fastx sample`

```bash
# single-end (e.g. long reads), sample by proportion or exact count
genorush fastx sample -i reads.fq.gz -p 0.1   -o sub.fq.gz -s 42
genorush fastx sample -i reads.fq.gz -n 50000 -o sub.fq.gz -s 42

# paired-end, sampled together in one pass — R1/R2 are always kept in sync
genorush fastx sample -i R1.fq.gz -I R2.fq.gz -o R1.sub.fq.gz -O R2.sub.fq.gz -p 0.1 -s 42
```

Unlike `seqkit sample`, which has no paired-end mode (you run it twice and
rely on passing the same seed to both invocations), this command reads
both mates in one process and samples pairs atomically, validating along
the way that R1/R2 read counts and IDs actually correspond. See
[`docs/en/sample.md`](docs/en/sample.md) for the full algorithm writeup
(deterministic parallel proportion sampling, single-pass reservoir
sampling for exact counts, and why each beats the naive approach).

### `fastx rescue`

```bash
# single-end: recover clean reads from a truncated/corrupted download
genorush fastx rescue -i reads.fq.gz -o rescued.fq.gz

# paired-end: recovers only pairs where both mates are intact and match
genorush fastx rescue -i R1.fq.gz -I R2.fq.gz -o R1.rescued.fq.gz -O R2.rescued.fq.gz
```

For interrupted downloads: everything decoded before the point of
corruption is still good data, and this command recovers exactly that,
stopping cleanly instead of erroring out. Exit code distinguishes a fully
clean read (`0`) from a partial rescue (`3`) from nothing salvageable
(`1`), so it composes into scripts. See
[`docs/en/rescue.md`](docs/en/rescue.md) for the full design writeup.

### `fastx interleave` / `fastx deinterleave`

```bash
genorush fastx interleave -i R1.fq.gz -I R2.fq.gz -o merged.fq.gz

# Splitting back apart. --layout auto (the default) probes the first few
# thousand records and picks a strategy; the other values skip the probe.
genorush fastx deinterleave -i merged.fq.gz -o R1.fq.gz -O R2.fq.gz
```

`--layout` in practice:

```bash
# by-suffix: route each record by the /1 or /2 marker in its own header,
# ignoring position entirely. This is the mode for files where position
# says nothing -- e.g. SRA-derived FASTQ whose mates come in several runs
# rather than two, and whose global read numbering gives the two mates of
# a pair different ids (@SRR17458599.1 1/2 mated to @SRR17458599.24174090
# 24174090/1). Single pass, constant memory.
genorush fastx deinterleave -i merged.fq.gz --layout by-suffix \
    -o R1.fq.gz -O R2.fq.gz -j 8

# interleaved: R1,R2,R1,R2,... Single pass. Each pair is verified to share
# a read id as it is written.
genorush fastx deinterleave -i merged.fq.gz --layout interleaved \
    -o R1.fq.gz -O R2.fq.gz

# concat: all R1 records, then all R2 records (`cat R1.fq R2.fq`). Needs
# the midpoint, so it counts records in a cheap first pass. Records that
# carry mate markers are checked against their half as they are written.
genorush fastx deinterleave -i merged.fq.gz --layout concat \
    -o R1.fq.gz -O R2.fq.gz

# --no-pair-check turns off those write-time checks for the positional
# modes. Reach for it only when your headers don't follow the standard
# /1+/2 or Illumina 1:...+2:... conventions and the checks cry wolf.
genorush fastx deinterleave -i merged.fq.gz --layout interleaved \
    --no-pair-check -o R1.fq.gz -O R2.fq.gz
```

A run that trips one of those checks stops with the offending record
number, names the mode that would handle the file, and exits non-zero; its
partial outputs must be discarded.

`fastx deinterleave` doesn't assume a merged file is properly interleaved:
`cat R1.fastq R2.fastq > merged.fastq` is common in the wild and is a
completely different byte layout that a naive splitter would silently get
wrong. Nor does it assume the file is positional at all — SRA-derived files
turn up with the mates in several runs rather than two, and with a global
read numbering that gives the two mates of a pair different ids, so nothing
but each record's own `/1`/`/2` marker can route them. `--layout auto`
probes the head of the file to pick a strategy, and whichever splitter runs
re-checks its assumption on every record as it writes. See
[`docs/en/interleave.md`](docs/en/interleave.md) for the detection
algorithm and the trade-offs.

### `fastx pair`

```bash
# R1 and R2 filtered independently are no longer aligned by position
genorush fastx pair -i R1.fq.gz -I R2.fq.gz \
    -o R1.paired.fq.gz -O R2.paired.fq.gz \
    -u R1.orphans.fq.gz -U R2.orphans.fq.gz -j 8
```

Quality control breaks pairing silently: a read is dropped from one mate file
and kept in the other, both files stay well-formed FASTQ, and everything
downstream that pairs by position starts aligning reads to the wrong mates.

Memory here is proportional to how far apart the files have drifted, not to
their size — a record is written out and released the moment its mate turns
up, so files that merely lost a few reads hold almost nothing. The comparable
tools index one whole file instead: `seqkit pair` has been reported taking
~380 GB of RAM on a 49 GB + 62 GB pair before dying. When the drift genuinely
is large (a reordered file), the join finishes through disk partitions under
`--max-memory` rather than growing without bound. See
[`docs/en/pair.md`](docs/en/pair.md).

### `fastx cat`

```bash
genorush fastx cat --r1 run1_R1.fq.gz --r1 run2_R1.fq.gz \
                    --r2 run1_R2.fq.gz --r2 run2_R2.fq.gz \
                    -o merged_R1.fq.gz -O merged_R2.fq.gz
```

For concatenating repeated sequencing runs of the same sample. Unlike
plain `cat`, this checks for duplicate read IDs across the inputs as it
streams through and aborts with the specific files/positions involved --
catching the realistic failure mode (the same file accidentally listed
twice) instead of silently doubling coverage. See
[`docs/en/cat.md`](docs/en/cat.md).

## Design notes for contributors

- `src/main.rs` wires a two-level `clap` command tree:
  `genorush <category> <action>`. Each category (`fastx/`, `gff/`, ...) is
  a module with a `mod.rs` that owns a `Subcommand` enum and a `run()`
  dispatcher; each action is its own file.
- `src/common/` holds logic shared across categories: `rename.rs` (the
  chunked-parallel line-transform engine), `fastq.rs` (a minimal FASTQ
  record model plus the concurrent-mate-reading/pairing infrastructure and
  parallel-block formatting shared by `sample`, `rescue`, `interleave`,
  `deinterleave`, and `cat`), `rng.rs` (a dependency-free SplitMix64 RNG,
  both a stateless index-keyed variant for parallel sampling and a
  stateful variant for sequential algorithms like reservoir sampling),
  `hash.rs` (a small FNV-1a hash used to compare/deduplicate read IDs
  across huge inputs without keeping full ID strings in memory --
  `deinterleave`'s layout detection and `cat`'s duplicate-ID check).
- `src/io_utils.rs` provides transparent gzip/bgzip-aware readers and
  writers used by every command — detect by magic bytes on read, by `.gz`
  extension on write. `BlockWriter` is the batch-oriented writer used by
  chunk-processing commands: it compresses multiple blocks into
  independent gzip members in parallel via `-j`/rayon (the same
  multi-member technique `pigz` uses), since standard gzip decompression
  can't be parallelized for a single stream but compressing data this tool
  generates itself can be.
- Every command ships with unit tests for its non-trivial shared logic
  (`cargo test`) and is clippy-clean (`cargo clippy --all-targets`).
- `tests/cli_pipes.rs` drives the built binary through real pipes, since
  argument defaults, exit codes and the data/log split only exist at that
  level. Nearly every case runs a command twice -- once to a file, once
  through a pipe -- and requires the two results to be byte-identical,
  because a pipeline that quietly differs from the file-based run is the
  failure worth catching. The helpers in `tests/common/` build fixtures and
  spawn processes in Rust rather than shelling out, so the suite runs on
  Windows as well as Linux and macOS.
- Full design rationale per command lives under `docs/en/` (English) and
  `docs/zh/` (Chinese, primary author's working language) — read those
  before extending a command, they document *why*, not just *what*.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).

## License

MIT, see [LICENSE](LICENSE).
