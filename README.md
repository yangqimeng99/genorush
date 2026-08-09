# GenoRush

[中文说明](README.zh.md)

A fast, natively multi-threaded, cross-platform command-line toolkit for
bioinformatics data wrangling, written in Rust. Built in the spirit of
[seqkit](https://github.com/shenwei356/seqkit): a single static binary,
one command per task, no runtime dependencies.

Status: early stage, actively growing. A few commands exist today; more are
planned under new categories (`vcf`, `sv`, ...) as needs come up.

## Install

Prebuilt binaries for Linux (static musl, works on any distro/glibc version), macOS (Intel and
Apple Silicon), and Windows are attached to
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

Every subcommand accepts a global `-j/--threads` flag (default: `1`; pass
`0` to use all logical cores).

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

## Design notes for contributors

- `src/main.rs` wires a two-level `clap` command tree:
  `genorush <category> <action>`. Each category (`fastx/`, `gff/`, ...) is
  a module with a `mod.rs` that owns a `Subcommand` enum and a `run()`
  dispatcher; each action is its own file.
- `src/common/` holds logic shared across categories: `rename.rs` (the
  chunked-parallel line-transform engine), `fastq.rs` (a minimal FASTQ
  record model plus the concurrent-mate-reading/pairing infrastructure
  shared by `sample` and `rescue`), `rng.rs` (a dependency-free SplitMix64
  RNG, both a stateless index-keyed variant for parallel sampling and a
  stateful variant for sequential algorithms like reservoir sampling).
- `src/io_utils.rs` provides transparent gzip/bgzip-aware readers and
  writers used by every command — detect by magic bytes on read, by `.gz`
  extension on write. `BlockWriter` is the batch-oriented writer used by
  chunk-processing commands (`fastx sample`): it compresses multiple
  blocks into independent gzip members in parallel via `-j`/rayon (the
  same multi-member technique `pigz` uses), since standard gzip
  decompression can't be parallelized for a single stream but compressing
  data this tool generates itself can be.
- Every command ships with unit tests for its non-trivial shared logic
  (`cargo test`) and is clippy-clean (`cargo clippy --all-targets`).
- Full design rationale per command lives under `docs/en/` (English) and
  `docs/zh/` (Chinese, primary author's working language) — read those
  before extending a command, they document *why*, not just *what*.

## License

MIT, see [LICENSE](LICENSE).
