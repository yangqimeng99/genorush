# Changelog

All notable changes to this project are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.3.0] - 2026-08-11

### Added

- `fastx interleave`: merge R1/R2 into a single standard interleaved FASTQ.
- `fastx deinterleave`: split a merged FASTQ back into R1/R2, auto-detecting
  whether the input is properly interleaved (`R1,R2,R1,R2,...`) or a naive
  `cat R1 R2`-style concatenation (all R1 records, then all R2) -- these
  are different byte layouts, not the same format, and a splitter that
  assumes one and gets the other silently fabricates broken pairs.
  Refuses to guess when neither hypothesis holds cleanly; `--layout` skips
  detection when the layout is already known. See
  `docs/en/interleave.md` / `docs/zh/interleave.md`.
- `fastx cat`: concatenate FASTQ files from repeated sequencing runs of the
  same sample (multiple lanes/flowcells), hashing every read ID as it
  streams through and aborting with the specific files/positions involved
  if a duplicate turns up -- catches the realistic failure mode (the same
  file accidentally included twice) that plain `cat` has no way to detect.
  Paired-end mode also re-verifies R1/R2 pairing within each source file
  pair. See `docs/en/cat.md` / `docs/zh/cat.md`.
- `common::hash`: a small, dependency-free FNV-1a hash, shared by
  `deinterleave`'s layout detection and `cat`'s duplicate-ID check for
  comparing read IDs across huge inputs without keeping every ID string in
  memory.
- `common::fastq::format_into_blocks` (moved out of `fastx sample`,
  unchanged behavior) is now shared by every command that buffers a batch
  of records for parallel-compressed `BlockWriter` output.

### Changed

- Release workflow no longer builds `macos-x86_64` (Intel, the `macos-13`
  runner): GitHub's shared runner pool for that specific image queued the
  job for 5+ minutes on two consecutive releases while every other
  platform (including `macos-latest`/arm64) started within seconds. Only
  Apple Silicon macOS binaries are published now; the v0.2.0 release may
  still end up with a bonus Intel build if that already-queued job
  eventually completes, but future releases won't wait on it.
- **`fastx rename`/`gff rename`/`fastx rescue` gzip output is now
  compressed in parallel** via `io_utils::BlockWriter`, closing the same
  single-threaded-compression gap found and fixed in `fastx sample`
  (v0.2.0) but never propagated to these three pre-existing commands --
  `-j` previously had no effect on their write side at all. `rescue` also
  gained `--chunk-records` (buffers confirmed-good records in batches
  instead of writing one at a time) as part of this migration, and no
  longer needs the explicit `drop(writer)` workaround for the
  `std::process::exit`-skips-destructors hazard documented in
  `docs/en/rescue.md` -- `BlockWriter` finishes each gzip member
  synchronously within `write_blocks` rather than deferring to `Drop`, so
  the hazard doesn't apply to it.
- **`fastx deinterleave` now reads via `common::fastq::spawn_reader`**
  (a background thread decompressing/parsing into a channel) instead of
  `io_utils::open_reader` directly, matching every other multi-file
  command. For `split_interleaved`/`split_concat` this is a real fix: the
  reader thread now keeps decompressing the next chunk while the current
  chunk's parallel compression runs, instead of the main thread sitting
  idle during compression. For `detect_layout`/`count_records` it's
  mainly a consistency/code-simplification change -- those two functions
  do so little non-reading work that the wall-clock benefit is negligible,
  since decompression was already the sole bottleneck either way.
- `io_utils::open_writer`/`write_lines` removed (dead code once the above
  three commands migrated off them).
- Every subcommand's `--help` (long form) now states a recommended
  `-j/--threads` value, derived from measuring this machine's actual
  single-thread decompression throughput (~159 MB/s) against gzip
  compression throughput at increasing thread counts (~12 MB/s at `-j 1`
  scaling to ~56 MB/s at `-j 8`, ~77 MB/s at `-j 12`, with sharply
  diminishing returns past `-j 8`). Since decompression is hard-capped at
  1-2 threads regardless of `-j` (a single gzip stream can't be
  parallelized) but compression scales with it, the useful ceiling is
  command-specific: `-j 8` for commands where write volume is close to
  read volume (`rename`/`gff rename`/`rescue`/`interleave`/`deinterleave`/
  `cat`, where compression remains the bottleneck even at high thread
  counts on the hardware measured), `-j 1-2` for `fastx sample` at
  typical light sampling proportions (write volume scales with `-p`/`-n`,
  so less compression work needs fewer threads to stop being the
  bottleneck), `-j 8` for `sample` when keeping most of the input.
  Exact numbers are hardware-dependent; the qualitative shape (capped
  read side, command-dependent write-side ceiling) is not.

## [0.2.0] - 2026-08-09

### Added

- `fastx rescue`: recovers the leading run of clean, well-formed FASTQ
  records (or read pairs) from a truncated/corrupted file, e.g. an
  interrupted download. Stops cleanly at the first read error, R1/R2
  count mismatch, or ID mismatch instead of erroring out, and reports
  the outcome via exit code (`0` = fully clean, `3` = partial rescue,
  `1` = nothing salvageable). See `docs/en/rescue.md` /
  `docs/zh/rescue.md`.
- `CONTRIBUTORS.md`.

### Changed

- **`-j/--threads` now defaults to `1`** instead of all logical cores,
  so parallelism is opt-in (`-j 0` for all cores, or `-j N`).
- **`fastx sample` gzip output is now compressed in parallel**
  (`io_utils::BlockWriter`, a multi-member gzip writer in the same
  spirit as `pigz`), scaling with `-j`. Fixes a real-world case where
  `-j` had no measurable effect on a 57+59 GB paired FASTQ dataset —
  root cause was single-threaded output compression alternating with
  (non-overlapping) input decompression as two non-concurrent phases.
  Benchmarked at a 4.1x wall-clock speedup (`-j 8` vs `-j 1`) on a
  compression-heavy synthetic workload; output is byte-identical
  regardless of thread count. See `docs/en/sample.md` /
  `docs/zh/sample.md` for the full writeup.
- `fastx sample`'s concurrent-mate-reading and pairing logic
  (`spawn_reader`, `PairStep`/`recv_pair_step`) moved into
  `common::fastq` so `fastx rescue` can reuse the same mechanism with
  the opposite error-handling policy (bail-on-desync vs.
  stop-and-keep-what-you-have).
- Project description broadened from "genomics/transcriptomics" to
  "bioinformatics" throughout `Cargo.toml`, `README.md`/`README.zh.md`,
  and `main.rs` — this is a general-purpose toolkit, not limited to
  those two subfields.
- Release workflow no longer builds a glibc Linux target: it links
  against the build runner's glibc version and fails to run on older
  systems with `GLIBC_2.XX not found`. Only the fully static `musl`
  build (verified with `ldd`: "statically linked") is shipped for
  Linux now — four platforms per release, not five.

### Fixed

- Release workflow's default `GITHUB_TOKEN` lacked `contents: write`
  permission, so the first tagged release build's asset-upload step
  failed with a 403 on every platform even though all five builds
  (including the since-removed glibc target) had already succeeded.

## [0.1.0] - 2026-08-09

Initial release.

### Added

- `fastx rename` / `gff rename`: reimplements
  [`ChangeChrNameInFaOrGff.py`](https://github.com/yangqimeng99/svlearn-paper-code/blob/main/scripts/ChangeChrNameInFaOrGff.py)
  from the SVLearn paper-code repository, verified byte-identical
  against the original script on FASTA/GFF test cases, with native
  gzip/bgzip support (by content, not extension) and rayon-parallel
  chunked line processing.
- `fastx sample`: FASTQ downsampling by proportion (parallel,
  deterministic SplitMix64-hash-based decisions) or exact count
  (single-pass reservoir sampling, O(k) memory). Paired-end R1/R2 are
  sampled together in one process with structural pairing guarantees
  (matching read counts, matching IDs), instead of relying on two
  separate same-seed `seqkit sample`-style invocations.
- CI (fmt/clippy/build/test across Linux/macOS/Windows) and a
  cross-platform release workflow.
