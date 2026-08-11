# `fastx rescue`: design and internals

Source: `src/fastx/rescue.rs`, shared infrastructure in `src/common/fastq.rs`
(`spawn_reader`, `PairStep`, `recv_pair_step`).

## Motivation

A very common real-world failure mode: a FASTQ download (single-end long
reads, or an Illumina R1/R2 pair) gets interrupted partway through, leaving
a `.fq.gz` whose gzip stream is truncated or has a bad trailing checksum.
The file isn't garbage — everything decoded *before* the truncation point
is perfectly good sequencing data, and for paired-end reads, every pair
where both mates were fully and correctly decoded is still usable for
alignment. Throwing the whole file away discards real, expensive-to-
regenerate data for no reason. `fastx rescue` recovers exactly the leading
run of clean, well-formed records (or read pairs) and writes them out,
stopping cleanly at the first sign of trouble instead of erroring out.

## Reusing `fastx sample`'s reading infrastructure, with different error semantics

`fastx sample` already had to solve "read two mate files concurrently and
detect when they go out of sync" (see `docs/en/sample.md`). `rescue` needs
the exact same mechanics — but where `sample` treats a desync as fatal
(abort, because sampling from misaligned pairs would silently corrupt the
output), `rescue`'s entire purpose is to keep going as far as possible and
then stop gracefully at that exact point. Rather than duplicate the
reading/pairing logic with inverted error handling, the shared primitives
were pulled out into `common::fastq`:

- **`spawn_reader(path) -> Receiver<Result<FastqRecord>>`**: decompresses
  and parses one file on its own OS thread, streaming records through a
  bounded channel. On a read/parse error, it sends that one `Err` and stops
  — no retry, no attempt to resync on later bytes (there's no reliable way
  to find the next record boundary in a stream that just proved untrustworthy).
  In paired-end mode, running one of these per mate file is what makes R1
  and R2 decompress concurrently rather than sequentially — the same I/O
  win described in `docs/en/sample.md`.

- **`PairStep` / `recv_pair_step(rx1, rx2)`**: pulls one record from each
  mate's channel and classifies the outcome into one of four cases —
  `Pair { r1, r2, ids_match }` (both sides produced a record; `ids_match`
  reports whether `FastqRecord::base_id()` agreed), `Eof` (both channels
  closed at the same step — a clean, synchronized end), `CountMismatch`
  (one side closed while the other still had data), or `ReadError` (either
  side hit a genuine read/parse failure). Deliberately, this function makes
  no judgment about which of these are "bad" — that's left entirely to the
  caller. `sample::recv_pair` (in `src/fastx/sample.rs`) turns
  `CountMismatch`, `ReadError`, and a false `ids_match` into hard errors via
  `bail!`. `rescue::run_pe` turns the exact same three outcomes into "stop
  the loop, log a warning, keep everything already written."

This is the same lesson `fastx rename`'s doc draws from a different angle:
when two commands need the same mechanism but different policies on top of
it, separate the *mechanism* (what happened) from the *policy* (what to do
about it) at the function-return level, rather than parameterizing one
function with a `bail_on_error: bool` flag that would need to grow a new
parameter for every future caller's variant policy.

## What counts as "the failure point"

Three distinct things can end a rescue run, and the log message
distinguishes them so a user (or a script parsing the log) knows what
actually happened:

1. **A read/parse error** — the gzip stream is corrupt or truncated, or a
   record is structurally broken (missing a line, a `+`-line that doesn't
   start with `+`, or a sequence/quality length mismatch). Surfaces via
   `PairStep::ReadError` / a direct channel `Err` in single-end mode.
2. **A read-count mismatch** (paired-end only) — one mate's stream ended
   (cleanly or not) before the other's. Even if the shorter file ended
   "cleanly," there's no valid pair to keep once one side runs out.
3. **An ID mismatch** (paired-end only, default on) — both mates still
   produced a record, but `FastqRecord::base_id()` disagrees. This usually
   *is* the corruption: garbage bytes downstream of a truncation point can
   still parse as syntactically valid FASTQ-looking lines by pure chance,
   and a base-ID check catches that even when the 4-line structure itself
   looks fine. `--no-pair-check` disables only this check (a
   `CountMismatch` or `ReadError` still stops the run) — use it only if
   your headers don't follow the `/1`+`/2` or Illumina `1:...`+`2:...`
   mate-suffix conventions `base_id()` relies on.

A clean run — reaching `PairStep::Eof` / a clean channel close on both
sides without ever hitting one of the above — is reported as "no
corruption detected," not as a rescue; the command works as a strict
pass-through validator in that case, at the cost of no more than reading
the file end to end once.

## Exit codes: making a partial rescue distinguishable from full success or total failure

A rescue command has three qualitatively different outcomes a caller might
want to branch on in a script, so it uses three exit codes instead of the
usual 0/1:

- **0** — read cleanly to the end; nothing was corrupted, nothing was lost.
- **3** — corruption was detected, but at least one record (or pair) was
  successfully rescued and written.
- **1** — corruption was detected and *nothing* was salvageable (e.g. the
  very first record fails), or the arguments/files were invalid to begin
  with. Nothing useful was produced.

This makes `if genorush fastx rescue ...; then ... elif [ $? -eq 3 ]; then
echo "partial, check logs"; else echo "unusable"; fi` a meaningful pattern
in a pipeline script.

### A subtlety: `std::process::exit` skips destructors (and how `BlockWriter` sidesteps it)

`finish()` in `rescue.rs` calls `std::process::exit(code)` directly for the
1/3 cases, since Rust's normal `main() -> Result<()>` exit-code convention
only distinguishes 0 from 1. `std::process::exit` terminates the process
immediately *without running destructors* — and the first version of this
command hit exactly the failure mode that matters here: it wrote through
`io_utils::open_writer`, whose gzip path was a single long-lived
`flate2::write::GzEncoder` that writes its trailer (CRC32 + uncompressed
size) in `Drop`, not on `flush()`. Calling `process::exit` while that
`GzEncoder` was still alive produced a `.gz` file with a valid header and
body but a missing trailer — rejected as truncated by every gzip reader,
even though the rescued *content* was byte-for-byte correct. The fix at
the time was an explicit `drop(writer)` before `finish()`, forcing the
trailer to be written while normal scope-based destruction was still in effect.

That workaround is gone now, not because the underlying hazard was patched
around again, but because it doesn't apply to the type in use anymore.
`rescue` was migrated to `io_utils::BlockWriter` (see `docs/en/sample.md`
for why it exists) to get parallel-compressed output the same as every
other chunk-processing command — confirmed-good records are buffered in
`--chunk-records`-sized batches and handed to `write_blocks`, flushing
whatever partial batch is left over the moment corruption is hit or the
input ends. `BlockWriter`'s gzip path creates a *fresh* `GzEncoder`
per block and calls `.finish()` on it synchronously, inside the same
`write_blocks` call that produced it — so every gzip member is fully
complete, trailer included, the instant `write_blocks` returns. There is
no persistent, Drop-dependent encoder state to lose to `process::exit`
anymore; `writer.flush()` (still called before `finish()`) only has to push
the underlying `BufWriter`'s buffered bytes to the OS, which is a
straightforward, non-hazardous flush. Worth remembering the general lesson
anyway — mixing a possibly-early-exiting code path with anything that
defers meaningful work to `Drop` is a trap — even though this specific
instance of it happens to no longer exist in this codebase.

## Limitations

This only detects *structural* corruption — a violation flate2 or
`read_fastq_record`'s own validation can actually observe (bad gzip data,
a truncated record, a header/plus-line/seq-qual-length violation). FASTQ
carries no per-record checksum, so bit-level corruption that leaves a
record's 4-line structure and lengths intact but scrambles its actual
sequence/quality content is invisible to this command, exactly as it would
be to any other FASTQ tool.

## Validated behavior

Tested against the same 200,000-pair / 20,000-long-read synthetic datasets
used for `fastx sample`, deliberately truncated: a single-end `.fq.gz` cut
to 60% of its size recovers a rescued file that `gzip -t` reports as fully
valid, and whose content is byte-for-byte identical to the corresponding
prefix of the original uncompressed file. A paired-end case with R1 intact
and R2 truncated to 70% recovers matching R1/R2 pair counts with zero ID
mismatches across the full rescued output. A plain-text (non-gzip) file
truncated mid-record is caught by the seq/qual-length structural check
rather than a decompression error, proving the detection isn't gzip-
specific. A fully clean file round-trips with exit code 0 and zero data
loss. A file corrupted from its very first record correctly exits 1 with
zero reads rescued. `--no-pair-check` was verified to suppress only the ID
check (an artificially ID-shifted but count-matched pair of files reads
through to completion under it) while still enforcing the count-mismatch
check.
