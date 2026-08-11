# `fastx rename` / `gff rename`: design and internals

Source: `src/common/rename.rs`, `src/fastx/rename.rs`, `src/fastx/mod.rs`, `src/gff/rename.rs`, `src/gff/mod.rs`, `src/io_utils.rs`.

## Origin

This command reimplements [`ChangeChrNameInFaOrGff.py`](https://github.com/yangqimeng99/svlearn-paper-code/blob/main/scripts/ChangeChrNameInFaOrGff.py)
from the SVLearn paper-code repository: given a two-column mapping table
(`new_name  old_name`), rewrite chromosome/contig/sequence names in a FASTA
or GFF file. It's a small tool, but a common first step in almost every
comparative-genomics or SV pipeline (aligning accession-style names like
`NC_019458.2` to short display names like `1`, `2`, `X`), so it's the
reference implementation for this project's overall architecture.

## CLI shape

```
genorush fastx rename <FILE> -n <NAME_MAP> -o <OUTPUT>
genorush gff   rename <FILE> -n <NAME_MAP> -o <OUTPUT>
```

The original Python script took a single `--fa`/`--gff` flag pair to select
the format. Here the format is encoded in the *category* (`fastx` vs `gff`)
instead: two thin leaf commands (`fastx::rename`, `gff::rename`) each supply
their own line-transform closure to one shared engine
(`common::rename::run`). This is the template every future `<category>
<action>` command in this project follows — see
[`../en/architecture.md`](../en/architecture.md) if present, or just read
`src/fastx/mod.rs` / `src/gff/mod.rs` for the pattern: a category module
owns a `clap::Subcommand` enum and a `run()` dispatcher; each leaf command
lives in its own file and calls into `common::` for anything reusable.

## Behavioral contract (ported from the Python script, verified byte-for-byte)

Mapping file (`-n/--name`): whitespace-separated, two columns per line,
`new_name old_name`. Loaded fully into a `HashMap<old_name, new_name>`.

FASTA transform: for a header line (`>...`), take the first
whitespace-delimited token, strip the leading `>`, look it up in the map.
If found, write `>{new_name}` — **note the sequence description after the
first token is dropped**, matching the original script's
`line.split()[0][1:]` exactly (a debatable design choice, but this is
byte-for-byte parity with the tool this replaces, not a novel design).
If not found, write `>{old_name}` (also dropping the description). Non-header
lines pass through unchanged (after trimming, see below).

GFF transform: comment lines (`#...`) pass through unchanged. Data lines are
split on the first tab; if column 1 is in the map, the line is rewritten as
`{new_name}\t{rest}`. This is always `{new}\t{rest}` — even when `rest` is
empty (a malformed/short line) — because that's what Python's
`'\t'.join(LineList[1:])` produces when `LineList` has only one element: an
empty string, joined after a tab that's still there. `src/gff/rename.rs`
keeps this quirk deliberately, with a comment explaining why, so a future
maintainer doesn't "fix" it and silently break parity.

Every line is trimmed with `.trim()` (both ends) before processing, mirroring
Python's `line.strip()`. This means leading whitespace on a sequence line,
or trailing `\r` from a CRLF file, is stripped on the way through — again,
inherited behavior, not a new decision.

## Where this implementation deliberately diverges

The original script has one silent bug and one platform-portability problem
that are not worth reproducing:

1. **Neither `--fa` nor `--gff` passed**: the Python script loads the input
   and then does nothing (neither `if` branch fires), silently producing an
   empty output file. `common::rename` doesn't have this ambiguity at all —
   the *category* (`fastx` vs `gff`) is the format selector, so there is no
   state where the format is unspecified.

2. **Gzip via `less`**: the original does
   `os.popen(f'less {input}').readlines()`. This only decompresses `.gz`
   input if the machine's `less` is wired up with `lesspipe`, is a shell
   command built from an unescaped filename (injection risk on adversarial
   paths), and simply doesn't exist as a strategy on Windows. `io_utils::open_reader`
   instead sniffs the first two bytes for the gzip magic number (`1f 8b`)
   and, if present, wraps the file in `flate2::read::MultiGzDecoder` —
   `Multi`, not the plain `GzDecoder`, because bgzip-compressed reference
   genomes are valid *concatenated multi-member* gzip streams, and a
   single-member decoder would silently truncate after the first block.
   Detection is by content, not by file extension, so a `.gz`-less but
   actually-gzipped stream still works. Output is gzip-compressed
   automatically when the output path ends in `.gz` — a feature the
   original didn't have at all (`click.File('w')` only ever writes plain
   text).

## Parallel processing model

Sequential per-line transforms are embarrassingly parallel here: each output
line depends only on its own input line and the (read-only, shared) name
map — never on neighboring lines or on any running state. `common::rename::run`
exploits this directly:

1. Read up to `--chunk-lines` (default 200,000) lines into a `Vec<String>`.
2. Hand the chunk to `rayon`'s `par_iter().map(transform).collect()` —
   `rayon::collect` on a `Vec` preserves input order regardless of how work
   is split across threads, so no explicit re-ordering step is needed.
3. Write the transformed chunk, repeat until EOF.

Chunking exists for one reason: **bounding memory**. A whole-genome FASTA
can be tens of millions of lines; without chunking, either the whole file
sits in memory as `Vec<String>`, or you pay the complexity of a fully
streaming parallel iterator. Chunking is the simple middle ground — peak
memory is `O(chunk_lines)`, independent of file size, and large enough
chunks keep rayon's per-task overhead negligible relative to per-line work.

`-j/--threads` (global flag, `src/main.rs`) sets rayon's global thread pool
size once at startup; `0` (the default) means "use all logical cores",
which is rayon's own default behavior.

### Parallel gzip output too

The above parallelizes the *transform* — computing each output line's
content. Writing that output, when the destination ends in `.gz`, is a
separate concern with its own single-threaded bottleneck if left
unaddressed: this was discovered the hard way in `fastx sample` (see
`docs/en/sample.md`'s `BlockWriter` section for the full story — a
57+59 GB real-world dataset where `-j` measurably did nothing, because
compression ran on one thread no matter the thread count). `common::rename::run`
was migrated to `io_utils::BlockWriter` for the same reason: each
transformed chunk's lines are split into `rayon::current_num_threads()`
groups (`lines_into_blocks`, `common/rename.rs`) and compressed into
independent gzip members in parallel, rather than streamed through one
long-lived `GzEncoder`. For plain-text output this changes nothing (no
compression to parallelize); for gzip output on a large reference genome,
`-j` now does real work on the write side, not just the transform side.

## Measured results

Validated against the real Python script (not a re-derivation of its logic —
the actual script, fetched and executed) on FASTA with header descriptions,
GFF with comment lines and unmapped contigs, and gzip input/output: output
is byte-identical (`diff` clean) in every case tested. On a 112 MB / 1.9M-line
synthetic FASTA (29 chromosomes), wall time was ~5.9 s for the Python script
vs. ~0.65–0.79 s here — the gap is dominated by Python's `os.popen`/`less`
subprocess overhead and interpreter line-processing cost, not by the
rename logic itself; this is I/O-bound work, so the multi-threaded chunking
buys comparatively little here (1 thread vs. 12 threads gave similar wall
time on that file). The parallel architecture starts paying off more
directly on commands with heavier per-record computation.

## Extending this pattern

To add a new `<category> rename`-like command:

1. If the transform is line-oriented and stateless like this one, write it
   as a `Fn(&str, &HashMap<String, String>) -> String` (or generalize
   `common::rename::run`'s signature further if the shared state isn't a
   name map) and reuse `common::rename::run` directly — see
   `src/fastx/rename.rs` and `src/gff/rename.rs` for the ~15-line pattern.
2. If the transform needs record-level (not line-level) structure — e.g.
   anything touching whole FASTA/FASTQ records — see `docs/en/sample.md`
   instead: `common::fastq` models 4-line FASTQ records, and the same
   chunk-then-`rayon::par_iter` strategy applies at the record granularity.
