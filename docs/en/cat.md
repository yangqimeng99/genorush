# `fastx cat`: design and internals

Source: `src/fastx/cat.rs`, `src/common/hash.rs`.

## Motivation

Concatenating raw FASTQ from multiple lanes/flowcells of the same
biological sample (topping up coverage across repeated sequencing runs) is
standard practice and normally safe — real Illumina read IDs encode
flowcell/lane/tile/coordinate information, so genuine cross-run ID
collisions aren't expected in practice. The realistic failure mode isn't
the ID scheme; it's operator error: the same file accidentally included
twice in a file list (a typo'd path, a copy-pasted glob that matched more
than intended). Plain `cat` gives that mistake zero visibility — it
silently inflates coverage and duplicates data flowing into downstream
alignment or variant calling, and there's nothing about the resulting file
that flags it as wrong. `fastx cat` does what `cat` does, but checks for
exactly this as it streams through, and aborts with the specific source
files and record positions involved the moment a duplicate ID shows up.

## Streaming duplicate-ID detection

As each source file streams through, every record's
`FastqRecord::base_id()` is hashed (`common::hash::fnv1a`, the same
function `fastx deinterleave` uses for layout detection — see
`docs/en/interleave.md` for why a hash instead of the full ID string) and
looked up in a running `HashMap<u64, (PathBuf, u64)>` keyed by that hash,
storing the source file and local record index where it was first seen. A
second occurrence of the same hash is reported immediately:

```
duplicate read ID "...": first seen in run1_R1.fq.gz (record #412),
again in run1_R1.fq.gz (record #412) -- did you accidentally include
the same file twice?
```

This fires mid-stream, not after the fact — no output past the duplicate
point has been committed to anything meaningful the caller would need to
clean up, and the specific files/positions in the message are exactly what
someone would need to go fix their file list. `--allow-duplicate-ids`
disables the check for cases where it's a false positive (a platform that
doesn't guarantee globally unique IDs).

## Paired-end mode checks two things at once

In paired-end mode (`--r2` given), each source pair is read concurrently
via `spawn_reader`/`recv_pair_step` — the same mechanism `fastx
sample`/`fastx rescue`/`fastx interleave` use — which means `fastx cat`
gets a second, independent check for free: R1/R2 pairing *within* each
source file pair. A source run where the mates have drifted out of sync
(different read counts, or IDs that don't correspond position-by-position)
is caught and reported before its data is concatenated in, not lumped in
with a downstream failure that would be much harder to trace back to
"which of the five input files was actually bad."

## Multiple sources, explicit order

Inputs are given as repeated `--r1`/`--r2` flags rather than a single
comma-separated list or a directory glob, so the concatenation order is
always exactly what's written on the command line, with no
platform-dependent glob-expansion ordering to reason about:

```
genorush fastx cat --r1 run1_R1.fq.gz --r1 run2_R1.fq.gz \
                    --r2 run1_R2.fq.gz --r2 run2_R2.fq.gz \
                    -o merged_R1.fq.gz -O merged_R2.fq.gz
```

`--r1`/`--r2` must be given the same number of times, in corresponding
order (source `i`'s R1 is `--r1`'s `i`-th occurrence, its R2 is `--r2`'s
`i`-th). Single-end mode is `--r1` only, `-O`/`--out2` omitted.

## Validated behavior

Tested with two constructed 500-pair "sequencing runs" (distinct read-ID
prefixes, as real different runs would have): concatenating both cleanly
succeeded with 1,000 pairs in the output, in source order, content and
R1/R2 pairing verified byte-identical/zero-mismatch against the expected
result. Deliberately passing the same source file twice was caught
immediately with a duplicate-ID error identifying the exact file and
record position; `--allow-duplicate-ids` was verified to suppress that
specific check while still processing the (now-doubled) input correctly.
