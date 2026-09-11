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
looked up in a running `HashSet<u64>` of everything seen so far. A second
occurrence of the same hash is reported immediately:

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

## What the set holds, and why it holds so little

This set is the command's entire memory footprint, and it is the one thing
here that grows with the input: one entry per read, kept until the run ends.
At the scale this tool is aimed at, what goes in each entry is not a detail.

An earlier version stored, per hash, the source path and record index where
it was first seen, so the error above could name both ends of a collision.
That is `HashMap<u64, (PathBuf, u64)>` — and a `path.to_path_buf()` heap
allocation *per read*, every one of which stays in the table. Measured on
2,000,000 reads, the check added 294 MiB over the same run with
`--allow-duplicate-ids`: **154 bytes per read**, once the table's 40-byte
entries, those allocations, and the transient doubling during a resize are
all counted. A 30x bovine WGS sample is around 270 million pairs, which puts
that at roughly 42 GB — on a shared machine, a job that dies hours in, after
writing most of its output.

Keeping only the hashes brings it to **28 bytes per read** measured the same
way (about 7.4 GB extrapolated), and made the check *faster*, not slower:

| | wall clock | peak RSS | check's own cost |
|---|---|---|---|
| `--allow-duplicate-ids` | 0.84 s | 42 MiB | — |
| storing paths and positions | 1.05 s | 336 MiB | 154 B/read, +25% time |
| hashes only | 0.92 s | 95 MiB | 28 B/read, +10% time |

Three things account for the speedup: the per-read allocation is gone, a
table a third the size keeps far more of itself in cache, and the keys are
hashed once instead of twice — `common::hash::BuildIdHasher` takes the FNV
value as given rather than running SipHash over it again, applying only
splitmix64's finalizer so the table still sees well-distributed bits (FNV's
low bits, which decide the bucket, are its weakest).

The 28 bytes are mostly not the entries themselves. A `HashSet<u64>` entry
is 8 bytes plus a 1-byte control tag, but capacity is a power of two held
below 7/8 full, and a resize briefly holds the old table alongside the new
one — the peak lands at roughly twice the steady state.

## Recovering the first occurrence

Storing only hashes costs the error message its "first seen in ... (record
#N)" half, so that half is recovered on demand: when a repeated hash turns
up, `locate_first` re-reads the inputs to find where the ID actually first
appeared. That pass is only ever taken once a repeat has been found, which
on real data means the command is about to stop anyway.

It also makes the check *exact*. FNV-1a is 64 bits and not
collision-proof, and the old version would abort on a collision between two
genuinely different IDs. The rescan distinguishes the two cases: if the only
record carrying that ID is the one in hand, nothing was repeated, and the run
continues (logged at debug level). Where the rescan is impossible — one of
the inputs is stdin, which cannot be read twice — the command reports the
duplicate without the earlier position and says why.

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
