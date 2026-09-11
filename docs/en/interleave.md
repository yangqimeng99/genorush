# `fastx interleave` / `fastx deinterleave`: design and internals

Source: `src/fastx/interleave.rs`, `src/fastx/deinterleave.rs`, `src/common/hash.rs`.

## Motivation

Two directions, two different problems. Merging R1/R2 into one interleaved
file (`fastx interleave`) is unambiguous — the output layout is ours to
choose, so it's always proper `R1,R2,R1,R2,...` interleaving, the
convention tools like `bwa mem -p` expect. Splitting a merged file back
apart (`fastx deinterleave`) is the hard direction, because "merged" isn't
one format in practice. Proper interleaved files alternate mates. But a lot
of files in circulation are just `cat R1.fastq R2.fastq > merged.fastq` —
every R1 record, then every R2 record, back to back. These are not the
same format wearing different clothes: a splitter that assumes interleaved
and gets a concatenated file silently produces garbage, pairing up two
unrelated reads as if they were mates, with no error at all. This surfaced
directly from real usage — see `docs/en/sample.md` and `docs/en/rescue.md`
for the same lesson learned twice already (pairing correctness has to be
verified, never assumed) — and `fastx deinterleave` applies it to a third
place it can go wrong.

## The third layout: when position means nothing

Both layouts above are *positional*: which mate a record belongs to is a
function of where it sits in the file. Real files exist where that premise
is simply false, and one turned up while validating this command against
the shell pipeline it replaces.

`SRR17458599.fastq.gz` (14.76 GB, 399,065,686 records) is laid out as:

```
records          1 …  24,174,089   all /2   (24,174,089 records)
records 24,174,090 … 223,706,932   all /1   (199,532,843 records)
records 223,706,933 … 399,065,686  all /2   (175,358,754 records)
```

Three runs, not two — so the concatenation hypothesis fails, and the runs
are long, so the interleaving hypothesis fails at the very first pair. Worse
for any id-based reasoning, the read numbering is *globally sequential*
across the whole file (verified: not one gap in 399 million records), which
means the two mates of a pair carry **different** ids: `@SRR17458599.1 1/2`
is mated to `@SRR17458599.24174090 24174090/1`. `base_id()` comparison
cannot pair this file up, and no positional rule can split it.

What does work is the one piece of information every record carries about
itself: the `/1` or `/2` marker in its header. That is what `--layout
by-suffix` routes on, and it is what the `awk` one-liners in circulation
have always done —

```awk
if (header ~ /\/1(\s|$)/)      print ... >> r1
else if (header ~ /\/2(\s|$)/) print ... >> r2
```

— which is why those pipelines handle files like this one when a positional
splitter cannot.

## Detecting the layout: probe the head, verify while writing

Testing the "concatenated" hypothesis requires knowing the *total record
count* before you can even state it: record `i` is hypothesized to match
record `i + n/2`, and `n/2` isn't known until the whole file has been read.
Earlier versions took that literally and always read the input twice — once
to hash every `base_id()` into memory and decide, once to split. On the file
above that cost 5 minutes and **3.05 GiB** of resident memory (8 bytes ×
399 M records) only to conclude that neither hypothesis held.

`--layout auto` now probes instead. It reads the first `PROBE_RECORDS`
(4096) records and asks two questions:

1. **Do adjacent records share a base id?** If every `(2i, 2i+1)` pair in
   the probe does, the file looks interleaved.
2. **Does every record carry a `/1` or `/2` marker?** If so — and the file
   is not interleaved — each record can be routed by its own header, no
   matter how the file is ordered: `by-suffix`.

Only when neither holds (no markers *and* no shared ids) does it fall back
to the full hash-everything scan, which remains the only way to recognise a
concatenation whose mates share ids and carry no markers. That is now the
sole code path whose memory grows with the input.

A probe is evidence, not proof, so the chosen splitter re-checks its own
assumption on every record as it writes: `split_interleaved` requires the
two records of each pair to share a base id, `split_by_suffix` requires
every record to carry a marker, and `split_concat` requires the marked
records in each half to agree with each other and differ across the
midpoint. Any violation stops the run with the offending record number and a
pointer to the mode that would handle it.

The trade-off is explicit: a violation is now caught *after* some output has
been written rather than before anything is written. The failing run exits
non-zero and says the outputs must be discarded. In exchange, the common
case stops reading multi-gigabyte inputs twice. Refusing to guess is still
the rule — what changed is that the refusal costs one probe instead of one
full pass, and that a file whose records identify themselves is no longer
treated as unguessable in the first place.

`--layout interleaved`/`--layout concat`/`--layout by-suffix` skip the probe
entirely. `interleaved` and `by-suffix` are single streaming passes with no
extra memory. `concat` still needs the midpoint up front, so it does a
first pass that only counts records — cheaper than the hash scan, but still
a real pre-pass, because the midpoint is unavoidably a function of the whole
file's length.

## Probing a stream you can't rewind

The probe reads records, and those records are still part of the input. An
earlier version simply re-opened the file and started over, which is free
enough on a file and impossible on a pipe.

Instead, the probed records are buffered (4096 of them, a few MB at most)
and handed back to the front of the stream: `Records` yields its replay
buffer before touching the channel again, so the splitter sees byte-for-byte
the same sequence of records it would have seen from a fresh read. This
makes `--layout auto` work on stdin, and removes the second open in the file
case as a side effect.

What the replay can't rescue is a layout that needs the *whole* input before
it can start:

- `--layout concat` needs the record count to locate the midpoint.
- `--layout auto` falls back to the hash-everything scan when the head is
  inconclusive, then splits.

Both read the input twice, so both are refused on stdin with a message
naming the single-pass alternatives (`interleaved`, `by-suffix`) rather than
consuming half the stream and failing somewhere less legible.

## Routing by header marker (`--layout by-suffix`)

`FastqRecord::mate_suffix()` scans the *entire* header for `/1` or `/2`
followed by whitespace or end-of-header, deliberately mirroring the `awk`
convention above. Scanning the whole header (rather than just the first
token, as `base_id()` does) is what makes the SRA layout work, where the
marker sits in a second field: `@SRR17458599.1 1/2`.

Two conventions are deliberately *not* treated as markers:

- Modern Illumina headers (`@ID 1:N:0:ATCG`). The mate number there is a
  field of a different grammar, and matching it would misread any header
  that happens to contain a similar-looking field. Files with those headers
  are interleaved or concatenated in practice, and the positional modes
  handle them.
- A marker glued to more text (`@READ/12`, `@READ/1x`). The marker must end
  the header or be followed by whitespace.

Where this mode is deliberately stricter than the `awk` it reproduces: a
record with no marker at all is a hard error, not a warning on stderr. The
shell version prints a warning and drops the record, which means a silent
loss of data in the middle of a multi-hour job that still exits 0. In a mode
whose entire premise is that position carries no information, there is
nothing to fall back on for that record.

Finally, `by-suffix` cannot guarantee what a positional split gets for free:
that R1 and R2 come out the same length. So the record counts are compared
at the end, and a mismatch is an error — the outputs are complete and
written, but they cannot be treated as positionally paired downstream, and
the run says so and exits non-zero.

## Why hashes, not full ID strings

This applies to the fallback detection path described above — the only one
that still holds per-record state.

`common::hash::fnv1a` is a small, dependency-free, non-cryptographic hash.
It is not collision-proof, and this is a deliberate, documented trade-off:
for real read IDs (structured, effectively always distinct across
unrelated reads), the odds of an accidental 64-bit collision are
astronomically small next to the actual failure modes this machinery
exists to catch — a file that isn't cleanly one layout or the other.
Storing full ID strings for hundreds of millions of records would cost
tens of gigabytes for no correctness benefit worth that cost; the same
hash is reused by `fastx cat` for its duplicate-ID check (see
`docs/en/cat.md`), on the same reasoning.

## Splitting without pre-aligning to chunk or midpoint boundaries

All three splitters (`split_interleaved`, `split_concat`,
`split_by_suffix`) process the input in `--chunk-records`-sized batches (parallel-compressed via
`common::fastq::format_into_blocks` / `io_utils::BlockWriter`, the same
machinery `fastx sample` uses — see `docs/en/sample.md`), but track a
*global* running index rather than resetting per chunk, so a chunk that
happens to straddle an odd/even boundary (interleaved) or the `n/2`
midpoint (concat) is still routed correctly record-by-record within that
chunk. Neither split function needs chunk boundaries to align with
anything meaningful in the data.

## `fastx interleave`: no detection needed, order survives arbitrary block splits

Since the output format is simply "whatever we choose to write," each
buffered chunk of pairs is flattened into `[r1, r2, r1, r2, ...]` before
being handed to `format_into_blocks`, which splits it into
`rayon::current_num_threads()` independently gzip-compressed blocks. This
works correctly even when a block boundary falls between an `r1` and its
`r2` — gzip decompression just concatenates all blocks back into one
continuous byte stream in order, and that stream is exactly the flattened
list regardless of where it was cut for parallel compression. Pairing
correctness during the *read* side reuses `spawn_reader`/`recv_pair_step`
and the same `--no-pair-check` escape hatch as `fastx sample`/`fastx rescue`.

## Validated behavior

Tested against three constructed 500-pair inputs: a genuine interleaved
file, a `cat R1 R2`-style concatenation, and a fully shuffled file that is
neither. `--layout auto` correctly identified the first two and refused to
guess on the third, reporting a specific first-mismatch index for both
hypotheses. Content from both correctly-detected cases was verified
byte-identical to the original R1/R2 after splitting. `--layout
interleaved`/`--layout concat` (explicit, detection skipped) were verified
against the same fixtures. A full round trip — `interleave` then
`deinterleave` — was verified to reproduce the original R1/R2 exactly.

Unit tests cover each splitter's routing and each guard directly: an
interleaved split, a rejected non-interleaved pair, `--no-pair-check`
overriding that rejection, a concat split (including one whose midpoint
falls inside a chunk), a rejected concat whose mate runs miss the midpoint,
`by-suffix` routing an out-of-order file and preserving input order within
each mate, a rejected unmarked record, and the four head-probe verdicts.

`by-suffix` was then validated against the pipeline it replaces, on the
real `SRR17458599` data described above:

- A 400,000-record slice was split by both the original `awk` and by
  `genorush` (once via `--layout auto`, once via explicit `--layout
  by-suffix`). All three R1 outputs are byte-identical (`cmp`, and md5
  `881b2b54…`), as are all three R2 outputs (md5 `486db68d…`).
- On that same slice, `--layout concat` splits at the midpoint and warns
  that `--out1` holds the R2 reads, since the file's first half is R2 — the
  silent mate swap that motivated the check.
- The full 14.76 GB file was then split end to end with `--layout auto`
  (which probed the head and chose `by-suffix`), producing 199,532,843
  records on each side. Decompressed, both outputs are md5-identical to the
  `awk` pipeline's original run: R1 `01b115c750218ce9ff7746fec545d0bd`,
  R2 `9a0643deb542c83d453176d6b419e840`.
- That run took 12m32s wall clock at `-j 12` with a peak RSS of **45 MiB**.
  For contrast, the same input under the old always-hash detection spent 5
  minutes and 3.05 GiB just to reach its "cannot determine the layout"
  error, before splitting anything.
