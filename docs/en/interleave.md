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

## Detecting the layout

The core difficulty: testing the "concatenated" hypothesis requires knowing
the *total record count* before you can even state it. Record `i` is
hypothesized to match record `i + n/2` — but you don't know `n/2` until
you've seen the whole file. There's no way around reading the entire input
before committing to an answer, so `--layout auto` (the default) does
exactly that:

1. Stream the whole file once via `read_fastq_record`, computing
   `common::hash::fnv1a(base_id())` for every record — an 8-byte hash, not
   the full ID string, so memory is `8 * n` bytes regardless of how long
   IDs are (a few GB even for a several-hundred-million-record file).
2. With every hash collected and the total count `n` known, test both
   hypotheses against the array in memory (cheap, `O(n)` comparisons, no
   further I/O):
   - interleaved: `hash[2*i] == hash[2*i+1]` for all `i`
   - concatenated: `hash[i] == hash[i + n/2]` for all `i`
3. Exactly one hypothesis holding cleanly is a confident answer. Both
   holding (a degenerate case, e.g. a tiny file) or neither holding is
   reported as a hard error with the exact first mismatch position for
   each hypothesis — this command refuses to guess. Guessing wrong here
   means silently corrupting every downstream analysis with fabricated
   pairs; an error asking for `--layout` to be passed explicitly is strictly
   preferable.
4. A second pass then performs the actual split using the confirmed layout.

This means `--layout auto` reads the input twice — once to hash and
decide, once to split — which is real, unavoidable cost given the
ambiguity, not a design shortcut. `--layout interleaved`/`--layout concat`
skip detection for a known layout. `interleaved` needs only a single
streaming pass with no extra memory (mate parity is determined purely by
position, checked as you go). `concat` still needs the midpoint up front,
so it does a first pass that only counts records (skipping the hash
computation `--layout auto` would do) — cheaper than full detection, but
still a real pre-pass, because the midpoint is unavoidably a function of
the whole file's length.

## Why hashes, not full ID strings

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

Both `split_interleaved` and `split_concat` process the input in
`--chunk-records`-sized batches (parallel-compressed via
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
