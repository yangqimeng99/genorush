# `fastx sample`: design and internals

Source: `src/fastx/sample.rs`, `src/common/fastq.rs`, `src/common/rng.rs`.

## Motivation

`seqkit sample` is the reference tool this command is modeled on and
deliberately improves on one specific gap: **it has no concept of
paired-end input**. To subsample an Illumina R1/R2 pair, you run it twice —
once per mate file — and rely on passing the *identical* `-s/--rand-seed`
both times so that two entirely independent invocations happen to make the
same keep/discard decision at every read index. `genorush fastx sample`
instead reads both mates in a single process and makes one decision per
read *pair*, so correct pairing is a property of the algorithm, not a
convention the caller has to uphold correctly on every invocation.

## Reading `seqkit sample`'s actual algorithm

Before designing around it, its source was read directly
(`shenwei356/seqkit`, `seqkit/cmd/sample.go`) rather than assumed from the
CLI docs. Two modes:

- **`-p` (proportion)**: single pass; for every record, draw
  `rand.Float64()` from a seeded `math/rand.Rand` and keep the record if the
  draw is `<= p`. One RNG draw per record, in file order, from one
  sequentially-advancing generator.
- **`-n` (number)**: either loads *all* records into memory
  (`fastx.GetSeqs`) and applies the same per-record Bernoulli trial with
  `proportion = n/len(records) * 1.1`, stopping once `n` are kept, or —
  with `--two-pass` — reads the file once just to count records
  (`fastx.GetSeqNumber`), then reads it again applying the same inflated-
  proportion trick. Either way it's an *approximate*-then-truncated scheme,
  not an exact single-pass sampling algorithm; `--two-pass` exists purely
  to avoid the full in-memory load, at the cost of reading the file twice.

Because `-p` mode consumes exactly one RNG draw per record in strict file
order, running it twice with the same seed on two files that have the same
number of records in the same order *does* produce consistent pairs — the
two RNG streams are byte-for-byte identical, draw for draw, so the i-th
decision matches in both runs, including where the loop breaks early in
`-n` mode (since both runs hit the target count at the same read index).
The scheme is correct under those preconditions; it's the *operational*
gap — two invocations, two full file reads, no cross-file validation, a
seed that has to be typed twice correctly, and total silence if a `-n`
seed happens to differ or record counts diverge — that this command closes.

## Two algorithmic changes, not just an API change

### 1. Proportion sampling: a stateless, parallel decision function

Go's `math/rand.Rand` is inherently sequential — call N advances the
generator to call N+1. `common::rng::deterministic_f64(seed, index)` is a
pure function instead: it runs the seed and the record's global index
through the SplitMix64 mixing function
(`z ^= z>>30; z *= C1; z ^= z>>27; z *= C2; z ^= z>>31`, the same finalizer
used inside `java.util.SplittableRandom` and Rust's own `rand` crate
internals) and returns a value that depends on `(seed, index)` alone —
computing draw #5,000,000 doesn't require having computed draws
#1..4,999,999 first. That's what makes it safe to batch records into
chunks and hand each chunk to `rayon::par_iter()`: every record's
keep/discard decision is independent of which thread computes it or in
what order, and the result is identical no matter how `--chunk-records` is
tuned or how many threads are running — the property `deterministic_f64`
is even unit-tested for directly (`common::rng::tests`). In paired-end mode
this composes for free: draw index `i` decides pair `i` as a whole (one
draw, not two), so there's no possibility of the mates disagreeing.

### 2. Exact-count sampling: single-pass reservoir sampling (Algorithm R)

`seqkit sample -n` either holds the whole input in memory or reads the
file twice. Neither is necessary: **reservoir sampling** solves "pick
exactly k items uniformly at random from a stream of unknown length" in a
single pass with `O(k)` memory, full stop — no need to know the total
count in advance, no second read. The classic algorithm (Vitter's
"Algorithm R"), implemented directly in `run_reservoir_se`/
`run_reservoir_pe`:

```text
for i, item in enumerate(stream):
    if i < k:
        reservoir[i] = item
    else:
        j = uniform_random(0, i)   # inclusive
        if j < k:
            reservoir[j] = item
```

Every item ends up in the reservoir with probability exactly `k/n` once the
full stream of length `n` has been consumed — a standard, well-known proof
by induction (the tool didn't invent the algorithm, just applied it in the
place `seqkit` doesn't). `common::rng::SplitMix64::next_below(bound)`
supplies the uniform-in-`[0, bound)` draw via Lemire's nearly-divisionless
unbiased method, avoiding the small modulo bias a naive `rng() % bound`
would introduce.

In paired-end mode the reservoir holds `(original_index, R1_record,
R2_record)` tuples — a pair is the atomic unit that gets kept or evicted
together, so exactness and pairing correctness both fall out of the same
data structure. Each kept item also remembers its original stream
position, so after the pass the reservoir is sorted back into input order
before writing (`chosen.sort_by_key(|(idx, ..)| *idx)`) — reservoir
sampling's replacement order is not the original order, and re-imposing it
is just a cheap `O(k log k)` sort, done purely for output readability/
determinism, not for correctness.

This means `-n` never needs a `--two-pass` flag at all: Algorithm R
*already* bounds memory to the sample size, in one pass, regardless of
input size. (Contrast with `seqkit`, where `--two-pass` is a real trade-off
the user has to opt into.)

## Concurrent mate reading

`spawn_reader()` puts each input file's decompression + FASTQ parsing on
its own OS thread, streaming parsed `FastqRecord`s to the main thread over
a bounded `std::sync::mpsc::sync_channel`. In single-end mode this mainly
overlaps I/O with the (cheap) sampling decision. In paired-end mode it's
the difference between decompressing R1 then R2 sequentially (as two
`seqkit sample` invocations would, back to back) and decompressing both
concurrently on separate threads — gzip decompression is CPU-bound, so
this is a genuine wall-clock win on multi-core machines, independent of
whatever thread count `-j` gives to rayon for the proportion-sampling
batches.

## Pairing correctness, made structural

`recv_pair()` (`src/fastx/sample.rs`) pulls one record from each mate's
channel per step and enforces two invariants before any sampling decision
happens:

1. **Synchronized EOF.** If one channel closes before the other, that's an
   error (`read 1/2 have different numbers of reads`) — `seqkit`'s
   twice-invocation approach has no way to detect this at all; each run
   only ever sees its own file.
2. **ID correspondence** (unless `--no-pair-check`). `FastqRecord::base_id()`
   strips a `@` prefix, keeps the first whitespace-delimited token, and
   strips a trailing `/1` or `/2` — covering both legacy Illumina headers
   (`@READ/1`, `@READ/2`) and modern ones (`@READ 1:N:0:...`,
   `@READ 2:N:0:...`). A mismatch at any position aborts immediately with
   the pair index and both IDs, rather than silently sampling a garbage
   pairing.

Both checks were exercised directly: a truncated R2 (one record short) is
caught at the exact pair index where the streams diverge; an artificially
shifted R2 (headers no longer aligned) is caught at pair `#0`;
`--no-pair-check` correctly disables only the ID check while the record-
count check still fires.

## CLI shape

```
genorush fastx sample -i reads.fq.gz -p 0.1 -o sub.fq.gz -s 42
genorush fastx sample -i reads.fq.gz -n 50000 -o sub.fq.gz -s 42
genorush fastx sample -i R1.fq.gz -I R2.fq.gz -o R1.sub.fq.gz -O R2.sub.fq.gz -p 0.1 -s 42
```

Flag names (`-i/--in1`, `-I/--in2`, `-o/--out1`, `-O/--out2`) intentionally
mirror `fastp`'s convention rather than `seqkit`'s, since paired-end I/O
flag naming is what most users in this domain already have muscle memory
for. `-p`/`-n` are `clap`-enforced mutually exclusive (`conflicts_with`);
"exactly one of them is required" is checked manually in `run()` rather
than via a `clap` `ArgGroup`, since both are plain `Option<T>` with no
default value and the manual check reads clearly. `-r/--non-deterministic`
intentionally has no `conflicts_with` on `-s/--seed` — `-s` carries a
`default_value_t`, which clap treats as "always present," so a `conflicts_with`
between the two would misfire on the common case of just passing `-r`
alone; the precedence (`-r` wins if both given) is instead resolved once in
`effective_seed()`.

## Validated behavior

Tested against 200,000 synthetic 150 bp read pairs and 20,000 synthetic
long reads (500–3,000 bp, single-end): proportion sampling lands within
sampling noise of the target (9.75% observed for `-p 0.1`, as expected);
`-n` reservoir sampling returns *exactly* the requested count; paired
proportion sampling produces identical R1/R2 counts with zero ID mismatches
across the full output; gzip round-trips cleanly; malformed/desynced
paired input is rejected with a clear, specific error and process exit
code 1. `common::rng` and `common::fastq` — the two modules every future
record-oriented command will likely reuse — carry unit tests covering
determinism, uniformity, and FASTQ structural-validation edge cases.
