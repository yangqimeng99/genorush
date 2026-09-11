# `fastx pair`: design and internals

Source: `src/fastx/pair.rs`.

## Motivation

Quality control breaks pairing. A trimmer drops a read from R1 for being too
short; its mate in R2 passes and stays. From that record on, the two files
disagree about which read sits at which position — and both files are still
perfectly well-formed FASTQ. Nothing errors. Every downstream tool that pairs
by position, which is most of them, then aligns reads against the wrong mates.

This is the same class of failure as the layouts `fastx deinterleave` refuses
to guess at, arriving from a different direction: not "which mate is this
record", but "does this record still have a mate at all".

## Why not index one file

The obvious implementation reads one file into a map keyed by read id, then
streams the other past it. It is what the available tools do, and it makes
memory a function of *input size*:

- `seqkit pair` has been reported consuming about 380 GB of RAM on a 49 GB +
  62 GB pair before dying, with 400 MB of output written
  ([shenwei356/seqkit#305](https://github.com/shenwei356/seqkit/issues/305),
  unanswered).
- BBMap's `repair.sh` needs a JVM and an up-front heap size.
- `fastq-pair` is small and fast but does not read gzip.

Real desync is almost never proportional to input size. A few percent of reads
get filtered; everything else stays in the same relative order. Paying for the
whole file to handle that is the wrong shape.

## The streaming join

Both files are read concurrently, one record at a time from each, strictly
alternating so the two read positions stay at the same index.

- A record whose mate has already been seen is written out **immediately**, as
  a pair, and neither is retained.
- A record whose mate has not been seen yet goes into `pending`.

Peak memory is therefore the number of records unmatched at the worst moment —
the true desync plus the orphans — and not the file size. Measured on a
2,000-record pair where each side is missing about 1% of the other's reads:

| input | peak records held |
|---|---|
| same order, a few reads missing from each side | 44 |
| one side fully reversed | 1,978 (one whole side) |
| one side fully reversed, spilled at an 8 KB budget | 84 |

The first row is the case this is built for. The second is what the design
must survive rather than optimise: with the mate of the first record sitting
at the far end of the other file, there is nothing to release until the very
end.

### Output alignment is free; input order is not

Because both records of a pair are written in the same step, output record `i`
of R1 and output record `i` of R2 are always the same pair. That alignment
costs nothing and holds in every regime — it is the whole point of the
command, not a constraint traded against memory.

Preserving the *input* order is a different property. It survives the
streaming join whenever the inputs were in the same relative order, because
matches then complete in order. On reordered input the output follows the
order in which pairs completed, and after a spill it follows the partitions.
The run says which happened rather than leaving it to be discovered.

### Why nothing is declared an orphan early

It is tempting to flush `pending` once the other stream has moved past where
a mate "should" have been, which would make memory proportional to desync
alone, orphans excluded. That is only sound if both files are in the same
relative order — exactly the assumption a re-pairing tool cannot make, since a
name-sorted or rebuilt file breaks it. A record flushed as an orphan whose
mate turns up later is silent data loss, the failure this command exists to
prevent. Orphanhood is decided at EOF, when it is provable.

## The on-disk fallback

When `pending` grows past `--max-memory`, holding it is not an option and
guessing is not either. The way out is the one databases use for this shape of
problem — a hash-partitioned (Grace) join.

Both remainders are written to `--partitions` files per side, chosen by
`hash(read id) % k`. Two mates hash identically, so **a pair can never be
split across partitions**, which means each partition joins on its own and
only one has to be resident at a time. Peak memory becomes a function of `k`
rather than of the input; partition files are deleted as the join consumes
them, so peak disk is the whole spill only at the moment the join starts.

Records are spilled in a small framing of their own — a `u64` original index,
then each of the four fields length-prefixed — rather than as FASTQ. The index
has to survive the trip so orphans come back out in input order, and lengths
are safer than delimiters when a header may contain anything, including tabs.

Spill files are written through the same `BlockWriter` every command uses, so
they are gzip-compressed in parallel, and read back through `open_reader`,
which detects that compression by content.

### Temporary files

The spill directory is created inside `--temp-dir` (default: the working
directory) and removed when the run ends, including on the error paths, via a
guard that runs on scope exit.

The directory is created before the free-space check, not after, and that
order is not incidental. Asking about a path that does not exist is not a
portable way to learn that it is unusable: POSIX `statvfs` reports `ENOENT`,
while Windows resolves any syntactically valid path to its volume root and
cheerfully answers about the volume. Creating the directory settles both
questions — does it exist, can we write there — on every platform, and the
guard removes it again if the space check then fails.

Free space is checked against the on-disk size of the inputs plus 10%. That
overshoots — only the unpaired remainder is ever spilled, compressed the way
the inputs are — and overshooting is the right direction: running out of
space halfway through a join leaves no answer and a directory full of
partitions.

One more platform difference shapes the join itself: Unix lets a file be
unlinked while it is still open, Windows does not. Partition files are
deleted as the join consumes them, so every writer and reader on a partition
has to be closed before its turn comes -- `Partitions::finish` consumes
itself to drop the writers, and each read phase is scoped.

## Validated behavior

`tests/pair_join.rs` runs the same inputs through all three regimes — in
memory and in order, in memory and reversed, and spilled to disk — and
requires the same set of pairs and the same two orphan lists from each, with
the outputs positionally aligned every time. A fallback that quietly paired
differently from the fast path would be worse than one that refused to run.

Also covered: that the fast path does not touch the disk, that the spill
directory is gone afterwards, that orphans with nowhere to go are still
counted and reported rather than silently dropped, that an unusable
`--temp-dir` is reported before any spilling starts, and that a `--max-memory`
which is not a size is rejected rather than read as zero.
