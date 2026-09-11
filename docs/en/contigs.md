# `check contigs`: design and internals

Source: `src/check/contigs.rs`, `src/common/contigs.rs`.

## Motivation

A BAM aligned against a reference whose contigs are called `1, 2, 3` and a VCF
annotated against one calling them `chr1, chr2, chr3` are both perfectly valid
files. Nothing about either is wrong. Tools that join them do not usually
crash — they find nothing for the contigs they cannot match, and report that
as zero, which is indistinguishable from a real answer.

Lengths are worse. Two builds of the same assembly can differ by a handful of
bases on one chromosome; every coordinate past that point is off by that much,
and nothing says so. GATK reports this class as "incompatible contigs" when it
happens to look, and the usual remedy is a round trip through Picard. This
command is the check on its own, ahead of time.

It is also where the project's format-reading layer starts, so the
command's other job is to prove that layer works on files the ecosystem
actually produces.

## What counts as a disagreement

The first file is the reference; every other is checked against it. That
asymmetry matters, because the two directions mean different things:

- **A contig the reference does not have** is fatal. Whatever consumes both
  files will silently find nothing there. The message names the convention
  when it can tell — `chr1 -- the reference calls it "1"` turns "not found"
  into a fix, and the `chr` prefix and the `MT`/`chrM` spelling account for
  most of these.
- **A length that differs** is fatal, and reported as what it is: two
  different genome builds.
- **A different order** is a remark. It matters to tools that require
  coordinate-sorted input against a fixed dictionary and is irrelevant to most
  others, so `--require-order` promotes it rather than the default failing.
- **Contigs the reference has that a file does not mention** is normal — a BED
  names only what it covers — so it is counted, not flagged. Unless that file
  *also* used names the reference lacks: then it is not a partial file, it is
  the naming problem seen from the other side, and calling it normal would
  mislead.

## Reading the contigs out of nine formats

`common::contigs` turns every supported format into the same shape: a name,
and a length where one is stated. Where the work happens differs sharply, and
so does the approach:

| format | source of truth | read with |
|---|---|---|
| `.fai`, `.chrom.sizes`, `.genome` | first two columns | directly |
| FASTA | the sibling `.fai` if present, else a full pass | directly |
| VCF | `##contig` header lines | `noodles-vcf` |
| BCF | the same, in the binary header | `noodles-bcf` |
| SAM | `@SQ` records | `noodles-sam` |
| BAM | the same, in the binary header | `noodles-bam` |
| CRAM | the same, in the file header container | `noodles-cram` |
| GFF/GTF | `##sequence-region` pragmas, else column 1 | directly |
| BED | column 1 | directly |

The split is deliberate. VCF, BCF, SAM, BAM and CRAM headers are real parsing
problems, and three of them are binary formats with their own compression
framing — that is what a library is for, and this is where `noodles` enters
the project. A GFF's contig names, on the other hand, are the first field of a
tab-separated line. Taking a version-coupled parser for `split('\t').next()`
would be a cost with no benefit; when a command needs to *manipulate* GFF
rather than glance at it, that is the time to add one.

Two details are easy to get wrong:

- **BCF, BAM and CRAM are opened as files**, not through the gzip-sniffing
  reader every text format here uses. They carry their own framing, and handing
  an already-decompressed stream to a reader that expects to do its own
  decompression gets nowhere.
- **CRAM's `read_header` expects the stream at the very start.** It checks the
  magic number and steps over the file definition itself, so calling
  `read_file_definition` first — which the neighbouring method invites —
  leaves it looking at the wrong bytes, and it reports an invalid header.

FASTA prefers a sibling `.fai` when one exists. That is what an index is for,
and it turns a full pass over a multi-gigabyte genome into reading a few
kilobytes.

## Validated behavior

`tests/check_contigs.rs` runs against fixtures written by the reference
implementations rather than by hand — `samtools view -b`, `samtools view -C`,
`bcftools view -Ob`, `samtools faidx` — because the point of these tests is
that the readers cope with real BAM, CRAM and BCF, and a fixture invented to
match the reader would prove nothing about that. Provenance is recorded in
`tests/data/contigs/README.md`.

Covered: every format read with the right contig count and lengths where they
exist; a `chr`-prefix mismatch failing with the naming hint and *not* also
being described as partial coverage; a length mismatch failing and naming the
cause; a reordering staying a remark until `--require-order`; partial coverage
staying a remark; and an unrecognised extension listing what is supported.
