Fixtures for `tests/check_contigs.rs`.

The binary ones were produced by the reference implementations rather than
written by hand, so the readers are tested against what the ecosystem actually
emits:

    samtools faidx ref.fa                       -> ref.fa.fai
    samtools view -b aln.sam                    -> aln.bam
    samtools view -C -T ref.fa aln.sam          -> aln.cram
    bcftools view -Ob -o calls.bcf calls.vcf    -> calls.bcf

with samtools 1.21 and bcftools 1.21. The reference genome is three short
contigs named `1`, `2` and `MT`; `chrprefix.vcf`, `wronglen.vcf` and
`reordered.vcf` each disagree with it in exactly one way.
