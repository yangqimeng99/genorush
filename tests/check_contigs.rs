//! `check contigs` against files the ecosystem actually produces.
//!
//! The fixtures under `tests/data/contigs/` were written by samtools and
//! bcftools, not by hand (see the README there). That matters more here than
//! anywhere else in this suite: the point of these tests is that the readers
//! cope with real BAM, CRAM and BCF, and a fixture invented to match the
//! reader would prove nothing about that.

mod common;

use common::*;

fn fixture(name: &str) -> String {
    data(&format!("contigs/{name}"))
        .to_str()
        .expect("fixture path")
        .to_string()
}

fn check(files: &[&str], extra: &[&str]) -> std::process::Output {
    let paths: Vec<String> = files.iter().map(|f| fixture(f)).collect();
    let mut args = vec!["check", "contigs"];
    args.extend(paths.iter().map(|s| s.as_str()));
    args.extend_from_slice(extra);
    run(&args)
}

#[test]
fn every_supported_format_is_read_and_agrees_with_the_reference() {
    let out = check(
        &[
            "ref.fa.fai",
            "aln.sam",
            "aln.bam",
            "aln.cram",
            "calls.vcf",
            "calls.bcf",
            "genes.gff3",
            "regions.bed",
        ],
        &[],
    );
    assert_ok(&out);
    let report = String::from_utf8_lossy(&out.stdout);

    // Each reader has to find all three contigs in the formats that state
    // them -- a header parsed as empty would otherwise look like agreement.
    for (file, kind, count) in [
        ("aln.sam", "SAM", 3),
        ("aln.bam", "BAM", 3),
        ("aln.cram", "CRAM", 3),
        ("calls.vcf", "VCF", 3),
        ("calls.bcf", "BCF", 3),
        ("genes.gff3", "GFF", 2),
        ("regions.bed", "BED", 2),
    ] {
        let line = report
            .lines()
            .find(|l| l.contains(file))
            .unwrap_or_else(|| panic!("{file} missing from the report:\n{report}"));
        assert!(line.contains(kind), "{file} read as the wrong kind: {line}");
        assert!(
            line.split_whitespace().any(|w| w == count.to_string()),
            "{file} should list {count} contigs: {line}"
        );
    }
    // Informational remarks are expected here (the GFF and BED cover part of
    // the genome); what must not appear is anything fatal.
    assert!(
        !report.contains("the reference does not have") && !report.contains("different length"),
        "no file should disagree with the reference:\n{report}"
    );
}

#[test]
fn lengths_come_through_for_the_formats_that_state_them() {
    let out = check(&["ref.fa.fai", "aln.bam", "calls.bcf", "regions.bed"], &[]);
    assert_ok(&out);
    let report = String::from_utf8_lossy(&out.stdout);

    for file in ["aln.bam", "calls.bcf"] {
        let line = report.lines().find(|l| l.contains(file)).unwrap();
        assert!(line.ends_with("yes"), "{file} should have lengths: {line}");
    }
    // A BED names contigs without sizing them, and saying "no" is different
    // from saying they disagree.
    let bed = report.lines().find(|l| l.contains("regions.bed")).unwrap();
    assert!(bed.contains("no"), "BED states no lengths: {bed}");
}

#[test]
fn a_chr_prefix_mismatch_fails_and_names_the_convention() {
    let out = check(&["ref.fa.fai", "chrprefix.vcf"], &[]);
    assert_refused(&out, "would change results silently");

    let report = String::from_utf8_lossy(&out.stdout);
    assert!(
        report.contains("the reference does not have"),
        "expected the unknown contigs to be listed:\n{report}"
    );
    // The hint is the useful part: "not found" sends someone hunting, "the
    // reference calls it 1" is a fix.
    assert!(
        report.contains(r#"chr1 -- the reference calls it "1""#),
        "expected a naming hint:\n{report}"
    );
    // And it must not also call this a partial file.
    assert!(
        !report.contains("normal for a file"),
        "a total naming mismatch is not a partial file:\n{report}"
    );
}

#[test]
fn a_length_mismatch_fails_and_says_it_is_a_different_build() {
    let out = check(&["ref.fa.fai", "wronglen.vcf"], &[]);
    assert_refused(&out, "would change results silently");

    let report = String::from_utf8_lossy(&out.stdout);
    assert!(
        report.contains("different genome builds"),
        "expected the build mismatch to be named:\n{report}"
    );
    assert!(
        report.contains("1501 here, 1500 in the reference"),
        "expected both lengths:\n{report}"
    );
}

#[test]
fn a_different_order_is_a_remark_until_it_is_asked_to_matter() {
    let relaxed = check(&["ref.fa.fai", "reordered.vcf"], &[]);
    assert_ok(&relaxed);
    let report = String::from_utf8_lossy(&relaxed.stdout);
    assert!(
        report.contains("different order") && report.contains("--require-order"),
        "expected a remark pointing at the flag:\n{report}"
    );

    let strict = check(&["ref.fa.fai", "reordered.vcf"], &["--require-order"]);
    assert_refused(&strict, "would change results silently");
}

#[test]
fn a_file_covering_part_of_the_genome_is_not_a_failure() {
    // A BED naming two of three contigs is the normal case, not a finding
    // worth failing on.
    let out = check(&["ref.fa.fai", "regions.bed"], &[]);
    assert_ok(&out);
    assert!(
        String::from_utf8_lossy(&out.stdout).contains("normal for a file"),
        "expected the partial coverage to be explained, not flagged"
    );
}

#[test]
fn an_unrecognised_extension_says_what_is_supported() {
    let dir = scratch("contigs_unknown");
    let odd = dir.join("mystery.txt");
    write(&odd, b"1\t2000\n");
    let out = run(&[
        "check",
        "contigs",
        &fixture("ref.fa.fai"),
        odd.to_str().unwrap(),
    ]);
    assert_refused(&out, "cannot tell what kind of file");
    assert_refused(&out, ".chrom.sizes");
}

#[test]
fn one_file_alone_is_not_a_comparison() {
    let out = run(&["check", "contigs", &fixture("ref.fa.fai")]);
    assert!(
        !out.status.success(),
        "comparing a file with nothing should be rejected"
    );
}
