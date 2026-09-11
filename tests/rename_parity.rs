//! Byte-for-byte parity with the Python script `fastx rename` / `gff rename`
//! replace.
//!
//! The strongest claim these commands make is that their output is identical
//! to [`ChangeChrNameInFaOrGff.py`][script] from the SVLearn paper-code
//! repository -- not similar, identical -- because pipelines already built on
//! that script have to be able to swap it out without anything downstream
//! noticing. Until now that claim rested on a one-off manual comparison.
//!
//! `tests/data/rename/expected.{fa,gff}` were produced by running that script
//! over `input.{fa,gff}` with `names.tsv`, and are committed as fixtures. The
//! tests here run the real binary over the same inputs and require the bytes
//! to match. Regenerating them, if the inputs ever change:
//!
//! ```text
//! python3 ChangeChrNameInFaOrGff.py -i tests/data/rename/input.fa  --fa \
//!     -n tests/data/rename/names.tsv -o tests/data/rename/expected.fa
//! python3 ChangeChrNameInFaOrGff.py -i tests/data/rename/input.gff --gff \
//!     -n tests/data/rename/names.tsv -o tests/data/rename/expected.gff
//! ```
//!
//! The inputs deliberately cover every rule in the contract (see
//! `docs/en/rename.md`), including the ones that look like bugs and are kept
//! anyway: sequence descriptions dropped, a single-column GFF line gaining a
//! trailing tab, whitespace trimmed off every line. A golden comparison alone
//! would catch a break but not name it, so each rule also has its own test
//! asserting the specific line it governs.
//!
//! [script]: https://github.com/yangqimeng99/svlearn-paper-code/blob/main/scripts/ChangeChrNameInFaOrGff.py

mod common;

use common::*;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    data(&format!("rename/{name}"))
}

/// Runs one of the rename commands over the committed input and returns what
/// it wrote. `extra` carries per-test flags such as a chunk size.
fn rename(category: &str, input: &str, out: &Path, extra: &[&str]) {
    let names = fixture("names.tsv");
    let input = fixture(input);
    let mut args = vec![
        category,
        "rename",
        input.to_str().unwrap(),
        "-n",
        names.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ];
    args.extend_from_slice(extra);
    assert_ok(&run(&args));
}

/// A distinct output path per call: several tests run the same command, and
/// the harness runs them on parallel threads.
fn out_path(stem: &str) -> PathBuf {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static N: AtomicUsize = AtomicUsize::new(0);
    scratch("rename_parity").join(format!(
        "{stem}.{}.{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ))
}

fn rename_fasta(extra: &[&str]) -> Vec<u8> {
    let out = out_path("out.fa");
    rename("fastx", "input.fa", &out, extra);
    read(&out)
}

fn rename_gff(extra: &[&str]) -> Vec<u8> {
    let out = out_path("out.gff");
    rename("gff", "input.gff", &out, extra);
    read(&out)
}

/// Compares as text so a failure prints something readable rather than two
/// byte arrays.
#[track_caller]
fn assert_same(actual: &[u8], expected: &[u8], what: &str) {
    let actual = String::from_utf8_lossy(actual);
    let expected = String::from_utf8_lossy(expected);
    assert_eq!(actual, expected, "{what} diverged from the Python script");
}

#[test]
fn fasta_output_matches_the_python_script_byte_for_byte() {
    assert_same(&rename_fasta(&[]), &read(&fixture("expected.fa")), "FASTA");
}

#[test]
fn gff_output_matches_the_python_script_byte_for_byte() {
    assert_same(&rename_gff(&[]), &read(&fixture("expected.gff")), "GFF");
}

#[test]
fn chunk_size_does_not_change_the_output() {
    // Lines are transformed in parallel batches; a batch boundary must not be
    // observable in the result.
    let expected_fa = read(&fixture("expected.fa"));
    let expected_gff = read(&fixture("expected.gff"));
    for chunk in ["1", "2", "3", "1000"] {
        assert_same(
            &rename_fasta(&["--chunk-lines", chunk]),
            &expected_fa,
            &format!("FASTA at --chunk-lines {chunk}"),
        );
        assert_same(
            &rename_gff(&["--chunk-lines", chunk]),
            &expected_gff,
            &format!("GFF at --chunk-lines {chunk}"),
        );
    }
}

#[test]
fn gzip_input_produces_the_same_output() {
    let dir = scratch("parity_gz");
    let gz_input = dir.join("input.fa.gz");
    write(&gz_input, &gzip(&read(&fixture("input.fa"))));

    let names = fixture("names.tsv");
    let out = dir.join("out.fa");
    assert_ok(&run(&[
        "fastx",
        "rename",
        gz_input.to_str().unwrap(),
        "-n",
        names.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));
    assert_same(&read(&out), &read(&fixture("expected.fa")), "gzip FASTA");
}

// ---------------------------------------------------------------------------
// The individual rules, so a break says which one
// ---------------------------------------------------------------------------

fn lines(bytes: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(bytes)
        .lines()
        .map(|s| s.to_string())
        .collect()
}

#[test]
fn a_mapped_sequence_name_loses_its_description() {
    // `line.split()[0][1:]` in the original: everything after the first token
    // is discarded. Debatable, and deliberately preserved.
    let out = lines(&rename_fasta(&[]));
    assert_eq!(out[0], ">1", "expected the description to be dropped");
}

#[test]
fn an_unmapped_sequence_name_also_loses_its_description() {
    let out = lines(&rename_fasta(&[]));
    assert!(
        out.contains(&">chrUnmapped".to_string()),
        "an unmapped name should pass through without its description: {out:?}"
    );
}

#[test]
fn a_header_of_nothing_but_a_marker_survives() {
    let out = lines(&rename_fasta(&[]));
    assert!(
        out.contains(&">".to_string()),
        "a bare '>' should still be written: {out:?}"
    );
}

#[test]
fn whitespace_is_trimmed_from_every_line() {
    let out = lines(&rename_fasta(&[]));
    assert!(
        out.contains(&"acgt with surrounding whitespace".to_string()),
        "leading and trailing whitespace should be trimmed: {out:?}"
    );
    assert!(
        out.iter().all(|l| l == l.trim()),
        "no line should keep surrounding whitespace: {out:?}"
    );
}

#[test]
fn a_carriage_return_is_stripped_with_the_rest() {
    let raw = rename_fasta(&[]);
    assert!(
        !raw.contains(&b'\r'),
        "CRLF input should come out with bare newlines"
    );
}

#[test]
fn gff_comments_pass_through_untouched() {
    let out = lines(&rename_gff(&[]));
    assert_eq!(out[0], "##gff-version 3");
    assert_eq!(out[1], "#!genome-build unit-test");
    assert!(
        out.contains(&"# a trailing comment".to_string()),
        "comments anywhere in the file should survive: {out:?}"
    );
}

#[test]
fn gff_renames_only_the_first_column() {
    let out = lines(&rename_gff(&[]));
    assert_eq!(out[2], "1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;Name=a");
}

#[test]
fn an_unmapped_gff_seqid_leaves_the_line_alone() {
    let out = lines(&rename_gff(&[]));
    assert!(
        out.contains(&"chrUnmapped\tsrc\tgene\t1\t10\t.\t+\t.\tID=g2".to_string()),
        "an unmapped seqid should pass the whole line through: {out:?}"
    );
}

#[test]
fn a_single_column_gff_line_keeps_its_trailing_tab() {
    // `'\t'.join(LineList[1:])` on a one-element list is the empty string,
    // written after a tab that is still there. It looks like a bug and it is
    // the contract: `src/gff/rename.rs` reproduces it on purpose.
    let out = lines(&rename_gff(&[]));
    assert!(
        out.contains(&"1\t".to_string()),
        "a mapped single-column line should become `1\\t`: {out:?}"
    );
}

#[test]
fn a_space_separated_mapping_row_is_honoured() {
    // The table is whitespace-separated, not tab-separated: `3   chr3`.
    let out = lines(&rename_fasta(&[]));
    assert!(
        out.contains(&">3".to_string()),
        "a space-separated mapping row should map like any other: {out:?}"
    );
}
