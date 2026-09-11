//! `fastx pair` across its three regimes.
//!
//! The command has one job -- decide which records have mates -- and three
//! very different ways of doing it depending on how far apart the inputs have
//! drifted: entirely in memory when they are merely missing reads, still in
//! memory but holding a whole side when they are reordered, and through disk
//! partitions when that no longer fits.
//!
//! The regimes exist for memory, not for semantics, so the tests here run the
//! same inputs through all three and require the same answer from each. A
//! fallback that quietly pairs differently from the fast path would be worse
//! than one that refuses to run.

mod common;

use common::*;
use std::path::Path;

/// Ids present in `keep`, as a mate file. Headers are SRA-shaped so the two
/// mates of a pair do not share their full header, only the id.
fn mate_file(dir: &Path, name: &str, ids: &[u32], mate: u8) -> std::path::PathBuf {
    let headers: Vec<String> = ids
        .iter()
        .map(|i| format!("@RD{i:06} {i}/{mate}"))
        .collect();
    let refs: Vec<&str> = headers.iter().map(|s| s.as_str()).collect();
    let path = dir.join(name);
    write(&path, &fastq_from(&refs));
    path
}

fn fastq_from(headers: &[&str]) -> Vec<u8> {
    headers
        .iter()
        .map(|h| format!("{h}\nACGTACGTAC\n+\nIIIIIIIIII\n"))
        .collect::<String>()
        .into_bytes()
}

fn ids_of(path: &Path) -> Vec<u32> {
    String::from_utf8_lossy(&read(path))
        .lines()
        .filter(|l| l.starts_with("@RD"))
        .map(|l| l[3..9].parse().expect("id"))
        .collect()
}

/// R1 keeps everything except multiples of 7; R2 everything except multiples
/// of 5. Neither side is a subset of the other, so both orphan lists matter.
struct Fixture {
    r1: Vec<u32>,
    r2: Vec<u32>,
    paired: Vec<u32>,
    orphan1: Vec<u32>,
    orphan2: Vec<u32>,
}

fn fixture_ids(n: u32) -> Fixture {
    Fixture {
        r1: (1..=n).filter(|i| i % 7 != 0).collect(),
        r2: (1..=n).filter(|i| i % 5 != 0).collect(),
        paired: (1..=n).filter(|i| i % 7 != 0 && i % 5 != 0).collect(),
        orphan1: (1..=n).filter(|i| i % 7 != 0 && i % 5 == 0).collect(),
        orphan2: (1..=n).filter(|i| i % 5 != 0 && i % 7 == 0).collect(),
    }
}

struct Outcome {
    pairs1: Vec<u32>,
    pairs2: Vec<u32>,
    orphans1: Vec<u32>,
    orphans2: Vec<u32>,
    stderr: String,
}

fn run_pair(dir: &Path, in1: &Path, in2: &Path, extra: &[&str]) -> Outcome {
    let (o1, o2) = (dir.join("p1.fq"), dir.join("p2.fq"));
    let (u1, u2) = (dir.join("u1.fq"), dir.join("u2.fq"));
    let mut args = vec![
        "fastx",
        "pair",
        "-i",
        in1.to_str().unwrap(),
        "-I",
        in2.to_str().unwrap(),
        "-o",
        o1.to_str().unwrap(),
        "-O",
        o2.to_str().unwrap(),
        "-u",
        u1.to_str().unwrap(),
        "-U",
        u2.to_str().unwrap(),
    ];
    args.extend_from_slice(extra);
    let out = run(&args);
    assert_ok(&out);
    Outcome {
        pairs1: ids_of(&o1),
        pairs2: ids_of(&o2),
        orphans1: ids_of(&u1),
        orphans2: ids_of(&u2),
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
    }
}

#[track_caller]
fn assert_correct(got: &Outcome, paired: &[u32], orphan1: &[u32], orphan2: &[u32]) {
    // The two outputs must line up record for record: that is the entire
    // point of re-pairing, and it holds however the answer was computed.
    assert_eq!(
        got.pairs1, got.pairs2,
        "outputs are not positionally aligned"
    );

    let mut sorted = got.pairs1.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, paired, "wrong set of pairs");

    let mut o1 = got.orphans1.clone();
    o1.sort_unstable();
    assert_eq!(o1, orphan1, "wrong read 1 orphans");
    let mut o2 = got.orphans2.clone();
    o2.sort_unstable();
    assert_eq!(o2, orphan2, "wrong read 2 orphans");
}

#[test]
fn same_order_inputs_pair_in_memory_and_keep_their_order() {
    let dir = scratch("pair_same_order");
    let fx = fixture_ids(300);
    let f1 = mate_file(&dir, "R1.fq", &fx.r1, 1);
    let f2 = mate_file(&dir, "R2.fq", &fx.r2, 2);

    let got = run_pair(&dir, &f1, &f2, &[]);
    assert_correct(&got, &fx.paired, &fx.orphan1, &fx.orphan2);

    // Nothing was reordered, so the output should still be ascending.
    let mut ascending = got.pairs1.clone();
    ascending.sort_unstable();
    assert_eq!(got.pairs1, ascending, "input order should survive");
    assert!(
        !got.stderr.contains("on disk"),
        "lightly desynced input must not need the disk: {}",
        got.stderr
    );
}

#[test]
fn a_reversed_mate_file_still_pairs_correctly_in_memory() {
    let dir = scratch("pair_reversed");
    let mut fx = fixture_ids(300);
    fx.r2.reverse();
    let f1 = mate_file(&dir, "R1.fq", &fx.r1, 1);
    let f2 = mate_file(&dir, "R2rev.fq", &fx.r2, 2);

    let got = run_pair(&dir, &f1, &f2, &[]);
    assert_correct(&got, &fx.paired, &fx.orphan1, &fx.orphan2);
}

#[test]
fn spilling_to_disk_reaches_the_same_answer() {
    let dir = scratch("pair_spilled");
    let mut fx = fixture_ids(300);
    fx.r2.reverse();
    let f1 = mate_file(&dir, "R1.fq", &fx.r1, 1);
    let f2 = mate_file(&dir, "R2rev.fq", &fx.r2, 2);

    // A budget far below one side forces the partitioned join; the partition
    // count keeps each partition inside that same budget.
    let got = run_pair(
        &dir,
        &f1,
        &f2,
        &[
            "--max-memory",
            "8K",
            "--partitions",
            "32",
            "--temp-dir",
            dir.to_str().unwrap(),
        ],
    );
    assert_correct(&got, &fx.paired, &fx.orphan1, &fx.orphan2);
    assert!(
        got.stderr.contains("Finishing the join on disk"),
        "the run should say it spilled: {}",
        got.stderr
    );

    // Deleting a file that is still open succeeds on Unix and fails on
    // Windows, so a handle left open by the join is invisible here unless the
    // warning it produces is treated as a failure.
    assert!(
        !got.stderr.contains("could not remove"),
        "spill files should all be deletable by the time the join is done: {}",
        got.stderr
    );

    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.file_name()
                .to_string_lossy()
                .starts_with("genorush-pair-")
        })
        .collect();
    assert!(
        leftovers.is_empty(),
        "the spill directory should be gone, found {} leftovers",
        leftovers.len()
    );
}

#[test]
fn orphans_are_counted_even_when_there_is_nowhere_to_put_them() {
    let dir = scratch("pair_dropped_orphans");
    let fx = fixture_ids(100);
    let f1 = mate_file(&dir, "R1.fq", &fx.r1, 1);
    let f2 = mate_file(&dir, "R2.fq", &fx.r2, 2);
    let (o1, o2) = (dir.join("p1.fq"), dir.join("p2.fq"));

    let out = run(&[
        "fastx",
        "pair",
        "-i",
        f1.to_str().unwrap(),
        "-I",
        f2.to_str().unwrap(),
        "-o",
        o1.to_str().unwrap(),
        "-O",
        o2.to_str().unwrap(),
    ]);
    assert_ok(&out);
    assert_eq!(ids_of(&o1).len(), fx.paired.len());

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("without a mate") && stderr.contains("were dropped"),
        "dropping orphans silently would hide data loss: {stderr}"
    );
}

#[test]
fn a_temp_dir_that_cannot_be_used_is_reported_before_anything_spills() {
    let dir = scratch("pair_bad_temp");
    let mut fx = fixture_ids(200);
    fx.r2.reverse();
    let f1 = mate_file(&dir, "R1.fq", &fx.r1, 1);
    let f2 = mate_file(&dir, "R2rev.fq", &fx.r2, 2);

    let out = run(&[
        "fastx",
        "pair",
        "-i",
        f1.to_str().unwrap(),
        "-I",
        f2.to_str().unwrap(),
        "-o",
        dir.join("p1.fq").to_str().unwrap(),
        "-O",
        dir.join("p2.fq").to_str().unwrap(),
        "--max-memory",
        "8K",
        "--temp-dir",
        dir.join("no/such/place").to_str().unwrap(),
    ]);
    assert_refused(&out, "free space");
}

#[test]
fn a_memory_budget_that_is_not_a_size_is_rejected() {
    let out = run(&[
        "fastx",
        "pair",
        "-i",
        "a.fq",
        "-I",
        "b.fq",
        "-o",
        "p1.fq",
        "-O",
        "p2.fq",
        "--max-memory",
        "plenty",
    ]);
    assert_refused(&out, "--max-memory");
}
