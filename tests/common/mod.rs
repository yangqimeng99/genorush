//! Shared plumbing for the CLI integration tests: build fixtures, run the
//! binary with real pipes attached, and read the results back.
//!
//! Everything here is written in Rust rather than shelling out, because CI
//! runs these tests on Windows as well as Linux and macOS -- `cat`, `zcat`
//! and friends are not available on all three, and their flags differ where
//! they are.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Output, Stdio};

/// A directory for one test's fixtures and outputs.
///
/// Tests run in parallel threads, so each one gets its own directory, named
/// after the test. Nothing is deleted: `target/` is disposable, files are
/// truncated on rewrite, and a test that fails leaves its inputs and outputs
/// behind to look at.
pub fn scratch(test: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("it-scratch")
        .join(test);
    std::fs::create_dir_all(&dir).expect("failed to create scratch dir");
    dir
}

pub fn write(path: &Path, bytes: &[u8]) {
    std::fs::write(path, bytes).expect("failed to write fixture");
}

pub fn read(path: &Path) -> Vec<u8> {
    std::fs::read(path).expect("failed to read output")
}

/// One FASTQ record, 100 bp, with a header of the given text.
pub fn record(header: &str) -> String {
    let seq: String = (0..100)
        .map(|i| ['A', 'C', 'G', 'T'][(i * 7 + header.len()) % 4])
        .collect();
    format!("@{header}\n{seq}\n+\n{}\n", "?".repeat(100))
}

/// A FASTQ built from headers, in order.
pub fn fastq(headers: &[String]) -> Vec<u8> {
    headers
        .iter()
        .map(|h| record(h))
        .collect::<String>()
        .into_bytes()
}

/// SRA-style headers: a globally sequential id, then the same number again
/// with the mate marker -- `@SRR000.7 7/1`. The two mates of a pair do not
/// share an id, which is the shape that forces marker-based routing.
pub fn sra_headers(range: std::ops::RangeInclusive<u32>, mate: u8) -> Vec<String> {
    range.map(|i| format!("SRR000.{i} {i}/{mate}")).collect()
}

/// A small FASTA with two sequences per chromosome name.
pub fn fasta(names: &[&str]) -> Vec<u8> {
    let mut out = String::new();
    for name in names {
        out.push_str(&format!(">{name} some description\n"));
        for _ in 0..4 {
            out.push_str("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        }
    }
    out.into_bytes()
}

pub fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    enc.write_all(bytes).expect("failed to gzip fixture");
    enc.finish().expect("failed to finish gzip fixture")
}

pub fn gunzip(bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(bytes)
        .read_to_end(&mut out)
        .expect("output was not readable gzip");
    out
}

pub fn is_gzip(bytes: &[u8]) -> bool {
    bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b
}

fn command(args: &[&str]) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_genorush"));
    cmd.args(args);
    cmd
}

/// Runs the binary with no stdin, capturing stdout and stderr.
pub fn run(args: &[&str]) -> Output {
    command(args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("failed to run genorush")
}

/// Runs the binary with `input` fed to its stdin.
///
/// stdin is written from a separate thread: a command that streams its
/// output would otherwise be able to fill the stdout pipe while this side is
/// still writing, and both processes would wait on each other forever.
pub fn run_with_stdin(args: &[&str], input: &[u8]) -> Output {
    let mut child = command(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn genorush");

    let mut sink = child.stdin.take().expect("stdin was piped");
    let payload = input.to_vec();
    let writer = std::thread::spawn(move || {
        // A command that stops early (a rejected argument combination) closes
        // stdin while this is still writing; that is the expected outcome of
        // those tests, not a failure of the harness.
        let _ = sink.write_all(&payload);
    });

    let out = child
        .wait_with_output()
        .expect("failed to wait for genorush");
    writer.join().expect("stdin writer panicked");
    out
}

/// Spawns the binary with its stdout piped, for tests that care about what
/// happens to the process rather than what it produced.
pub fn spawn_piped(args: &[&str]) -> Child {
    command(args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("failed to spawn genorush")
}

#[track_caller]
pub fn assert_ok(out: &Output) {
    assert!(
        out.status.success(),
        "expected success, got {}\nstderr:\n{}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
}

/// Asserts the run failed and said why, with `needle` somewhere in stderr.
#[track_caller]
pub fn assert_refused(out: &Output, needle: &str) {
    assert!(
        !out.status.success(),
        "expected a refusal, but the command succeeded\nstdout:\n{}",
        String::from_utf8_lossy(&out.stdout)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains(needle),
        "refused, but not for the expected reason.\nlooking for: {needle}\nstderr:\n{stderr}"
    );
}
