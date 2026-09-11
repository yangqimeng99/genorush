//! Helpers shared by the in-crate unit tests.
//!
//! Test-only: the module is `#[cfg(test)]`, so none of this is compiled into
//! the binary.

use std::fs;
use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::common::fastq::read_fastq_record;

/// A unique path for a test artifact.
///
/// Scratch space lives under the build directory rather than the system temp
/// dir: it is already git-ignored, it is guaranteed writable wherever `cargo`
/// itself can write, and it keeps test artifacts next to the build they came
/// from instead of scattered in /tmp. The counter and pid keep tests that run
/// in parallel -- or two `cargo test` runs at once -- off each other's files.
pub fn scratch_path(stem: &str) -> PathBuf {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("test-scratch");
    fs::create_dir_all(&dir).expect("failed to create test scratch dir");
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    dir.join(format!("{stem}.{n}.{}.fastq", std::process::id()))
}

/// Writes a FASTQ whose records carry the given headers, with a fixed body so
/// tests can assert on routing and ordering alone.
pub fn write_fastq(path: &Path, headers: &[&str]) {
    let mut body = String::new();
    for h in headers {
        body.push_str(&format!("{h}\nACGT\n+\nIIII\n"));
    }
    fs::write(path, body).expect("failed to write test fastq");
}

/// The headers of every record in a FASTQ, parsed back through the same
/// reader the commands use, so a structurally broken output fails here rather
/// than silently comparing equal.
pub fn headers_of(path: &Path) -> Vec<String> {
    let body = fs::read_to_string(path).expect("failed to read output");
    let mut reader = Cursor::new(body.into_bytes());
    let mut out = Vec::new();
    let mut line_no = 1;
    while let Some(rec) = read_fastq_record(&mut reader, line_no).expect("malformed output") {
        out.push(rec.header);
        line_no += 4;
    }
    out
}
