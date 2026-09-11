//! End-to-end tests for reading stdin and writing stdout.
//!
//! These drive the built binary with real pipes rather than calling into the
//! code, because what they protect only exists at that level: argument
//! defaults, which stream a path resolves to, exit codes, and the split
//! between data on stdout and logs on stderr.
//!
//! The shape of almost every test is the same. Run the command the ordinary
//! way, to a file, and then run it again through a pipe; the two results must
//! be byte-identical. A pipeline that quietly produces something *slightly*
//! different from the file-based run is the failure worth guarding against,
//! so nothing here asserts on record counts or headers alone.
//!
//! Not covered here: refusing gzip output to a terminal. Attaching a pty
//! needs `script`, whose flags differ between Linux and macOS and which does
//! not exist on Windows, and CI runs all three. The rule is a pure function
//! of (stdout-bound, gzip, is-a-tty) and is unit-tested as such in
//! `io_utils`; only the wiring to the real `isatty` is untested.

mod common;

use common::*;

// ---------------------------------------------------------------------------
// A pipe must produce exactly what a file produces
// ---------------------------------------------------------------------------

#[test]
fn rename_through_a_pipe_matches_the_file_run() {
    let dir = scratch("rename_pipe");
    let input = dir.join("genome.fa");
    let map = dir.join("map.tsv");
    let out = dir.join("ref.fa");
    write(&input, &fasta(&["chr1", "chr2", "chr3"]));
    write(&map, b"1\tchr1\n2\tchr2\n3\tchr3\n");

    assert_ok(&run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--chunk-lines",
        "3",
    ]));

    let piped = run_with_stdin(
        &[
            "fastx",
            "rename",
            "-",
            "-n",
            map.to_str().unwrap(),
            "-o",
            "-",
            "--chunk-lines",
            "3",
        ],
        &fasta(&["chr1", "chr2", "chr3"]),
    );
    assert_ok(&piped);
    assert_eq!(piped.stdout, read(&out));
}

#[test]
fn gzip_input_on_a_pipe_is_detected_and_z_compresses_the_way_out() {
    let dir = scratch("rename_gzip_pipe");
    let map = dir.join("map.tsv");
    let out = dir.join("ref.fa");
    let body = fasta(&["chr1", "chr2"]);
    write(&map, b"1\tchr1\n2\tchr2\n");

    let input = dir.join("genome.fa");
    write(&input, &body);
    assert_ok(&run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));

    // Nothing on the command line says the input is compressed: it is
    // recognised from the bytes, which is the only option on a stream.
    let piped = run_with_stdin(
        &[
            "-z",
            "fastx",
            "rename",
            "-",
            "-n",
            map.to_str().unwrap(),
            "-o",
            "-",
        ],
        &gzip(&body),
    );
    assert_ok(&piped);
    assert!(is_gzip(&piped.stdout), "-z must compress stdout");
    assert_eq!(gunzip(&piped.stdout), read(&out));
}

/// Every BGZF stream ends with this: an empty member saying the file is
/// whole. Its absence is how a reader detects truncation, and its presence is
/// what separates an indexable `.gz` from one that merely decompresses.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

#[test]
fn bgzf_output_is_indexable_gzip_with_the_same_contents() {
    let dir = scratch("bgzf_flag");
    let input = dir.join("genome.fa");
    let map = dir.join("map.tsv");
    let plain = dir.join("plain.fa");
    let bgz = dir.join("indexed.fa.gz");
    write(&input, &fasta(&["chr1", "chr2"]));
    write(&map, b"1\tchr1\n2\tchr2\n");

    assert_ok(&run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        plain.to_str().unwrap(),
    ]));
    assert_ok(&run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        bgz.to_str().unwrap(),
        "--bgzf",
    ]));

    let raw = read(&bgz);
    assert!(is_gzip(&raw), "BGZF is gzip");
    assert_eq!(
        raw[3] & 0x04,
        0x04,
        "the extra field flag carries the block size"
    );
    assert_eq!(&raw[12..14], b"BC", "BGZF's own subfield should be there");
    assert_eq!(
        &raw[raw.len() - BGZF_EOF.len()..],
        &BGZF_EOF,
        "without the EOF block the file reads as truncated"
    );

    // Same bytes either way: choosing BGZF is a framing decision, not a
    // content one.
    assert_eq!(gunzip(&raw), read(&plain));

    // And the ordinary gzip path must not be quietly producing BGZF.
    let plain_gz = dir.join("plain2.fa.gz");
    assert_ok(&run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        plain_gz.to_str().unwrap(),
    ]));
    assert!(!read(&plain_gz).ends_with(&BGZF_EOF));
}

#[test]
fn gzip_flag_compresses_a_file_whose_name_does_not_end_in_gz() {
    let dir = scratch("forced_gzip_file");
    let input = dir.join("genome.fa");
    let map = dir.join("map.tsv");
    let plain = dir.join("plain.fa");
    let forced = dir.join("compressed.bin");
    write(&input, &fasta(&["chr1"]));
    write(&map, b"1\tchr1\n");

    assert_ok(&run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        plain.to_str().unwrap(),
    ]));
    assert_ok(&run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        forced.to_str().unwrap(),
        "-z",
    ]));

    let bytes = read(&forced);
    assert!(is_gzip(&bytes), "-z must compress regardless of the name");
    assert_eq!(gunzip(&bytes), read(&plain));
}

#[test]
fn sampling_by_proportion_keeps_the_same_reads_through_a_pipe() {
    let dir = scratch("sample_proportion");
    let input = dir.join("reads.fq");
    let out = dir.join("ref.fq");
    let body = fastq(&sra_headers(1..=500, 1));
    write(&input, &body);

    let common: Vec<&str> = vec!["fastx", "sample", "-p", "0.25", "-s", "42"];
    let mut file_args = common.clone();
    file_args.extend(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()]);
    assert_ok(&run(&file_args));

    let mut pipe_args = common;
    pipe_args.extend(["-i", "-", "-o", "-"]);
    let piped = run_with_stdin(&pipe_args, &body);
    assert_ok(&piped);
    // The keep/discard decision is a function of (seed, global index), so
    // the stream and the file must select exactly the same reads.
    assert_eq!(piped.stdout, read(&out));
    assert!(!piped.stdout.is_empty(), "0.25 of 500 reads should be some");
}

#[test]
fn sampling_an_exact_count_matches_through_a_pipe() {
    let dir = scratch("sample_number");
    let input = dir.join("reads.fq");
    let out = dir.join("ref.fq");
    let body = fastq(&sra_headers(1..=500, 1));
    write(&input, &body);

    assert_ok(&run(&[
        "fastx",
        "sample",
        "-n",
        "37",
        "-s",
        "7",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));

    let piped = run_with_stdin(
        &[
            "fastx", "sample", "-n", "37", "-s", "7", "-i", "-", "-o", "-",
        ],
        &body,
    );
    assert_ok(&piped);
    assert_eq!(piped.stdout, read(&out));
}

#[test]
fn interleave_takes_one_mate_from_stdin() {
    let dir = scratch("interleave_pipe");
    let r1_path = dir.join("R1.fq");
    let r2_path = dir.join("R2.fq");
    let out = dir.join("ref.fq");
    let r1 = fastq(&(1..=100).map(|i| format!("read{i}/1")).collect::<Vec<_>>());
    let r2 = fastq(&(1..=100).map(|i| format!("read{i}/2")).collect::<Vec<_>>());
    write(&r1_path, &r1);
    write(&r2_path, &r2);

    assert_ok(&run(&[
        "fastx",
        "interleave",
        "-i",
        r1_path.to_str().unwrap(),
        "-I",
        r2_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--chunk-records",
        "7",
    ]));

    let piped = run_with_stdin(
        &[
            "fastx",
            "interleave",
            "-i",
            "-",
            "-I",
            r2_path.to_str().unwrap(),
            "-o",
            "-",
            "--chunk-records",
            "7",
        ],
        &r1,
    );
    assert_ok(&piped);
    assert_eq!(piped.stdout, read(&out));
}

/// The layout that motivated `by-suffix`, in miniature: a run of R2, then
/// all of R1, then the rest of R2, with globally sequential ids so the two
/// mates of a pair never share one.
fn three_run_fixture() -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let head = sra_headers(1..=20, 2);
    let middle = sra_headers(21..=70, 1);
    let tail = sra_headers(71..=100, 2);

    let mut merged = head.clone();
    merged.extend(middle.clone());
    merged.extend(tail.clone());

    let mut r2 = head;
    r2.extend(tail);

    (fastq(&merged), fastq(&middle), fastq(&r2))
}

#[test]
fn deinterleave_auto_works_on_a_pipe_and_agrees_with_the_file_run() {
    let dir = scratch("deinterleave_pipe");
    let merged_path = dir.join("merged.fq");
    let (merged, want_r1, want_r2) = three_run_fixture();
    write(&merged_path, &merged);

    let file_r1 = dir.join("file_1.fq");
    let file_r2 = dir.join("file_2.fq");
    assert_ok(&run(&[
        "fastx",
        "deinterleave",
        "-i",
        merged_path.to_str().unwrap(),
        "-o",
        file_r1.to_str().unwrap(),
        "-O",
        file_r2.to_str().unwrap(),
        "--chunk-records",
        "9",
    ]));
    assert_eq!(read(&file_r1), want_r1, "file run: read 1");
    assert_eq!(read(&file_r2), want_r2, "file run: read 2");

    // Same input as a stream. `auto` has to probe the head and then keep
    // reading the same stream, since a pipe cannot be reopened.
    let pipe_r1 = dir.join("pipe_1.fq");
    let pipe_r2 = dir.join("pipe_2.fq");
    let piped = run_with_stdin(
        &[
            "fastx",
            "deinterleave",
            "-i",
            "-",
            "-o",
            pipe_r1.to_str().unwrap(),
            "-O",
            pipe_r2.to_str().unwrap(),
            "--chunk-records",
            "9",
        ],
        &merged,
    );
    assert_ok(&piped);
    assert_eq!(read(&pipe_r1), want_r1, "piped run: read 1");
    assert_eq!(read(&pipe_r2), want_r2, "piped run: read 2");
    assert!(
        String::from_utf8_lossy(&piped.stderr).contains("routing by marker"),
        "the probe should have chosen marker routing for this layout"
    );
}

#[test]
fn deinterleave_streams_one_mate_while_writing_the_other() {
    let dir = scratch("deinterleave_split_destinations");
    let merged_path = dir.join("merged.fq");
    let r2_out = dir.join("R2.fq");
    let (merged, want_r1, want_r2) = three_run_fixture();
    write(&merged_path, &merged);

    let out = run(&[
        "fastx",
        "deinterleave",
        "-i",
        merged_path.to_str().unwrap(),
        "--layout",
        "by-suffix",
        "-o",
        "-",
        "-O",
        r2_out.to_str().unwrap(),
    ]);
    assert_ok(&out);
    assert_eq!(out.stdout, want_r1);
    assert_eq!(read(&r2_out), want_r2);
}

#[test]
fn cat_takes_one_source_from_stdin() {
    let dir = scratch("cat_pipe");
    let a = dir.join("run1.fq");
    let b = dir.join("run2.fq");
    let out = dir.join("ref.fq");
    let first = fastq(&sra_headers(1..=50, 1));
    let second = fastq(&sra_headers(51..=100, 1));
    write(&a, &first);
    write(&b, &second);

    assert_ok(&run(&[
        "fastx",
        "cat",
        "--r1",
        a.to_str().unwrap(),
        "--r1",
        b.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));

    let piped = run_with_stdin(
        &[
            "fastx",
            "cat",
            "--r1",
            "-",
            "--r1",
            b.to_str().unwrap(),
            "-o",
            "-",
        ],
        &first,
    );
    assert_ok(&piped);
    assert_eq!(piped.stdout, read(&out));
}

#[test]
fn rescue_reads_a_truncated_stream() {
    let dir = scratch("rescue_pipe");
    let input = dir.join("trunc.fq");
    let out = dir.join("ref.fq");
    let mut body = fastq(&sra_headers(1..=100, 1));
    body.truncate(body.len() - 150); // cut through the final record

    write(&input, &body);
    // A partial rescue exits 3, by design, so success is not asserted here.
    let file_run = run(&[
        "fastx",
        "rescue",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_eq!(file_run.status.code(), Some(3), "expected a partial rescue");

    let piped = run_with_stdin(&["fastx", "rescue", "-i", "-", "-o", "-"], &body);
    assert_eq!(piped.status.code(), Some(3), "same outcome on a stream");
    assert_eq!(piped.stdout, read(&out));
    assert!(!piped.stdout.is_empty(), "most records were still intact");
}

// ---------------------------------------------------------------------------
// Omitting the output means stdout
// ---------------------------------------------------------------------------

#[test]
fn rename_without_an_output_writes_stdout() {
    let dir = scratch("default_rename");
    let input = dir.join("genome.fa");
    let map = dir.join("map.tsv");
    let out = dir.join("ref.fa");
    write(&input, &fasta(&["chr1", "chr2"]));
    write(&map, b"1\tchr1\n2\tchr2\n");

    assert_ok(&run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));
    let defaulted = run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
    ]);
    assert_ok(&defaulted);
    assert_eq!(defaulted.stdout, read(&out));
}

#[test]
fn sample_without_an_output_writes_stdout() {
    let dir = scratch("default_sample");
    let input = dir.join("reads.fq");
    let out = dir.join("ref.fq");
    write(&input, &fastq(&sra_headers(1..=200, 1)));

    assert_ok(&run(&[
        "fastx",
        "sample",
        "-i",
        input.to_str().unwrap(),
        "-p",
        "0.5",
        "-s",
        "1",
        "-o",
        out.to_str().unwrap(),
    ]));
    let defaulted = run(&[
        "fastx",
        "sample",
        "-i",
        input.to_str().unwrap(),
        "-p",
        "0.5",
        "-s",
        "1",
    ]);
    assert_ok(&defaulted);
    assert_eq!(defaulted.stdout, read(&out));
}

#[test]
fn interleave_without_an_output_writes_stdout() {
    let dir = scratch("default_interleave");
    let r1 = dir.join("R1.fq");
    let r2 = dir.join("R2.fq");
    let out = dir.join("ref.fq");
    write(
        &r1,
        &fastq(&(1..=20).map(|i| format!("r{i}/1")).collect::<Vec<_>>()),
    );
    write(
        &r2,
        &fastq(&(1..=20).map(|i| format!("r{i}/2")).collect::<Vec<_>>()),
    );

    assert_ok(&run(&[
        "fastx",
        "interleave",
        "-i",
        r1.to_str().unwrap(),
        "-I",
        r2.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));
    let defaulted = run(&[
        "fastx",
        "interleave",
        "-i",
        r1.to_str().unwrap(),
        "-I",
        r2.to_str().unwrap(),
    ]);
    assert_ok(&defaulted);
    assert_eq!(defaulted.stdout, read(&out));
}

#[test]
fn deinterleave_without_an_output_streams_read_one() {
    let dir = scratch("default_deinterleave");
    let merged_path = dir.join("merged.fq");
    let r2_out = dir.join("R2.fq");
    let (merged, want_r1, want_r2) = three_run_fixture();
    write(&merged_path, &merged);

    // Only -O is required: read 1 streams out, read 2 goes to the file.
    let out = run(&[
        "fastx",
        "deinterleave",
        "-i",
        merged_path.to_str().unwrap(),
        "-O",
        r2_out.to_str().unwrap(),
    ]);
    assert_ok(&out);
    assert_eq!(out.stdout, want_r1);
    assert_eq!(read(&r2_out), want_r2);
}

#[test]
fn cat_without_an_output_writes_stdout() {
    let dir = scratch("default_cat");
    let a = dir.join("run1.fq");
    let out = dir.join("ref.fq");
    write(&a, &fastq(&sra_headers(1..=40, 1)));

    assert_ok(&run(&[
        "fastx",
        "cat",
        "--r1",
        a.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]));
    let defaulted = run(&["fastx", "cat", "--r1", a.to_str().unwrap()]);
    assert_ok(&defaulted);
    assert_eq!(defaulted.stdout, read(&out));
}

// ---------------------------------------------------------------------------
// Combinations that cannot work are refused, not half-done
// ---------------------------------------------------------------------------

#[test]
fn refuses_two_inputs_reading_the_same_stdin() {
    let out = run(&[
        "fastx", "sample", "-i", "-", "-I", "-", "-o", "a.fq", "-O", "b.fq", "-p", "0.5",
    ]);
    assert_refused(&out, "only one stdin/stdout");
}

#[test]
fn refuses_two_outputs_sharing_one_pipe() {
    let dir = scratch("reject_two_stdout");
    let r1 = dir.join("R1.fq");
    let r2 = dir.join("R2.fq");
    write(
        &r1,
        &fastq(&(1..=4).map(|i| format!("r{i}/1")).collect::<Vec<_>>()),
    );
    write(
        &r2,
        &fastq(&(1..=4).map(|i| format!("r{i}/2")).collect::<Vec<_>>()),
    );

    let out = run(&[
        "fastx",
        "sample",
        "-i",
        r1.to_str().unwrap(),
        "-I",
        r2.to_str().unwrap(),
        "-o",
        "-",
        "-O",
        "-",
        "-p",
        "0.5",
    ]);
    assert_refused(&out, "only one stdin/stdout");
}

#[test]
fn refuses_to_send_both_mates_of_a_split_to_stdout() {
    let dir = scratch("reject_deinterleave_stdout");
    let merged_path = dir.join("merged.fq");
    let (merged, _, _) = three_run_fixture();
    write(&merged_path, &merged);

    let out = run(&[
        "fastx",
        "deinterleave",
        "-i",
        merged_path.to_str().unwrap(),
        "-o",
        "-",
        "-O",
        "-",
    ]);
    assert_refused(&out, "only one stdin/stdout");
}

#[test]
fn refuses_two_stdin_sources_in_cat() {
    let out = run(&[
        "fastx", "cat", "--r1", "-", "--r2", "-", "-o", "a.fq", "-O", "b.fq",
    ]);
    assert_refused(&out, "only one stdin/stdout");
}

#[test]
fn refuses_stdin_for_both_the_input_and_the_mapping_table() {
    let out = run(&["fastx", "rename", "-", "-n", "-", "-o", "out.fa"]);
    assert_refused(&out, "only one stdin/stdout");
}

#[test]
fn refuses_a_layout_that_would_need_to_read_the_pipe_twice() {
    let dir = scratch("reject_concat_on_pipe");
    let r2_out = dir.join("R2.fq");
    let (merged, _, _) = three_run_fixture();

    let out = run_with_stdin(
        &[
            "fastx",
            "deinterleave",
            "-i",
            "-",
            "--layout",
            "concat",
            "-O",
            r2_out.to_str().unwrap(),
        ],
        &merged,
    );
    assert_refused(&out, "a pipe can only be read once");
    // The message has to leave the user somewhere to go.
    assert_refused(&out, "by-suffix");
}

#[test]
fn an_unopenable_input_is_still_an_error() {
    let dir = scratch("missing_input");
    let map = dir.join("map.tsv");
    write(&map, b"1\tchr1\n");
    let missing = dir.join("does-not-exist.fa");

    let out = run(&[
        "fastx",
        "rename",
        missing.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
        "-o",
        "-",
    ]);
    assert_refused(&out, "failed to open input file");
}

// ---------------------------------------------------------------------------
// Behaving like a well-mannered pipeline citizen
// ---------------------------------------------------------------------------

#[test]
fn closing_the_output_pipe_early_ends_the_run_cleanly() {
    use std::io::Read;

    let dir = scratch("broken_pipe");
    let input = dir.join("genome.fa");
    let map = dir.join("map.tsv");
    // Big enough that the child is still writing when the reader goes away:
    // it fills the pipe buffer and blocks long before it is finished.
    let names: Vec<String> = (1..=4000).map(|i| format!("chr{i}")).collect();
    let refs: Vec<&str> = names.iter().map(|s| s.as_str()).collect();
    write(&input, &fasta(&refs));
    write(&map, b"1\tchr1\n");

    let mut child = spawn_piped(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
    ]);

    let mut stdout = child.stdout.take().expect("stdout was piped");
    let mut head = [0u8; 64];
    stdout
        .read_exact(&mut head)
        .expect("failed to read the head");
    drop(stdout); // this is what `| head` does

    let status = child.wait().expect("failed to wait for genorush");
    assert!(
        status.success(),
        "a reader that stops early is not a failure, got {status}"
    );
}

#[test]
fn stdout_carries_data_only_and_logs_go_to_stderr() {
    let dir = scratch("stream_separation");
    let input = dir.join("genome.fa");
    let map = dir.join("map.tsv");
    write(&input, &fasta(&["chr1"]));
    write(&map, b"1\tchr1\n");

    let out = run(&[
        "fastx",
        "rename",
        input.to_str().unwrap(),
        "-n",
        map.to_str().unwrap(),
    ]);
    assert_ok(&out);

    let stdout = String::from_utf8(out.stdout).expect("output should be text here");
    assert!(
        stdout.starts_with(">1"),
        "stdout should open with the renamed record, got: {:?}",
        &stdout[..stdout.len().min(40)]
    );
    assert!(
        !stdout.contains("INFO"),
        "log lines must not contaminate piped data:\n{stdout}"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("INFO") && stderr.contains("lines processed"),
        "the run's logs should still be there, on stderr:\n{stderr}"
    );
}
