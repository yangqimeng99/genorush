//! Minimal FASTQ record model: read/write a 4-line record verbatim, with just
//! enough structural validation to catch truncated or corrupt input early.
//!
//! Also hosts the concurrent-mate-reading infrastructure shared by every
//! command that processes single- or paired-end FASTQ (`fastx sample`,
//! `fastx rescue`): `spawn_reader` decompresses and parses one file on its
//! own thread, and `recv_pair_step` consumes two such threads in lockstep,
//! classifying every possible outcome (a good pair, synchronized EOF, a
//! count mismatch, an ID mismatch, or a read error) so callers can decide
//! for themselves whether a given outcome is fatal (`sample`, which bails)
//! or just the point where recovery should stop (`rescue`, which doesn't).

use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver};
use std::thread;

use anyhow::{bail, Result};

use crate::io_utils::open_reader;

#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub header: String, // includes leading '@'
    pub seq: String,
    pub plus: String, // includes leading '+'
    pub qual: String,
}

impl FastqRecord {
    pub fn write_to(&self, w: &mut dyn Write) -> std::io::Result<()> {
        writeln!(w, "{}", self.header)?;
        writeln!(w, "{}", self.seq)?;
        writeln!(w, "{}", self.plus)?;
        writeln!(w, "{}", self.qual)
    }

    /// The read/pair identifier: the first whitespace-delimited token after
    /// `@`, with a trailing `/1` or `/2` mate suffix stripped if present.
    /// Covers both legacy Illumina (`@ID/1`) and modern (`@ID 1:N:0:...`)
    /// header conventions.
    pub fn base_id(&self) -> &str {
        let id = self.header.strip_prefix('@').unwrap_or(&self.header);
        let first = id.split_whitespace().next().unwrap_or(id);
        first
            .strip_suffix("/1")
            .or_else(|| first.strip_suffix("/2"))
            .unwrap_or(first)
    }
}

/// Reads one 4-line FASTQ record. Returns `Ok(None)` at a clean EOF (nothing
/// read before hitting end of stream), or an error on a truncated record or
/// a header/plus-line/seq-qual-length structural violation.
pub fn read_fastq_record(reader: &mut dyn BufRead, line_no: u64) -> Result<Option<FastqRecord>> {
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(None);
    }
    let header = header.trim_end_matches(['\n', '\r']).to_string();
    if !header.starts_with('@') {
        bail!("malformed FASTQ at line {line_no}: header does not start with '@': {header:?}");
    }

    let mut seq = String::new();
    if reader.read_line(&mut seq)? == 0 {
        bail!(
            "truncated FASTQ record at line {}: missing sequence line",
            line_no + 1
        );
    }
    let seq = seq.trim_end_matches(['\n', '\r']).to_string();

    let mut plus = String::new();
    if reader.read_line(&mut plus)? == 0 {
        bail!(
            "truncated FASTQ record at line {}: missing '+' line",
            line_no + 2
        );
    }
    let plus = plus.trim_end_matches(['\n', '\r']).to_string();
    if !plus.starts_with('+') {
        bail!(
            "malformed FASTQ at line {}: expected '+' line, got: {plus:?}",
            line_no + 2
        );
    }

    let mut qual = String::new();
    if reader.read_line(&mut qual)? == 0 {
        bail!(
            "truncated FASTQ record at line {}: missing quality line",
            line_no + 3
        );
    }
    let qual = qual.trim_end_matches(['\n', '\r']).to_string();

    if seq.len() != qual.len() {
        bail!(
            "malformed FASTQ record ending at line {}: sequence length ({}) != quality length ({}) for read {:?}",
            line_no + 3,
            seq.len(),
            qual.len(),
            header
        );
    }

    Ok(Some(FastqRecord {
        header,
        seq,
        plus,
        qual,
    }))
}

/// Spawns a background thread that decompresses and parses `path`, streaming
/// parsed records out through a bounded channel. Putting this on its own
/// thread is what lets two mate files decompress concurrently instead of
/// sequentially. On a read/parse error the thread sends that one `Err` and
/// then stops — it does not retry or skip ahead, since a caller that wants
/// to keep reading past a bad record wouldn't know which bytes to resync on.
pub fn spawn_reader(path: PathBuf) -> Result<Receiver<Result<FastqRecord>>> {
    let mut reader = open_reader(&path)?;
    let (tx, rx) = mpsc::sync_channel::<Result<FastqRecord>>(4096);
    thread::spawn(move || {
        let mut line_no: u64 = 1;
        loop {
            match read_fastq_record(reader.as_mut(), line_no) {
                Ok(Some(rec)) => {
                    line_no += 4;
                    if tx.send(Ok(rec)).is_err() {
                        break;
                    }
                }
                Ok(None) => break,
                Err(e) => {
                    let _ = tx.send(Err(e));
                    break;
                }
            }
        }
    });
    Ok(rx)
}

/// Every possible outcome of pulling one record from each of two mate
/// channels in lockstep. Deliberately has no notion of "fatal" — that's a
/// judgment call left to the caller (see `sample::recv_pair`, which turns
/// `CountMismatch`/`ReadError`/a false `ids_match` into hard errors, vs.
/// `rescue`, which turns the same outcomes into "stop here, keep what we
/// already have").
pub enum PairStep {
    /// Both channels produced a record. `ids_match` reports whether
    /// `FastqRecord::base_id()` agreed between the two.
    Pair {
        r1: FastqRecord,
        r2: FastqRecord,
        ids_match: bool,
    },
    /// Both channels closed at the same step: a clean, synchronized end of input.
    Eof,
    /// One channel closed while the other still had a record: R1/R2 have a
    /// different number of reads.
    CountMismatch,
    /// One (or both) of the underlying readers hit a read/parse error —
    /// truncated data, corrupt gzip, or a malformed record.
    ReadError(anyhow::Error),
}

pub fn recv_pair_step(
    rx1: &Receiver<Result<FastqRecord>>,
    rx2: &Receiver<Result<FastqRecord>>,
) -> PairStep {
    match (rx1.recv(), rx2.recv()) {
        (Ok(Err(e)), _) => PairStep::ReadError(e),
        (_, Ok(Err(e))) => PairStep::ReadError(e),
        (Ok(Ok(r1)), Ok(Ok(r2))) => {
            let ids_match = r1.base_id() == r2.base_id();
            PairStep::Pair { r1, r2, ids_match }
        }
        (Err(_), Err(_)) => PairStep::Eof,
        _ => PairStep::CountMismatch,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn rec(header: &str) -> FastqRecord {
        FastqRecord {
            header: header.to_string(),
            seq: "ACGT".to_string(),
            plus: "+".to_string(),
            qual: "IIII".to_string(),
        }
    }

    #[test]
    fn reads_one_record_then_eof() {
        let mut c = Cursor::new(b"@r1 1:N:0:AT\nACGT\n+\nIIII\n".to_vec());
        let rec = read_fastq_record(&mut c, 1).unwrap().unwrap();
        assert_eq!(rec.header, "@r1 1:N:0:AT");
        assert_eq!(rec.seq, "ACGT");
        assert_eq!(rec.plus, "+");
        assert_eq!(rec.qual, "IIII");
        assert!(read_fastq_record(&mut c, 5).unwrap().is_none());
    }

    #[test]
    fn reads_multiple_records_in_sequence() {
        let mut c = Cursor::new(b"@a\nAC\n+\nII\n@b\nGT\n+\nJJ\n".to_vec());
        let r1 = read_fastq_record(&mut c, 1).unwrap().unwrap();
        let r2 = read_fastq_record(&mut c, 5).unwrap().unwrap();
        assert_eq!(r1.header, "@a");
        assert_eq!(r2.header, "@b");
        assert!(read_fastq_record(&mut c, 9).unwrap().is_none());
    }

    #[test]
    fn rejects_missing_at_prefix() {
        let mut c = Cursor::new(b"not-a-header\nACGT\n+\nIIII\n".to_vec());
        assert!(read_fastq_record(&mut c, 1).is_err());
    }

    #[test]
    fn rejects_missing_plus_prefix() {
        let mut c = Cursor::new(b"@r\nACGT\nXXXX\nIIII\n".to_vec());
        assert!(read_fastq_record(&mut c, 1).is_err());
    }

    #[test]
    fn rejects_seq_qual_length_mismatch() {
        let mut c = Cursor::new(b"@r\nACGT\n+\nIII\n".to_vec());
        assert!(read_fastq_record(&mut c, 1).is_err());
    }

    #[test]
    fn rejects_truncated_record() {
        let mut c = Cursor::new(b"@r\nACGT\n".to_vec());
        assert!(read_fastq_record(&mut c, 1).is_err());
    }

    #[test]
    fn base_id_strips_legacy_mate_suffix() {
        let r = FastqRecord {
            header: "@READ_1/1".into(),
            seq: String::new(),
            plus: "+".into(),
            qual: String::new(),
        };
        assert_eq!(r.base_id(), "READ_1");
    }

    #[test]
    fn base_id_strips_illumina_style_suffix() {
        let r1 = FastqRecord {
            header: "@READ_1 1:N:0:ATCG".into(),
            seq: String::new(),
            plus: "+".into(),
            qual: String::new(),
        };
        let r2 = FastqRecord {
            header: "@READ_1 2:N:0:ATCG".into(),
            seq: String::new(),
            plus: "+".into(),
            qual: String::new(),
        };
        assert_eq!(r1.base_id(), "READ_1");
        assert_eq!(r1.base_id(), r2.base_id());
    }

    #[test]
    fn recv_pair_step_matches_when_ids_agree() {
        let (tx1, rx1) = mpsc::sync_channel(1);
        let (tx2, rx2) = mpsc::sync_channel(1);
        tx1.send(Ok(rec("@r 1:N:0:A"))).unwrap();
        tx2.send(Ok(rec("@r 2:N:0:A"))).unwrap();
        match recv_pair_step(&rx1, &rx2) {
            PairStep::Pair { ids_match, .. } => assert!(ids_match),
            _ => panic!("expected Pair"),
        }
    }

    #[test]
    fn recv_pair_step_flags_id_mismatch_without_failing() {
        let (tx1, rx1) = mpsc::sync_channel(1);
        let (tx2, rx2) = mpsc::sync_channel(1);
        tx1.send(Ok(rec("@r1"))).unwrap();
        tx2.send(Ok(rec("@r2"))).unwrap();
        match recv_pair_step(&rx1, &rx2) {
            PairStep::Pair { ids_match, r1, r2 } => {
                assert!(!ids_match);
                assert_eq!(r1.base_id(), "r1");
                assert_eq!(r2.base_id(), "r2");
            }
            _ => panic!("expected Pair with ids_match = false"),
        }
    }

    #[test]
    fn recv_pair_step_detects_synchronized_eof() {
        let (tx1, rx1) = mpsc::sync_channel::<Result<FastqRecord>>(1);
        let (tx2, rx2) = mpsc::sync_channel::<Result<FastqRecord>>(1);
        drop(tx1);
        drop(tx2);
        assert!(matches!(recv_pair_step(&rx1, &rx2), PairStep::Eof));
    }

    #[test]
    fn recv_pair_step_detects_count_mismatch() {
        let (tx1, rx1) = mpsc::sync_channel(1);
        let (tx2, rx2) = mpsc::sync_channel::<Result<FastqRecord>>(1);
        tx1.send(Ok(rec("@r1"))).unwrap();
        drop(tx1);
        drop(tx2);
        assert!(matches!(
            recv_pair_step(&rx1, &rx2),
            PairStep::CountMismatch
        ));
    }

    #[test]
    fn recv_pair_step_surfaces_read_error() {
        let (tx1, rx1) = mpsc::sync_channel::<Result<FastqRecord>>(1);
        let (tx2, rx2) = mpsc::sync_channel(1);
        tx1.send(Err(anyhow::anyhow!("boom"))).unwrap();
        tx2.send(Ok(rec("@r2"))).unwrap();
        assert!(matches!(recv_pair_step(&rx1, &rx2), PairStep::ReadError(_)));
    }
}
