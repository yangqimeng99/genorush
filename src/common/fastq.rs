//! Minimal FASTQ record model: read/write a 4-line record verbatim, with just
//! enough structural validation to catch truncated or corrupt input early.

use std::io::{BufRead, Write};

use anyhow::{bail, Result};

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
        bail!("truncated FASTQ record at line {}: missing sequence line", line_no + 1);
    }
    let seq = seq.trim_end_matches(['\n', '\r']).to_string();

    let mut plus = String::new();
    if reader.read_line(&mut plus)? == 0 {
        bail!("truncated FASTQ record at line {}: missing '+' line", line_no + 2);
    }
    let plus = plus.trim_end_matches(['\n', '\r']).to_string();
    if !plus.starts_with('+') {
        bail!("malformed FASTQ at line {}: expected '+' line, got: {plus:?}", line_no + 2);
    }

    let mut qual = String::new();
    if reader.read_line(&mut qual)? == 0 {
        bail!("truncated FASTQ record at line {}: missing quality line", line_no + 3);
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

    Ok(Some(FastqRecord { header, seq, plus, qual }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

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
}
