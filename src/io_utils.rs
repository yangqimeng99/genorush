use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Opens `path` for reading, transparently decompressing gzip/bgzip input.
/// Detection is by magic bytes rather than file extension, so gzip data
/// piped through a renamed file still works. `MultiGzDecoder` is required
/// (not `GzDecoder`) because bgzip-compressed genome references are valid
/// concatenated multi-member gzip streams.
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let mut file = File::open(path)
        .with_context(|| format!("failed to open input file: {}", path.display()))?;

    let mut magic = [0u8; 2];
    let read_n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;

    if read_n == 2 && magic == GZIP_MAGIC {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Opens `path` for writing. If the path ends in `.gz`, output is gzip
/// compressed on the fly; otherwise it is written as plain text.
pub fn open_writer(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path)
        .with_context(|| format!("failed to create output file: {}", path.display()))?;
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    if is_gz {
        Ok(Box::new(GzEncoder::new(
            BufWriter::new(file),
            Compression::default(),
        )))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

/// Reads at most `max_lines` lines from `reader` into `out`, trimming
/// leading/trailing whitespace from each line (mirrors Python's `str.strip()`
/// semantics used by the original script). Returns the number of lines read;
/// 0 means EOF.
pub fn read_line_chunk(
    reader: &mut dyn BufRead,
    out: &mut Vec<String>,
    max_lines: usize,
) -> Result<usize> {
    out.clear();
    let mut buf = String::new();
    let mut n_read = 0;
    while n_read < max_lines {
        buf.clear();
        let n = reader.read_line(&mut buf).context("failed reading line")?;
        if n == 0 {
            break;
        }
        out.push(buf.trim().to_string());
        n_read += 1;
    }
    Ok(n_read)
}

pub fn write_lines(writer: &mut dyn Write, lines: &[String]) -> io::Result<()> {
    for line in lines {
        writer.write_all(line.as_bytes())?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
