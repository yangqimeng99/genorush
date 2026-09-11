//! Reading the contig list out of whatever genomics file is at hand.
//!
//! Almost every format carries some idea of "the sequences this file is about":
//! a FASTA index has names and lengths, a VCF header has `##contig` lines, a
//! BAM header has `@SQ` records, a BED has names in its first column and
//! nothing else. They agree until they don't, and when they don't the symptom
//! is rarely an error -- a tool that cannot find a contig usually just reports
//! nothing for it.
//!
//! This module turns all of them into one shape so they can be compared. It is
//! also the place where this project first reaches for `noodles`: VCF, BCF,
//! SAM, BAM and CRAM headers are real parsing problems, and BCF/BAM/CRAM are
//! binary formats with their own compression framing. GFF, BED and FASTA
//! indexes are not -- their contig names are the first column of a text file --
//! so they are read here directly. Taking a version-coupled parser for a
//! `split('\t').next()` would be a cost without a benefit; when a command needs
//! to *manipulate* GFF rather than glance at it, that is the time to add one.

use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::io_utils::open_reader;

/// One sequence, as some file described it. `length` is `None` for formats
/// that name contigs without sizing them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub length: Option<u64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// `.fai`, `.chrom.sizes`, `.genome` -- name and length, two columns.
    Sizes,
    Fasta,
    Vcf,
    Bcf,
    Sam,
    Bam,
    Cram,
    Gff,
    Bed,
}

impl Kind {
    pub fn label(self) -> &'static str {
        match self {
            Kind::Sizes => "sizes",
            Kind::Fasta => "FASTA",
            Kind::Vcf => "VCF",
            Kind::Bcf => "BCF",
            Kind::Sam => "SAM",
            Kind::Bam => "BAM",
            Kind::Cram => "CRAM",
            Kind::Gff => "GFF",
            Kind::Bed => "BED",
        }
    }
}

/// The contigs one file describes, in the order that file lists them.
pub struct ContigSet {
    pub path: PathBuf,
    pub kind: Kind,
    pub contigs: Vec<Contig>,
}

impl ContigSet {
    pub fn get(&self, name: &str) -> Option<&Contig> {
        self.contigs.iter().find(|c| c.name == name)
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.contigs.iter().map(|c| c.name.as_str())
    }
}

/// Picks a format from the file name, ignoring a trailing `.gz`.
///
/// Content sniffing would be more robust for the binary formats, but the text
/// ones are genuinely ambiguous -- a `.gff` and a `.bed` differ only by column
/// count -- so a name that says what it is beats a guess that might not.
pub fn detect_kind(path: &Path) -> Result<Kind> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    let stem = name.strip_suffix(".gz").unwrap_or(&name);
    let stem = stem.strip_suffix(".bgz").unwrap_or(stem);

    let kind = if stem.ends_with(".fai")
        || stem.ends_with(".chrom.sizes")
        || stem.ends_with(".genome")
        || stem.ends_with(".sizes")
    {
        Kind::Sizes
    } else if [".fa", ".fasta", ".fna", ".fas"]
        .iter()
        .any(|e| stem.ends_with(e))
    {
        Kind::Fasta
    } else if stem.ends_with(".vcf") {
        Kind::Vcf
    } else if stem.ends_with(".bcf") {
        Kind::Bcf
    } else if stem.ends_with(".sam") {
        Kind::Sam
    } else if stem.ends_with(".bam") {
        Kind::Bam
    } else if stem.ends_with(".cram") {
        Kind::Cram
    } else if [".gff", ".gff3", ".gtf"].iter().any(|e| stem.ends_with(e)) {
        Kind::Gff
    } else if stem.ends_with(".bed") {
        Kind::Bed
    } else {
        bail!(
            "cannot tell what kind of file {} is from its name. Recognised: \
             .fai/.chrom.sizes/.genome, .fa/.fasta/.fna, .vcf, .bcf, .sam, .bam, \
             .cram, .gff/.gff3/.gtf, .bed -- each optionally .gz",
            path.display()
        );
    };
    Ok(kind)
}

pub fn read_contigs(path: &Path) -> Result<ContigSet> {
    let kind = detect_kind(path)?;
    let contigs = match kind {
        Kind::Sizes => read_sizes(path)?,
        Kind::Fasta => read_fasta(path)?,
        Kind::Vcf => read_vcf(path)?,
        Kind::Bcf => read_bcf(path)?,
        Kind::Sam => read_sam(path)?,
        Kind::Bam => read_bam(path)?,
        Kind::Cram => read_cram(path)?,
        Kind::Gff => read_gff(path)?,
        Kind::Bed => read_bed(path)?,
    };
    Ok(ContigSet {
        path: path.to_path_buf(),
        kind,
        contigs,
    })
}

/// `.fai` and `.chrom.sizes` are the same two leading columns; a `.fai` simply
/// has three more that nothing here needs.
fn read_sizes(path: &Path) -> Result<Vec<Contig>> {
    let reader = open_reader(path)?;
    let mut out = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed reading {}", path.display()))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split_whitespace();
        let name = fields
            .next()
            .with_context(|| format!("{}: line {} has no name", path.display(), i + 1))?;
        let length = fields.next().and_then(|v| v.parse::<u64>().ok());
        out.push(Contig {
            name: name.to_string(),
            length,
        });
    }
    Ok(out)
}

/// Reads names and lengths from the FASTA itself, unless an index sits beside
/// it -- which is the whole point of an index, and turns a full pass over a
/// multi-gigabyte genome into reading a few kilobytes.
fn read_fasta(path: &Path) -> Result<Vec<Contig>> {
    let fai = PathBuf::from(format!("{}.fai", path.display()));
    if fai.is_file() {
        log::info!("using the index beside {}", path.display());
        return read_sizes(&fai);
    }

    let reader = open_reader(path)?;
    let mut out: Vec<Contig> = Vec::new();
    let mut length = 0u64;
    for line in reader.lines() {
        let line = line.with_context(|| format!("failed reading {}", path.display()))?;
        let line = line.trim();
        if let Some(header) = line.strip_prefix('>') {
            if let Some(last) = out.last_mut() {
                last.length = Some(length);
            }
            length = 0;
            out.push(Contig {
                name: header.split_whitespace().next().unwrap_or("").to_string(),
                length: None,
            });
        } else {
            length += line.len() as u64;
        }
    }
    if let Some(last) = out.last_mut() {
        last.length = Some(length);
    }
    Ok(out)
}

fn read_vcf(path: &Path) -> Result<Vec<Contig>> {
    let mut reader = noodles_vcf::io::Reader::new(open_reader(path)?);
    let header = reader
        .read_header()
        .with_context(|| format!("failed reading the VCF header of {}", path.display()))?;
    Ok(from_vcf_header(&header))
}

fn from_vcf_header(header: &noodles_vcf::Header) -> Vec<Contig> {
    header
        .contigs()
        .iter()
        .map(|(name, contig)| Contig {
            name: name.to_string(),
            length: contig.length().map(|l| l as u64),
        })
        .collect()
}

fn from_sam_header(header: &noodles_sam::Header) -> Vec<Contig> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, seq)| Contig {
            name: String::from_utf8_lossy(name.as_ref()).into_owned(),
            length: Some(usize::from(seq.length()) as u64),
        })
        .collect()
}

/// BCF, BAM and CRAM carry their own compression framing, so they are opened
/// as files rather than through the gzip-sniffing reader every text format
/// uses: handing an already-decompressed stream to a reader that expects to do
/// its own decompression gets nowhere.
fn open_binary(path: &Path) -> Result<File> {
    File::open(path).with_context(|| format!("failed to open {}", path.display()))
}

fn read_bcf(path: &Path) -> Result<Vec<Contig>> {
    let mut reader = noodles_bcf::io::Reader::new(open_binary(path)?);
    let header = reader
        .read_header()
        .with_context(|| format!("failed reading the BCF header of {}", path.display()))?;
    Ok(from_vcf_header(&header))
}

fn read_sam(path: &Path) -> Result<Vec<Contig>> {
    let mut reader = noodles_sam::io::Reader::new(open_reader(path)?);
    let header = reader
        .read_header()
        .with_context(|| format!("failed reading the SAM header of {}", path.display()))?;
    Ok(from_sam_header(&header))
}

fn read_bam(path: &Path) -> Result<Vec<Contig>> {
    let mut reader = noodles_bam::io::Reader::new(open_binary(path)?);
    let header = reader
        .read_header()
        .with_context(|| format!("failed reading the BAM header of {}", path.display()))?;
    Ok(from_sam_header(&header))
}

fn read_cram(path: &Path) -> Result<Vec<Contig>> {
    // `read_header` expects the stream at the very start: it checks the magic
    // number and steps over the file definition itself. Reading the definition
    // first, as the sibling method invites, leaves it looking at the wrong
    // bytes and reports an invalid header.
    let mut reader = noodles_cram::io::Reader::new(open_binary(path)?);
    let header = reader
        .read_header()
        .with_context(|| format!("failed reading the CRAM header of {}", path.display()))?;
    Ok(from_sam_header(&header))
}

/// A GFF states its contigs properly only if it carries `##sequence-region`
/// pragmas. Most do not, so the fallback is the first column of the data --
/// which finds every contig the file actually refers to, which is the question
/// being asked.
fn read_gff(path: &Path) -> Result<Vec<Contig>> {
    let reader = open_reader(path)?;
    let mut out: Vec<Contig> = Vec::new();
    let mut seen_pragma = false;
    for line in reader.lines() {
        let line = line.with_context(|| format!("failed reading {}", path.display()))?;
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("##sequence-region") {
            let mut f = rest.split_whitespace();
            if let Some(name) = f.next() {
                let start = f.next().and_then(|v| v.parse::<u64>().ok());
                let end = f.next().and_then(|v| v.parse::<u64>().ok());
                if !seen_pragma {
                    out.clear();
                    seen_pragma = true;
                }
                out.push(Contig {
                    name: name.to_string(),
                    // The pragma gives a span, which is the length when it
                    // starts at 1 -- and it always does in practice.
                    length: match (start, end) {
                        (Some(s), Some(e)) if e >= s => Some(e - s + 1),
                        _ => None,
                    },
                });
            }
            continue;
        }
        if seen_pragma || line.is_empty() || line.starts_with('#') {
            continue;
        }
        let name = line.split('\t').next().unwrap_or("");
        if !name.is_empty() && !out.iter().any(|c| c.name == name) {
            out.push(Contig {
                name: name.to_string(),
                length: None,
            });
        }
    }
    Ok(out)
}

fn read_bed(path: &Path) -> Result<Vec<Contig>> {
    let reader = open_reader(path)?;
    let mut out: Vec<Contig> = Vec::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("failed reading {}", path.display()))?;
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let name = line.split('\t').next().unwrap_or("");
        if !name.is_empty() && !out.iter().any(|c| c.name == name) {
            out.push(Contig {
                name: name.to_string(),
                length: None,
            });
        }
    }
    Ok(out)
}
