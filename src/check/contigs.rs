//! `check contigs`: do these files agree about the genome they describe?
//!
//! The failure this exists for is quiet. A BAM aligned against a reference
//! whose contigs are called `1, 2, 3` and a VCF annotated against one calling
//! them `chr1, chr2, chr3` are both perfectly valid files. Tools that join
//! them do not usually crash -- they find nothing for the contigs they cannot
//! match and report that as zero, which is indistinguishable from a real
//! answer. Lengths differing by a few bases are worse still: everything lines
//! up until it doesn't, somewhere in the middle of a chromosome.
//!
//! GATK reports these as "incompatible contigs" when it happens to look, and
//! the usual remedy is Picard. This is the check on its own, ahead of time,
//! reading headers rather than data.

use std::collections::BTreeSet;
use std::path::PathBuf;

use anyhow::{bail, Result};
use clap::Args;

use crate::common::contigs::{read_contigs, ContigSet};
use crate::io_utils::OutputOpts;

#[derive(Args, Debug)]
pub struct ContigsArgs {
    /// Files to compare. The first is the reference the others are checked
    /// against -- usually the genome FASTA, or its `.fai`.
    #[arg(value_name = "FILE", required = true, num_args = 2..)]
    files: Vec<PathBuf>,

    /// Treat a different contig order as a failure, not a remark. Order
    /// matters to tools that require coordinate-sorted input against a fixed
    /// dictionary; it is irrelevant to most others.
    #[arg(long)]
    require_order: bool,
}

/// Something worth telling the user about, and whether it should stop a
/// pipeline.
struct Finding {
    fatal: bool,
    text: String,
}

pub fn run(args: ContigsArgs, _opts: OutputOpts) -> Result<()> {
    let sets: Vec<ContigSet> = args
        .files
        .iter()
        .map(|p| read_contigs(p))
        .collect::<Result<_>>()?;

    let (reference, others) = sets.split_first().expect("at least two files required");
    print_summary(reference, others);

    let mut findings = Vec::new();
    for set in others {
        findings.extend(compare(reference, set, args.require_order));
    }

    if findings.is_empty() {
        println!("\nAll files agree with {}.", reference.path.display());
        return Ok(());
    }

    println!();
    for f in &findings {
        println!("{}", f.text);
    }

    let fatal = findings.iter().filter(|f| f.fatal).count();
    if fatal > 0 {
        bail!(
            "{fatal} of {} findings would change results silently rather than fail loudly",
            findings.len()
        );
    }
    Ok(())
}

fn print_summary(reference: &ContigSet, others: &[ContigSet]) {
    let width = std::iter::once(reference)
        .chain(others)
        .map(|s| s.path.display().to_string().len())
        .max()
        .unwrap_or(10)
        .max(6);

    println!(
        "{:<width$}  {:<6}  {:>7}  lengths",
        "file", "kind", "contigs"
    );
    for (i, set) in std::iter::once(reference).chain(others).enumerate() {
        let sized = set.contigs.iter().filter(|c| c.length.is_some()).count();
        let lengths = if set.contigs.is_empty() {
            "-".to_string()
        } else if sized == set.contigs.len() {
            "yes".to_string()
        } else if sized == 0 {
            "no".to_string()
        } else {
            format!("{sized}/{}", set.contigs.len())
        };
        println!(
            "{:<width$}  {:<6}  {:>7}  {}{}",
            set.path.display().to_string(),
            set.kind.label(),
            set.contigs.len(),
            lengths,
            if i == 0 { "   (reference)" } else { "" }
        );
    }
}

/// The most common reason a name is missing is a naming convention, not a
/// missing sequence. Saying which convention turns "not found" into a fix.
fn alias_hint(name: &str, reference: &ContigSet) -> Option<String> {
    let candidates: Vec<String> = match name.strip_prefix("chr") {
        Some(bare) => vec![bare.to_string()],
        None => vec![format!("chr{name}")],
    };
    let mito = match name {
        "chrM" | "chrMT" => vec!["MT".to_string(), "M".to_string()],
        "MT" | "M" => vec!["chrM".to_string(), "chrMT".to_string()],
        _ => Vec::new(),
    };
    candidates
        .into_iter()
        .chain(mito)
        .find(|c| reference.get(c).is_some())
        .map(|c| format!("the reference calls it {c:?}"))
}

fn compare(reference: &ContigSet, other: &ContigSet, require_order: bool) -> Vec<Finding> {
    let mut findings = Vec::new();
    let label = other.path.display().to_string();

    // A name the reference has never heard of is the loud version of the
    // silent failure: whatever consumes both will simply find nothing there.
    let mut unknown = Vec::new();
    for contig in &other.contigs {
        if reference.get(&contig.name).is_none() {
            let hint = alias_hint(&contig.name, reference)
                .map(|h| format!(" -- {h}"))
                .unwrap_or_default();
            unknown.push(format!("      {}{}", contig.name, hint));
        }
    }
    if !unknown.is_empty() {
        findings.push(Finding {
            fatal: true,
            text: format!(
                "  {label}: {} contig(s) the reference does not have\n{}",
                unknown.len(),
                unknown.join("\n")
            ),
        });
    }

    // Two files that both state a length and disagree are describing
    // different genome builds, whatever their names suggest.
    let mut conflicts = Vec::new();
    for contig in &other.contigs {
        let (Some(theirs), Some(ours)) = (
            contig.length,
            reference.get(&contig.name).and_then(|c| c.length),
        ) else {
            continue;
        };
        if theirs != ours {
            conflicts.push(format!(
                "      {}: {theirs} here, {ours} in the reference",
                contig.name
            ));
        }
    }
    if !conflicts.is_empty() {
        findings.push(Finding {
            fatal: true,
            text: format!(
                "  {label}: {} contig(s) with a different length -- these are different \
                 genome builds\n{}",
                conflicts.len(),
                conflicts.join("\n")
            ),
        });
    }

    // Order only matters to some tools, so it is a remark unless asked for.
    let shared_here: Vec<&str> = other
        .names()
        .filter(|n| reference.get(n).is_some())
        .collect();
    let shared_there: Vec<&str> = reference
        .names()
        .filter(|n| other.get(n).is_some())
        .collect();
    if shared_here != shared_there {
        findings.push(Finding {
            fatal: require_order,
            text: format!(
                "  {label}: lists its contigs in a different order than the reference{}",
                if require_order {
                    ""
                } else {
                    " (pass --require-order if that matters here)"
                }
            ),
        });
    }

    // Contigs the reference has and this file does not mention are normal --
    // a BED only names what it covers -- so they are counted, not flagged.
    // Unless the file also used names the reference lacks: then this is not a
    // partial file, it is the same naming problem seen from the other side,
    // and calling it normal would be actively misleading.
    let names_disagree = findings.iter().any(|f| f.fatal);
    let missing: BTreeSet<&str> = reference
        .names()
        .filter(|n| other.get(n).is_none())
        .collect();
    if !missing.is_empty() && !names_disagree {
        findings.push(Finding {
            fatal: false,
            text: format!(
                "  {label}: does not mention {} of the reference's contigs (normal for a \
                 file that only covers part of the genome)",
                missing.len()
            ),
        });
    }

    findings
}
