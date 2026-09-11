pub mod contigs;

use anyhow::Result;
use clap::{Args, Subcommand};

use crate::io_utils::OutputOpts;

/// Consistency checks across files that are supposed to describe the same data
#[derive(Args, Debug)]
pub struct CheckCli {
    #[command(subcommand)]
    command: CheckCommand,
}

#[derive(Subcommand, Debug)]
enum CheckCommand {
    /// Compare the contigs several files describe, against the first of them
    ///
    /// Reads only headers and indexes where it can, so this is fast even on
    /// whole genomes and alignments: a FASTA is read in full only when no
    /// `.fai` sits beside it.
    Contigs(contigs::ContigsArgs),
}

pub fn run(cli: CheckCli, opts: OutputOpts) -> Result<()> {
    match cli.command {
        CheckCommand::Contigs(args) => contigs::run(args, opts),
    }
}
