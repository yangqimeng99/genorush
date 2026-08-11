pub mod rename;

use anyhow::Result;
use clap::{Args, Subcommand};

/// GFF/GTF annotation utilities
#[derive(Args, Debug)]
pub struct GffCli {
    #[command(subcommand)]
    command: GffCommand,
}

#[derive(Subcommand, Debug)]
enum GffCommand {
    /// Rename seqid (column 1) in a GFF/GTF file via a mapping table
    ///
    /// -j/--threads: output volume is close to input volume here (a
    /// near-1:1 rewrite), so compression -- not the single-threaded read
    /// side -- is the bottleneck and benefits from most cores you can
    /// spare. Recommended -j 8; more helps further on machines with more
    /// cores, with fast-diminishing returns per added thread past that.
    Rename(rename::RenameArgs),
}

pub fn run(cli: GffCli) -> Result<()> {
    match cli.command {
        GffCommand::Rename(args) => rename::run(args),
    }
}
