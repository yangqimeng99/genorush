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
    Rename(rename::RenameArgs),
}

pub fn run(cli: GffCli) -> Result<()> {
    match cli.command {
        GffCommand::Rename(args) => rename::run(args),
    }
}
