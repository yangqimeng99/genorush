pub mod rename;
pub mod sample;

use anyhow::Result;
use clap::{Args, Subcommand};

/// FASTA/FASTQ sequence utilities
#[derive(Args, Debug)]
pub struct FastxCli {
    #[command(subcommand)]
    command: FastxCommand,
}

#[derive(Subcommand, Debug)]
enum FastxCommand {
    /// Rename sequence names in a FASTA file via a mapping table
    Rename(rename::RenameArgs),
    /// Downsample FASTQ reads by proportion or exact count (single- or paired-end)
    Sample(sample::SampleArgs),
}

pub fn run(cli: FastxCli) -> Result<()> {
    match cli.command {
        FastxCommand::Rename(args) => rename::run(args),
        FastxCommand::Sample(args) => sample::run(args),
    }
}
