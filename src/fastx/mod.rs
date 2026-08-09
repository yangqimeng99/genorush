pub mod rename;
pub mod rescue;
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
    /// Recover the leading run of clean reads from a truncated/corrupted FASTQ (single- or paired-end)
    Rescue(rescue::RescueArgs),
}

pub fn run(cli: FastxCli) -> Result<()> {
    match cli.command {
        FastxCommand::Rename(args) => rename::run(args),
        FastxCommand::Sample(args) => sample::run(args),
        FastxCommand::Rescue(args) => rescue::run(args),
    }
}
