pub mod cat;
pub mod deinterleave;
pub mod interleave;
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
    /// Split a merged paired-end FASTQ (interleaved or R1-then-R2 concatenated) back into R1/R2
    Deinterleave(deinterleave::DeinterleaveArgs),
    /// Merge R1/R2 into a single standard interleaved FASTQ
    Interleave(interleave::InterleaveArgs),
    /// Concatenate FASTQ files from repeated sequencing runs, checking for duplicate read IDs
    Cat(cat::CatArgs),
}

pub fn run(cli: FastxCli) -> Result<()> {
    match cli.command {
        FastxCommand::Rename(args) => rename::run(args),
        FastxCommand::Sample(args) => sample::run(args),
        FastxCommand::Rescue(args) => rescue::run(args),
        FastxCommand::Deinterleave(args) => deinterleave::run(args),
        FastxCommand::Interleave(args) => interleave::run(args),
        FastxCommand::Cat(args) => cat::run(args),
    }
}
