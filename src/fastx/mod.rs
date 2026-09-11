pub mod cat;
pub mod deinterleave;
pub mod interleave;
pub mod rename;
pub mod rescue;
pub mod sample;

use anyhow::Result;
use clap::{Args, Subcommand};

use crate::io_utils::OutputOpts;

/// FASTA/FASTQ sequence utilities
#[derive(Args, Debug)]
pub struct FastxCli {
    #[command(subcommand)]
    command: FastxCommand,
}

#[derive(Subcommand, Debug)]
enum FastxCommand {
    /// Rename sequence names in a FASTA file via a mapping table
    ///
    /// -j/--threads: output volume is close to input volume here (a
    /// near-1:1 rewrite), so compression -- not the single-threaded read
    /// side -- is the bottleneck and benefits from most cores you can
    /// spare. Recommended -j 8; more helps further on machines with more
    /// cores, with fast-diminishing returns per added thread past that.
    Rename(rename::RenameArgs),
    /// Downsample FASTQ reads by proportion or exact count (single- or paired-end)
    ///
    /// -j/--threads: write volume (and thus how much -j helps) scales with
    /// how much you keep. Recommended -j 1-2 for typical light sampling
    /// (-p well under 1, or a small -n); -j 8 if keeping most of the input
    /// (-p close to 1, or a large -n). See docs/en/sample.md.
    Sample(sample::SampleArgs),
    /// Recover the leading run of clean reads from a truncated/corrupted FASTQ (single- or paired-end)
    ///
    /// -j/--threads: on typical input (corruption near the end, so most of
    /// the file gets rescued), write volume is close to read volume and
    /// compression dominates. Recommended -j 8; less matters more if
    /// corruption is near the very start (little gets written either way).
    Rescue(rescue::RescueArgs),
    /// Split a merged paired-end FASTQ back into R1/R2, by position or by /1 and /2 header markers
    ///
    /// Handles both positional layouts (interleaved, or R1-then-R2
    /// concatenated) and files where position carries no information at all
    /// and each record's own `/1`/`/2` marker is the only thing that says
    /// which mate it is -- see --layout.
    ///
    /// -j/--threads: nothing is discarded, just split across two files, so
    /// total write volume equals read volume and compression is the
    /// bottleneck. Recommended -j 8.
    Deinterleave(deinterleave::DeinterleaveArgs),
    /// Merge R1/R2 into a single standard interleaved FASTQ
    ///
    /// -j/--threads: write volume equals read volume (nothing discarded)
    /// and two mates' worth of data needs compressing, so this benefits
    /// from threads more than most commands here. Recommended -j 8.
    Interleave(interleave::InterleaveArgs),
    /// Concatenate FASTQ files from repeated sequencing runs, checking for duplicate read IDs
    ///
    /// -j/--threads: a straight concatenation, so write volume equals read
    /// volume and compression is the bottleneck. Recommended -j 8.
    Cat(cat::CatArgs),
}

pub fn run(cli: FastxCli, opts: OutputOpts) -> Result<()> {
    match cli.command {
        FastxCommand::Rename(args) => rename::run(args, opts),
        FastxCommand::Sample(args) => sample::run(args, opts),
        FastxCommand::Rescue(args) => rescue::run(args, opts),
        FastxCommand::Deinterleave(args) => deinterleave::run(args, opts),
        FastxCommand::Interleave(args) => interleave::run(args, opts),
        FastxCommand::Cat(args) => cat::run(args, opts),
    }
}
