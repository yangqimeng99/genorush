mod common;
mod fastx;
mod gff;
mod io_utils;

use anyhow::Result;
use clap::{Parser, Subcommand};

/// GenoRush - a fast, parallel, cross-platform CLI toolkit for bioinformatics data
#[derive(Parser, Debug)]
#[command(name = "genorush", version, propagate_version = true)]
struct Cli {
    /// Worker threads for parallel processing (default: 1; pass 0 to use all
    /// logical cores). Only speeds up compressing output, not decompressing
    /// input (a single gzip stream can't be decompressed in parallel), so
    /// the useful ceiling depends on the command and how much it writes
    /// relative to what it reads -- see `<category> <action> --help` for a
    /// per-command recommendation.
    #[arg(short = 'j', long, global = true, default_value_t = 1)]
    threads: usize,

    #[command(subcommand)]
    command: TopCommand,
}

#[derive(Subcommand, Debug)]
enum TopCommand {
    Fastx(fastx::FastxCli),
    Gff(gff::GffCli),
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .expect("failed to initialize thread pool");
    }
    log::info!("using {} worker thread(s)", rayon::current_num_threads());

    match cli.command {
        TopCommand::Fastx(c) => fastx::run(c),
        TopCommand::Gff(c) => gff::run(c),
    }
}
