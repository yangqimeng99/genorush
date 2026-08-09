mod common;
mod fastx;
mod gff;
mod io_utils;

use anyhow::Result;
use clap::{Parser, Subcommand};

/// GenoRush - a fast, parallel, cross-platform CLI toolkit for genomics/transcriptomics data
#[derive(Parser, Debug)]
#[command(name = "genorush", version, propagate_version = true)]
struct Cli {
    /// Worker threads for parallel processing (0 = use all logical cores)
    #[arg(short = 'j', long, global = true, default_value_t = 0)]
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
