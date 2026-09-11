mod common;
mod fastx;
mod gff;
mod io_utils;

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

    /// Gzip-compress every output, whatever the output path looks like.
    ///
    /// Outputs are normally compressed when their path ends in `.gz`, which
    /// stdout has no way to express: `-o -` is a pipe, not a filename. This
    /// flag is how a pipeline asks for compressed output. It is a no-op for
    /// a path that already ends in `.gz`.
    #[arg(short = 'z', long, global = true)]
    gzip: bool,

    #[command(subcommand)]
    command: TopCommand,
}

#[derive(Subcommand, Debug)]
enum TopCommand {
    Fastx(fastx::FastxCli),
    Gff(gff::GffCli),
}

/// True if the failure is really "the process on the other end of our stdout
/// stopped reading" -- `genorush ... -o - | head` being the everyday case.
///
/// Rust ignores SIGPIPE, so instead of dying quietly like a C tool would,
/// the write comes back as `ErrorKind::BrokenPipe` and would otherwise be
/// reported as a command failure: an error message the user did not cause
/// and a non-zero status that trips `set -o pipefail` in any wrapping
/// script.
fn is_broken_pipe(err: &anyhow::Error) -> bool {
    err.chain().any(|cause| {
        cause
            .downcast_ref::<std::io::Error>()
            .is_some_and(|e| e.kind() == std::io::ErrorKind::BrokenPipe)
    })
}

fn main() -> std::process::ExitCode {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .expect("failed to initialize thread pool");
    }
    log::info!("using {} worker thread(s)", rayon::current_num_threads());

    let opts = io_utils::OutputOpts { gzip: cli.gzip };

    let result = match cli.command {
        TopCommand::Fastx(c) => fastx::run(c, opts),
        TopCommand::Gff(c) => gff::run(c, opts),
    };

    match result {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) if is_broken_pipe(&e) => {
            log::info!("output pipe closed by the reader; stopping early");
            std::process::ExitCode::SUCCESS
        }
        Err(e) => {
            // Matches what `fn main() -> Result<()>` used to print, so error
            // output is unchanged for every failure that isn't a broken pipe.
            eprintln!("Error: {e:?}");
            std::process::ExitCode::FAILURE
        }
    }
}
