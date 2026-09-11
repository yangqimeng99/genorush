use std::collections::HashMap;

use anyhow::Result;
use clap::Args;

use crate::common::rename::{self, RenameCommonArgs};
use crate::io_utils::OutputOpts;

#[derive(Args, Debug)]
pub struct RenameArgs {
    #[command(flatten)]
    common: RenameCommonArgs,
}

fn fasta_line(line: &str, dict: &HashMap<String, String>, out: &mut Vec<u8>) {
    match line.strip_prefix('>') {
        Some(rest) => {
            // A header is rewritten either way: mapped or not, the original
            // script keeps only the first token, so the description goes.
            let old_name = rest.split_whitespace().next().unwrap_or("");
            let name = dict.get(old_name).map_or(old_name, |n| n.as_str());
            out.push(b'>');
            out.extend_from_slice(name.as_bytes());
        }
        // Sequence lines -- nearly the whole file -- are copied straight
        // through.
        None => out.extend_from_slice(line.as_bytes()),
    }
}

pub fn run(args: RenameArgs, opts: OutputOpts) -> Result<()> {
    rename::run(&args.common, fasta_line, opts)
}
