use std::collections::HashMap;

use anyhow::Result;
use clap::Args;

use crate::common::rename::{self, RenameCommonArgs};

#[derive(Args, Debug)]
pub struct RenameArgs {
    #[command(flatten)]
    common: RenameCommonArgs,
}

fn fasta_line(line: &str, dict: &HashMap<String, String>) -> String {
    match line.strip_prefix('>') {
        Some(rest) => {
            let old_name = rest.split_whitespace().next().unwrap_or("");
            match dict.get(old_name) {
                Some(new_name) => format!(">{new_name}"),
                None => format!(">{old_name}"),
            }
        }
        None => line.to_string(),
    }
}

pub fn run(args: RenameArgs) -> Result<()> {
    rename::run(&args.common, fasta_line)
}
