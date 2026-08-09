use std::collections::HashMap;

use anyhow::Result;
use clap::Args;

use crate::common::rename::{self, RenameCommonArgs};

#[derive(Args, Debug)]
pub struct RenameArgs {
    #[command(flatten)]
    common: RenameCommonArgs,
}

fn gff_line(line: &str, dict: &HashMap<String, String>) -> String {
    if line.starts_with('#') {
        return line.to_string();
    }
    let mut fields = line.splitn(2, '\t');
    let seqid = fields.next().unwrap_or("");
    match dict.get(seqid) {
        // Always emit the tab even when `rest` is empty: matches the
        // original Python script's `'\t'.join(...)` behavior exactly, so
        // byte-for-byte output parity holds on malformed/short lines too.
        Some(new_name) => {
            let rest = fields.next().unwrap_or("");
            format!("{new_name}\t{rest}")
        }
        None => line.to_string(),
    }
}

pub fn run(args: RenameArgs) -> Result<()> {
    rename::run(&args.common, gff_line)
}
