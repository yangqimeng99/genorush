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

fn gff_line(line: &str, dict: &HashMap<String, String>, out: &mut Vec<u8>) {
    if line.starts_with('#') {
        out.extend_from_slice(line.as_bytes());
        return;
    }
    let mut fields = line.splitn(2, '\t');
    let seqid = fields.next().unwrap_or("");
    match dict.get(seqid) {
        // Always emit the tab even when `rest` is empty: matches the
        // original Python script's `'\t'.join(...)` behavior exactly, so
        // byte-for-byte output parity holds on malformed/short lines too.
        Some(new_name) => {
            out.extend_from_slice(new_name.as_bytes());
            out.push(b'\t');
            out.extend_from_slice(fields.next().unwrap_or("").as_bytes());
        }
        // A seqid that is not in the mapping leaves the line untouched.
        None => out.extend_from_slice(line.as_bytes()),
    }
}

pub fn run(args: RenameArgs, opts: OutputOpts) -> Result<()> {
    rename::run(&args.common, gff_line, opts)
}
