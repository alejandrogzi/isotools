use anyhow::Result;
use config::{COLORIZE, OVERLAP};
// use log::error;
use packbed::packbed;
// use rayon::prelude::*;

use crate::cli::Args;

pub fn detect_intron_retentions(args: Args) -> Result<()> {
    let tracks = packbed(args.refs, args.query, OVERLAP, COLORIZE)?;

    Ok(())
}
