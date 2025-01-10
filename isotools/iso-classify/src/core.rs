use std::collections::{BTreeMap, BTreeSet};

use anyhow::Result;
use packbed::{packbed, PackMode};

use config::{OVERLAP_CDS, OVERLAP_EXON};

use crate::cli::IntronArgs as Args;

pub fn classify_introns(args: Args) -> Result<()> {
    let isoseqs = packbed(
        args.iso,
        args.toga,
        OVERLAP_CDS,
        OVERLAP_EXON,
        PackMode::Intron,
    )?;

    // dbg!(isoseqs);

    Ok(())
}
