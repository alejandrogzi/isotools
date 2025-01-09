use std::collections::{BTreeMap, BTreeSet};

use anyhow::Result;
use packbed::packbed;

use config::{Strand, OVERLAP_CDS, OVERLAP_EXON};

use crate::cli::IntronArgs as Args;

pub struct IntronPred {
    pub chrom: String,
    pub strand: Strand,
    pub introns: BTreeMap<(u64, u64), IntronPredDescriptor>,
}

pub fn classify_introns(args: Args) -> Result<()> {
    let isoseqs = packbed(args.iso, None, OVERLAP_CDS, OVERLAP_EXON, PackMode::Intron)?;
    // we need to convert isoseqs into a DashMap<String, Vec<IntronPred>>

    Ok(())
}
