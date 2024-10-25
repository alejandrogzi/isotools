use anyhow::Result;
use config::{COLORIZE, OVERLAP};
use hashbrown::HashSet;
use packbed::{packbed, GenePred};
use rayon::prelude::*;
use std::sync::Arc;
// use log::error;

use crate::cli::Args;
use crate::utils::{unpack_blacklist, Bed4};

pub fn detect_intron_retentions(args: Args) -> Result<()> {
    let tracks = packbed(args.refs, args.query, OVERLAP, COLORIZE)?;
    let blacklist = unpack_blacklist(args.blacklist)?;

    tracks.par_iter().for_each(|(chr, buckets)| {
        let binding = HashSet::new();
        let banned = blacklist.get(chr).unwrap_or(&binding);
        buckets.par_iter().for_each(|bucket| {
            let rs = process_bucket(bucket, banned);
        });
    });

    Ok(())
}

/// Process a bucket of overlapping transcripts
///
/// # Arguments
///     bucket: &Vec<Arc<GenePred>> - list of overlapping transcripts
///     ban: &HashSet<(u64, u64)> - list of blacklisted introns
///
/// # Returns
///    Vec<Arc<GenePred>> - list of retained introns
///
/// # Example
/// ```rust
/// let rs = process_bucket(bucket, banned);
/// ```
///
/// # Description
///
/// This function processes a bucket of overlapping transcripts to
/// identify retained introns. In brief, loops over all query transcripts
/// in the bucket and check if any of their exons completely overlap with
/// any reference intron. If so, the transcript is considered to have a
/// retained intron.
fn process_bucket(bucket: &Vec<Arc<GenePred>>, ban: &HashSet<(u64, u64)>) {
    for query in bucket {
        if query.is_ref() {
            continue;
        }

        for ref_ in bucket {
            if !ref_.is_ref() {
                continue;
            }

            let exons = &query.exons;
            let introns = &ref_.introns;
        }
    }
}
