use std::collections::BTreeSet;

use anyhow::Result;
use config::{
    get_progress_bar, write_objs, HIT, OVERLAP_CDS, OVERLAP_EXON, PASS,
    TRUNCATION_RECOVERY_THRESHOLD, TRUNCATION_THRESHOLD,
};
use dashmap::DashSet;
use hashbrown::HashSet;
use log::{info, warn};
use packbed::{packbed, GenePred, RefGenePred};
use rayon::prelude::*;

use crate::cli::Args;
use crate::utils::unpack_blacklist;

pub fn detect_truncations(args: Args) -> Result<()> {
    info!("Detecting 5'end truncations...");

    let tracks = packbed(args.refs, args.query, OVERLAP_CDS, OVERLAP_EXON)?;
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let hit_acc: DashSet<String> = DashSet::new();
    let pass_acc: DashSet<String> = DashSet::new();
    // let dirt_acc: DashSet<String> = DashSet::new();

    let pb = get_progress_bar(tracks.len() as u64, "Processing...");
    tracks.par_iter().for_each(|bucket| {
        let components = bucket.value().to_owned();

        components.into_par_iter().for_each(|comp| {
            let (hits, pass) = process_component(comp, &blacklist, args.recover);

            hits.into_iter().for_each(|hit| {
                hit_acc.insert(hit);
            });
            pass.into_iter().for_each(|p| {
                pass_acc.insert(p);
            });
        });

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!("Reads with 5'end truncations: {}", hit_acc.len());

    [hit_acc, pass_acc]
        .par_iter()
        .zip([HIT, PASS].par_iter())
        .for_each(|(rx, path)| write_objs(&rx, path));

    Ok(())
}

#[inline(always)]
pub fn process_component(
    comp: (RefGenePred, Vec<GenePred>),
    ban: &HashSet<String>,
    recover: bool,
) -> (Vec<String>, Vec<String>) {
    let mut truncations = Vec::new();
    let mut pass = Vec::new();

    let mut tmp_dirt = Vec::new();
    let mut owners = BTreeSet::new();

    let refs = comp.0;
    let queries = comp.1;

    let ref_starts = &refs.starts;
    let ref_middles = &refs.middles;

    let (mut t, totals) = (0_f32, queries.len() as f32);

    for query in queries.iter() {
        if ban.contains(query.name()) {
            continue;
        }

        let (query_start, query_end) = query.get_first_exon();

        let is_complete = ref_starts.iter().any(|(s, e)| {
            if query_end < *s {
                return false;
            }

            (query_start >= *s) && (query_start < *e)
        });

        if is_complete {
            // still checks if read start is inside any middle boundaries
            if ref_middles.iter().any(|(s, e)| {
                if (query_start >= *s) && (query_start < *e) {
                    owners.insert((s, e));
                    tmp_dirt.push(query);
                    return true;
                } else {
                    return false;
                }
            }) {
                let line = query.line().to_owned();
                truncations.push(line);

                t += 1.0;
            } else {
                let line = query.line().to_owned();
                pass.push(line);
            }
        } else {
            // we do not have any overlap with consensus starts.
            // we need to see if we overlap any middle exons, if so read
            // is truncated, otherwise it is a novel start.
            let is_truncated = ref_middles.iter().any(|(mid_exon_start, mid_exon_end)| {
                if query_end < *mid_exon_start {
                    return false;
                }

                if (query_start >= *mid_exon_start && query_start < *mid_exon_end)
                    || (query_end > *mid_exon_start && query_end <= *mid_exon_end)
                    || (query_start < *mid_exon_start && query_end > *mid_exon_end)
                {
                    owners.insert((mid_exon_start, mid_exon_end));
                    tmp_dirt.push(query);
                    return true;
                } else {
                    return false;
                }
            });

            if is_truncated {
                let line = query.line().to_owned();
                truncations.push(line);

                t += 1.0;
            } else {
                let line = query.line().to_owned();
                pass.push(line);
            }
        }
    }

    // after classying reads, we check bucket frequencies
    // if the number of truncated reads is greater than 50%
    // of the total reads in the bucket, we consider the bucket
    // to be dirty and return the reads for recovery if args.recover
    if recover {
        if (t / totals) >= TRUNCATION_THRESHOLD {
            warn!("Bucket {:?} is dirty -> {}", refs.reads, t / totals);
            let dirt = tmp_dirt;

            let new_passes = recover_from_dirt(dirt, owners, &refs);
            new_passes.iter().for_each(|p| {
                truncations.retain(|x| x != p);
                pass.push(p.to_owned());
            });
        }
    } else {
        drop(tmp_dirt)
    }

    return (truncations, pass);
}

pub fn recover_from_dirt(
    mut dirt: Vec<&GenePred>,
    owners: BTreeSet<(&u64, &u64)>,
    refs: &RefGenePred,
) -> Vec<String> {
    // 1. see how many of each owner we have in refs.reads

    let mut local_passes = vec![];
    let background = refs.reads.len() as f32;

    for owner in owners.iter() {
        let mut count = 0.0;
        for read in refs.reads.iter() {
            let ref_exons = read.get_middle_exons();
            let (s, e) = owner;

            if ref_exons.contains(&(**s, **e)) {
                count += 1.0;
            }

            // ref_exons.iter().for_each(|(start, end)| {
            //     if (start >= *s) && (*e <= end) {
            //          count += 1;
            //     }
            // });
        }

        // 2. if the owner has more than 50% support, we consider it
        //  a valid owner and keep reads truncated by that owner as
        //  truncated reads; otherwise, we consider the owner to be
        //  a weak ownner and send all truncated reads to the pass
        //  bucket
        let ratio = count / background;
        if ratio < TRUNCATION_RECOVERY_THRESHOLD {
            // send reads truncated by this owner to pass
            for read in dirt.clone().iter_mut() {
                let (owner_start, owner_end) = *owner;
                let (query_start, query_end) = read.get_first_exon();

                if (query_start >= *owner_start && query_start < *owner_end)
                    || (query_end > *owner_start && query_end <= *owner_end)
                    || (query_start < *owner_start && query_end > *owner_end)
                {
                    // send to pass and remove from dirt
                    let line = read.line().to_owned();
                    local_passes.push(line);
                    dirt.retain(|x| x.name() != read.name());
                }
            }
        }
    }

    local_passes
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile;

    #[test]
    fn test_detection_fn_with_tempfile() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        let path = file.path().to_path_buf();
        write!(
            file,
            "s1/t778845/t804002/tm54164U_210309_085211/65275776/ccs_PerID0.996_5Clip0_3Clip0_PolyA274_PolyARead275/t60/t-/t778845/t804002/t255,0,0/t14/t1268,142,88,200,253,203,254,167,120,142,218,197,132,302/t0,14396,14736,14943,16461,16846,17598,18073,18511,18890,20330,21290,22627,24855"
        ).unwrap();

        write!(
            file,
            "s1/t778870/t803968/tm54164U_210309_085211/92276372/ccs_PerID1.000_5Clip0_3Clip0_PolyA71_PolyARead72/t60/t-/t778870/t803968/t255,0,0/t11/t1243,142,88,200,253,1006,167,218,197,132,268/t0,14371,14711,14918,16436,16821,18048,20305,21265,22602,24830"
        ).unwrap();

        let args = Args {
            refs: [path.clone()].to_vec(),
            query: [path].to_vec(),
            threads: 1,
            blacklist: Vec::new(),
            recover: false,
            skip_exon: false,
        };

        assert!(detect_truncations(args).is_ok());
    }
}
