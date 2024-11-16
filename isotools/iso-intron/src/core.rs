use anyhow::Result;
use config::{get_progress_bar, write_objs, BED3, HIT, OVERLAP_CDS, OVERLAP_EXON, PASS, SCALE};
use dashmap::DashSet;
use hashbrown::HashSet;
use log::info;
use packbed::{packbed, GenePred, RefGenePred};
use rayon::prelude::*;

use crate::cli::Args;
use crate::utils::{unpack_blacklist, Bed4};

/// Detects intron retentions in an query set of reads
/// based on a provided reference set.
///
/// # Example
/// ```ignore
/// use iso_intron::cli::Args;
/// use iso_intron::core::detect_intron_retentions;
///
/// let args = Args {
///   refs: vec![],
///   query: vec![],
///   blacklist: vec![],
///   plot: true,
///   threads: 4,
/// };
///
/// detect_intron_retentions(args);
/// ```
///
/// # Description
///
/// Entry point for detecting intron retentions in a query set of reads.
/// The function packs the reference and query transcripts into a collection
/// and processes each bucket of overlapping transcripts to identify
/// retained introns. Here we interpret intron retentions as the complete
/// exon overlap with a reference intron:
///
/// ```text
///     5'                          3'
///     XXXXX-----XXXXXX--------XXXXXX
///                     ^^^^^^^^
///     XXXXX-----XXXXXXXXXXXXXXXXXXXX
/// ```
pub fn detect_intron_retentions(args: Args) -> Result<()> {
    info!("Detecting intron retentions...");

    let tracks = packbed(args.refs, args.query, OVERLAP_CDS, OVERLAP_EXON)?;
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let hit_acc: DashSet<String> = DashSet::new();
    let pass_acc: DashSet<String> = DashSet::new();
    let misc_acc: DashSet<String> = DashSet::new();

    let pb = get_progress_bar(tracks.len() as u64, "Processing...");
    tracks.par_iter().for_each(|bucket| {
        let chr = bucket.key();
        let components = bucket.value().to_owned();

        let binding = HashSet::new();
        let banned = blacklist.get(chr).unwrap_or(&binding);

        components.into_par_iter().for_each(|comp| {
            let (hits, pass, blocks) = process_component(comp, banned, args.plot);

            hits.into_iter().for_each(|hit| {
                hit_acc.insert(hit);
            });
            pass.into_iter().for_each(|p| {
                pass_acc.insert(p);
            });
            if let Some(b) = blocks {
                misc_acc.insert(b);
            }
        });

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!("Reads with retained introns: {}", hit_acc.len());

    [hit_acc, pass_acc]
        .par_iter()
        .zip([HIT, PASS].par_iter())
        .for_each(|(rx, path)| write_objs(&rx, path));

    if args.plot {
        write_objs(&misc_acc, BED3);
    }

    Ok(())
}

/// Process a bucket of overlapping transcripts
///
/// # Example
/// ```rust
/// use iso_intron::core::process_bucket;
/// use packbed::Bed12;
/// use hashbrown::HashSet;
/// use std::sync::Arc;
///
/// let line = "chr1\t100\t200\ttest\t0\t+\t100\t200\t0\t1\t1,1\t0,100";
/// let bucket = vec![Bed12::parse(line, false, false).unwrap().into()];
/// let banned = HashSet::from_iter(vec![(100, 200)]);
/// let plot = false;
///
/// let rs = process_bucket(&bucket, &banned, plot);
/// ```
///
/// # Description
///
/// This function processes a bucket of overlapping transcripts to
/// identify retained introns. In brief, loops over all consensus query
/// transcripts in the bucket and check if any of their exons completely
/// overlap with any reference intron. If so, the transcript is considered
/// to have a retained intron.
#[inline(always)]
pub fn process_component(
    comp: (RefGenePred, Vec<GenePred>),
    ban: &HashSet<(u64, u64)>,
    plot: bool,
) -> (Vec<String>, Vec<String>, Option<String>) {
    let mut hits = Vec::new();
    let mut pass = Vec::new();
    let mut blocks = if plot { Some(String::new()) } else { None };

    let refs = comp.0;
    let queries = comp.1;

    for query in queries {
        let mut hit: bool = false;

        for exon in &query.exons {
            for intron in &refs.introns {
                if exon.1 < intron.0 {
                    break;
                }

                if ban.contains(intron) {
                    continue;
                }

                if exon.0 < intron.0 && exon.1 > intron.1 {
                    hit = true;
                    if plot {
                        match query.strand {
                            '+' => Bed4::from(&query.chrom, intron.0, intron.1, &query.name)
                                .send(&mut blocks.as_mut().unwrap()),
                            '-' => Bed4::from(
                                &query.chrom,
                                SCALE - intron.1,
                                SCALE - intron.0,
                                &query.name,
                            )
                            .send(&mut blocks.as_mut().unwrap()),
                            _ => panic!("Invalid strand"),
                        }

                        continue;
                    }

                    // implement the logic to change the name based on UTR overlaping + in-frame

                    break;
                } else {
                    continue;
                }
            }
        }

        if hit {
            let line = query.line().to_owned();
            hits.push(line);
        } else {
            let line = query.line().to_owned();
            pass.push(line);
        }
    }

    (hits, pass, blocks)
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
            plot: false,
        };

        assert!(detect_intron_retentions(args).is_ok());
    }
}
