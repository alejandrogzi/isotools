//! Core module for detecting intron retentions in a query set of reads
//!
//!
//! expected output type: HashMap/DashMap<ReadName, RetentionDescriptor>
//! expected output structure: ReadName -> <T> where T: ModuleDescriptor
//!     - ReadName is the name of the read
//!     - <T> is a struct implementing ModuleDescriptor trait
//!
//! expected output example:
//!     {
//!        "read1": {
//!           "intron_retention": true,
//!           "retains_rt_intron": true,
//!           "is_retention_supported": false,
//!           "is_retention_supported_map": [false, false],
//!           "retention_support_ratio": [0.23, 0.4],
//!           "exon_support_ratio": [0.65, 0.43],
//!           "number_of_retentions": 2,
//!           "number_of_unsupported_retentions": 2,
//!           "number_of_support_retentions": 0,
//!           "location_of_retention": ["chr1:1000-2000", "chr1:3000-4000"],
//!           "retention_in_cds": [false, true],
//!           "retention_in_utr": [true, false],
//!         }
//!     }

use std::collections::BTreeSet;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Arc;

use anyhow::Result;
use config::{
    get_progress_bar, write_objs, IntronRetentionValue, ModuleDescriptor, ModuleMap, ModuleType,
    BED3, INTRON_RETENTIONS, INTRON_RETENTION_FREE, INTRON_RETENTION_RECOVERY_THRESHOLD,
    OVERLAP_CDS, OVERLAP_EXON, RETENTION_RATIO_THRESHOLD, SCALE,
};
use dashmap::DashSet;
use hashbrown::{HashMap, HashSet};
use log::{info, warn};
use packbed::{packbed, GenePred, RefGenePred};
use rayon::prelude::*;
use serde_json::Value;

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

    let tracks = packbed(args.refs, Some(args.query), OVERLAP_CDS, OVERLAP_EXON)?;
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let hit_acc: DashSet<String> = DashSet::new();
    let pass_acc: DashSet<String> = DashSet::new();
    let misc_acc: DashSet<String> = DashSet::new();

    let pb = get_progress_bar(tracks.len() as u64, "Processing...");
    let dirty_count = AtomicU32::new(0);
    let n_comps = AtomicU32::new(0);

    tracks.par_iter().for_each(|bucket| {
        let chr = bucket.key();
        let components = bucket.value().to_owned();
        n_comps.fetch_add(components.len() as u32, Ordering::Relaxed);

        let binding = HashSet::new();
        let banned = blacklist.get(chr).unwrap_or(&binding);

        components.into_par_iter().for_each(|comp| {
            let (hits, pass, blocks, descriptor, is_dirty) =
                process_component(comp, banned, args.plot, args.recover);

            hits.into_iter().for_each(|hit| {
                hit_acc.insert(hit);
            });
            pass.into_iter().for_each(|p| {
                pass_acc.insert(p);
            });
            if let Some(b) = blocks {
                misc_acc.insert(b);
            }

            if is_dirty {
                dirty_count.fetch_add(1, Ordering::Relaxed);
            }
        });

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!("Reads with retained introns: {}", hit_acc.len());

    if args.recover {
        warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            dirty_count,
            dirty_count.load(Ordering::Relaxed) as f64 / n_comps.load(Ordering::Relaxed) as f64
                * 100.0
        );
    }

    [hit_acc, pass_acc]
        .par_iter()
        .zip([INTRON_RETENTIONS, INTRON_RETENTION_FREE].par_iter())
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
/// use iso_intron::core::process_component;
/// use packbed::{Bed12, GenePred, RefGenePred};
/// use hashbrown::HashSet;
/// use std::sync::Arc;
///
/// let line = "chr1\t100\t200\ttest\t0\t+\t100\t200\t0\t1\t1,1\t0,100";
/// let gp = Bed12::parse(line, true, true).unwrap();
/// let comp = (RefGenePred::from(vec![gp.clone()]), vec![gp]);
/// let banned = HashSet::from_iter(vec![(100, 200)]);
/// let plot = false;
///
/// let rs = process_component(comp, &banned, plot);
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
    recover: bool,
) -> (
    Vec<String>,
    Vec<String>,
    Option<String>,
    HashMap<String, Box<dyn ModuleMap>>,
    bool,
) {
    let mut hits = Vec::new();
    let mut pass = Vec::new();
    let mut blocks = if plot { Some(String::new()) } else { None };

    let mut tmp_dirt = Vec::new();
    let mut exon_owners = BTreeSet::new();
    let mut introns_owned = BTreeSet::new();

    let mut descriptor = HashMap::new();

    let refs = comp.0;
    let queries = comp.1;

    let comp_size = queries.len() + refs.reads.len();
    let (mut count, totals) = (0_f32, queries.len() as f32);

    for query in queries.iter() {
        let query = Arc::new(query);
        let mut hit: bool = false;
        let five_utr = query.get_five_utr();
        let three_utr = query.get_three_utr();

        descriptor.insert(
            query.name.clone(),
            ModuleDescriptor::with_schema(ModuleType::IntronRetention),
        );

        let mut number_of_retentions = 0_u32;
        let mut location_of_retentions = Vec::new();
        let mut retention_in_cds = Vec::new();
        let mut retention_in_utr = Vec::new();
        let mut retention_in_frame = Vec::new();

        for exon in &query.exons {
            for intron in &refs.introns {
                if exon.1 < intron.0 {
                    break;
                }

                if ban.contains(intron) {
                    continue;
                }

                if exon.0 < intron.0 && exon.1 > intron.1 {
                    if !hit {
                        count += 1.0;
                        tmp_dirt.push(query.clone());
                    }

                    hit = true;
                    exon_owners.insert((exon.0, exon.1));
                    introns_owned.insert((intron.0, intron.1));

                    number_of_retentions += 1;

                    if (five_utr.0 <= intron.0 && intron.0 <= five_utr.1)
                        || (three_utr.0 <= intron.1 && intron.1 <= three_utr.1)
                    {
                        retention_in_cds.push(false);
                        retention_in_utr.push(true);
                    } else {
                        retention_in_cds.push(true);
                        retention_in_utr.push(false);
                    }

                    if (intron.1 - intron.0) % 3 == 0 {
                        retention_in_frame.push(true);
                    } else {
                        retention_in_frame.push(false);
                    }

                    match query.strand {
                        '+' => {
                            location_of_retentions
                                .push(format!("{}:{}-{}", query.chrom, intron.0, intron.1,));

                            if plot {
                                Bed4::from(&query.chrom, intron.0 - 1, intron.1 + 1, &query.name)
                                    .send(&mut blocks.as_mut().unwrap())
                            }
                        }
                        '-' => {
                            location_of_retentions.push(format!(
                                "{}:{}-{}",
                                query.chrom,
                                SCALE - intron.1,
                                SCALE - intron.0,
                            ));

                            if plot {
                                Bed4::from(
                                    &query.chrom,
                                    SCALE - intron.1 - 1,
                                    SCALE - intron.0 + 1,
                                    &query.name,
                                )
                                .send(&mut blocks.as_mut().unwrap())
                            }
                        }
                        _ => panic!("Invalid strand"),
                    }

                    continue;
                } else {
                    continue;
                }
            }
        }

        let handle = descriptor.get_mut(&query.name).unwrap();

        handle
            .set_value(
                Box::new(IntronRetentionValue::ComponentSize),
                Value::Number(comp_size.into()),
            )
            .ok();
        handle
            .set_value(
                Box::new(IntronRetentionValue::RefComponentSize),
                Value::Number(comp_size.into()),
            )
            .ok();
        handle
            .set_value(
                Box::new(IntronRetentionValue::QueryComponentSize),
                Value::Number(comp_size.into()),
            )
            .ok();

        if hit {
            let line = query.line().to_owned();
            hits.push(line);

            handle
                .set_value(
                    Box::new(IntronRetentionValue::IsIntronRetention),
                    Value::Bool(true),
                )
                .ok();
            handle
                .set_value(
                    Box::new(IntronRetentionValue::NumberOfRetentions),
                    Value::Number(number_of_retentions.into()),
                )
                .ok();
            handle
                .set_value(
                    Box::new(IntronRetentionValue::RetentionLocation),
                    Value::Array(
                        location_of_retentions
                            .into_iter()
                            .map(Value::String)
                            .collect(),
                    ),
                )
                .ok();
            handle
                .set_value(
                    Box::new(IntronRetentionValue::IsRetentionInCds),
                    Value::Array(retention_in_cds.into_iter().map(Value::Bool).collect()),
                )
                .ok();
            handle
                .set_value(
                    Box::new(IntronRetentionValue::IsRetentionInUtr),
                    Value::Array(retention_in_utr.into_iter().map(Value::Bool).collect()),
                )
                .ok();
            handle
                .set_value(
                    Box::new(IntronRetentionValue::IsIntronRetainedInFrame),
                    Value::Array(retention_in_frame.into_iter().map(Value::Bool).collect()),
                )
                .ok();
        } else {
            let line = query.line().to_owned();
            pass.push(line);
        }
    }

    // after classying reads, we check bucket frequencies
    // if the number of truncated reads is greater than 50%
    // of the total reads in the bucket, we consider the bucket
    // to be dirty and return the reads for recovery if args.recover
    let ratio = count / totals;
    let mut is_dirty = false;
    if recover {
        if ratio >= RETENTION_RATIO_THRESHOLD {
            // log::warn!("Bucket {:?} is dirty -> {}", refs.reads, ratio);
            let dirt = tmp_dirt;
            is_dirty = true;

            let new_passes = recover_from_dirt(
                dirt,
                exon_owners,
                introns_owned,
                &refs,
                &mut descriptor,
                ratio,
            );
            new_passes.iter().for_each(|p| {
                hits.retain(|x| x != p);
                pass.push(p.to_owned());
            });
        } else {
            for query in queries.iter() {
                let handle = descriptor.get_mut(&query.name).unwrap();

                handle
                    .set_value(
                        Box::new(IntronRetentionValue::ComponentRetentionRatio),
                        serde_json::json!(ratio),
                    )
                    .ok();
            }
        }
    }

    // dbg!(&descriptor);

    (hits, pass, blocks, descriptor, is_dirty)
}

pub fn recover_from_dirt(
    dirt: Vec<Arc<&GenePred>>,
    exon_owners: BTreeSet<(u64, u64)>,
    introns_owned: BTreeSet<(u64, u64)>,
    refs: &RefGenePred,
    descriptor: &mut HashMap<String, Box<dyn ModuleMap>>,
    ratio: f32,
) -> Vec<String> {
    let mut local_passes = vec![];
    let mut owner_ratios = HashMap::new();
    let mut owned_ratios = HashMap::new();

    // exons: 1) build ratios
    for owner in exon_owners.iter() {
        let mut background = 0_f32;
        let mut count = 0_f32;
        for read in refs.reads.iter() {
            let ref_exons = &read.exons;
            let (owner_start, owner_end) = owner;

            if ref_exons.contains(&(*owner_start, *owner_end)) {
                count += 1.0;
            }

            // any read that spans the owner is considered background
            if (read.start < *owner_start) && (*owner_end < read.end) {
                background += 1.0;
            }
        }

        let ratio = count / background;
        owner_ratios.insert(owner, ratio);
    }

    // introns: 1) build ratios and 2) remove introns with low ratios
    let mut supported_introns = vec![];
    for owned in introns_owned.iter() {
        let mut background = 0.0;
        let mut count = 0.0;
        for read in refs.reads.iter() {
            let ref_introns = read.get_introns();
            let (owned_start, owned_end) = owned;

            if ref_introns.contains(&(*owned_start, *owned_end)) {
                count += 1.0;
            }

            // any read that spans the owned intron is considered background
            if (read.start < *owned_start) && (*owned_end < read.end) {
                background += 1.0;
            }
        }

        let ratio = count / background;
        owned_ratios.insert(owned, ratio);

        if ratio >= INTRON_RETENTION_RECOVERY_THRESHOLD {
            supported_introns.push(owned);
        }
    }

    for read in dirt.clone().iter_mut() {
        let handle = descriptor.get_mut(&read.name).unwrap();

        handle
            .set_value(
                Box::new(IntronRetentionValue::IsDirtyComponent),
                Value::Bool(true),
            )
            .ok();
        handle
            .set_value(
                Box::new(IntronRetentionValue::ComponentRetentionRatio),
                serde_json::json!(ratio),
            )
            .ok();

        let mut number_of_unsupported_retentions = 0;
        let mut number_of_supported_retentions = 0;
        let mut retention_support_ratio = vec![];
        let mut exon_support_ratio = vec![];
        let mut is_retention_supported_map = vec![];

        for exon in &read.exons {
            for intron in introns_owned.iter() {
                if exon.1 < intron.0 {
                    break;
                }

                if exon.0 < intron.0 && exon.1 > intron.1 {
                    retention_support_ratio.push(
                        owned_ratios
                            .get(intron)
                            .expect("ERROR: Intron not found in intron_owned_ratios"),
                    );
                    exon_support_ratio.push(
                        owner_ratios
                            .get(exon)
                            .expect("ERROR: Exon not found in exon_owner_ratios"),
                    );

                    if supported_introns.contains(&intron) {
                        // exon eats a 'likely' real intron, fill up the following descriptors:
                        // unrecover,
                        // is_retention_supported_map
                        number_of_supported_retentions += 1;
                        is_retention_supported_map.push(true);
                    } else {
                        // exon eats a 'likely' false intron, fill up the following descriptors:
                        // recover,
                        // is_retention_supported_map
                        number_of_unsupported_retentions += 1;
                        is_retention_supported_map.push(false);
                    }
                }
            }
        }

        if number_of_supported_retentions < 1
            && number_of_unsupported_retentions + number_of_supported_retentions > 0
        {
            local_passes.push(read.line().to_owned());

            handle
                .set_value(
                    Box::new(IntronRetentionValue::IsRetentionSupported),
                    serde_json::json!(false),
                )
                .ok();
        } else {
            handle
                .set_value(
                    Box::new(IntronRetentionValue::IsRetentionSupported),
                    serde_json::json!(true),
                )
                .ok();
        }

        handle
            .set_value(
                Box::new(IntronRetentionValue::NumberOfUnrecovers),
                serde_json::json!(number_of_supported_retentions),
            )
            .ok();
        handle
            .set_value(
                Box::new(IntronRetentionValue::NumberOfRecovers),
                serde_json::json!(number_of_unsupported_retentions),
            )
            .ok();
        handle
            .set_value(
                Box::new(IntronRetentionValue::RetentionSupportRatio),
                serde_json::json!(retention_support_ratio),
            )
            .ok();
        handle
            .set_value(
                Box::new(IntronRetentionValue::ExonSupportRatio),
                serde_json::json!(exon_support_ratio),
            )
            .ok();
        handle
            .set_value(
                Box::new(IntronRetentionValue::IsRetentionSupportedMap),
                serde_json::json!(is_retention_supported_map),
            )
            .ok();
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
            plot: false,
            recover: false,
        };

        assert!(detect_intron_retentions(args).is_ok());
    }
}
