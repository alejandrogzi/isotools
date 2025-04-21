use std::collections::BTreeSet;

use anyhow::Result;
use config::{
    get_progress_bar, write_descriptor, ModuleDescriptor, ModuleMap, ModuleType, OverlapType,
    StartTruncationValue, TRUNCATION_DESCRIPTOR, TRUNCATION_RECOVERY_THRESHOLD,
    TRUNCATION_THRESHOLD,
};
use dashmap::DashMap;
use hashbrown::{HashMap, HashSet};
use log::{info, warn};
use packbed::{packbed, BedPackage, GenePred, RefGenePred};
use rayon::prelude::*;
use serde_json::Value;

use crate::cli::Args;
use crate::utils::*;

pub fn detect_truncations(args: Args) -> Result<DashMap<String, Box<dyn ModuleMap>>> {
    info!("Detecting 5'end truncations...");

    let tracks = packbed(
        args.refs,
        Some(args.query),
        OverlapType::CDSBound,
        packbed::PackMode::Default,
    )?;
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let accumulator = ParallelAccumulator::default();
    let counter = ParallelCounter::default();

    let pb = get_progress_bar(tracks.len() as u64, "Processing...");
    tracks.into_par_iter().for_each(|bucket| {
        let components = bucket.1;
        counter.inc_components(components.len() as u32);

        process_components(components, &blacklist, &accumulator, &counter, args.recover);

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!(
        "Reads with 5'end truncations: {}",
        accumulator.truncations.len()
    );

    if args.recover {
        let (count, ratio) = counter.get_stat();
        warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            count, ratio
        );
    }

    if !args.in_memory {
        info!("INFO: Writing results...");

        write_results(&accumulator);
        write_descriptor(&accumulator.descriptor, TRUNCATION_DESCRIPTOR);
    }

    Ok(accumulator.descriptor)
}

pub fn process_components(
    components: Vec<Box<dyn BedPackage>>,
    banned: &HashSet<String>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
    recover: bool,
) {
    components.into_par_iter().for_each(|comp| {
        let comp = comp
            .as_any_owned()
            .downcast::<(RefGenePred, Vec<GenePred>)>()
            .expect("ERROR: Failed to downcast to RefGenePred");

        let (truncations, no_truncations, descriptor, is_dirty) =
            process_component(comp, &banned, recover);

        truncations.into_iter().for_each(|hit| {
            accumulator.truncations.insert(hit);
        });
        no_truncations.into_iter().for_each(|pass| {
            accumulator.no_truncations.insert(pass);
        });

        descriptor.into_iter().for_each(|(name, desc)| {
            accumulator.descriptor.insert(name, desc);
        });

        if is_dirty {
            counter.inc_dirty();
        }
    });
}

#[inline(always)]
pub fn process_component(
    comp: Box<(RefGenePred, Vec<GenePred>)>,
    ban: &HashSet<String>,
    recover: bool,
) -> (
    Vec<String>,
    Vec<String>,
    HashMap<String, Box<dyn ModuleMap>>,
    bool,
) {
    let mut truncations = Vec::new();
    let mut pass = Vec::new();

    let mut tmp_dirt = Vec::new();
    let mut owners = BTreeSet::new();

    let mut descriptor = HashMap::new();

    let refs = comp.0;
    let queries = comp.1;

    let ref_starts = &refs.starts;
    let ref_middles = &refs.middles;

    let comp_size = queries.len() + refs.reads.len();
    let (mut t, totals) = (0_f32, queries.len() as f32);

    for query in queries.iter() {
        if ban.contains(query.name()) {
            continue;
        }

        descriptor.insert(
            query.name.clone(),
            ModuleDescriptor::with_schema(ModuleType::StartTruncation),
        );
        let handle = descriptor.get_mut(&query.name).unwrap();

        handle
            .set_value(
                Box::new(StartTruncationValue::ComponentSize),
                serde_json::json!(comp_size),
            )
            .ok();
        handle
            .set_value(
                Box::new(StartTruncationValue::RefComponentSize),
                serde_json::json!(refs.reads.len()),
            )
            .ok();
        handle
            .set_value(
                Box::new(StartTruncationValue::QueryComponentSize),
                serde_json::json!(queries.len()),
            )
            .ok();

        let (query_start, query_end) = query.get_first_exon();

        let is_complete = ref_starts.iter().any(|(s, e)| {
            if query_end < *s {
                return false;
            }

            (query_start >= *s) && (query_start < *e)
        });

        if is_complete {
            handle
                .set_value(
                    Box::new(StartTruncationValue::IsNovelStart),
                    Value::Bool(false),
                )
                .ok();

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

                handle
                    .set_value(
                        Box::new(StartTruncationValue::IsReadTruncated),
                        Value::Bool(true),
                    )
                    .ok();

                t += 1.0;
            } else {
                let line = query.line().to_owned();
                pass.push(line);

                handle
                    .set_value(
                        Box::new(StartTruncationValue::IsReadTruncated),
                        Value::Bool(false),
                    )
                    .ok();
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

                handle
                    .set_value(
                        Box::new(StartTruncationValue::IsReadTruncated),
                        Value::Bool(true),
                    )
                    .ok();

                t += 1.0;
            } else {
                let line = query.line().to_owned();
                pass.push(line);

                handle
                    .set_value(
                        Box::new(StartTruncationValue::IsNovelStart),
                        Value::Bool(true),
                    )
                    .ok();
            }
        }
    }

    // after classying reads, we check bucket frequencies
    // if the number of truncated reads is greater than 50%
    // of the total reads in the bucket, we consider the bucket
    // to be dirty and return the reads for recovery if args.recover
    let ratio = t / totals;
    let mut is_dirty = false;
    if recover {
        if ratio >= TRUNCATION_THRESHOLD {
            // warn!("Bucket {:?} is dirty -> {}", refs.reads, t / totals);
            is_dirty = true;
            let dirt = tmp_dirt;

            let new_passes = recover_from_dirt(dirt, owners, &refs, &mut descriptor, ratio);
            new_passes.iter().for_each(|p| {
                truncations.retain(|x| x != p);
                pass.push(p.to_owned());
            });
        }
    } else {
        drop(tmp_dirt)
    }

    // dbg!(&descriptor);

    return (truncations, pass, descriptor, is_dirty);
}

pub fn recover_from_dirt(
    mut dirt: Vec<&GenePred>,
    owners: BTreeSet<(&u64, &u64)>,
    refs: &RefGenePred,
    descriptor: &mut HashMap<String, Box<dyn ModuleMap>>,
    ratio: f32,
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

        // 2. if the owner (middle exon) has more than 50% support, we consider it
        //  a valid owner and keep reads truncated by that owner as
        //  truncated reads; otherwise, we consider the owner to be
        //  a weak ownner and send all truncated reads to the pass
        //  bucket
        let owner_ratio = count / background;
        if owner_ratio < TRUNCATION_RECOVERY_THRESHOLD {
            // send reads truncated by this owner to pass
            for read in dirt.clone().iter_mut() {
                let handle = descriptor.get_mut(&read.name).unwrap();

                handle
                    .set_value(
                        Box::new(StartTruncationValue::IsDirtyComponent),
                        Value::Bool(true),
                    )
                    .ok();

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

                    handle
                        .set_value(
                            Box::new(StartTruncationValue::TruncationSupportRatio),
                            serde_json::json!(owner_ratio),
                        )
                        .ok();
                    handle
                        .set_value(
                            Box::new(StartTruncationValue::IsTruncationSupported),
                            serde_json::json!(false),
                        )
                        .ok();
                    handle
                        .set_value(
                            Box::new(StartTruncationValue::ComponentTruncationRatio),
                            serde_json::json!(ratio),
                        )
                        .ok();
                }
            }
        } else {
            for read in dirt.iter() {
                let handle = descriptor.get_mut(&read.name).unwrap();

                handle
                    .set_value(
                        Box::new(StartTruncationValue::IsDirtyComponent),
                        Value::Bool(true),
                    )
                    .ok();

                let (owner_start, owner_end) = *owner;
                let (query_start, query_end) = read.get_first_exon();

                if (query_start >= *owner_start && query_start < *owner_end)
                    || (query_end > *owner_start && query_end <= *owner_end)
                    || (query_start < *owner_start && query_end > *owner_end)
                {
                    handle
                        .set_value(
                            Box::new(StartTruncationValue::IsTruncationSupported),
                            serde_json::json!(true),
                        )
                        .ok();
                    handle
                        .set_value(
                            Box::new(StartTruncationValue::ComponentTruncationRatio),
                            serde_json::json!(ratio),
                        )
                        .ok();
                    handle
                        .set_value(
                            Box::new(StartTruncationValue::TruncationSupportRatio),
                            serde_json::json!(owner_ratio),
                        )
                        .ok();
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
            in_memory: true,
        };

        assert!(detect_truncations(args).is_ok());
    }
}
