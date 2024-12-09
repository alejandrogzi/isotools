use anyhow::Result;
use config::{
    get_progress_bar, write_objs, FusionDetectionValue, ModuleDescriptor, ModuleMap, ModuleType,
    FUSIONS, FUSION_FREE, FUSION_RATIO_THRESHOLD, FUSION_REVIEW, OVERLAP_CDS, OVERLAP_EXON,
};
use dashmap::DashSet;
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{buckerize, combine, parse_tracks, unpack, GenePred, RefGenePred};
use rayon::prelude::*;
use serde_json::Value;

use std::sync::atomic::{AtomicU32, Ordering};

use crate::cli::Args;
use crate::utils::{exonic_overlap, prepare_refs, unpack_blacklist};

pub fn detect_fusions(args: Args) -> Result<()> {
    info!("Preparing files for fusion detection...");

    let ref_str = prepare_refs(args.refs)?;
    let refs = parse_tracks(ref_str.as_str(), OVERLAP_CDS, true)?;
    let query = unpack(args.query, OVERLAP_CDS, false)?;

    let (tracks, n) = combine(refs, query);
    let buckets = buckerize(tracks, OVERLAP_CDS, OVERLAP_EXON, n);

    let blacklist = if !args.blacklist.is_empty() {
        unpack_blacklist(args.blacklist, true, false).ok()
    } else {
        None
    };
    let pb = get_progress_bar(buckets.len() as u64, "Detecting fusions...");

    let n_comps = AtomicU32::new(0);
    let dirty_count = AtomicU32::new(0);

    let hits: DashSet<String> = DashSet::new();
    let passes: DashSet<String> = DashSet::new();
    let reviews: DashSet<String> = DashSet::new();

    buckets.par_iter().for_each(|bucket| {
        let chr = bucket.key();
        let components = bucket.value().to_owned();
        n_comps.fetch_add(components.len() as u32, Ordering::Relaxed);

        let binding = HashSet::new();
        let banned = if let Some(bl) = blacklist.as_ref() {
            bl.get(chr).unwrap_or(&binding)
        } else {
            &binding
        };

        components.into_par_iter().for_each(|comp| {
            let (fusions, no_fusions, review, descriptor, is_dirty) =
                process_component(comp, banned, args.recover);

            fusions.into_iter().for_each(|hit| {
                hits.insert(hit);
            });
            no_fusions.into_iter().for_each(|p| {
                passes.insert(p);
            });
            if let Some(r) = review {
                r.into_iter().for_each(|r| {
                    reviews.insert(r);
                });
            }

            if is_dirty {
                dirty_count.fetch_add(1, Ordering::Relaxed);
            }
        });

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!("Detected fusions: {}", hits.len());
    info!("Fusion-free reads: {}", passes.len());

    if args.recover {
        log::warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            dirty_count,
            dirty_count.load(Ordering::Relaxed) as f64 / n_comps.load(Ordering::Relaxed) as f64
                * 100.0
        );

        if reviews.len() > 0 {
            write_objs(&reviews, FUSION_REVIEW);
        }
    }

    [&hits, &passes]
        .par_iter()
        .zip([FUSIONS, FUSION_FREE].par_iter())
        .for_each(|(rx, path)| write_objs(&rx, path));

    Ok(())
}

fn process_component(
    component: (RefGenePred, Vec<GenePred>),
    banned: &HashSet<String>,
    recover: bool,
) -> (
    Vec<String>,
    Vec<String>,
    Option<Vec<String>>,
    HashMap<String, Box<dyn ModuleMap>>,
    bool,
) {
    let mut descriptor = HashMap::new();

    let mut fusions = Vec::new();
    let mut no_fusions = Vec::new();

    let mut is_dirty = false;

    let refs = component.0;
    let queries = component.1;

    let genes = refs.get_names();

    let comp_size = queries.len() + genes.len();
    let query_size = queries.len();
    let (mut fusion_count, totals) = (0_f32, queries.len() as f32);

    if genes.len() > 1 {
        // fusion
        // log::warn!("Multiple ref genes in component: {:?}", genes);
        let refs = refs.smash_exons_by_name();

        queries.iter().for_each(|query| {
            if banned.contains(&query.name) {
                no_fusions.push(query.line.clone());
                return;
            }

            descriptor.insert(
                query.name.clone(),
                ModuleDescriptor::with_schema(ModuleType::FusionDetection),
            );
            let handle = descriptor.get_mut(&query.name).unwrap();

            let mut count = 0;

            // check if query read overlaps any of the genes
            for ref_exons in refs.iter() {
                if exonic_overlap(&query.exons, &ref_exons) {
                    count += 1;
                }
            }

            if count > 1 {
                fusions.push(query.line.clone());
                fusion_count += 1.0;

                handle
                    .set_value(
                        Box::new(FusionDetectionValue::IsFusedRead),
                        Value::Bool(true),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::LocationOfFusion),
                        Value::String(format!("{}:{}-{}", query.chrom, query.start, query.end)),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::FusionInFrame),
                        Value::Bool((query.cds_end - query.cds_start) % 3 == 0),
                    )
                    .ok();
            } else {
                no_fusions.push(query.line.clone());

                handle
                    .set_value(
                        Box::new(FusionDetectionValue::IsFusedRead),
                        Value::Bool(false),
                    )
                    .ok();
            }

            handle
                .set_value(
                    Box::new(FusionDetectionValue::ComponentSize),
                    Value::Number(comp_size.into()),
                )
                .ok();
            handle
                .set_value(
                    Box::new(FusionDetectionValue::RefComponentSize),
                    Value::Number(genes.len().into()),
                )
                .ok();
            handle
                .set_value(
                    Box::new(FusionDetectionValue::QueryComponentSize),
                    Value::Number(query_size.into()),
                )
                .ok();
        });
    } else if genes.is_empty() {
        // species-specific gene | intergenic region | missing gene in refs
        // log::warn!("No ref genes in component, check loci: {:?}", queries);
        queries.iter().for_each(|query| {
            no_fusions.push(query.line.clone());
        });
    } else {
        // all good, no fusions here
        queries.iter().for_each(|query| {
            no_fusions.push(query.line.clone());
        });
    }

    let ratio = fusion_count / totals;
    if recover {
        if ratio >= FUSION_RATIO_THRESHOLD {
            is_dirty = true;
            let mut review = vec![];

            for query in queries.iter() {
                review.push(query.line.clone());
                let handle = descriptor.get_mut(&query.name).unwrap();

                handle
                    .set_value(
                        Box::new(FusionDetectionValue::ComponentFusionRatio),
                        serde_json::json!(ratio),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::IsDirtyComponent),
                        Value::Bool(true),
                    )
                    .ok();
            }

            return (vec![], vec![], Some(review), descriptor, is_dirty);
        } else {
            for query in queries.iter() {
                let handle = descriptor.get_mut(&query.name).unwrap();

                handle
                    .set_value(
                        Box::new(FusionDetectionValue::ComponentFusionRatio),
                        serde_json::json!(ratio),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::IsDirtyComponent),
                        Value::Bool(false),
                    )
                    .ok();
            }
        }
    } else {
        for query in queries.iter() {
            let handle = descriptor.get_mut(&query.name).unwrap();

            handle
                .set_value(
                    Box::new(FusionDetectionValue::ComponentFusionRatio),
                    serde_json::json!(ratio),
                )
                .ok();
            handle
                .set_value(
                    Box::new(FusionDetectionValue::IsDirtyComponent),
                    Value::Bool(false),
                )
                .ok();
        }
    }

    // dbg!(&descriptor);

    (fusions, no_fusions, None, descriptor, is_dirty)
}
