use anyhow::Result;
use config::{
    exonic_overlap, get_progress_bar, write_objs, FusionDetectionValue, ModuleDescriptor,
    ModuleMap, ModuleType, FUSIONS, FUSION_FREE, FUSION_RATIO_THRESHOLD, FUSION_REVIEW,
};
use dashmap::DashSet;
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{buckerize, combine, parse_tracks, unpack, GenePred, RefGenePred};
use rayon::prelude::*;
use serde_json::Value;

use std::sync::atomic::{AtomicU32, Ordering};

use crate::cli::Args;
use crate::utils::{prepare_refs, unpack_blacklist};

pub fn detect_fusions(args: Args) -> Result<()> {
    info!("Preparing files for fusion detection...");

    let ref_str = prepare_refs(args.refs)?;
    let refs = parse_tracks::<GenePred>(ref_str.as_str(), config::OverlapType::Exon, true)?;
    let query = unpack(args.query, config::OverlapType::CDS, false)?;

    let (tracks, n) = combine(refs, query);
    let buckets = buckerize(
        tracks,
        config::OverlapType::CDS,
        n,
        packbed::PackMode::Default,
    );

    let blacklist = if !args.blacklist.is_empty() {
        unpack_blacklist(args.blacklist, config::OverlapType::Exon, false).ok()
    } else {
        None
    };
    let pb = get_progress_bar(buckets.len() as u64, "Detecting fusions...");

    let counter = ParallelCounter::default();
    let acc = ParallelAccumulator::default();

    buckets.par_iter().for_each(|bucket| {
        let chr = bucket.key();
        let components = bucket.value().to_owned();
        counter.inc_comp(components.len() as u32);

        let binding = HashSet::new();
        let banned = if let Some(bl) = blacklist.as_ref() {
            bl.get(chr).unwrap_or(&binding)
        } else {
            &binding
        };

        components.into_par_iter().for_each(|comp| {
            let comp = comp
                .as_any()
                .downcast_ref::<(RefGenePred, Vec<GenePred>)>()
                .expect("ERROR: Failed to downcast to RefGenePred");

            let (fusions, no_fusions, review, descriptor, is_dirty) =
                process_component(comp, banned, args.recover);

            acc.add(fusions, no_fusions, review);

            if is_dirty {
                counter.inc_dirty(1);
            }
        });

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!("Detected fusions: {}", acc.fusions.len());
    info!("Fusion-free reads: {}", acc.passes.len());

    if args.recover {
        log::warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            counter.num_of_dirty,
            counter.load_ratio()
        );

        if acc.review.len() > 0 {
            write_objs(&acc.review, FUSION_REVIEW);
        }
    }

    [&acc.fusions, &acc.passes]
        .par_iter()
        .zip([FUSIONS, FUSION_FREE].par_iter())
        .for_each(|(rx, path)| write_objs(&rx, path));

    if args.recover {
        write_objs(&acc.review, FUSION_REVIEW);
    }

    Ok(())
}

struct ParallelAccumulator {
    fusions: DashSet<String>,
    review: DashSet<String>,
    passes: DashSet<String>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            fusions: DashSet::new(),
            review: DashSet::new(),
            passes: DashSet::new(),
        }
    }
}

impl ParallelAccumulator {
    fn add(&self, fusions: Vec<String>, passes: Vec<String>, review: Option<Vec<String>>) {
        fusions.into_iter().for_each(|fusion| {
            self.fusions.insert(fusion);
        });
        passes.into_iter().for_each(|pass| {
            self.passes.insert(pass);
        });
        if let Some(r) = review {
            r.into_iter().for_each(|r| {
                self.review.insert(r);
            });
        }
    }
}

struct ParallelCounter {
    num_of_comps: AtomicU32,
    num_of_dirty: AtomicU32,
}

impl Default for ParallelCounter {
    fn default() -> Self {
        Self {
            num_of_comps: AtomicU32::new(0),
            num_of_dirty: AtomicU32::new(0),
        }
    }
}

impl ParallelCounter {
    fn inc_comp(&self, value: u32) {
        self.num_of_comps.fetch_add(value, Ordering::Relaxed);
    }

    fn inc_dirty(&self, value: u32) {
        self.num_of_dirty.fetch_add(value, Ordering::Relaxed);
    }

    fn load_ratio(&self) -> f64 {
        self.num_of_dirty.load(Ordering::Relaxed) as f64
            / self.num_of_comps.load(Ordering::Relaxed) as f64
    }
}

fn process_component(
    component: &(RefGenePred, Vec<GenePred>),
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

    let refs = &component.0;
    let queries = &component.1;

    let genes = refs.get_names_split();

    let comp_size = queries.len() + genes.len();
    let query_size = queries.len();
    let (mut fusion_count, totals) = (0_f32, queries.len() as f32);

    if genes.len() > 1 {
        // INFO: fusion loci
        // INFO: we create a per-gene collection of exons
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

            // check if query read overlaps any of the gene exon collections
            // WARN: we are not keeping track of gene names here
            for ref_exons in refs.iter() {
                if exonic_overlap(&query.exons, &ref_exons) {
                    count += 1;
                }
            }

            if count > 1 {
                // query read overlaps more than one gene exon collection
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

            descriptor.insert(
                query.name.clone(),
                ModuleDescriptor::with_schema(ModuleType::FusionDetection),
            );
            let handle = descriptor.get_mut(&query.name).unwrap();

            handle
                .set_value(
                    Box::new(FusionDetectionValue::IsFusedRead),
                    Value::Bool(false),
                )
                .ok();
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
    } else {
        // all good, no fusions here
        queries.iter().for_each(|query| {
            no_fusions.push(query.line.clone());

            descriptor.insert(
                query.name.clone(),
                ModuleDescriptor::with_schema(ModuleType::FusionDetection),
            );
            let handle = descriptor.get_mut(&query.name).unwrap();

            handle
                .set_value(
                    Box::new(FusionDetectionValue::IsFusedRead),
                    Value::Bool(false),
                )
                .ok();
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
    }

    let ratio = fusion_count / totals;
    if recover {
        // if the fusion ratio in the component is above the threshold,
        // mark all queries as dirty and submit them for revie
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
