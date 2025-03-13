use anyhow::Result;
use config::{
    exonic_overlap, get_progress_bar, splice_site_overlap, tsv_to_map, write_objs,
    FusionDetectionValue, MatchType, ModuleDescriptor, ModuleMap, ModuleType, FUSIONS,
    FUSION_FAKES, FUSION_FREE, FUSION_RATIO_THRESHOLD, FUSION_REVIEW, SCALE,
};
use dashmap::DashSet;
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{
    buckerize, combine, packbed, par_reader, parse_tracks, unpack, GenePred, RefGenePred,
};
use rayon::prelude::*;
use serde_json::Value;

use std::sync::{
    atomic::{AtomicU32, Ordering},
    Arc,
};

use crate::cli::Args;
use crate::utils::{prepare_refs, unpack_blacklist, IsoformParser};

pub fn detect_fusions(args: Args) -> Result<()> {
    info!("Preparing files for fusion detection...");

    let ref_str = prepare_refs(args.refs)?;
    let refs = parse_tracks::<GenePred>(ref_str.as_str(), config::OverlapType::Exon, true)?; // INFO: OverlapType does not matter with TOGA
    let query = unpack(args.query, config::OverlapType::Exon, false)?;

    let (tracks, n) = combine(refs, query);
    let buckets = buckerize(
        tracks,
        config::OverlapType::Exon, // INFO: needs to be exonic to account for fused UTRs
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

    let match_type = MatchType::from(args.intron_match);

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

            let (fusions, no_fusions, fake_fusions, review, _, is_dirty) =
                process_component(comp, banned, args.recover, match_type).unwrap_or_default();

            acc.add(fusions, no_fusions, review, fake_fusions);

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

    if !acc.fakes.is_empty() {
        write_objs(&acc.fakes, FUSION_FAKES);
    }

    Ok(())
}

struct ParallelAccumulator {
    fusions: DashSet<String>,
    review: DashSet<String>,
    passes: DashSet<String>,
    fakes: DashSet<String>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            fusions: DashSet::new(),
            review: DashSet::new(),
            passes: DashSet::new(),
            fakes: DashSet::new(),
        }
    }
}

impl ParallelAccumulator {
    fn add(
        &self,
        fusions: Vec<String>,
        passes: Vec<String>,
        review: Option<Vec<String>>,
        fakes: Vec<String>,
    ) {
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

        if !fakes.is_empty() {
            fakes.into_iter().for_each(|fake| {
                self.fakes.insert(fake);
            });
        }
    }

    fn add_fusions(&self, fusions: Vec<String>) {
        fusions.into_iter().for_each(|fusion| {
            self.fusions.insert(fusion);
        });
    }

    fn add_passes(&self, passes: Vec<String>) {
        passes.into_iter().for_each(|pass| {
            self.passes.insert(pass);
        });
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
    match_type: MatchType,
) -> Option<(
    Vec<String>,
    Vec<String>,
    Vec<String>,
    Option<Vec<String>>,
    HashMap<String, Box<dyn ModuleMap>>,
    bool,
)> {
    if component.1.is_empty() {
        return None;
    }

    let mut descriptor = HashMap::new();

    let mut fusions = Vec::new();
    let mut fake_fusions = Vec::new();
    let mut no_fusions = Vec::new();

    let mut is_dirty = false;

    let refs = &component.0;
    let queries = &component.1;

    let genes = refs.get_names_split();

    let comp_size = queries.len() + genes.len();
    let query_size = queries.len();
    let (mut real_fusion_count, mut fake_fusion_count, totals) =
        (0_f32, 0_f32, queries.len() as f32);

    if genes.len() > 1 {
        // INFO: fusion loci [more than one gene in component]
        // INFO: we create a per-gene collection of exons
        let ref_exons = refs.smash_exons_by_name();
        let ref_introns = refs.smash_introns_by_name();

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

            // INFO: check if query read overlaps any of the gene exon collections
            // WARN: we are not keeping track of gene names here
            for r_exons in ref_exons.iter() {
                if exonic_overlap(&query.exons, r_exons) {
                    count += 1;
                }
            }

            // INFO: query read overlaps more than one gene exon collection
            if count > 1 {
                // INFO: we need to pass a 2nd check following this logic:
                //
                //  gene1:  XXXX---XXX--XXX
                //  gene2:                   XXX---XXX---XXX
                //  normal: XXXX---XXX--XXX-XXXXXXXXXXXX [does not follow ex-in struct]
                //  fusion: XXXX---XXX--XXX--XXX---XXX---XXX [follows ex-in struct]
                //
                // INFO: where not all reads that overlap 2 different genes
                // INFO: should be catalogued as fusions!
                let mut splicing_overlaps = 0_f32;
                for r_introns in ref_introns.iter() {
                    if splice_site_overlap(&query.introns, r_introns, match_type) {
                        splicing_overlaps += 1.0;
                    }
                }

                if splicing_overlaps > 1.0 {
                    real_fusion_count += 1.0;
                    fusions.push(query.line.clone());

                    handle
                        .set_value(
                            Box::new(FusionDetectionValue::IsFusedRead),
                            Value::Bool(true),
                        )
                        .ok();
                    handle
                        .set_value(
                            Box::new(FusionDetectionValue::FusionInFrame),
                            Value::Bool((query.cds_end - query.cds_start) % 3 == 0),
                        )
                        .ok();

                    let location = match query.strand {
                        config::Strand::Forward => {
                            format!("{}:{}-{}", query.chrom, query.start, query.end)
                        }
                        config::Strand::Reverse => {
                            format!(
                                "{}:{}-{}",
                                query.chrom,
                                SCALE - query.end,
                                SCALE - query.start
                            )
                        }
                    };
                    handle
                        .set_value(
                            Box::new(FusionDetectionValue::LocationOfFusion),
                            Value::String(location),
                        )
                        .ok();
                } else {
                    fake_fusion_count += 1.0;
                    fake_fusions.push(query.line.clone());
                    handle
                        .set_value(
                            Box::new(FusionDetectionValue::IsFusedRead),
                            Value::Bool(false),
                        )
                        .ok();
                }
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
        // INFO: species-specific gene | intergenic region | missing gene in refs
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
        // INFO: all good, no fusions here
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

    let whole_ratio = (real_fusion_count + fake_fusion_count) / totals;
    let real_ratio = real_fusion_count / totals;
    let fake_ratio = fake_fusion_count / totals;

    if recover {
        // INFO: if the fusion ratio in the component is above the threshold,
        // INFO: mark all queries as dirty and submit them for revie
        if real_ratio >= FUSION_RATIO_THRESHOLD {
            is_dirty = true;
            let mut review = vec![];

            for query in queries.iter() {
                review.push(query.line.clone());
                let handle = descriptor.get_mut(&query.name).unwrap();

                handle
                    .set_value(
                        Box::new(FusionDetectionValue::WholeComponentFusionRatio),
                        serde_json::json!(whole_ratio),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::RealComponentFusionRatio),
                        serde_json::json!(real_ratio),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::FakeComponentFusionRatio),
                        serde_json::json!(fake_ratio),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::IsDirtyComponent),
                        Value::Bool(true),
                    )
                    .ok();
            }

            return Some((vec![], vec![], vec![], Some(review), descriptor, is_dirty));
        } else {
            for query in queries.iter() {
                let handle = descriptor.get_mut(&query.name).unwrap();

                handle
                    .set_value(
                        Box::new(FusionDetectionValue::WholeComponentFusionRatio),
                        serde_json::json!(whole_ratio),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::RealComponentFusionRatio),
                        serde_json::json!(real_ratio),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(FusionDetectionValue::FakeComponentFusionRatio),
                        serde_json::json!(fake_ratio),
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
                    Box::new(FusionDetectionValue::WholeComponentFusionRatio),
                    serde_json::json!(whole_ratio),
                )
                .ok();
            handle
                .set_value(
                    Box::new(FusionDetectionValue::RealComponentFusionRatio),
                    serde_json::json!(real_ratio),
                )
                .ok();
            handle
                .set_value(
                    Box::new(FusionDetectionValue::FakeComponentFusionRatio),
                    serde_json::json!(fake_ratio),
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

    Some((
        fusions,
        no_fusions,
        fake_fusions,
        None,
        descriptor,
        is_dirty,
    ))
}

pub fn detect_fusions_with_mapping(args: Args) -> Result<()> {
    let contents = par_reader(args.refs).expect("ERROR: Failed to read isoform file(s)!");
    let refs = tsv_to_map::<IsoformParser, String>(Arc::new(contents), 0, 1)
        .expect("ERROR: Failed to parse isoform file(s)!");

    let buckets = packbed(
        args.query,
        None,
        args.overlap_type,
        packbed::PackMode::Paired,
    )
    .expect("ERROR: Failed to pack query reads!");

    let pb = get_progress_bar(buckets.len() as u64, "Detecting fusions in mapping mode...");

    let counter = ParallelCounter::default();
    let acc = ParallelAccumulator::default();

    buckets.par_iter().for_each(|bucket| {
        let _ = bucket.key();
        let components = bucket.value().to_owned();
        counter.inc_comp(components.len() as u32);

        components.into_par_iter().for_each(|comp| {
            let comp = comp
                .as_any()
                .downcast_ref::<(Vec<GenePred>, Vec<GenePred>)>()
                .expect("ERROR: Failed to downcast to RefGenePred");

            let (fusions, no_fusions) = process_mapping_component(comp, &refs);

            acc.add_passes(no_fusions);
            if let Some(f) = fusions {
                acc.add_fusions(f);
                // counter.inc_dirty(1);
            }
        });

        pb.inc(1);
    });

    pb.finish_and_clear();

    info!("Detected fusions: {}", acc.fusions.len());
    info!("Fusion-free reads: {}", acc.passes.len());

    [&acc.fusions, &acc.passes]
        .par_iter()
        .zip([FUSIONS, FUSION_FREE].par_iter())
        .for_each(|(rx, path)| write_objs(&rx, path));

    Ok(())
}

fn process_mapping_component(
    component: &(Vec<GenePred>, Vec<GenePred>),
    isoforms: &HashMap<String, Vec<String>>,
) -> (Option<Vec<String>>, Vec<String>) {
    let query: &Vec<GenePred> = component.0.as_ref(); // INFO: treating refs back to queries, discarding fake queries

    let mut genes = HashMap::new();

    query.iter().for_each(|query| {
        let tx_to_gene = isoforms
            .get(&query.name)
            .expect("ERROR: Failed to get isoforms!");

        // INFO: creates fills a hashmap with genes as keys and exons as values
        tx_to_gene.iter().for_each(|gene| {
            let exons = genes.entry(gene).or_insert_with(HashSet::new);
            exons.extend(query.exons.iter().cloned());
        });
    });

    // INFO: if at any point names is > 1, we have a fusion loci
    if genes.len() > 1 {
        let mut fusions = vec![];
        let mut no_fusions = vec![];

        query.iter().for_each(|query| {
            let mut count = 0;

            for (_, exons) in genes.iter() {
                if exonic_overlap(&query.exons, exons) {
                    count += 1;
                }
            }

            if count > 1 {
                fusions.push(query.line.clone());
            } else {
                no_fusions.push(query.line.clone());
            }
        });

        return (Some(fusions), no_fusions);
    }

    return (None, query.iter().map(|q| q.line.clone()).collect());
}
