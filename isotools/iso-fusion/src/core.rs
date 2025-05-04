//! Core module for detecting fusions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting fusions
//! and processing the components of reads and introns in parallel.
//!
//! In short, each read is checked for the presence of fusions using
//! a set of reference coordinates. If a read has a true fusion, it
//! is discarded. If a read has an RT intron, it is also discarded.
//! The veracity of a 'fusion' is determined by the exact match of at
//! least 1 splice site in the reference. If a read has a fake fusion, it is
//! marked as such. The results are written to a set of files.

use anyhow::Result;
use config::{
    exonic_overlap, get_progress_bar, par_write_results, splice_site_overlap, tsv_to_map,
    write_descriptor, write_objs, FusionDetectionValue, MatchType, ModuleDescriptor, ModuleMap,
    ModuleType, FUSIONS, FUSION_DESCRIPTOR, FUSION_FAKES, FUSION_FREE, FUSION_RATIO_THRESHOLD,
    FUSION_REVIEW, SCALE,
};
use dashmap::DashMap;
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{
    buckerize, combine, packbed, par_reader, parse_tracks, unpack, BedPackage, GenePred,
    RefGenePred,
};
use rayon::prelude::*;
use serde_json::Value;

use std::collections::BTreeSet;
use std::path::PathBuf;
use std::sync::Arc;

use crate::cli::Args;
use crate::utils::{
    prepare_refs, unpack_blacklist, IsoformParser, ParallelAccumulator, ParallelCounter,
};

/// Detects fusions in a query set of reads
///
/// # Arguments
///
/// * `args` - The command line arguments
///
/// # Returns
///
/// * `Result<()>` - The result of the operation
///
/// # Example
///
/// ```rust, no_run
/// let args = Args::new();
/// detect_fusions(args).unwrap();
/// ```
pub fn detect_fusions(args: Args) -> Result<DashMap<String, Box<dyn ModuleMap>>> {
    info!("Preparing files for fusion detection...");

    let buckets = make_buckets(args.refs, args.query).expect("ERROR: Failed to make buckets!");
    let blacklist =
        unpack_blacklist(args.blacklist, config::OverlapType::Exon, false).unwrap_or_default();
    let match_type = MatchType::from(args.intron_match);

    let pb = get_progress_bar(buckets.len() as u64, "Detecting fusions...");

    let counter = ParallelCounter::default();
    let accumulator = ParallelAccumulator::default();

    buckets.into_par_iter().for_each(|bucket| {
        let chr = bucket.0;
        let components = bucket.1;

        counter.inc_comp(components.len() as u32);

        let binding = HashSet::new();
        let banned = blacklist.get(&chr).unwrap_or(&binding);

        process_components(
            components,
            banned,
            &accumulator,
            &counter,
            args.recover,
            match_type,
        );

        pb.inc(1);
    });

    pb.finish_and_clear();

    info!("Detected fusions: {}", accumulator.fusions.len());
    info!("Fusion-free reads: {}", accumulator.passes.len());

    if args.recover {
        log::warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            counter.num_of_dirty,
            counter.load_ratio()
        );
    }

    if !args.in_memory {
        write_descriptor(&accumulator.descriptor, FUSION_DESCRIPTOR);
        par_write_results(
            &accumulator,
            vec![
                args.prefix.join(FUSIONS),
                args.prefix.join(FUSION_FREE),
                args.prefix.join(FUSION_REVIEW),
                args.prefix.join(FUSION_FAKES),
            ],
            None,
        );
    }

    Ok(accumulator.descriptor)
}

/// Create a per-chromosome bucket of tracks for fusion detection
///
/// # Arguments
///
/// * `refs` - A vector of PathBuf containing the reference files
/// * `query` - A vector of PathBuf containing the query files
///
/// # Example
///
/// ```rust, no_run
/// use iso_fusion::core::make_buckets;
///
/// let refs = vec![PathBuf::from("path/to/refs.bed")];
/// let query = vec![PathBuf::from("path/to/query.bed")];
///
/// let buckets = make_buckets(refs, query).expect("ERROR: Failed to make buckets!");
///
/// assert_eq!(buckets.len(), 1);
/// ```
fn make_buckets(
    refs: Vec<PathBuf>,
    query: Vec<PathBuf>,
) -> Result<DashMap<String, Vec<Box<dyn BedPackage>>>> {
    let ref_str = prepare_refs(refs)?;
    let refs = parse_tracks::<GenePred>(ref_str.as_str(), config::OverlapType::Exon, true)?; // INFO: OverlapType does not matter with TOGA
    let query = unpack(query, config::OverlapType::Exon, false)?;

    let (tracks, n) = combine(refs, query);
    let buckets = buckerize(
        tracks,
        config::OverlapType::Exon, // INFO: needs to be exonic to account for fused UTRs
        n,
        packbed::PackMode::Default,
    );

    Ok(buckets)
}

/// Processes the components of reads in parallel
///
/// # Arguments
///
/// * `components` - The components to process per chromosome
/// * `banned` - the set of banned introns for the chromosome
/// * `accumulator` - the accumulator to use
/// * `counter` - The counter to use
/// * `recover` - Indicates if the component should be recovered
/// * `match_type` - The type of match to use
///
/// # Returns
///
/// * None
///
/// # Example
///
/// ```rust, no_run
/// let components = vec![];
/// let banned = HashSet::new();
/// let accumulator = ParallelAccumulator::default();
/// let counter = ParallelCounter::default();
/// let recover = false;
/// let match_type = MatchType::SpliceSite;
///
/// process_components(components, &banned, &accumulator, &counter, recover, match_type);
///
/// assert_eq!(accumulator.num_retentions(), 0);
/// assert_eq!(counter.num_components(), 0);
/// ```
fn process_components(
    components: Vec<Box<dyn BedPackage>>,
    banned: &HashSet<String>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
    recover: bool,
    match_type: MatchType,
) {
    components.into_par_iter().for_each(|comp| {
        let comp = comp
            .as_any()
            .downcast_ref::<(RefGenePred, Vec<GenePred>)>()
            .expect("ERROR: Failed to downcast to RefGenePred");

        let (fusions, no_fusions, fake_fusions, review, descriptor, is_dirty) =
            process_component(comp, banned, recover, match_type).unwrap_or_default();

        accumulator.add(fusions, no_fusions, review, fake_fusions, descriptor);

        if is_dirty {
            counter.inc_dirty(1);
        }
    });
}

/// FusionSchema struct
///
///
/// This struct is used to store the fusion schema for a read.
/// Mimics FusionDetectionValue enum and fields of FusionDetectionDescriptor
///
/// # Fields
///
/// * `is_fused_read` - Indicates if the read is a fusion read
/// * `is_fusion_supported` - Indicates if the fusion is supported
/// * `component_size` - The size of the component
/// * `ref_component_size` - The size of the reference component
/// * `query_component_size` - The size of the query component
/// * `whole_component_fusion_ratio` - The ratio of the whole component that is fused
/// * `real_component_fusion_ratio` - The ratio of the real component that is fused
/// * `fake_component_fusion_ratio` - The ratio of the fake component that is fused
/// * `is_dirty_component` - Indicates if the component is dirty
/// * `location_of_fusion` - The location of the fusion
/// * `fusion_in_frame` - Indicates if the fusion is in frame
///
pub struct FusionSchema {
    is_fused_read: Value,
    ref_component_size: Value,
    query_component_size: Value,
    whole_component_fusion_ratio: Value,
    real_component_fusion_ratio: Value,
    fake_component_fusion_ratio: Value,
    is_dirty_component: Value,
    location_of_fusion: Value,
    fusion_in_frame: Value,
}

impl FusionSchema {
    /// Fills up the descriptor with the values from the schema
    ///
    /// # Arguments
    ///
    /// * `descriptor` - The descriptor to fill up
    /// * `read` - The read to fill up the descriptor with
    ///
    /// # Returns
    ///
    /// * None
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut descriptor = HashMap::new();
    /// let read = GenePred::new();
    ///
    /// schema.difuse(&mut descriptor, &read);
    /// ```
    fn diffuse(&self, descriptor: &mut HashMap<String, Box<dyn ModuleMap>>, read: &GenePred) {
        descriptor.insert(
            read.name.clone(),
            ModuleDescriptor::with_schema(ModuleType::FusionDetection),
        );
        let handle = descriptor.get_mut(&read.name).unwrap();

        for (field, value) in self.as_vals() {
            handle.set_value(Box::new(field), value).ok();
        }
    }

    /// Converts the schema into an iterable collection
    ///
    /// # Returns
    ///
    /// * Vec<(FusionDetectionValue, Value)> - The vector of tuples
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let schema = FusionSchema::default();
    ///
    /// let vals = schema.as_vals();
    ///
    /// for (field, value) in vals {
    ///     println!("{}: {}", field, value);
    /// }
    /// ```
    fn as_vals(&self) -> Vec<(FusionDetectionValue, Value)> {
        vec![
            (
                FusionDetectionValue::IsFusedRead,
                self.is_fused_read.clone(),
            ),
            (
                FusionDetectionValue::RefFusionComponentSize,
                self.ref_component_size.clone(),
            ),
            (
                FusionDetectionValue::QueryFusionComponentSize,
                self.query_component_size.clone(),
            ),
            (
                FusionDetectionValue::WholeComponentFusionRatio,
                self.whole_component_fusion_ratio.clone(),
            ),
            (
                FusionDetectionValue::RealComponentFusionRatio,
                self.real_component_fusion_ratio.clone(),
            ),
            (
                FusionDetectionValue::FakeComponentFusionRatio,
                self.fake_component_fusion_ratio.clone(),
            ),
            (
                FusionDetectionValue::IsDirtyFusionComponent,
                self.is_dirty_component.clone(),
            ),
            (
                FusionDetectionValue::LocationOfFusion,
                self.location_of_fusion.clone(),
            ),
            (
                FusionDetectionValue::FusionInFrame,
                self.fusion_in_frame.clone(),
            ),
        ]
    }
}

/// Implements the Default trait for FusionSchema
///
/// # Example
///
/// ```rust, no_run
/// let schema = FusionSchema::default();
///
/// assert_eq!(schema.is_fused_read, Value::Null);
/// assert_eq!(schema.is_fusion_supported, Value::Null);
/// assert_eq!(schema.ref_component_size, Value::Null);
/// ```
impl Default for FusionSchema {
    fn default() -> Self {
        FusionSchema {
            is_fused_read: Value::Bool(false),
            ref_component_size: Value::Number(0.into()),
            query_component_size: Value::Number(0.into()),
            whole_component_fusion_ratio: Value::Null,
            real_component_fusion_ratio: Value::Null,
            fake_component_fusion_ratio: Value::Null,
            is_dirty_component: Value::Bool(false),
            location_of_fusion: Value::Null,
            fusion_in_frame: Value::Bool(false),
        }
    }
}

/// LocalCounter struct
///
/// This struct is used to count the number
/// of real and fake fusions
///
/// # Fields
///
/// * `real_fusion_count` - The count of real fusions
/// * `fake_fusion_count` - The count of fake fusions
/// * `totals` - The total number of reads
///
/// # Example
///
/// ```rust, no_run
/// let counter = LocalCounter::new(100);
///
/// assert_eq!(counter.real_fusion_count, 0.0);
/// assert_eq!(counter.fake_fusion_count, 0.0);
/// assert_eq!(counter.totals, 100.0);
/// ```
struct LocalCounter {
    real_fusion_count: f32,
    fake_fusion_count: f32,
    totals: f32,
}

impl LocalCounter {
    /// Creates a new LocalCounter
    ///
    /// # Arguments
    ///
    /// * `totals` - The total number of reads
    ///
    /// # Returns
    ///
    /// * Self - The new LocalCounter
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = LocalCounter::new(100);
    ///
    /// assert_eq!(counter.real_fusion_count, 0.0);
    /// assert_eq!(counter.fake_fusion_count, 0.0);
    /// assert_eq!(counter.totals, 100.0);
    /// ```
    fn new(totals: usize) -> Self {
        LocalCounter {
            real_fusion_count: 0.0,
            fake_fusion_count: 0.0,
            totals: totals as f32,
        }
    }

    /// Increments the real fusion count
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut counter = LocalCounter::new(100);
    ///
    /// counter.inc_real();
    /// assert_eq!(counter.real_fusion_count, 1.0);
    /// ```
    fn inc_real(&mut self) {
        self.real_fusion_count += 1.0;
    }

    /// Increments the fake fusion count
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut counter = LocalCounter::new(100);
    ///
    /// counter.inc_fake();
    /// assert_eq!(counter.fake_fusion_count, 1.0);
    /// ```
    fn inc_fake(&mut self) {
        self.fake_fusion_count += 1.0;
    }

    /// Gets the ratios of fusions
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = LocalCounter::new(100);
    ///
    /// let ratios = counter.get_ratios();
    /// assert_eq!(ratios.0, 0.0);
    /// assert_eq!(ratios.1, 0.0);
    /// assert_eq!(ratios.2, 0.0);
    ///
    /// counter.inc_real();
    /// counter.inc_fake();
    ///
    /// let ratios = counter.get_ratios();
    /// assert_eq!(ratios.0, 0.01);
    /// assert_eq!(ratios.1, 0.01);
    /// assert_eq!(ratios.2, 0.0);
    /// ```
    fn get_ratios(&self) -> (f32, f32, f32) {
        let whole_ratio = (self.real_fusion_count + self.fake_fusion_count) / self.totals;
        let real_ratio = self.real_fusion_count / self.totals;
        let fake_ratio = self.fake_fusion_count / self.totals;

        (whole_ratio, real_ratio, fake_ratio)
    }

    /// Gets the ratio of real fusions
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = LocalCounter::new(100);
    ///
    /// counter.inc_real();
    ///
    /// let ratio = counter.get_real_ratio();
    /// assert_eq!(ratio, 0.01);
    /// ```
    fn get_real_ratio(&self) -> f32 {
        self.real_fusion_count / self.totals
    }
}

/// Processes a component of reads
///
/// # Arguments
///
/// * `component` - A tuple containing a RefGenePred and a vector of GenePred
/// * `banned` - A HashSet of banned genes
/// * `recover` - A boolean indicating if the component should be recovered
/// * `match_type` - A MatchType enum
///
/// # Returns
///
/// * An Option containing a tuple of fusions, no_fusions, fake_fusions, review, descriptor, and is_dirty
///
/// # Example
///
/// ```rust, no_run
/// use iso_fusion::core::process_component;
///
/// let component = (
///    RefGenePred::default(),
///   vec![GenePred::default()]
/// );
///
/// let banned = HashSet::new();
/// let recover = false;
/// let match_type = MatchType::SpliceSite;
///
/// process_component(
///     &component,
///     &banned,
///     recover,
///     match_type
/// );
/// ```
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
    // INFO: skipping early if no reads in component
    if component.1.is_empty() {
        return None;
    }

    let mut descriptor = HashMap::new();

    let mut fusions = Vec::new();
    let mut fake_fusions = Vec::new();
    let mut no_fusions = Vec::new();

    let refs = &component.0;
    let queries = &component.1;

    let genes = refs.get_names_split();

    let mut counter = LocalCounter::new(queries.len());

    // INFO: if genes is empty the follwing occurs
    // INFO: species-specific gene | intergenic region | missing gene in refs
    // INFO: if genes len = 1, we have a single gene in the component
    if genes.len() > 1 {
        identify_fusions(
            refs,
            &queries,
            &mut descriptor,
            &mut counter,
            &mut no_fusions,
            &mut fusions,
            &mut fake_fusions,
            banned,
            match_type,
        );
    }

    // INFO: second pass for seen reads
    // INFO: first pass for non-fusions
    fill_schema(
        &queries,
        &mut descriptor,
        &mut no_fusions,
        genes.len(),
        &mut counter,
    );

    if recover {
        // INFO: if the fusion ratio in the component is above the threshold,
        // INFO: mark all queries as dirty and submit them for revie
        if counter.get_real_ratio() >= FUSION_RATIO_THRESHOLD {
            let review = recover_component(&queries, &mut descriptor);
            return Some((vec![], vec![], vec![], Some(review), descriptor, true));
        }
    }

    // dbg!(&descriptor);

    Some((fusions, no_fusions, fake_fusions, None, descriptor, false))
}

/// Recover reads from a fusion component
///
/// # Arguments
///
/// * `queries` - A vector of GenePred structs
/// * `descriptor` - A mutable reference to a HashMap containing the descriptor
///
/// # Returns
///
/// A vector of strings containing the review
///
/// # Example
///
/// ```
/// let queries = vec![GenePred { name: "query1".to_string(), line: "line1".to_string() }];
/// let mut descriptor = HashMap::new();
/// let review = recover_component(&queries, &mut descriptor);
/// assert_eq!(review, vec!["line1"]);
/// ```
fn recover_component(
    queries: &Vec<GenePred>,
    descriptor: &mut HashMap<String, Box<dyn ModuleMap>>,
) -> Vec<String> {
    let mut review = vec![];

    for query in queries.iter() {
        review.push(query.line.clone());
        let handle = descriptor.get_mut(&query.name).unwrap();

        handle
            .set_value(
                Box::new(FusionDetectionValue::IsDirtyFusionComponent),
                Value::Bool(true),
            )
            .ok();
    }

    review
}

/// Fills the schema for seen and unseen reads
///
/// # Arguments
///
/// * `reads` - A vector of GenePred structs
/// * `descriptor` - A mutable reference to a HashMap containing the descriptor
/// * `no_fusions` - A mutable reference to a vector of strings containing no fusions
/// * `genes` - The number of genes
/// * `counter` - A mutable reference to a LocalCounter struct
///
/// # Example
///
/// ```rust, no_run
/// let reads = vec![GenePred::default()];
///
/// let mut descriptor = HashMap::new();
/// let mut no_fusions = vec![];
/// let genes = 1;
/// let mut counter = LocalCounter::new(100);
///
/// fill_schema(&reads, &mut descriptor, &mut no_fusions, genes, &mut counter);
///
/// assert_eq!(no_fusions.len(), 1);
/// assert_eq!(descriptor.len(), 1);
/// ```
fn fill_schema(
    reads: &Vec<GenePred>,
    descriptor: &mut HashMap<String, Box<dyn ModuleMap>>,
    no_fusions: &mut Vec<String>,
    genes: usize,
    counter: &mut LocalCounter,
) {
    reads.iter().for_each(|read| {
        let (whole, real, fake) = counter.get_ratios();

        // INFO: read already exists in the descriptor -> fill descriptor directly
        if descriptor.contains_key(&read.name) {
            let handle = descriptor
                .get_mut(&read.name)
                .expect("ERROR: Failed to get descriptor, this is a bug!");

            handle
                .set_value(
                    Box::new(FusionDetectionValue::WholeComponentFusionRatio),
                    serde_json::json!(whole),
                )
                .ok();
            handle
                .set_value(
                    Box::new(FusionDetectionValue::RealComponentFusionRatio),
                    serde_json::json!(real),
                )
                .ok();
            handle
                .set_value(
                    Box::new(FusionDetectionValue::FakeComponentFusionRatio),
                    serde_json::json!(fake),
                )
                .ok();
        } else {
            // INFO: read does not exist in the descriptor
            // INFO: and it is not a fusion
            no_fusions.push(read.line.clone());

            let mut schema = FusionSchema::default();

            schema.is_fused_read = Value::Bool(false);
            schema.ref_component_size = Value::Number(genes.into());
            schema.query_component_size = Value::Number(reads.len().into());
            schema.whole_component_fusion_ratio = serde_json::json!(whole);
            schema.real_component_fusion_ratio = serde_json::json!(real);
            schema.fake_component_fusion_ratio = serde_json::json!(fake);
            schema.is_dirty_component = Value::Bool(false);

            schema.diffuse(descriptor, read);
        }
    });
}

/// Identifies fusions in the reads
///
/// # Arguments
///
/// * `refs` - A reference to a RefGenePred struct
/// * `reads` - A vector of GenePred structs
/// * `descriptor` - A mutable reference to a HashMap containing the descriptor
/// * `counter` - A mutable reference to a LocalCounter struct
/// * `no_fusions` - A mutable reference to a vector of strings containing no fusions
/// * `fusions` - A mutable reference to a vector of strings containing fusions
/// * `fake_fusions` - A mutable reference to a vector of strings containing fake fusions
/// * `banned` - A HashSet of banned genes
/// * `match_type` - A MatchType enum
///
/// # Example
///
/// ```rust, no_run
/// let refs = RefGenePred::default();
/// let reads = vec![GenePred::default()];
///
/// let mut descriptor = HashMap::new();
/// let mut counter = LocalCounter::new(100);
/// let mut no_fusions = vec![];
/// let mut fusions = vec![];
/// let mut fake_fusions = vec![];
/// let banned = HashSet::new();
/// let match_type = MatchType::Exact;
///
/// identify_fusions(
///     &refs,
///     &reads,
///     &mut descriptor,
///     &mut counter,
///     &mut no_fusions,
///     &mut fusions,
///     &mut fake_fusions,
///     &banned,
///     match_type
/// );
///
/// assert_eq!(no_fusions.len(), 1);
/// assert_eq!(fusions.len(), 0);
/// assert_eq!(fake_fusions.len(), 0);
/// ```
fn identify_fusions(
    refs: &RefGenePred,
    reads: &Vec<GenePred>,
    descriptor: &mut HashMap<String, Box<dyn ModuleMap>>,
    counter: &mut LocalCounter,
    no_fusions: &mut Vec<String>,
    fusions: &mut Vec<String>,
    fake_fusions: &mut Vec<String>,
    banned: &HashSet<String>,
    match_type: MatchType,
) {
    // INFO: fusion loci [more than one gene in component]
    // INFO: we create a per-gene collection of exons
    let ref_exons = refs.smash_exons_by_name();
    let ref_introns = refs.smash_introns_by_name();

    reads.iter().for_each(|query| {
        if banned.contains(&query.name) {
            no_fusions.push(query.line.clone());
            return;
        }

        let mut schema = FusionSchema::default();
        let mut count = 0;

        // INFO: check if query read overlaps any of the gene exon collections
        // WARN: we are not keeping track of gene names here
        for r_exons in ref_exons.iter() {
            if exonic_overlap(&query.exons, r_exons) {
                count += 1;
            }
        }

        schema.is_dirty_component = Value::Bool(false);

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

                // counterreal_fusion_count += 1.0;
                counter.inc_real();
                fusions.push(query.line.clone());

                schema.is_fused_read = Value::Bool(true);
                schema.fusion_in_frame = Value::Bool((query.cds_end - query.cds_start) % 3 == 0);
                schema.location_of_fusion = Value::String(location);
            } else {
                counter.inc_fake();
                // fake_fusion_count += 1.0;
                fake_fusions.push(query.line.clone());
                schema.is_fused_read = Value::Bool(false);
            }
        } else {
            no_fusions.push(query.line.clone());
            schema.is_fused_read = Value::Bool(false);
        }

        schema.diffuse(descriptor, query);
    });
}

/// Detect fusions using mapping data
///
/// # Arguments
///
/// * `args` - A struct containing the command line arguments
///
/// # Example
///
/// ```rust, no_run
/// use iso_fusion::cli::Args;
/// use iso_fusion::core::detect_fusions_with_mapping;
///
/// let args = Args {
///     refs: vec!["path/to/isoforms.tsv".to_string()],
///     query: "path/to/mapping.bed".to_string(),
///     overlap_type: "exon".to_string(),
/// };
///
/// detect_fusions_with_mapping(args);
/// ```
pub fn detect_fusions_with_mapping(args: Args) -> Result<()> {
    let contents = par_reader(args.refs).expect("ERROR: Failed to read isoform file(s)!");
    let refs = tsv_to_map::<IsoformParser, String>(Arc::new(contents), 1, 0) // INFO: gene/ttranscript to transcript->gene Hash
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

    let match_type = MatchType::from(args.intron_match);

    buckets.par_iter().for_each(|bucket| {
        let _ = bucket.key();
        let components = bucket.value().to_owned();
        counter.inc_comp(components.len() as u32);

        components.into_par_iter().for_each(|comp| {
            let comp = comp
                .as_any()
                .downcast_ref::<(Vec<GenePred>, Vec<GenePred>)>()
                .expect("ERROR: Failed to downcast to RefGenePred");

            let (fusions, no_fusions) = process_mapping_component(comp, &refs, match_type);

            acc.add_passes(no_fusions);

            if let Some(f) = fusions {
                acc.add_fusions(f);
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

/// Receives a component of mapping data and returns fusions and non-fusions
///
/// # Arguments
///
/// * `component` - A tuple containing two vectors of GenePred structs
/// * `isoforms` - A HashMap containing isoform data
/// * `match_type` - A MatchType enum
///
/// # Example
///
/// ```rust, no_run
/// use iso_fusion::core::process_mapping_component;
/// use iso_fusion::config::MatchType;
/// use std::collections::HashMap;
/// use packbed::GenePred;
///
/// let component = (
///    vec![GenePred::default()],
///    vec![GenePred::default()]
/// );
///
/// let isoforms = HashMap::new();
///
/// let match_type = MatchType::Exact;
///
/// process_mapping_component(
///    &component,
///    &isoforms,
///    match_type
/// );
/// ```
fn process_mapping_component(
    component: &(Vec<GenePred>, Vec<GenePred>),
    isoforms: &HashMap<String, Vec<String>>,
    match_type: MatchType,
) -> (Option<Vec<String>>, Vec<String>) {
    let query: &Vec<GenePred> = component.0.as_ref(); // INFO: treating refs back to queries, discarding fake queries

    let mut genes = HashMap::new();

    query.iter().for_each(|q| {
        // INFO: getting gene name from transcript name
        let tx_to_gene = isoforms
            .get(&q.name)
            .expect("ERROR: Failed to get isoforms!")
            .get(0) // WARN: should be safe to get the first element
            .expect("ERROR: Failed to get gene name!");

        // INFO: fills a hashmap with genes as keys and exons as values
        let exons = genes.entry(tx_to_gene).or_insert_with(BTreeSet::new);
        exons.extend(q.exons.iter().cloned());
    });

    return get_fusions_from_component(query, genes, match_type);
}

/// Single-component handler for fusions
///
/// # Arguments
///
/// * `component` - A vector of GenePred structs
/// * `genes` - A HashMap containing gene names and exons
/// * `match_type` - A MatchType enum
///
/// # Example
///
/// ```rust, no_run
/// use iso_fusion::core::get_fusions_from_component;
/// use iso_fusion::config::MatchType;
/// use std::collections::{BTreeSet, HashMap};
/// use packbed::GenePred;
///
/// let component = vec![GenePred::default()];
///
/// let genes = HashMap::new();
///
/// let match_type = MatchType::Exact;
///
/// get_fusions_from_component(
///   &component,
///   genes,
///   match_type
/// );
/// ```
fn get_fusions_from_component(
    component: &Vec<GenePred>,
    genes: HashMap<&String, BTreeSet<(u64, u64)>>,
    match_type: MatchType,
) -> (Option<Vec<String>>, Vec<String>) {
    // INFO: if at any point names is > 1, we have a fusion loci
    if genes.len() > 1 {
        let mut fusions = vec![];
        let mut no_fusions = vec![];

        let ref_introns = get_intron_coords_from_gene_map(&genes);

        component.iter().for_each(|query| {
            let mut count = 0;

            for (_, exons) in genes.iter() {
                if exonic_overlap(&query.exons, exons) {
                    count += 1;
                }
            }

            // INFO: query read overlaps more than one gene exon collection
            // INFO: we need to pass again to prove that query shares splice sites
            // INFO: with more than one gene
            if count > 1 {
                let mut splicing_overlaps = 0_f32;
                for r_introns in ref_introns.iter() {
                    if splice_site_overlap(&query.introns, r_introns, match_type) {
                        splicing_overlaps += 1.0;
                    }
                }

                if splicing_overlaps > 1.0 {
                    fusions.push(query.line.clone());
                }
            } else {
                no_fusions.push(query.line.clone());
            }
        });

        return (Some(fusions), no_fusions);
    }

    return (None, component.iter().map(|q| q.line.clone()).collect());
}

/// Get intron coordinates from gene map
///
/// # Arguments
///
/// * `genes` - A HashMap containing gene names and exons
///
/// # Example
///
/// ```rust, no_run
/// use iso_fusion::core::get_intron_coords_from_gene_map;
/// use std::collections::HashMap;
///
/// let genes = HashMap::new();
///
/// get_intron_coords_from_gene_map(&genes);
/// ```
fn get_intron_coords_from_gene_map(
    genes: &HashMap<&String, BTreeSet<(u64, u64)>>,
) -> Vec<HashSet<(u64, u64)>> {
    let mut introns = vec![];

    for (_, exons) in genes.iter() {
        let mut intron = HashSet::new();
        let mut exons = exons.iter().peekable();

        while let Some(exon) = exons.next() {
            if let Some(next_exon) = exons.peek() {
                intron.insert((exon.1 + 1, next_exon.0 - 1));
            }
        }

        introns.push(intron);
    }

    introns
}
