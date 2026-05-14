// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

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
use dashmap::DashMap;
use genepred::{Bed12, GenePred};
use log::info;
use packbed::{pack, OverlapType, Role};
use rayon::prelude::*;

use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::PathBuf;

use crate::cli::Args;
use crate::utils::*;

pub const FUSIONS: &str = "fusions.bed";
pub const FUSION_FREE: &str = "fusions.free.bed";
pub const FUSION_REVIEW: &str = "fusions.review.bed";
pub const FUSION_FAKES: &str = "fusions.fakes.bed";
pub const FUSION_DESCRIPTOR: &str = "fusions.tsv";

// tag separators
pub const SEP: &str = "#"; // WARN: should change to ':' -> breaks ORF caller [translationai hdf5]
pub const BIG_SEP: &str = "__"; // WARN: should change to '::' -> breaks ORF caller [translationai hdf5]

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
pub fn detect_fusions(args: Args) {
    info!("Preparing files for fusion detection...");

    let (buckets, record_to_parent, parent_info) =
        make_buckets(args.refs, args.query).expect("ERROR: Failed to make buckets!");
    let match_type = MatchType::from(args.intron_match);

    info!("Detected {} reference transcripts", record_to_parent.len());

    let counter = ParallelCounter::default();
    let accumulator = ParallelAccumulator::default();

    buckets.into_par_iter().for_each(|(_chr, components)| {
        counter.inc_comp(components.len() as u32);

        process_components(
            components,
            &record_to_parent,
            &parent_info,
            &accumulator,
            &counter,
            args.recover,
            match_type,
            args.threshold,
            args.tag,
        );
    });

    info!("Detected fusions: {}", accumulator.fusions.len());
    info!("Fusion-free reads: {}", accumulator.passes.len());

    if args.recover {
        log::warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            counter.num_of_dirty,
            counter.load_ratio()
        );
    }

    if args.descriptor {
        write_descriptor(&accumulator.descriptor, args.prefix.join(FUSION_DESCRIPTOR));
    }

    par_write_results(
        &accumulator,
        vec![
            args.prefix.join(FUSIONS),
            args.prefix.join(FUSION_FREE),
            args.prefix.join(FUSION_REVIEW),
        ],
        None,
    );
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
    mut refs: Vec<PathBuf>,
    query: Vec<PathBuf>,
) -> Result<(
    DashMap<String, Vec<Vec<GenePred>>>,
    DashMap<Vec<u8>, Vec<u8>>,
    DashMap<Vec<u8>, GroupData>,
)> {
    let (record_to_parent, parent_info) = create_reference_map(refs.clone())?;
    log::info!(
        "Created reference map with {} entries",
        record_to_parent.len()
    );
    log::info!("Created parent info with {} entries", parent_info.len());

    let mut modes = Vec::new();
    modes.extend(std::iter::repeat(Role::Reference).take(refs.len()));
    modes.extend(std::iter::repeat(Role::Query).take(query.len()));

    refs.extend(query);
    let tracks = pack(refs, modes, OverlapType::Exon)
        .unwrap_or_else(|_| panic!("ERROR: Failed to pack reference transcripts!"));

    Ok((tracks, record_to_parent, parent_info))
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
#[allow(clippy::too_many_arguments)]
fn process_components(
    components: Vec<Vec<GenePred>>,
    record_to_parent: &DashMap<Vec<u8>, Vec<u8>>,
    parent_info: &DashMap<Vec<u8>, GroupData>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
    recover: bool,
    match_type: MatchType,
    threshold: f32,
    tag: bool,
) {
    components.into_iter().for_each(|comp| {
        let (reference_regions, query_regions): (Vec<_>, Vec<_>) =
            comp.into_iter().partition(|region| {
                let role = region
                    .get_extra("role".as_bytes())
                    .unwrap_or_else(|| panic!("ERROR: role not found for region: {region:?}"))
                    .clone()
                    .into_scalar()
                    .unwrap_or_else(|| panic!("ERROR: role not found for region: {region:?}"));

                matches!(role.as_slice(), b"reference")
            });

        let (fusions, no_fusions, review, descriptor, is_dirty) = process_component(
            reference_regions,
            query_regions,
            record_to_parent,
            parent_info,
            recover,
            match_type,
            threshold,
            tag,
        )
        .unwrap_or_default();

        accumulator.add(fusions, no_fusions, review, descriptor);

        if is_dirty {
            counter.inc_dirty(1);
        }
    });
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
    reference_regions: Vec<GenePred>,
    mut query_regions: Vec<GenePred>,
    record_to_parent: &DashMap<Vec<u8>, Vec<u8>>,
    parent_info: &DashMap<Vec<u8>, GroupData>,
    recover: bool,
    match_type: MatchType,
    threshold: f32,
    tag: bool,
) -> Option<(
    Vec<String>,
    Vec<String>,
    Option<Vec<String>>,
    HashMap<String, FusionSchema>,
    bool,
)> {
    // INFO: skipping early if no reads in component
    if query_regions.is_empty() {
        return None;
    }

    let mut descriptor = HashMap::new();

    let mut fusions = Vec::new();
    let mut no_fusions = Vec::new();

    let mut parents = HashSet::new();
    for r in reference_regions.iter() {
        let parent = record_to_parent.get(r.name().unwrap()).unwrap_or_else(|| {
            panic!(
                "ERROR: Failed to get record_to_parent for record: {:?}",
                r.name()
            )
        });

        if !parents.contains(parent.value()) {
            parents.insert(parent.clone());
        }
    }

    let mut counter = LocalCounter::new(query_regions.len());

    // INFO: if genes is empty the follwing occurs
    // INFO: species-specific gene | intergenic region | missing gene in refs
    // INFO: if genes len = 1, we have a single gene in the component
    if parents.len() > 1 {
        let query_len = query_regions.len();
        identify_fusions(
            &mut query_regions,
            query_len,
            &parents,
            parent_info,
            &mut descriptor,
            &mut counter,
            &mut no_fusions,
            &mut fusions,
            match_type,
            tag,
        );
    }

    // INFO: second pass for seen reads
    // INFO: first pass for non-fusions
    fill_schema(
        &query_regions,
        &mut descriptor,
        &mut no_fusions,
        parents.len(),
        &mut counter,
    );

    if recover {
        // INFO: if the fusion ratio in the component is above the threshold,
        // INFO: mark all queries as dirty and submit them for review
        if counter.get_real_ratio() >= threshold {
            let review = recover_component(&mut query_regions, &mut descriptor);
            return Some((vec![], vec![], Some(review), descriptor, true));
        }
    }

    Some((fusions, no_fusions, None, descriptor, false))
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
    queries: &mut [GenePred],
    descriptor: &mut HashMap<String, FusionSchema>,
) -> Vec<String> {
    let mut review = vec![];

    for query in queries.iter_mut() {
        // INFO: append :RW tag to read name
        let query_name =
            unsafe { std::str::from_utf8_unchecked(query.name().unwrap()) }.to_string();
        query.set_name(Some(format!("{}{SEP}RW", query_name).as_bytes().to_vec()));
        let query_line =
            unsafe { std::str::from_utf8_unchecked(&query.to_bed::<Bed12>()) }.to_string();

        review.push(query_line);
        let handle = descriptor.get_mut(&query_name).unwrap_or_else(|| {
            panic!("ERROR: Failed to get descriptor for read: {query_name:?}, this is a bug!")
        });

        handle.is_dirty_component = true;
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
    reads: &[GenePred],
    descriptor: &mut HashMap<String, FusionSchema>,
    no_fusions: &mut Vec<String>,
    genes: usize,
    counter: &mut LocalCounter,
) {
    reads.iter().for_each(|read| {
        let (whole, real, fake) = counter.get_ratios();
        let read_name = unsafe { std::str::from_utf8_unchecked(read.name().unwrap()) }.to_string();

        // INFO: read already exists in the descriptor -> fill descriptor directly
        if descriptor.contains_key(&read_name) {
            let handle = descriptor.get_mut(&read_name).unwrap_or_else(|| {
                panic!("ERROR: Failed to get descriptor for read: {read_name:?}, this is a bug!")
            });

            handle.whole_component_fusion_ratio = whole;
            handle.real_component_fusion_ratio = real;
            handle.fake_component_fusion_ratio = fake;
        } else {
            // INFO: read does not exist in the descriptor
            // INFO: and it is not a fusion
            let read_line =
                unsafe { std::str::from_utf8_unchecked(&read.to_bed::<Bed12>()) }.to_string();
            no_fusions.push(read_line);

            let mut schema = FusionSchema::default();

            schema.is_fused_read = false;
            schema.ref_component_size = genes as u32;
            schema.query_component_size = reads.len() as u32;
            schema.whole_component_fusion_ratio = whole;
            schema.real_component_fusion_ratio = real;
            schema.fake_component_fusion_ratio = fake;
            schema.is_dirty_component = false;

            descriptor.insert(read_name, schema);
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
#[allow(clippy::too_many_arguments)]
fn identify_fusions(
    query_regions: &mut Vec<GenePred>,
    query_len: usize,
    parents: &HashSet<Vec<u8>>,
    parent_info: &DashMap<Vec<u8>, GroupData>,
    descriptor: &mut HashMap<String, FusionSchema>,
    counter: &mut LocalCounter,
    no_fusions: &mut Vec<String>,
    fusions: &mut Vec<String>,
    match_type: MatchType,
    _tag: bool,
) {
    let reference_exons: Vec<BTreeSet<(u64, u64)>> = parents
        .iter()
        .map(|p| parent_info.get(p).unwrap().exons.clone())
        .collect();
    let reference_introns: Vec<BTreeSet<(u64, u64)>> = parents
        .iter()
        .map(|p| parent_info.get(p).unwrap().introns.clone())
        .collect();

    query_regions.iter_mut().for_each(|query| {
        let mut query_name =
            unsafe { std::str::from_utf8_unchecked(query.name().unwrap()) }.to_string();
        let query_exons = query.exons().clone();
        let query_introns = query.introns().clone();

        let mut schema = FusionSchema::default();
        let mut count = 0;

        // INFO: check if query read overlaps any of the gene exon collections
        // WARN: we are not keeping track of gene names here
        for r_exons in reference_exons.iter() {
            if exonic_overlap(&query_exons, r_exons) {
                count += 1;
            }
        }

        schema.is_dirty_component = false;
        schema.ref_component_size = reference_exons.len() as u32;
        schema.query_component_size = query_len as u32;

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
            for r_introns in reference_introns.iter() {
                if splice_site_overlap(&query_introns, r_introns, match_type) {
                    splicing_overlaps += 1.0;
                }
            }

            if splicing_overlaps > 1.0 {
                // counter.real_fusion_count += 1.0;
                counter.inc_real();
                let query_line =
                    unsafe { std::str::from_utf8_unchecked(&query.to_bed::<Bed12>()) }.to_string();
                fusions.push(query_line);

                schema.is_fused_read = true;
            } else {
                counter.inc_fake();
                // fake_fusion_count += 1.0;
                // INFO: append fake tag to read {:FK}

                // if tag {
                query_name = format!("{}{SEP}FK", query_name);
                query.set_name(Some(query_name.as_bytes().to_vec()));
                // }

                let query_line =
                    unsafe { std::str::from_utf8_unchecked(&query.to_bed::<Bed12>()) }.to_string();

                no_fusions.push(query_line);
                schema.is_fused_read = false;

                descriptor.insert(query_name, schema);
                return;
            }
        } else {
            let query_line =
                unsafe { std::str::from_utf8_unchecked(&query.to_bed::<Bed12>()) }.to_string();
            no_fusions.push(query_line);
            schema.is_fused_read = false;
        }

        descriptor.insert(query_name, schema);
    });
}
