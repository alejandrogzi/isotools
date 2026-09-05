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

use std::{
    borrow::Borrow,
    collections::{BTreeSet, HashMap, HashSet},
    fmt::Debug,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicU32, Ordering},
};

use anyhow::{Ok, Result};
use dashmap::{DashMap, DashSet};
use genepred::GenePred;
use log::info;
use num_traits::{Num, NumCast};
use packbed::{pack, OverlapType, Role};
use petgraph::unionfind::UnionFind;
use rayon::prelude::*;
use std::collections::hash_map::Entry;

/// Prepare reference transcripts for fusions
///
/// # Parameters
///
/// - `refs`: A vector of paths to reference files.
///
/// # Returns
///
/// - A Result containing a String with the prepared reference transcripts.
///
/// # Example
///
/// ```rust, no_run
/// let refs = vec![PathBuf::from("path/to/ref1.bed"), PathBuf::from("path/to/ref2.bed")];
/// let result = prepare_refs(refs);
///
/// assert!(result.is_ok());
/// ```
pub fn create_reference_map<P: AsRef<Path> + Debug + Sync + Send>(
    refs: Vec<P>,
) -> Result<(DashMap<Vec<u8>, Vec<u8>>, DashMap<Vec<u8>, GroupData>), anyhow::Error> {
    let mut modes = Vec::new();
    modes.extend(std::iter::repeat(Role::Reference).take(refs.len()));
    let tracks = pack(refs, modes, OverlapType::Exon)
        .unwrap_or_else(|e| panic!("ERROR: Failed to pack reference transcripts {:?}", e));

    let record_to_parent = DashMap::new();
    let parent_info = DashMap::new();

    tracks.into_par_iter().for_each(|(_chrom, components)| {
        let (local_record_to_parent, local_parent_info) = group_components(&components);

        local_record_to_parent.into_iter().for_each(|(k, v)| {
            record_to_parent.insert(k, v);
        });

        local_parent_info.into_iter().for_each(|(k, v)| {
            parent_info.insert(k, v);
        });
    });

    info!("Reference transcript map created!");
    Ok((record_to_parent, parent_info))
}

/// Parallel accumulator for the processing function
///
/// # Fields
///
/// - `fusions`: A set of fusion names.
/// - `review`: A set of review names.
/// - `passes`: A set of pass names.
/// - `fakes`: A set of fake names.
/// - `descriptor`: A map of strings to boxed `ModuleMap` trait objects.
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.fusions.len(), 0);
/// ```
pub struct ParallelAccumulator {
    pub fusions: DashSet<String>,
    pub review: DashSet<String>,
    pub passes: DashSet<String>,
    pub descriptor: DashMap<String, FusionSchema>,
}

/// ParallelAccumulator constructor
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.fusions.len(), 0);
/// assert_eq!(accumulator.review.len(), 0);
/// assert_eq!(accumulator.passes.len(), 0);
/// assert_eq!(accumulator.fakes.len(), 0);
/// assert_eq!(accumulator.descriptor.len(), 0);
/// ```
impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            fusions: DashSet::new(),
            review: DashSet::new(),
            passes: DashSet::new(),
            descriptor: DashMap::new(),
        }
    }
}

pub trait ParallelCollector {
    fn len(&self) -> usize;
    fn get_collections(&self) -> Result<Vec<&DashSet<String>>, Box<dyn std::error::Error>>;
}

/// ParallelCollector trait for ParallelAccumulator
impl ParallelCollector for ParallelAccumulator {
    /// Get the number of fields in the accumulator
    fn len(&self) -> usize {
        ParallelAccumulator::NUM_FIELDS
    }

    /// Get the a collection of items from the accumulator
    fn get_collections(&self) -> Result<Vec<&DashSet<String>>, Box<dyn std::error::Error>> {
        let mut collections = Vec::with_capacity(ParallelAccumulator::NUM_FIELDS);

        collections.push(&self.fusions);
        collections.push(&self.passes);
        collections.push(&self.review);

        std::result::Result::Ok(collections)
    }
}

impl ParallelAccumulator {
    /// Number of fields in the accumulator of type DashSet<String>
    pub const NUM_FIELDS: usize = 3;

    /// Add items to the accumulator
    ///
    /// # Parameters
    ///
    /// - `fusions`: A vector of fusion names.
    /// - `passes`: A vector of pass names.
    /// - `review`: An optional vector of review names.
    /// - `fakes`: A vector of fake names.
    /// - `descriptor`: A map of strings to boxed `ModuleMap` trait objects.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let accumulator = ParallelAccumulator::default();
    ///
    /// accumulator.add(
    ///     vec!["fusion1".to_string(), "fusion2".to_string()],
    ///     vec!["pass1".to_string(), "pass2".to_string()],
    ///     Some(vec!["review1".to_string()]),
    ///     vec!["fake1".to_string(), "fake2".to_string()],
    ///     HashMap::new(),
    /// );
    ///
    /// assert_eq!(accumulator.fusions.len(), 2);
    /// assert_eq!(accumulator.passes.len(), 2);
    /// assert_eq!(accumulator.review.len(), 1);
    /// assert_eq!(accumulator.fakes.len(), 2);
    /// assert_eq!(accumulator.descriptor.len(), 0);
    /// ```
    pub fn add(
        &self,
        fusions: Vec<String>,
        passes: Vec<String>,
        review: Option<Vec<String>>,
        descriptor: HashMap<String, FusionSchema>,
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

        descriptor.into_iter().for_each(|(k, v)| {
            self.descriptor.insert(k, v);
        });
    }

    /// Add fusions to the accumulator
    ///
    /// # Parameters
    ///
    /// - `fusions`: A vector of fusion names.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let accumulator = ParallelAccumulator::default();
    /// accumulator.add_fusions(vec!["fusion1".to_string(), "fusion2".to_string()]);
    ///
    /// assert_eq!(accumulator.fusions.len(), 2);
    /// ```
    pub fn add_fusions(&self, fusions: Vec<String>) {
        fusions.into_iter().for_each(|fusion| {
            self.fusions.insert(fusion);
        });
    }

    /// Add review items to the accumulator
    ///
    /// # Parameters
    ///
    /// - `review`: A vector of review names.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let accumulator = ParallelAccumulator::default();
    /// accumulator.add_review(vec!["review1".to_string(), "review2".to_string()]);
    ///
    /// assert_eq!(accumulator.review.len(), 2);
    /// ```
    pub fn add_passes(&self, passes: Vec<String>) {
        passes.into_iter().for_each(|pass| {
            self.passes.insert(pass);
        });
    }
}

/// Parallel counter for the processing function
///
/// # Fields
///
/// - `num_of_comps`: Number of components processed.
/// - `num_of_dirty`: Number of dirty components processed.
///
/// # Example
///
/// ```rust, no_run
/// let counter = ParallelCounter::default();
///
/// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 0);
/// ```
pub struct ParallelCounter {
    pub num_of_comps: AtomicU32,
    pub num_of_dirty: AtomicU32,
}

/// Default implementation for ParallelCounter
///
/// # Example
///
/// ```rust, no_run
/// let counter = ParallelCounter::default();
///
/// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 0);
/// assert_eq!(counter.num_of_dirty.load(Ordering::Relaxed), 0);
/// ```
impl Default for ParallelCounter {
    fn default() -> Self {
        Self {
            num_of_comps: AtomicU32::new(0),
            num_of_dirty: AtomicU32::new(0),
        }
    }
}

impl ParallelCounter {
    /// Increment the number of components
    ///
    /// # Parameters
    ///
    /// - `value`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::default();
    /// counter.inc_comp(5);
    ///
    /// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 5);
    /// ```
    pub fn inc_comp(&self, value: u32) {
        self.num_of_comps.fetch_add(value, Ordering::Relaxed);
    }

    /// Increment the number of dirty components
    ///
    /// # Parameters
    ///
    /// - `value`: The number of dirty components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::default();
    /// counter.inc_dirty(3);
    ///
    /// assert_eq!(counter.num_of_dirty.load(Ordering::Relaxed), 3);
    /// ```
    pub fn inc_dirty(&self, value: u32) {
        self.num_of_dirty.fetch_add(value, Ordering::Relaxed);
    }

    /// Get the component ratio (dirty/total)
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::default();
    /// counter.inc_comp(10);
    /// counter.inc_dirty(3);
    ///
    /// assert_eq!(counter.load_ratio(), 0.3);
    /// ```
    pub fn load_ratio(&self) -> f64 {
        self.num_of_dirty.load(Ordering::Relaxed) as f64
            / self.num_of_comps.load(Ordering::Relaxed) as f64
    }
}

/// Splice site overlap between two sets of introns.
///
/// This function is used to determine if two sets of
/// introns overlap by comparing the start and end
/// positions of each intron. The latter will depend
/// on the match type.
///
/// # Arguments
/// * `introns_a` - a collection of introns
/// * `introns_b` - a collection of introns
/// * `match_type` - the match type [Intron, SpliceSite]
///
/// # Returns
/// * `bool` - true if there is an overlap, false otherwise
///
/// # Example
/// ```rust
/// use isotools::config::fns::splice_site_overlap;
/// use hashbrown::HashSet;
///
/// let introns_a = vec![(1, 10), (20, 30)];
/// let introns_b = [(5, 15), (25, 35)].iter().cloned().collect::<HashSet<_>>();
/// assert_eq!(splice_site_overlap(&introns_a, &introns_b), true);
/// ```
#[inline(always)]
pub fn splice_site_overlap<'a, N, I>(
    introns_a: &[(N, N)],
    introns_b: I,
    match_type: MatchType,
) -> bool
where
    N: Num + NumCast + Copy + PartialOrd + Eq + std::hash::Hash,
    I: Clone + IntoIterator<Item = &'a (N, N)>,
    N: 'a,
{
    match match_type {
        MatchType::Intron => {
            let b: HashSet<(N, N)> = introns_b.clone().into_iter().copied().collect();
            introns_a.iter().any(|a| b.contains(a))
        }
        MatchType::SpliceSite => {
            let splice_sites: HashSet<N> =
                introns_b.into_iter().flat_map(|&(x, y)| [x, y]).collect();

            introns_a
                .iter()
                .any(|&(x, y)| splice_sites.contains(&x) || splice_sites.contains(&y))
        }
    }
}

/// Match type
///
/// This enum is used to store the type of match.
///
/// # Example
///
/// ```rust, no_run
/// use iso::MatchType;
///
/// let intron = MatchType::Intron;
/// let splice_site = MatchType::SpliceSite;
/// ```
#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub enum MatchType {
    Intron,
    #[default]
    SpliceSite,
}

impl From<bool> for MatchType {
    fn from(value: bool) -> Self {
        match value {
            true => MatchType::Intron,
            false => MatchType::SpliceSite,
        }
    }
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
///
pub struct FusionSchema {
    pub is_fused_read: bool,
    pub ref_component_size: u32,
    pub query_component_size: u32,
    pub whole_component_fusion_ratio: f32,
    pub real_component_fusion_ratio: f32,
    pub fake_component_fusion_ratio: f32,
    pub is_dirty_component: bool,
    pub is_fuzzy_masked: bool,
}

impl FusionSchema {}

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
            is_fused_read: false,
            ref_component_size: 0,
            query_component_size: 0,
            whole_component_fusion_ratio: 0.0,
            real_component_fusion_ratio: 0.0,
            fake_component_fusion_ratio: 0.0,
            is_dirty_component: false,
            is_fuzzy_masked: false,
        }
    }
}

impl std::fmt::Display for FusionSchema {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.is_fused_read,
            self.ref_component_size,
            self.query_component_size,
            self.whole_component_fusion_ratio,
            self.real_component_fusion_ratio,
            self.fake_component_fusion_ratio,
            self.is_dirty_component,
            self.is_fuzzy_masked
        )
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
pub struct LocalCounter {
    pub real_fusion_count: f32,
    pub fake_fusion_count: f32,
    pub totals: f32,
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
    pub fn new(totals: usize) -> Self {
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
    pub fn inc_real(&mut self) {
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
    pub fn inc_fake(&mut self) {
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
    pub fn get_ratios(&self) -> (f32, f32, f32) {
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
    pub fn get_real_ratio(&self) -> f32 {
        self.real_fusion_count / self.totals
    }
}

/// Exonic overlap between two sets of exons.
///
/// This function is used to determine if two
/// sets of exons overlap by comparing the start
/// and end positions of each exon. The sets of
/// exons could be any type of collection.
///
/// # Arguments
/// * `exons_a` - a collection of exons
/// * `exons_b` - a collection of exons
///
/// # Returns
/// * `bool` - true if there is an overlap, false otherwise
///
/// # Example
/// ```rust
/// use isotools::config::fns::exonic_overlap;
///
/// let exons_a = vec![(1, 10), (20, 30)];
/// let exons_b = vec![(5, 15), (25, 35)];
/// assert_eq!(exonic_overlap(&exons_a, &exons_b), true);
/// ```
#[inline(always)]
pub fn exonic_overlap<N, I1, I2>(exons_a: &I1, exons_b: &I2) -> bool
where
    N: Num + NumCast + Copy + PartialOrd,
    I1: IntoIterator,
    I2: IntoIterator,
    I1::Item: Borrow<(N, N)>,
    I2::Item: Borrow<(N, N)>,
    for<'a> &'a I1: IntoIterator<Item = &'a I1::Item>,
    for<'a> &'a I2: IntoIterator<Item = &'a I2::Item>,
{
    let mut iter_a = exons_a.into_iter();
    let mut iter_b = exons_b.into_iter();

    let mut exon_a = iter_a.next();
    let mut exon_b = iter_b.next();

    loop {
        match (exon_a, exon_b) {
            (Some(start_end_a), Some(start_end_b)) => {
                let (start_a, end_a) = start_end_a.borrow();
                let (start_b, end_b) = start_end_b.borrow();

                if *start_a < *end_b && *start_b < *end_a {
                    return true;
                }

                if *end_a < *end_b {
                    exon_a = iter_a.next();
                } else {
                    exon_b = iter_b.next();
                }
            }
            _ => break,
        }
    }

    false
}

/// Write results from a ParallelAccumulator to files in parallel
///
/// # Arguments
///
/// * `accumulator` - The accumulator containing the results
/// * `filenames` - A vector of PathBufs representing the filenames
/// * `outdir` - An optional PathBuf representing the output directory
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::par_write_results;
///
/// let accumulator = ParallelAccumulator::default();
/// let filenames = vec![PathBuf::from("file1.txt"), PathBuf::from("file2.txt")];
/// let outdir = Some(PathBuf::from("output"));
///
/// par_write_results(accumulator, filenames, outdir);
/// ```
pub fn par_write_results<K: ParallelCollector>(
    accumulator: &K,
    filenames: Vec<PathBuf>,
    outdir: Option<PathBuf>,
) {
    if accumulator.len() != filenames.len() {
        log::error!("ERROR: Number of filenames does not match the number of results!");
        std::process::exit(1);
    }

    let collections = accumulator
        .get_collections()
        .expect("ERROR: Could not get collections!");

    if let Some(outdir) = outdir {
        let files = filenames
            .iter()
            .map(|filename| outdir.join(filename))
            .collect::<Vec<_>>();

        write_pairs(collections, files);
    } else {
        write_pairs(collections, filenames);
    }
}

/// Inner function to write pairs of collections to files
///
/// # Arguments
///
/// * `collections` - A vector of DashSet<String> representing the collections
/// * `filenames` - A vector of PathBufs representing the filenames
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::write_pairs;
///
/// let collections = vec![DashSet::new(), DashSet::new()];
/// let filenames = vec![PathBuf::from("file1.txt"), PathBuf::from("file2.txt")];
///
/// write_pairs(collections, filenames);
/// ```
fn write_pairs(collections: Vec<&DashSet<String>>, filenames: Vec<PathBuf>) {
    collections
        .par_iter()
        .zip(filenames.par_iter())
        .for_each(|(collection, path)| {
            if collection.is_empty() {
                log::warn!("WARN: A collection from the accumulator is empty! Skipping...");
                return;
            }

            write_objs(
                collection,
                path.to_str().expect("ERROR: Invalid path to write!"),
            );
        });
}

/// Write a DashSet to a file
///
/// # Arguments
///
/// * `data` - The DashSet to write
/// * `fname` - The name of the file to write to
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::write_objs;
///
/// let data = DashSet::new();
/// data.insert("line1".to_string());
/// data.insert("line2".to_string());
///
/// let fname = "output.txt";
///
/// write_objs(&data, fname);
/// ```
pub fn write_objs<T>(data: &DashSet<T>, fname: &str)
where
    T: AsRef<[u8]> + Sync + Send + Eq + std::hash::Hash,
{
    log::info!("Reads in {}: {:?}. Writing...", fname, data.len());
    let f = match File::create(fname) {
        std::result::Result::Ok(f) => f,
        Err(e) => panic!("Error creating file: {}", e),
    };
    let mut writer = BufWriter::new(f);

    for line in data.iter() {
        let bytes = line.as_ref();
        writer
            .write_all(bytes)
            .unwrap_or_else(|e| panic!("ERROR: Error writing to file -> {e}"));
        writer
            .write_all(b"\n")
            .unwrap_or_else(|e| panic!("ERROR: Error newline to file -> {e}"));
    }
}

/// Write a HashMap descriptor to a file
pub fn write_descriptor(data: &DashMap<String, FusionSchema>, path: PathBuf) {
    log::info!("Reads in {}: {:?}. Writing...", path.display(), data.len());
    let f = match File::create(path) {
        std::result::Result::Ok(f) => f,
        Err(e) => panic!("Error creating file: {}", e),
    };
    let mut writer = BufWriter::new(f);

    data.iter().for_each(|entry| {
        let (key, value) = entry.pair();
        let line = format!("{}\t{}", key, value);
        let bytes = line.as_bytes();
        writer
            .write_all(bytes)
            .unwrap_or_else(|e| panic!("ERROR: Error writing to file -> {e}"));
        writer
            .write_all(b"\n")
            .unwrap_or_else(|e| panic!("ERROR: Error newline to file -> {e}"));
    })
}

type Exon = (u64, u64);

/// GroupData struct
///
/// This struct is used to store the exons and introns of a group of records.
#[derive(Debug, Clone, Default)]
pub struct GroupData {
    pub exons: BTreeSet<Exon>,
    pub introns: BTreeSet<Exon>,
    pub cds: BTreeSet<Exon>,
    pub genes: HashSet<Vec<u8>>,
}

/// Parse the `#GENE` tag from a record name (`TX#GENE#N`).
/// Returns the lowercased gene symbol, or `None` if absent.
pub fn parse_gene_tag(name: &[u8]) -> Option<Vec<u8>> {
    let mut parts = name.split(|b| *b == b'#');
    // part[0] = tx id, part[1] = gene symbol
    let _tx = parts.next()?;
    let gene = parts.next()?;
    if gene.is_empty() {
        return None;
    }
    let mut lower = gene.to_vec();
    lower.make_ascii_lowercase();
    Some(lower)
}

#[inline]
fn intervals_overlap(a: (u64, u64), b: (u64, u64)) -> bool {
    a.0 < b.1 && b.0 < a.1
}

fn cds_overlap(a: &BTreeSet<Exon>, b: &BTreeSet<Exon>) -> bool {
    for x in a.iter() {
        for y in b.iter() {
            if intervals_overlap(*x, *y) {
                return true;
            }
        }
    }
    false
}

fn shared_splice_sites(a: &BTreeSet<Exon>, b: &BTreeSet<Exon>) -> usize {
    use std::collections::HashSet as StdHashSet;
    let mut sa = StdHashSet::new();
    for (s, e) in a.iter() {
        sa.insert(*s);
        sa.insert(*e);
    }
    let mut shared = StdHashSet::new();
    for (s, e) in b.iter() {
        if sa.contains(s) {
            shared.insert(*s);
        }
        if sa.contains(e) {
            shared.insert(*e);
        }
    }
    shared.len()
}

/// Build the set of non-fusible (masked) parent pairs.
///
/// A pair is masked when both parents look like the same locus split by
/// heterogeneous annotation (e.g. truncated CDS / liftover shifts):
/// - same `#GENE` tag (case-insensitive), OR
/// - >=1 exact shared exon, OR
/// - >=1 exact shared intron (full intron match), OR
/// - CDS overlap (any bp) + >=1 shared splice site.
///
/// Pairs are stored canonically ordered; lookup must use
/// [`masked_pair`] / [`is_masked`].
pub fn build_exclusion_mask(
    parents: &HashSet<Vec<u8>>,
    parent_info: &DashMap<Vec<u8>, GroupData>,
) -> HashSet<(Vec<u8>, Vec<u8>)> {
    let list: Vec<Vec<u8>> = parents.iter().cloned().collect();
    let mut mask = HashSet::new();
    for i in 0..list.len() {
        for j in (i + 1)..list.len() {
            let (a, b) = (&list[i], &list[j]);
            let (da, db) = match (parent_info.get(a), parent_info.get(b)) {
                (Some(x), Some(y)) => (x, y),
                _ => continue,
            };
            let same_gene = !da.genes.is_empty()
                && !db.genes.is_empty()
                && da.genes.intersection(&db.genes).next().is_some();
            let exact_exon = !da.exons.is_disjoint(&db.exons);
            let full_intron = !da.introns.is_disjoint(&db.introns);
            let fuzzy_cds =
                cds_overlap(&da.cds, &db.cds) && shared_splice_sites(&da.introns, &db.introns) >= 1;
            if same_gene || exact_exon || full_intron || fuzzy_cds {
                let (x, y) = if a <= b {
                    (a.clone(), b.clone())
                } else {
                    (b.clone(), a.clone())
                };
                mask.insert((x, y));
            }
        }
    }
    mask
}

#[inline]
fn canonical_pair(a: &[u8], b: &[u8]) -> (Vec<u8>, Vec<u8>) {
    if a <= b {
        (a.to_vec(), b.to_vec())
    } else {
        (b.to_vec(), a.to_vec())
    }
}

/// Check whether two parents are masked as non-fusible.
pub fn is_masked(mask: &HashSet<(Vec<u8>, Vec<u8>)>, a: &[u8], b: &[u8]) -> bool {
    mask.contains(&canonical_pair(a, b))
}

/// Group components into records and their exons/introns by their parent.
///
/// # Arguments
///
/// * `components` - A vector of vectors of GenePred structs
///
/// # Returns
///
/// * A tuple containing a DashMap of record names to parent names and a DashMap of parent names to GroupData structs
///
/// # Example
///
/// ```rust, no_run
/// use isotools::utils::group_components;
///
/// let components = vec![vec![GenePred::default()]];
/// let (record_to_parent, parent_info) = group_components(&components);
///
/// assert_eq!(record_to_parent.len(), 1);
/// assert_eq!(parent_info.len(), 1);
/// ```
pub fn group_components(
    components: &[Vec<GenePred>],
) -> (DashMap<Vec<u8>, Vec<u8>>, DashMap<Vec<u8>, GroupData>) {
    let record_to_parent: DashMap<Vec<u8>, Vec<u8>> = DashMap::new();
    let parent_info: DashMap<Vec<u8>, GroupData> = DashMap::new();

    components
        .par_iter()
        .enumerate()
        .for_each(|(comp_idx, component)| {
            if component.is_empty() {
                return;
            }
            let n = component.len();
            let mut strand = '0';

            // Cache CDS once — coding_exons() likely allocates.
            let coding: Vec<Vec<Exon>> = component
                .iter()
                .map(|r| {
                    if strand == '0' {
                        strand = match r.strand() {
                            Some(genepred::Strand::Reverse) => '-',
                            Some(genepred::Strand::Forward) => '+',
                            _ => panic!("ERROR: Unexpected strand: {:?}", r.strand()),
                        }
                    }
                    r.coding_exons()
                })
                .collect();

            // 1. Union-Find: for each CDS exon, union every record that contains it.
            //    We only need to remember the *first* record that ever saw an exon;
            //    later records get unioned to it. Transitivity falls out automatically.
            let mut uf = UnionFind::<usize>::new(n);
            let mut exon_to_first: HashMap<(char, Exon), usize> = HashMap::new();

            for (i, cds) in coding.iter().enumerate() {
                for &exon in cds {
                    match exon_to_first.entry((strand, exon)) {
                        Entry::Occupied(e) => {
                            uf.union(*e.get(), i);
                        }
                        Entry::Vacant(v) => {
                            v.insert(i);
                        }
                    }
                }
            }

            // 2. Bucket records by their UF root. Skip non-coding records.
            let mut groups: HashMap<usize, Vec<usize>> = HashMap::with_capacity(8);
            for i in 0..n {
                if coding[i].is_empty() {
                    continue;
                }
                groups.entry(uf.find_mut(i)).or_default().push(i);
            }
            if groups.is_empty() {
                return;
            }

            // 3. Deterministic sub_idx: sort groups by their lowest member index.
            let chrom = &component[0].chrom;
            let mut entries: Vec<Vec<usize>> = groups.into_values().collect();
            for g in entries.iter_mut() {
                g.sort_unstable();
            }
            entries.sort_unstable_by_key(|g| g[0]);

            // 4. Aggregate exons/introns per group, write into the DashMaps.
            for (sub_idx, idxs) in entries.into_iter().enumerate() {
                let mut group_name = Vec::with_capacity(chrom.len() + 24);
                group_name.extend_from_slice(b"COMP_");
                group_name.extend_from_slice(chrom);
                write!(&mut group_name, "_{}_{}_{}", strand, comp_idx, sub_idx).unwrap();

                let mut exons: BTreeSet<Exon> = BTreeSet::new();
                let mut introns: BTreeSet<Exon> = BTreeSet::new();
                let mut cds: BTreeSet<Exon> = BTreeSet::new();
                let mut genes: HashSet<Vec<u8>> = HashSet::new();

                for &i in &idxs {
                    let rec = &component[i];
                    exons.extend(rec.exons());
                    introns.extend(rec.introns());
                    cds.extend(coding[i].iter().copied());
                    if let Some(name) = &rec.name {
                        record_to_parent.insert(name.clone(), group_name.clone());
                        if let Some(g) = parse_gene_tag(name) {
                            genes.insert(g);
                        }
                    }
                }

                parent_info.insert(
                    group_name,
                    GroupData {
                        exons,
                        introns,
                        cds,
                        genes,
                    },
                );
            }
        });

    (record_to_parent, parent_info)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn group(exons: &[(u64, u64)], introns: &[(u64, u64)], cds: &[(u64, u64)], genes: &[&str]) -> GroupData {
        GroupData {
            exons: exons.iter().copied().collect(),
            introns: introns.iter().copied().collect(),
            cds: cds.iter().copied().collect(),
            genes: genes
                .iter()
                .map(|g| g.as_bytes().to_ascii_lowercase())
                .collect(),
        }
    }

    #[test]
    fn gene_tag_is_case_insensitive() {
        assert_eq!(
            parse_gene_tag(b"HLbosTau10.ENSBTAT00000064442#ARPC5#34"),
            Some(b"arpc5".to_vec())
        );
        assert_eq!(
            parse_gene_tag(b"mm10.ENSMUST00000097536#Arpc5#16"),
            Some(b"arpc5".to_vec())
        );
        assert_eq!(parse_gene_tag(b"HAP1_SUPER_2-g175.t1"), None);
        assert_eq!(parse_gene_tag(b"TX##1"), None);
    }

    #[test]
    fn arpc5_truncated_cds_parents_are_masked() {
        // Main ARPC5 group vs truncated-CDS HLbosTau10.ENSBTAT group.
        let main = group(
            &[(48579271, 48579538), (48581899, 48581972), (48584202, 48584379)],
            &[(48579538, 48581899), (48581972, 48584202)],
            &[(48579395, 48579538), (48581899, 48581972)],
            &["arpc5"],
        );
        let trunc = group(
            &[(48579090, 48580309), (48581719, 48581972), (48584202, 48584379)],
            &[(48580309, 48581719), (48581972, 48584202)],
            &[(48579395, 48579894)],
            &["arpc5"],
        );
        let info = DashMap::new();
        info.insert(b"A".to_vec(), main);
        info.insert(b"B".to_vec(), trunc);
        let parents: HashSet<Vec<u8>> =
            [b"A".to_vec(), b"B".to_vec()].into_iter().collect();
        let mask = build_exclusion_mask(&parents, &info);
        assert_eq!(mask.len(), 1);
        assert!(is_masked(&mask, b"A", b"B"));
        assert!(is_masked(&mask, b"B", b"A"));
    }

    #[test]
    fn distant_genes_are_not_masked() {
        let g1 = group(
            &[(1000, 1200), (1500, 1700)],
            &[(1200, 1500)],
            &[(1000, 1700)],
            &["gene1"],
        );
        let g2 = group(
            &[(10000, 10200), (10500, 10700)],
            &[(10200, 10500)],
            &[(10000, 10700)],
            &["gene2"],
        );
        let info = DashMap::new();
        info.insert(b"A".to_vec(), g1);
        info.insert(b"B".to_vec(), g2);
        let parents: HashSet<Vec<u8>> =
            [b"A".to_vec(), b"B".to_vec()].into_iter().collect();
        let mask = build_exclusion_mask(&parents, &info);
        assert!(mask.is_empty());
    }

    #[test]
    fn liftover_shift_covered_by_cds_plus_splice() {
        // No exact exon (3bp shift), but CDS overlaps + 1 shared splice site.
        let a = group(
            &[(1000, 1100)],
            &[(1100, 1500)],
            &[(1000, 1100)],
            &[],
        );
        let b = group(
            &[(1000, 1103)],
            &[(1103, 1500)],
            &[(1000, 1103)],
            &[],
        );
        let info = DashMap::new();
        info.insert(b"A".to_vec(), a);
        info.insert(b"B".to_vec(), b);
        let parents: HashSet<Vec<u8>> =
            [b"A".to_vec(), b"B".to_vec()].into_iter().collect();
        let mask = build_exclusion_mask(&parents, &info);
        assert_eq!(mask.len(), 1);
    }
}
