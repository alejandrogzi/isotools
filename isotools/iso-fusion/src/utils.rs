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

use std::fmt::Debug;
use std::path::Path;
use std::str::FromStr;
use std::sync::atomic::{AtomicU32, Ordering};

use anyhow::{Ok, Result};
use dashmap::{DashMap, DashSet};
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{packbed, par_reader, Bed12, GenePred, RefGenePred};
use rayon::prelude::*;

use config::{get_progress_bar, ModuleMap, OverlapType, ParallelCollector, TsvParser};

/// Unpack blacklist from a vector of .bed paths
///
/// # Parameters
///
/// - `paths`: A vector of PathBuf representing the paths to the .bed files.
///
/// # Returns
///
/// - An Option containing a HashMap of strings to HashSets of Strings if the paths are not empty.
///
/// # Example
///
/// ```rust, no_run
/// let paths = vec![PathBuf::from("path/to/bed1.bed"), PathBuf::from("path/to/bed2.bed")];
/// let result = unpack_blacklist(paths);
///
/// assert!(result.is_some());
/// ```
pub fn unpack_blacklist<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
    cds_overlap: OverlapType,
    is_ref: bool,
) -> Result<HashMap<String, HashSet<String>>, anyhow::Error> {
    let contents = par_reader(files)?;
    let tracks = parse_tracks(&contents, cds_overlap, is_ref)?;

    Ok(tracks)
}

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
pub fn prepare_refs<P: AsRef<Path> + Debug + Sync + Send>(
    refs: Vec<P>,
) -> Result<String, anyhow::Error> {
    let tracks = packbed(
        refs,
        None,
        config::OverlapType::Exon,
        packbed::PackMode::Default,
    )?;

    let results: Vec<String> = tracks
        .into_par_iter()
        .map(|bucket| {
            let components = bucket.1;
            let mut acc = String::new();

            for comp in components {
                let refs = comp
                    .as_any()
                    .downcast_ref::<(RefGenePred, Vec<GenePred>)>()
                    .expect("ERROR: Failed to downcast to GenePred")
                    .0
                    .clone();

                let loci = refs.merge_names();
                refs.reads.into_iter().for_each(|mut record| {
                    let mut line = record.line.clone();

                    // if loci has more than 1 gene name -> fusion
                    if loci.split('.').count() > 1 {
                        line = record.mut_name_from_line(&loci);
                    }

                    acc += format!("{}\n", line).as_str();
                });
            }

            acc
        })
        .collect();

    info!("Reference transcripts fixed for fusions!");
    Ok(results.concat())
}

/// Parse tracks from a string [iso-fusion specific]
///
/// # Parameters
///
/// - `contents`: A string slice containing the contents to parse.
/// - `cds_overlap`: An `OverlapType` enum value.
/// - `is_ref`: A boolean indicating if the contents are reference.
///
/// # Returns
///
/// - A Result containing a HashMap of strings to HashSets of strings.
///
/// # Example
///
/// ```rust, no_run
/// let contents = "chr1\t100\t200\tgene1\nchr2\t150\t250\tgene2";
/// let cds_overlap = OverlapType::Exon;
/// let is_ref = true;
///
/// let result = parse_tracks(contents, cds_overlap, is_ref);
///
/// assert!(result.is_ok());
/// ```
fn parse_tracks<'a>(
    contents: &'a str,
    cds_overlap: OverlapType,
    is_ref: bool,
) -> Result<HashMap<String, HashSet<String>>, anyhow::Error> {
    let pb = get_progress_bar(
        contents.lines().count() as u64,
        "Parsing BED12 blacklist...",
    );
    let tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|row| Bed12::read(row, cds_overlap, is_ref).ok())
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, HashSet<String>>, record| {
                acc.entry(record.chrom.clone())
                    .or_default()
                    .insert(record.name);
                pb.inc(1);
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (k, v) in map {
                    let acc_v = acc.entry(k).or_insert(HashSet::new());
                    acc_v.extend(v);
                }
                acc
            },
        );

    pb.finish_and_clear();
    info!("Records parsed: {}", tracks.values().flatten().count());

    Ok(tracks)
}

/// Isoform parser to parse data usin tsv_to_map function
#[derive(Debug)]
pub struct IsoformParser {
    fields: Vec<String>,
}

/// TsvParser trait for IsoformParser
impl TsvParser for IsoformParser {
    fn parse(line: &str) -> Result<Self, anyhow::Error> {
        Ok(Self {
            fields: line.split('\t').map(|s| s.to_string()).collect(),
        })
    }

    fn key(&self, index: usize) -> &str {
        &self.fields[index]
    }

    fn value<V: FromStr>(&self, index: usize) -> Result<V, V::Err> {
        self.fields[index].parse::<V>()
    }
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
    pub fakes: DashSet<String>,
    pub descriptor: DashMap<String, Box<dyn ModuleMap>>,
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
            fakes: DashSet::new(),
            descriptor: DashMap::new(),
        }
    }
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
        collections.push(&self.fakes);

        std::result::Result::Ok(collections)
    }
}

impl ParallelAccumulator {
    /// Number of fields in the accumulator of type DashSet<String>
    pub const NUM_FIELDS: usize = 4;

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
        fakes: Vec<String>,
        descriptor: HashMap<String, Box<dyn ModuleMap>>,
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
        };

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
