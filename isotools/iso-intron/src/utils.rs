//! Core module for detecting intron retentions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting intron retentions
//! and processing the components of reads and introns in parallel.
//!
//! In short, each read is checked for the presence of intron retentions
//! or RT introns. If a read has a true intron retention, it is discarded.
//! If a read has an RT intron, it is also discarded. The veracity of an
//! 'intron' is determined by 'iso-classify', using machine-learning models,
//! ab initio gene prediction, and other heuristics. The process is heavily
//! parallelized to offer fast performance on large datasets.

use dashmap::{DashMap, DashSet};
use hashbrown::{HashMap, HashSet};
use packbed::{par_reader, record::Bed4};

use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Arc;

use config::{bed_to_map, CoordType, ModuleMap, ParallelCollector};

/// Parallel counter for the processing function
///
/// # Fields
///
/// - `dirties`: Atomic counter for dirty items.
/// - `components`: Atomic counter for components.
/// - `retentions`: Atomic counter for retentions.
///
/// # Example
///
/// ```rust, no_run
/// let counter = ParallelCounter::default();
///
/// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 0);
/// ```
pub struct ParallelCounter {
    pub dirties: AtomicU32,
    pub components: AtomicU32,
    pub retentions: AtomicU32,
}

impl ParallelCounter {
    /// Constructor for ParallelCounter
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    ///
    /// assert_eq!(counter.dirties.load(Ordering::Relaxed), 0);
    /// assert_eq!(counter.components.load(Ordering::Relaxed), 0);
    /// assert_eq!(counter.retentions.load(Ordering::Relaxed), 0);
    /// ```
    fn new() -> Self {
        Self {
            dirties: AtomicU32::new(0),
            components: AtomicU32::new(0),
            retentions: AtomicU32::new(0),
        }
    }

    /// Increment the number of components
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_components(5);
    ///
    /// assert_eq!(counter.components.load(Ordering::Relaxed), 5);
    /// ```
    pub fn inc_components(&self, count: u32) {
        self.components.fetch_add(count, Ordering::Relaxed);
    }

    /// Increment the number of dirty items
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_dirty();
    ///
    /// assert_eq!(counter.dirties.load(Ordering::Relaxed), 1);
    /// ```
    pub fn inc_dirty(&self) {
        self.dirties.fetch_add(1, Ordering::Relaxed);
    }

    /// Get the number of dirty items
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_dirty();
    ///
    /// assert_eq!(counter.get_dirty(), 1);
    /// ```
    pub fn get_counters(&self) -> (f64, f64) {
        (
            self.dirties.load(Ordering::Relaxed) as f64,
            self.components.load(Ordering::Relaxed) as f64,
        )
    }

    /// Get the number of dirty items
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_dirty();
    ///
    /// assert_eq!(counter.get_dirty(), 1);
    /// ```
    pub fn get_stat(&self) -> (f64, f64) {
        let (dirties, components) = self.get_counters();
        (dirties, (dirties / components) * 100.0)
    }

    /// Get the number of retentions
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_retentions();
    ///
    /// assert_eq!(counter.get_retentions(), 1);
    /// ```
    pub fn inc_retentions(&self) {
        self.retentions.fetch_add(1, Ordering::Relaxed);
    }
}

/// Default implementation for ParallelCounter
///
/// # Example
///
/// ```rust, no_run
/// let counter = ParallelCounter::default();
///
/// assert_eq!(counter.dirties.load(Ordering::Relaxed), 0);
/// assert_eq!(counter.components.load(Ordering::Relaxed), 0);
/// assert_eq!(counter.retentions.load(Ordering::Relaxed), 0);
/// ```
impl Default for ParallelCounter {
    fn default() -> Self {
        Self::new()
    }
}

/// Parallel accumulator for the processing function
///
/// # Fields
///
/// - `retentions`: A set of strings representing the retentions.
/// - `non_retentions`: A set of strings representing the non-retentions.
/// - `miscellaneous`: A set of strings representing miscellaneous items.
/// - `descriptor`: A map of strings to boxed `ModuleMap` trait objects.
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.pass.len(), 0);
/// ```
pub struct ParallelAccumulator {
    pub retentions: DashSet<String>,
    pub non_retentions: DashSet<String>,
    pub miscellaneous: DashSet<String>,
    pub descriptor: DashMap<String, Box<dyn ModuleMap>>,
}

/// ParallelAccumulator constructor
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.retentions.len(), 0);
/// assert_eq!(accumulator.non_retentions.len(), 0);
/// assert_eq!(accumulator.miscellaneous.len(), 0);
/// assert_eq!(accumulator.descriptor.len(), 0);
/// ```
impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            retentions: DashSet::new(),
            non_retentions: DashSet::new(),
            miscellaneous: DashSet::new(),
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

        collections.push(&self.retentions);
        collections.push(&self.non_retentions);
        collections.push(&self.miscellaneous);

        Ok(collections)
    }
}

impl ParallelAccumulator {
    /// Number of fields in the accumulator of type DashSet<String>
    pub const NUM_FIELDS: usize = 3;

    /// Get the number of retentions
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let accumulator = ParallelAccumulator::default();
    /// assert_eq!(accumulator.num_retentions(), 0);
    /// ```
    pub fn num_retentions(&self) -> usize {
        self.retentions.len()
    }

    /// Add items to the accumulator
    ///
    /// # Parameters
    ///
    /// - `keep`: A vector of strings to be retained.
    /// - `discard`: A vector of strings to be discarded.
    /// - `descriptor`: A HashMap of strings to boxed `ModuleMap` trait objects.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut accumulator = ParallelAccumulator::default();
    /// accumulator.add(vec!["item1".to_string()], vec!["item2".to_string()], HashMap::new());
    ///
    /// assert_eq!(accumulator.num_retentions(), 1);
    /// ```
    pub fn add(
        &self,
        keep: Vec<String>,
        discard: Vec<String>,
        review: Option<Vec<String>>,
        descriptor: HashMap<String, Box<dyn ModuleMap>>,
    ) {
        for item in keep {
            self.non_retentions.insert(item);
        }
        for item in discard {
            self.retentions.insert(item);
        }
        for (key, value) in descriptor {
            self.descriptor.insert(key, value);
        }

        if let Some(review) = review {
            for item in review {
                self.miscellaneous.insert(item);
            }
        }
    }
}

/// Unpack blacklist from a vector of .bed paths
///
/// # Parameters
///
/// - `paths`: A vector of PathBuf representing the paths to the .bed files.
///
/// # Returns
///
/// - An Option containing a HashMap of strings to HashSets of tuples (u64, u64) if the paths are not empty.
///
/// # Example
///
/// ```rust, no_run
/// let paths = vec![PathBuf::from("path/to/bed1.bed"), PathBuf::from("path/to/bed2.bed")];
/// let result = unpack_blacklist(paths);
///
/// assert!(result.is_some());
/// ```
pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    if paths.is_empty() {
        return None;
    }

    let contents = Arc::new(par_reader(paths).unwrap());
    let tracks = bed_to_map::<Bed4>(contents, CoordType::Bounds).unwrap();

    Some(tracks)
}
