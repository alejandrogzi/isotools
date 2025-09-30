//! Core module for detecting orphans in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the core functions for detecting orphans in a query set of reads
//! and processing the components in a parallel fashion.
//!
//! In short, each component of reads is subjected to any of the two modes: guided or
//! self-guided. Guided mode is the default mode and is used when the user provides a TOGA file
//! as input. Self-guided mode is used when the user does not provide a TOGA file as input.
//! Both, guided and self-guided, cover an extensive amount of curated oprhan cases under
//! the assumption that they do not represent a valid source of evidence for transcription.
//! The process is heavily parallellized to offer fast performance on large datasets.

use std::sync::atomic::{AtomicU32, Ordering};

use config::ParallelCollector;
use dashmap::DashSet;

/// Parallel accumulator for the processing function
///
/// # Fields
///
/// * `keep` - Keeps the lines of the input file that are kept
/// * `orphans` - Keeps the lines of the input file that are orphans
///
/// # Example
///
/// ```
/// use isotools::utils::ParallelAccumulator;
///
/// let accumulator = ParallelAccumulator::default();
/// ```
pub struct ParallelAccumulator {
    pub keep: DashSet<String>,
    pub orphans: DashSet<String>,
}

/// Default implementation for the parallel accumulator
///
/// # Example
///
/// ```
/// use isotools::utils::ParallelAccumulator;
///
/// let accumulator = ParallelAccumulator::default();
/// ```
impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            keep: DashSet::new(),
            orphans: DashSet::new(),
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

        collections.push(&self.keep);
        collections.push(&self.orphans);

        Ok(collections)
    }
}

impl ParallelAccumulator {
    /// Number of fields in the accumulator of type DashSet<String>
    pub const NUM_FIELDS: usize = 2;

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
    pub fn add(&self, keep: Vec<String>, discard: Vec<String>) {
        for item in keep {
            self.keep.insert(item);
        }
        for item in discard {
            self.orphans.insert(item);
        }
    }

    /// Get the number of orphans
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let accumulator = ParallelAccumulator::default();
    /// assert_eq!(accumulator.oprhans(), 0);
    /// ```
    pub fn num_orphans(&self) -> usize {
        self.orphans.len()
    }
}

/// Parallel counter for the processing function
pub struct ParallelCounter {
    // general
    pub components: AtomicU32,

    // single-exon categories
    pub read_se_mc_toga_se: AtomicU32, // multi-comp TOGA single-exon
    pub read_se_mc_toga_me: AtomicU32, // multi-comp TOGA multi-exon
    pub read_se_sc_toga_se: AtomicU32, // single-read comp TOGA single-exon
    pub read_se_sc_toga_me: AtomicU32, // single-read comp TOGA multi-exon

    // multi-exon categories
    pub read_me_mc_toga_se: AtomicU32, // multi-comp TOGA single-exon
    pub read_me_sc_toga_se: AtomicU32, // single-comp TOGA single-exon
    pub read_me_sc_toga_me: AtomicU32, // single-comp TOGA multi-exon

    // de-novo categories
    pub read_se_sc_no_toga: AtomicU32, // single-exon comp
    pub component_less_than_threshold: AtomicU32, // less than treshold
    pub read_de_novo_unsupported: AtomicU32, // unsupported exon
}

/// Default implementation for the parallel counter
///
/// # Example
///
/// ```
/// use isotools::utils::ParallelCounter;
///
/// let counter = ParallelCounter::default();
/// ```
impl Default for ParallelCounter {
    fn default() -> Self {
        Self {
            components: AtomicU32::new(0),
            read_se_mc_toga_se: AtomicU32::new(0),
            read_se_mc_toga_me: AtomicU32::new(0),
            read_se_sc_toga_se: AtomicU32::new(0),
            read_se_sc_toga_me: AtomicU32::new(0),
            read_me_mc_toga_se: AtomicU32::new(0),
            read_me_sc_toga_se: AtomicU32::new(0),
            read_me_sc_toga_me: AtomicU32::new(0),
            read_se_sc_no_toga: AtomicU32::new(0),
            component_less_than_threshold: AtomicU32::new(0),
            read_de_novo_unsupported: AtomicU32::new(0),
        }
    }
}

impl ParallelCounter {
    /// Create a new parallel counter
    ///
    /// # Example
    ///
    /// ```
    /// use isotools::utils::ParallelCounter;
    ///
    /// let counter = ParallelCounter::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
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
    /// assert_eq!(counter.components(), 5);
    /// ```
    pub fn inc_components(&self, count: u32) {
        self.components.fetch_add(count, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon multi-component TOGA single-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_mc_toga_se(5);
    ///
    /// assert_eq!(counter.read_se_mc_toga_se(), 5);
    /// ```
    pub fn inc_read_se_mc_toga_se(&self) {
        self.read_se_mc_toga_se.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon multi-component TOGA multi-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_mc_toga_me(5);
    ///
    /// assert_eq!(counter.read_se_mc_toga_me(), 5);
    /// ```
    pub fn inc_read_se_mc_toga_me(&self) {
        self.read_se_mc_toga_me.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon single-component TOGA single-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_sc_toga_se(5);
    ///
    /// assert_eq!(counter.read_se_sc_toga_se(), 5);
    /// ```
    pub fn inc_read_se_sc_toga_se(&self) {
        self.read_se_sc_toga_se.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon single-component TOGA multi-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_sc_toga_me(5);
    ///
    /// assert_eq!(counter.read_se_sc_toga_me(), 5);
    /// ```
    pub fn inc_read_se_sc_toga_me(&self) {
        self.read_se_sc_toga_me.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with multi-exon multi-component TOGA single-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_me_mc_toga_se(5);
    ///
    /// assert_eq!(counter.read_me_mc_toga_se(), 5);
    /// ```
    pub fn inc_read_me_mc_toga_se(&self) {
        self.read_me_mc_toga_se.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with multi-exon single-component TOGA single-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_me_sc_toga_se(5);
    ///
    /// assert_eq!(counter.read_me_sc_toga_se(), 5);
    /// ```
    pub fn inc_read_me_sc_toga_se(&self) {
        self.read_me_sc_toga_se.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with multi-exon single-component TOGA multi-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_me_sc_toga_me(5);
    ///
    /// assert_eq!(counter.read_me_sc_toga_me(), 5);
    /// ```
    pub fn inc_read_me_sc_toga_me(&self) {
        self.read_me_sc_toga_me.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon single-component TOGA no TOGA
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_sc_no_toga(5);
    ///
    /// assert_eq!(counter.read_se_sc_no_toga(), 5);
    /// ```
    pub fn inc_read_se_sc_no_toga(&self) {
        self.read_se_sc_no_toga.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with less than threshold
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_component_less_than_threshold(5);
    ///
    /// assert_eq!(counter.component_less_than_threshold(), 5);
    /// ```
    pub fn inc_component_less_than_threshold(&self) {
        self.component_less_than_threshold
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with unsupported exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_de_novo_unsupported(5);
    ///
    /// assert_eq!(counter.read_de_novo_unsupported(), 5);
    /// ```
    pub fn inc_read_de_novo_unsupported(&self) {
        self.read_de_novo_unsupported
            .fetch_add(1, Ordering::Relaxed);
    }
}
