// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for detecting orphans in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the core functions for detecting orphans in a query set of reads
//! and processing the components in a parallel fashion.
//!
//! In short, each component of reads is subjected to any of the two modes: guided or
//! self-guided. Guided mode is the default mode and is used when the user provides a reference file
//! as input. Self-guided mode is used when the user does not provide a reference file as input.
//! Both, guided and self-guided, cover an extensive amount of curated oprhan cases under
//! the assumption that they do not represent a valid source of evidence for transcription.
//! The process is heavily parallellized to offer fast performance on large datasets.

use dashmap::DashSet;
use std::sync::atomic::{AtomicU32, Ordering};

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
    pub keep: DashSet<Vec<u8>>,
    pub orphans: DashSet<Vec<u8>>,
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
    pub fn add(&self, keep: Vec<Vec<u8>>, discard: Vec<Vec<u8>>) {
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
    pub read_se_mc_reference_se: AtomicU32, // multi-comp reference single-exon
    pub read_se_mc_reference_me: AtomicU32, // multi-comp reference multi-exon
    pub read_se_sc_reference_se: AtomicU32, // single-read comp reference single-exon
    pub read_se_sc_reference_me: AtomicU32, // single-read comp reference multi-exon

    // multi-exon categories
    pub read_me_mc_reference_se: AtomicU32, // multi-comp reference single-exon
    pub read_me_mc_reference_me: AtomicU32, // multi-comp reference multi-exon
    pub read_me_sc_reference_se: AtomicU32, // single-comp reference single-exon
    pub read_me_sc_reference_me: AtomicU32, // single-comp reference multi-exon

    // de-novo categories
    pub read_se_sc_no_reference: AtomicU32, // single-exon comp
    pub component_less_than_threshold: AtomicU32, // less than treshold
    pub read_de_novo_unsupported: AtomicU32, // unsupported exon

    // rescue categories
    pub rescue: AtomicU32,                   // rescue
    pub read_no_splice_match: AtomicU32,     // no splice match
    pub component_above_discards: AtomicU32, // component above discards
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
            read_se_mc_reference_se: AtomicU32::new(0),
            read_se_mc_reference_me: AtomicU32::new(0),
            read_se_sc_reference_se: AtomicU32::new(0),
            read_se_sc_reference_me: AtomicU32::new(0),
            read_me_mc_reference_se: AtomicU32::new(0),
            read_me_mc_reference_me: AtomicU32::new(0),
            read_me_sc_reference_se: AtomicU32::new(0),
            read_me_sc_reference_me: AtomicU32::new(0),
            read_se_sc_no_reference: AtomicU32::new(0),
            component_less_than_threshold: AtomicU32::new(0),
            read_de_novo_unsupported: AtomicU32::new(0),
            read_no_splice_match: AtomicU32::new(0),
            rescue: AtomicU32::new(0),
            component_above_discards: AtomicU32::new(0),
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

    /// Return the number of components
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// assert_eq!(counter.num_components(), 0);
    /// ```
    pub fn num_components(&self) -> u32 {
        self.components.load(Ordering::Relaxed)
    }

    /// Increment the number of reads with single-exon multi-component reference single-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_mc_reference_se();
    ///
    /// assert_eq!(counter.read_se_mc_reference_se(), 1);
    /// ```
    pub fn inc_read_se_mc_reference_se(&self) {
        self.read_se_mc_reference_se.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon multi-component reference multi-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_mc_reference_me();
    ///
    /// assert_eq!(counter.read_se_mc_reference_me(), 1);
    /// ```
    pub fn inc_read_se_mc_reference_me(&self) {
        self.read_se_mc_reference_me.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon single-component reference single-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_sc_reference_se();
    ///
    /// assert_eq!(counter.read_se_sc_reference_se(), 1);
    /// ```
    pub fn inc_read_se_sc_reference_se(&self) {
        self.read_se_sc_reference_se.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon single-component reference multi-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_sc_reference_me();
    ///
    /// assert_eq!(counter.read_se_sc_reference_me(), 1);
    /// ```
    pub fn inc_read_se_sc_reference_me(&self) {
        self.read_se_sc_reference_me.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with multi-exon multi-component reference single-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_me_mc_reference_se();
    ///
    /// assert_eq!(counter.read_me_mc_reference_se(), 1);
    /// ```
    pub fn inc_read_me_mc_reference_se(&self) {
        self.read_me_mc_reference_se.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with multi-exon multi-component reference multi-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_me_mc_reference_me();
    ///
    /// assert_eq!(counter.read_me_mc_reference_me(), 1);
    /// ```
    pub fn inc_read_me_mc_reference_me(&self) {
        self.read_me_mc_reference_me.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with multi-exon single-component reference single-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_me_sc_reference_se();
    ///
    /// assert_eq!(counter.read_me_sc_reference_se(), 1);
    /// ```
    pub fn inc_read_me_sc_reference_se(&self) {
        self.read_me_sc_reference_se.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with multi-exon single-component reference multi-exon
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_me_sc_reference_me();
    ///
    /// assert_eq!(counter.read_me_sc_reference_me(), 1);
    /// ```
    pub fn inc_read_me_sc_reference_me(&self) {
        self.read_me_sc_reference_me.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with single-exon single-component reference no reference
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_se_sc_no_reference();
    ///
    /// assert_eq!(counter.read_se_sc_no_reference(), 1);
    /// ```
    pub fn inc_read_se_sc_no_reference(&self) {
        self.read_se_sc_no_reference.fetch_add(1, Ordering::Relaxed);
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
    /// counter.inc_component_less_than_threshold();
    ///
    /// assert_eq!(counter.component_less_than_threshold(), 1);
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
    /// counter.inc_read_de_novo_unsupported();
    ///
    /// assert_eq!(counter.read_de_novo_unsupported(), 1);
    /// ```
    pub fn inc_read_de_novo_unsupported(&self) {
        self.read_de_novo_unsupported
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of rescues
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_rescue();
    ///
    /// assert_eq!(counter.rescue(), 1);
    /// ```
    pub fn inc_rescue(&self) {
        self.rescue.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of reads with no splice match
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_read_no_splice_match();
    ///
    /// assert_eq!(counter.read_no_splice_match(), 1);
    /// ```
    pub fn inc_read_no_splice_match(&self) {
        self.read_no_splice_match.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the number of components above discards
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_component_above_discards();
    ///
    /// assert_eq!(counter.component_above_discards(), 1);
    /// ```
    pub fn inc_component_above_discards(&self) {
        self.component_above_discards
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Summarize the counter
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    ///
    /// counter.inc_read_se_mc_reference_se();
    /// counter.inc_read_se_mc_reference_me();
    /// counter.inc_read_se_sc_reference_se();
    /// counter.inc_read_se_sc_reference_me();
    /// counter.inc_read_me_mc_reference_se();
    /// counter.inc_read_me_sc_reference_se();
    /// counter.inc_read_me_sc_reference_me();
    /// counter.inc_read_se_sc_no_reference();
    /// counter.inc_component_less_than_threshold();
    /// counter.inc_read_de_novo_unsupported();
    /// counter.inc_read_no_splice_match();
    /// counter.inc_rescue();
    /// counter.inc_component_above_discards();
    //
    /// let report = counter.report_stats();
    /// ```
    pub fn report_stats(&self) -> String {
        let mut report = String::new();

        report.push_str(&format!(
            "Total number of components: {}\n",
            self.num_components()
        ));

        let inc_read_se_mc_reference_se = self.read_se_mc_reference_se.load(Ordering::Relaxed);
        let inc_read_se_mc_reference_me = self.read_se_mc_reference_me.load(Ordering::Relaxed);
        let inc_read_se_sc_reference_se = self.read_se_sc_reference_se.load(Ordering::Relaxed);
        let inc_read_se_sc_reference_me = self.read_se_sc_reference_me.load(Ordering::Relaxed);
        let inc_read_me_mc_reference_se = self.read_me_mc_reference_se.load(Ordering::Relaxed);
        let inc_read_me_mc_reference_me = self.read_me_mc_reference_me.load(Ordering::Relaxed);
        let inc_read_me_sc_reference_se = self.read_me_sc_reference_se.load(Ordering::Relaxed);
        let inc_read_me_sc_reference_me = self.read_me_sc_reference_me.load(Ordering::Relaxed);
        let inc_read_se_sc_no_reference = self.read_se_sc_no_reference.load(Ordering::Relaxed);
        let inc_component_less_than_threshold =
            self.component_less_than_threshold.load(Ordering::Relaxed);
        let inc_read_de_novo_unsupported = self.read_de_novo_unsupported.load(Ordering::Relaxed);
        let inc_read_no_splice_match = self.read_no_splice_match.load(Ordering::Relaxed);
        let inc_rescue = self.rescue.load(Ordering::Relaxed);
        let inc_component_above_discards = self.component_above_discards.load(Ordering::Relaxed);

        report.push_str(&format!(
            "Total number of reads with single-exon multi-component reference single-exon: {}\n",
            inc_read_se_mc_reference_se
        ));
        report.push_str(&format!(
            "Total number of reads with single-exon multi-component reference multi-exon: {}\n",
            inc_read_se_mc_reference_me
        ));
        report.push_str(&format!(
            "Total number of reads with single-exon single-component reference single-exon: {}\n",
            inc_read_se_sc_reference_se
        ));
        report.push_str(&format!(
            "Total number of reads with single-exon single-component reference multi-exon: {}\n",
            inc_read_se_sc_reference_me
        ));
        report.push_str(&format!(
            "Total number of reads with multi-exon multi-component reference single-exon: {}\n",
            inc_read_me_mc_reference_se
        ));
        report.push_str(&format!(
            "Total number of reads with multi-exon multi-component reference multi-exon: {}\n",
            inc_read_me_mc_reference_me
        ));
        report.push_str(&format!(
            "Total number of reads with multi-exon single-component reference single-exon: {}\n",
            inc_read_me_sc_reference_se
        ));
        report.push_str(&format!(
            "Total number of reads with multi-exon single-component reference multi-exon: {}\n",
            inc_read_me_sc_reference_me
        ));
        report.push_str(&format!(
            "Total number of reads with single-exon single-component reference no reference: {}\n",
            inc_read_se_sc_no_reference
        ));
        report.push_str(&format!(
            "Total number of denovo reads with less than min_read_num_denovo: {}\n",
            inc_component_less_than_threshold
        ));
        report.push_str(&format!(
            "Total number of reads with unsupported exon: {}\n",
            inc_read_de_novo_unsupported
        ));
        report.push_str(&format!(
            "Total number of components above min_discard_percentage: {}\n",
            inc_component_above_discards
        ));
        report.push_str(&format!(
            "Total number of reads with no splice match against reference: {}\n",
            inc_read_no_splice_match
        ));
        report.push_str(&format!("Total number of rescues: {}\n", inc_rescue));

        report
    }
}

/// Report the stats of the counter to the log
///
/// # Arguments
///
/// * `counter` - The counter to report
///
/// # Example
///
/// ```
/// use isotools::utils::ParallelCounter;
///
/// let counter = ParallelCounter::new();
/// counter.inc_read_se_mc_reference_se(5);
/// counter.inc_read_se_mc_reference_me(5);
/// counter.inc_read_se_sc_reference_se(5);
/// counter.inc_read_se_sc_reference_me(5);
/// counter.inc_read_me_mc_reference_se(5);
/// counter.inc_read_me_sc_reference_se(5);
/// counter.inc_read_me_sc_reference_me(5);
/// counter.inc_read_se_sc_no_reference(5);
/// counter.inc_component_less_than_threshold(5);
/// counter.inc_read_de_novo_unsupported(5);
/// counter.inc_read_no_splice_match(5);
/// counter.inc_rescue(5);
/// counter.inc_component_above_discards(5);
///
/// __report_stats(&counter);
/// ```
pub fn __report_stats(counter: &ParallelCounter) {
    let stats = counter.report_stats();
    log::info!("INFO: Reporting stats:\n{}", stats);
}
