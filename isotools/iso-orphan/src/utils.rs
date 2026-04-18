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

/// Parallel counter for the processing function.
///
/// Tracks classification decisions across guided and de novo modes.
pub struct ParallelCounter {
    // general
    pub components: AtomicU32,

    // guided mode — multi-exon reads
    pub guided_me_junction_keep: AtomicU32,
    pub guided_me_boundary_keep: AtomicU32,
    pub guided_me_splice_keep: AtomicU32,
    pub guided_me_orphan: AtomicU32,

    // guided mode — single-exon reads
    pub guided_se_keep: AtomicU32,
    pub guided_se_orphan: AtomicU32,

    // guided mode — out-of-bounds redirect
    pub guided_oob_denovo: AtomicU32,

    // de novo mode
    pub denovo_single_read: AtomicU32,
    pub denovo_small_component: AtomicU32,
    pub denovo_me_cluster_keep: AtomicU32,
    pub denovo_me_cluster_orphan: AtomicU32,
    pub denovo_me_intron_keep: AtomicU32,
    pub denovo_me_splice_orphan: AtomicU32,
    pub denovo_se_cluster_keep: AtomicU32,
    pub denovo_se_cluster_orphan: AtomicU32,
}

impl Default for ParallelCounter {
    fn default() -> Self {
        Self {
            components: AtomicU32::new(0),
            guided_me_junction_keep: AtomicU32::new(0),
            guided_me_boundary_keep: AtomicU32::new(0),
            guided_me_splice_keep: AtomicU32::new(0),
            guided_me_orphan: AtomicU32::new(0),
            guided_se_keep: AtomicU32::new(0),
            guided_se_orphan: AtomicU32::new(0),
            guided_oob_denovo: AtomicU32::new(0),
            denovo_single_read: AtomicU32::new(0),
            denovo_small_component: AtomicU32::new(0),
            denovo_me_cluster_keep: AtomicU32::new(0),
            denovo_me_cluster_orphan: AtomicU32::new(0),
            denovo_me_intron_keep: AtomicU32::new(0),
            denovo_me_splice_orphan: AtomicU32::new(0),
            denovo_se_cluster_keep: AtomicU32::new(0),
            denovo_se_cluster_orphan: AtomicU32::new(0),
        }
    }
}

impl ParallelCounter {
    pub fn new() -> Self {
        Self::default()
    }

    // -- general --

    pub fn inc_components(&self, count: u32) {
        self.components.fetch_add(count, Ordering::Relaxed);
    }

    pub fn num_components(&self) -> u32 {
        self.components.load(Ordering::Relaxed)
    }

    // -- guided mode: multi-exon reads --

    /// Multi-exon read kept by junction scoring against individual reference transcripts.
    pub fn inc_guided_me_junction_keep(&self) {
        self.guided_me_junction_keep.fetch_add(1, Ordering::Relaxed);
    }

    /// Multi-exon read rescued by exon boundary matching (fallback when junction score is low).
    pub fn inc_guided_me_boundary_keep(&self) {
        self.guided_me_boundary_keep.fetch_add(1, Ordering::Relaxed);
    }

    /// Multi-exon read rescued by splice site score (fallback when junction score is low).
    pub fn inc_guided_me_splice_keep(&self) {
        self.guided_me_splice_keep.fetch_add(1, Ordering::Relaxed);
    }

    /// Multi-exon read orphaned (no junction support, no boundary match).
    pub fn inc_guided_me_orphan(&self) {
        self.guided_me_orphan.fetch_add(1, Ordering::Relaxed);
    }

    // -- guided mode: single-exon reads --

    /// Single-exon read kept by reciprocal overlap with reference exon.
    pub fn inc_guided_se_keep(&self) {
        self.guided_se_keep.fetch_add(1, Ordering::Relaxed);
    }

    /// Single-exon read orphaned (insufficient overlap with any reference exon).
    pub fn inc_guided_se_orphan(&self) {
        self.guided_se_orphan.fetch_add(1, Ordering::Relaxed);
    }

    // -- guided mode: out-of-bounds redirect --

    /// Out-of-bounds reads redirected from guided to de novo.
    pub fn inc_guided_oob_denovo(&self, count: u32) {
        self.guided_oob_denovo.fetch_add(count, Ordering::Relaxed);
    }

    // -- de novo mode --

    /// Single-read component orphaned (insufficient evidence without support).
    pub fn inc_denovo_single_read(&self) {
        self.denovo_single_read.fetch_add(1, Ordering::Relaxed);
    }

    /// Component smaller than min_cluster_support → all reads orphaned.
    pub fn inc_denovo_small_component(&self) {
        self.denovo_small_component.fetch_add(1, Ordering::Relaxed);
    }

    /// Multi-exon read kept (in a supported intron-chain cluster).
    pub fn inc_denovo_me_cluster_keep(&self) {
        self.denovo_me_cluster_keep.fetch_add(1, Ordering::Relaxed);
    }

    /// Multi-exon read orphaned (in an unsupported intron-chain cluster, low intron support).
    pub fn inc_denovo_me_cluster_orphan(&self) {
        self.denovo_me_cluster_orphan
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Multi-exon read kept by per-intron support (non-dominant cluster).
    pub fn inc_denovo_me_intron_keep(&self) {
        self.denovo_me_intron_keep.fetch_add(1, Ordering::Relaxed);
    }

    /// Multi-exon read orphaned by low splice-site score (would have passed intron support).
    pub fn inc_denovo_me_splice_orphan(&self) {
        self.denovo_me_splice_orphan.fetch_add(1, Ordering::Relaxed);
    }

    /// Single-exon read kept (in a supported overlap cluster).
    pub fn inc_denovo_se_cluster_keep(&self) {
        self.denovo_se_cluster_keep.fetch_add(1, Ordering::Relaxed);
    }

    /// Single-exon read orphaned (in an unsupported overlap cluster).
    pub fn inc_denovo_se_cluster_orphan(&self) {
        self.denovo_se_cluster_orphan
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Generate a formatted stats report.
    pub fn report_stats(&self) -> String {
        let mut report = String::new();

        report.push_str(&format!(
            "Total number of components: {}\n",
            self.num_components()
        ));

        // guided mode stats
        let me_jk = self.guided_me_junction_keep.load(Ordering::Relaxed);
        let me_bk = self.guided_me_boundary_keep.load(Ordering::Relaxed);
        let me_sk = self.guided_me_splice_keep.load(Ordering::Relaxed);
        let me_or = self.guided_me_orphan.load(Ordering::Relaxed);
        let se_k = self.guided_se_keep.load(Ordering::Relaxed);
        let se_o = self.guided_se_orphan.load(Ordering::Relaxed);

        report.push_str(&format!(
            "Guided multi-exon: junction-keep={}, boundary-keep={}, splice-keep={}, orphan={}\n",
            me_jk, me_bk, me_sk, me_or
        ));
        report.push_str(&format!(
            "Guided single-exon: keep={}, orphan={}\n",
            se_k, se_o
        ));
        let oob = self.guided_oob_denovo.load(Ordering::Relaxed);
        report.push_str(&format!("Guided out-of-bounds -> de novo: {}\n", oob));

        // de novo stats
        let dn_sr = self.denovo_single_read.load(Ordering::Relaxed);
        let dn_sc = self.denovo_small_component.load(Ordering::Relaxed);
        let dn_mek = self.denovo_me_cluster_keep.load(Ordering::Relaxed);
        let dn_meo = self.denovo_me_cluster_orphan.load(Ordering::Relaxed);
        let dn_mik = self.denovo_me_intron_keep.load(Ordering::Relaxed);
        let dn_mso = self.denovo_me_splice_orphan.load(Ordering::Relaxed);
        let dn_sek = self.denovo_se_cluster_keep.load(Ordering::Relaxed);
        let dn_seo = self.denovo_se_cluster_orphan.load(Ordering::Relaxed);

        report.push_str(&format!(
            "De novo: single-read-orphan={}, small-component-orphan={}\n",
            dn_sr, dn_sc
        ));
        report.push_str(&format!(
            "De novo multi-exon: cluster-keep={}, intron-support-keep={}, cluster-orphan={}, splice-orphan={}\n",
            dn_mek, dn_mik, dn_meo, dn_mso
        ));
        report.push_str(&format!(
            "De novo single-exon clusters: keep={}, orphan={}\n",
            dn_sek, dn_seo
        ));

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
/// ```rust, no_run
/// use isotools::utils::ParallelCounter;
///
/// let counter = ParallelCounter::new();
/// __report_stats(&counter);
/// ```
pub fn __report_stats(counter: &ParallelCounter) {
    let stats = counter.report_stats();
    log::info!("INFO: Reporting stats:\n{}", stats);
}
