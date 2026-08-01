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
use std::sync::Mutex;

use crate::report::{ClassifiedRead, FinalDecision, ReadReport};

/// Parallel accumulator for the processing function
///
/// Every classified read is routed by its own [`FinalDecision`], so the BED files and
/// the TSV report always describe the same decision.
///
/// # Fields
///
/// * `keep` - Keeps the lines of the input file that are kept
/// * `orphans` - Keeps the lines of the input file that are orphans
/// * `reports` - Keeps one classification report per processed query record
///
/// # Example
///
/// ```
/// use iso_orphan::utils::ParallelAccumulator;
///
/// let accumulator = ParallelAccumulator::default();
/// assert_eq!(accumulator.num_reports(), 0);
/// ```
pub struct ParallelAccumulator {
    pub keep: DashSet<Vec<u8>>,
    pub orphans: DashSet<Vec<u8>>,
    pub reports: Mutex<Vec<ReadReport>>,
}

/// Default implementation for the parallel accumulator
///
/// # Example
///
/// ```
/// use iso_orphan::utils::ParallelAccumulator;
///
/// let accumulator = ParallelAccumulator::default();
/// assert!(accumulator.keep.is_empty());
/// ```
impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            keep: DashSet::new(),
            orphans: DashSet::new(),
            reports: Mutex::new(Vec::new()),
        }
    }
}

impl ParallelAccumulator {
    /// Add classified reads to the accumulator
    ///
    /// The BED line of each read is inserted into the retention or the scrap set
    /// according to its decision, while its report is appended verbatim. BED sets
    /// deduplicate identical lines; reports never do, so every processed query record
    /// keeps exactly one row.
    ///
    /// # Parameters
    ///
    /// - `classified`: The classified reads of one component.
    ///
    /// # Example
    ///
    /// ```
    /// use iso_orphan::report::ClassifiedRead;
    /// use iso_orphan::utils::ParallelAccumulator;
    ///
    /// let accumulator = ParallelAccumulator::default();
    /// let classified: Vec<ClassifiedRead> = Vec::new();
    /// accumulator.add(classified);
    ///
    /// assert_eq!(accumulator.num_reports(), 0);
    /// ```
    pub fn add(&self, classified: Vec<ClassifiedRead>) {
        let mut reports = Vec::with_capacity(classified.len());

        for read in classified {
            let ClassifiedRead { bed_line, report } = read;

            match report.decision {
                FinalDecision::Pass => self.keep.insert(bed_line),
                FinalDecision::Scrap => self.orphans.insert(bed_line),
            };

            reports.push(report);
        }

        self.reports
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .extend(reports);
    }

    /// Get the number of orphans
    ///
    /// # Example
    ///
    /// ```
    /// use iso_orphan::utils::ParallelAccumulator;
    ///
    /// let accumulator = ParallelAccumulator::default();
    /// assert_eq!(accumulator.num_orphans(), 0);
    /// ```
    pub fn num_orphans(&self) -> usize {
        self.orphans.len()
    }

    /// Get the number of accumulated reports
    ///
    /// # Example
    ///
    /// ```
    /// use iso_orphan::utils::ParallelAccumulator;
    ///
    /// let accumulator = ParallelAccumulator::default();
    /// assert_eq!(accumulator.num_reports(), 0);
    /// ```
    pub fn num_reports(&self) -> usize {
        self.reports
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .len()
    }

    /// Take the accumulated reports out of the accumulator
    ///
    /// Rows are returned in accumulation order, which depends on thread scheduling;
    /// callers must sort them with `crate::report::sort_reports` before writing.
    ///
    /// # Example
    ///
    /// ```
    /// use iso_orphan::utils::ParallelAccumulator;
    ///
    /// let accumulator = ParallelAccumulator::default();
    /// let reports = accumulator.take_reports();
    ///
    /// assert!(reports.is_empty());
    /// ```
    pub fn take_reports(&self) -> Vec<ReadReport> {
        std::mem::take(
            &mut *self
                .reports
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner()),
        )
    }
}

/// Parallel counter for the processing function.
///
/// Tracks classification decisions across guided and de novo modes.
pub struct ParallelCounter {
    // general
    pub components: AtomicU32,

    // guided mode — reads retained on reference evidence
    pub guided_me_junction_keep: AtomicU32,
    pub guided_me_boundary_keep: AtomicU32,
    pub guided_se_keep: AtomicU32,

    // guided mode — reads that continue to de novo evaluation
    pub guided_oob_denovo: AtomicU32,
    pub guided_fallback_denovo: AtomicU32,

    // de novo mode — component context, reported but never terminal
    pub reads_in_single_read_component: AtomicU32,
    pub reads_in_low_support_component: AtomicU32,

    // de novo mode — outcomes
    pub splice_score_unavailable: AtomicU32,
    pub denovo_me_cluster_keep: AtomicU32,
    pub denovo_me_cluster_orphan: AtomicU32,
    pub denovo_me_intron_keep: AtomicU32,
    pub denovo_me_splice_keep: AtomicU32,
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
            guided_se_keep: AtomicU32::new(0),
            guided_oob_denovo: AtomicU32::new(0),
            guided_fallback_denovo: AtomicU32::new(0),
            reads_in_single_read_component: AtomicU32::new(0),
            reads_in_low_support_component: AtomicU32::new(0),
            splice_score_unavailable: AtomicU32::new(0),
            denovo_me_cluster_keep: AtomicU32::new(0),
            denovo_me_cluster_orphan: AtomicU32::new(0),
            denovo_me_intron_keep: AtomicU32::new(0),
            denovo_me_splice_keep: AtomicU32::new(0),
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

    /// Single-exon read kept by reciprocal overlap with reference exon.
    pub fn inc_guided_se_keep(&self) {
        self.guided_se_keep.fetch_add(1, Ordering::Relaxed);
    }

    // -- guided mode: reads continuing to de novo evaluation --

    /// Reads with no reference-exon overlap, evaluated de novo.
    pub fn inc_guided_oob_denovo(&self, count: u32) {
        self.guided_oob_denovo.fetch_add(count, Ordering::Relaxed);
    }

    /// Reads that overlap a reference exon but whose reference evidence was
    /// insufficient, continuing to de novo evaluation.
    pub fn inc_guided_fallback_denovo(&self) {
        self.guided_fallback_denovo.fetch_add(1, Ordering::Relaxed);
    }

    // -- de novo mode --

    /// Read sitting alone in its component. Context only, never terminal.
    pub fn inc_reads_in_single_read_component(&self) {
        self.reads_in_single_read_component
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Read in a component smaller than min_cluster_support. Context only.
    pub fn inc_reads_in_low_support_component(&self) {
        self.reads_in_low_support_component
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Read for which no splice-site evidence could be computed.
    ///
    /// Counted whether or not scores were supplied: an unstranded read, or one outside
    /// BigWig coverage, cannot be rescued or vetoed by splice evidence, and a large count
    /// here explains an otherwise surprising number of scraps.
    pub fn inc_splice_score_unavailable(&self) {
        self.splice_score_unavailable
            .fetch_add(1, Ordering::Relaxed);
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

    /// Multi-exon read kept by splice-site evidence alone (no structural support).
    pub fn inc_denovo_me_splice_keep(&self) {
        self.denovo_me_splice_keep.fetch_add(1, Ordering::Relaxed);
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
        let se_k = self.guided_se_keep.load(Ordering::Relaxed);

        report.push_str(&format!(
            "Reference-supported: multi-exon junction={}, multi-exon boundary={}, single-exon overlap={}\n",
            me_jk, me_bk, se_k
        ));
        let oob = self.guided_oob_denovo.load(Ordering::Relaxed);
        let fallback = self.guided_fallback_denovo.load(Ordering::Relaxed);
        report.push_str(&format!(
            "Guided -> de novo: no-reference-exon-overlap={}, insufficient-reference-evidence={}\n",
            oob, fallback
        ));

        // de novo stats
        let dn_sr = self.reads_in_single_read_component.load(Ordering::Relaxed);
        let dn_sc = self.reads_in_low_support_component.load(Ordering::Relaxed);
        let dn_mek = self.denovo_me_cluster_keep.load(Ordering::Relaxed);
        let dn_meo = self.denovo_me_cluster_orphan.load(Ordering::Relaxed);
        let dn_mik = self.denovo_me_intron_keep.load(Ordering::Relaxed);
        let dn_msk = self.denovo_me_splice_keep.load(Ordering::Relaxed);
        let dn_mso = self.denovo_me_splice_orphan.load(Ordering::Relaxed);
        let dn_sek = self.denovo_se_cluster_keep.load(Ordering::Relaxed);
        let dn_seo = self.denovo_se_cluster_orphan.load(Ordering::Relaxed);

        report.push_str(&format!(
            "De novo context: reads alone in a component={}, reads in a low-support component={}\n",
            dn_sr, dn_sc
        ));
        report.push_str(&format!(
            "De novo multi-exon: cluster-keep={}, intron-support-keep={}, splice-keep={}, unsupported-orphan={}, splice-orphan={}\n",
            dn_mek, dn_mik, dn_msk, dn_meo, dn_mso
        ));
        report.push_str(&format!(
            "De novo single-exon clusters: keep={}, orphan={}\n",
            dn_sek, dn_seo
        ));
        report.push_str(&format!(
            "Reads without splice-site evidence: {}\n",
            self.splice_score_unavailable.load(Ordering::Relaxed)
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
/// ```
/// use iso_orphan::utils::{__report_stats, ParallelCounter};
///
/// let counter = ParallelCounter::new();
/// __report_stats(&counter);
/// ```
pub fn __report_stats(counter: &ParallelCounter) {
    let stats = counter.report_stats();
    log::info!("INFO: Reporting stats:\n{}", stats);
}
