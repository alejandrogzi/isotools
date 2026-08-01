// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Per-read classification report for iso-orphan.
//!
//! This module owns the structured classification result of the program. Every valid
//! query BED12 record processed by iso-orphan produces exactly one [`ClassifiedRead`],
//! which pairs the record's BED line with the [`ReadReport`] describing how it was
//! classified. The report is the single source of truth for:
//!
//! 1. the final [`FinalDecision`] (`PASS` or `SCRAP`),
//! 2. the [`DecisionReason`] codes explaining the decision path,
//! 3. the feature values calculated while classifying the read,
//! 4. the selection of the record for the HQ or the scraps BED file, and
//! 5. the rows of `<outdir>/orphans/<prefix>.report.tsv`.
//!
//! Because the BED selection and the TSV row are derived from the same value, the two
//! outputs cannot disagree.
//!
//! # Formatting contract
//!
//! * `.` marks an unavailable or non-applicable scalar value.
//! * Booleans are lowercase `true` / `false`.
//! * Floating-point features and thresholds use three decimal places.
//! * Reason codes are joined with `/` in evaluation order.
//! * Applied thresholds are joined with `;`.
//! * Textual BED fields are escaped so a single record cannot corrupt the TSV layout.
//!
//! # Example
//!
//! ```
//! use iso_orphan::report::{report_header, write_reports, REPORT_COLUMNS};
//!
//! let mut out = Vec::new();
//! write_reports(&mut out, &[]).unwrap();
//!
//! let text = String::from_utf8(out).unwrap();
//! assert_eq!(text, format!("{}\n", report_header()));
//! assert_eq!(REPORT_COLUMNS[REPORT_COLUMNS.len() - 2], "reasons");
//! assert_eq!(REPORT_COLUMNS[REPORT_COLUMNS.len() - 1], "decision");
//! ```

use std::io::{self, Write};

use genepred::{Bed12, GenePred, Strand};

use crate::scoring::BoundaryEvidence;

/// Placeholder written for unavailable or non-applicable values.
pub const MISSING: &str = ".";

/// Column order of `<prefix>.report.tsv`.
///
/// The last two columns are always `reasons` and `decision`.
pub const REPORT_COLUMNS: [&str; 26] = [
    "read_id",
    "chrom",
    "start",
    "end",
    "strand",
    "exon_count",
    "read_type",
    "mode",
    "evaluation_path",
    "group_size",
    "reference_count",
    "cluster_size",
    "overlaps_reference_exon",
    "best_reference_id",
    "best_reference_overlap",
    "junction_matches",
    "junction_count",
    "junction_fraction",
    "boundary_match",
    "boundary_evidence",
    "intron_support_fraction",
    "median_splice_score",
    "weakest_splice_score",
    "applied_thresholds",
    "reasons",
    "decision",
];

/// Final retention decision of a query read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum FinalDecision {
    /// The read is written to `<prefix>.hq.bed`.
    Pass,
    /// The read is written to `<prefix>.scraps.bed`.
    Scrap,
}

impl FinalDecision {
    /// Stable textual representation used in the `decision` column.
    pub fn as_str(&self) -> &'static str {
        match self {
            FinalDecision::Pass => "PASS",
            FinalDecision::Scrap => "SCRAP",
        }
    }

    /// Whether the read was retained.
    pub fn is_pass(&self) -> bool {
        matches!(self, FinalDecision::Pass)
    }
}

/// Stable reason codes describing the evaluated criteria of a read.
///
/// A read carries every applicable code in evaluation order, including the failed
/// criteria that motivated a fallback or a rescue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecisionReason {
    // -- guided single-exon reads --
    /// Reciprocal overlap with a reference exon met `min_overlap_frac`.
    ReferenceOverlapPass,
    /// Reciprocal overlap with a reference exon was below `min_overlap_frac`.
    LowReferenceOverlap,

    // -- guided multi-exon reads --
    /// Query-junction match fraction met `min_junction_frac`.
    JunctionSupportPass,
    /// Query-junction match fraction was below `min_junction_frac`.
    LowJunctionSupport,
    /// No multi-exon reference was available for junction comparison.
    NoComparableReferenceJunctions,

    // -- boundary matching --
    /// An exon boundary shared with a reference rescued the read.
    BoundaryMatchRescue,
    /// No exon boundary was shared with any reference.
    NoBoundaryMatch,

    // -- splice-site scores --
    /// The median splice score rescued a read that failed boundary matching.
    SpliceScoreRescue,
    /// The median splice score met `min_splice_score`.
    SpliceScorePass,
    /// The median splice score was below `min_splice_score`.
    LowSpliceScore,
    /// No usable median splice score existed for the read.
    SpliceScoreUnavailable,

    // -- guided to de novo redirection --
    /// The read shares no exon with any reference and was evaluated de novo.
    NoReferenceExonOverlap,

    // -- de novo component checks --
    /// The read is alone in its component.
    SingleReadComponent,
    /// The component holds fewer reads than `min_cluster_support`.
    LowComponentSupport,

    // -- de novo multi-exon reads --
    /// The read belongs to an intron-chain cluster of at least `min_cluster_support` reads.
    DominantIntronChainCluster,
    /// The read's intron-chain cluster is smaller than `min_cluster_support`.
    LowIntronChainClusterSupport,
    /// Per-intron support across the component rescued the read.
    IntronSupportRescue,
    /// Per-intron support across the component was below `min_intron_support_frac`.
    LowIntronSupport,
    /// No other read of the component carried a junction, so support is unmeasurable.
    IntronSupportUnavailable,

    // -- de novo single-exon reads --
    /// The read belongs to a reciprocal-overlap cluster of at least `min_cluster_support` reads.
    SupportedSingleExonCluster,
    /// The read's reciprocal-overlap cluster is smaller than `min_cluster_support`.
    LowSingleExonClusterSupport,
}

impl DecisionReason {
    /// Stable uppercase snake-case code used in the `reasons` column.
    pub fn as_str(&self) -> &'static str {
        match self {
            DecisionReason::ReferenceOverlapPass => "REFERENCE_OVERLAP_PASS",
            DecisionReason::LowReferenceOverlap => "LOW_REFERENCE_OVERLAP",
            DecisionReason::JunctionSupportPass => "JUNCTION_SUPPORT_PASS",
            DecisionReason::LowJunctionSupport => "LOW_JUNCTION_SUPPORT",
            DecisionReason::NoComparableReferenceJunctions => "NO_COMPARABLE_REFERENCE_JUNCTIONS",
            DecisionReason::BoundaryMatchRescue => "BOUNDARY_MATCH_RESCUE",
            DecisionReason::NoBoundaryMatch => "NO_BOUNDARY_MATCH",
            DecisionReason::SpliceScoreRescue => "SPLICE_SCORE_RESCUE",
            DecisionReason::SpliceScorePass => "SPLICE_SCORE_PASS",
            DecisionReason::LowSpliceScore => "LOW_SPLICE_SCORE",
            DecisionReason::SpliceScoreUnavailable => "SPLICE_SCORE_UNAVAILABLE",
            DecisionReason::NoReferenceExonOverlap => "NO_REFERENCE_EXON_OVERLAP",
            DecisionReason::SingleReadComponent => "SINGLE_READ_COMPONENT",
            DecisionReason::LowComponentSupport => "LOW_COMPONENT_SUPPORT",
            DecisionReason::DominantIntronChainCluster => "DOMINANT_INTRON_CHAIN_CLUSTER",
            DecisionReason::LowIntronChainClusterSupport => "LOW_INTRON_CHAIN_CLUSTER_SUPPORT",
            DecisionReason::IntronSupportRescue => "INTRON_SUPPORT_RESCUE",
            DecisionReason::LowIntronSupport => "LOW_INTRON_SUPPORT",
            DecisionReason::IntronSupportUnavailable => "INTRON_SUPPORT_UNAVAILABLE",
            DecisionReason::SupportedSingleExonCluster => "SUPPORTED_SINGLE_EXON_CLUSTER",
            DecisionReason::LowSingleExonClusterSupport => "LOW_SINGLE_EXON_CLUSTER_SUPPORT",
        }
    }
}

/// Exonic architecture of a query read, derived from its intron chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadType {
    /// The read has no intron.
    SingleExon,
    /// The read has at least one intron.
    MultiExon,
}

impl ReadType {
    /// Classify a record by the presence of introns, matching the classifier routing.
    pub fn from_record(record: &GenePred) -> Self {
        if record.introns().is_empty() {
            ReadType::SingleExon
        } else {
            ReadType::MultiExon
        }
    }

    /// Stable textual representation used in the `read_type` column.
    pub fn as_str(&self) -> &'static str {
        match self {
            ReadType::SingleExon => "SINGLE_EXON",
            ReadType::MultiExon => "MULTI_EXON",
        }
    }
}

/// Initial program mode, as selected by the CLI arguments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    /// Reference-guided mode (`--all` or `--overlapping`).
    Guided,
    /// Self-guided mode (`--non-overlapping`).
    DeNovo,
}

impl RunMode {
    /// Stable textual representation used in the `mode` column.
    pub fn as_str(&self) -> &'static str {
        match self {
            RunMode::Guided => "GUIDED",
            RunMode::DeNovo => "DE_NOVO",
        }
    }
}

/// Classifier that produced the final decision of a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvaluationPath {
    /// Scored against the reference transcripts of its component.
    GuidedReference,
    /// Guided read redirected to the de novo classifier for lacking reference-exon overlap.
    GuidedToDeNovo,
    /// Classified by the de novo classifier.
    DeNovo,
}

impl EvaluationPath {
    /// Stable textual representation used in the `evaluation_path` column.
    pub fn as_str(&self) -> &'static str {
        match self {
            EvaluationPath::GuidedReference => "GUIDED_REFERENCE",
            EvaluationPath::GuidedToDeNovo => "GUIDED_TO_DE_NOVO",
            EvaluationPath::DeNovo => "DE_NOVO",
        }
    }
}

/// A threshold that was actually compared against a feature of the read.
///
/// Only the thresholds of the read's own evaluation path are recorded, so the
/// `applied_thresholds` column stays a faithful trace of the applied criteria.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AppliedThreshold {
    /// `--min-overlap-frac`
    MinOverlapFrac(f64),
    /// `--min-junction-frac`
    MinJunctionFrac(f64),
    /// `--junction-tolerance`
    JunctionTolerance(u64),
    /// `--end-tolerance`
    EndTolerance(u64),
    /// `--min-read-num-denovo`
    MinClusterSupport(usize),
    /// `--min-single-exon-support`
    MinSingleExonSupport(usize),
    /// `--min-intron-support-frac`
    MinIntronSupportFrac(f64),
    /// `--intron-support-threshold`
    IntronSupportThreshold(f64),
    /// `--min-splice-score`
    MinSpliceScore(f32),
}

impl AppliedThreshold {
    /// Render as `name=value`, with three decimal places for fractional thresholds.
    pub fn render(&self) -> String {
        match self {
            AppliedThreshold::MinOverlapFrac(v) => format!("min_overlap_frac={}", fmt_f64(*v)),
            AppliedThreshold::MinJunctionFrac(v) => format!("min_junction_frac={}", fmt_f64(*v)),
            AppliedThreshold::JunctionTolerance(v) => format!("junction_tolerance={}", v),
            AppliedThreshold::EndTolerance(v) => format!("end_tolerance={}", v),
            AppliedThreshold::MinClusterSupport(v) => format!("min_cluster_support={}", v),
            AppliedThreshold::MinSingleExonSupport(v) => {
                format!("min_single_exon_support={}", v)
            }
            AppliedThreshold::MinIntronSupportFrac(v) => {
                format!("min_intron_support_frac={}", fmt_f64(*v))
            }
            AppliedThreshold::IntronSupportThreshold(v) => {
                format!("intron_support_threshold={}", fmt_f64(*v))
            }
            AppliedThreshold::MinSpliceScore(v) => {
                format!("min_splice_score={}", fmt_f64(*v as f64))
            }
        }
    }
}

/// Structured classification result of a single query read.
///
/// Identity and context fields are always present; evidence fields are optional and
/// are rendered as `.` when the corresponding check was not part of the read's path.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadReport {
    // -- identity --
    /// BED name field of the query record.
    pub read_id: Option<Vec<u8>>,
    /// BED chromosome field.
    pub chrom: Vec<u8>,
    /// BED start (0-based, inclusive).
    pub start: u64,
    /// BED end (0-based, exclusive).
    pub end: u64,
    /// BED strand field.
    pub strand: Option<Strand>,
    /// Number of exons of the query record.
    pub exon_count: usize,
    /// Exonic architecture of the query record.
    pub read_type: ReadType,

    // -- context --
    /// Initial program mode.
    pub mode: RunMode,
    /// Classifier that produced the decision.
    pub evaluation_path: EvaluationPath,
    /// Number of query reads in the group used for the decision.
    pub group_size: usize,
    /// Number of reference transcripts in the original packed component.
    pub reference_count: usize,
    /// Size of the cluster assigned to the read, when clustering was performed.
    pub cluster_size: Option<usize>,

    // -- calculated evidence --
    /// Whether the read shares an exon with any reference transcript.
    pub overlaps_reference_exon: Option<bool>,
    /// Reference transcript that produced the best applicable score.
    pub best_reference_id: Option<Vec<u8>>,
    /// Best reciprocal overlap used for guided single-exon classification.
    pub best_reference_overlap: Option<f64>,
    /// Matched query junctions in the best reference comparison.
    pub junction_matches: Option<usize>,
    /// Total number of query junctions.
    pub junction_count: Option<usize>,
    /// Best query-junction match fraction.
    pub junction_fraction: Option<f64>,
    /// Whether a coherent structural feature was shared with a single reference, when
    /// boundary rescue was evaluated.
    pub boundary_match: Option<bool>,
    /// Which coherent structural feature carried the boundary rescue.
    pub boundary_evidence: Option<BoundaryEvidence>,
    /// Fraction of the read's introns supported across the component.
    pub intron_support_fraction: Option<f64>,
    /// Median splice-site score across the read's splice sites.
    ///
    /// This is the value compared against `min_splice_score`.
    pub median_splice_score: Option<f32>,
    /// Lowest splice-site score across the read's splice sites.
    ///
    /// Diagnostic only: it gates nothing, but it shows whether a passing median hides an
    /// uncovered or genuinely weak site.
    pub weakest_splice_score: Option<f32>,
    /// Thresholds compared during the read's evaluation.
    pub applied_thresholds: Vec<AppliedThreshold>,

    // -- outcome --
    /// Reason codes in evaluation order.
    pub reasons: Vec<DecisionReason>,
    /// Final decision.
    pub decision: FinalDecision,
}

impl ReadReport {
    /// Start a report for `record`, carrying identity and context only.
    ///
    /// Evidence fields are filled in by the classifier as checks are performed, and the
    /// outcome is sealed with [`ReadReport::finish`].
    pub fn new(
        record: &GenePred,
        mode: RunMode,
        evaluation_path: EvaluationPath,
        group_size: usize,
        reference_count: usize,
    ) -> Self {
        Self {
            read_id: record.name().map(|name| name.to_vec()),
            chrom: record.chrom().to_vec(),
            start: record.start(),
            end: record.end(),
            strand: record.strand(),
            exon_count: record.exons().len(),
            read_type: ReadType::from_record(record),
            mode,
            evaluation_path,
            group_size,
            reference_count,
            cluster_size: None,
            overlaps_reference_exon: None,
            best_reference_id: None,
            best_reference_overlap: None,
            junction_matches: None,
            junction_count: None,
            junction_fraction: None,
            boundary_match: None,
            boundary_evidence: None,
            intron_support_fraction: None,
            median_splice_score: None,
            weakest_splice_score: None,
            applied_thresholds: Vec::new(),
            reasons: Vec::new(),
            decision: FinalDecision::Scrap,
        }
    }

    /// Seal the report with its ordered reason codes and its final decision.
    pub fn finish(mut self, reasons: Vec<DecisionReason>, decision: FinalDecision) -> Self {
        self.reasons = reasons;
        self.decision = decision;
        self
    }

    /// Record a threshold that was compared during evaluation.
    pub fn apply(&mut self, threshold: AppliedThreshold) {
        self.applied_thresholds.push(threshold);
    }

    /// Render the `reasons` column.
    pub fn reasons_field(&self) -> String {
        if self.reasons.is_empty() {
            return MISSING.to_string();
        }
        self.reasons
            .iter()
            .map(|reason| reason.as_str())
            .collect::<Vec<_>>()
            .join("/")
    }

    /// Render the `applied_thresholds` column.
    pub fn thresholds_field(&self) -> String {
        if self.applied_thresholds.is_empty() {
            return MISSING.to_string();
        }
        self.applied_thresholds
            .iter()
            .map(|threshold| threshold.render())
            .collect::<Vec<_>>()
            .join(";")
    }

    /// Render the report as a single TSV data row, without a trailing newline.
    pub fn to_tsv_row(&self) -> String {
        let fields: [String; REPORT_COLUMNS.len()] = [
            self.read_id
                .as_deref()
                .map(sanitize)
                .unwrap_or_else(|| MISSING.to_string()),
            sanitize(&self.chrom),
            self.start.to_string(),
            self.end.to_string(),
            self.strand
                .map(|strand| strand.to_string())
                .unwrap_or_else(|| MISSING.to_string()),
            self.exon_count.to_string(),
            self.read_type.as_str().to_string(),
            self.mode.as_str().to_string(),
            self.evaluation_path.as_str().to_string(),
            self.group_size.to_string(),
            self.reference_count.to_string(),
            fmt_opt_usize(self.cluster_size),
            fmt_opt_bool(self.overlaps_reference_exon),
            self.best_reference_id
                .as_deref()
                .map(sanitize)
                .unwrap_or_else(|| MISSING.to_string()),
            fmt_opt_f64(self.best_reference_overlap),
            fmt_opt_usize(self.junction_matches),
            fmt_opt_usize(self.junction_count),
            fmt_opt_f64(self.junction_fraction),
            fmt_opt_bool(self.boundary_match),
            self.boundary_evidence
                .map(|evidence| evidence.as_str().to_string())
                .unwrap_or_else(|| MISSING.to_string()),
            fmt_opt_f64(self.intron_support_fraction),
            fmt_opt_f64(self.median_splice_score.map(f64::from)),
            fmt_opt_f64(self.weakest_splice_score.map(f64::from)),
            self.thresholds_field(),
            self.reasons_field(),
            self.decision.as_str().to_string(),
        ];

        fields.join("\t")
    }
}

/// A query read paired with the classification result that produced its BED line.
#[derive(Debug, Clone, PartialEq)]
pub struct ClassifiedRead {
    /// BED12 line of the query record.
    pub bed_line: Vec<u8>,
    /// Classification result of the query record.
    pub report: ReadReport,
}

impl ClassifiedRead {
    /// Pair `record`'s BED12 line with its sealed report.
    pub fn new(record: &GenePred, report: ReadReport) -> Self {
        Self {
            bed_line: record.to_bed::<Bed12>(),
            report,
        }
    }

    /// Final decision of the read.
    pub fn decision(&self) -> FinalDecision {
        self.report.decision
    }
}

/// Header line of `<prefix>.report.tsv`.
pub fn report_header() -> String {
    REPORT_COLUMNS.join("\t")
}

/// Sort report rows deterministically.
///
/// The primary key is `chrom, start, end, read_id, decision, reasons`. Remaining ties
/// are broken by the rendered row, so thread scheduling cannot change the file. Rows
/// that are equal in every column stay duplicated.
pub fn sort_reports(reports: &mut [ReadReport]) {
    reports.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
            .then_with(|| a.read_id.cmp(&b.read_id))
            .then_with(|| a.decision.as_str().cmp(b.decision.as_str()))
            .then_with(|| {
                a.reasons
                    .iter()
                    .map(|reason| reason.as_str())
                    .cmp(b.reasons.iter().map(|reason| reason.as_str()))
            })
            .then_with(|| a.to_tsv_row().cmp(&b.to_tsv_row()))
    });
}

/// Write the header and one row per report, in the given order.
pub fn write_reports<W: Write>(writer: &mut W, reports: &[ReadReport]) -> io::Result<()> {
    writer.write_all(report_header().as_bytes())?;
    writer.write_all(b"\n")?;

    for report in reports {
        writer.write_all(report.to_tsv_row().as_bytes())?;
        writer.write_all(b"\n")?;
    }

    writer.flush()
}

/// Escape a textual BED field so it cannot break the TSV layout.
///
/// Tabs, carriage returns and newlines become their escaped two-character form.
/// Invalid UTF-8 is replaced lossily, and an empty field becomes [`MISSING`].
fn sanitize(field: &[u8]) -> String {
    let text = String::from_utf8_lossy(field);
    let mut out = String::with_capacity(text.len());

    for ch in text.chars() {
        match ch {
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            '\n' => out.push_str("\\n"),
            _ => out.push(ch),
        }
    }

    if out.is_empty() {
        MISSING.to_string()
    } else {
        out
    }
}

/// Render a float with three decimal places, or [`MISSING`] when not finite.
fn fmt_f64(value: f64) -> String {
    if value.is_finite() {
        format!("{:.3}", value)
    } else {
        MISSING.to_string()
    }
}

/// Render an optional float, or [`MISSING`] when absent.
fn fmt_opt_f64(value: Option<f64>) -> String {
    value.map(fmt_f64).unwrap_or_else(|| MISSING.to_string())
}

/// Render an optional count, or [`MISSING`] when absent.
fn fmt_opt_usize(value: Option<usize>) -> String {
    value
        .map(|v| v.to_string())
        .unwrap_or_else(|| MISSING.to_string())
}

/// Render an optional boolean in lowercase, or [`MISSING`] when absent.
fn fmt_opt_bool(value: Option<bool>) -> String {
    match value {
        Some(true) => "true".to_string(),
        Some(false) => "false".to_string(),
        None => MISSING.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// A report with no evidence at all, built without touching the classifier.
    fn blank_report() -> ReadReport {
        ReadReport {
            read_id: Some(b"read1".to_vec()),
            chrom: b"chr1".to_vec(),
            start: 100,
            end: 600,
            strand: Some(Strand::Forward),
            exon_count: 3,
            read_type: ReadType::MultiExon,
            mode: RunMode::Guided,
            evaluation_path: EvaluationPath::GuidedReference,
            group_size: 8,
            reference_count: 3,
            cluster_size: None,
            overlaps_reference_exon: None,
            best_reference_id: None,
            best_reference_overlap: None,
            junction_matches: None,
            junction_count: None,
            junction_fraction: None,
            boundary_match: None,
            boundary_evidence: None,
            intron_support_fraction: None,
            median_splice_score: None,
            weakest_splice_score: None,
            applied_thresholds: Vec::new(),
            reasons: Vec::new(),
            decision: FinalDecision::Scrap,
        }
    }

    fn make_record(name: &[u8], start: u64, end: u64) -> GenePred {
        GenePred {
            chrom: b"chr1".to_vec(),
            start,
            end,
            name: Some(name.to_vec()),
            strand: Some(Strand::Reverse),
            thick_start: Some(start),
            thick_end: Some(end),
            block_count: Some(2),
            block_starts: Some(vec![start, start + 200]),
            block_ends: Some(vec![start + 100, end]),
            extras: HashMap::new(),
        }
    }

    // -----------------------------------------------------------------------
    // Header
    // -----------------------------------------------------------------------

    #[test]
    fn test_header_column_order() {
        let header = report_header();
        let columns: Vec<&str> = header.split('\t').collect();

        assert_eq!(columns.len(), 26);
        assert_eq!(
            columns,
            vec![
                "read_id",
                "chrom",
                "start",
                "end",
                "strand",
                "exon_count",
                "read_type",
                "mode",
                "evaluation_path",
                "group_size",
                "reference_count",
                "cluster_size",
                "overlaps_reference_exon",
                "best_reference_id",
                "best_reference_overlap",
                "junction_matches",
                "junction_count",
                "junction_fraction",
                "boundary_match",
                "boundary_evidence",
                "intron_support_fraction",
                "median_splice_score",
                "weakest_splice_score",
                "applied_thresholds",
                "reasons",
                "decision",
            ]
        );
    }

    #[test]
    fn test_header_final_columns_are_reasons_and_decision() {
        let header = report_header();
        let columns: Vec<&str> = header.split('\t').collect();

        assert_eq!(columns[columns.len() - 2], "reasons");
        assert_eq!(columns[columns.len() - 1], "decision");
    }

    #[test]
    fn test_header_has_no_comment_prefix() {
        let header = report_header();
        assert!(!header.starts_with('#'));
    }

    // -----------------------------------------------------------------------
    // Row rendering
    // -----------------------------------------------------------------------

    #[test]
    fn test_row_has_one_field_per_column() {
        let row = blank_report().to_tsv_row();
        assert_eq!(row.split('\t').count(), REPORT_COLUMNS.len());
    }

    #[test]
    fn test_missing_values_written_as_dot() {
        let report = blank_report().finish(
            vec![DecisionReason::NoComparableReferenceJunctions],
            FinalDecision::Scrap,
        );
        let row = report.to_tsv_row();
        let fields: Vec<&str> = row.split('\t').collect();

        // every optional evidence column is unavailable here
        for column in [
            "cluster_size",
            "overlaps_reference_exon",
            "best_reference_id",
            "best_reference_overlap",
            "junction_matches",
            "junction_count",
            "junction_fraction",
            "boundary_match",
            "intron_support_fraction",
            "median_splice_score",
            "applied_thresholds",
        ] {
            let idx = REPORT_COLUMNS.iter().position(|c| *c == column).unwrap();
            assert_eq!(fields[idx], ".", "column {} must be '.'", column);
        }
    }

    #[test]
    fn test_missing_read_id_written_as_dot() {
        let mut report = blank_report();
        report.read_id = None;
        assert_eq!(report.to_tsv_row().split('\t').next().unwrap(), ".");

        report.read_id = Some(Vec::new());
        assert_eq!(report.to_tsv_row().split('\t').next().unwrap(), ".");
    }

    #[test]
    fn test_strand_rendering() {
        let mut report = blank_report();
        let strand_idx = REPORT_COLUMNS.iter().position(|c| *c == "strand").unwrap();

        for (strand, expected) in [
            (Some(Strand::Forward), "+"),
            (Some(Strand::Reverse), "-"),
            (Some(Strand::Unknown), "."),
            (None, "."),
        ] {
            report.strand = strand;
            let row = report.to_tsv_row();
            let fields: Vec<&str> = row.split('\t').collect();
            assert_eq!(fields[strand_idx], expected);
        }
    }

    #[test]
    fn test_float_formatting_is_three_decimals() {
        let mut report = blank_report();
        report.best_reference_overlap = Some(0.5);
        report.junction_fraction = Some(2.0 / 3.0);
        report.intron_support_fraction = Some(0.0);
        report.median_splice_score = Some(0.12);

        let row = report.to_tsv_row();
        let fields: Vec<&str> = row.split('\t').collect();
        let at = |name: &str| {
            let idx = REPORT_COLUMNS.iter().position(|c| *c == name).unwrap();
            fields[idx]
        };

        assert_eq!(at("best_reference_overlap"), "0.500");
        assert_eq!(at("junction_fraction"), "0.667");
        assert_eq!(at("intron_support_fraction"), "0.000");
        assert_eq!(at("median_splice_score"), "0.120");
    }

    #[test]
    fn test_non_finite_floats_are_not_written_as_nan() {
        let mut report = blank_report();
        report.junction_fraction = Some(f64::NAN);
        report.intron_support_fraction = Some(f64::INFINITY);
        report.median_splice_score = Some(f32::NAN);

        let row = report.to_tsv_row();
        assert!(!row.contains("NaN"));
        assert!(!row.contains("inf"));
    }

    #[test]
    fn test_booleans_are_lowercase() {
        let mut report = blank_report();
        report.overlaps_reference_exon = Some(true);
        report.boundary_match = Some(false);

        let row = report.to_tsv_row();
        assert!(row.contains("\ttrue\t"));
        assert!(row.contains("\tfalse\t"));
        assert!(!row.contains("True"));
        assert!(!row.contains("False"));
    }

    #[test]
    fn test_decision_column_only_pass_or_scrap() {
        for (decision, expected) in [
            (FinalDecision::Pass, "PASS"),
            (FinalDecision::Scrap, "SCRAP"),
        ] {
            let report = blank_report().finish(vec![DecisionReason::LowIntronSupport], decision);
            let row = report.to_tsv_row();
            assert_eq!(row.rsplit('\t').next().unwrap(), expected);
        }
    }

    #[test]
    fn test_reasons_are_slash_delimited_in_order() {
        let report = blank_report().finish(
            vec![
                DecisionReason::LowIntronChainClusterSupport,
                DecisionReason::IntronSupportRescue,
                DecisionReason::SpliceScorePass,
            ],
            FinalDecision::Pass,
        );

        assert_eq!(
            report.reasons_field(),
            "LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_RESCUE/SPLICE_SCORE_PASS"
        );
        assert!(!report.reasons_field().contains(' '));
    }

    #[test]
    fn test_thresholds_are_semicolon_delimited() {
        let mut report = blank_report();
        report.apply(AppliedThreshold::MinIntronSupportFrac(0.5));
        report.apply(AppliedThreshold::IntronSupportThreshold(0.5));
        report.apply(AppliedThreshold::MinSpliceScore(0.5));

        assert_eq!(
            report.thresholds_field(),
            "min_intron_support_frac=0.500;intron_support_threshold=0.500;min_splice_score=0.500"
        );
        assert!(!report.thresholds_field().contains(' '));
    }

    #[test]
    fn test_integer_thresholds_render_without_decimals() {
        assert_eq!(
            AppliedThreshold::MinClusterSupport(5).render(),
            "min_cluster_support=5"
        );
        assert_eq!(
            AppliedThreshold::EndTolerance(0).render(),
            "end_tolerance=0"
        );
        assert_eq!(
            AppliedThreshold::JunctionTolerance(3).render(),
            "junction_tolerance=3"
        );
        assert_eq!(
            AppliedThreshold::MinOverlapFrac(0.5).render(),
            "min_overlap_frac=0.500"
        );
    }

    #[test]
    fn test_textual_fields_are_sanitized() {
        let mut report = blank_report();
        report.read_id = Some(b"bad\tname\nwith\rbreaks".to_vec());
        report.best_reference_id = Some(b"ref\t1".to_vec());

        let row = report.to_tsv_row();
        assert_eq!(row.split('\t').count(), REPORT_COLUMNS.len());
        assert!(!row.contains('\n'));
        assert!(!row.contains('\r'));
        assert!(row.starts_with("bad\\tname\\nwith\\rbreaks\t"));
        assert!(row.contains("ref\\t1"));
    }

    // -----------------------------------------------------------------------
    // Writing
    // -----------------------------------------------------------------------

    #[test]
    fn test_write_reports_emits_header_then_rows() {
        let reports = vec![
            blank_report().finish(
                vec![DecisionReason::JunctionSupportPass],
                FinalDecision::Pass,
            ),
            blank_report().finish(
                vec![DecisionReason::LowReferenceOverlap],
                FinalDecision::Scrap,
            ),
        ];

        let mut buffer = Vec::new();
        write_reports(&mut buffer, &reports).unwrap();
        let text = String::from_utf8(buffer).unwrap();
        let lines: Vec<&str> = text.lines().collect();

        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], report_header());
        assert!(lines[1].ends_with("\tPASS"));
        assert!(lines[2].ends_with("\tSCRAP"));
        assert!(text.ends_with('\n'));
    }

    #[test]
    fn test_write_reports_without_rows_emits_only_header() {
        let mut buffer = Vec::new();
        write_reports(&mut buffer, &[]).unwrap();

        assert_eq!(
            String::from_utf8(buffer).unwrap(),
            format!("{}\n", report_header())
        );
    }

    // -----------------------------------------------------------------------
    // Sorting
    // -----------------------------------------------------------------------

    #[test]
    fn test_sort_reports_orders_by_documented_key() {
        let mut a = blank_report();
        a.chrom = b"chr2".to_vec();
        a.read_id = Some(b"z".to_vec());

        let mut b = blank_report();
        b.start = 50;
        b.read_id = Some(b"b".to_vec());

        let mut c = blank_report();
        c.read_id = Some(b"a".to_vec());

        let mut reports = vec![a, b, c];
        sort_reports(&mut reports);

        let ids: Vec<Vec<u8>> = reports.iter().map(|r| r.read_id.clone().unwrap()).collect();
        assert_eq!(ids, vec![b"b".to_vec(), b"a".to_vec(), b"z".to_vec()]);
    }

    #[test]
    fn test_sort_reports_breaks_ties_on_decision_and_reasons() {
        let pass = blank_report().finish(
            vec![DecisionReason::JunctionSupportPass],
            FinalDecision::Pass,
        );
        let scrap = blank_report().finish(
            vec![
                DecisionReason::LowJunctionSupport,
                DecisionReason::NoBoundaryMatch,
            ],
            FinalDecision::Scrap,
        );

        let mut reports = vec![scrap.clone(), pass.clone()];
        sort_reports(&mut reports);

        assert_eq!(reports[0].decision, FinalDecision::Pass);
        assert_eq!(reports[1].decision, FinalDecision::Scrap);

        let mut same_decision = vec![
            blank_report().finish(
                vec![
                    DecisionReason::LowJunctionSupport,
                    DecisionReason::NoBoundaryMatch,
                ],
                FinalDecision::Scrap,
            ),
            blank_report().finish(
                vec![DecisionReason::LowReferenceOverlap],
                FinalDecision::Scrap,
            ),
        ];
        sort_reports(&mut same_decision);
        assert_eq!(
            same_decision[0].reasons_field(),
            "LOW_JUNCTION_SUPPORT/NO_BOUNDARY_MATCH"
        );
        assert_eq!(same_decision[1].reasons_field(), "LOW_REFERENCE_OVERLAP");
    }

    #[test]
    fn test_sort_reports_keeps_identical_rows_duplicated() {
        let report = blank_report().finish(
            vec![DecisionReason::JunctionSupportPass],
            FinalDecision::Pass,
        );
        let mut reports = vec![report.clone(), report.clone(), report];

        sort_reports(&mut reports);

        assert_eq!(reports.len(), 3);
        let rows: Vec<String> = reports.iter().map(|r| r.to_tsv_row()).collect();
        assert_eq!(rows[0], rows[1]);
        assert_eq!(rows[1], rows[2]);
    }

    #[test]
    fn test_sort_reports_is_order_independent() {
        let mut base = Vec::new();
        for i in 0..8u64 {
            let mut report = blank_report();
            report.read_id = Some(format!("read{}", i % 3).into_bytes());
            report.start = 100 + (i % 4) * 10;
            report.end = report.start + 500;
            base.push(report.finish(
                vec![DecisionReason::DominantIntronChainCluster],
                if i % 2 == 0 {
                    FinalDecision::Pass
                } else {
                    FinalDecision::Scrap
                },
            ));
        }

        let mut forward = base.clone();
        let mut reversed = base.clone();
        reversed.reverse();

        sort_reports(&mut forward);
        sort_reports(&mut reversed);

        let rows = |reports: &[ReadReport]| -> Vec<String> {
            reports.iter().map(|r| r.to_tsv_row()).collect()
        };
        assert_eq!(rows(&forward), rows(&reversed));
    }

    // -----------------------------------------------------------------------
    // Construction from records
    // -----------------------------------------------------------------------

    #[test]
    fn test_new_captures_record_identity() {
        let record = make_record(b"q1", 100, 600);
        let report = ReadReport::new(&record, RunMode::DeNovo, EvaluationPath::DeNovo, 6, 0);

        assert_eq!(report.read_id.as_deref(), Some(b"q1".as_slice()));
        assert_eq!(report.chrom, b"chr1".to_vec());
        assert_eq!(report.start, 100);
        assert_eq!(report.end, 600);
        assert_eq!(report.exon_count, 2);
        assert_eq!(report.read_type, ReadType::MultiExon);
        assert_eq!(report.group_size, 6);
        assert_eq!(report.reference_count, 0);
        assert!(report.reasons.is_empty());
    }

    #[test]
    fn test_classified_read_pairs_bed_line_with_report() {
        let record = make_record(b"q1", 100, 600);
        let report = ReadReport::new(&record, RunMode::DeNovo, EvaluationPath::DeNovo, 6, 0)
            .finish(
                vec![DecisionReason::SingleReadComponent],
                FinalDecision::Scrap,
            );
        let classified = ClassifiedRead::new(&record, report);

        assert_eq!(classified.bed_line, record.to_bed::<Bed12>());
        assert_eq!(classified.decision(), FinalDecision::Scrap);
        assert!(!classified.decision().is_pass());
    }
}
