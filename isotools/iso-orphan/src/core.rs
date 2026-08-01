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
//! Both, guided and self-guided, cover an extensive amount of curated orphan cases under
//! the assumption that they do not represent a valid source of evidence for transcription.
//! The process is heavily parallellized to offer fast performance on large datasets.

use std::path::PathBuf;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use crate::{
    cli::{Args, Mode},
    denovo::{cluster_multi_exon, cluster_single_exon},
    report::{
        sort_reports, write_reports, AppliedThreshold, ClassifiedRead, DecisionReason,
        EvaluationPath, FinalDecision, ReadReport, RunMode,
    },
    scoring::{
        best_boundary_match, best_junction_match, best_single_exon_match,
        overlaps_any_reference_exon, IntronSupport, ScoringParams,
    },
    splice::{SpliceEvidence, SpliceScoreProvider},
    utils::*,
};

use dashmap::DashMap;
use genepred::GenePred;
use packbed::Role;
use rayon::prelude::*;

pub type Components = DashMap<String, Vec<Vec<GenePred>>>;

/// Score below which BigWig splice values are not kept in memory.
///
/// Sites under this value cannot clear any sensible `min_splice_score`, so dropping them
/// costs nothing. The effective prefilter is capped by `min_splice_score` itself.
const SPLICE_SCORE_PREFILTER: f32 = 0.01;

/// Detects orphans in a dataset of reads
///
/// # Arguments
///
/// * `args` - The arguments passed to the program
///
/// # Returns
///
/// None
///
/// # Examples
///
/// ```text
/// use isotools::iso_orphan::__detect_orphans;
/// use isotools::cli::Args;
///
/// let args = Args {
///     bed: vec!["tests/data/mm39.bed".into()],
///     ref:  None,
///     all: false,
///     overlapping: false,
///     non_overlapping: false,
///     keep_orphans: false,
///     min_read_num_denovo: 5,
///     outdir: "tests/data".into(),
///     name: "test".into(),
///     threads: 1,
///     level: log::Level::Info,
/// };
///
/// __detect_orphans(args);
/// ```
pub fn __detect_orphans(args: Args) {
    log::info!("INFO: Detecting orphans in dataset -> {:?}", args.query);

    let mode = Mode::from(&args);

    let outdir = args.outdir.join("orphans");
    std::fs::create_dir_all(&outdir).unwrap_or_else(|_| {
        log::error!("ERROR: Could not create output directory -> {:?}", outdir);
        std::process::exit(1);
    });

    let overlap_mode = args
        .overlap_mode
        .parse::<packbed::OverlapType>()
        .unwrap_or_else(|_| {
            log::error!(
                "ERROR: Could not parse overlap mode -> {:?}",
                args.overlap_mode
            );
            std::process::exit(1);
        });

    let params = ScoringParams {
        junction_tolerance: args.junction_tolerance,
        min_junction_frac: args.min_junction_frac,
        min_overlap_frac: args.min_overlap_frac,
        min_cluster_support: args.min_read_num_denovo,
        min_single_exon_support: args
            .min_single_exon_support
            .unwrap_or(args.min_read_num_denovo),
        end_tolerance: args.end_tolerance,
        min_intron_support_frac: args.min_intron_support_frac,
        intron_support_threshold: args.intron_support_threshold,
        min_splice_score: args.min_splice_score,
    };

    // Load optional splice-site scores from BigWig files
    //
    // INFO: the load prefilter must stay at or below min_splice_score, otherwise a site
    // scoring just under the prefilter would be dropped and later read as 0.0
    let splice_prefilter = SPLICE_SCORE_PREFILTER.min(args.min_splice_score);
    let splice_scores = args.splicing_scores.as_ref().map(|dir| {
        SpliceScoreProvider::from_dir(dir, splice_prefilter).unwrap_or_else(|e| {
            log::error!("ERROR: Failed to load splice scores: {}", e);
            std::process::exit(1);
        })
    });

    match mode {
        Mode::Guided => {
            log::info!("INFO: Running using guided mode");

            let mut inputs = args.refs.unwrap_or_else(|| {
                log::error!("ERROR: No reference file provided while using reference guided mode!");
                std::process::exit(1);
            });

            let mut modes = std::iter::repeat_n(Role::Reference, inputs.len()).collect::<Vec<_>>();
            modes.extend(vec![Role::Query]);
            inputs.extend(vec![args.query]);

            // CASE: single-exon-component reads with non-single-exon reference refs [SOLVED]
            // INFO:  the following case is presented:
            //
            // [SOLVED: CDS overlap mode]
            // ref:   xxx------------XXXX---
            //        |||
            // read1: xxx------XXXXx
            //
            // The above case would illustrate a component where a non-single-exon
            // component read could repesent real transcription but not at a significant
            // level.
            //
            // For a clearer illustration go to mm39 chr8:71,358,478-71,360,769
            //
            // Note that the same case but with CDS overlap is not a problem and
            // represents a case of a real transcription of a shorter isoform or a
            // truncated one.
            let tracks = packbed::pack(inputs, modes, overlap_mode).unwrap_or_else(|e| {
                log::error!("ERROR: Failed to packed reads: {:?}", e);
                std::process::exit(1);
            });

            __process(
                tracks,
                &mode,
                &params,
                outdir,
                args.prefix,
                splice_scores.as_ref(),
            );
        }
        Mode::DeNovo => {
            log::info!("INFO: Running using non-guided mode");

            // CASE: single-exon overlap non-supported UTR exon [SOLVED: CDS overlap mode]
            // INFO: the following case is presented:
            //
            // read1:       xxxXXXXxxxxxxxxxx
            // read2:      xxxxxxxxxxXXXXXxxxx
            // read3: --------------xxx----------
            // read4: ---------------------------
            // read5: ---------------------------
            //
            // The above case would illustrate a component where single-exon
            // reads do belong to the component by touching a non-supported
            // additional UTR-exon. This does not represent reliable transcription.
            //
            // SOLVED: CDS overlap mode would make isolate read1 and read2 to their
            // own single-exon components that, by extension, should be discarded.
            //
            // For a clearer illustration go to HLpteAle1A HAP1_SUPER_1:113,320,574-113,325,638
            // or mm39 chr8:73,246,564-73,253,282 or mm39 chr8:73,318,599-73,322,867
            //
            // CASE: single-exon read overlaps UTR [SOLVED: CDS overlap mode]
            // INFO: the following case is presented:
            //
            // read1: xxXXXXxxxx
            // read2: -------xxxxXXXXX------XXX--
            // read3: ---------xxXXXXX------XXX--
            // read4: ---------xxxxxxxxxXX--XXX--
            // read5: ----------xXXXXX-----------
            //                   |||||      |||
            //
            // The above case would illustrate a component where a single-exon read
            // overlaps with a non-single-exon read UTR.
            //
            // For a clearer illustration go to mm39 chr8:72,735,773-72,743,045
            // or mm39 chr8:73,437,570-73,445,000
            //
            // SOLVED: CDS overlap mode would make isolate read1, then will be discarded
            // because it is a single-exon component.
            //
            // or its variant:
            //
            // multi-exon read (likely orphan) matches exact UTR
            //
            // read1: xxxxx------xXXXXxx
            // read2: xxxxx----------------XXXX-
            // read3: xxxxx----------------XXXX-
            // read4: xxxxx----------------XXXX-
            // read5: xxxxx----------------XXXX-
            //
            // For a clearer illustration go to mm39 chr8:73,373,118-73,377,702
            //
            // SOLVED: CDS overlap mode would make isolate read1, then will be discarded
            // because it is a single-exon component.
            let tracks = packbed::pack(vec![args.query], vec![Role::Query], overlap_mode)
                .unwrap_or_else(|e| {
                    log::error!("ERROR: Failed to packed reads: {:?}", e);
                    std::process::exit(1);
                });

            __process(
                tracks,
                &mode,
                &params,
                outdir,
                args.prefix,
                splice_scores.as_ref(),
            );
        }
    };
}

/// Processes the components of reads in parallel
///
/// Writes three files to `outdir`: `<filename>.hq.bed`, `<filename>.scraps.bed` and
/// `<filename>.report.tsv`. The BED files and the report are both derived from the
/// same [`ClassifiedRead`] values, so they cannot disagree.
///
/// # Arguments
///
/// * `tracks` - A vector of vectors of GenePred structs
/// * `mode` - A reference to a Mode enum
/// * `params` - A reference to a ScoringParams struct
/// * `outdir` - A PathBuf struct
/// * `filename` - A String struct
/// * `splice_scores` - A reference to a SpliceScoreProvider struct
///
/// # Example
///
/// ```text
/// let tracks = vec![vec![GenePred::default()]];
/// let mode = Mode::Guided;
/// let params = ScoringParams::default();
/// let outdir = PathBuf::default();
/// let filename = String::default();
/// let splice_scores = None;
///
/// __process(tracks, &mode, &params, outdir, filename, splice_scores);
/// ```
fn __process(
    tracks: Components,
    mode: &Mode,
    params: &ScoringParams,
    outdir: PathBuf,
    filename: String,
    splice_scores: Option<&SpliceScoreProvider>,
) {
    let accumulator = ParallelAccumulator::default();
    let counter = ParallelCounter::default();

    tracks.into_par_iter().for_each(|bucket| {
        let chr = bucket.0;
        let components = bucket.1;

        log::debug!(
            "DEBUG: {} components in bucket -> {:?}",
            components.len(),
            chr
        );

        counter.inc_components(components.len() as u32);

        __process_components(
            components,
            &accumulator,
            &counter,
            mode,
            params,
            splice_scores,
        );
    });

    log::info!("INFO: Orphans found: {}", accumulator.num_orphans());

    __report_stats(&counter);

    let pass = format!("{}.hq.bed", filename);
    let orphans = format!("{}.scraps.bed", &filename);
    let report = format!("{}.report.tsv", &filename);

    let mut p_writer = BufWriter::new(File::create(outdir.join(pass)).unwrap());
    let mut o_writer = BufWriter::new(File::create(outdir.join(orphans)).unwrap());

    // INFO: one row per processed query record; never deduplicated, since query names
    // are not unique and identical records must still be reported separately
    let mut reports = accumulator.take_reports();
    sort_reports(&mut reports);

    let mut r_writer = BufWriter::new(File::create(outdir.join(report)).unwrap());
    write_reports(&mut r_writer, &reports).unwrap();

    log::info!("INFO: Reads classified in report: {}", reports.len());

    accumulator.keep.into_iter().for_each(|line| {
        p_writer.write_all(&line).unwrap();
        p_writer.write_all(b"\n".as_slice()).unwrap();
    });

    accumulator.orphans.into_iter().for_each(|line| {
        o_writer.write_all(&line).unwrap();
        o_writer.write_all(b"\n".as_slice()).unwrap();
    });
}

/// Processes a component of reads
///
/// Reference and query records are separated, the shared de novo evidence of the
/// component is built once over **all** of its query reads, and every query read is then
/// classified against the reference first and against that evidence second.
///
/// # Arguments
///
/// * `components` - A vector of vectors of GenePred structs
/// * `accumulator` - A reference to a ParallelAccumulator struct
/// * `counter` - A reference to a ParallelCounter struct
/// * `mode` - A reference to a Mode enum
/// * `params` - A reference to a ScoringParams struct
/// * `splice_scores` - A reference to a SpliceScoreProvider struct
///
/// # Example
///
/// ```text
/// let components = vec![vec![GenePred::default()]];
/// let accumulator = ParallelAccumulator::default();
/// let counter = ParallelCounter::default();
/// let mode = Mode::Guided;
/// let params = ScoringParams::default();
/// let splice_scores = None;
///
/// __process_components(
///     components,
///     &accumulator,
///     &counter,
///     &mode,
///     &params,
///     splice_scores
/// );
/// ```
fn __process_components(
    components: Vec<Vec<GenePred>>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
    mode: &Mode,
    params: &ScoringParams,
    splice_scores: Option<&SpliceScoreProvider>,
) {
    components.into_par_iter().for_each(|component| {
        let (references, queries): (Vec<GenePred>, Vec<GenePred>) =
            component.into_iter().partition(is_reference);

        if queries.is_empty() {
            let projections: Vec<Option<&[u8]>> = references.iter().map(|p| p.name()).collect();
            log::trace!(
                "DEBUG: reference projection without reads -> {:?}",
                projections
            );
            return;
        }

        let run_mode = match mode {
            Mode::Guided => RunMode::Guided,
            Mode::DeNovo => RunMode::DeNovo,
        };

        log::debug!(
            "DEBUG: processing {} mode; component of {} queries and {} references",
            run_mode.as_str(),
            queries.len(),
            references.len()
        );

        // INFO: in self-guided mode the references are not consulted even when the input
        // happens to carry some
        let references: &[GenePred] = match mode {
            Mode::Guided => &references,
            Mode::DeNovo => &[],
        };

        accumulator.add(classify_component(
            references,
            &queries,
            run_mode,
            counter,
            params,
            splice_scores,
        ));
    });
}

/// Whether a packed record carries the reference role.
fn is_reference(record: &GenePred) -> bool {
    let Some(role) = record.get_extra(b"role") else {
        log::error!(
            "ERROR: Could not get role from record {:?}! Was the input packed by packbed?",
            record.name()
        );
        std::process::exit(1);
    };

    role.clone().into_scalar() == Some(b"reference".to_vec())
}

/// De novo evidence shared by every query read of one component.
///
/// Building this once per component, over all of its query reads, is what lets a read
/// that fails reference evaluation be judged against the same evidence as a read that
/// never had a reference to compare with. Restricting the evidence to the subset of reads
/// that failed a routing test would make a read's support depend on how many *other*
/// reads happened to touch a reference exon.
struct ComponentEvidence {
    /// Number of query reads in the component.
    component_size: usize,
    /// Indices, into the component's query reads, of the multi-exon reads.
    multi_exon: Vec<usize>,
    /// Intron-chain cluster size of each multi-exon read, indexed as `multi_exon`.
    chain_cluster_size: Vec<usize>,
    /// Junction support across the multi-exon reads of the component.
    intron_support: IntronSupport,
    /// Reciprocal-overlap cluster size of each single-exon read.
    ///
    /// Indexed as `single_exon`.
    overlap_cluster_size: Vec<usize>,
    /// Indices, into the component's query reads, of the single-exon reads.
    single_exon: Vec<usize>,
    /// Position of each query read within `multi_exon` or `single_exon`.
    local_index: Vec<usize>,
}

impl ComponentEvidence {
    /// Cluster the query reads of a component and count their junction support.
    fn build(queries: &[GenePred], params: &ScoringParams) -> Self {
        let (multi_exon, single_exon): (Vec<usize>, Vec<usize>) =
            (0..queries.len()).partition(|&idx| !queries[idx].introns().is_empty());

        let mut local_index = vec![0usize; queries.len()];
        for (local, &global) in multi_exon.iter().enumerate() {
            local_index[global] = local;
        }
        for (local, &global) in single_exon.iter().enumerate() {
            local_index[global] = local;
        }

        let me_reads: Vec<&GenePred> = multi_exon.iter().map(|&idx| &queries[idx]).collect();
        let se_reads: Vec<&GenePred> = single_exon.iter().map(|&idx| &queries[idx]).collect();

        let intron_support = IntronSupport::build(&me_reads, params.junction_tolerance);

        let mut chain_cluster_size = vec![0usize; me_reads.len()];
        for members in cluster_multi_exon(&me_reads, params.junction_tolerance).values() {
            for &local in members {
                chain_cluster_size[local] = members.len();
            }
        }

        let mut overlap_cluster_size = vec![0usize; se_reads.len()];
        for members in cluster_single_exon(&se_reads, params.min_overlap_frac).values() {
            for &local in members {
                overlap_cluster_size[local] = members.len();
            }
        }

        Self {
            component_size: queries.len(),
            multi_exon,
            chain_cluster_size,
            intron_support,
            overlap_cluster_size,
            single_exon,
            local_index,
        }
    }

    /// Size of the cluster assigned to query read `idx`.
    fn cluster_size(&self, idx: usize, is_multi_exon: bool) -> usize {
        let local = self.local_index[idx];
        if is_multi_exon {
            self.chain_cluster_size[local]
        } else {
            self.overlap_cluster_size[local]
        }
    }

    /// Fraction of query read `idx`'s junctions supported by the other reads.
    fn support_fraction(&self, idx: usize, threshold: f64) -> Option<f64> {
        self.intron_support
            .support_fraction(self.local_index[idx], threshold)
    }

    /// Number of multi-exon query reads in the component.
    fn multi_exon_count(&self) -> usize {
        self.multi_exon.len()
    }

    /// Number of single-exon query reads in the component.
    fn single_exon_count(&self) -> usize {
        self.single_exon.len()
    }
}

/// Outcome of comparing one query read against the references of its component.
enum ReferenceOutcome {
    /// The reference evidence is coherent and sufficient; the read is retained.
    Supported(Vec<DecisionReason>),
    /// The reference evidence is insufficient. The read is not discarded here: it
    /// continues to de novo evaluation carrying these reasons as context.
    Insufficient(Vec<DecisionReason>),
}

/// Classifies every query read of one component.
///
/// Reference support and de novo support are applied **sequentially** rather than as
/// alternative routes:
///
/// 1. coherent reference support → `PASS`
/// 2. otherwise, independent support from the other query reads → `PASS`
/// 3. otherwise → `SCRAP`
///
/// Reference overlap is therefore evidence, not a permanent routing gate: a novel isoform
/// at an annotated locus can still be retained on the strength of the other reads, and a
/// read with no reference overlap is not treated as a different kind of object.
///
/// # Arguments
///
/// * `references` - Reference transcripts in the component, empty in self-guided mode
/// * `queries` - Query reads to classify
/// * `mode` - Initial program mode, reported as-is
/// * `counter` - Stats counter
/// * `params` - Scoring parameters
/// * `splice_scores` - Optional splice-site score provider
fn classify_component(
    references: &[GenePred],
    queries: &[GenePred],
    mode: RunMode,
    counter: &ParallelCounter,
    params: &ScoringParams,
    splice_scores: Option<&SpliceScoreProvider>,
) -> Vec<ClassifiedRead> {
    let evidence = ComponentEvidence::build(queries, params);
    let mut classified: Vec<ClassifiedRead> = Vec::with_capacity(queries.len());

    log::debug!(
        "DEBUG: component of {} queries -> {} multi-exon, {} single-exon, {} junctions",
        evidence.component_size,
        evidence.multi_exon_count(),
        evidence.single_exon_count(),
        evidence.intron_support.junction_count(),
    );

    for (idx, read) in queries.iter().enumerate() {
        let is_multi_exon = !read.introns().is_empty();

        let mut report = ReadReport::new(
            read,
            mode,
            EvaluationPath::GuidedReference,
            evidence.component_size,
            references.len(),
        );

        // -------------------------------------------------------------------
        // 1. reference support
        // -------------------------------------------------------------------
        let mut context: Vec<DecisionReason> = Vec::new();

        if !references.is_empty() {
            let overlaps = overlaps_any_reference_exon(read, references);
            report.overlaps_reference_exon = Some(overlaps);

            if overlaps {
                match evaluate_reference(read, references, is_multi_exon, &mut report, params) {
                    ReferenceOutcome::Supported(reasons) => {
                        count_reference_pass(counter, is_multi_exon, &reasons);
                        classified.push(ClassifiedRead::new(
                            read,
                            report.finish(reasons, FinalDecision::Pass),
                        ));
                        continue;
                    }
                    ReferenceOutcome::Insufficient(reasons) => {
                        counter.inc_guided_fallback_denovo();
                        context = reasons;
                    }
                }
            } else {
                counter.inc_guided_oob_denovo(1);
                context.push(DecisionReason::NoReferenceExonOverlap);
            }

            report.evaluation_path = EvaluationPath::GuidedToDeNovo;
        } else {
            report.evaluation_path = EvaluationPath::DeNovo;
        }

        // -------------------------------------------------------------------
        // 2. independent support from the other query reads of the component
        // -------------------------------------------------------------------
        let (reasons, decision) = evaluate_de_novo(
            read,
            idx,
            is_multi_exon,
            &evidence,
            &mut report,
            context,
            counter,
            params,
            splice_scores,
        );

        classified.push(ClassifiedRead::new(read, report.finish(reasons, decision)));
    }

    assert_eq!(
        classified.len(),
        queries.len(),
        "Partition violation: classified({}) != input({})",
        classified.len(),
        queries.len()
    );

    classified
}

/// Evaluates one query read against the reference transcripts of its component.
///
/// Records every feature it computes on `report`, and returns whether the reference
/// evidence is sufficient to retain the read on its own.
fn evaluate_reference(
    read: &GenePred,
    references: &[GenePred],
    is_multi_exon: bool,
    report: &mut ReadReport,
    params: &ScoringParams,
) -> ReferenceOutcome {
    if !is_multi_exon {
        // Single-exon reads: score by reciprocal overlap with the best reference exon.
        //
        // CASE: single-exon read(s) with single-exon reference refs
        //
        // ref:  xxxxxxxxxxxxXXXXXxxxxxxx
        // read1:     xxxxxxxXXXXXxxxx
        // read2:  xxxxxxxxxxXXXXXx
        // read3:          xxXXXXXxxxx
        // read5: xxxxxXXXXXXXXXXXxxxxx
        //             ^^^^^^|||||
        //
        // Here, we would argue that all of the reads (being single-exon) match
        // at least one CDS coordinate + all of them are within reference boundaries
        //
        // For a clearer illustration go to mm39 chr3:61,269,596-61,278,844
        // or mm39 chr3:65,014,198-65,019,415 or mm39 chr3:117,366,942-117,374,421
        let best = best_single_exon_match(read, references);
        let overlap = best.reciprocal_overlap;

        report.best_reference_id = best.reference_id;
        report.best_reference_overlap = Some(overlap);
        report.apply(AppliedThreshold::MinOverlapFrac(params.min_overlap_frac));

        if overlap >= params.min_overlap_frac {
            log::debug!(
                "DEBUG: read {:?} single-exon overlap={:.3} >= {:.3} -> keep!",
                read.name(),
                overlap,
                params.min_overlap_frac,
            );
            return ReferenceOutcome::Supported(vec![DecisionReason::ReferenceOverlapPass]);
        }

        log::debug!(
            "DEBUG: read {:?} single-exon overlap={:.3} < {:.3} -> de novo!",
            read.name(),
            overlap,
            params.min_overlap_frac,
        );
        return ReferenceOutcome::Insufficient(vec![DecisionReason::LowReferenceOverlap]);
    }

    // Multi-exon reads: score by splice junction agreement with individual
    // reference transcripts.
    //
    // CASE: multi-exon read with multi-exon reference
    //
    // ref:      xxxxxxXXXXXXX----XXXXX----
    //               ^^||||^^^
    // read1: xxxxxxxXXXXXXXXX----XXXXXxx -> shorter isoform
    // read2:       xxxXXXXXXX----XXXXX----
    // read3:     xxxxxXXXXXXX----XXXXX----
    //
    // For a clearer illustration go to mm39 chr3:65,435,178-65,440,499
    let Some(best) = best_junction_match(read, references, params.junction_tolerance) else {
        // INFO: every reference in the component is single-exon, so there is no junction
        // to compare against. A single shared coordinate is not transcript support for a
        // spliced read, so it is evaluated de novo instead
        log::debug!(
            "DEBUG: read {:?} multi-exon vs single-exon refs only -> de novo!",
            read.name(),
        );
        report.junction_count = Some(read.introns().len());
        return ReferenceOutcome::Insufficient(vec![
            DecisionReason::NoComparableReferenceJunctions,
        ]);
    };

    let score = best.score;
    report.best_reference_id = best.reference_id;
    report.junction_matches = Some(score.junction_matches);
    report.junction_count = Some(score.query_junction_count);
    report.junction_fraction = Some(score.query_junction_frac);
    report.apply(AppliedThreshold::MinJunctionFrac(params.min_junction_frac));

    if score.passes(params.min_junction_frac) {
        log::debug!(
            "DEBUG: read {:?} multi-exon junction_frac={:.3} ({}/{} matches) -> keep!",
            read.name(),
            score.query_junction_frac,
            score.junction_matches,
            score.query_junction_count,
        );
        return ReferenceOutcome::Supported(vec![DecisionReason::JunctionSupportPass]);
    }

    // The junction fraction is too low. A rescue is still possible, but only on a
    // complete structural feature shared with one single reference transcript.
    report.apply(AppliedThreshold::JunctionTolerance(
        params.junction_tolerance,
    ));
    report.apply(AppliedThreshold::EndTolerance(params.end_tolerance));

    let rescue = best_boundary_match(
        read,
        references,
        params.junction_tolerance,
        params.end_tolerance,
    );
    report.boundary_match = Some(rescue.is_some());

    match rescue {
        Some(rescue) => {
            log::debug!(
                "DEBUG: read {:?} multi-exon junction_frac={:.3} low, but coherent {} with {:?} -> boundary keep!",
                read.name(),
                score.query_junction_frac,
                rescue.evidence.as_str(),
                rescue.reference_id.as_ref().map(|id| String::from_utf8_lossy(id).to_string()),
            );
            report.boundary_evidence = Some(rescue.evidence);
            report.best_reference_id = rescue.reference_id;
            ReferenceOutcome::Supported(vec![
                DecisionReason::LowJunctionSupport,
                DecisionReason::BoundaryMatchRescue,
            ])
        }
        None => {
            log::debug!(
                "DEBUG: read {:?} multi-exon junction_frac={:.3}, no coherent reference feature -> de novo!",
                read.name(),
                score.query_junction_frac,
            );
            ReferenceOutcome::Insufficient(vec![
                DecisionReason::LowJunctionSupport,
                DecisionReason::NoBoundaryMatch,
            ])
        }
    }
}

/// Evaluates one query read against the independent evidence of its component.
///
/// Multi-exon and single-exon reads are treated differently on purpose. A multi-exon read
/// carries junctions, which are structural claims that other reads and the genome
/// sequence can corroborate, so low abundance means it needs stronger evidence rather
/// than automatic rejection. A single-exon read carries no junction at all, so abundance
/// remains its only evidence and stays terminal.
///
/// `context` holds the reasons already established upstream, such as a failed reference
/// comparison, and is prepended to the reasons produced here.
#[allow(clippy::too_many_arguments)]
fn evaluate_de_novo(
    read: &GenePred,
    idx: usize,
    is_multi_exon: bool,
    evidence: &ComponentEvidence,
    report: &mut ReadReport,
    context: Vec<DecisionReason>,
    counter: &ParallelCounter,
    params: &ScoringParams,
    splice_scores: Option<&SpliceScoreProvider>,
) -> (Vec<DecisionReason>, FinalDecision) {
    let mut reasons = context;

    // INFO: component size is a grouping property, not evidence for a transcript. It is
    // reported as context so a low-abundance locus stays visible, but it never decides
    // the outcome on its own
    if evidence.component_size <= 1 {
        counter.inc_reads_in_single_read_component();
        reasons.push(DecisionReason::SingleReadComponent);
    } else if evidence.component_size < params.min_cluster_support {
        counter.inc_reads_in_low_support_component();
        reasons.push(DecisionReason::LowComponentSupport);
    }

    let cluster_size = evidence.cluster_size(idx, is_multi_exon);
    report.cluster_size = Some(cluster_size);

    if !is_multi_exon {
        // -------------------------------------------------------------------
        // Single-exon reads: reciprocal-overlap cluster abundance only
        // -------------------------------------------------------------------
        //
        // CASE: single-exon reads following same/diff exonic pattern
        //
        // read1:       xxxXXXXXXXxxx
        // read2:      xxxxxxxxxxXXXXXxxxx
        //
        // For a clearer illustration go to HLpteAle1A HAP1_SUPER_1:112,960,536-112,964,453
        // or mm39 chr8:72,042,766-72,051,539 or mm39 chr8:73,174,008-73,179,856
        report.apply(AppliedThreshold::MinSingleExonSupport(
            params.min_single_exon_support,
        ));
        report.apply(AppliedThreshold::MinOverlapFrac(params.min_overlap_frac));

        if cluster_size >= params.min_single_exon_support {
            counter.inc_denovo_se_cluster_keep();
            reasons.push(DecisionReason::SupportedSingleExonCluster);
            return (reasons, FinalDecision::Pass);
        }

        counter.inc_denovo_se_cluster_orphan();
        reasons.push(DecisionReason::LowSingleExonClusterSupport);
        return (reasons, FinalDecision::Scrap);
    }

    // -----------------------------------------------------------------------
    // Multi-exon reads: intron chain -> junction support -> splice evidence
    // -----------------------------------------------------------------------
    report.apply(AppliedThreshold::MinClusterSupport(
        params.min_cluster_support,
    ));

    if cluster_size >= params.min_cluster_support {
        // Strong evidence: dominant intron-chain cluster → unconditional keep
        counter.inc_denovo_me_cluster_keep();
        reasons.push(DecisionReason::DominantIntronChainCluster);
        return (reasons, FinalDecision::Pass);
    }

    reasons.push(DecisionReason::LowIntronChainClusterSupport);

    let support = evidence.support_fraction(idx, params.intron_support_threshold);
    report.intron_support_fraction = support;
    report.apply(AppliedThreshold::MinIntronSupportFrac(
        params.min_intron_support_frac,
    ));
    report.apply(AppliedThreshold::IntronSupportThreshold(
        params.intron_support_threshold,
    ));

    match support {
        Some(fraction) if fraction >= params.min_intron_support_frac => {
            // The other reads of the component corroborate this read's junctions. A
            // splice-site model may still veto it, but an unavailable score does not:
            // absent evidence is not evidence of absence.
            reasons.push(DecisionReason::IntronSupportRescue);

            match splice_evidence_of(read, report, counter, params, splice_scores) {
                Some(evidence) if evidence.passes(params.min_splice_score) => {
                    counter.inc_denovo_me_intron_keep();
                    reasons.push(DecisionReason::SpliceScorePass);
                    (reasons, FinalDecision::Pass)
                }
                Some(_) => {
                    counter.inc_denovo_me_splice_orphan();
                    reasons.push(DecisionReason::LowSpliceScore);
                    (reasons, FinalDecision::Scrap)
                }
                None => {
                    counter.inc_denovo_me_intron_keep();
                    reasons.push(DecisionReason::SpliceScoreUnavailable);
                    (reasons, FinalDecision::Pass)
                }
            }
        }
        support => {
            // Structural support is missing, either because the other reads do not
            // corroborate the junctions or because there is no other read to ask. The
            // read is not discarded yet: independent sequence evidence can still carry
            // it, which is what keeps a low-abundance locus recoverable.
            reasons.push(match support {
                Some(_) => DecisionReason::LowIntronSupport,
                None => DecisionReason::IntronSupportUnavailable,
            });

            match splice_evidence_of(read, report, counter, params, splice_scores) {
                Some(evidence) if evidence.passes(params.min_splice_score) => {
                    counter.inc_denovo_me_splice_keep();
                    reasons.push(DecisionReason::SpliceScoreRescue);
                    (reasons, FinalDecision::Pass)
                }
                Some(_) => {
                    counter.inc_denovo_me_splice_orphan();
                    reasons.push(DecisionReason::LowSpliceScore);
                    (reasons, FinalDecision::Scrap)
                }
                None => {
                    counter.inc_denovo_me_cluster_orphan();
                    reasons.push(DecisionReason::SpliceScoreUnavailable);
                    (reasons, FinalDecision::Scrap)
                }
            }
        }
    }
}

/// Splice-site evidence for `read`, recording it on `report`.
fn splice_evidence_of(
    read: &GenePred,
    report: &mut ReadReport,
    counter: &ParallelCounter,
    params: &ScoringParams,
    splice_scores: Option<&SpliceScoreProvider>,
) -> Option<SpliceEvidence> {
    report.apply(AppliedThreshold::MinSpliceScore(params.min_splice_score));

    let evidence = splice_scores.and_then(|scores| scores.splice_evidence(read));

    match evidence {
        Some(evidence) => {
            report.median_splice_score = Some(evidence.median);
            report.weakest_splice_score = Some(evidence.weakest);
        }
        // INFO: an unstranded read, or one outside BigWig coverage, can never be rescued
        // or vetoed by splice evidence; counting these explains scrap totals
        None => counter.inc_splice_score_unavailable(),
    }

    evidence
}

/// Increments the counter matching a reference-supported outcome.
fn count_reference_pass(
    counter: &ParallelCounter,
    is_multi_exon: bool,
    reasons: &[DecisionReason],
) {
    if !is_multi_exon {
        counter.inc_guided_se_keep();
    } else if reasons.contains(&DecisionReason::BoundaryMatchRescue) {
        counter.inc_guided_me_boundary_keep();
    } else {
        counter.inc_guided_me_junction_keep();
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::ReadType;
    use crate::scoring::{BoundaryEvidence, ScoringParams};
    use crate::splice::SpliceScoreProvider;
    use crate::utils::{ParallelAccumulator, ParallelCounter};
    use std::collections::HashMap;

    fn make_record(
        name: &[u8],
        start: u64,
        end: u64,
        block_starts: Vec<u64>,
        block_ends: Vec<u64>,
    ) -> GenePred {
        GenePred {
            chrom: b"chr1".to_vec(),
            start,
            end,
            name: Some(name.to_vec()),
            strand: None,
            thick_start: Some(start),
            thick_end: Some(end),
            block_count: Some(block_starts.len() as u32),
            block_starts: Some(block_starts),
            block_ends: Some(block_ends),
            extras: HashMap::new(),
        }
    }

    fn make_stranded_record(
        name: &[u8],
        start: u64,
        end: u64,
        block_starts: Vec<u64>,
        block_ends: Vec<u64>,
        strand: genepred::Strand,
    ) -> GenePred {
        GenePred {
            chrom: b"chr1".to_vec(),
            start,
            end,
            name: Some(name.to_vec()),
            strand: Some(strand),
            thick_start: Some(start),
            thick_end: Some(end),
            block_count: Some(block_starts.len() as u32),
            block_starts: Some(block_starts),
            block_ends: Some(block_ends),
            extras: HashMap::new(),
        }
    }

    /// Classify a guided component.
    fn guided(
        references: Vec<GenePred>,
        queries: Vec<GenePred>,
        counter: &ParallelCounter,
        params: &ScoringParams,
        splice_scores: Option<&SpliceScoreProvider>,
    ) -> Vec<ClassifiedRead> {
        classify_component(
            &references,
            &queries,
            RunMode::Guided,
            counter,
            params,
            splice_scores,
        )
    }

    /// Classify a self-guided component.
    fn de_novo(
        queries: &[GenePred],
        counter: &ParallelCounter,
        params: &ScoringParams,
        splice_scores: Option<&SpliceScoreProvider>,
    ) -> Vec<ClassifiedRead> {
        classify_component(
            &[],
            queries,
            RunMode::DeNovo,
            counter,
            params,
            splice_scores,
        )
    }

    /// Split classified reads into the BED lines of the two outputs, exactly as the
    /// accumulator does.
    fn split(classified: &[ClassifiedRead]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
        let mut keep = Vec::new();
        let mut orphans = Vec::new();

        for read in classified {
            if read.decision().is_pass() {
                keep.push(read.bed_line.clone());
            } else {
                orphans.push(read.bed_line.clone());
            }
        }

        (keep, orphans)
    }

    /// The single report of the read named `name`.
    fn report_of<'a>(classified: &'a [ClassifiedRead], name: &[u8]) -> &'a ReadReport {
        let mut found = classified
            .iter()
            .map(|read| &read.report)
            .filter(|report| report.read_id.as_deref() == Some(name));

        let report = found.next().unwrap_or_else(|| {
            panic!(
                "no report for read {:?}",
                String::from_utf8_lossy(name).to_string()
            )
        });
        assert!(
            found.next().is_none(),
            "more than one report for read {:?}",
            String::from_utf8_lossy(name).to_string()
        );

        report
    }

    /// A splice-score provider scoring every splice site of `introns` with `score`.
    ///
    /// BigWig scores are 0-based per-base values: the donor of an intron `(start, end)`
    /// sits on its first base, `start`, and the acceptor on its last base, `end - 1`.
    fn provider_for(introns: &[(u64, u64)], score: f32) -> SpliceScoreProvider {
        use dashmap::DashMap;

        let donor = DashMap::new();
        let acceptor = DashMap::new();
        let mut d = HashMap::new();
        let mut a = HashMap::new();

        for &(start, end) in introns {
            d.insert(start, score);
            a.insert(end - 1, score);
        }

        donor.insert("chr1".to_string(), d);
        acceptor.insert("chr1".to_string(), a);

        SpliceScoreProvider::from_maps(donor, DashMap::new(), acceptor, DashMap::new())
    }

    /// Reference with exons [100,200), [300,400), [500,600); introns (200,300), (400,500).
    fn reference() -> GenePred {
        make_record(b"ref1", 100, 600, vec![100, 300, 500], vec![200, 400, 600])
    }

    /// `min_cluster_support` reads reproducing the reference intron chain.
    fn dominant_reads(count: u64) -> Vec<GenePred> {
        (0..count)
            .map(|i| {
                make_record(
                    format!("dom{}", i).as_bytes(),
                    100 + i * 2,
                    600,
                    vec![100 + i * 2, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect()
    }

    // -----------------------------------------------------------------------
    // Guided reference evaluation
    // -----------------------------------------------------------------------

    #[test]
    fn test_guided_single_exon_overlap_pass() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let ref1 = make_record(b"ref1", 100, 200, vec![100], vec![200]);
        let query = make_record(b"q1", 100, 200, vec![100], vec![200]);

        let classified = guided(vec![ref1], vec![query], &counter, &params, None);
        let report = report_of(&classified, b"q1");

        assert_eq!(report.decision, FinalDecision::Pass);
        assert_eq!(report.reasons_field(), "REFERENCE_OVERLAP_PASS");
        assert_eq!(report.read_type, ReadType::SingleExon);
        assert_eq!(report.evaluation_path, EvaluationPath::GuidedReference);
        assert_eq!(report.overlaps_reference_exon, Some(true));
        assert_eq!(
            report.best_reference_id.as_deref(),
            Some(b"ref1".as_slice())
        );
        assert_eq!(report.best_reference_overlap, Some(1.0));
        assert_eq!(report.thresholds_field(), "min_overlap_frac=0.500");
    }

    #[test]
    fn test_guided_multi_exon_junction_pass() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let classified = guided(vec![reference()], vec![query], &counter, &params, None);
        let report = report_of(&classified, b"q1");

        assert_eq!(report.decision, FinalDecision::Pass);
        assert_eq!(report.reasons_field(), "JUNCTION_SUPPORT_PASS");
        assert_eq!(report.junction_matches, Some(2));
        assert_eq!(report.junction_count, Some(2));
        assert_eq!(report.junction_fraction, Some(1.0));
        assert_eq!(report.boundary_match, None);
        assert_eq!(report.boundary_evidence, None);
        assert_eq!(report.thresholds_field(), "min_junction_frac=0.500");
    }

    // -----------------------------------------------------------------------
    // Coherent boundary rescue
    // -----------------------------------------------------------------------

    #[test]
    fn test_guided_rescue_requires_a_complete_shared_junction() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            // two query junctions, only one of which matches -> 0.5 is not enough
            min_junction_frac: 0.9,
            ..ScoringParams::default()
        };

        // shares the complete junction (200,300) with ref1 and adds a novel one
        let query = make_record(b"q1", 100, 700, vec![100, 300, 650], vec![200, 400, 700]);

        let classified = guided(vec![reference()], vec![query], &counter, &params, None);
        let report = report_of(&classified, b"q1");

        assert_eq!(report.decision, FinalDecision::Pass);
        assert_eq!(
            report.reasons_field(),
            "LOW_JUNCTION_SUPPORT/BOUNDARY_MATCH_RESCUE"
        );
        assert_eq!(report.boundary_match, Some(true));
        assert_eq!(
            report.boundary_evidence,
            Some(BoundaryEvidence::SpliceJunction)
        );
        assert_eq!(
            report.best_reference_id.as_deref(),
            Some(b"ref1".as_slice())
        );
        assert_eq!(
            report.thresholds_field(),
            "min_junction_frac=0.900;junction_tolerance=0;end_tolerance=0"
        );
    }

    #[test]
    fn test_guided_rescue_on_a_whole_shared_exon() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        // reuses the whole internal exon [300,400) but shares no junction
        let query = make_record(b"q1", 300, 5100, vec![300, 5000], vec![400, 5100]);

        let classified = guided(vec![reference()], vec![query], &counter, &params, None);
        let report = report_of(&classified, b"q1");

        assert_eq!(report.decision, FinalDecision::Pass);
        assert_eq!(report.boundary_evidence, Some(BoundaryEvidence::ExonPair));
    }

    #[test]
    fn test_guided_rescue_on_both_transcript_ends() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            end_tolerance: 20,
            ..ScoringParams::default()
        };

        // same locus as ref1 within tolerance, but an entirely novel intron chain
        let query = make_record(b"q1", 110, 590, vec![110, 350, 450], vec![250, 420, 590]);

        let classified = guided(vec![reference()], vec![query], &counter, &params, None);
        let report = report_of(&classified, b"q1");

        assert_eq!(report.decision, FinalDecision::Pass);
        assert_eq!(
            report.boundary_evidence,
            Some(BoundaryEvidence::BothTranscriptEnds)
        );
    }

    #[test]
    fn test_guided_single_shared_coordinate_is_not_a_rescue() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        // one coordinate, 300, coincides with a reference exon start; nothing else does
        let query = make_record(b"q1", 300, 5100, vec![300, 5000], vec![450, 5100]);

        let classified = guided(vec![reference()], vec![query], &counter, &params, None);
        let report = report_of(&classified, b"q1");

        assert_eq!(report.boundary_match, Some(false));
        assert_eq!(report.boundary_evidence, None);
        assert!(report.reasons.contains(&DecisionReason::NoBoundaryMatch));
        assert_eq!(
            report.evaluation_path,
            EvaluationPath::GuidedToDeNovo,
            "an insufficient reference comparison continues to de novo evaluation"
        );
    }

    #[test]
    fn test_guided_rescue_never_pools_incompatible_references() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        // ref_a supplies the query's first junction, ref_b the second; neither reference
        // supplies both, and each ref contradicts the other
        let ref_a = make_record(b"ref_a", 100, 400, vec![100, 300], vec![200, 400]);
        let ref_b = make_record(b"ref_b", 100, 900, vec![100, 800], vec![700, 900]);
        // query junctions: (200,300) from ref_a and (700,800) from ref_b
        let query = make_record(b"q1", 100, 900, vec![100, 300, 800], vec![200, 700, 900]);

        let classified = guided(vec![ref_a, ref_b], vec![query], &counter, &params, None);
        let report = report_of(&classified, b"q1");

        // each reference on its own contributes exactly one of the two query junctions,
        // which is a coherent partial-transcript relationship with that reference
        assert_eq!(report.junction_matches, Some(1));
        assert_eq!(
            report.junction_fraction,
            Some(0.5),
            "no single reference explains both junctions"
        );
    }

    // -----------------------------------------------------------------------
    // Multi-exon query against single-exon references only
    // -----------------------------------------------------------------------

    #[test]
    fn test_guided_multi_exon_vs_single_exon_refs_goes_to_de_novo() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let se_ref = make_record(b"ref1", 1000, 2000, vec![1000], vec![2000]);
        // shares the 1000 boundary but no junction: not transcript support
        let query = make_stranded_record(
            b"q1",
            1000,
            1400,
            vec![1000, 1300],
            vec![1200, 1400],
            genepred::Strand::Forward,
        );

        let classified = guided(vec![se_ref], vec![query], &counter, &params, None);
        let report = report_of(&classified, b"q1");

        assert_eq!(report.evaluation_path, EvaluationPath::GuidedToDeNovo);
        assert!(report
            .reasons
            .contains(&DecisionReason::NoComparableReferenceJunctions));
        assert_eq!(report.boundary_match, None);
        assert_eq!(
            report.decision,
            FinalDecision::Scrap,
            "a lone read with no other evidence is still discarded"
        );
    }

    #[test]
    fn test_guided_multi_exon_vs_single_exon_refs_rescued_by_splice_scores() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let se_ref = make_record(b"ref1", 1000, 2000, vec![1000], vec![2000]);
        let query = make_stranded_record(
            b"q1",
            1000,
            1400,
            vec![1000, 1300],
            vec![1200, 1400],
            genepred::Strand::Forward,
        );
        let provider = provider_for(&[(1200, 1300)], 0.9);

        let classified = guided(
            vec![se_ref],
            vec![query],
            &counter,
            &params,
            Some(&provider),
        );
        let report = report_of(&classified, b"q1");

        assert_eq!(report.decision, FinalDecision::Pass);
        assert_eq!(
            report.reasons_field(),
            "NO_COMPARABLE_REFERENCE_JUNCTIONS/SINGLE_READ_COMPONENT/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_UNAVAILABLE/SPLICE_SCORE_RESCUE"
        );
        assert_eq!(report.weakest_splice_score, Some(0.9));
    }

    // -----------------------------------------------------------------------
    // Guided failure falls through to de novo evidence
    // -----------------------------------------------------------------------

    #[test]
    fn test_guided_failure_rescued_by_other_query_reads() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 3,
            ..ScoringParams::default()
        };

        // three reads sharing a novel intron chain inside the reference locus; they
        // overlap reference exons, so the old routing sent them to guided-only evaluation.
        // Neither transcript end, nor any whole exon, agrees with the reference, so no
        // coherent reference rescue applies
        let novel: Vec<GenePred> = (0..3u64)
            .map(|i| {
                make_record(
                    format!("novel{}", i).as_bytes(),
                    150 + i,
                    580,
                    vec![150 + i, 350, 550],
                    vec![250, 450, 580],
                )
            })
            .collect();

        let classified = guided(vec![reference()], novel, &counter, &params, None);
        let report = report_of(&classified, b"novel0");

        assert_eq!(
            report.decision,
            FinalDecision::Pass,
            "a novel isoform at an annotated locus survives on query support"
        );
        assert_eq!(
            report.reasons_field(),
            "LOW_JUNCTION_SUPPORT/NO_BOUNDARY_MATCH/DOMINANT_INTRON_CHAIN_CLUSTER"
        );
        assert_eq!(report.evaluation_path, EvaluationPath::GuidedToDeNovo);
        assert_eq!(report.overlaps_reference_exon, Some(true));
        assert_eq!(report.cluster_size, Some(3));
    }

    #[test]
    fn test_guided_single_exon_failure_falls_through_to_de_novo() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 4,
            min_single_exon_support: 3,
            ..ScoringParams::default()
        };

        // four single-exon reads that overlap a reference exon too weakly to pass, but
        // form a strong reciprocal-overlap cluster among themselves
        let queries: Vec<GenePred> = (0..4u64)
            .map(|i| {
                make_record(
                    format!("se{}", i).as_bytes(),
                    150 + i,
                    400 + i,
                    vec![150 + i],
                    vec![400 + i],
                )
            })
            .collect();

        let classified = guided(vec![reference()], queries, &counter, &params, None);
        let report = report_of(&classified, b"se0");

        assert_eq!(report.decision, FinalDecision::Pass);
        assert_eq!(
            report.reasons_field(),
            "LOW_REFERENCE_OVERLAP/SUPPORTED_SINGLE_EXON_CLUSTER"
        );
        assert_eq!(report.cluster_size, Some(4));
    }

    #[test]
    fn test_guided_no_reference_exon_overlap_keeps_its_context() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        // two identical chains sitting inside the reference intron (200,300)
        let queries = vec![
            make_record(b"oob0", 210, 290, vec![210, 260], vec![240, 290]),
            make_record(b"oob1", 210, 290, vec![210, 260], vec![240, 290]),
        ];

        let classified = guided(vec![reference()], queries, &counter, &params, None);
        let report = report_of(&classified, b"oob0");

        assert_eq!(report.overlaps_reference_exon, Some(false));
        assert_eq!(report.evaluation_path, EvaluationPath::GuidedToDeNovo);
        assert_eq!(
            report.reasons_field(),
            "NO_REFERENCE_EXON_OVERLAP/DOMINANT_INTRON_CHAIN_CLUSTER"
        );
        assert_eq!(report.decision, FinalDecision::Pass);
    }

    #[test]
    fn test_de_novo_evidence_spans_the_whole_component() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 4,
            ..ScoringParams::default()
        };

        // three reads overlap reference exons and one does not, but all four share the
        // same novel intron chain; the evidence base must include all of them
        let mut queries: Vec<GenePred> = (0..3u64)
            .map(|i| {
                make_record(
                    format!("in{}", i).as_bytes(),
                    150 + i,
                    580,
                    vec![150 + i, 350, 550],
                    vec![250, 450, 580],
                )
            })
            .collect();
        queries.push(make_record(
            b"oob",
            210,
            580,
            vec![210, 350, 550],
            vec![250, 450, 580],
        ));

        let classified = guided(vec![reference()], queries, &counter, &params, None);

        for name in [b"in0".as_slice(), b"oob".as_slice()] {
            let report = report_of(&classified, name);
            assert_eq!(
                report.cluster_size,
                Some(4),
                "{:?} sees the whole component",
                String::from_utf8_lossy(name).to_string()
            );
            assert_eq!(report.decision, FinalDecision::Pass);
        }
    }

    // -----------------------------------------------------------------------
    // Component size is no longer a terminal classifier for multi-exon reads
    // -----------------------------------------------------------------------

    #[test]
    fn test_small_component_multi_exon_reads_can_be_retained() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            ..ScoringParams::default()
        };

        // three reads with a coherent shared structure in a component of three
        let reads = dominant_reads(3);
        let classified = de_novo(&reads, &counter, &params, None);
        let (keep, orphans) = split(&classified);

        assert_eq!(
            keep.len(),
            3,
            "coherent structure is retained below min_cluster_support"
        );
        assert!(orphans.is_empty());

        let report = report_of(&classified, b"dom0");
        assert_eq!(
            report.reasons_field(),
            "LOW_COMPONENT_SUPPORT/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_RESCUE/SPLICE_SCORE_UNAVAILABLE"
        );
        assert_eq!(
            report.intron_support_fraction,
            Some(1.0),
            "the other two reads corroborate both junctions"
        );
        assert_eq!(report.group_size, 3);
    }

    #[test]
    fn test_lone_multi_exon_read_rescued_by_splice_evidence() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let read = make_stranded_record(
            b"lone",
            100,
            600,
            vec![100, 300],
            vec![200, 600],
            genepred::Strand::Forward,
        );
        let provider = provider_for(&[(200, 300)], 0.8);

        let classified = de_novo(&[read], &counter, &params, Some(&provider));
        let report = report_of(&classified, b"lone");

        assert_eq!(report.decision, FinalDecision::Pass);
        assert_eq!(
            report.reasons_field(),
            "SINGLE_READ_COMPONENT/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_UNAVAILABLE/SPLICE_SCORE_RESCUE"
        );
        assert_eq!(
            report.intron_support_fraction, None,
            "there is no other read to draw support from"
        );
        assert_eq!(report.weakest_splice_score, Some(0.8));
    }

    #[test]
    fn test_lone_multi_exon_read_without_splice_evidence_is_scrapped() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let read = make_record(b"lone", 100, 600, vec![100, 300], vec![200, 600]);
        let classified = de_novo(&[read], &counter, &params, None);
        let report = report_of(&classified, b"lone");

        assert_eq!(report.decision, FinalDecision::Scrap);
        assert_eq!(
            report.reasons_field(),
            "SINGLE_READ_COMPONENT/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_UNAVAILABLE/SPLICE_SCORE_UNAVAILABLE"
        );
    }

    #[test]
    fn test_lone_multi_exon_read_with_weak_splice_evidence_is_scrapped() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let read = make_stranded_record(
            b"lone",
            100,
            600,
            vec![100, 300],
            vec![200, 600],
            genepred::Strand::Forward,
        );
        let provider = provider_for(&[(200, 300)], 0.1);

        let classified = de_novo(&[read], &counter, &params, Some(&provider));
        let report = report_of(&classified, b"lone");

        assert_eq!(report.decision, FinalDecision::Scrap);
        assert!(report.reasons_field().ends_with("LOW_SPLICE_SCORE"));
        assert_eq!(report.weakest_splice_score, Some(0.1));
    }

    #[test]
    fn test_single_exon_reads_keep_a_terminal_abundance_requirement() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_single_exon_support: 5,
            ..ScoringParams::default()
        };

        // three overlapping single-exon reads: no junction, so no orthogonal evidence
        let reads: Vec<GenePred> = (0..3u64)
            .map(|i| {
                make_record(
                    format!("se{}", i).as_bytes(),
                    100 + i,
                    300 + i,
                    vec![100 + i],
                    vec![300 + i],
                )
            })
            .collect();

        let classified = de_novo(&reads, &counter, &params, None);
        let (keep, orphans) = split(&classified);

        assert!(keep.is_empty());
        assert_eq!(orphans.len(), 3);

        let report = report_of(&classified, b"se0");
        assert_eq!(
            report.reasons_field(),
            "LOW_COMPONENT_SUPPORT/LOW_SINGLE_EXON_CLUSTER_SUPPORT"
        );
        assert_eq!(
            report.thresholds_field(),
            "min_single_exon_support=5;min_overlap_frac=0.500"
        );
    }

    #[test]
    fn test_single_exon_support_is_configurable_apart_from_cluster_support() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            min_single_exon_support: 3,
            ..ScoringParams::default()
        };

        let reads: Vec<GenePred> = (0..3u64)
            .map(|i| {
                make_record(
                    format!("se{}", i).as_bytes(),
                    100 + i,
                    300 + i,
                    vec![100 + i],
                    vec![300 + i],
                )
            })
            .collect();

        let classified = de_novo(&reads, &counter, &params, None);
        let (keep, _) = split(&classified);

        assert_eq!(keep.len(), 3);
        assert_eq!(
            report_of(&classified, b"se0").reasons_field(),
            "LOW_COMPONENT_SUPPORT/SUPPORTED_SINGLE_EXON_CLUSTER"
        );
    }

    // -----------------------------------------------------------------------
    // Leave-one-out intron support
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_cannot_support_its_own_junctions() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            ..ScoringParams::default()
        };

        // two multi-exon reads with mutually inconsistent chains, plus single-exon reads
        // to make the component large. Counting a read's own junctions would give each of
        // them a support fraction of 1.000
        let mut reads = vec![
            make_record(b"me_a", 1000, 1400, vec![1000, 1300], vec![1100, 1400]),
            make_record(b"me_b", 1000, 1600, vec![1000, 1500], vec![1100, 1600]),
        ];
        for i in 0..8u64 {
            reads.push(make_record(
                format!("se{}", i).as_bytes(),
                1000,
                1600,
                vec![1000],
                vec![1600],
            ));
        }

        let classified = de_novo(&reads, &counter, &params, None);

        for name in [b"me_a".as_slice(), b"me_b".as_slice()] {
            let report = report_of(&classified, name);
            assert_eq!(
                report.intron_support_fraction,
                Some(0.0),
                "{:?} has no independent support",
                String::from_utf8_lossy(name).to_string()
            );
            assert_eq!(report.decision, FinalDecision::Scrap);
            assert!(report.reasons.contains(&DecisionReason::LowIntronSupport));
        }
    }

    #[test]
    fn test_intron_support_counts_only_other_reads() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            ..ScoringParams::default()
        };

        // six reads share (200,300) and (400,500); the variant shares only the first
        let mut reads = dominant_reads(6);
        reads.push(make_record(
            b"variant",
            100,
            650,
            vec![100, 300, 550],
            vec![200, 400, 650],
        ));

        let classified = de_novo(&reads, &counter, &params, None);
        let report = report_of(&classified, b"variant");

        // (200,300) is carried by 6 of the 6 other reads, (400,550) by none
        assert_eq!(report.intron_support_fraction, Some(0.5));
        assert_eq!(report.decision, FinalDecision::Pass);
        assert!(report
            .reasons
            .contains(&DecisionReason::IntronSupportRescue));
    }

    #[test]
    fn test_intron_support_uses_the_junction_tolerance() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 99,
            junction_tolerance: 5,
            ..ScoringParams::default()
        };

        // five reads whose junctions differ by a few bases: within the tolerance that
        // clusters intron chains, so they must also support each other's junctions
        let reads: Vec<GenePred> = (0..5u64)
            .map(|i| {
                make_record(
                    format!("r{}", i).as_bytes(),
                    100,
                    600,
                    vec![100, 300 + i],
                    vec![200 + i, 600],
                )
            })
            .collect();

        let classified = de_novo(&reads, &counter, &params, None);
        let report = report_of(&classified, b"r0");

        assert_eq!(
            report.intron_support_fraction,
            Some(1.0),
            "junctions within junction_tolerance corroborate each other"
        );
        assert_eq!(report.decision, FinalDecision::Pass);
    }

    #[test]
    fn test_exact_intron_support_without_tolerance() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 99,
            junction_tolerance: 0,
            ..ScoringParams::default()
        };

        let reads: Vec<GenePred> = (0..5u64)
            .map(|i| {
                make_record(
                    format!("r{}", i).as_bytes(),
                    100,
                    600,
                    vec![100, 300 + i],
                    vec![200 + i, 600],
                )
            })
            .collect();

        let classified = de_novo(&reads, &counter, &params, None);
        let report = report_of(&classified, b"r0");

        assert_eq!(
            report.intron_support_fraction,
            Some(0.0),
            "with an exact tolerance no read corroborates another"
        );
        assert_eq!(report.decision, FinalDecision::Scrap);
    }

    // -----------------------------------------------------------------------
    // Dominant clusters and splice vetoes
    // -----------------------------------------------------------------------

    #[test]
    fn test_dominant_intron_chain_cluster_passes_unconditionally() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 3,
            ..ScoringParams::default()
        };

        let reads: Vec<GenePred> = dominant_reads(5)
            .into_iter()
            .map(|mut read| {
                read.strand = Some(genepred::Strand::Forward);
                read
            })
            .collect();
        // low scores everywhere: a dominant cluster is not vetoed
        let provider = provider_for(&[(200, 300), (400, 500)], 0.01);

        let classified = de_novo(&reads, &counter, &params, Some(&provider));
        let (keep, orphans) = split(&classified);

        assert_eq!(keep.len(), 5);
        assert!(orphans.is_empty());

        let report = report_of(&classified, b"dom0");
        assert_eq!(report.reasons_field(), "DOMINANT_INTRON_CHAIN_CLUSTER");
        assert_eq!(report.cluster_size, Some(5));
        assert_eq!(
            report.median_splice_score, None,
            "splice scores are not consulted for a dominant cluster"
        );
    }

    #[test]
    fn test_splice_veto_after_intron_support_rescue() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            ..ScoringParams::default()
        };

        let mut reads: Vec<GenePred> = (0..6u64)
            .map(|i| {
                make_stranded_record(
                    format!("dom{}", i).as_bytes(),
                    100 + i * 2,
                    600,
                    vec![100 + i * 2, 300, 500],
                    vec![200, 400, 600],
                    genepred::Strand::Forward,
                )
            })
            .collect();
        reads.push(make_stranded_record(
            b"variant",
            100,
            650,
            vec![100, 300, 550],
            vec![200, 400, 650],
            genepred::Strand::Forward,
        ));

        let provider = provider_for(&[(200, 300), (400, 500), (400, 550)], 0.1);
        let classified = de_novo(&reads, &counter, &params, Some(&provider));
        let report = report_of(&classified, b"variant");

        assert_eq!(report.decision, FinalDecision::Scrap);
        assert_eq!(
            report.reasons_field(),
            "LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_RESCUE/LOW_SPLICE_SCORE"
        );
    }

    #[test]
    fn test_median_splice_score_decides_not_the_weakest_site() {
        use dashmap::DashMap;

        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        // three of four sites are strong and one is outside BigWig coverage: the read is
        // rescued on its median, and its weakest site is reported without gating anything
        let donor = DashMap::new();
        let acceptor = DashMap::new();
        let mut d = HashMap::new();
        let mut a = HashMap::new();
        d.insert(200u64, 0.9f32);
        d.insert(400u64, 0.9f32);
        a.insert(299u64, 0.9f32);
        // acceptor of intron (400,500) is missing
        donor.insert("chr1".to_string(), d);
        acceptor.insert("chr1".to_string(), a);
        let provider =
            SpliceScoreProvider::from_maps(donor, DashMap::new(), acceptor, DashMap::new());

        let read = make_stranded_record(
            b"lone",
            100,
            600,
            vec![100, 300, 500],
            vec![200, 400, 600],
            genepred::Strand::Forward,
        );

        let classified = de_novo(&[read], &counter, &params, Some(&provider));
        let report = report_of(&classified, b"lone");

        assert_eq!(report.median_splice_score, Some(0.9));
        assert_eq!(report.weakest_splice_score, Some(0.0));
        assert_eq!(
            report.decision,
            FinalDecision::Pass,
            "one uncovered splice site must not sink an otherwise supported read"
        );
        assert!(report.reasons.contains(&DecisionReason::SpliceScoreRescue));
    }

    #[test]
    fn test_unstranded_read_has_no_splice_evidence() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let read = make_record(b"lone", 100, 600, vec![100, 300], vec![200, 600]);
        let provider = provider_for(&[(200, 300)], 0.9);

        let classified = de_novo(&[read], &counter, &params, Some(&provider));
        let report = report_of(&classified, b"lone");

        assert_eq!(report.median_splice_score, None);
        assert_eq!(report.weakest_splice_score, None);
        assert_eq!(report.decision, FinalDecision::Scrap);
        assert!(report
            .reasons
            .contains(&DecisionReason::SpliceScoreUnavailable));
    }

    #[test]
    fn test_reverse_strand_uses_reverse_tracks() {
        use dashmap::DashMap;

        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        // on the minus strand the donor of intron (200,300) is its last base
        let donor_rev = DashMap::new();
        let acceptor_rev = DashMap::new();
        let mut d = HashMap::new();
        let mut a = HashMap::new();
        d.insert(299u64, 0.8f32);
        a.insert(200u64, 0.8f32);
        donor_rev.insert("chr1".to_string(), d);
        acceptor_rev.insert("chr1".to_string(), a);
        let provider =
            SpliceScoreProvider::from_maps(DashMap::new(), donor_rev, DashMap::new(), acceptor_rev);

        let read = make_stranded_record(
            b"lone",
            100,
            600,
            vec![100, 300],
            vec![200, 600],
            genepred::Strand::Reverse,
        );

        let classified = de_novo(&[read], &counter, &params, Some(&provider));
        let report = report_of(&classified, b"lone");

        assert_eq!(report.weakest_splice_score, Some(0.8));
        assert_eq!(report.decision, FinalDecision::Pass);
    }

    // -----------------------------------------------------------------------
    // Report integrity
    // -----------------------------------------------------------------------

    #[test]
    fn test_one_row_per_query_record() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        let queries = vec![
            make_record(b"in1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]),
            make_record(b"in2", 150, 250, vec![150], vec![250]),
            make_record(b"oob1", 220, 280, vec![220], vec![280]),
            make_record(b"oob2", 700, 900, vec![700, 850], vec![800, 900]),
        ];
        let expected = queries.len();

        let classified = guided(vec![reference()], queries, &counter, &params, None);

        assert_eq!(classified.len(), expected);
        for read in &classified {
            assert!(matches!(
                read.report.decision,
                FinalDecision::Pass | FinalDecision::Scrap
            ));
            assert!(!read.report.reasons.is_empty());
        }
    }

    #[test]
    fn test_duplicate_query_names_keep_separate_rows() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        let reads: Vec<GenePred> = (0..6u64)
            .map(|i| {
                make_record(
                    b"dup",
                    100 + i * 2,
                    600,
                    vec![100 + i * 2, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect();

        let classified = de_novo(&reads, &counter, &params, None);
        assert_eq!(classified.len(), 6);

        let mut rows: Vec<ReadReport> = classified.into_iter().map(|read| read.report).collect();
        sort_reports(&mut rows);
        assert_eq!(rows.len(), 6);
    }

    #[test]
    fn test_identical_duplicate_records_keep_separate_rows() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        let reads: Vec<GenePred> = (0..6)
            .map(|_| make_record(b"same", 100, 600, vec![100, 300, 500], vec![200, 400, 600]))
            .collect();

        let classified = de_novo(&reads, &counter, &params, None);
        assert_eq!(classified.len(), 6);

        let accumulator = ParallelAccumulator::default();
        accumulator.add(classified);

        assert_eq!(
            accumulator.num_reports(),
            6,
            "the report must not deduplicate records"
        );
        assert_eq!(
            accumulator.keep.len(),
            1,
            "the BED output keeps its existing deduplication"
        );
    }

    #[test]
    fn test_report_decisions_agree_with_bed_outputs() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        let queries = vec![
            make_record(b"in1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]),
            make_record(b"in2", 150, 400, vec![150], vec![400]),
            make_record(b"oob1", 220, 280, vec![220], vec![280]),
        ];

        let classified = guided(vec![reference()], queries, &counter, &params, None);
        let expected_keep: Vec<Vec<u8>> = classified
            .iter()
            .filter(|read| read.report.decision.is_pass())
            .map(|read| read.bed_line.clone())
            .collect();
        let expected_orphans: Vec<Vec<u8>> = classified
            .iter()
            .filter(|read| !read.report.decision.is_pass())
            .map(|read| read.bed_line.clone())
            .collect();

        let accumulator = ParallelAccumulator::default();
        accumulator.add(classified);

        assert_eq!(accumulator.keep.len(), expected_keep.len());
        assert_eq!(accumulator.orphans.len(), expected_orphans.len());
        for line in &expected_keep {
            assert!(accumulator.keep.contains(line));
            assert!(!accumulator.orphans.contains(line));
        }
        for line in &expected_orphans {
            assert!(accumulator.orphans.contains(line));
            assert!(!accumulator.keep.contains(line));
        }
    }

    #[test]
    fn test_deterministic_across_runs() {
        let params = ScoringParams {
            min_cluster_support: 3,
            ..ScoringParams::default()
        };

        let make_reads = || -> Vec<GenePred> {
            let mut reads = dominant_reads(5);
            reads.push(make_record(
                b"variant",
                100,
                650,
                vec![100, 300, 550],
                vec![200, 400, 650],
            ));
            reads
        };

        let counter1 = ParallelCounter::default();
        let first = de_novo(&make_reads(), &counter1, &params, None);
        let counter2 = ParallelCounter::default();
        let second = de_novo(&make_reads(), &counter2, &params, None);

        let mut rows1: Vec<ReadReport> = first.into_iter().map(|read| read.report).collect();
        let mut rows2: Vec<ReadReport> = second.into_iter().map(|read| read.report).collect();
        sort_reports(&mut rows1);
        sort_reports(&mut rows2);

        let rendered =
            |rows: &[ReadReport]| -> Vec<String> { rows.iter().map(|r| r.to_tsv_row()).collect() };
        assert_eq!(rendered(&rows1), rendered(&rows2));
    }

    #[test]
    fn test_bed_line_matches_reported_identity() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let classified = guided(vec![reference()], vec![query], &counter, &params, None);
        let read = &classified[0];

        let line = String::from_utf8(read.bed_line.clone()).unwrap();
        let bed: Vec<&str> = line.split('\t').collect();

        assert_eq!(bed[0].as_bytes(), read.report.chrom.as_slice());
        assert_eq!(bed[1], read.report.start.to_string());
        assert_eq!(bed[2], read.report.end.to_string());
        assert_eq!(bed[3].as_bytes(), read.report.read_id.as_deref().unwrap());
    }

    #[test]
    fn test_self_guided_mode_ignores_references() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        let queries = vec![
            make_record(b"q0", 100, 600, vec![100, 300, 500], vec![200, 400, 600]),
            make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]),
        ];

        let classified = de_novo(&queries, &counter, &params, None);
        let report = report_of(&classified, b"q0");

        assert_eq!(report.mode, RunMode::DeNovo);
        assert_eq!(report.evaluation_path, EvaluationPath::DeNovo);
        assert_eq!(report.overlaps_reference_exon, None);
        assert_eq!(report.reference_count, 0);
    }
}
