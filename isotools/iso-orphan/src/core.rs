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
    scoring::{
        best_junction_score, best_single_exon_overlap, compute_intron_support, has_boundary_match,
        intron_support_fraction, overlaps_any_reference_exon, ScoringParams,
    },
    splice::SpliceScoreProvider,
    utils::*,
};

use dashmap::DashMap;
use genepred::{Bed12, GenePred};
use packbed::Role;
use rayon::prelude::*;

pub type Components = DashMap<String, Vec<Vec<GenePred>>>;

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
/// ```
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
        end_tolerance: args.end_tolerance,
        min_intron_support_frac: args.min_intron_support_frac,
        intron_support_threshold: args.intron_support_threshold,
        min_splice_score: args.min_splice_score,
    };

    // Load optional splice-site scores from BigWig files
    let splice_scores = args.splicing_scores.as_ref().map(|dir| {
        SpliceScoreProvider::from_dir(dir, 0.01).unwrap_or_else(|e| {
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

            let mut modes = std::iter::repeat(Role::Reference)
                .take(inputs.len())
                .collect::<Vec<_>>();
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
/// ```rust, no_run
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

    let mut p_writer = BufWriter::new(File::create(outdir.join(pass)).unwrap());
    let mut o_writer = BufWriter::new(File::create(outdir.join(orphans)).unwrap());

    accumulator.keep.into_iter().for_each(|line| {
        p_writer.write_all(&line).unwrap();
        p_writer.write_all(&b"\n".as_slice()).unwrap();
    });

    accumulator.orphans.into_iter().for_each(|line| {
        o_writer.write_all(&line).unwrap();
        o_writer.write_all(&b"\n".as_slice()).unwrap();
    });
}

/// Processes a component of reads
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
/// ```rust, no_run
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
    components.into_par_iter().for_each(|component| match mode {
        Mode::Guided => {
            let (refs, queries): (Vec<GenePred>, Vec<GenePred>) =
                component.into_iter().partition(|record| {
                    let role = record
                        .get_extra(b"role")
                        .unwrap_or_else(|| panic!("ERROR: Could not get role from record!"))
                        .clone()
                        .into_scalar();

                    role == Some(b"reference".to_vec())
                });

            if !queries.is_empty() {
                log::debug!(
                    "DEBUG: processing guided mode; component of size {}",
                    queries.len()
                );
            } else {
                log::trace!("DEBUG: component of size 0 -> skipping");
                return;
            }

            let (keep, orphans) = guided(refs, queries, counter, params, splice_scores);
            accumulator.add(keep, orphans);
        }
        Mode::DeNovo => {
            let (_refs, queries) = component.into_iter().partition(|record| {
                let role = record
                    .get_extra(b"role")
                    .unwrap_or_else(|| panic!("ERROR: Could not get role from record!"))
                    .clone()
                    .into_scalar();

                role == Some(b"reference".to_vec())
            });

            let (keep, orphans) = self_guided(queries, counter, params, splice_scores);
            accumulator.add(keep, orphans);
        }
    });
}

/// Guided mode: score each query read against individual reference transcripts.
///
/// Reads within reference bounds are scored by junction/overlap agreement.
/// Reads extending beyond reference bounds are redirected to de novo evaluation.
///
/// # Arguments
///
/// * `references` - Reference transcripts in the component
/// * `queries` - Query reads to classify
/// * `counter` - Stats counter
/// * `params` - Scoring parameters
/// * `splice_scores` - Optional splice-site score provider
fn guided(
    references: Vec<GenePred>,
    queries: Vec<GenePred>,
    counter: &ParallelCounter,
    params: &ScoringParams,
    splice_scores: Option<&SpliceScoreProvider>,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut keep = Vec::new();
    let mut orphans = Vec::new();

    // INFO: single component reads -> no reference
    if references.is_empty() {
        // CASE: single-exon-component reads [PARTIALLY SOLVED]
        // INFO: the following case is presented:
        //
        // [SOLVED: single-component reads will be just discarded]
        // ref:   XX-----------XXX--XX---XXXX
        // read1:    xxXXXXxx  |||  ||   ||||
        // read2: XX-----------XXX--XX---XXXX
        //
        // Here, read1 is a single-exon component read that can be easily
        // discarded because it represents a common case of background noise.
        // For a clearer illustration go to mm39 chr8:71,358,478-71,360,769,
        // or mm39 chr8:72,647,487-72,650,648.
        //
        // or its variant [UNSOLVED]
        //
        // ref:   XX-----------XXX--XX---XXXX
        // read1:    xxXXXXxx
        // read2:    xxxxXXXxx
        // read3:     xXXXXXXxx
        // read4: XX-----------XXX--XX---XXXX
        //
        // Here, read1-3 are single-reads in a component that does not necessarilly
        // represent a single transcript. Here, we should try CDS exact matches. If
        // the reads do not match exactly, we should discard them.
        //
        // For a clearer illustration go to HLpteAle1A HAP1_SUPER_1:112,985,006-112,991,219
        // or mm39 chr3:118,192,628-118,376,044 or mm39 chr3:141,166,307-141,389,231
        return do_self_guided_check(&queries, counter, params, splice_scores);
    }

    let total_queries = queries.len();

    if queries.is_empty() {
        let projections: Vec<Option<&[u8]>> = references.iter().map(|p| p.name()).collect();
        log::trace!(
            "DEBUG: reference projection without reads -> {:?}",
            projections
        );
        return (keep, orphans);
    }

    // -----------------------------------------------------------------------
    // Partition queries by exon overlap with references
    // -----------------------------------------------------------------------
    //
    // CASE: overlapping vs non-overlapping reads
    //
    // ref:    xxxxXXXX----XXXXX----XXXxxx
    //
    // read1:  xxxxXXXX----XXXXX              → overlaps ref exons → guided
    // read2:        xxxXXXXxxx               → overlaps ref exon  → guided
    // read3:              xxxxxxxxxxxx       → sits in ref intron  → de novo
    // read4: xxxx                            → no exon overlap     → de novo
    //
    // Reads whose exons overlap at least one reference exon are scored in
    // guided mode. Reads with zero exon overlap (intronic, intergenic,
    // chimeric with no reference support) are redirected to de novo.

    let (in_bounds, oob): (Vec<GenePred>, Vec<GenePred>) = queries
        .into_iter()
        .partition(|read| overlaps_any_reference_exon(read, &references));

    log::debug!(
        "DEBUG: {} in bound reads vs {} out of bound in component of size {}",
        in_bounds.len(),
        oob.len(),
        total_queries
    );

    // Score in-bounds reads against references
    let has_multi_exon_refs = references.iter().any(|r| !r.introns().is_empty());
    log::debug!("DEBUG: has_multi_exon_refs={}", has_multi_exon_refs);

    for read in &in_bounds {
        let is_single_exon_read = read.introns().is_empty();
        let line = read.to_bed::<Bed12>();

        if is_single_exon_read {
            // Single-exon reads: score by reciprocal overlap with best reference exon.
            //
            // CASE: single-exon read(s) with single-exon reference refs
            // INFO: the following case is presented:
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
            //
            // CASE: single-exon fuzzy overlap with reference CDS
            // INFO: the following case is presented:
            //
            // ref:      xxxxxxXXXXXXX----------
            //               ^^||||^^^
            // read1: xxxxxxxXXXXXXXxxxx
            // read2:       xxxXXXXXXX----------
            // read3:     xxxxxXXXXXXX----------
            //
            // For a clearer illustration go to mm39 chr3:65,014,198-65,019,415
            // or mm39 chr3:117,366,942-117,374,421
            let overlap = best_single_exon_overlap(read, &references);

            if overlap >= params.min_overlap_frac {
                log::debug!(
                    "DEBUG: read {:?} single-exon overlap={:.3} >= {:.3} -> keep!",
                    read.name(),
                    overlap,
                    params.min_overlap_frac,
                );
                counter.inc_guided_se_keep();
                keep.push(line);
            } else {
                log::debug!(
                    "DEBUG: read {:?} single-exon overlap={:.3} < {:.3} -> orphan!",
                    read.name(),
                    overlap,
                    params.min_overlap_frac,
                );
                counter.inc_guided_se_orphan();
                orphans.push(line);
            }
        } else {
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
            if has_multi_exon_refs {
                if let Some(best) =
                    best_junction_score(read, &references, params.junction_tolerance)
                {
                    if best.passes(params.min_junction_frac) {
                        log::debug!(
                            "DEBUG: read {:?} multi-exon junction_frac={:.3} ({}/{} matches) -> keep!",
                            read.name(),
                            best.query_junction_frac,
                            best.junction_matches,
                            best.query_junction_count,
                        );
                        counter.inc_guided_me_junction_keep();
                        keep.push(line);
                    } else if has_boundary_match(read, &references, params.end_tolerance) {
                        log::debug!(
                            "DEBUG: read {:?} multi-exon junction_frac={:.3} low, but has boundary match -> boundary keep!",
                            read.name(),
                            best.query_junction_frac,
                        );
                        counter.inc_guided_me_boundary_keep();
                        keep.push(line);
                    } else {
                        log::debug!(
                            "DEBUG: read {:?} multi-exon junction_frac={:.3}, no boundary match -> orphan!",
                            read.name(),
                            best.query_junction_frac,
                        );
                        counter.inc_guided_me_orphan();
                        orphans.push(line);
                    }
                } else {
                    counter.inc_guided_me_orphan();
                    orphans.push(line);
                }
            } else {
                // All references are single-exon; use boundary matching or test splice sites if
                // scores provide
                log::debug!("DEBUG: multi-exon read vs single-exon refs, boundary match or splice site score");

                if has_boundary_match(read, &references, params.junction_tolerance) {
                    log::debug!(
                        "DEBUG: read {:?} multi-exon vs single-exon refs, boundary match -> keep!",
                        read.name(),
                    );
                    counter.inc_guided_me_boundary_keep();
                    keep.push(line);
                } else if splice_scores
                    .and_then(|scores| scores.median_splice_score(read))
                    .map(|median| {
                        log::debug!("DEBUG: splice score median={}", median);
                        median >= params.min_splice_score
                    })
                    .unwrap_or(false)
                // No scores or no strand → pass
                {
                    log::debug!(
                            "DEBUG: read {:?} multi-exon vs single-exon refs, no boundary match but splice site score above threshold -> keep!", read.name()
                        );

                    counter.inc_guided_me_splice_keep();
                    keep.push(line);
                } else {
                    log::debug!(
                        "DEBUG: read {:?} multi-exon vs single-exon refs, no boundary match nor splice site score -> orphan!",
                        read.name(),
                    );
                    counter.inc_guided_me_orphan();
                    orphans.push(line);
                }
            }
        }
    }

    // Redirect out-of-bounds reads to de novo evaluation
    if !oob.is_empty() {
        log::debug!(
            "DEBUG: {} reads extend beyond reference bounds -> de novo",
            oob.len()
        );
        counter.inc_guided_oob_denovo(oob.len() as u32);
        let (oob_keep, oob_orphans) = do_self_guided_check(&oob, counter, params, splice_scores);
        keep.extend(oob_keep);
        orphans.extend(oob_orphans);
    }

    debug_assert_eq!(
        keep.len() + orphans.len(),
        total_queries,
        "Partition violation in guided mode: keep({}) + orphan({}) != input({})",
        keep.len(),
        orphans.len(),
        total_queries
    );

    (keep, orphans)
}

/// Self-guided mode: score each query read against individual reference transcripts.
///
/// Reads within reference bounds are scored by junction/overlap agreement.
/// Reads extending beyond reference bounds are redirected to de novo evaluation.
///
/// # Arguments
///
/// * `queries` - A vector of GenePred structs
/// * `counter` - A reference to a ParallelCounter struct
/// * `params` - A reference to a ScoringParams struct
/// * `splice_scores` - A reference to a SpliceScoreProvider struct
///
/// # Example
///
/// ```rust, no_run
/// let queries = vec![GenePred::default()];
/// let counter = ParallelCounter::default();
/// let params = ScoringParams::default();
/// let splice_scores = None;
///
/// self_guided(queries, &counter, &params, splice_scores);
/// ```
fn self_guided(
    queries: Vec<GenePred>,
    counter: &ParallelCounter,
    params: &ScoringParams,
    splice_scores: Option<&SpliceScoreProvider>,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    do_self_guided_check(&queries, counter, params, splice_scores)
}

/// Self-guided mode executor: classifies reads within a component using de novo clustering
/// with per-intron support evaluation and optional splice-site score filtering.
///
/// Decision flow for multi-exon reads:
/// 1. Dominant intron-chain cluster (≥ min_cluster_support) → KEEP
/// 2. Per-intron support ≥ threshold → tentative KEEP
///    a. If splice scores available AND median < min_splice_score → ORPHAN
///    b. Else → KEEP
/// 3. Insufficient intron support → ORPHAN
fn do_self_guided_check(
    reads: &[GenePred],
    counter: &ParallelCounter,
    params: &ScoringParams,
    splice_scores: Option<&SpliceScoreProvider>,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut keep = Vec::new();
    let mut orphans = Vec::new();

    // CASE: single-component single-exon
    // INFO: the following case is presented:
    //
    // [SOLVED: discarding single-component reads]
    // ref:             XXX--XX---XXXX
    // read1: xxXXXXxx  |||  ||   ||||
    // read2:           XXX--XX---XXXX
    //
    // or its variants:
    //
    // [SOLVED: discarding single-component reads]
    // ref:   XX-----------XXX--XX---XXXX
    // read1:    xxXXXXxx  |||  ||   ||||
    // read2: XX-----------XXX--XX---XXXX
    //
    // Here, read1 is a single-exon component read that can be easily
    // discarded because it represents a common case of background noise.
    // For a clearer illustration go to mm39 chr8:71,358,478-71,360,769,
    // or mm39 chr8:72,647,487-72,650,648.
    //
    // ref:                  XXX--XX---XXXX
    // read1: xxXXXX--XXXxx  |||  ||   ||||
    // read2:                XXX--XX---XXXX
    //
    // The above case would illustrate a component where a non-single-exon
    // component read could repesent real transcription but not at a significant
    // level.
    //
    // For a clearer illustration go to mm39 chr8:71,358,478-71,360,769
    if reads.len() <= 1 {
        log::trace!("DEBUG: single-component read with no reference refs -> orphan!");
        counter.inc_denovo_single_read();

        for read in reads {
            orphans.push(read.to_bed::<Bed12>());
        }

        return (keep, orphans);
    }

    // CASE: multi-exon components
    //
    // [PARTIALLY SOLVED: min_read_num_denovo]
    // ref:   XX-----------XXX--XX---XXXX
    // read1:    xxXXXXxx
    // read2:    xxxxXXXxx
    // read3: XX-----------XXX--XX---XXXX
    //
    // For a clearer illustration go to mm39 chr3:97,596,920-97,624,824
    // or mm39 chr3:103,646,891-103,652,710
    if reads.len() < params.min_cluster_support {
        log::debug!(
            "DEBUG: component with {} reads < min_cluster_support={} -> orphan!",
            reads.len(),
            params.min_cluster_support,
        );
        counter.inc_denovo_small_component();

        for read in reads {
            orphans.push(read.to_bed::<Bed12>());
        }

        return (keep, orphans);
    }

    // -----------------------------------------------------------------------
    // Partition reads into multi-exon and single-exon groups
    // -----------------------------------------------------------------------
    //
    // CASE: single-exon reads following same/diff exonic pattern
    //
    // read1:       xxxXXXXXXXxxx
    // read2:      xxxxxxxxxxXXXXXxxxx
    //
    // For a clearer illustration go to HLpteAle1A HAP1_SUPER_1:112,960,536-112,964,453
    //
    // CASE: single-exon read partially overlap CDS
    //
    // read1: xxxxxxxxxXXXXxxx
    // read2: -------xxxxXXXXX------XXX--
    // read3: ---------xxXXXXX------XXX--
    // read4: ---------xxxxxxxxxXX--XXX--
    // read5: ----------xXXXXX-----------
    //                   |||||      |||
    //
    // For a clearer illustration go to mm39 chr8:72,042,766-72,051,539
    // or mm39 chr8:73,174,008-73,179,856
    let (multi_exon_indexed, single_exon_indexed): (Vec<_>, Vec<_>) = reads
        .iter()
        .enumerate()
        .partition(|(_, r)| !r.introns().is_empty());

    // -----------------------------------------------------------------------
    // Multi-exon reads: cluster by intron chain + per-intron support + splice veto
    // -----------------------------------------------------------------------
    let me_reads: Vec<&GenePred> = multi_exon_indexed.iter().map(|(_, r)| *r).collect();
    let me_indices: Vec<usize> = multi_exon_indexed.iter().map(|(i, _)| *i).collect();

    // Build component-wide intron support map
    let intron_support = compute_intron_support(&me_reads);
    let total_me = me_reads.len();

    // Cluster by intron chain
    let me_clusters = cluster_multi_exon(&me_reads, params.junction_tolerance);

    for (_cid, members) in &me_clusters {
        let is_dominant = members.len() >= params.min_cluster_support;

        for &local_idx in members {
            let global_idx = me_indices[local_idx];
            let read = &reads[global_idx];
            let line = read.to_bed::<Bed12>();

            if is_dominant {
                // Strong evidence: dominant intron-chain cluster → unconditional keep
                keep.push(line);
                counter.inc_denovo_me_cluster_keep();
            } else {
                // Check per-intron support across the component
                let support_frac = intron_support_fraction(
                    read,
                    &intron_support,
                    total_me,
                    params.intron_support_threshold,
                );

                if support_frac >= params.min_intron_support_frac {
                    // Intron support passes; apply splice-score veto if available
                    let splice_pass = splice_scores
                        .and_then(|scores| scores.median_splice_score(read))
                        .map(|median| median >= params.min_splice_score)
                        .unwrap_or(true); // No scores or no strand → pass

                    if splice_pass {
                        keep.push(line);
                        counter.inc_denovo_me_intron_keep();
                    } else {
                        orphans.push(line);
                        counter.inc_denovo_me_splice_orphan();
                    }
                } else {
                    orphans.push(line);
                    counter.inc_denovo_me_cluster_orphan();
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Single-exon reads: cluster by reciprocal overlap
    // -----------------------------------------------------------------------
    let se_reads: Vec<&GenePred> = single_exon_indexed.iter().map(|(_, r)| *r).collect();
    let se_indices: Vec<usize> = single_exon_indexed.iter().map(|(i, _)| *i).collect();

    let se_clusters = cluster_single_exon(&se_reads, params.min_overlap_frac);
    for (_cid, members) in &se_clusters {
        if members.len() >= params.min_cluster_support {
            for &local_idx in members {
                let global_idx = se_indices[local_idx];
                keep.push(reads[global_idx].to_bed::<Bed12>());
                counter.inc_denovo_se_cluster_keep();
            }
        } else {
            for &local_idx in members {
                let global_idx = se_indices[local_idx];
                orphans.push(reads[global_idx].to_bed::<Bed12>());
                counter.inc_denovo_se_cluster_orphan();
            }
        }
    }

    debug_assert_eq!(
        keep.len() + orphans.len(),
        reads.len(),
        "Partition violation in de novo mode: keep({}) + orphan({}) != input({})",
        keep.len(),
        orphans.len(),
        reads.len()
    );

    (keep, orphans)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringParams;
    use crate::utils::ParallelCounter;
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

    // -----------------------------------------------------------------------
    // Guided: normal in-bounds case
    // -----------------------------------------------------------------------

    #[test]
    fn test_guided_normal_in_bounds() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let ref1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let queries: Vec<GenePred> = (0..5)
            .map(|i| {
                make_record(
                    format!("q{}", i).as_bytes(),
                    100 + i * 10,
                    600,
                    vec![100 + i * 10, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect();

        let n = queries.len();
        let (keep, orphans) = guided(vec![ref1], queries, &counter, &params, None);

        assert_eq!(keep.len() + orphans.len(), n, "Partition violation");
        assert!(
            keep.len() >= 4,
            "Most in-bounds reads with matching junctions should be kept"
        );
    }

    // -----------------------------------------------------------------------
    // Guided: bridging / out-of-bounds redirect
    // -----------------------------------------------------------------------

    #[test]
    fn test_guided_oob_redirected_to_denovo() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        // Reference: exons at [100,200), [300,400), [500,600)
        let ref1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        // 2 queries overlapping ref exons → guided
        // 2 queries with NO exon overlap → de novo
        let queries = vec![
            // Overlaps ref exon [300,400)
            make_record(b"in1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]),
            // Overlaps ref exon [100,200) partially
            make_record(b"in2", 150, 250, vec![150], vec![250]),
            // Sits in ref intron [200,300) — no exon overlap
            make_record(b"oob1", 220, 280, vec![220], vec![280]),
            // Entirely beyond ref — no exon overlap
            make_record(b"oob2", 700, 900, vec![700, 850], vec![800, 900]),
        ];

        let n = queries.len();
        let (keep, orphans) = guided(vec![ref1], queries, &counter, &params, None);

        assert_eq!(keep.len() + orphans.len(), n, "Partition violation");
        let oob_count = counter
            .guided_oob_denovo
            .load(std::sync::atomic::Ordering::Relaxed);
        assert_eq!(oob_count, 2, "2 reads have no exon overlap with reference");
    }

    // -----------------------------------------------------------------------
    // Guided: subset redirected to de novo (no refs)
    // -----------------------------------------------------------------------

    #[test]
    fn test_guided_no_refs_falls_to_denovo() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        let queries: Vec<GenePred> = (0..5)
            .map(|i| {
                make_record(
                    format!("q{}", i).as_bytes(),
                    100,
                    600,
                    vec![100, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect();

        let n = queries.len();
        let (keep, orphans) = guided(vec![], queries, &counter, &params, None);

        assert_eq!(keep.len() + orphans.len(), n, "Partition violation");
    }

    // -----------------------------------------------------------------------
    // De novo: exact intron-chain cluster keep
    // -----------------------------------------------------------------------

    #[test]
    fn test_denovo_exact_intron_chain_keep() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 3,
            ..ScoringParams::default()
        };

        // 5 reads with exact same intron chain → dominant cluster → all kept
        let reads: Vec<GenePred> = (0..5u64)
            .map(|i| {
                make_record(
                    format!("r{}", i).as_bytes(),
                    100 + i * 3,
                    600,
                    vec![100 + i * 3, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect();

        let (keep, orphans) = do_self_guided_check(&reads, &counter, &params, None);

        assert_eq!(
            keep.len(),
            5,
            "All reads in dominant cluster should be kept"
        );
        assert_eq!(orphans.len(), 0);
    }

    // -----------------------------------------------------------------------
    // De novo: 9/10 supported introns → keep via intron support
    // -----------------------------------------------------------------------

    #[test]
    fn test_denovo_9_of_10_introns_supported_keep() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            min_intron_support_frac: 0.5,
            intron_support_threshold: 0.5,
            ..ScoringParams::default()
        };

        // 6 reads form the dominant cluster with chain [(200,300), (400,500)]
        let mut reads: Vec<GenePred> = (0..6u64)
            .map(|i| {
                make_record(
                    format!("dom{}", i).as_bytes(),
                    100 + i * 2,
                    600,
                    vec![100 + i * 2, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect();

        // 1 variant read: shares intron (200,300) but has a different second intron (400,550)
        // This read's cluster has only 1 member (< min_cluster_support=5).
        // But intron (200,300) is supported (7/7 reads have it).
        // Intron (400,550) is NOT supported (1/7 reads).
        // So 1/2 = 0.5 of this read's introns are supported → passes 0.5 threshold
        reads.push(make_record(
            b"variant",
            100,
            650,
            vec![100, 300, 550],
            vec![200, 400, 650],
        ));

        let n = reads.len();
        let (keep, orphans) = do_self_guided_check(&reads, &counter, &params, None);

        assert_eq!(keep.len() + orphans.len(), n, "Partition violation");
        // The variant read should be kept by intron support
        assert_eq!(keep.len(), 7, "6 dominant + 1 intron-support-kept");
    }

    // -----------------------------------------------------------------------
    // De novo: weak intron support → discard
    // -----------------------------------------------------------------------

    #[test]
    fn test_denovo_weak_intron_support_discard() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            min_intron_support_frac: 0.5,
            intron_support_threshold: 0.5,
            ..ScoringParams::default()
        };

        // 6 reads with chain [(200,300), (400,500)]
        let mut reads: Vec<GenePred> = (0..6u64)
            .map(|i| {
                make_record(
                    format!("dom{}", i).as_bytes(),
                    100 + i * 2,
                    600,
                    vec![100 + i * 2, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect();

        // Outlier: completely different intron chain (no introns supported)
        reads.push(make_record(
            b"outlier",
            700,
            1200,
            vec![700, 900, 1100],
            vec![800, 1000, 1200],
        ));

        let n = reads.len();
        let (keep, orphans) = do_self_guided_check(&reads, &counter, &params, None);

        assert_eq!(keep.len() + orphans.len(), n, "Partition violation");
        assert_eq!(keep.len(), 6, "Only dominant cluster kept");
        assert_eq!(orphans.len(), 1, "Outlier orphaned");
    }

    // -----------------------------------------------------------------------
    // Splice score: keep when above threshold
    // -----------------------------------------------------------------------

    #[test]
    fn test_denovo_splice_score_keep() {
        use crate::splice::SpliceScoreProvider;
        use dashmap::DashMap;

        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            min_intron_support_frac: 0.5,
            intron_support_threshold: 0.5,
            min_splice_score: 0.5,
            ..ScoringParams::default()
        };

        // Build splice scores: high scores at intron boundaries
        let donor_fwd = DashMap::new();
        let mut d = HashMap::new();
        d.insert(200u64, 0.9f32);
        d.insert(400u64, 0.8f32);
        d.insert(450u64, 0.85f32); // for variant read
        donor_fwd.insert("chr1".to_string(), d);

        let acceptor_fwd = DashMap::new();
        let mut a = HashMap::new();
        a.insert(300u64, 0.85f32);
        a.insert(500u64, 0.9f32);
        a.insert(550u64, 0.8f32);
        acceptor_fwd.insert("chr1".to_string(), a);

        let provider =
            SpliceScoreProvider::from_maps(donor_fwd, DashMap::new(), acceptor_fwd, DashMap::new());

        // 6 dominant reads + 1 variant with good splice scores
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

        // Variant: shares (200,300) with dominant, novel (450,550)
        reads.push(make_stranded_record(
            b"variant",
            100,
            650,
            vec![100, 300, 550],
            vec![200, 450, 650],
            genepred::Strand::Forward,
        ));

        let n = reads.len();
        let (keep, orphans) = do_self_guided_check(&reads, &counter, &params, Some(&provider));

        assert_eq!(keep.len() + orphans.len(), n, "Partition violation");
        // Variant has good splice scores → kept
        assert_eq!(keep.len(), 7);
    }

    // -----------------------------------------------------------------------
    // Splice score: discard when below threshold
    // -----------------------------------------------------------------------

    #[test]
    fn test_denovo_splice_score_discard() {
        use crate::splice::SpliceScoreProvider;
        use dashmap::DashMap;

        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            min_intron_support_frac: 0.5,
            intron_support_threshold: 0.3,
            min_splice_score: 0.5,
            ..ScoringParams::default()
        };

        // Build splice scores: LOW scores everywhere
        let donor_fwd = DashMap::new();
        let mut d = HashMap::new();
        d.insert(200u64, 0.1f32);
        d.insert(400u64, 0.1f32);
        d.insert(450u64, 0.1f32);
        donor_fwd.insert("chr1".to_string(), d);

        let acceptor_fwd = DashMap::new();
        let mut a = HashMap::new();
        a.insert(300u64, 0.1f32);
        a.insert(500u64, 0.1f32);
        a.insert(550u64, 0.1f32);
        acceptor_fwd.insert("chr1".to_string(), a);

        let provider =
            SpliceScoreProvider::from_maps(donor_fwd, DashMap::new(), acceptor_fwd, DashMap::new());

        // 6 dominant reads (kept regardless of splice scores) + 1 variant
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

        // Variant: shares (200,300), novel intron — passes intron support but LOW splice scores
        reads.push(make_stranded_record(
            b"variant",
            100,
            650,
            vec![100, 300, 550],
            vec![200, 450, 650],
            genepred::Strand::Forward,
        ));

        let n = reads.len();
        let (keep, orphans) = do_self_guided_check(&reads, &counter, &params, Some(&provider));

        assert_eq!(keep.len() + orphans.len(), n, "Partition violation");
        // Dominant 6 kept, variant orphaned by splice score veto
        assert_eq!(keep.len(), 6);
        assert_eq!(orphans.len(), 1);
        let splice_orphan_count = counter
            .denovo_me_splice_orphan
            .load(std::sync::atomic::Ordering::Relaxed);
        assert_eq!(splice_orphan_count, 1);
    }

    // -----------------------------------------------------------------------
    // Partition correctness
    // -----------------------------------------------------------------------

    #[test]
    fn test_partition_correctness_guided_single_exon() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let ref1 = make_record(b"r1", 100, 200, vec![100], vec![200]);
        let queries = vec![
            make_record(b"q1", 100, 200, vec![100], vec![200]),
            make_record(b"q2", 500, 600, vec![500], vec![600]),
            make_record(b"q3", 110, 190, vec![110], vec![190]),
        ];

        let n = queries.len();
        let (keep, orphans) = guided(vec![ref1], queries, &counter, &params, None);

        assert_eq!(keep.len() + orphans.len(), n);
        assert!(keep.len() >= 1, "Exact match should be kept");
    }

    #[test]
    fn test_denovo_partition_mixed() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 2,
            ..ScoringParams::default()
        };

        let mut reads: Vec<GenePred> = Vec::new();
        for i in 0..3u64 {
            reads.push(make_record(
                format!("me{}", i).as_bytes(),
                100 + i * 5,
                600,
                vec![100 + i * 5, 300, 500],
                vec![200, 400, 600],
            ));
        }
        for i in 0..3u64 {
            reads.push(make_record(
                format!("se{}", i).as_bytes(),
                100 + i * 10,
                300 + i * 10,
                vec![100 + i * 10],
                vec![300 + i * 10],
            ));
        }

        let n = reads.len();
        let (keep, orphans) = do_self_guided_check(&reads, &counter, &params, None);

        assert_eq!(keep.len() + orphans.len(), n, "Partition violation");
    }

    #[test]
    fn test_keep_orphan_disjoint() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let queries: Vec<GenePred> = (0..10u64)
            .map(|i| {
                if i < 5 {
                    make_record(
                        format!("q{}", i).as_bytes(),
                        100 + i * 3,
                        600,
                        vec![100 + i * 3, 300, 500],
                        vec![200, 400, 600],
                    )
                } else {
                    make_record(
                        format!("q{}", i).as_bytes(),
                        700 + i * 10,
                        1200 + i * 10,
                        vec![700 + i * 10, 900 + i * 10, 1100 + i * 10],
                        vec![800 + i * 10, 1000 + i * 10, 1200 + i * 10],
                    )
                }
            })
            .collect();

        let n = queries.len();
        let (keep, orphans) = guided(vec![reference], queries, &counter, &params, None);

        assert_eq!(keep.len() + orphans.len(), n, "All reads classified");
        for k in &keep {
            assert!(!orphans.contains(k), "Read appears in both keep and orphan");
        }
    }

    // -----------------------------------------------------------------------
    // Deterministic behavior
    // -----------------------------------------------------------------------

    #[test]
    fn test_deterministic_across_runs() {
        let params = ScoringParams {
            min_cluster_support: 3,
            ..ScoringParams::default()
        };

        let make_reads = || -> Vec<GenePred> {
            let mut reads = Vec::new();
            for i in 0..5u64 {
                reads.push(make_record(
                    format!("dom{}", i).as_bytes(),
                    100 + i * 2,
                    600,
                    vec![100 + i * 2, 300, 500],
                    vec![200, 400, 600],
                ));
            }
            reads.push(make_record(
                b"variant",
                100,
                650,
                vec![100, 300, 550],
                vec![200, 400, 650],
            ));
            reads
        };

        // Run twice and verify same result
        let counter1 = ParallelCounter::default();
        let (keep1, orphans1) = do_self_guided_check(&make_reads(), &counter1, &params, None);

        let counter2 = ParallelCounter::default();
        let (keep2, orphans2) = do_self_guided_check(&make_reads(), &counter2, &params, None);

        assert_eq!(keep1.len(), keep2.len());
        assert_eq!(orphans1.len(), orphans2.len());
    }

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn test_denovo_single_read_orphaned() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let reads = vec![make_record(b"lone", 100, 200, vec![100], vec![200])];
        let (keep, orphans) = do_self_guided_check(&reads, &counter, &params, None);

        assert!(keep.is_empty());
        assert_eq!(orphans.len(), 1);
    }

    #[test]
    fn test_denovo_small_component_orphaned() {
        let counter = ParallelCounter::default();
        let params = ScoringParams {
            min_cluster_support: 5,
            ..ScoringParams::default()
        };

        let reads: Vec<GenePred> = (0..3u64)
            .map(|i| {
                make_record(
                    format!("r{}", i).as_bytes(),
                    100 + i,
                    200 + i,
                    vec![100 + i],
                    vec![200 + i],
                )
            })
            .collect();

        let (keep, orphans) = do_self_guided_check(&reads, &counter, &params, None);

        assert!(keep.is_empty());
        assert_eq!(orphans.len(), 3);
    }

    #[test]
    fn test_multi_isoform_reference_coherence() {
        let counter = ParallelCounter::default();
        let params = ScoringParams::default();

        let ref1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let ref2 = make_record(b"r2", 100, 800, vec![100, 350, 700], vec![300, 600, 800]);
        let query = make_record(b"q1", 100, 800, vec![100, 350, 700], vec![300, 600, 800]);

        let (keep, orphans) = guided(vec![ref1, ref2], vec![query], &counter, &params, None);

        assert_eq!(keep.len(), 1);
        assert_eq!(orphans.len(), 0);
    }
}
