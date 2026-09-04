// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for intron classification
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module provides a comprehensive sub-pipeline for classifying
//! introns within genomic sequences. It different data sources to
//! categorize introns based on their predicted splicing potential
//! and structural characteristics.
//!
//! In essence, this module identifies and characterizes introns from
//! input long-read sequencing data. It performs data integration,
//! collecting splice site prediction scores (from tools like SpliceAI
//! and MaxEntScan), analyzing genomic sequence context, and detecting
//! specific sequence patterns such as RT repeats and NULLG motifs. Through
//! a parallel processing approach, each intron is evaluated to determine
//! its "support type", indicating whether it is likely to be a genuine
//! spliced intron, an RT-driven event, or an unclear case requiring
//! further investigation. The final output is a detailed, classified list
//! of introns, enabling deeper insights into alternative splicing and
//! RNULL processing.

use std::fmt::Debug;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::Result;
use dashmap::DashSet;
use genepred::{GenePred, Strand};
use hashbrown::hash_map::Entry;
use hashbrown::{HashMap, HashSet};
use packbed::{pack, Role};
use rayon::prelude::*;
use rust_lapper::Lapper;

use crate::cli::IntronArgs as Args;
use crate::utils::*;

const WINDOW_SIZE: usize = 8;
const RT_REPEAT: usize = 12;
const MISMATCHES: u32 = 1;
const MAX_ENT_DONOR_MIN_SIZE: usize = 9;
const WIGGLE_SWITCH: [usize; 2] = [2, 4];
const NULLG_PATTERNS: [&[u8]; 3] = [b"CAG", b"TAG", b"AAG"];

type ScanScores = Option<(SpliceScoreMap, SpliceScoreMap)>;

const INTRON_SIZE_THRESHOLD: u64 = 50000; // 50kb
const MINIMUM_SPLICEAI_SS_SIGNAL: f32 = 0.001;

/// Classifies introns based on various criteria,
/// including splice site scores, genomic context,
/// and repeat sequences.
///
/// This function orchestrates the entire intron classification pipeline. It loads data,
/// parallelizes the processing by chromosome, and writes the final classified introns to a file.
///
/// # Arguments
///
/// * `args`: An `Args` struct containing all necessary configuration and input paths.
///
/// # Returns
///
/// * A `Result<PathBuf>` which is the path to the output file if the classification is successful.
///
/// # Note
///
/// The final output has the following format, since the header is not included:
/// # Output columns with descriptions:
/// # chrom: Chromosome identifier
/// # start: Start position of the intron (1-based)
/// # end: End position of the intron (1-based)
/// # strand: Genomic strand (+ or -)
/// # seen: Frequency of how many reads contain this intron
/// # spanned: Frequency of how many reads span this intron
/// # splice_ai_donor: SpliceAI score for the donor site (probability of being a real donor)
/// # splice_ai_acceptor: SpliceAI score for the acceptor site (probability of being a real acceptor)
/// # max_ent_donor: MaxEntScan score for the donor site (strength prediction)
/// # max_ent_acceptor: MaxEntScan score for the acceptor site (strength prediction)
/// # donor_sequence: Nucleotide sequence around the donor site
/// # acceptor_sequence: Nucleotide sequence around the acceptor site
/// # donor_context: MaxEntScan 9-mer donor context sequence
/// # acceptor_context: MaxEntScan 23-mer acceptor context sequence
/// # intron_position: Classification of the intron's position according to a reference set
/// # is_reference_supported: Boolean indicating if the intron is supported by a reference set
/// # is_in_frame: Boolean indicating if the intron maintains the reading frame
/// # donor_rt_context: RT-switch context sequence for the donor site
/// # acceptor_rt_context: RT-switch context sequence for the acceptor site
/// # is_rt_intron: Boolean indicating if the intron is an RT-switch intron
/// # is_nag_intron: Boolean indicating if the intron is a reference-nag intron
/// # support: Classification of the intron's support type
///
/// # Example
///
/// ```rust, ignore
/// // Assume 'args' is a properly constructed Args struct
/// let output_path = classify_introns(args).expect("Failed to classify introns");
/// println!("Intron classification complete. Output at: {:?}", output_path);
/// ```
pub fn classify_introns(args: Args) -> Result<()> {
    log::info!("INFO: Classifying introns...");

    let spliceosome = if let Some(iic) = args.iic {
        log::info!("INFO: Parsing intronIC...");
        read_iic(iic)
    } else {
        log::warn!("WARN: No intronIC provided, using empty spliceosome!");
        HashMap::new()
    };

    let isoseqs = if let Some(mut reference) = args.reference {
        let mut modes = std::iter::repeat(Role::Reference)
            .take(reference.len())
            .collect::<Vec<_>>();
        modes.extend(vec![Role::Query]);
        reference.extend(vec![args.input]);

        packbed::pack(reference, modes, packbed::OverlapType::Exon)?
    } else {
        pack(
            vec![args.input],
            vec![Role::Query],
            packbed::OverlapType::Exon,
        )?
    };

    let mut chrs = Vec::with_capacity(isoseqs.len());
    isoseqs.iter().for_each(|bucket| {
        log::debug!(
            "DEBUG: Bucket key: {}, value size: {}",
            bucket.key(),
            bucket.value().len()
        );

        if !bucket.value().is_empty() {
            let bind = bucket.key().clone();
            let chr = bind
                .split(':')
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: Could not get chromosome from key -> {:?}",
                        bucket.key()
                    )
                })
                .to_string();

            chrs.push(chr);
        }
    });

    log::info!(
        "INFO: Collecting data for {} chromosome keys...",
        chrs.len()
    );

    let (splice_plus, splice_minus) = get_splice_scores(args.spliceai, chrs);

    let genome = if let Some(sequence) = args.sequence {
        get_sequences(sequence)
    } else {
        HashMap::new()
    };

    // if args.scan is true, we need to load both databases
    let scan_scores = if args.scan { load_scan_scores() } else { None };
    let repeats = if let Some(path) = args.repeats {
        load_repeats(path)
    } else {
        log::warn!("WARN: No repeats provided, skipping repeats overlapping...");
        None
    };

    let accumulator = ParallelAccumulator::default();

    isoseqs.into_par_iter().for_each(|bucket| {
        log::debug!("INFO: Processing bucket key -> {:?}", bucket.0);

        let chr = bucket.0.split(':').next().unwrap_or_else(|| {
            panic!(
                "ERROR: Could not get chromosome from key -> {:?}",
                bucket.0.clone()
            )
        });
        let components = bucket.1;

        let splice_map = create_splice_map(chr, &splice_plus, &splice_minus);

        let binding = HashMap::new();
        let repeats = repeats.as_ref().unwrap_or(&binding).get(chr.as_bytes());

        distribute(
            components,
            &splice_map,
            &scan_scores,
            &genome,
            &accumulator,
            args.nag,
            args.rt_freq_threshold,
            args.spliceai_min_ss_signal,
            args.intron_freq_threshold,
            args.maxent_min_ss_signal,
            &spliceosome,
            repeats,
            args.do_intron_track,
        );
    });
    log::info!("INFO: Classified {} introns", accumulator.introns.len());

    let filename = format!(
        "{}.{}",
        args.prefix.clone().unwrap_or("isotools".to_string()),
        INTRON_CLASSIFICATION
    );
    std::fs::create_dir_all(&args.outdir)
        .unwrap_or_else(|e| panic!("ERROR: Could not create output directory -> {e}!"));
    let output = args.outdir.join(filename);

    log::info!("INFO: Writing output to {}", output.display());
    let mut writer = BufWriter::new(
        File::create(output)
            .unwrap_or_else(|e| panic!("ERROR: Could not create output file -> {e}!")),
    );

    accumulator.introns.iter().for_each(|descriptor| {
        writer
            .write_all(&descriptor)
            .unwrap_or_else(|e| panic!("ERROR: Could not write intron -> {e}!"));
        writer
            .write_all(b"\n")
            .unwrap_or_else(|e| panic!("ERROR: Could not write newline -> {e}!"));
    });

    if args.do_intron_track {
        log::info!("INFO: Writing intron track");

        let mut writer = BufWriter::new(
            File::create(format!(
                "{}.introns_track,bed",
                args.prefix.unwrap_or("isotools".to_string())
            ))
            .unwrap_or_else(|e| panic!("ERROR: Could not create output file -> {e}!")),
        );

        accumulator.track.iter().for_each(|track| {
            writer
                .write_all(&track)
                .unwrap_or_else(|e| panic!("ERROR: Could not write intron -> {e}!"));
            writer
                .write_all(b"\n")
                .unwrap_or_else(|e| panic!("ERROR: Could not write newline -> {e}!"));
        });
    }

    Ok(())
}

struct ParallelAccumulator {
    introns: DashSet<Vec<u8>>,
    track: DashSet<Vec<u8>>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            introns: DashSet::new(),
            track: DashSet::new(),
        }
    }
}

/// Distributes intron buckets for parallel processing.
///
/// This function takes a vector of `BedPackage` components (intron buckets) and processes them in parallel.
/// It applies various classification rules to each intron, collecting the results in a shared accumulator.
///
/// # Arguments
///
/// * `components`: A `Vec` of `Box<dyn BedPackage>` representing intron data for a specific chromosome.
/// * `banned`: A `HashSet` of `(u64, u64)` tuples representing blacklisted introns.
/// * `splice_map`: A tuple of `SharedSpliceMap` containing SpliceAI scores for plus and minus strands.
/// * `scan_scores`: An `Option` containing `ScanScores` for MaxEnt scanning.
/// * `genome`: An `Option` containing the genome sequence data.
/// * `accumulator`: A `ParallelAccumulator` for thread-safe collection of results.
/// * `nag`: A `bool` indicating whether to scan for NULLG patterns.
///
/// # Example
///
/// ```rust, ignore
/// distribute(components, banned, &splice_map, &scan_scores, &genome, &accumulator, nag);
/// ```
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn distribute(
    components: Vec<Vec<GenePred>>,
    splice_map: &(SharedSpliceMap, SharedSpliceMap),
    scan_scores: &ScanScores,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    accumulator: &ParallelAccumulator,
    nag: bool,
    rt_frequency_threshold: f32,
    spliceai_min_ss_signal: f32,
    intron_frequency_threshold: f32,
    maxent_min_ss_signal: f32,
    spliceosome: &HashMap<String, USpliceType>,
    repeats: Option<&Lapper<u64, ()>>,
    do_intron_track: bool,
) {
    components.into_par_iter().for_each(|comp| {
        let (refs, queries) = comp.into_iter().partition(|record| {
            let role = record
                .get_extra(b"role")
                .unwrap_or_else(|| panic!("ERROR: Could not get role from record!"))
                .clone()
                .into_scalar();

            role == Some(b"reference".to_vec())
        });

        let result = process_component(
            refs,
            queries,
            splice_map,
            scan_scores,
            genome,
            nag,
            rt_frequency_threshold,
            spliceai_min_ss_signal,
            intron_frequency_threshold,
            maxent_min_ss_signal,
            spliceosome,
            repeats,
            do_intron_track,
        );

        let (lines, track_lines) = result;

        lines.into_iter().for_each(|descriptor| {
            if !descriptor.is_empty() {
                accumulator.introns.insert(descriptor);
            }
        });

        if do_intron_track {
            track_lines.into_iter().for_each(|track| {
                if !track.is_empty() {
                    accumulator.track.insert(track);
                }
            });
        }
    });
}

/// Processes a single `IntronBucket` and classifies its introns.
///
/// This function iterates through all introns in a given `IntronBucket`, applies classification logic
/// (e.g., getting splice site context, SpliceAI scores, MaxEnt scores, and checking for RT and NULLG repeats),
/// and determines the support type for each intron.
///
/// # Arguments
///
/// * `component`: A mutable reference to an `IntronBucket`.
/// * `banned`: A `HashSet` of `(u64, u64)` tuples representing blacklisted introns.
/// * `splice_map`: A tuple of `SharedSpliceMap` containing SpliceAI scores.
/// * `scan_scores`: An `Option` containing `ScanScores` for MaxEnt scanning.
/// * `genome`: An `Option` containing the genome sequence data.
/// * `nag`: A `bool` indicating whether to scan for NULLG patterns.
///
/// # Returns
///
/// * A `HashMap<(u64, u64), String>`
///   where the key is the intron coordinates and the value is a formatted string of its classification descriptor.
///
/// # Example
///
/// ```rust, ignore
/// let classified_introns = process_component(&mut component, &banned, &splice_map, &scan_scores, &genome, nag);
/// ```
#[inline(always)]
#[allow(clippy::collapsible_else_if)]
#[allow(clippy::if_same_then_else)]
#[allow(clippy::too_many_arguments)]
fn process_component(
    refs: Vec<GenePred>,
    queries: Vec<GenePred>,
    splice_map: &(SharedSpliceMap, SharedSpliceMap),
    scan_scores: &ScanScores,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    nag: bool,
    rt_frequency_threshold: f32,
    spliceai_min_ss_signal: f32,
    intron_frequency_threshold: f32,
    maxent_min_ss_signal: f32,
    spliceosome: &HashMap<String, USpliceType>,
    repeats: Option<&Lapper<u64, ()>>,
    do_intron_track: bool,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    if queries.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // INFO: on queries in case of no references -> panics
    let chr = queries[0].chrom();
    let strand = queries[0]
        .strand()
        .unwrap_or_else(|| panic!("ERROR: Could not get strand from {}!", queries[0]));

    // INFO: getting stranded spliceAi scores
    let splice_scores = match strand {
        genepred::Strand::Forward => Some(&splice_map.0),
        genepred::Strand::Reverse => Some(&splice_map.1),
        _ => None,
    };

    // INFO: collecting query introns and reference introns
    let mut ref_introns = HashSet::new();
    let mut ref_cds_introns = HashSet::new();
    let mut ref_utr_introns = HashSet::new();
    refs.iter().for_each(|ref_| {
        ref_.introns().iter().for_each(|(start, end)| {
            if !ref_introns.contains(&(*start, *end)) {
                ref_introns.insert((*start, *end));

                if *start >= ref_.thick_start().unwrap() && *end <= ref_.thick_end().unwrap() {
                    if !ref_cds_introns.contains(&(*start, *end)) {
                        ref_cds_introns.insert((*start, *end));
                    }

                    // INFO: we remove any intron that was already in the UTR
                    // INFO: rule -> any intron in CDS must be considered as CDS intron
                    if ref_utr_introns.contains(&(*start, *end)) {
                        ref_utr_introns.remove(&(*start, *end));
                    }
                } else {
                    if !ref_utr_introns.contains(&(*start, *end)) {
                        ref_utr_introns.insert((*start, *end));
                    }
                }
            }
        });
    });

    let mut spans = HashSet::new();
    queries.iter().for_each(|query| {
        spans.insert((query.start(), query.end()));
    });

    let mut nag_introns: HashMap<(u64, u64), Intron> = HashMap::new(); // collect NULLG hits
    let mut query_introns = HashMap::new();
    queries.iter().for_each(|query| {
        // INFO: determine if query is in frame
        let is_in_frame = query.cds_length() % 3 == 0;

        query.introns().iter_mut().for_each(|(start, end)| {
            match query_introns.entry((*start, *end)) {
                Entry::Occupied(mut entry) => {
                    let stats: &mut Intron = entry.get_mut();
                    stats.seen += 1;
                }
                Entry::Vacant(entry) => {
                    let key = format!(
                        "{}:{}-{}({})",
                        std::str::from_utf8(chr).unwrap(),
                        *start - 1,
                        *end + 1,
                        strand
                    );

                    let mut stats = Intron::new();
                    stats.chrom = chr.to_vec();
                    stats.start = *start;
                    stats.end = *end;
                    stats.strand = strand;
                    stats.seen = 1;
                    stats.is_in_frame = is_in_frame;
                    stats.intron_len = *end - *start;

                    // INFO: ask if intron is within how many spans
                    for (span_start, span_end) in &spans {
                        if *start >= *span_start && *end <= *span_end {
                            stats.spanned += 1;
                        }
                    }

                    // INFO: if intron is within reference introns, set is_reference_supported
                    if ref_introns.contains(&(*start, *end)) {
                        stats.is_reference_supported = true;
                    }

                    // INFO: get splice context and spliceAI scores
                    get_sj_context(
                        &(*start + 1, *end - 1),
                        &mut stats,
                        &strand,
                        chr,
                        genome,
                        scan_scores,
                    );
                    get_sj_ai_scores(&(*start, *end), &mut stats, splice_scores, &strand);

                    // INFO: get splice type from intronIC
                    let sup_key = format!(
                        "{}:{}-{}({})",
                        std::str::from_utf8(chr).unwrap(),
                        *start,
                        *end,
                        strand
                    );
                    let splice_u_type = spliceosome.get(&sup_key).unwrap_or_else(|| {
                        if stats.is_reference_supported {
                            &USpliceType::Unknown
                        } else {
                            log::debug!(
                                "DEBUG: Couldn't find splice type for {:?} using {:?} with metadata -> {:?}. Setting to UNKNOWN",
                                (*start, *end),
                                key,
                                stats
                            );
                            &USpliceType::Unknown
                        }
                    });
                    stats.splice_u_type = splice_u_type.clone();

                    // INFO: determine position in reference set
                    if ref_cds_introns.contains(&(*start, *end)) {
                        stats.intron_position = Position::CDS;
                    } else if ref_utr_introns.contains(&(*start, *end)) {
                        stats.intron_position = Position::UTR;
                    } else {
                        stats.intron_position = Position::Unknown;
                    }

                    // INFO: if repeats are provided, check if intron is within a repeat
                    if let Some(repeats) = repeats {
                        if is_intron_within_repeat(repeats, (*start, *end)) {
                            stats.within_repeat = SpanRepeat::Within;
                        } else if is_intron_covered_by_repeat(repeats, (*start, *end), 0.5) {
                            stats.within_repeat = SpanRepeat::Spans;
                        } else {
                            stats.within_repeat = SpanRepeat::Null;
                        }
                    }

                    // INFO: if --nag set and acceptor is AG and reference supported, look for NULLGNULLG
                    if nag && stats.acceptor_sequence == b"AG" && stats.is_reference_supported {
                        scan_nag_repeats(
                            &(*start, *end),
                            &mut stats,
                            &strand,
                            chr,
                            genome,
                            scan_scores,
                            splice_scores,
                            &ref_cds_introns,
                            &ref_utr_introns,
                            &mut nag_introns,
                        );
                    }

                    entry.insert(stats);
                }
            }
        });
    });

    // INFO: merging NULLG introns in after all entries are inserted
    query_introns.extend(nag_introns);

    log::debug!("DEBUG: Second pass over introns...");
    __veredict(
        &mut query_introns,
        rt_frequency_threshold,
        spliceai_min_ss_signal,
        intron_frequency_threshold,
        maxent_min_ss_signal,
    );

    let mut rows = Vec::new();
    let mut track_rows = Vec::new();

    query_introns.values().for_each(|descriptor| {
        rows.push(descriptor.as_bytes());

        if do_intron_track {
            match descriptor.support {
                SupportType::Artifact
                | SupportType::Unclear
                | SupportType::WeakRT
                | SupportType::StrongRT => {
                    track_rows.push(descriptor.as_track_bytes());
                }
                _ => {
                    log::debug!(
                        "DEBUG: Intron not classified as artifacts/unclear/RT: {:?}",
                        descriptor
                    );
                }
            }
        }
    });

    (rows, track_rows)
}

/// Determines the support type for each intron.
///
/// This function takes a vector of `Intron` descriptors and determines the support type for each intron.
/// It uses the `Intron` descriptors to determine the support type based on various criteria, such as
/// splice site context, MaxEnt scores, and intron frequency.
///
/// # Arguments
///
/// * `query_introns`: A mutable reference to a `HashMap` containing the intron descriptors.
/// * `rt_frequency_threshold`: A float representing the RT intron frequency threshold.
/// * `spliceai_min_ss_signal`: A float representing the spliceAI minimum splice signal threshold.
/// * `intron_frequency_threshold`: A float representing the intron frequency threshold.
/// * `maxent_min_ss_signal`: A float representing the MaxEnt minimum splice signal threshold.
///
/// # Example
///
/// ```rust, ignore
/// __veredict(
///     &mut query_introns,
///     rt_frequency_threshold,
///     spliceai_min_ss_signal,
///     intron_frequency_threshold,
///     maxent_min_ss_signal,
/// );
/// ```
#[allow(clippy::collapsible_else_if)]
pub fn __veredict(
    query_introns: &mut HashMap<(u64, u64), Intron>,
    rt_frequency_threshold: f32,
    spliceai_min_ss_signal: f32,
    intron_frequency_threshold: f32,
    maxent_min_ss_signal: f32,
) {
    query_introns
        .iter_mut()
        .for_each(|((_intron_start, _intron_end), descriptor)| {
            // INFO: not including NAG-derived introns bc they are already Splicing
            // WARN: strange TOGA supported RT introns -> s14	9965456	9968743
            if descriptor.is_rt_intron {
                if descriptor.is_reference_supported {
                    if descriptor.intron_position == Position::UTR {
                        // INFO: TOGA2 UTRs are not 100% reliable yet
                        match descriptor.within_repeat {
                            SpanRepeat::Within | SpanRepeat::Spans => {
                                descriptor.support = SupportType::StrongRT;
                            }
                            SpanRepeat::Null => {
                                // INFO: UTR TOGA2 support with small repeats but not within
                                // INFO: test spliceAi score and decide
                                if descriptor.splice_ai_donor >= spliceai_min_ss_signal
                                    && descriptor.splice_ai_acceptor >= spliceai_min_ss_signal
                                {
                                    descriptor.support = SupportType::Splicing;
                                } else {
                                    descriptor.support = SupportType::WeakRT;
                                }
                            }
                        }
                    } else {
                        // INFO: TOGA2 intron is splicing
                        descriptor.support = SupportType::Splicing;
                    }
                } else {
                    // INFO: intron is RT but not TOGA2 supported
                    match descriptor.within_repeat {
                        // INFO: within/spans repeats -> sure RT no matter frequency
                        SpanRepeat::Within | SpanRepeat::Spans => {
                            descriptor.support = SupportType::StrongRT;
                        }
                        SpanRepeat::Null => {
                            // INFO: completey new intron with small repeats and not within db
                            // INFO: if intron > 50% -> really unclear
                            if (descriptor.seen as f32 / descriptor.spanned as f32)
                                >= rt_frequency_threshold
                            {
                                // INFO: second assertion due to 1/1 clear artifacts passing as
                                // unclear:
                                // scaffold_1	120379762	120496991	+	1	1	0	0	0	-33.18959	tt	TG
                                // taagcaacattttcagt	GCAGCAGAGAGAGCAACATGGGT	NEW	NO_REFERENCE_SUPPORT
                                // OUT_OF_FRAME	taagcaacattt	AGAGCAACATGG	RT_INTRON	NOT_NAG_SS
                                // U2	NO_REPEAT	UNCLEAR
                                if descriptor.splice_ai_donor > 0.0
                                    && descriptor.splice_ai_acceptor > 0.0
                                    && (descriptor.max_ent_donor > 0.0
                                        || descriptor.max_ent_acceptor > 0.0)
                                {
                                    descriptor.support = SupportType::Unclear;
                                } else {
                                    descriptor.support = SupportType::Artifact;
                                }
                            } else {
                                // INFO: low frequency RT intron -> spliceAi scores decide
                                // strong/weak
                                if descriptor.splice_ai_donor >= spliceai_min_ss_signal
                                    && descriptor.splice_ai_acceptor >= spliceai_min_ss_signal
                                {
                                    descriptor.support = SupportType::WeakRT;
                                } else {
                                    descriptor.support = SupportType::StrongRT;
                                }
                            }
                        }
                    }
                }
            } else {
                let has_signal = (descriptor.max_ent_donor >= maxent_min_ss_signal
                    && descriptor.max_ent_acceptor >= maxent_min_ss_signal)
                    && (descriptor.splice_ai_donor > 0.01 && descriptor.splice_ai_acceptor > 0.01);
                let has_minimum_signal = (descriptor.max_ent_donor >= 0.0
                    && descriptor.max_ent_acceptor >= 0.0)
                    && (descriptor.splice_ai_donor >= MINIMUM_SPLICEAI_SS_SIGNAL
                        && descriptor.splice_ai_acceptor >= MINIMUM_SPLICEAI_SS_SIGNAL);

                if descriptor.is_reference_supported
                    || (descriptor.splice_ai_donor >= spliceai_min_ss_signal
                        && descriptor.splice_ai_acceptor >= spliceai_min_ss_signal)
                {
                    descriptor.support = SupportType::Splicing;
                } else if ((descriptor.seen as f32 / descriptor.spanned as f32)
                    >= intron_frequency_threshold)
                    || has_signal
                {
                    // INFO: new branch testing size threshold
                    // INFO: intron bigger than 50kb without any signal is likely an artifact
                    if descriptor.intron_len >= INTRON_SIZE_THRESHOLD {
                        if descriptor.splice_u_type == USpliceType::U12
                            || (descriptor.splice_ai_donor > 0.01
                                && descriptor.splice_ai_acceptor > 0.01)
                        {
                            descriptor.support = SupportType::Unclear;
                        } else {
                            descriptor.support = SupportType::Artifact;
                        }
                    // WARN: how do we deal with 1/2 (50%) cases?
                    // INFO: new branch for MaxEnt only -> not trusting it alone
                    // INFO: here we test if maxEnt is significant + if there is
                    // INFO: spliceAi signal [> 0.0]
                    } else if descriptor.splice_u_type == USpliceType::U12 {
                        descriptor.support = SupportType::Splicing;
                    } else if has_minimum_signal {
                        descriptor.support = SupportType::Unclear;
                    } else {
                        // INFO: intron type is unknown but is high frequency OR has ss signal -> unclear
                        descriptor.support = SupportType::Artifact;
                    }
                } else {
                    // INFO: if intron is U12 -> unclear directly
                    if descriptor.splice_u_type == USpliceType::U12 || has_minimum_signal {
                        descriptor.support = SupportType::Unclear;
                    } else {
                        // INFO: intron is U2 or unknown + low freq + low ss signal
                        descriptor.support = SupportType::Artifact
                    }
                }
            }
        });
}

/// Extracts the genomic sequence context around the splice junction.
///
/// This function retrieves the donor and acceptor splice site sequences and their surrounding context from the genome.
/// It also handles reverse strand sequences by taking the reverse complement. It then calls functions to calculate
/// MaxEnt scores and scan for RT repeats.
///
/// # Arguments
///
/// * `intron`: A tuple `(u64, u64)` representing the start and end coordinates of the intron.
/// * `descriptor`: A mutable reference to an `Intron` struct to store the extracted sequences and scores.
/// * `strand`: The `Strand` of the intron (`Forward` or `Reverse`).
/// * `chr`: The chromosome name as a `String`.
/// * `genome`: An `Option` containing the genome sequence data.
/// * `scan_scores`: An `Option` containing `ScanScores` for MaxEnt scanning.
///
/// # Example
///
/// ```rust, ignore
/// get_sj_context(&intron, &mut descriptor, &strand, &chr, &genome, &scan_scores);
/// ```
fn get_sj_context(
    intron: &(u64, u64),
    descriptor: &mut Intron,
    strand: &genepred::Strand,
    chr: &[u8],
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    scan_scores: &ScanScores,
) {
    let intron_start = intron.0;
    let intron_end = intron.1;

    if !genome.is_empty() {
        match strand {
            genepred::Strand::Forward => {
                // INFO: For TOGA nag we need +3 upstream and +3 downstream
                // INFO: For RT repeats we need 12-mers (11 exon + 1 intron) at 5'
                // INFO: and 11 intron + 1 exon at 3'

                let donor_context = Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read donor context!")
                        [intron_start as usize - 12..intron_start as usize + 5]
                        .as_ref(),
                );
                let donor_seq = donor_context.slice(11, 13);
                let donor_rt_context = donor_context.slice(0, 12);

                let acceptor_context = Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read acceptor context!")
                        [intron_end as usize - 19..intron_end as usize + 4]
                        .as_ref(),
                );
                let acceptor_seq = acceptor_context.slice(18, 20);
                let acceptor_rt_context = acceptor_context.slice(9, 21);

                descriptor.donor_sequence = donor_seq;
                descriptor.donor_context = donor_context;

                descriptor.acceptor_sequence = acceptor_seq;
                descriptor.acceptor_context = acceptor_context;

                descriptor.donor_rt_context = donor_rt_context;
                descriptor.acceptor_rt_context = acceptor_rt_context;
            }
            genepred::Strand::Reverse => {
                let donor_context = Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read donor context!")
                        [intron_end as usize - 5..intron_end as usize + 12]
                        .as_ref(),
                )
                .reverse_complement();
                let donor_seq = donor_context.slice(11, 13);
                let donor_rt_context = donor_context.slice(0, 12);

                let acceptor_context = Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read acceptor context!")
                        [intron_start as usize - 4..intron_start as usize + 19]
                        .as_ref(),
                )
                .reverse_complement();
                let acceptor_seq = acceptor_context.slice(18, 20);
                let acceptor_rt_context = acceptor_context.slice(9, 21);

                descriptor.donor_sequence = donor_seq;
                descriptor.donor_context = donor_context;

                descriptor.acceptor_sequence = acceptor_seq;
                descriptor.acceptor_context = acceptor_context;

                descriptor.donor_rt_context = donor_rt_context;
                descriptor.acceptor_rt_context = acceptor_rt_context;
            }
            _ => panic!("ERROR: Unknown strand {:?}!", strand),
        }

        get_sj_max_entropy(descriptor, scan_scores);
        scan_rt_repeats(descriptor);
    }
}

/// Calculates the MaxEnt scores for the donor and acceptor splice sites.
///
/// This function looks up the pre-calculated MaxEnt scores for the splice site contexts stored in the
/// `scan_scores` database and updates the intron descriptor.
///
/// # Arguments
///
/// * `descriptor`: A mutable reference to an `Intron` struct.
/// * `scan_scores`: An `Option` containing `ScanScores` for MaxEnt lookup.
///
/// # Example
///
/// ```rust, ignore
/// get_sj_max_entropy(&mut descriptor, &scan_scores);
/// ```
fn get_sj_max_entropy(descriptor: &mut Intron, scan_scores: &ScanScores) {
    if let Some(scan_scores) = scan_scores {
        let (donor_score_map, acceptor_score_map): &(SpliceScoreMap, SpliceScoreMap) = scan_scores;

        let donor_max_ent_context = &descriptor.donor_context.slice(8, 17);

        if donor_max_ent_context.len() != MAX_ENT_DONOR_MIN_SIZE {
            eprintln!("ERROR: Donor context is not 9 bases long!");
        }

        let donor_score = donor_score_map
            .get(donor_max_ent_context)
            .and_then(|r| r.first())
            .unwrap_or(&0.0);

        descriptor.max_ent_donor = *donor_score as f32;

        let acceptor_score =
            calculate_acceptor_score(&descriptor.acceptor_context, acceptor_score_map);

        descriptor.max_ent_acceptor = acceptor_score as f32;
    }
}

/// Retrieves the SpliceAI scores for the donor and acceptor splice sites.
///
/// This function queries the SpliceAI score maps for the given intron's
/// donor and acceptor coordinates and updates the intron descriptor
/// with the found scores.
///
/// # Arguments
///
/// * `intron`: A tuple `(u64, u64)` representing the start and end coordinates of the intron.
/// * `descriptor`: A mutable reference to an `Intron` struct.
/// * `splice_scores`: An `Option` containing a reference to `SharedSpliceMap`.
/// * `strand`: The `Strand` of the intron.
///
/// # Example
///
/// ```rust, ignore
/// get_sj_ai_scores(&intron, &mut descriptor, splice_scores, &strand);
/// ```
fn get_sj_ai_scores(
    intron: &(u64, u64),
    descriptor: &mut Intron,
    splice_scores: Option<&SharedSpliceMap>,
    strand: &Strand,
) {
    let intron_start = intron.0;
    let intron_end = intron.1;

    if let Some(splice_scores) = splice_scores {
        if let Some(donor_score_map) = splice_scores.0.as_ref() {
            let acceptor_score_map = splice_scores
                .1
                .as_ref()
                .expect("ERROR: Acceptor score map is None, this is a bug!");

            // donor(+)/acceptor(-) [-1 to match bigtools coords]
            let (intron_donor, intron_acceptor) = match strand {
                genepred::Strand::Forward => (intron_start as usize, intron_end as usize - 1),
                genepred::Strand::Reverse => (intron_end as usize - 1, intron_start as usize),
                _ => panic!("ERROR: Unknown strand {:?}!", strand),
            };

            let (donor_score, acceptor_score) = (
                donor_score_map
                    .get(&intron_donor)
                    .map(|r| *r)
                    .unwrap_or(0.0),
                acceptor_score_map
                    .get(&intron_acceptor)
                    .map(|r| *r)
                    .unwrap_or(0.0),
            );

            descriptor.splice_ai_donor = donor_score;
            descriptor.splice_ai_acceptor = acceptor_score;
        }
    } else {
        log::warn!("WARN: No spliceAI scores to grab values from!");
    }
}

/// Processes a potential NULLG-derived intron.
///
/// This function creates a new intron descriptor for a NULLG-derived
/// intron, retrieves its genomic context and splice site scores,
/// and formats a new descriptor string. This new intron is considered
/// a splicing intron.
///
/// # Arguments
///
/// * `base_intron`: A tuple `(u64, u64)` of the original intron coordinates.
/// * `offset`: An `i64` representing the coordinate shift for the new intron.
/// * `strand`: The `Strand` of the intron.
/// * `chr`: The chromosome name.
/// * `genome`: An `Option` containing the genome sequence data.
/// * `scan_scores`: An `Option` containing `ScanScores`.
/// * `splice_scores`: An `Option` containing a reference to `SharedSpliceMap`.
///
/// # Returns
///
/// * A `String` containing the formatted descriptor for the new NULLG intron.
///
/// # Example
///
/// ```rust, ignore
/// let nag_descriptor = process_nag_pattern(&intron, -3, &strand, &chr, &genome, &scan_scores, splice_scores);
/// ```
#[allow(clippy::too_many_arguments)]
fn process_nag_pattern(
    base_intron: &(u64, u64),
    offset: i64,
    strand: &Strand,
    chr: &[u8],
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    scan_scores: &ScanScores,
    splice_scores: Option<&SharedSpliceMap>,
    ref_cds_introns: &HashSet<(u64, u64)>,
    ref_utr_introns: &HashSet<(u64, u64)>,
    acc: &mut HashMap<(u64, u64), Intron>,
) {
    // INFO: move the acceptor position by offset [look for NULLG]
    let intron = match strand {
        Strand::Forward => (base_intron.0, (base_intron.1 as i64 + offset) as u64),
        Strand::Reverse => (base_intron.0 - (offset as u64), base_intron.1),
        _ => panic!("ERROR: Unknown strand {:?}!", strand),
    };

    if acc.contains_key(&intron) {
        log::debug!(
            "DEBUG: NULLG pattern already in accumulator -> {:?}",
            intron
        );

        let descriptor = acc.get_mut(&intron).unwrap_or_else(|| {
            panic!(
                "ERROR: could not retrieve descriptor from accumulator -> {:?}",
                intron
            )
        });

        descriptor.is_nag_intron = true;

        log::debug!(
            "DEBUG: Will replace seen intron with NAG flag -> {:?}",
            intron
        );
    } else {
        log::debug!("DEBUG: Found and insert potential NAG {:?}!", intron);

        let mut new_descriptor = Intron::new();

        new_descriptor.chrom = chr.to_vec();
        new_descriptor.start = intron.0;
        new_descriptor.end = intron.1;
        new_descriptor.strand = *strand;
        new_descriptor.intron_len = intron.1 - intron.0;

        get_sj_context(
            &(intron.0 + 1, intron.1 - 1),
            &mut new_descriptor,
            strand,
            chr,
            genome,
            scan_scores,
        );
        get_sj_ai_scores(&intron, &mut new_descriptor, splice_scores, strand);
        new_descriptor.is_nag_intron = true;

        // INFO: determine position in reference set
        if ref_cds_introns.contains(&intron) {
            new_descriptor.intron_position = Position::CDS;
        } else if ref_utr_introns.contains(&intron) {
            new_descriptor.intron_position = Position::UTR;
        } else {
            new_descriptor.intron_position = Position::Unknown;
        }

        // WARN: NAG-derived introns will be always SupportType::Splicing
        new_descriptor.support = SupportType::Splicing;
        acc.insert(intron, new_descriptor);
    }
}

/// Scans for NULLG patterns in the acceptor context.
///
/// If the pre- or post-acceptor sequence matches one of the NULLG
/// patterns (`CAG`, `TAG`, `AAG`), this function calls `process_nag_pattern`
/// to create and classify new introns at the wiggled splice sites.
///
/// # Arguments
///
/// * `intron`: A tuple `(u64, u64)` of the original intron coordinates.
/// * `descriptor`: A mutable reference to the `Intron` struct.
/// * `strand`: The `Strand` of the intron.
/// * `chr`: The chromosome name.
/// * `genome`: An `Option` containing the genome sequence data.
/// * `scan_scores`: An `Option` containing `ScanScores`.
/// * `splice_scores`: An `Option` containing a reference to `SharedSpliceMap`.
/// * `acc`: A mutable `HashMap` to store the new NULLG introns.
///
/// # Example
///
/// ```rust, ignore
/// scan_nag_repeats(&intron, &mut descriptor, &strand, &chr, &genome, &scan_scores, splice_scores, &mut acc);
/// ```
#[allow(clippy::too_many_arguments)]
fn scan_nag_repeats(
    intron: &(u64, u64),
    descriptor: &mut Intron,
    strand: &Strand,
    chr: &[u8],
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    scan_scores: &ScanScores,
    splice_scores: Option<&SharedSpliceMap>,
    ref_cds_introns: &HashSet<(u64, u64)>,
    ref_utr_introns: &HashSet<(u64, u64)>,
    acc: &mut HashMap<(u64, u64), Intron>,
) {
    let pre_acceptor = &descriptor.acceptor_context.slice(14, 17);
    let post_acceptor = &descriptor.acceptor_context.slice(20, 23);

    if NULLG_PATTERNS.contains(&pre_acceptor.as_slice()) {
        process_nag_pattern(
            intron,
            -3,
            strand,
            chr,
            genome,
            scan_scores,
            splice_scores,
            ref_cds_introns,
            ref_utr_introns,
            acc,
        );
    }

    if NULLG_PATTERNS.contains(&post_acceptor.as_slice()) {
        process_nag_pattern(
            intron,
            3,
            strand,
            chr,
            genome,
            scan_scores,
            splice_scores,
            ref_cds_introns,
            ref_utr_introns,
            acc,
        );
    }
}

/// Scans for RT repeats between the donor and acceptor splice sites.
///
/// This function checks if there are repeated sequences between the donor and acceptor contexts,
/// which is an indicator of an RT-driven intron.
///
/// # Arguments
///
/// * `descriptor`: A mutable reference to an `Intron` struct.
///
/// # Example
///
/// ```rust, ignore
/// scan_rt_repeats(&mut descriptor);
/// ```
fn scan_rt_repeats(descriptor: &mut Intron) {
    unsafe {
        scan_sequence(descriptor);
    }
}

/// Performs the core logic for detecting RT repeats.
///
/// This unsafe function compares k-mers (8-mers in this case) from the donor and acceptor contexts.
/// If a match with a Hamming distance less than or equal to `MISMATCHES` is found, it sets the `is_rt_intron`
/// flag to true and returns.
///
/// # Arguments
///
/// * `descriptor`: A mutable reference to an `Intron` struct.
///
/// # Example
///
/// ```rust, ignore
/// unsafe {
///     scan_sequence(&mut descriptor);
/// }
/// ```
#[inline(always)]
unsafe fn scan_sequence(descriptor: &mut Intron) {
    let donor = &descriptor.donor_rt_context;
    let acceptor = &descriptor.acceptor_rt_context;

    if donor.len() != RT_REPEAT || acceptor.len() != RT_REPEAT {
        log::error!("ERROR: RT context is not 12 bases long, this is a bug!");
        std::process::exit(1);
    }

    #[inline(always)]
    fn base_to_bits(base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 0,
            _ => panic!("ERROR: Invalid base -> {base}"),
        }
    }

    // INFO: convert 8-mer to u16 (2 bits per base = 16 bits total)
    #[inline(always)]
    fn window_to_int(window: &[u8]) -> u16 {
        let mut value: u16 = 0;
        for &base in window {
            value = (value << 2) | (base_to_bits(base) as u16);
        }
        value
    }

    #[inline(always)]
    fn hamming_distance(a: u16, b: u16) -> u32 {
        (a ^ b).count_ones()
    }

    if donor.len() < WINDOW_SIZE || acceptor.len() < WINDOW_SIZE {
        return;
    }

    // INFO: build a set of all 8-mers in seq1
    let mut kmers = HashSet::with_capacity(donor.len() - WINDOW_SIZE + 1);
    for idx in 0..=donor.len() - WINDOW_SIZE {
        let kmer = window_to_int(&donor[idx..idx + WINDOW_SIZE]);
        kmers.insert((kmer, idx));
    }

    // INFO: compare k-mers from seq2 against the set
    for idx in 0..=acceptor.len() - WINDOW_SIZE {
        let kmer2 = window_to_int(&acceptor[idx..idx + WINDOW_SIZE]);

        for &(kmer1, _) in &kmers {
            if hamming_distance(kmer1, kmer2) <= MISMATCHES {
                // INFO: returning bool and skipping early
                // repeats.push((i, j));
                // dbg!("Match found: seq1[{}..{}], seq2[{}..{}]", _, _ + 8, idx, idx + 8);
                descriptor.is_rt_intron = true;
                return;
            }
        }
    }

    descriptor.is_rt_intron = false;
}

/// Wiggles splice sites and checks for support in a reference set.
///
/// This function is currently unused (`#[allow(unused)]`) but is intended to check if slightly shifted splice
/// sites are present in a set of reference introns. This could potentially be used to recover or reclassify
/// introns that are just a few bases off from a known splicing event.
///
/// # Arguments
///
/// * `acc`: A mutable `HashMap` to store and update intron classification results.
/// * `intron`: A tuple `(u64, u64)` of the intron coordinates.
/// * `ref_introns`: A `HashSet` of reference intron coordinates.
/// * `strand`: The `Strand` of the intron.
///
/// # Example
///
/// ```rust, ignore
/// wiggle_splice_sites(&mut acc, &intron, &ref_introns, &strand);
/// ```
#[allow(unused)]
#[deprecated(note = "This function is no longer used")]
fn wiggle_splice_sites(
    acc: &mut HashMap<(u64, u64), String>,
    intron: &(u64, u64),
    ref_introns: &HashSet<(u64, u64)>,
    strand: &Strand,
) {
    fn process_intron(
        acc: &mut HashMap<(u64, u64), String>,
        ref_introns: &HashSet<(u64, u64)>,
        intron_start: u64,
        intron_end: u64,
        strand: &Strand,
    ) {
        // INFO: look if swithed intron is inside ref_introns
        if ref_introns.contains(&(intron_start, intron_end)) {
            let (intron_start, intron_end) = match strand {
                Strand::Forward => (intron_start, intron_end),
                Strand::Reverse => (intron_end, intron_start),
                _ => panic!("ERROR: Unknown strand {:?}!", strand),
            };

            // INFO: if its already classified, change support and toga_support
            if acc.contains_key(&(intron_start, intron_end)) {
                let mut fields = acc
                    .get(&(intron_start, intron_end))
                    .unwrap()
                    .split('\t')
                    .collect::<Vec<_>>();

                fields[15] = "true";
                fields[21] = "SPLICED";

                acc.insert((intron_start, intron_end), fields.join("\t"));
            } else {
                // WARN: if its not classified, we need to insert it, leaving the rest of the fields empty
                // WARN: the String inserted here should be unreacheable!
                acc.insert((intron_start, intron_end), String::new());
            }
        }
    }

    for wiggle in WIGGLE_SWITCH {
        // INFO: backward
        process_intron(
            acc,
            ref_introns,
            intron.0 - wiggle as u64,
            intron.1 - wiggle as u64,
            strand,
        );

        // INFO: forward
        process_intron(
            acc,
            ref_introns,
            intron.0 + wiggle as u64,
            intron.1 + wiggle as u64,
            strand,
        );
    }
}

/// Reads intronIC output and returns a `HashMap` of intron IDs to `USpliceType`.
///
/// # Arguments
///
/// * `path`: A `PathBuf` representing the path to the intronIC output file.
///
/// # Returns
///
/// * A `HashMap` containing the intron IDs as keys and `USpliceType` as values.
///
/// # Example
///
/// ```rust, ignore
/// let intron_types = read_iic(path);
/// ```
pub fn read_iic(path: PathBuf) -> HashMap<String, USpliceType> {
    // INFO: format intronIC output:
    // chr10:75482935-75483890(-)  -90.0   GT-AG ...  .   .   u2  .
    // INFO: we are interested in the first and n-1 columns [coords, type]
    let reader = BufReader::new(
        File::open(&path)
            .unwrap_or_else(|e| panic!("ERROR: Could not open intronIC output file -> {e}!")),
    );

    log::debug!("DEBUG: opened intronIC output file -> {path:?}!");

    let mut table = HashMap::new();

    for line in reader.lines().skip(1) {
        let record =
            line.unwrap_or_else(|e| panic!("ERROR: Could not read intronIC output line -> {e}!",));

        let mut fields = record.split('\t');

        let id = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not read intronIC id from {:?}", &record))
            .to_string();
        let splice_type = fields
            .nth(11)
            .unwrap_or_else(|| panic!("ERROR: Could not read intronIC type from {:?}!", &record));

        let splice_type = match splice_type {
            "u2" => USpliceType::U2,
            "u12" => USpliceType::U12,
            // INFO: SUPER_6:7202593-7202620(-);[o:s] NA tt-AC ... NA NA NA omitted_short
            "NA" => USpliceType::Unknown,
            _ => panic!("ERROR: Unknown intronIC type -> {splice_type}!"),
        };

        table.insert(id, splice_type);
    }

    log::info!(
        "INFO: Read {} intronIC types from {}",
        table.len(),
        path.display()
    );

    table
}

/// Splice site type (U2, U12, Unknown)
///
/// This enum is used to store the type of splice site.
///
/// # Example
///
/// ```rust, no_run
/// use iso::SpliceSite;
///
/// let donor = SpliceSite::Donor;
/// let acceptor = SpliceSite::Acceptor;
/// ```
#[derive(Debug, Eq, PartialEq, Clone)]
pub enum USpliceType {
    U2,
    U12,
    Unknown,
}

/// Implements std::fmt::Display for USpliceType
impl std::fmt::Display for USpliceType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            USpliceType::U2 => write!(f, "U2"),
            USpliceType::U12 => write!(f, "U12"),
            USpliceType::Unknown => write!(f, "UNKNOWN"),
        }
    }
}

/// Implements std::str::FromStr for USpliceType
impl std::str::FromStr for USpliceType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "U2" => Ok(USpliceType::U2),
            "U12" => Ok(USpliceType::U12),
            "Unknown" => Ok(USpliceType::Unknown),
            _ => Err(format!("ERROR: Unknown splice type -> {s}")),
        }
    }
}

/// Holds a variety of statistical and contextual data for a predicted intron.
///
/// This struct contains metrics from different prediction tools and classification
/// systems, used to evaluate the quality and nature of a predicted intron.
#[derive(Debug, PartialEq, Clone)]
pub struct Intron {
    /// Chromosome name
    pub chrom: Vec<u8>,
    /// Start position of the intron (1-based)
    pub start: u64,
    /// End position of the intron (1-based)
    pub end: u64,
    /// Strand of the intron
    pub strand: genepred::Strand,
    /// The frequency of how many reads contain this intron.
    pub seen: usize,
    /// The frequency of how many reads span this intron.
    pub spanned: usize,
    /// SpliceAI score for the donor site.
    pub splice_ai_donor: f32,
    /// SpliceAI score for the acceptor site.
    pub splice_ai_acceptor: f32,
    /// MaxEntScan score for the donor site.
    pub max_ent_donor: f32,
    /// MaxEntScan score for the acceptor site.
    pub max_ent_acceptor: f32,
    /// The sequence around the donor site.
    pub donor_sequence: Vec<u8>,
    /// The sequence around the acceptor site.
    pub acceptor_sequence: Vec<u8>,
    /// The MaxEntScan 9-mer donor context sequence.
    pub donor_context: Sequence,
    /// The MaxEntScan 23-mer acceptor context sequence.
    pub acceptor_context: Sequence,
    /// The classification of the intron's position according to reference set.
    pub intron_position: Position,
    /// A boolean indicating if the intron is supported by reference set.
    pub is_reference_supported: bool,
    /// A boolean indicating if the intron maintains the reading frame.
    pub is_in_frame: bool,
    /// The RT-switch context sequence for the donor site.
    pub donor_rt_context: Vec<u8>,
    /// The RT-switch context sequence for the acceptor site.
    pub acceptor_rt_context: Vec<u8>,
    /// A boolean indicating if the intron is an RT-switch intron.
    pub is_rt_intron: bool,
    /// A boolean indicating if the intron is a reference-nag intron.
    pub is_nag_intron: bool,
    /// A classification of the intron's splice type.
    pub splice_u_type: USpliceType,
    /// A boolean indicating if the intron is within a repeat.
    pub within_repeat: SpanRepeat,
    /// Length of the intron
    pub intron_len: u64,
    /// A classification of the intron's support type.
    pub support: SupportType,
}

impl Intron {
    /// Creates a new `Intron` instance with all fields initialized to default values.
    ///
    /// This is a convenience constructor for creating a blank statistics object.
    ///
    /// # Returns
    ///
    /// A new `Intron` instance with default values.
    ///
    pub fn new() -> Self {
        Self {
            chrom: Vec::new(),
            start: 0,
            end: 0,
            strand: genepred::Strand::Unknown,
            seen: 0,
            spanned: 0,
            splice_ai_donor: 0.0,
            splice_ai_acceptor: 0.0,
            max_ent_donor: 0.0,
            max_ent_acceptor: 0.0,
            donor_sequence: Vec::new(),
            acceptor_sequence: Vec::new(),
            donor_context: Sequence::new(&[]),
            acceptor_context: Sequence::new(&[]),
            intron_position: Position::Unknown,
            is_reference_supported: false,
            is_in_frame: false,
            donor_rt_context: Vec::new(),
            acceptor_rt_context: Vec::new(),
            is_rt_intron: false,
            is_nag_intron: false,
            splice_u_type: USpliceType::Unknown,
            within_repeat: SpanRepeat::Null,
            intron_len: 0,
            support: SupportType::Unclear,
        }
    }

    pub fn as_bytes(&self) -> Vec<u8> {
        format!("{}", self).as_bytes().to_vec()
    }

    pub fn as_track_bytes(&self) -> Vec<u8> {
        format!(
            "{}\t{}\t{}\t{}",
            std::str::from_utf8(&self.chrom).unwrap_or("NULL"),
            self.start,
            self.end,
            self.support
        )
        .as_bytes()
        .to_vec()
    }
}

/// Implements std::fmt::Display for Intron
impl std::fmt::Display for Intron {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            std::str::from_utf8(&self.chrom).unwrap_or("NULL"),
            self.start,
            self.end,
            self.strand,
            self.seen,
            self.spanned,
            self.splice_ai_donor,
            self.splice_ai_acceptor,
            self.max_ent_donor,
            self.max_ent_acceptor,
            std::str::from_utf8(&self.donor_sequence).unwrap_or("NULL"),
            std::str::from_utf8(&self.acceptor_sequence).unwrap_or("NULL"),
            std::str::from_utf8(&self.donor_context.seq).unwrap_or("NULL"),
            std::str::from_utf8(&self.acceptor_context.seq).unwrap_or("NULL"),
            self.intron_position,
            if self.is_reference_supported { "REFERENCE_SUPPORT" } else { "NO_REFERENCE_SUPPORT" },
            if self.is_in_frame { "IN_FRAME" } else { "OUT_OF_FRAME" },
            std::str::from_utf8(&self.donor_rt_context).unwrap_or("NULL"),
            std::str::from_utf8(&self.acceptor_rt_context).unwrap_or("NULL"),
            if self.is_rt_intron { "RT_INTRON" } else { "NOT_RT_INTRON" },
            if self.is_nag_intron { "NAG_SS" } else { "NOT_NAG_SS" },
            self.splice_u_type,
            self.within_repeat,
            self.intron_len,
            self.support
        )
    }
}

/// Position of the intron in the reference
#[derive(Debug, Eq, PartialEq, Clone)]
#[allow(clippy::upper_case_acronyms)]
pub enum Position {
    CDS,
    UTR,
    Unknown,
}

/// Implements std::fmt::Display for Position
impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Position::CDS => write!(f, "CDS"),
            Position::UTR => write!(f, "UTR"),
            Position::Unknown => write!(f, "NEW"),
        }
    }
}

/// Support type
///
/// This enum is used to store the type of support.
///
/// # Example
///
/// ```rust, no_run
/// use iso::SupportType;
///
/// let spliced = SupportType::Splicing;
/// let rt = SupportType::RT;
/// let unclear = SupportType::Unclear;
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum SupportType {
    Splicing,
    StrongRT,
    WeakRT,
    Unclear,
    Artifact,
}

impl std::fmt::Display for SupportType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SupportType::Splicing => write!(f, "SPLICED"),
            SupportType::StrongRT => write!(f, "STRONG_RT"),
            SupportType::WeakRT => write!(f, "WEAK_RT"),
            SupportType::Unclear => write!(f, "UNCLEAR"),
            SupportType::Artifact => write!(f, "ARTIFACT"),
        }
    }
}

impl std::str::FromStr for SupportType {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "SPLICED" => Ok(SupportType::Splicing),
            "STRONG_RT" => Ok(SupportType::StrongRT),
            "WEAK_RT" => Ok(SupportType::WeakRT),
            "UNCLEAR" => Ok(SupportType::Unclear),
            "ARTIFACT" => Ok(SupportType::Artifact),
            _ => Err("ERROR: Cannot parse support type!".into()),
        }
    }
}

/// Span repeat type
///
/// This enum is used to store the type of span repeat.
///
/// # Example
///
/// ```rust, no_run
/// use iso::SpanRepeat;
///
/// let within = SpanRepeat::Within;
/// let spans = SpanRepeat::Spans;
/// let null = SpanRepeat::Null;
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum SpanRepeat {
    Within,
    Spans,
    Null,
}

/// Implements std::fmt::Display for SpanRepeat
impl std::fmt::Display for SpanRepeat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpanRepeat::Within => write!(f, "WITHIN_REPEAT"),
            SpanRepeat::Spans => write!(f, "SPANS_REPEAT"),
            SpanRepeat::Null => write!(f, "NO_REPEAT"),
        }
    }
}
