//! Core module for detecting orphans in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the core functions for detecting orphans in a query set of reads
//! and processing the components in a parallel fashion.
//!
//! In short, each component of reads is subjected to any of the two modes: guided or
//! self-guided. Guided mode is the default mode and is used when the user provides a TOGA file
//! as input. Self-guided mode is used when the user does not provide a TOGA file as input.
//! Both, guided and self-guided, cover an extensive amount of curated orphan cases under
//! the assumption that they do not represent a valid source of evidence for transcription.
//! The process is heavily parallellized to offer fast performance on large datasets.

use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use crate::{
    cli::{Args, Mode},
    utils::*,
};

use dashmap::DashMap;
use genepred::{Bed12, GenePred};
use packbed::{OverlapType, Role};
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
///     toga: None,
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

    match mode {
        Mode::Guided => {
            log::info!("INFO: Running using guided mode");

            let mut inputs = args.refs.unwrap_or_else(|| {
                log::error!("ERROR: No TOGA file provided while using reference guided mode!");
                std::process::exit(1);
            });

            let mut modes = std::iter::repeat(Role::Reference)
                .take(inputs.len())
                .collect::<Vec<_>>();
            modes.extend(vec![Role::Query]);
            inputs.extend(vec![args.query]);

            // CASE: single-exon-component reads with non-single-exon TOGA refs [SOLVED]
            // INFO:  the following case is presented:
            //
            // [SOLVED: CDS overlap mode]
            // toga:  xxx------------XXXX---
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
            let tracks =
                packbed::pack(inputs, vec![Role::Reference, Role::Query], OverlapType::CDS)
                    .unwrap_or_else(|e| {
                        log::error!("ERROR: Failed to packed reads: {:?}", e);
                        std::process::exit(1);
                    });

            __process(
                tracks,
                &mode,
                args.min_read_num_denovo,
                outdir,
                args.prefix,
                args.min_discard_percentage,
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
            let tracks = packbed::pack(vec![args.query], vec![Role::Query], OverlapType::CDS)
                .unwrap_or_else(|e| {
                    log::error!("ERROR: Failed to packed reads: {:?}", e);
                    std::process::exit(1);
                });

            __process(
                tracks,
                &mode,
                args.min_read_num_denovo,
                outdir,
                args.prefix,
                args.min_discard_percentage,
            );
        }
    };
}

/// Processes the components and writes the results to a file
///
/// # Arguments
///
/// * `tracks` - The components to process
/// * `mode` - The mode to use
/// * `min_read_num_denovo` - The min_read_num_denovo to use
/// * `outdir` - The output directory
/// * `filename` - The filename to use
///
/// # Returns
///
/// None
///
/// # Examples
///
/// ```
/// use isotools::iso_orphan::__process;
/// use isotools::cli::Mode;
/// use isotools::utils::Components;
///
/// let tracks: Components = DashMap::new();
/// let mode = Mode::Guided;
/// let min_read_num_denovo = 5;
/// let outdir = "tests/data".into();
/// let filename = "test".into();
///
/// __process(tracks, &mode, min_read_num_denovo, outdir, filename);
/// ```
fn __process(
    tracks: Components,
    mode: &Mode,
    min_read_num_denovo: usize,
    outdir: PathBuf,
    filename: String,
    min_discard_percentage: f32,
) {
    let accumulator = ParallelAccumulator::default();
    let counter = ParallelCounter::default();

    tracks.into_par_iter().for_each(|bucket| {
        let _chr = bucket.0;
        let components = bucket.1;

        counter.inc_components(components.len() as u32);

        __process_components(
            components,
            &accumulator,
            &counter,
            mode,
            min_read_num_denovo,
            min_discard_percentage,
        );
    });

    log::info!("INFO: Orphans found: {}", accumulator.num_orphans());

    __report_stats(&counter);

    let pass = format!("{}.orphan_free.bed", filename);
    let orphans = format!("{}.orphans.bed", &filename);

    let mut p_writer = BufWriter::new(File::create(outdir.join(pass)).unwrap());
    let mut o_writer = BufWriter::new(File::create(outdir.join(orphans)).unwrap());

    accumulator.keep.into_iter().for_each(|line| {
        p_writer.write_all(&line).unwrap();
    });

    accumulator.orphans.into_iter().for_each(|line| {
        o_writer.write_all(&line).unwrap();
    });
}

/// Parallel processing of components
///
/// # Arguments
///
/// * `components` - The components to process
/// * `accumulator` - The accumulator to use
/// * `counter` - The counter to use
/// * `mode` - The mode to use
/// * `min_read_num_denovo` - The min_read_num_denovo to use
///
/// # Returns
///
/// None
///
/// # Examples
///
/// ```
/// use isotools::iso_orphan::__process_components;
/// use isotools::cli::Mode;
/// use isotools::utils::{Components, ParallelAccumulator, ParallelCounter};
///
/// let components: Vec<Box<dyn BedPackage>> = vec![];
/// let accumulator = ParallelAccumulator::default();
/// let counter = ParallelCounter::default();
/// let mode = Mode::Guided;
/// let min_read_num_denovo = 5;
///
/// __process_components(components, &accumulator, &counter, &mode, min_read_num_denovo);
/// ```
fn __process_components(
    components: Vec<Vec<GenePred>>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
    mode: &Mode,
    min_read_num_denovo: usize,
    min_discard_percentage: f32,
) {
    components.into_par_iter().for_each(|component| match mode {
        Mode::Guided => {
            let (refs, queries) = component.into_iter().partition(|record| {
                let role = record
                    .get_extra(b"role")
                    .unwrap_or_else(|| panic!("ERROR: Could not get role from record!"))
                    .clone()
                    .into_scalar();

                role == Some(b"reference".to_vec())
            });

            let (keep, orphans) = guided(
                refs,
                queries,
                counter,
                min_read_num_denovo,
                min_discard_percentage,
            );
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

            let (keep, orphans) = self_guided(queries, counter, min_read_num_denovo);
            accumulator.add(keep, orphans);
        }
    });
}

/// Guided mode
///
/// # Arguments
///
/// * `components` - The components to process
/// * `counter` - The counter to use
/// * `min_read_num_denovo` - The min_read_num_denovo to use
///
/// # Returns
///
/// A tuple containing the keep and orphans
///
/// # Examples
///
/// ```
/// use isotools::iso_orphan::guided;
/// use isotools::utils::ParallelCounter;
///
/// let components = Box::new((RefGenePred::default(), vec![]));
/// let counter = ParallelCounter::default();
/// let min_read_num_denovo = 5;
///
/// let (keep, orphans) = guided(components, &counter, min_read_num_denovo);
/// ```
#[allow(clippy::boxed_local)]
fn guided(
    references: Vec<GenePred>,
    mut queries: Vec<GenePred>,
    counter: &ParallelCounter,
    min_read_num_denovo: usize,
    min_discard_percentage: f32,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut keep = Vec::new();
    let mut orphans = Vec::new();

    // INFO: single component reads -> no TOGA
    if references.is_empty() {
        // CASE: single-exon-component reads [PARTIALLY SOLVED]
        // INFO: the following case is presented:
        //
        // [SOLVED: single-component reads will be just discarded]
        // toga:  XX-----------XXX--XX---XXXX
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
        // toga:  XX-----------XXX--XX---XXXX
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
        do_self_guided_check(&mut queries, counter, min_read_num_denovo);
    }

    // INFO: ask if ANY TOGA is a single-exon
    // INFO: single-exon signifies no introns
    let is_reference_single_exon = references
        .iter()
        .any(|projection| projection.introns().is_empty());

    let reference_exons = references
        .iter()
        .flat_map(|record| record.exons())
        .collect::<Vec<(u64, u64)>>();

    if queries.len() > 1 {
        let mut discards = 0; // INFO: only in common branch
        queries.sort_by(|a, b| a.exonic_length().cmp(&b.exonic_length()));

        for read in queries.iter() {
            let is_single_exon_read = read.introns().is_empty();
            let line = read.to_bed::<Bed12>();

            if is_single_exon_read {
                if is_reference_single_exon {
                    counter.inc_read_se_mc_toga_se();

                    // CASE: single-exon read(s) with single-exon TOGA refs
                    // INFO: the following case is presented:
                    //
                    // toga: xxxxxxxxxxxxXXXXXxxxxxxx
                    // read1:     xxxxxxxXXXXXxxxx
                    // read2:  xxxxxxxxxxXXXXXx
                    // read3:          xxXXXXXxxxx
                    // read5: xxxxxXXXXXXXXXXXxxxxx
                    //             ^^^^^^|||||
                    //
                    // Here, we would argue that all of the reads (being single-exon) match
                    // at least one CDS coordinate + all of them are within TOGA boundaries
                    //
                    // For a clearer illustration go to mm39 chr3:61,269,596-61,278,844
                    let is_match = read
                        .exons()
                        .iter()
                        .any(|exon| reference_exons.contains(exon));

                    if is_match {
                        log::debug!(
                            "INFO: read: {:?} single-exon in a multi-read component overlaps with TOGA single-exon CDS -> keep!",
                            &read,
                        );
                        keep.push(line);
                    } else {
                        log::debug!( "INFO: read: {:?} single-exon in a multi-read component does not overlap with TOGA single-exon CDS -> orphan!", &read);
                        discards += 1;
                    }
                } else {
                    counter.inc_read_se_mc_toga_me();

                    // INFO: main question: your CDS matches exactly or partially a TOGA CDS?
                    // CASE: single-exon fuzzy overlap with TOGA CDS
                    // INFO: the following case is presented:
                    //
                    // toga:     xxxxxxXXXXXXX----------
                    //               ^^||||^^^
                    // read1: xxxxxxxXXXXXXXxxxx
                    // read2:       xxxXXXXXXX----------
                    // read3:     xxxxxXXXXXXX----------
                    //
                    // Here, read1 overlaps at the CDS level and thus is part of the component;
                    // however, is not an exact CDS match.
                    //
                    // For a clearer illustration go to mm39 chr3:65,014,198-65,019,415
                    // or mm39 chr3:117,366,942-117,374,421
                    let is_match = read
                        .exons()
                        .iter()
                        .any(|exon| reference_exons.contains(exon));

                    if is_match {
                        log::debug!(
                            "INFO: read: {:?} single-exon in a multi-read component overlaps with TOGA multi-exon CDS -> keep!",
                            &read,
                        );
                        keep.push(read.to_bed::<Bed12>());
                    } else {
                        log::debug!("INFO: read: {:?} single-exon in a multi-read component does not overlap with TOGA multi-exon CDS -> orphan!", &read);
                        discards += 1;
                    }
                }
            } else {
                if is_reference_single_exon {
                    counter.inc_read_me_mc_toga_se();
                } else {
                    counter.inc_read_me_mc_toga_me();
                }

                // INFO: means CDS TOGA overlap and not single-exon
                // INFO: branch where most of the cases will be in
                // INFO: could present the following case:
                //
                // toga:     xxxxxxXXXXXXX----XXXXX----
                //               ^^||||^^^
                // read1: xxxxxxxXXXXXXXXX----XXXXXxx -> shorter isoform
                // read2:       xxxXXXXXXX----XXXXX----
                // read3:     xxxxxXXXXXXX----XXXXX----
                //
                // or its variants:
                //
                // toga:     xxxxxxXXXXXXX----------
                //               ^^||||^^^
                // read1: xxxxxxxXXXXXXXXX----xxxxx
                // read2:       xxxXXXXXXX----------
                // read3:     xxxxxXXXXXXX----------
                //
                // or:
                //
                // toga:  xxxxXXXXX-------XXXXXXX----
                //                               ^^
                // read1: xxxxXXXXXxxxxxxxXXXXXXXXXxx
                // read2: xxxxXXXXX-------XXXXXXX----
                //
                // Here, we will evaluate if the present read has more than 1
                // exact CDS match.
                //
                // For a clearer illustration go to mm39 chr3:65,435,178-65,440,499

                // WARN: not making is_reference_single_exon check because it is not needed
                let mut matches = 0;
                for read_exon in read.exons().iter() {
                    if reference_exons.contains(read_exon) {
                        matches += 1;
                    }

                    if matches > 1 {
                        log::debug!(
                            "DEBUG: read {:?} multi-exon in a mult-read component has more than 1 exact CDS matches -> keep!",
                            &read
                        );

                        keep.push(line);
                        break;
                    }
                }

                if matches <= 1 {
                    discards += 1;

                    log::debug!(
                        "DEBUG: read {:?} multi-exon in a mult-read component has less than 1 exact CDS matches -> increase discards!",
                        &read
                    );
                    discards += 1;
                }
            }
        }

        // INFO: if we are discarding more than 50% of reads, perform splice site matching
        if (discards as f32) / queries.len() as f32 >= min_discard_percentage {
            counter.inc_component_above_discards();

            let reference_exons_flat: BTreeSet<u64> =
                reference_exons.iter().flat_map(|(s, e)| [*s, *e]).collect();

            for read in queries.iter() {
                let is_splice_match = read.exons().iter().any(|(exon_start, exon_end)| {
                    reference_exons_flat.contains(exon_start)
                        || reference_exons_flat.contains(exon_end)
                });

                let line = read.to_bed::<Bed12>();

                if is_splice_match {
                    log::debug!(
                        "DEBUG: read {:?} from any category that was being discarded has at least 1 splice site match with reference -> rescue!",
                        &read
                    );

                    counter.inc_rescue();
                    if !keep.contains(&line) {
                        keep.push(line);
                    }
                } else {
                    log::debug!(
                        "DEBUG: read {:?} multi-exon in a mult-read component has no splice site matches with reference -> orphan confirmed!",
                        &read
                    );
                    counter.inc_read_no_splice_match();

                    if !orphans.contains(&line) {
                        orphans.push(line);
                    }
                }
            }
        } else {
            // INFO: likely real drop
            for read in queries.iter() {
                let line = read.to_bed::<Bed12>();
                if !orphans.contains(&line) && !keep.contains(&line) {
                    orphans.push(line);
                    counter.inc_read_no_splice_match();
                }
            }
        }
    } else {
        if queries.is_empty() {
            // INFO: TOGA projection without reads
            let projections: Vec<Option<&[u8]>> = references.iter().map(|p| p.name()).collect();
            log::trace!("DEBUG: TOGA projection without reads -> {:?}", projections);
        }

        // INFO: weird case of TOGA-single-exon and single-read component
        if is_reference_single_exon {
            for read in queries {
                let is_single_exon_read = read.introns().is_empty();
                let line = read.to_bed::<Bed12>();

                if is_single_exon_read {
                    counter.inc_read_se_sc_toga_se();

                    // INFO: CDS must matche once
                    // INFO: branch -> CDS-overlap with TOGA but single-read component
                    //
                    // ideal case:
                    //
                    // toga: xxxxXXXXXXxxxx
                    //           ||||||
                    // read1:  xxXXXXXXxxxx
                    //
                    // its counterpart:
                    //
                    // toga: xxxxXXXXXXxxxx
                    //         ^^|||
                    // read1: xXXXXXxxx <- likely not a supporting transcript
                    let is_match = read
                        .exons()
                        .iter()
                        .any(|exon| reference_exons.contains(exon));

                    if is_match {
                        log::debug!(
                            "DEBUG: read: {:?} single-exon and match exactly with TOGA single-exon -> keep!",
                            &read,
                        );

                        keep.push(line);
                    } else {
                        log::debug!(
                            "DEBUG: read: {:?} single-exon and does not match exactly with TOGA single-exon -> orphan!",
                            &read,
                        );
                        orphans.push(line);
                    }
                } else {
                    counter.inc_read_me_sc_toga_se();

                    // INFO: TOGA-single-exon CDS overlap with multi-exon single-read component
                    // INFO: at least one splice site should match
                    //
                    // ideal case:
                    //
                    // toga: xxxxXXXXXXxxxx
                    //           ||||||
                    // read1:  xxXXXXXXxxxx-----xxxxxx
                    //
                    // its counterpart:
                    //
                    // toga: xxxxXXXXXXxxxx
                    //         ^^|||
                    // read1: xXXXXXxxx-------xxxxx <- likely not a supporting transcript

                    let reference_exons_flat: BTreeSet<u64> =
                        reference_exons.iter().flat_map(|(s, e)| [*s, *e]).collect();

                    let is_splice_match = read.exons().iter().any(|(exon_start, exon_end)| {
                        reference_exons_flat.contains(exon_start)
                            || reference_exons_flat.contains(exon_end)
                    });

                    if is_splice_match {
                        log::debug!(
                            "DEBUG: read: {:?} multi-exon matches at least 1 splice site with TOGA single-exon -> keep!",
                            &read,
                        );

                        keep.push(line);
                    } else {
                        log::debug!(
                            "DEBUG: read: {:?} multi-exon does not match any splice site with TOGA single-exon -> orphan!",
                            &read,
                        );

                        orphans.push(line);
                    }
                }
            }
        } else {
            // INFO: TOGA multi-exon, single-read component
            // INFO: ask if CDS matches are more than 1 OR is within TOGA boundaries

            for read in queries {
                let is_single_exon_read = read.introns().is_empty();
                let line = read.to_bed::<Bed12>();

                if is_single_exon_read {
                    counter.inc_read_se_sc_toga_me();

                    // CASE: single-exon read component with multi-exon TOGA refs
                    // INFO: the following case is presented:
                    //
                    // toga:  xxxXXXX----XXXXX----XXXxxx
                    //           ||||
                    // read1: xxxXXXX
                    //
                    // or any of its  variants:
                    //
                    // toga:  xxxXXXX----XXXXX----XXXxxx
                    // read1: xxxxxxxxxxxxxxXXXXxxxx
                    //
                    // Here, we ask for specific exon match otherwise would be hard
                    // to distinguish with a single-exon read component
                    let is_match = read
                        .exons()
                        .iter()
                        .any(|exon| reference_exons.contains(exon));

                    if is_match {
                        log::debug!(
                            "DEBUG: read: {:?} single-exon in single-read component and match exactly with TOGA multi-exon -> keep!",
                            &read,
                        );
                        keep.push(line);
                    } else {
                        log::debug!(
                            "DEBUG: read: {:?} single-exon in single-read component and does not match exactly with TOGA multi-exon -> orphan!",
                            &read,
                        );
                        orphans.push(line);
                    }
                } else {
                    counter.inc_read_me_sc_toga_me();

                    // CASE: multi-exon single-read component with multi-exon TOGA refs
                    //
                    // Here, we ask for at least one specific exon match otherwise would be hard
                    // to distinguish with a single-read component
                    let is_match = read
                        .exons()
                        .iter()
                        .any(|exon| reference_exons.contains(exon));

                    if is_match {
                        log::debug!(
                                "DEBUG: read: {:?} multi-exon in single-read component and match more than once exactly with TOGA multi-exon -> keep!",
                                &read,
                            );

                        keep.push(line);
                    } else {
                        log::debug!(
                                    "DEBUG: read: {:?} multi-exon in single-read component and does not match more than once exactly with TOGA multi-exon -> orphan!",
                                    &read,
                                );
                        orphans.push(line);
                    }
                }
            }
        }
    }

    (keep, orphans)
}

/// Self-guided mode
///
/// # Arguments
///
/// * `components` - The components to process
/// * `counter` - The counter to use
/// * `min_read_num_denovo` - The min_read_num_denovo to use
///
/// # Returns
///
/// A tuple containing the keep and orphans
///
/// # Examples
///
/// ```
/// use isotools::iso_orphan::self_guided;
/// use isotools::utils::ParallelCounter;
///
/// let components = Box::new((vec![], vec![]));
/// let counter = ParallelCounter::default();
/// let min_read_num_denovo = 5;
///
/// let (keep, orphans) = self_guided(components, &counter, min_read_num_denovo);
/// ```
#[allow(clippy::boxed_local)]
fn self_guided(
    mut queries: Vec<GenePred>,
    counter: &ParallelCounter,
    min_read_num_denovo: usize,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    do_self_guided_check(&mut queries, counter, min_read_num_denovo)
}

/// Self-guided mode executor
///
/// # Arguments
///
/// * `reads` - The reads to process
/// * `counter` - The counter to use
/// * `min_read_num_denovo` - The min_read_num_denovo to use
///
/// # Returns
///
/// A tuple containing the keep and orphans
///
/// # Examples
///
/// ```
/// use isotools::iso_orphan::do_self_guided_check;
/// use isotools::utils::ParallelCounter;
///
/// let reads = vec![];
/// let counter = ParallelCounter::default();
/// let min_read_num_denovo = 5;
///
/// let (keep, orphans) = do_self_guided_check(&mut reads, &counter, min_read_num_denovo);
/// ```
fn do_self_guided_check(
    reads: &mut Vec<GenePred>,
    counter: &ParallelCounter,
    min_read_num_denovo: usize,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut keep = Vec::new();
    let mut orphans = Vec::new();

    let is_single_component_read = reads.len() == 1;

    // CASE: single-component single-exon
    // INFO: the following case is presented:
    //
    // [SOLVED: discarding single-component reads]
    // toga:            XXX--XX---XXXX
    // read1: xxXXXXxx  |||  ||   ||||
    // read2:           XXX--XX---XXXX
    //
    // or its variants:
    //
    // [SOLVED: discarding single-component reads]
    // toga:  XX-----------XXX--XX---XXXX
    // read1:    xxXXXXxx  |||  ||   ||||
    // read2: XX-----------XXX--XX---XXXX
    //
    //
    // Here, read1 is a single-exon component read that can be easily
    // discarded because it represents a common case of background noise.
    // For a clearer illustration go to mm39 chr8:71,358,478-71,360,769,
    // or mm39 chr8:72,647,487-72,650,648.
    //
    // toga:                 XXX--XX---XXXX
    // read1: xxXXXX--XXXxx  |||  ||   ||||
    // read2:                XXX--XX---XXXX
    //
    // The above case would illustrate a component where a non-single-exon
    // component read could repesent real transcription but not at a significant
    // level.
    //
    // For a clearer illustration go to mm39 chr8:71,358,478-71,360,769
    if is_single_component_read {
        log::trace!(
            "DEBUG: single-component single-exon read with no TOGA refs -> {:?} -> orphan!",
            reads
        );
        counter.inc_read_se_sc_no_toga();

        for read in reads {
            let line = read.to_bed::<Bed12>();
            orphans.push(line);
        }

        return (keep, orphans);
    }

    // CASE: multi-exon components
    // INFO: the following case is presented:
    //
    // [PARTIALLY SOLVED: min_read_num_denovo]
    // toga:  XX-----------XXX--XX---XXXX
    // read1:    xxXXXXxx
    // read2:    xxxxXXXxx
    // read3: XX-----------XXX--XX---XXXX
    //
    // Here, read1 and read2 overlap and form a unique non-TOGA component
    // with more than 1 read. Even though a min_read_num_denovo would be too relative
    // to label a component as background transcription, we force any isolated
    // component to have more than 5 reads.
    //
    // For a clearer illustration go to mm39 chr3:97,596,920-97,624,824
    // or mm39 chr3:103,646,891-103,652,710
    if reads.len() < min_read_num_denovo {
        log::debug!(
            "DEBUG: component with less than {} reads min_read_num_denovo -> {:?} -> orphan!",
            min_read_num_denovo,
            reads
        );
        counter.inc_component_less_than_threshold();

        for read in reads {
            let line = read.to_bed::<Bed12>();
            orphans.push(line);
        }

        return (keep, orphans);
    }

    // INFO: sorting reads by absolute exonic_len to allow group single-exons
    // INFO: while collection information from bigger reads
    reads.sort_by(|a, b| a.exonic_length().cmp(&b.exonic_length()));

    // INFO: establishing exonic matches from reads
    // INFO: { exon -> matches }, then rank by matches and keep if reads with exon > 50% ocurrence
    //
    // CASE: single-exon reads following same/diff exonic pattern
    // INFO: the following case is presented:
    //
    // read1:       xxxXXXXXXXxxx
    // read2:      xxxxxxxxxxXXXXXxxxx
    //
    // Here, both reads would belong to the same component because exonic
    // overlap (accounting for UTRs); however, their CDS structure is not consistent.
    //
    // read1:       xxxXXXXXXXxxx
    // read2:      xxxxXXXXXXXxxxxxxx
    //
    // In contrast, the above case would illustrate a component where single-exon
    // reads follow the same exonic pattern. This likely repesents uniform transcription.
    //
    // For a clearer illustration go to HLpteAle1A HAP1_SUPER_1:112,960,536-112,964,453
    //
    // CASE: single-exon read partially overlap CDS
    // INFO: the following case is presented:
    //
    // read1: xxxxxxxxxXXXXxxx
    // read2: -------xxxxXXXXX------XXX--
    // read3: ---------xxXXXXX------XXX--
    // read4: ---------xxxxxxxxxXX--XXX--
    // read5: ----------xXXXXX-----------
    //                   |||||      |||
    //
    // The above case would illustrate a component where a single-exon read
    // CDS overlaps partially with a non-single-exon read CDS. We would argue
    // that this read - even if it matches at the CDS level - does not follow
    // the same transcription pattern as the other reads.
    //
    // For a clearer illustration go to mm39 chr8:72,042,766-72,051,539
    // or mm39 chr8:73,174,008-73,179,856
    let mut rank = HashMap::new();
    for read in reads.iter() {
        for exon in read.exons().iter() {
            rank.entry(exon.clone())
                .and_modify(|e| *e += 1)
                .or_insert(1);
        }
    }

    let unique_exons = rank.len();
    for (exon, matches) in rank.iter() {
        log::debug!("DEBUG: exon: {:?} matches: {}", exon, matches);

        if (*matches as f32 / unique_exons as f32) > 0.5 {
            log::debug!(
                "DEBUG: supported exon: {:?} with matches: {} in set of reads: {:?}",
                exon,
                matches,
                reads
            );

            for read in reads.iter() {
                let line = read.to_bed::<Bed12>();
                if read.exons().contains(exon) {
                    log::debug!("DEBUG: read: {:?} has exon: {:?} -> keep!", read, exon);

                    if !keep.contains(&line) {
                        keep.push(line);
                    }
                } else {
                    log::debug!(
                        "DEBUG: read: {:?} does not have exon: {:?} -> orphan!",
                        read,
                        exon
                    );

                    if !orphans.contains(&line) {
                        counter.inc_read_de_novo_unsupported();
                        orphans.push(line);
                    }
                }
            }
        } else {
            log::debug!(
                "DEBUG: non-supported exon: {:?} with matches: {} in set of reads: {:?}",
                exon,
                matches,
                reads
            );
        }
    }

    (keep, orphans)
}
