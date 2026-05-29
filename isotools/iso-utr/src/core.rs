// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for detecting intron retentions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2026
//!
//! This module contains the main algorithm for detecting truncations
//! in a query set of reads.
//!
//! In short it takes a set of query reads and a set of reference reads and
//! detects truncations in the query reads. It does this by checking if the
//! query reads overlap with any middle exon from the reference set of reads.
//! A recovery step can be performed by evaluating the support of the middle
//! exon in the reference set of reads.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufWriter, Write},
};

use genepred::{GenePred, Strand};
use hashbrown::HashMap;
use log::{info, warn};
use packbed::{OverlapType, Role};
use rayon::prelude::*;

use crate::cli::Args;
use crate::utils::*;

pub const SCALE: u64 = 100000000000; // 100Gb

pub fn detect_truncations(mut args: Args) -> Result<(), Box<dyn std::error::Error>> {
    info!("Detecting 5'end truncations...");

    let mut modes = std::iter::repeat(Role::Reference)
        .take(args.refs.len())
        .collect::<Vec<_>>();
    modes.extend(vec![Role::Query]);
    args.refs.extend(args.query);

    log::info!("INFO: Packing BED12 files...");

    let tracks = packbed::pack(args.refs, modes, OverlapType::Exon).unwrap_or_else(|e| {
        log::error!("Failed to pack beds: {:?}", e);
        std::process::exit(1);
    });

    let accumulator = ParallelAccumulator::default();
    let counter = ParallelCounter::default();

    log::info!("INFO: Processing components...");
    tracks.into_par_iter().for_each(|bucket| {
        let components = bucket.1;
        counter.inc_components(components.len() as u32);

        process_components(
            components,
            &accumulator,
            &counter,
            args.recover,
            args.recovery_threshold,
            args.exon_recovery_threshold,
        );
    });

    log::info!("INFO: Processed {} components", counter.load_components());

    if args.recover {
        let (count, ratio) = counter.get_stat();
        warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            count, ratio
        );
    }

    let file = args.outdir.join(format!("{}.truncation.tsv", args.prefix));
    let mut writer = BufWriter::new(File::create(file).unwrap());

    accumulator.lines.into_iter().for_each(|line| {
        writer
            .write_all(&line)
            .unwrap_or_else(|e| panic!("ERROR: Could not write schema -> {e}!"));
    });

    Ok(())
}

/// Processes a component of reads and detects truncations
pub fn process_components(
    components: Vec<Vec<GenePred>>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
    recover: bool,
    recovery_threshold: f32,
    exon_recovery_threshold: f32,
) {
    components.into_par_iter().for_each(|component| {
        counter.inc_components(1);
        let (refs, queries) = component.into_iter().partition(|record| {
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
            recover,
            recovery_threshold,
            exon_recovery_threshold,
            counter,
        );

        result.into_iter().for_each(|descriptor| {
            if !descriptor.is_empty() {
                accumulator.lines.insert(descriptor);
            }
        });
    });
}

/// Returns the first and last exon of a record.
///
/// # Arguments
///
/// * `record` - The record to get the exons from.
///
/// # Returns
///
/// A tuple containing the first and last exon of the record.
///
/// # Example
///
/// ```ignore
/// let record = GenePred::new(b"chr1", 100, 200, b"+", b"exon1", b"exon2");
/// let (first, last) = terminal_exon(&record);
/// assert_eq!(first, 100);
/// assert_eq!(last, 200);
/// ```
#[inline(always)]
fn terminal_exon(record: &GenePred) -> (u64, u64) {
    let exons = record.exons();

    match record.strand() {
        Some(Strand::Forward) => exons.first().copied().unwrap_or_else(|| {
            log::error!(
                "ERROR: Could not get first exon from record -> {:?}",
                record.name()
            );
            std::process::exit(1);
        }),
        Some(Strand::Reverse) => {
            let mut terminal_exon = exons.last().copied().unwrap_or_else(|| {
                log::error!(
                    "ERROR: Could not get last exon from record -> {:?}",
                    record.name()
                );
                std::process::exit(1);
            });
            terminal_exon = (SCALE - terminal_exon.1, SCALE - terminal_exon.0);
            terminal_exon
        }
        Some(Strand::Unknown) | None => {
            log::error!(
                "ERROR: Could not get strand from record -> {:?}",
                record.name()
            );
            std::process::exit(1);
        }
    }
}

/// Processes a component of reads and detects truncations
///
/// # Arguments
///
/// * `refs` - A `Vec<GenePred>` to process
/// * `queries` - A `Vec<GenePred>` to process
/// * `recover` - A `bool` to recover reads from
/// * `recovery_threshold` - A `f32` to recover reads from
/// * `exon_recovery_threshold` - A `f32` to recover reads from
///
/// # Example
///
/// ```ignore
/// process_component(
///     refs,
///     queries,
///     recover,
///     recovery_threshold,
///     exon_recovery_threshold,
/// );
/// ```
#[inline(always)]
pub fn process_component(
    refs: Vec<GenePred>,
    queries: Vec<GenePred>,
    recover: bool,
    recovery_threshold: f32,
    exon_recovery_threshold: f32,
    counter: &ParallelCounter,
) -> Vec<Vec<u8>> {
    let mut descriptor: HashMap<Option<&[u8]>, Schema<'_>> = HashMap::new();
    let mut accumulator = Vec::new();

    let mut reference_starts: BTreeSet<(u64, u64)> = BTreeSet::new();
    let mut reference_middle_exons: BTreeSet<(u64, u64)> = BTreeSet::new();
    let mut reference_middle_exons_support: BTreeMap<(u64, u64), f32> = BTreeMap::new();

    log::debug!("DEBUG: Getting reference exons...");
    refs.iter().for_each(|record| match record.strand() {
        Some(Strand::Forward) => {
            let exons = record.exons();
            reference_starts.insert(exons.first().copied().unwrap_or_else(|| {
                log::error!(
                    "ERROR: Could not get first exon from record -> {:?}",
                    record.name()
                );
                std::process::exit(1);
            }));

            exons[1..].into_iter().for_each(|exon| {
                reference_middle_exons.insert(*exon);
                reference_middle_exons_support
                    .entry(*exon)
                    .and_modify(|x| *x += 1.0)
                    .or_insert(1.0);
            });
        }
        Some(Strand::Reverse) => {
            let exons = record.exons();

            let mut terminal_exon = exons.last().copied().unwrap_or_else(|| {
                log::error!(
                    "ERROR: Could not get last exon from record -> {:?}",
                    record.name()
                );
                std::process::exit(1);
            });
            terminal_exon = (SCALE - terminal_exon.1, SCALE - terminal_exon.0);
            reference_starts.insert(terminal_exon);

            exons[..exons.len() - 1].into_iter().for_each(|exon| {
                let exon = (SCALE - exon.1, SCALE - exon.0);

                reference_middle_exons.insert(exon);
                reference_middle_exons_support
                    .entry(exon)
                    .and_modify(|x| *x += 1.0)
                    .or_insert(1.0);
            });
        }
        Some(Strand::Unknown) | None => {
            log::error!(
                "ERROR: Could not get strand from record -> {:?}",
                record.name()
            );
            std::process::exit(1);
        }
    });
    log::debug!("DEBUG: Reference exons -> {:?}", reference_starts);
    log::debug!("DEBUG: Middle exons -> {:?}", reference_middle_exons);

    let (mut truncations, totals) = (0_f32, queries.len() as f32);

    for query in queries.iter() {
        log::debug!("DEBUG: Processing query -> {:?}", query.name());
        let mut schema = Schema::default();
        let query_name = query.name().unwrap_or_else(|| {
            log::error!("ERROR: Could not get name from record -> {:?}", query);
            std::process::exit(1);
        });
        schema.id = query_name;

        let (query_first_exon_start, query_first_exon_end) = terminal_exon(query);
        log::debug!(
            "DEBUG: Query exons -> {:?}",
            (query_first_exon_start, query_first_exon_end)
        );

        let is_complete = reference_starts.iter().any(|&(s, e)| {
            if query_first_exon_end < s {
                schema.is_novel_start = b"NOVEL_START";
                return false;
            }

            (query_first_exon_start >= s) && (query_first_exon_start < e)
        });

        if is_complete {
            log::debug!("DEBUG: Query is complete -> {:?}", query.name());
            if schema.is_novel_start != b"NOVEL_START" {
                schema.is_novel_start = b"NOT_NOVEL_START";
            }

            // still checks if read start is inside any middle boundaries
            log::debug!(
                "DEBUG: Checking middle exons -> {:?}",
                reference_middle_exons
            );
            reference_middle_exons
                .iter()
                .for_each(|&(mid_exon_start, mid_exon_end)| {
                    if query_first_exon_end < mid_exon_start {
                        log::debug!("DEBUG: No overlap with middle exon -> {:?}", query.name());
                        return;
                    }

                    if (query_first_exon_start >= mid_exon_start)
                        && (query_first_exon_start < mid_exon_end)
                    {
                        log::debug!("DEBUG: Truncation found -> {:?}", query.name());
                        schema.status = b"TRUNCATED";
                        schema.code = b"T";

                        match query.strand() {
                            Some(Strand::Forward) => {
                                schema.exon = Vec::new();
                                schema.exon.extend(query.chrom());
                                schema.exon.push(b':');
                                schema.exon.extend(mid_exon_start.to_string().as_bytes());
                                schema.exon.push(b'-');
                                schema.exon.extend(mid_exon_end.to_string().as_bytes());
                            }
                            Some(Strand::Reverse) => {
                                schema.exon = Vec::new();
                                schema.exon.extend(query.chrom());
                                schema.exon.push(b':');
                                schema
                                    .exon
                                    .extend((SCALE - mid_exon_end).to_string().as_bytes());
                                schema.exon.push(b'-');
                                schema
                                    .exon
                                    .extend((SCALE - mid_exon_start).to_string().as_bytes());
                            }
                            _ => {
                                log::error!(
                                    "ERROR: Could not get strand from record -> {:?}",
                                    query.name()
                                );
                                std::process::exit(1);
                            }
                        }

                        truncations += 1.0;
                        return;
                    } else {
                        if schema.code == b"T" {
                            log::debug!("DEBUG: Read already marked as truncated, overlapping another middle exon -> {:?}", query.name());
                            return;
                        }

                        log::debug!("DEBUG: Pass found -> {:?}", query.name());
                        schema.status = b"PASS";
                        schema.code = b"A";
                        return;
                    }
                });

            // INFO: doest not overlap any middle exon
            if schema.status.is_empty() {
                log::debug!(
                    "DEBUG: No overlap with ANY middle exon -> {:?}",
                    query.name()
                );
                schema.status = b"PASS";
                schema.code = b"A";
            };
        } else {
            // we do not have any overlap with consensus starts.
            // we need to see if we overlap any middle exons, if so read
            // is truncated, otherwise it is a novel start.
            // let is_truncated =
            reference_middle_exons
                .iter()
                .for_each(|&(mid_exon_start, mid_exon_end)| {
                    if query_first_exon_end < mid_exon_start {
                        return;
                    }

                    if (query_first_exon_start >= mid_exon_start
                        && query_first_exon_start < mid_exon_end)
                        || (query_first_exon_end > mid_exon_start
                            && query_first_exon_end <= mid_exon_end)
                        || (query_first_exon_start < mid_exon_start
                            && query_first_exon_end > mid_exon_end)
                    {
                        schema.status = b"TRUNCATED";
                        schema.code = b"T";

                        match query.strand() {
                            Some(Strand::Forward) => {
                                schema.exon = Vec::new();
                                schema.exon.extend(query.chrom());
                                schema.exon.push(b':');
                                schema.exon.extend(mid_exon_start.to_string().as_bytes());
                                schema.exon.push(b'-');
                                schema.exon.extend(mid_exon_end.to_string().as_bytes());
                            }
                            Some(Strand::Reverse) => {
                                schema.exon = Vec::new();
                                schema.exon.extend(query.chrom());
                                schema.exon.push(b':');
                                schema
                                    .exon
                                    .extend((SCALE - mid_exon_end).to_string().as_bytes());
                                schema.exon.push(b'-');
                                schema
                                    .exon
                                    .extend((SCALE - mid_exon_start).to_string().as_bytes());
                            }
                            _ => {
                                log::error!(
                                    "ERROR: Could not get strand from record -> {:?}",
                                    query.name()
                                );
                                std::process::exit(1);
                            }
                        }

                        truncations += 1.0;
                    } else {
                        schema.status = b"PASS";
                        schema.code = b"A";
                        schema.is_novel_start = b"NOVEL_START";
                    }
                });

            // INFO: doest not overlap any middle exon
            if schema.status.is_empty() {
                log::debug!(
                    "DEBUG: No overlap with ANY middle exon -> {:?}",
                    query.name()
                );
                schema.status = b"PASS";
                schema.code = b"A";
            };
        }

        descriptor.insert(Some(query_name), schema);
    }

    // INFO: after classying reads, we check bucket frequencies
    // if the number of truncated reads is greater than 50%
    // of the total reads in the bucket, we consider the bucket
    // to be dirty and return the reads for recovery if args.recover
    let ratio = truncations / totals;
    if recover {
        if ratio >= recovery_threshold {
            counter.inc_dirty();
            log::debug!("Bucket {:?} is dirty -> {}", queries, ratio);
            reference_middle_exons_support
                .iter_mut()
                .for_each(|(_k, v)| {
                    *v = *v as f32 / refs.len() as f32;
                });

            recover_reads(
                &queries,
                reference_middle_exons_support,
                &mut descriptor,
                exon_recovery_threshold,
                ratio,
                &mut accumulator,
            );
        } else {
            push_descriptors(&queries, &mut descriptor, ratio, &mut accumulator);
        }
    } else {
        push_descriptors(&queries, &mut descriptor, ratio, &mut accumulator);
    }

    accumulator
}

fn push_descriptors<'a>(
    queries: &'a [GenePred],
    descriptor: &mut HashMap<Option<&'a [u8]>, Schema<'a>>,
    ratio: f32,
    accumulator: &mut Vec<Vec<u8>>,
) {
    for query in queries.iter() {
        let query_key = Some(query.name().unwrap_or_else(|| {
            log::error!("ERROR: Could not get name from record -> {:?}", query);
            std::process::exit(1);
        }));

        let schema = descriptor.get_mut(&query_key).unwrap_or_else(|| {
            log::error!("ERROR: Could not get handle for query {:?}", query.name());
            std::process::exit(1);
        });

        schema.ratio = ratio;
        accumulator.push(schema.to_line());
    }
}

/// Recovers reads from the Bucket
///
/// # Arguments
///
/// * `queries` - A slice of `GenePred` to recover reads from
/// * `reference_exon_support` - A `BTreeMap` of `(u64, u64)` to `f32` to recover reads from
/// * `descriptor` - A `HashMap` of `Option<&[u8]>` to `Schema`
/// * `recovery_threshold` - A `f32` to recover reads from
/// * `ratio` - A `f32` to recover reads from
/// * `accumulator` - A `Vec<Vec<u8>>` to recover reads from
///
/// # Example
///
/// ```ignore
/// recover_reads(
///     &queries,
///     reference_exon_support,
///     &mut descriptor,
///     exon_recovery_threshold,
///     ratio,
///     &mut accumulator,
/// );
/// ```
fn recover_reads<'a>(
    queries: &'a [GenePred],
    reference_exon_support: BTreeMap<(u64, u64), f32>,
    descriptor: &mut HashMap<Option<&'a [u8]>, Schema<'a>>,
    recovery_threshold: f32,
    ratio: f32,
    accumulator: &mut Vec<Vec<u8>>,
) {
    for (ref_exon, support) in reference_exon_support.iter() {
        // 2. if the owner (middle exon) has more than 50% support, we consider it
        //  a valid owner and keep reads truncated by that owner as
        //  truncated reads; otherwise, we consider the owner to be
        //  a weak ownner and send all truncated reads to the pass
        //  bucket

        let (owner_start, owner_end) = *ref_exon;

        // send reads truncated by this owner to pass
        for query in queries.iter() {
            let query_key = Some(query.name().unwrap_or_else(|| {
                log::error!("ERROR: Could not get name from record -> {:?}", query);
                std::process::exit(1);
            }));

            let schema = descriptor.get_mut(&query_key).unwrap_or_else(|| {
                log::error!("ERROR: Could not get handle for query {:?}", query.name());
                std::process::exit(1);
            });

            schema.ratio = ratio;

            // INFO: if read was already a pass, leave it as pass
            if schema.status == b"PASS" {
                // INFO: marker for reviewed reads
                schema.code = b"AW";
                accumulator.push(schema.to_line());
                continue;
            }

            // INFO: from this point we only have truncated reads to evaluate
            let (query_start, query_end) = terminal_exon(query);

            if (query_start >= owner_start && query_start < owner_end)
                || (query_end > owner_start && query_end <= owner_end)
                || (query_start < owner_start && query_end > owner_end)
            {
                // INFO: read is truncated by this owner and owner is NOT supported
                // INFO: modifying read as a forced pass
                if *support < recovery_threshold {
                    log::debug!(
                        "DEBUG: Owner {:?} has less than threshold {:?}",
                        ref_exon,
                        support
                    );

                    schema.status = b"PASS";
                    schema.code = b"FW";
                } else {
                    // INFO: read is truncated by this owner and owner is supported
                    log::debug!(
                        "DEBUG: Owner {:?} has more than threshold {:?}",
                        ref_exon,
                        support
                    );
                    schema.status = b"TRUNCATED";
                    schema.code = b"TW";
                }
            }

            accumulator.push(schema.to_line());
        }
    }
}

/// Schema
#[derive(Debug, Default, Clone)]
struct Schema<'a> {
    id: &'a [u8],
    status: &'a [u8],
    code: &'a [u8],
    ratio: f32,
    is_novel_start: &'a [u8],
    exon: Vec<u8>,
}

impl Schema<'_> {
    /// Formats the schema into a line of TSV
    ///
    /// # Returns
    ///
    /// A `Vec<u8>` containing the formatted line of TSV
    pub fn to_line(&self) -> Vec<u8> {
        let mut body = String::new();

        body.push_str("<h2>Truncation</h2><br>");
        body.push_str(&format!(
            "- status: {}<br>",
            std::str::from_utf8(&self.status).unwrap_or("NULL")
        ));
        body.push_str(&format!(
            "- code: {} [T: truncation, A: no events, F: forced, W: reviewed]<br>",
            std::str::from_utf8(&self.code).unwrap_or("NULL")
        ));
        body.push_str(&format!(
            "- support: {}<br>",
            std::str::from_utf8(&self.is_novel_start).unwrap_or("NULL")
        ));
        body.push_str(&format!("- ratio: {}<br>", self.ratio));
        body.push_str(&format!(
            "- exon: {}<br>",
            std::str::from_utf8(&self.exon).unwrap_or("NULL")
        ));

        // fmt: id\tstatus\tcode\thtml  — all on one line
        let mut out = Vec::new();
        out.extend_from_slice(&self.id);
        out.push(b'\t');
        out.extend_from_slice(&self.status);
        out.push(b'\t');
        out.extend_from_slice(&self.code);
        out.push(b'\t');
        out.extend_from_slice(body.as_bytes());
        out.push(b'\n'); // single trailing newline for the TSV row

        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile;

    #[test]
    fn test_detection_fn_with_tempfile() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        let _path = file.path().to_path_buf();
        write!(
            file,
            "s1/t778845/t804002/tm54164U_210309_085211/65275776/ccs_PerID0.996_5Clip0_3Clip0_PolyA274_PolyARead275/t60/t-/t778845/t804002/t255,0,0/t14/t1268,142,88,200,253,203,254,167,120,142,218,197,132,302/t0,14396,14736,14943,16461,16846,17598,18073,18511,18890,20330,21290,22627,24855"
        ).unwrap();

        write!(
            file,
            "s1/t778870/t803968/tm54164U_210309_085211/92276372/ccs_PerID1.000_5Clip0_3Clip0_PolyA71_PolyARead72/t60/t-/t778870/t803968/t255,0,0/t11/t1243,142,88,200,253,1006,167,218,197,132,268/t0,14371,14711,14918,16436,16821,18048,20305,21265,22602,24830"
        ).unwrap();
    }

    #[test]
    fn recover_mode_emits_non_dirty_component_descriptors() {
        let reference = record(b"ref", &[(100, 150), (200, 250)]);
        let query_a = record(b"query_a", &[(110, 145), (200, 250)]);
        let query_b = record(b"query_b", &[(120, 140), (205, 240)]);
        let counter = ParallelCounter::default();

        let lines = process_component(
            vec![reference],
            vec![query_a, query_b],
            true,
            0.5,
            0.5,
            &counter,
        );

        assert_eq!(lines.len(), 2);

        let lines = lines
            .into_iter()
            .map(|line| String::from_utf8(line).unwrap())
            .collect::<Vec<_>>();
        assert!(lines
            .iter()
            .any(|line| line.starts_with("query_a\tPASS\tA\t")));
        assert!(lines
            .iter()
            .any(|line| line.starts_with("query_b\tPASS\tA\t")));
    }

    fn record(name: &[u8], exons: &[(u64, u64)]) -> GenePred {
        let mut record = GenePred::from_coords(
            b"chr1".to_vec(),
            exons.first().unwrap().0,
            exons.last().unwrap().1,
            Default::default(),
        );
        record.set_name(Some(name.to_vec()));
        record.set_strand(Some(Strand::Forward));
        record.set_block_count(Some(exons.len() as u32));
        record.set_block_starts(Some(exons.iter().map(|(start, _)| *start).collect()));
        record.set_block_ends(Some(exons.iter().map(|(_, end)| *end).collect()));
        record
    }
}
