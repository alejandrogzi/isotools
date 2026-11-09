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
use log::{info, warn};
use packbed::Role;
use rayon::prelude::*;

use crate::cli::Args;
use crate::utils::*;

pub const SCALE: u64 = 100000000000; // 100Gb

/// Donor of a middle exon k -> (end of exon k-1, exon k), in transcript orientation
type MiddleDonors = BTreeMap<u64, Vec<(u64, (u64, u64))>>;

pub fn detect_truncations(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    // INFO: main() built a pool and dropped it, so --threads was ignored
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?
        .install(|| truncations(args))
        .map_err(|e| -> Box<dyn std::error::Error> { e })
}

fn truncations(mut args: Args) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    info!("Detecting 5'end truncations...");

    let mut modes = std::iter::repeat(Role::Reference)
        .take(args.refs.len())
        .collect::<Vec<_>>();
    modes.extend(vec![Role::Query]);
    args.refs.extend(args.query);

    log::info!("INFO: Packing BED12 files...");

    let tracks = packbed::pack(args.refs, modes, args.overlap_type.into()).unwrap_or_else(|e| {
        log::error!("Failed to pack beds: {:?}", e);
        std::process::exit(1);
    });

    let counter = ParallelCounter::default();

    log::info!("INFO: Processing components...");
    let mut lines: Vec<Vec<u8>> = tracks
        .into_par_iter()
        .flat_map_iter(|bucket| {
            process_components(
                bucket.1,
                &counter,
                args.recover,
                args.recovery_threshold,
                args.exon_recovery_threshold,
                args.start_junction,
            )
        })
        .collect();

    // INFO: rows were written in thread-completion order [DashSet], so two runs gave the
    // INFO: same rows in a different order; sorted by id, exact duplicates dropped
    lines.par_sort_unstable();
    lines.dedup();

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

    for line in lines.iter() {
        writer
            .write_all(line)
            .unwrap_or_else(|e| panic!("ERROR: Could not write schema -> {e}!"));
    }

    Ok(())
}

/// Processes the components of a chrom:strand bucket and returns their rows
pub fn process_components(
    components: Vec<Vec<GenePred>>,
    counter: &ParallelCounter,
    recover: bool,
    recovery_threshold: f32,
    exon_recovery_threshold: f32,
    start_junction: bool,
) -> Vec<Vec<u8>> {
    components
        .into_par_iter()
        .flat_map_iter(|component| {
            counter.inc_components(1);
            let (refs, queries) = component.into_iter().partition(|record| {
                let role = record
                    .get_extra(b"role")
                    .unwrap_or_else(|| panic!("ERROR: Could not get role from record!"))
                    .clone()
                    .into_scalar();

                role == Some(b"reference".to_vec())
            });

            process_component(
                refs,
                queries,
                recover,
                recovery_threshold,
                exon_recovery_threshold,
                start_junction,
                counter,
            )
            .into_iter()
            .filter(|line| !line.is_empty())
        })
        .collect()
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

/// Returns the coordinates of the first coding exon of a record.
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
/// let (first, last) = terminal_cds_exon(&record);
/// assert_eq!(first, 100);
/// assert_eq!(last, 200);
/// ```
#[inline(always)]
fn terminal_cds_exon(record: &GenePred) -> (u64, u64) {
    let exons = record.coding_exons();

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
/// * `start_junction` - Decide truncations by the query's first junction (`--start-junction`)
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
///     start_junction,
///     counter,
/// );
/// ```
#[inline(always)]
pub fn process_component(
    refs: Vec<GenePred>,
    queries: Vec<GenePred>,
    recover: bool,
    recovery_threshold: f32,
    exon_recovery_threshold: f32,
    start_junction: bool,
    counter: &ParallelCounter,
) -> Vec<Vec<u8>> {
    let mut accumulator = Vec::new();

    // INFO: every coordinate below is in transcript orientation [SCALE-flipped on '-']
    // INFO: reference first exons, with the name of the record that starts there
    let mut reference_starts: Vec<((u64, u64), &[u8])> = Vec::new();
    // INFO: non-first exons -> (record, exon index) owning them, so a truncation can be
    // INFO: checked against the intron chain of the record that contains its start
    let mut containers: BTreeMap<(u64, u64), Vec<(usize, usize)>> = BTreeMap::new();
    let mut reference_introns: Vec<Vec<(u64, u64)>> = Vec::with_capacity(refs.len());
    // INFO: donor of every middle exon k [neither first nor last] -> (end of exon k-1,
    // INFO: exon k), for first exons that start in intron k-1 and share that donor
    let mut middle_donors: MiddleDonors = BTreeMap::new();
    let mut reference_middle_exons_support: BTreeMap<(u64, u64), f32> = BTreeMap::new();

    log::debug!("DEBUG: Getting reference exons...");
    refs.iter().enumerate().for_each(|(idx, record)| {
        if matches!(record.strand(), Some(Strand::Unknown) | None) {
            log::error!(
                "ERROR: Could not get strand from record -> {:?}",
                record.name()
            );
            std::process::exit(1);
        }

        let exons = oriented_exons(record);
        let weight = cn_weight(record) as f32;

        reference_starts.push((
            exons.first().copied().unwrap_or_else(|| {
                log::error!(
                    "ERROR: Could not get first exon from record -> {:?}",
                    record.name()
                );
                std::process::exit(1);
            }),
            record.name().unwrap_or_default(),
        ));

        exons[1..].iter().enumerate().for_each(|(k, exon)| {
            containers.entry(*exon).or_default().push((idx, k + 1));
            reference_middle_exons_support
                .entry(*exon)
                .and_modify(|x| *x += weight)
                .or_insert(weight);
        });

        for k in 1..exons.len().saturating_sub(1) {
            middle_donors
                .entry(exons[k].1)
                .or_default()
                .push((exons[k - 1].1, exons[k]));
        }

        reference_introns.push(oriented_introns(&exons));
    });
    log::debug!("DEBUG: Reference starts -> {:?}", reference_starts);
    log::debug!("DEBUG: Middle exons -> {:?}", containers.keys());

    // INFO: a `#CN<n>` record stands for n identical reads
    let (mut truncations, totals) = (0_f32, queries.iter().map(cn_weight).sum::<usize>() as f32);

    // INFO: first pass -> classify every query on its own; the CDS guard below needs the
    // INFO: full PASS set, so applying it while iterating made results depend on the
    // INFO: (unordered) query order [Q4 5b]
    let mut pass_cds_start_exons = BTreeSet::new();
    let mut classified = Vec::with_capacity(queries.len());

    for query in queries.iter() {
        log::debug!("DEBUG: Processing query -> {:?}", query.name());
        let mut schema = Schema::default();
        let query_name = query.name().unwrap_or_else(|| {
            log::error!("ERROR: Could not get name from record -> {:?}", query);
            std::process::exit(1);
        });
        schema.id = query_name;

        let (query_first_exon_start, query_first_exon_end) = terminal_exon(query);
        let query_cds_start_exon = terminal_cds_exon(query);
        let query_introns = oriented_introns(&oriented_exons(query));

        // INFO: a start is novel when no OTHER record starts there; with --ref == --query
        // INFO: every query is its own reference start, which made this label unreachable
        let is_complete = reference_starts.iter().any(|&((s, e), name)| {
            name != query_name && query_first_exon_start >= s && query_first_exon_start < e
        });
        schema.is_novel_start = if is_complete {
            b"NOT_NOVEL_START"
        } else {
            b"NOVEL_START"
        };

        // INFO: truncated -> the start lies inside a non-first exon of a container AND the
        // INFO: query's intron chain continues the container's chain from that exon on; an
        // INFO: alternative first exon inside another isoform's exon is not a truncation
        // INFO: [Q4 3d: SNTG1]. Smallest container exon first, so the label is deterministic.
        // INFO: --start-junction only asks the query's first intron to be the intron leaving
        // INFO: exon k (a 5' fragment whose chain diverges downstream is still truncated);
        // INFO: an alternative first exon splices elsewhere, so it still passes
        let first_intron = query_introns.first();
        let continues = |introns: &[(u64, u64)], k: usize| -> bool {
            match start_junction {
                true => first_intron.map_or(true, |first| introns.get(k) == Some(first)),
                false => is_sub_chain(&query_introns, introns, k),
            }
        };
        let mut container = containers
            .range(..=(query_first_exon_start, u64::MAX))
            .filter(|((mid_exon_start, mid_exon_end), _)| {
                query_first_exon_start >= *mid_exon_start && query_first_exon_start < *mid_exon_end
            })
            .find_map(|(exon, owners)| {
                owners
                    .iter()
                    .any(|&(idx, k)| continues(&reference_introns[idx], k))
                    .then_some(*exon)
            });

        // INFO: both modes: a first exon that ends at the donor of a middle exon k of another
        // INFO: model and starts inside intron k-1, without reaching exon k-1 [no retention
        // INFO: of intron k-1], is exon k extended into its upstream intron. Only the donor
        // INFO: must match: the query's acceptor and downstream chain may differ. Smallest
        // INFO: such exon, so the label is deterministic
        if container.is_none() && first_intron.is_some() {
            container = middle_donors
                .get(&query_first_exon_end)
                .and_then(|candidates| {
                    candidates
                        .iter()
                        .filter(|&&(upstream_exon_end, (exon_start, _))| {
                            query_first_exon_start >= upstream_exon_end
                                && query_first_exon_start < exon_start
                        })
                        .map(|&(_, exon)| exon)
                        .min()
                });
        }

        match container {
            Some(exon) => {
                log::debug!("DEBUG: Truncation found -> {:?}", query.name());
                schema.exon = exon_label(query, exon);
            }
            None => {
                schema.status = b"PASS";
                schema.code = b"A";
                pass_cds_start_exons.insert(query_cds_start_exon);
            }
        }

        let weight = cn_weight(query) as f32;
        classified.push((
            query_name,
            weight,
            query_cds_start_exon,
            container.is_some(),
            schema,
        ));
    }

    // INFO: second pass -> a truncated query that shares its first coding exon with any
    // INFO: PASS query keeps a complete CDS start, so it passes [order-free guard].
    // INFO: Schemas follow the query order, not a name-keyed map: two records with the
    // INFO: same name used to collapse into whichever one was inserted last
    let mut schemas: Vec<Schema<'_>> = Vec::with_capacity(classified.len());
    for (query_name, weight, query_cds_start_exon, is_truncated, mut schema) in classified {
        if is_truncated {
            if pass_cds_start_exons.contains(&query_cds_start_exon) {
                log::debug!(
                    "DEBUG: Query CDS start exon is in PASS CDS start exons -> {:?}",
                    query_name
                );
                schema.status = b"PASS";
                schema.code = b"A";
                schema.exon = Vec::new();
            } else {
                schema.status = b"TRUNCATED";
                schema.code = b"T";
                truncations += weight;
            }
        }

        schemas.push(schema);
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
            let ref_weight = refs.iter().map(cn_weight).sum::<usize>() as f32;
            reference_middle_exons_support
                .iter_mut()
                .for_each(|(_k, v)| {
                    *v /= ref_weight;
                });

            recover_reads(
                &queries,
                reference_middle_exons_support,
                &mut schemas,
                exon_recovery_threshold,
                ratio,
                &mut accumulator,
            );
        } else {
            push_descriptors(&mut schemas, ratio, &mut accumulator);
        }
    } else {
        push_descriptors(&mut schemas, ratio, &mut accumulator);
    }

    accumulator
}

/// Exons of a record in transcript orientation [SCALE-flipped and reversed on '-']
fn oriented_exons(record: &GenePred) -> Vec<(u64, u64)> {
    let exons = record.exons();
    match record.strand() {
        Some(Strand::Reverse) => exons
            .iter()
            .rev()
            .map(|(start, end)| (SCALE - end, SCALE - start))
            .collect(),
        _ => exons,
    }
}

/// Introns between consecutive oriented exons
fn oriented_introns(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    exons
        .windows(2)
        .map(|pair| (pair[0].1, pair[1].0))
        .collect()
}

/// True if `query` is the container's intron chain starting right after its exon `k`
///
/// A mono-exonic query (empty chain) is a sub-chain of any container.
///
/// # Example
///
/// ```ignore
/// // container introns [(200, 300), (400, 500)], query starts in exon 1
/// assert!(is_sub_chain(&[(400, 500)], &[(200, 300), (400, 500)], 1));
/// assert!(!is_sub_chain(&[(350, 500)], &[(200, 300), (400, 500)], 1));
/// ```
fn is_sub_chain(query: &[(u64, u64)], container: &[(u64, u64)], k: usize) -> bool {
    container.len() >= k + query.len() && container[k..k + query.len()] == *query
}

/// `chrom:start-end` of an oriented exon, in genomic coordinates
fn exon_label(query: &GenePred, (start, end): (u64, u64)) -> Vec<u8> {
    let (start, end) = match query.strand() {
        Some(Strand::Reverse) => (SCALE - end, SCALE - start),
        _ => (start, end),
    };

    let mut label = Vec::new();
    label.extend(query.chrom());
    label.push(b':');
    label.extend(start.to_string().as_bytes());
    label.push(b'-');
    label.extend(end.to_string().as_bytes());
    label
}

/// Sets the component ratio on every schema and pushes its row
fn push_descriptors(schemas: &mut [Schema<'_>], ratio: f32, accumulator: &mut Vec<Vec<u8>>) {
    for schema in schemas.iter_mut() {
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
/// * `schemas` - One `Schema` per query, in query order
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
///     &mut schemas,
///     exon_recovery_threshold,
///     ratio,
///     &mut accumulator,
/// );
/// ```
fn recover_reads(
    queries: &[GenePred],
    reference_exon_support: BTreeMap<(u64, u64), f32>,
    schemas: &mut [Schema<'_>],
    recovery_threshold: f32,
    ratio: f32,
    accumulator: &mut Vec<Vec<u8>>,
) {
    debug_assert_eq!(queries.len(), schemas.len());
    for (query, schema) in queries.iter().zip(schemas.iter_mut()) {
        schema.ratio = ratio;

        // INFO: if read was already a pass, leave it as pass
        if schema.status == b"PASS" {
            // INFO: marker for reviewed reads
            schema.code = b"AW";
            accumulator.push(schema.to_line());
            continue;
        }

        let (query_start, query_end) = terminal_exon(query);
        let mut has_weak_owner = false;
        let mut has_supported_owner = false;

        for (ref_exon, support) in reference_exon_support.iter() {
            let (owner_start, owner_end) = *ref_exon;
            if !overlaps(query_start, query_end, owner_start, owner_end) {
                continue;
            }

            if *support < recovery_threshold {
                log::debug!(
                    "DEBUG: Owner {:?} has less than threshold {:?}",
                    ref_exon,
                    support
                );
                has_weak_owner = true;
            } else {
                log::debug!(
                    "DEBUG: Owner {:?} has more than threshold {:?}",
                    ref_exon,
                    support
                );
                has_supported_owner = true;
                break;
            }
        }

        if has_supported_owner {
            schema.status = b"TRUNCATED";
            schema.code = b"TW";
        } else if has_weak_owner {
            schema.status = b"PASS";
            schema.code = b"FW";
        }

        accumulator.push(schema.to_line());
    }
}

/// Support of a record: the `CN<n>` tag of its name (`name__FC0#TC0#...#CN12`), 1 when absent.
fn cn_weight(record: &GenePred) -> usize {
    let name = String::from_utf8_lossy(record.name().unwrap_or_default());
    let tags = name.split_once("__").map_or("", |(_, tags)| tags);
    tags.split('#')
        .find_map(|tag| tag.strip_prefix("CN")?.parse().ok())
        .unwrap_or(1)
}

#[inline(always)]
fn overlaps(query_start: u64, query_end: u64, owner_start: u64, owner_end: u64) -> bool {
    (query_start >= owner_start && query_start < owner_end)
        || (query_end > owner_start && query_end <= owner_end)
        || (query_start < owner_start && query_end > owner_end)
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
    use clap::Parser;
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
            false,
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

    #[test]
    fn recovery_emits_once_per_query_after_evaluating_all_owners() {
        let queries = vec![
            record(b"supported", &[(130, 140)]),
            record(b"weak", &[(210, 220)]),
            record(b"already_pass", &[(310, 320)]),
            record(b"unchanged", &[(410, 420)]),
        ];
        let mut schemas: Vec<Schema<'_>> = Vec::new();
        for query in queries.iter() {
            let name = query.name().unwrap();
            let mut schema = Schema {
                id: name,
                status: b"TRUNCATED",
                code: b"T",
                ..Default::default()
            };

            if name == b"already_pass" {
                schema.status = b"PASS";
                schema.code = b"A";
            }

            schemas.push(schema);
        }

        let mut reference_exon_support = BTreeMap::new();
        reference_exon_support.insert((100, 150), 0.25);
        reference_exon_support.insert((120, 170), 0.75);
        reference_exon_support.insert((200, 250), 0.25);

        let mut lines = Vec::new();
        recover_reads(
            &queries,
            reference_exon_support,
            &mut schemas,
            0.5,
            0.75,
            &mut lines,
        );

        assert_eq!(lines.len(), queries.len());

        let lines = lines
            .into_iter()
            .map(|line| String::from_utf8(line).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(count_rows_for(&lines, "supported"), 1);
        assert_eq!(count_rows_for(&lines, "weak"), 1);
        assert_eq!(count_rows_for(&lines, "already_pass"), 1);
        assert_eq!(count_rows_for(&lines, "unchanged"), 1);
        assert!(lines
            .iter()
            .any(|line| line.starts_with("supported\tTRUNCATED\tTW\t")));
        assert!(lines
            .iter()
            .any(|line| line.starts_with("weak\tPASS\tFW\t")));
        assert!(lines
            .iter()
            .any(|line| line.starts_with("already_pass\tPASS\tAW\t")));
        assert!(lines
            .iter()
            .any(|line| line.starts_with("unchanged\tTRUNCATED\tT\t")));
    }

    #[test]
    fn cn_tag_weighs_a_record_like_that_many_identical_records() {
        let counter = ParallelCounter::default();
        let (first, e1, e2) = ((100, 150), (200, 250), (300, 350));
        // references share a first exon, then carry the middle exon e1, e2 or none
        let refs = [
            ("a", vec![first, e1], 3),
            ("b", vec![first, e1], 2),
            ("c", vec![first, e2], 2),
            ("d", vec![first], 3),
        ];
        // one passing query, one inside e1 and one inside e2
        let queries = [
            ("p", vec![(110, 150)], 5),
            ("t1", vec![(210, 240)], 3),
            ("t2", vec![(310, 340)], 2),
        ];
        // a CN-n record is kept as is (`weighted`) or becomes n identical untagged records
        let build = |records: &[(&str, Vec<(u64, u64)>, usize)], weighted: bool| {
            let mut out = Vec::new();
            for (name, exons, cn) in records {
                match weighted {
                    true => out.push(record(format!("{name}__FC0#CN{cn}").as_bytes(), exons)),
                    false => out
                        .extend((0..*cn).map(|k| record(format!("{name}_{k}").as_bytes(), exons))),
                }
            }
            out
        };
        let run = |weighted| -> Vec<String> {
            let (refs, queries) = (build(&refs, weighted), build(&queries, weighted));
            process_component(refs, queries, true, 0.5, 0.5, false, &counter)
                .into_iter()
                .map(|line| String::from_utf8(line).unwrap())
                .collect()
        };
        let (weighted, identical) = (run(true), run(false));

        // 5 of 10 queries are truncated, so the component is dirty right on the 0.5 threshold;
        // support is 5/10 for e1 (recovered, on the threshold) and 2/10 for e2 (not recovered)
        assert!(weighted[0].contains("\tPASS\tAW\t"), "{}", weighted[0]);
        assert!(weighted[1].contains("\tTRUNCATED\tTW\t"), "{}", weighted[1]);
        assert!(weighted[2].contains("\tPASS\tFW\t"), "{}", weighted[2]);

        // every tagged record reads like each of its identical copies
        let tail = |line: &str| line.split_once('\t').unwrap().1.to_string();
        let mut copies = identical.iter();
        for (line, (_, _, cn)) in weighted.iter().zip(&queries) {
            for _ in 0..*cn {
                assert_eq!(tail(line), tail(copies.next().unwrap()));
            }
        }
    }

    fn run_same_file(records: Vec<GenePred>) -> Vec<String> {
        run_same_file_with(records, false)
    }

    fn run_same_file_with(records: Vec<GenePred>, start_junction: bool) -> Vec<String> {
        // INFO: ARK passes the same BED as --ref and --query
        let counter = ParallelCounter::default();
        let mut lines: Vec<String> = process_component(
            records.clone(),
            records,
            false,
            0.5,
            0.5,
            start_junction,
            &counter,
        )
        .into_iter()
        .map(|line| String::from_utf8(line).unwrap())
        .collect();
        lines.sort();
        lines
    }

    #[test]
    fn start_junction_catches_diverging_and_intron_starting_fragments() {
        // container exons 0..3; its intron leaving exon 1 is [400, 500)
        let container = record(
            b"container",
            &[(100, 200), (300, 400), (500, 600), (700, 800)],
        );
        // starts in exon 1, leaves it through [400, 500), then diverges (new last exon)
        let diverging = record(b"diverging", &[(330, 400), (500, 560), (900, 950)]);
        // starts in intron 0 [200, 300) and runs into exon 1, then continues the chain
        let intronic = record(b"intronic", &[(250, 400), (500, 600), (700, 800)]);
        // starts in exon 1 but splices at 380: an alternative first exon, never truncated
        let alternative = record(b"alternative", &[(320, 380), (500, 600)]);
        // starts in intron 0 but its first intron [390, 500) is not the container's
        let intronic_other = record(b"intronic_other", &[(260, 390), (500, 600)]);
        let records = vec![container, diverging, intronic, alternative, intronic_other];

        // default: the whole chain from exon k must match; a start in intron 0 that ends at
        // exon 1's donor [400] is a truncation in both modes, one ending at 390 is not
        let lines = run_same_file(records.clone());
        for id in ["container", "diverging", "alternative", "intronic_other"] {
            assert!(
                line_for(&lines, id).contains("\tPASS\tA\t"),
                "{id}: {lines:?}"
            );
        }
        let line = line_for(&lines, "intronic");
        assert!(line.contains("\tTRUNCATED\tT\t"), "{line}");
        assert!(line.contains("exon: chr1:300-400"), "{line}");

        let lines = run_same_file_with(records, true);
        for id in ["diverging", "intronic"] {
            let line = line_for(&lines, id);
            assert!(line.contains("\tTRUNCATED\tT\t"), "{id}: {line}");
            assert!(line.contains("exon: chr1:300-400"), "{id}: {line}");
        }
        for id in ["container", "alternative", "intronic_other"] {
            assert!(
                line_for(&lines, id).contains("\tPASS\tA\t"),
                "{id}: {lines:?}"
            );
        }
    }

    #[test]
    fn start_junction_intron_start_on_the_minus_strand() {
        // '-' strand: transcript order is [500, 600) -> [300, 400) -> [100, 200); the intron
        // upstream of exon [300, 400) is [400, 500) and the one leaving it is [200, 300)
        let container = stranded(
            b"container",
            &[(100, 200), (300, 400), (500, 600)],
            Strand::Reverse,
        );
        // 5' end at 450 (inside [400, 500)), first intron [200, 300): the container's
        let intronic = stranded(b"intronic", &[(100, 200), (300, 450)], Strand::Reverse);

        for start_junction in [false, true] {
            let lines =
                run_same_file_with(vec![container.clone(), intronic.clone()], start_junction);
            let line = line_for(&lines, "intronic");
            assert!(line.contains("\tTRUNCATED\tT\t"), "{line}");
            assert!(line.contains("exon: chr1:300-400"), "{line}");
            assert!(
                line_for(&lines, "container").contains("\tPASS\tA\t"),
                "{lines:?}"
            );
        }
    }

    #[test]
    fn first_exon_sharing_a_middle_donor_and_starting_in_its_upstream_intron_is_truncated() {
        // container exons 0..3 on '+': [100, 200) [300, 400) [500, 600) [700, 800); exon 1
        // is a middle exon with donor 400, its upstream intron is [200, 300)
        let container = record(
            b"container",
            &[(100, 200), (300, 400), (500, 600), (700, 800)],
        );
        // starts in intron 0, ends at donor 400, then the container's chain
        let same_chain = record(b"same_chain", &[(250, 400), (500, 600), (700, 800)]);
        // donor 400 shared, but an alternative acceptor [550] and a new last exon: only the
        // donor has to match
        let other_acceptor = record(b"other_acceptor", &[(260, 400), (550, 650)]);
        // starts on the first base of intron 0 [the end of exon 0]: still no retention
        let intron_edge = record(b"intron_edge", &[(200, 400), (500, 600)]);
        // starts inside exon 0 and spans intron 0: a retention, not this truncation
        let retaining = record(b"retaining", &[(150, 400), (500, 600)]);
        // starts in intron 0 but its donor [390] is not exon 1's
        let other_donor = record(b"other_donor", &[(270, 390), (500, 600)]);
        // mono-exonic inside intron 0 + exon 1: no donor to share
        let mono = record(b"mono", &[(280, 400)]);
        let records = vec![
            container,
            same_chain,
            other_acceptor,
            intron_edge,
            retaining,
            other_donor,
            mono,
        ];

        for start_junction in [false, true] {
            let lines = run_same_file_with(records.clone(), start_junction);
            for id in ["same_chain", "other_acceptor", "intron_edge"] {
                let line = line_for(&lines, id);
                assert!(line.contains("\tTRUNCATED\tT\t"), "{id}: {line}");
                assert!(line.contains("exon: chr1:300-400"), "{id}: {line}");
            }
            for id in ["container", "retaining", "other_donor", "mono"] {
                assert!(
                    line_for(&lines, id).contains("\tPASS\tA\t"),
                    "{id} [-s={start_junction}]: {lines:?}"
                );
            }
        }
    }

    #[test]
    fn shared_donor_on_the_minus_strand_and_last_exons_have_no_donor() {
        // '-' strand, genomic exons [100, 200) [300, 400) [500, 600) [700, 800): transcript
        // order is [700, 800) -> [500, 600) -> [300, 400) -> [100, 200). The middle exon
        // [500, 600) has its donor at 500 and its upstream intron is [600, 700)
        let container = stranded(
            b"container",
            &[(100, 200), (300, 400), (500, 600), (700, 800)],
            Strand::Reverse,
        );
        // 5' end at 650 [inside intron 600-700], first exon [500, 650): donor 500 shared
        let truncated = stranded(b"truncated", &[(300, 400), (500, 650)], Strand::Reverse);
        // 5' end at 750 [inside exon 700-800]: retains intron [600, 700)
        let retaining = stranded(b"retaining", &[(300, 400), (500, 750)], Strand::Reverse);
        // ends at the 3' end of the container's last exon [100, 200): on '-' its donor would
        // sit at 100, but the container's last exon has no donor
        let last = stranded(b"last", &[(20, 50), (100, 250)], Strand::Reverse);

        let lines = run_same_file(vec![container, truncated, retaining, last]);
        let line = line_for(&lines, "truncated");
        assert!(line.contains("\tTRUNCATED\tT\t"), "{line}");
        assert!(line.contains("exon: chr1:500-600"), "{line}");
        for id in ["container", "retaining", "last"] {
            assert!(
                line_for(&lines, id).contains("\tPASS\tA\t"),
                "{id}: {lines:?}"
            );
        }
    }

    #[test]
    fn records_sharing_a_name_keep_their_own_rows() {
        // two models with one name: the name-keyed map kept whichever came last, so the
        // row depended on the record order; each record now keeps its own row
        let container = record(b"container", &[(100, 200), (300, 400), (500, 600)]);
        let truncated = record(b"dup", &[(320, 400), (500, 600)]);
        let pass = record(b"dup", &[(20, 60), (500, 600)]);

        let forward = run_same_file(vec![container.clone(), truncated.clone(), pass.clone()]);
        let backward = run_same_file(vec![pass, truncated, container]);
        assert_eq!(forward, backward);
        let dups: Vec<&String> = forward.iter().filter(|l| l.starts_with("dup\t")).collect();
        assert_eq!(dups.len(), 2, "{forward:?}");
        assert!(
            dups.iter().any(|l| l.contains("\tTRUNCATED\tT\t")),
            "{dups:?}"
        );
        assert!(dups.iter().any(|l| l.contains("\tPASS\tA\t")), "{dups:?}");
    }

    #[test]
    fn output_is_sorted_and_identical_across_runs() {
        let dir = tempfile::tempdir().unwrap();
        let bed = dir.path().join("models.bed");
        let mut file = File::create(&bed).unwrap();
        // three loci on both strands; one exact duplicate row
        let rows = [
            "chr1\t100\t600\tc\t0\t+\t100\t600\t0,0,0\t3\t100,100,100,\t0,200,400,",
            "chr1\t250\t600\tt\t0\t+\t250\t600\t0,0,0\t2\t150,100,\t0,250,",
            "chr1\t250\t600\tt\t0\t+\t250\t600\t0,0,0\t2\t150,100,\t0,250,",
            "chr2\t100\t800\tm\t0\t-\t100\t800\t0,0,0\t4\t100,100,100,100,\t0,200,400,600,",
            "chr2\t300\t650\tn\t0\t-\t300\t650\t0,0,0\t2\t100,150,\t0,200,",
            "chr3\t100\t400\ta\t0\t+\t100\t400\t0,0,0\t2\t100,100,\t0,200,",
        ];
        writeln!(file, "{}", rows.join("\n")).unwrap();

        let run = |prefix: &str, threads: &str| -> String {
            let args = Args::parse_from([
                "iso-utr",
                "--ref",
                bed.to_str().unwrap(),
                "--query",
                bed.to_str().unwrap(),
                "--threads",
                threads,
                "--prefix",
                prefix,
                "--outdir",
                dir.path().to_str().unwrap(),
            ]);
            detect_truncations(args).unwrap();
            std::fs::read_to_string(dir.path().join(format!("{prefix}.truncation.tsv"))).unwrap()
        };

        let first = run("a", "1");
        let ids: Vec<&str> = first
            .lines()
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        assert_eq!(ids, ["a", "c", "m", "n", "t"], "{first}");
        assert!(first.contains("t\tTRUNCATED\tT\t"), "{first}");
        assert!(first.contains("n\tTRUNCATED\tT\t"), "{first}");
        for (k, threads) in ["2", "4", "8"].iter().enumerate() {
            assert_eq!(run(&format!("r{k}"), threads), first);
        }
    }

    fn line_for<'a>(lines: &'a [String], id: &str) -> &'a str {
        let prefix = format!("{id}\t");
        lines.iter().find(|line| line.starts_with(&prefix)).unwrap()
    }

    #[test]
    fn alternative_first_exon_inside_another_exon_is_not_a_truncation() {
        // Q4 3d: SNTG1-like, '-' strand. The container's middle exon [300, 400) holds the
        // query's TSS (399), but the query splices at 350, not at the container's 300
        let container = stranded(
            b"container",
            &[(100, 200), (300, 400), (500, 600)],
            Strand::Reverse,
        );
        let alternative = stranded(b"alternative", &[(100, 200), (350, 399)], Strand::Reverse);
        // a true 5' truncation continues the container's chain
        let truncated = stranded(b"truncated", &[(100, 200), (300, 380)], Strand::Reverse);

        let lines = run_same_file(vec![container, alternative, truncated]);
        assert!(
            line_for(&lines, "alternative").contains("\tPASS\tA\t"),
            "{lines:?}"
        );
        assert!(
            line_for(&lines, "container").contains("\tPASS\tA\t"),
            "{lines:?}"
        );
        let line = line_for(&lines, "truncated");
        assert!(line.contains("\tTRUNCATED\tT\t"), "{line}");
        assert!(line.contains("exon: chr1:300-400"), "{line}");
    }

    #[test]
    fn shorter_and_mono_exonic_fragments_of_a_container_are_truncations() {
        let container = record(
            b"container",
            &[(100, 200), (300, 400), (500, 600), (700, 800)],
        );
        // starts in exon 1, continues the chain but ends before the container's last exon
        let shorter = record(b"shorter", &[(320, 400), (500, 560)]);
        // mono-exonic, inside exon 2 [empty chain]
        let mono = record(b"mono", &[(520, 580)]);
        // starts in exon 1 but skips exon 2: not the container's chain. Its own start
        // (330) differs from `shorter`'s: with the same first coding exon, this PASS would
        // rescue `shorter` through the CDS guard
        let skipping = record(b"skipping", &[(330, 400), (700, 800)]);

        let lines = run_same_file(vec![container, shorter, mono, skipping]);
        assert!(
            line_for(&lines, "shorter").contains("\tTRUNCATED\tT\t"),
            "{lines:?}"
        );
        assert!(
            line_for(&lines, "mono").contains("\tTRUNCATED\tT\t"),
            "{lines:?}"
        );
        assert!(
            line_for(&lines, "skipping").contains("\tPASS\tA\t"),
            "{lines:?}"
        );
    }

    #[test]
    fn cds_guard_does_not_depend_on_query_order() {
        // `truncated` starts inside the container's exon 2 and continues its chain, but
        // shares its first coding exon [320, 400) with `pass`, a PASS model
        let container = with_cds(
            record(b"container", &[(100, 200), (300, 400), (500, 600)]),
            100,
        );
        let truncated = with_cds(record(b"truncated", &[(320, 400), (500, 600)]), 320);
        let pass = with_cds(record(b"pass", &[(50, 80), (320, 400), (500, 600)]), 320);
        let records = [container, truncated, pass];

        let orders = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        let results: Vec<Vec<String>> = orders
            .iter()
            .map(|order| run_same_file(order.iter().map(|&i| records[i].clone()).collect()))
            .collect();

        // the old guard passed `truncated` only when `pass` came first
        for lines in &results {
            assert_eq!(lines, &results[0]);
        }
        assert!(
            line_for(&results[0], "truncated").contains("\tPASS\tA\t"),
            "{:?}",
            results[0]
        );
    }

    #[test]
    fn novel_start_is_reachable_when_ref_is_query() {
        let container = record(b"container", &[(100, 200), (300, 400)]);
        let shared = record(b"shared", &[(150, 200), (300, 400)]);
        let novel = record(b"novel", &[(20, 60), (300, 400)]);

        let lines = run_same_file(vec![container, shared, novel]);
        assert!(
            line_for(&lines, "novel").contains("support: NOVEL_START"),
            "{lines:?}"
        );
        assert!(
            line_for(&lines, "shared").contains("support: NOT_NOVEL_START"),
            "{lines:?}"
        );
    }

    fn stranded(name: &[u8], exons: &[(u64, u64)], strand: Strand) -> GenePred {
        let mut record = record(name, exons);
        record.set_strand(Some(strand));
        record
    }

    fn with_cds(mut record: GenePred, thick_start: u64) -> GenePred {
        record.set_thick_start(Some(thick_start));
        record
    }

    fn count_rows_for(lines: &[String], id: &str) -> usize {
        let prefix = format!("{id}\t");
        lines
            .iter()
            .filter(|line| line.starts_with(&prefix))
            .count()
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
        // INFO: a CDS spanning the record, which `process_component` needs for every query
        record.set_thick_start(Some(exons.first().unwrap().0));
        record.set_thick_end(Some(exons.last().unwrap().1));
        record.set_block_count(Some(exons.len() as u32));
        record.set_block_starts(Some(exons.iter().map(|(start, _)| *start).collect()));
        record.set_block_ends(Some(exons.iter().map(|(_, end)| *end).collect()));
        record
    }
}
