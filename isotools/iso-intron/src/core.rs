// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for detecting intron retentions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting intron retentions
//! and processing the components of reads and introns in parallel.
//!
//! In short, each read is checked for the presence of intron retentions
//! or RT introns. If a read has a true intron retention, it is discarded.
//! If a read has an RT intron, it is also discarded. The veracity of an
//! 'intron' is determined by 'iso-classify', using machine-learning models,
//! ab initio gene prediction, and other heuristics. The process is heavily
//! parallelized to offer fast performance on large datasets.

use dashmap::DashSet;
use genepred::GenePred;
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{OverlapType, Role};
use rayon::prelude::*;
use rust_lapper::Lapper;

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::cli::*;
use crate::utils::*;

/// Detects intron retentions in a query set of reads
///
/// # Arguments
///
/// * `args` - The command line arguments
///
/// # Returns
///
/// * `Result<()>` - The result of the operation
///
/// # Example
///
/// ```rust, no_run
/// let args = Args::new();
/// detect_intron_retentions(args).unwrap();
/// ```
pub fn detect_intron_retentions(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    info!("INFO: Detecting intron retentions...");

    let tracks = packbed::pack(vec![args.query], vec![Role::Query], OverlapType::Exon)
        .unwrap_or_else(|e| panic!("ERROR: Could not pack query -> {e}!"));
    let (index, reference_introns) = load_introns(&args.introns)
        .unwrap_or_else(|| panic!("ERROR: Could not load introns from {:?}!", args.introns));
    // INFO: --retention-source query -> a record can only retain an intron that another query
    // INFO: record splices; the table also holds introns of reads that never reach the query set
    let query_source = args.retention_source == RetentionSource::Query;
    let index = if query_source {
        query_retention_index(&tracks, &reference_introns)
    } else {
        index
    };
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let accumulator = DashSet::new();
    let counter = ParallelCounter::default();
    let allow_missing = args.allow_missing;
    let recover = args.recover;
    let ignore_utr = args.ignore_utr;

    tracks.into_par_iter().for_each(|bucket| {
        let chr = bucket.0;
        let components = bucket.1;

        counter.inc_components(components.len() as u32);

        let binding = HashSet::new();
        let banned = blacklist.get(chr.as_bytes()).unwrap_or(&binding);

        let lapper = Lapper::new(vec![Iv {
            start: 0,
            stop: 0,
            val: (),
        }]);

        // INFO: necessary to keep it strand-aware -> avoids retentions on the wrong strand
        let local_index = index.get(chr.as_bytes()).unwrap_or_else(|| {
            // INFO: check if intron collection is empty -> input was empty
            if reference_introns.is_empty() {
                log::warn!(
                    "WARN: No introns found in input -> {}. All reads will be free of IRs!",
                    args.introns.display()
                );
                &lapper
            } else {
                // INFO: weird edge case where all chr components are made up by single exons only
                // INFO: test if anything in components has introns > 0, if not, return empty index
                // INFO: with --retention-source query a strand without indexed introns has no query
                // INFO: intron in the table, so nothing on it can be retained
                if components
                    .iter()
                    .all(|comp| comp.iter().all(|read| read.introns().is_empty()))
                    || args.allow_missing
                    || query_source
                {
                    log::warn!(
                        "WARN: No introns found for {} in input -> {}. This is an edge case a whole chromsome of single exons. All reads will be free of IRs!",
                        chr,
                        args.introns.display()
                    );
                    return &lapper;
                }

                log::error!("ERROR: Could not find introns for chromosome -> {chr:?}!");
                std::process::exit(1);
            }
        });

        process_components(
            components,
            &reference_introns,
            local_index,
            banned,
            &accumulator,
            &counter,
            recover,
            ignore_utr,
            allow_missing,
        );
    });

    let file = args.outdir.join(format!("{}.retentions.tsv", args.prefix));
    let mut writer = BufWriter::new(File::create(file).unwrap());
    for schema in accumulator.into_iter() {
        writer
            .write_all(&schema)
            .unwrap_or_else(|e| panic!("ERROR: Could not write schema -> {e}!"));
    }

    info!("Reads with retained introns: {}", counter.num_retentions());
    info!("Reads with ignored introns: {}", counter.num_ignored());
    Ok(())
}

/// Builds the retention index from the table introns that at least one query record splices
///
/// # Arguments
///
/// * `tracks` - The packed query records, keyed by `chr:strand`
/// * `reference_introns` - The iso-classify table, keyed by `chr:start-end(strand)`
///
/// # Returns
///
/// * `IntronIndex` - One Lapper per `chr:strand`, holding only the query-spliced table introns
fn query_retention_index(
    tracks: &packbed::Map<String, Vec<Vec<GenePred>>>,
    reference_introns: &HashMap<Vec<u8>, Intron>,
) -> IntronIndex {
    let mut spliced: HashSet<Vec<u8>> = HashSet::new();
    for entry in tracks.iter() {
        for read in entry.value().iter().flatten() {
            let Some(strand) = read.strand() else {
                continue;
            };
            let strand = strand.to_string();
            for (start, end) in read.introns() {
                spliced.insert(intron_key(&read.chrom, start, end, &strand));
            }
        }
    }

    let mut collector: HashMap<Vec<u8>, Vec<Iv>> = HashMap::new();
    for (key, intron) in reference_introns.iter() {
        if !spliced.contains(key) {
            continue;
        }

        let mut chr_stranded = intron.chrom.clone();
        chr_stranded.push(b':');
        chr_stranded.extend_from_slice(intron.strand.to_string().as_bytes());
        collector.entry(chr_stranded).or_default().push(Iv {
            start: intron.start,
            stop: intron.end,
            val: (),
        });
    }

    let kept: usize = collector.values().map(Vec::len).sum();
    info!(
        "INFO: --retention-source query keeps {} of {} table introns as retention candidates",
        kept,
        reference_introns.len()
    );

    collector
        .into_iter()
        .map(|(chr_stranded, ivs)| (chr_stranded, Lapper::new(ivs)))
        .collect()
}

/// Lookup key of the intron table: `chr:start-end(strand)`
fn intron_key(chrom: &[u8], start: u64, end: u64, strand: &str) -> Vec<u8> {
    let mut key = Vec::with_capacity(chrom.len() + 32);
    key.extend_from_slice(chrom);
    key.push(b':');
    key.extend_from_slice(start.to_string().as_bytes());
    key.push(b'-');
    key.extend_from_slice(end.to_string().as_bytes());
    key.push(b'(');
    key.extend_from_slice(strand.as_bytes());
    key.push(b')');
    key
}

/// Processes the components of reads and introns in parallel
///
/// # Arguments
///
/// * `components` - The components to process per chromosome
/// * `banned` - the set of banned introns for the chromosome
/// * `accumulator` - the accumulator to use
/// * `counter` - The counter to use
/// * `recover` - Indicates if the component should be recovered
///
/// # Returns
///
/// * None
///
/// # Example
///
/// ```rust, no_run
/// let components = vec![];
/// let banned = HashSet::new();
/// let accumulator = ParallelAccumulator::default();
/// let counter = ParallelCounter::default();
///
/// process_components(components, &banned, &accumulator, &counter, false);
///
/// assert_eq!(accumulator.num_retentions(), 0);
/// assert_eq!(counter.num_components(), 0);
/// ```
#[inline(always)]
fn process_components(
    components: Vec<Vec<GenePred>>,
    reference_introns: &HashMap<Vec<u8>, Intron>,
    index: &Lapper<u64, ()>,
    banned: &HashSet<(u64, u64)>,
    accumulator: &DashSet<Vec<u8>>,
    counter: &ParallelCounter,
    recover: bool,
    ignore_utr: bool,
    allow_missing: bool,
) {
    components.into_iter().for_each(|component| {
        let local_collector = process_component(
            component,
            reference_introns,
            index,
            banned,
            counter,
            recover,
            ignore_utr,
            allow_missing,
        );

        local_collector.into_iter().for_each(|item| {
            accumulator.insert(item);
        });
    });
}

/// Processes a component of reads and introns
///
/// # Arguments
///
/// * `comp` - The component to process
/// * `ban` - The set of banned introns
/// * `counter` - The counter to use
/// * `recover` - Indicates if the component should be recovered
///
/// # Returns
///
/// * `Vec<String>` - The vector of reads to keep
/// * `Vec<String>` - The vector of reads to discard
/// * `Option<Vec<String>>` - The vector of reads to review
/// * `HashMap<String, Box<dyn ModuleMap>>` - The descriptor to fill up
///
/// # Example
///
/// ```rust, no_run
/// let mut comp = (vec![], vec![]);
/// let ban = HashSet::new();
/// let counter = ParallelCounter::default();
/// let recover = false;
///
/// let (keep, discard, review, descriptor) = process_component(&mut comp, &ban, &counter, recover);
///
/// assert_eq!(keep.len(), 0);
/// assert_eq!(discard.len(), 0);
/// assert_eq!(review, None);
/// assert_eq!(descriptor.len(), 0);
/// ```
#[inline(always)]
pub fn process_component<'a>(
    component: Vec<GenePred>,
    reference_introns: &'a HashMap<Vec<u8>, Intron>,
    index: &Lapper<u64, ()>,
    ban: &HashSet<(u64, u64)>,
    counter: &ParallelCounter,
    recover: bool,
    ignore_utr: bool,
    allow_missing: bool,
) -> Vec<Vec<u8>> {
    let mut descriptor: HashMap<&[u8], Schema<'a>> = HashMap::new();
    let mut accumulator: Vec<Vec<u8>> = Vec::new();
    let (mut count, totals) = (0_f32, component.len() as f32);

    // INFO: for every read -> see which introns are retained and which the read has!
    for read in component.iter() {
        let mut schema = Schema::default();

        schema.id = read
            .name()
            .unwrap_or_else(|| panic!("ERROR: Read has no name!"))
            .to_vec();
        schema.size = component.len() as u32;

        detect_rt_intron(read, ban, reference_introns, &mut schema, allow_missing);
        detect_retention(
            read,
            reference_introns,
            index,
            &mut schema,
            ignore_utr,
            allow_missing,
        );

        if schema.events > 0 {
            count += 1.0;
            counter.inc_retentions();
        }

        if !schema.ig_html.is_empty() {
            counter.inc_ignored();
        }

        descriptor.insert(read.name().unwrap(), schema.clone());
    }

    let ratio = count / totals;

    // INFO: second pass to fill up the descriptor with comp ratio
    for read in component.iter() {
        let schema = descriptor
            .get_mut(&read.name().unwrap())
            .unwrap_or_else(|| {
                panic!(
                    "ERROR: Read not found in component, this is likely a bug! -> {}",
                    std::str::from_utf8(read.name().unwrap()).unwrap()
                );
            });

        schema.ratio = ratio;
        schema.component_events = count as u32;
        schema.component_size = component.len() as u32;

        if recover && ratio > RETENTION_RATIO_THRESHOLD {
            // INFO: set component status to 'REVIEW'
            // INFO: check current class and append REVIEW not replace
            let mut new_status = schema.status.clone();
            new_status.extend_from_slice(b"/REVIEW");
            schema.status = new_status;
        }

        accumulator.push(schema.to_line());
    }

    accumulator
}

/// Schema struct
///
/// # Fields
///
/// * `status`: Vec<u8> - The status of the component
/// * `ratio`: f32 - The ratio of the component     
/// * `events`: u32 - The number of events in the component     
/// * `code`: Vec<u8> - The code of the component
/// * `ir_html`: Vec<Intron> - The html of the intron retentions
/// * `fr_html`: Vec<Intron> - The html of the false retentions
/// * `rt_html`: Vec<Intron> - The html of the true retentions
#[derive(Debug, Clone, PartialEq, Default)]
struct Schema<'a> {
    pub id: Vec<u8>,
    pub status: Vec<u8>,
    pub ratio: f32,
    pub component_events: u32,
    pub component_size: u32,
    pub events: u32,
    pub code: Vec<u8>,
    pub size: u32,
    pub ir_html: Vec<&'a Intron>,
    pub fr_html: Vec<&'a Intron>,
    pub rt_html: Vec<&'a Intron>,
    pub ig_html: Vec<&'a Intron>,
}

impl Schema<'_> {
    /// INFO: fmt -> id\tstatus\tcode\thtml
    ///
    /// # Arguments
    ///
    /// * `&self` - This function takes a reference to a `Schema` struct.
    ///
    /// # Returns
    ///
    /// * `Vec<u8>` - A vector of bytes representing the RT introns in the schema.
    ///
    /// # Notes
    ///
    /// html format is:
    ///
    /// /h2 Intron retentions
    /// - status: {value}
    /// - code: {value} [R: retention, T: RT retention, X: has RT intron, A: no events]
    /// - events: {value}
    /// - ratio: {value}
    ///
    /// /h3 True retentions:
    ///
    /// {table}
    ///
    /// /h3 RT retentions:
    ///
    /// {table}
    ///
    /// /h3 Spliced RT introns in read:
    ///
    /// {table}
    ///
    pub fn to_line(&self) -> Vec<u8> {
        let mut body = String::new();

        body.push_str("<h2>Intron retentions</h2><br>");
        body.push_str(&format!(
            "- status: {}<br>",
            std::str::from_utf8(&self.status).unwrap_or("NULL")
        ));
        body.push_str(&format!(
            "- code: {} [R: retention, K: RT retention, X: has RT intron, G: retains artifact, A: no events, Q: has artifact, W: ignored retention]<br>",
            std::str::from_utf8(&self.code).unwrap_or("NULL")
        ));
        body.push_str(&format!("- events: {}<br>", self.events));
        body.push_str(&format!(
            "- ratio: {} ({}/{})<br>",
            self.ratio, self.component_events, self.component_size
        ));

        body.push_str("<h3>True retentions:</h3><br>");
        if self.ir_html.is_empty() {
            body.push_str("none<br>");
        } else {
            for intron in &self.ir_html {
                body.push_str(&intron_html(intron));
            }
        }

        body.push_str("<h3>RT/Artifact retentions:</h3><br>");
        if self.fr_html.is_empty() {
            body.push_str("none<br>");
        } else {
            for intron in &self.fr_html {
                body.push_str(&intron_html(intron));
            }
        }

        body.push_str("<h3>Ignored retentions:</h3><br>");
        if self.ig_html.is_empty() {
            body.push_str("none<br>");
        } else {
            for intron in &self.ig_html {
                body.push_str(&intron_html(intron));
            }
        }

        body.push_str("<h3>Spliced RT introns in read:</h3><br>");
        if self.rt_html.is_empty() {
            body.push_str("none<br>");
        } else {
            for intron in &self.rt_html {
                body.push_str(&intron_html(intron));
            }
        }

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

/// One intron as an html table row of the descriptor column
///
/// The Display impl is the tab-separated iso-classify row: its tabs would split the
/// html into extra TSV/BED columns downstream (POLISH BEDs with 22-94 columns), so
/// tabs become spaces and newlines are dropped; the descriptor stays one column.
fn intron_html(intron: &Intron) -> String {
    format!("{}<br>", intron)
        .replace('\n', "")
        .replace('\t', " ")
}

/// Detects if a read has any RT introns
///
/// # Arguments
///
/// * `read` - The read to check for intron retentions
/// * `reference_introns` - The reference introns
/// * `index` - The index of the reference introns
/// * `schema` - The schema to set the values in
fn detect_rt_intron<'a>(
    read: &GenePred,
    ban: &HashSet<(u64, u64)>,
    reference_introns: &'a HashMap<Vec<u8>, Intron>,
    schema: &mut Schema<'a>,
    allow_missing: bool,
) {
    let read_introns = read.introns();
    let mut rt = false;

    for intron in read_introns {
        if ban.contains(&(intron.0, intron.1)) {
            continue;
        }

        // INFO: build lookup key
        let mut key = Vec::new();
        key.extend_from_slice(&read.chrom);
        key.extend_from_slice(b":");
        key.extend_from_slice(intron.0.to_string().as_bytes());
        key.extend_from_slice(b"-");
        key.extend_from_slice(intron.1.to_string().as_bytes());
        key.extend_from_slice(b"(");
        key.extend_from_slice(read.strand().unwrap().to_string().as_bytes());
        key.extend_from_slice(b")");

        let Some(info) = reference_introns.get(&key) else {
            if allow_missing {
                log::warn!(
                    "WARN: Intron not found but allowed because of --allow-missing! -> {:?}",
                    key
                );
                continue;
            } else {
                panic!(
                    "ERROR: Intron not found, this is likely a bug! -> {:?}",
                    key
                );
            }
        };

        match info.support {
            SupportType::StrongRT => {
                if !rt {
                    rt = true;
                    // INFO: add 'X' to code representing that read has RT intron
                    if !schema.code.contains(&b'X') {
                        schema.code.push(b'X');
                    }

                    if schema.status != b"RETENTION" && schema.status != b"HAS_STRONG_RT" {
                        schema.status = b"HAS_STRONG_RT".to_vec();
                    }
                }

                schema.rt_html.push(info);
            }
            SupportType::WeakRT => {
                if !rt {
                    rt = true;
                    // INFO: add 'X' to code representing that read has RT intron
                    if !schema.code.contains(&b'X') {
                        schema.code.push(b'X');
                    }

                    if schema.status != b"RETENTION"
                        && schema.status != b"HAS_STRONG_RT"
                        && schema.status != b"HAS_WEAK_RT"
                    {
                        schema.status = b"HAS_WEAK_RT".to_vec();
                    }
                }

                schema.rt_html.push(info);
            }
            SupportType::Artifact => {
                if !rt {
                    rt = true;
                    // INFO: add 'Q' to code representing that read has an artifact
                    if !schema.code.contains(&b'Q') {
                        schema.code.push(b'Q');
                    }

                    if schema.status != b"RETENTION"
                        && schema.status != b"HAS_STRONG_RT"
                        && schema.status != b"HAS_ARTIFACT"
                    {
                        schema.status = b"HAS_ARTIFACT".to_vec();
                    }
                }

                schema.rt_html.push(info);
            }
            _ => {}
        }
    }
}

/// Determine if any read exons retain any intron
///
/// # Arguments
///
/// * `read` - The read to check for intron retention
/// * `ref_introns` - The reference introns
/// * `schema` - The schema to set the values in
///
/// # Returns
///
/// * None -> the schema is modified in place
///
/// # Example
///
/// ```rust, no_run
/// let mut schema = Schema::default();
/// detect_retention(&read, &ref_introns, &mut schema);
///
/// match schema.intron_retention {
///     true => println!("Intron retention detected"),
///     false => println!("No intron retention detected"),
/// }
/// ```
fn detect_retention<'a>(
    read: &GenePred,
    reference_introns: &'a HashMap<Vec<u8>, Intron>,
    index: &Lapper<u64, ()>,
    schema: &mut Schema<'a>,
    ignore_utr: bool,
    allow_missing: bool,
) {
    let read_exons = read.exons();
    let mut flawed = schema.code.contains(&b'Q') || schema.code.contains(&b'X');

    for exon in read_exons {
        for (retention_start, retention_stop) in intron_retentions(index, exon) {
            let mut key = Vec::new();
            key.extend_from_slice(&read.chrom);
            key.extend_from_slice(b":");
            key.extend_from_slice(retention_start.to_string().as_bytes());
            key.extend_from_slice(b"-");
            key.extend_from_slice(retention_stop.to_string().as_bytes());
            key.extend_from_slice(b"(");
            key.extend_from_slice(read.strand().unwrap().to_string().as_bytes());
            key.extend_from_slice(b")");

            let Some(info) = reference_introns.get(&key) else {
                if allow_missing {
                    log::warn!(
                        "WARN: Intron not found but allowed because of --allow-missing! -> {} from {}",
                        std::str::from_utf8(&key).unwrap(),
                        read
                    );
                    continue;
                } else {
                    panic!(
                        "ERROR: Intron not found, this is likely a bug! -> {} from {}",
                        std::str::from_utf8(&key).unwrap(),
                        read
                    );
                }
            };

            match info.support {
                // INFO: add 'K' to code representing that read retains true RT intron
                SupportType::StrongRT | SupportType::WeakRT => {
                    if !schema.code.contains(&b'K') {
                        schema.code.push(b'K');
                    }
                    schema.fr_html.push(info);
                }
                // INFO: add 'G' to code representing that read retains artifact intron
                SupportType::Artifact => {
                    if !schema.code.contains(&b'G') {
                        schema.code.push(b'G');
                    }
                    schema.fr_html.push(info);
                }
                SupportType::Unclear | SupportType::Splicing => {
                    if !ignore_utr {
                        if schema.events == 0 {
                            // INFO: add 'R' to code representing that reads retains a true intron
                            schema.code.push(b'R');
                            flawed = true;

                            if schema.status != b"RETENTION" {
                                schema.status = b"RETENTION".to_vec();
                            }
                        }

                        schema.events += 1;
                        schema.ir_html.push(info);
                    } else {
                        match info.intron_position {
                            Position::UTR => {
                                // INFO: we must ignore this retention but still count it
                                schema.events += 1;
                                schema.ig_html.push(info);

                                if !schema.code.contains(&b'W') {
                                    schema.code.push(b'W');
                                }
                            }
                            Position::CDS | Position::Unknown => {
                                // INFO: true retention
                                if schema.events == 0 {
                                    // INFO: add 'R' to code representing that reads retains a true intron
                                    schema.code.push(b'R');
                                    flawed = true;

                                    if schema.status != b"RETENTION" {
                                        schema.status = b"RETENTION".to_vec();
                                    }
                                }

                                schema.events += 1;
                                schema.ir_html.push(info);
                            }
                        }
                    }
                }
            }
        }
    }

    if !flawed {
        schema.code.push(b'A');
        schema.status = b"PASS".to_vec();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // iso-classify row: 25 tab-separated columns, `seen` is column 5
    fn intron_row(start: u64, end: u64, seen: usize) -> Intron {
        let line = format!(
            "chr1\t{start}\t{end}\t+\t{seen}\t{seen}\t0.9\t0.9\t8.0\t8.0\tGT\tAG\tCAGGTAAGT\tTTTTTTTTTTTTTTTTTTTCAGG\tNEW\tNO_REFERENCE_SUPPORT\tIN_FRAME\tAAAAAAAAAAAA\tCCCCCCCCCCCC\tNOT_RT_INTRON\tNOT_NAG_SS\tU2\tNO_REPEAT\t{}\tSPLICED",
            end - start
        );
        Intron::from(&line).unwrap()
    }

    fn run(model: &str, seen: usize) -> String {
        let intron = intron_row(2000, 3000, seen);
        let mut introns = HashMap::new();
        introns.insert(b"chr1:2000-3000(+)".to_vec(), intron);
        let index = Lapper::new(vec![Iv {
            start: 2000,
            stop: 3000,
            val: (),
        }]);

        // single-exon model over the intron
        let mut read = GenePred::from_coords(b"chr1".to_vec(), 1000, 5000, Default::default());
        read.set_name(Some(model.as_bytes().to_vec()));
        read.set_strand(Some(genepred::Strand::Forward));

        let counter = ParallelCounter::default();
        let lines = process_component(
            vec![read],
            &introns,
            &index,
            &HashSet::new(),
            &counter,
            false,
            false,
            false,
        );
        String::from_utf8(lines[0].clone()).unwrap()
    }

    #[test]
    fn retained_spliced_intron_is_a_retention_regardless_of_support() {
        // owner rule: a retention of a true intron is a retention; read support is not
        // required (every intron is validated separately downstream)
        for seen in [1, 40] {
            let line = run("R632@54__TC0#CN150", seen);
            assert!(line.contains("\tRETENTION\tR\t"), "{line}");
        }
    }

    #[test]
    fn descriptor_html_stays_one_column() {
        // Q4 5a: tabs of the iso-classify row split the html into 22-94 BED columns
        for seen in [1, 40] {
            let line = run("R632@54__TC0#CN150", seen);
            assert_eq!(line.trim_end().split('\t').count(), 4, "{line}");
            assert!(line.contains("chr1 2000 3000 +"), "{line}");
        }
    }

    /// Runs the binary entry point on a two-record query and returns the `wide` model's row
    fn run_source(source: &str) -> String {
        use clap::Parser;

        let dir = std::env::temp_dir().join(format!(
            "iso_intron_source_{}_{}",
            source,
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();

        // `wide` is one exon over both table introns; `spliced` splices only 2000-3000.
        // 3500-4000 comes from a read outside the query set (e.g. one that never got an ORF)
        let query = dir.join("query.bed");
        std::fs::write(
            &query,
            "chr1\t1000\t5000\twide\t0\t+\t1000\t5000\t0,0,0\t1\t4000,\t0,\n\
             chr1\t1000\t5000\tspliced\t0\t+\t1000\t5000\t0,0,0\t2\t1000,2000,\t0,2000,\n",
        )
        .unwrap();
        let introns = dir.join("introns.tsv");
        let rows: Vec<String> = [(2000, 3000), (3500, 4000)]
            .iter()
            .map(|(start, end)| {
                format!(
                    "chr1\t{start}\t{end}\t+\t5\t5\t0.9\t0.9\t8.0\t8.0\tGT\tAG\tCAGGTAAGT\tTTTTTTTTTTTTTTTTTTTCAGG\tNEW\tNO_REFERENCE_SUPPORT\tIN_FRAME\tAAAAAAAAAAAA\tCCCCCCCCCCCC\tNOT_RT_INTRON\tNOT_NAG_SS\tU2\tNO_REPEAT\t{}\tSPLICED\n",
                    end - start
                )
            })
            .collect();
        std::fs::write(&introns, rows.concat()).unwrap();

        let args = Args::parse_from([
            "iso-intron",
            "--introns",
            introns.to_str().unwrap(),
            "--query",
            query.to_str().unwrap(),
            "--outdir",
            dir.to_str().unwrap(),
            "--prefix",
            "t",
            "--threads",
            "1",
            "--retention-source",
            source,
        ]);
        detect_intron_retentions(args).unwrap();

        let out = std::fs::read_to_string(dir.join("t.retentions.tsv")).unwrap();
        std::fs::remove_dir_all(&dir).ok();
        out.lines()
            .find(|line| line.starts_with("wide\t"))
            .unwrap()
            .to_string()
    }

    #[test]
    fn retention_source_query_ignores_introns_no_query_record_splices() {
        let all = run_source("all");
        assert!(all.contains("- events: 2<br>"), "{all}");
        assert!(all.contains("chr1 3500 4000 +"), "{all}");

        let query = run_source("query");
        assert!(query.contains("\tRETENTION\tR\t"), "{query}");
        assert!(query.contains("- events: 1<br>"), "{query}");
        assert!(query.contains("chr1 2000 3000 +"), "{query}");
        assert!(!query.contains("chr1 3500 4000 +"), "{query}");
    }
}
