// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for segmenting reads based on their polyA features
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for grouping reads
//! and processing components based on polyA features in parallel.
//!
//! In short, this modules provides iso-segment. The segment module
//! filters reads based on alignment quality and predicts the polyA
//! tail using a two-state HMM model.

use std::{
    collections::HashMap,
    fmt::Write as FmtWrite,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{
        mpsc::{sync_channel, Receiver, SyncSender},
        Arc,
    },
    thread::JoinHandle,
};

use crate::{cli::Args, constants::*};

use bio::stats::hmm::{discrete_emission::Model, viterbi, State};
use ndarray::{array, Array2};
use noodles_bam::{io::writer::Builder, Record};
use noodles_core::{region::Interval, Position};
use noodles_sam::{
    alignment::{
        io::Write as BamWrite,
        record::{cigar::op::Kind, data::field::Tag, Cigar},
        RecordBuf,
    },
    Header,
};
use rayon::prelude::*;

const TAG_CORRECTED: Tag = Tag::new(b'x', b'c');

/// Segment reads based on their polyA features
///
/// # Arguments
///
/// * `args` - The command-line arguments.
///
/// # Returns
///
/// * `Result<(), String>` - An error message if the segmentation fails.
///
/// # Example
///
/// ```rust, no_run
/// use iso_segment::segment;
///
/// let args = Args::default();
/// segment(args).unwrap();
/// ```
pub fn segment(args: Args) -> Result<(), String> {
    let (chroms, refs, header) = read_bam(&args.bam);
    let (outputs, handles) = spawn_output_workers(&args, &header, Arc::clone(&refs));

    // INFO: process each reference region in parallel
    chroms.par_iter().for_each(|chr| {
        // INFO: thread-local copy of the reference sequences
        let references = Arc::clone(&refs);
        let outputs = outputs.clone();

        let mut reader = noodles_bam::io::indexed_reader::Builder::default()
            .build_from_path(&args.bam)
            .unwrap_or_else(|e| panic!("ERROR: could not open file: {e}"));

        // INFO: must re-read header per thread since IndexedReader owns its header
        let header = reader
            .read_header()
            .unwrap_or_else(|e| panic!("ERROR: could not open file: {e}"));
        let hmm = HMM::init(args.p2p, args.emit_a);

        let mut track = 0;

        let end = Position::new(references[chr.as_bytes()].length().get()).unwrap();
        let bounds = Interval::from(Position::MIN..=end);
        let region = noodles_core::Region::new(chr.as_str(), bounds);

        reader
            .query(&header, &region)
            .unwrap_or_else(|e| panic!("ERROR: could not query genomic region: {e}"))
            .for_each(|result| {
                if let Ok(record) = result {
                    track += 1;
                    process_record(
                        record,
                        &header,
                        track,
                        chr.as_str(),
                        &args,
                        &hmm,
                        &outputs,
                        args.singleton,
                        args.fragmented,
                    );
                }
            });
    });

    drop(outputs);

    for (name, handle) in handles {
        handle
            .join()
            .map_err(|_| format!("ERROR: could not join {name}"))?;
    }

    Ok(())
}

/// Container for output channels to accept/reject writers
///
/// Provides synchronized senders for distributing processed records
/// to either high-quality (accepted) or low-quality (rejected) outputs.
///
/// # Example
///
/// ```rust, no_run
/// use std::sync::mpsc::sync_channel;
///
/// let (accept, _) = sync_channel(100);
/// let (reject, _) = sync_channel(100);
/// let senders = OutputSenders { accept, reject };
/// ```
#[derive(Clone)]
struct OutputSenders {
    accept: SyncSender<RecordBuf>,
    reject: SyncSender<RecordBuf>,
}

impl OutputSenders {
    /// Send record to appropriate output channel
    ///
    /// # Arguments
    ///
    /// * `accepted` - Whether record passed quality filters
    /// * `record` - The record to send
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use noodles_sam::alignment::record::RecordBuf;
    ///
    /// let senders = OutputSenders { accept, reject };
    /// let record = RecordBuf::new();
    /// senders.send(true, record);
    /// ```
    fn send(&self, accepted: bool, record: RecordBuf) {
        let sender = if accepted { &self.accept } else { &self.reject };
        sender
            .send(record)
            .expect("ERROR: failed to send record to output writer");
    }
}

/// Build output file path based on parameters
///
/// # Arguments
///
/// * `outdir` - Output directory path
/// * `prefix` - File name prefix
/// * `split` - Whether to split by chromosome
/// * `delimiter` - Delimiter for chromosome prefix
/// * `chromosome` - Chromosome name (if splitting)
/// * `quality_suffix` - Quality prefix (hq/lq)
/// * `extension` - File extension
///
/// # Returns
///
/// * `PathBuf` - Full output path
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
///
/// let path = build_output_path(
///     &PathBuf::from("out"),
///     &PathBuf::from("sample"),
///     true,
///     "@",
///     Some("chr1"),
///     "hq",
///     "bam",
/// );
/// ```
#[inline(always)]
fn build_output_path(
    outdir: &Path,
    prefix: &Path,
    split: bool,
    delimiter: &str,
    chromosome: Option<&str>,
    quality_suffix: &str,
    extension: &str,
) -> PathBuf {
    let stem = if split {
        chromosome
            .map(|chr| format!("{chr}{delimiter}{}", prefix.display()))
            .unwrap_or_else(|| prefix.display().to_string())
    } else {
        prefix.display().to_string()
    };

    outdir.join(format!("{stem}.{quality_suffix}.{extension}"))
}

/// Spawn output worker threads for BAM/BED writing
///
/// Creates threaded workers that write records to output files
/// based on the configured output format (BAM or BED).
///
/// # Arguments
///
/// * `args` - Command-line arguments
/// * `header` - BAM header
/// * `refs` - Reference sequences
///
/// # Returns
///
/// * Tuple of (OutputSenders, handles)
///
/// # Example
///
/// ```rust, no_run
/// use iso_segment::cli::Args;
/// use noodles_sam::Header;
///
/// let args = Args::default();
/// let header = Header::default();
/// let refs = std::sync::Arc::new(noodles_sam::header::ReferenceSequences::new());
/// let (senders, handles) = spawn_output_workers(&args, &header, refs);
/// ```
fn spawn_output_workers(
    args: &Args,
    header: &Header,
    refs: Arc<noodles_sam::header::ReferenceSequences>,
) -> (OutputSenders, Vec<(&'static str, JoinHandle<()>)>) {
    let capacity = (args.threads.max(1) * 32).clamp(128, 1024);

    if args.bed {
        let (accept, accept_handle) = spawn_bed_output_worker(
            args.clone(),
            Arc::clone(&refs),
            PASS_PREFIX,
            RGB_ACCEPT,
            capacity,
        );
        let (reject, reject_handle) =
            spawn_bed_output_worker(args.clone(), refs, FAIL_PREFIX, RGB_REJECT, capacity);

        (
            OutputSenders { accept, reject },
            vec![
                ("accept BED writer", accept_handle),
                ("reject BED writer", reject_handle),
            ],
        )
    } else {
        let header = Arc::new(header.clone());
        let (accept, accept_handle) =
            spawn_bam_output_worker(args.clone(), Arc::clone(&header), PASS_PREFIX, capacity);
        let (reject, reject_handle) =
            spawn_bam_output_worker(args.clone(), header, FAIL_PREFIX, capacity);

        (
            OutputSenders { accept, reject },
            vec![
                ("accept BAM writer", accept_handle),
                ("reject BAM writer", reject_handle),
            ],
        )
    }
}

/// Spawn a BAM output worker thread
///
/// Creates a threaded worker for writing BAM records to disk.
///
/// # Arguments
///
/// * `args` - Command-line arguments
/// * `header` - BAM header
/// * `quality_suffix` - Quality prefix (hq/lq)
/// * `capacity` - Channel capacity
///
/// # Returns
///
/// * Tuple of (sender, handle)
///
/// # Example
///
/// ```rust, no_run
/// use iso_segment::cli::Args;
/// use noodles_sam::Header;
///
/// let args = Args::default();
/// let header = Arc::new(Header::default());
/// let (sender, handle) = spawn_bam_output_worker(args, header, "hq", 128);
/// ```
fn spawn_bam_output_worker(
    args: Args,
    header: Arc<Header>,
    quality_suffix: &'static str,
    capacity: usize,
) -> (SyncSender<RecordBuf>, JoinHandle<()>) {
    let (sender, receiver) = sync_channel(capacity);
    let handle =
        std::thread::spawn(move || write_bam_records(receiver, args, header, quality_suffix));
    (sender, handle)
}

/// Spawn a BED output worker thread
///
/// Creates a threaded worker for writing BED records to disk.
///
/// # Arguments
///
/// * `args` - Command-line arguments
/// * `refs` - Reference sequences
/// * `quality_suffix` - Quality prefix (hq/lq)
/// * `rgb` - RGB color string
/// * `capacity` - Channel capacity
///
/// # Returns
///
/// * Tuple of (sender, handle)
///
/// # Example
///
/// ```rust, no_run
/// use iso_segment::cli::Args;
/// use std::sync::Arc;
///
/// let args = Args::default();
/// let refs = Arc::new(noodles_sam::header::ReferenceSequences::new());
/// let (sender, handle) = spawn_bed_output_worker(args, refs, "hq", "43,118,219", 128);
/// ```
fn spawn_bed_output_worker(
    args: Args,
    refs: Arc<noodles_sam::header::ReferenceSequences>,
    quality_suffix: &'static str,
    rgb: &'static str,
    capacity: usize,
) -> (SyncSender<RecordBuf>, JoinHandle<()>) {
    let (sender, receiver) = sync_channel(capacity);
    let handle =
        std::thread::spawn(move || write_bed_records(receiver, args, refs, quality_suffix, rgb));
    (sender, handle)
}

/// Write BAM records to output files
///
/// Continuously receives records from channel and writes to BAM files.
/// Supports both split (by chromosome) and unified output modes.
///
/// # Arguments
///
/// * `receiver` - Channel receiver for records
/// * `args` - Command-line arguments
/// * `header` - BAM header
/// * `quality_suffix` - Quality prefix (hq/lq)
///
/// # Example
///
/// ```rust, no_run
/// use std::sync::mpsc::sync_channel;
/// use iso_segment::cli::Args;
///
/// let (tx, rx) = sync_channel(100);
/// let args = Args::default();
/// // write_bam_records(rx, args, header, "hq");
/// ```
fn write_bam_records(
    receiver: Receiver<RecordBuf>,
    args: Args,
    header: Arc<Header>,
    quality_suffix: &'static str,
) {
    log::info!("INFO: Writing streamed BAM records to {quality_suffix} output");

    if args.split {
        let refs = header.reference_sequences();
        let mut writers = HashMap::new();

        while let Ok(record) = receiver.recv() {
            let chromosome = resolve_record_chromosome(&record, refs);
            let writer = writers.entry(chromosome.clone()).or_insert_with(|| {
                let output = build_output_path(
                    &args.outdir,
                    &args.prefix,
                    args.split,
                    &args.delimiter,
                    Some(&chromosome),
                    quality_suffix,
                    "bam",
                );

                let mut writer = Builder.build_from_path(&output).unwrap_or_else(|_| {
                    panic!("ERROR: could not create file: {}", output.display())
                });

                writer
                    .write_header(&header)
                    .expect("ERROR: failed to write header");

                writer
            });

            writer
                .write_alignment_record(&header, &record)
                .expect("ERROR: failed to write record");
        }

        for (_, mut writer) in writers {
            writer
                .finish(&header)
                .expect("ERROR: failed to finish writing BAM file");
        }
    } else {
        let output = build_output_path(
            &args.outdir,
            &args.prefix,
            args.split,
            &args.delimiter,
            None,
            quality_suffix,
            "bam",
        );

        let mut writer = Builder
            .build_from_path(&output)
            .unwrap_or_else(|_| panic!("ERROR: could not create file: {}", output.display()));

        writer
            .write_header(&header)
            .expect("ERROR: failed to write header");

        while let Ok(record) = receiver.recv() {
            writer
                .write_alignment_record(&header, &record)
                .expect("ERROR: failed to write record");
        }

        writer
            .finish(&header)
            .expect("ERROR: failed to finish writing BAM file");
    }
}

/// Write BED records to output files
///
/// Continuously receives records from channel and writes to BED files.
/// Supports both split (by chromosome) and unified output modes.
///
/// # Arguments
///
/// * `receiver` - Channel receiver for records
/// * `args` - Command-line arguments
/// * `refs` - Reference sequences
/// * `quality_suffix` - Quality prefix (hq/lq)
/// * `rgb` - RGB color string
///
/// # Example
///
/// ```rust, no_run
/// use std::sync::mpsc::sync_channel;
/// use iso_segment::cli::Args;
/// use std::sync::Arc;
///
/// let (tx, rx) = sync_channel(100);
/// let args = Args::default();
/// let refs = Arc::new(noodles_sam::header::ReferenceSequences::new());
/// // write_bed_records(rx, args, refs, "hq", "43,118,219");
/// ```
fn write_bed_records(
    receiver: Receiver<RecordBuf>,
    args: Args,
    refs: Arc<noodles_sam::header::ReferenceSequences>,
    quality_suffix: &'static str,
    rgb: &'static str,
) {
    log::info!("INFO: Writing streamed BED records to {quality_suffix} output");

    if args.split {
        let mut writers = HashMap::new();

        while let Ok(record) = receiver.recv() {
            let Some((chromosome, line)) = record_to_bed(&record, refs.as_ref(), rgb) else {
                continue;
            };

            let writer = writers.entry(chromosome.clone()).or_insert_with(|| {
                let output = build_output_path(
                    &args.outdir,
                    &args.prefix,
                    args.split,
                    &args.delimiter,
                    Some(&chromosome),
                    quality_suffix,
                    "bed",
                );

                BufWriter::new(
                    File::create(&output)
                        .unwrap_or_else(|e| panic!("ERROR: could not create file: {e}")),
                )
            });

            writeln!(writer, "{line}").expect("ERROR: failed to write BED line");
        }

        for (_, mut writer) in writers {
            writer.flush().expect("ERROR: failed to flush BED file");
        }
    } else {
        let output = build_output_path(
            &args.outdir,
            &args.prefix,
            args.split,
            &args.delimiter,
            None,
            quality_suffix,
            "bed",
        );

        let mut writer = BufWriter::new(
            File::create(&output).unwrap_or_else(|e| panic!("ERROR: could not create file: {e}")),
        );

        while let Ok(record) = receiver.recv() {
            let Some((_, line)) = record_to_bed(&record, refs.as_ref(), rgb) else {
                continue;
            };

            writeln!(writer, "{line}").expect("ERROR: failed to write BED line");
        }

        writer.flush().expect("ERROR: failed to flush BED file");
    }
}

/// Resolve chromosome name from a BAM record
///
/// # Arguments
///
/// * `record` - BAM record
/// * `refs` - Reference sequences
///
/// # Returns
///
/// * `String` - Chromosome name or "unmapped"
///
/// # Example
///
/// ```rust, no_run
/// use iso_segment::cli::Args;
/// // let chrom = resolve_record_chromosome(&record, &refs);
/// ```
fn resolve_record_chromosome(
    record: &RecordBuf,
    refs: &noodles_sam::header::ReferenceSequences,
) -> String {
    record
        .reference_sequence_id()
        .and_then(|id| refs.get_index(id).map(|(name, _)| name.to_string()))
        .unwrap_or_else(|| "unmapped".to_string())
}

/// Convert BAM record to BED format
///
/// # Arguments
///
/// * `record` - BAM record
/// * `refs` - Reference sequences
/// * `rgb` - RGB color string
///
/// # Returns
///
/// * Option of (chromosome, BED line)
///
/// # Example
///
/// ```rust, no_run
/// // let bed = record_to_bed(&record, &refs, "43,118,219");
/// ```
fn record_to_bed(
    record: &RecordBuf,
    refs: &noodles_sam::header::ReferenceSequences,
    rgb: &str,
) -> Option<(String, String)> {
    let chromosome = resolve_record_chromosome(record, refs);
    let line = record_to_bed_line(record, &chromosome, rgb)?;
    Some((chromosome, line))
}

/// Generate BED line from BAM record
///
/// # Arguments
///
/// * `record` - BAM record
/// * `chr` - Chromosome name
/// * `rgb` - RGB color string
///
/// # Returns
///
/// * Option of BED line string
///
/// # Example
///
/// ```rust, no_run
/// // let line = record_to_bed_line(&record, "chr1", "43,118,219");
/// ```
#[inline(always)]
fn record_to_bed_line(record: &RecordBuf, chr: &str, rgb: &str) -> Option<String> {
    let start = record.alignment_start()?.get() - 1;
    let score = record
        .mapping_quality()
        .map(|q| q.get())
        .unwrap_or(0)
        .to_string();

    let mut name = record.name().unwrap().to_string();
    let strand = if record.flags().is_reverse_complemented() {
        '-'
    } else {
        '+'
    };

    let is_corrected = record.data().get(&TAG_CORRECTED).is_some();
    if is_corrected {
        name = format!("C{name}").into();
    }

    let mut blocks = Vec::new();
    let mut ref_pos = start;
    let mut block_start = ref_pos;
    let mut block_len = 0;

    let cigar = record.cigar();
    for c in cigar.iter() {
        let c = c.unwrap();
        let len = c.len();
        match c.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                block_len += len;
                ref_pos += len;
            }
            Kind::Skip => {
                if block_len > 0 {
                    blocks.push((block_start, block_len));
                }

                ref_pos += len;
                block_start = ref_pos;
                block_len = 0;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => continue,
        }
    }

    if block_len > 0 {
        blocks.push((block_start, block_len));
    }

    if blocks.is_empty() {
        return None;
    }

    let chrom_end = blocks.last().map(|(s, l)| s + l).unwrap_or(start);
    let thick_start = start;
    let thick_end = chrom_end;

    let mut block_sizes = String::with_capacity(16 * blocks.len());
    let mut block_starts = String::with_capacity(16 * blocks.len());

    for (s, l) in &blocks {
        let _ = write!(block_sizes, "{},", l);
        let _ = write!(block_starts, "{},", s - start);
    }
    block_sizes.pop();
    block_starts.pop();

    let mut line = String::with_capacity(256);
    write!(
        &mut line,
        "{chr}\t{start}\t{chrom_end}\t{name}\t{score}\t{strand}\t\
         {thick_start}\t{thick_end}\t{rgb}\t{}\t{block_sizes}\t{block_starts}",
        blocks.len()
    )
    .unwrap();

    Some(line)
}

/// Reads a BAM file and returns the chromosome names,
/// reference sequences, and header
///
/// # Arguments
///
/// * `bam` - Path to the BAM file
///
/// # Returns
///
/// * A tuple containing:
///     * `Vec<String>` - Vector of chromosome names
///     * `Arc<noodles_sam::header::ReferenceSequences>` - Reference sequences
///     * `noodles_sam::header::Header` - Header of the BAM file
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
///
/// let bam_path = PathBuf::from("path/to/bam/file.bam");
/// let (chroms, refs, header) = read_bam(&bam_path);
///
/// assert_eq!(chroms, vec!["chr1", "chr2"]);
/// ```
fn read_bam(
    bam: &PathBuf,
) -> (
    Vec<String>,
    Arc<noodles_sam::header::ReferenceSequences>,
    Header,
) {
    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .build_from_path(bam)
        .unwrap_or_else(|e| panic!("ERROR: could not read .bam file: {e}"));

    let header = reader.read_header().unwrap();
    let refs = Arc::new(header.reference_sequences().clone()); // clone for sharing

    let chroms: Vec<String> = refs.iter().map(|(name, _)| name.to_string()).collect();

    (chroms, refs, header)
}

/// Processes a single alignment record to analyze poly-A tails and apply quality filters
///
/// This function performs several operations:
/// 1. Filters out low-quality reads based on minimum identity threshold
/// 2. Predicts poly-A tail length using HMM models
/// 3. Optionally tags reads with additional information
/// 4. Distributes reads to accept/reject accumulators based on quality criteria
///
/// # Arguments
///
/// * `record` - The alignment record to process
/// * `header` - SAM/BAM header reference
/// * `track` - Track identifier used for read tagging
/// * `chr` - Chromosome name string reference
/// * `args` - Configuration parameters for processing
/// * `hmm` - Shared HMM model for the worker processing this chromosome
/// * `outputs` - Streaming output senders for accepted/rejected reads
///
/// # Behavior
///
/// - Reads failing the minimum identity threshold (`args.min_identity`) are immediately discarded
/// - For reads with 3' end clipping:
///   - Predicts poly-A tail length using HMM if not hard-clipped
///   - Optionally adjusts prediction using suffix analysis if configured
/// - Updates the record's clipping information based on poly-A prediction
/// - Tags the read name if configured (`args.tag`)
/// - Distributes to accept/reject accumulators based on:
///   - Identity threshold (`args.identity`)
///   - Maximum 5' clip length (`args.max_clip_five`)
///   - Maximum 3' clip length (`args.max_clip_three`)
///
/// # Example
///
/// ```rust, no_run
/// use std::sync::Arc;
/// use noodles_sam::header::Header;
///
/// let record = get_alignment_record(); // Assume this gets a record
/// let header = Header::default();
/// let args = Args::default();
/// let accumulator = ParallelAccumulator::new();
///
/// process_record(
///     record,
///     &header,
///     1,
///     &"chr1".to_string(),
///     &args,
///     &accumulator
/// );
/// ```
#[allow(clippy::too_many_arguments)]
fn process_record(
    record: Record,
    header: &Header,
    track: u64,
    chr: &str,
    args: &Args,
    hmm: &HMM,
    outputs: &OutputSenders,
    singleton: bool,
    fragmented: bool,
) {
    let mut read = Read::from_mapping_record(&record);

    // INFO: enforces minimum identity to filter
    // out low quality reads from the analysis
    if read.identity < args.min_identity {
        return;
    }

    if read.three_clip > 0 || args.tail_suffix > 0 {
        read.set_sequence(&record);

        if read.three_clip > 0 && !read.has_hard_clip_three {
            let clipped_seq = read.get_rev_clipped_seq_sized();
            read.set_polya_len(predict_tail(clipped_seq, hmm));
        }

        if args.tail_suffix > 0 {
            let suffix_len = read.three_clip + args.tail_suffix;
            let mut polya_suffix_len =
                predict_tail_with_suffix(suffix_len, &read, hmm, args.suffix_step_size);

            if read.polya_len > polya_suffix_len {
                polya_suffix_len = read.polya_len;
            }

            read.set_polya_read_len(polya_suffix_len);
        }

        read.set_three_clip(read.three_clip.saturating_sub(read.polya_len));
    }

    let accepted = read.identity >= args.identity
        && read.five_clip <= args.max_clip_five
        && read.three_clip <= args.max_clip_three;

    let mut record = RecordBuf::try_from_alignment_record(header, &record)
        .unwrap_or_else(|err| panic!("ERROR: failed to convert record: {err:?}"));

    if args.tag {
        *record.name_mut() = Some(
            read.tag_read(track, chr, &args.batch, singleton, fragmented)
                .into(),
        );
    }

    outputs.send(accepted, record);
}

/// Predict the length of the polyA tail using a suffix
///
/// # Arguments
///
/// * `suffix_len` - Length of the suffix
/// * `read` - A reference to the Read object
/// * `hmm` - A reference to the HMM object
/// * `step_size` - Step size for the suffix length
///
/// # Returns
///
/// * `usize` - The length of the tail
///
/// # Example
///
/// ```rust, no_run
/// use iso::HMM;
///
/// let suffix_len = 10;
/// let read = Read::new();
/// let hmm = HMM::init(0.9, 0.99);
/// let step_size = 5;
/// let polya_len = predict_tail_with_suffix(suffix_len, &read, &hmm, step_size);
///
/// assert_eq!(polya_len, 10);
/// ```
fn predict_tail_with_suffix(
    mut suffix_len: usize,
    read: &Read,
    hmm: &HMM,
    step_size: usize,
) -> usize {
    loop {
        let suffix = read.get_rev_seq_suffix_sized(suffix_len);
        let polya_len = predict_tail(suffix, hmm);

        if polya_len == suffix_len {
            suffix_len += step_size;
        } else {
            return polya_len;
        }
    }
}

/// Predict the length of the polyA tail
/// using a two-state Hidden Markov Model (HMM)
///
/// # Arguments
///
/// * `sequence` - A vector of usize representing the sequence
///
/// # Returns
///
/// * `usize` - The length of the tail
///
/// # Example
///
/// ```rust, no_run
/// use iso::HMM;
///
/// let sequence = vec![0, 1, 2, 3, 4];
/// let polya_len = predict_tail(sequence);
///
/// assert_eq!(polya_len, 5);
/// ```
fn predict_tail(sequence: Vec<usize>, hmm: &HMM) -> usize {
    hmm.segment(&sequence)
}

/// A two-state Hidden Markov Model (HMM)
/// for polyA tail detection
///
/// # Example
///
/// ```rust, no_run
/// use iso::HMM;
///
/// let hmm = HMM::init(0.9, 0.99);
/// let sequence = vec![0, 1, 2, 3, 4];
/// let polya_len = hmm.segment(&sequence);
///
/// assert_eq!(polya_len, 5);
/// ```
#[allow(clippy::upper_case_acronyms)]
struct HMM {
    model: Model,
}

impl HMM {
    /// Initialize a two-state HMM for polyA tail detection
    ///
    /// # Arguments
    ///
    /// * `forward` - Probability of staying in state P
    /// * `weight` - Probability of emitting 'A' in state P
    ///
    /// # Returns
    ///
    /// * `HMM` - A new instance of the HMM
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let hmm = HMM::init(0.9, 0.99);
    /// let sequence = vec![0, 1, 2, 3, 4];
    /// let polya_len = hmm.segment(&sequence);
    ///
    /// assert_eq!(polya_len, 5);
    /// ```
    fn init(forward: f64, weight: f64) -> Self {
        let transition = HMM::__get_transition(forward);
        let emission = HMM::__get_emission(weight);
        let initial = array![0.5, 0.5];

        let model = Model::with_float(&transition, &emission, &initial)
            .unwrap_or_else(|err| panic!("ERROR: failed to build HMM model: {err:?}"));

        Self { model }
    }

    /// Get the transition matrix for the HMM
    ///
    /// # Arguments
    ///
    /// * `forward` - Probability of staying in state P
    ///
    /// # Returns
    ///
    /// * `Array2<f64>` - Transition matrix
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let forward = 0.9;
    /// let transition = HMM::__get_transition(forward);
    ///
    /// assert_eq!(transition, array![[0.9, 0.1], [0.0, 1.0]]);
    /// ```
    fn __get_transition(forward: f64) -> Array2<f64> {
        array![
            [forward, 1.0 - forward], // INFO: from P to P and P to S
            [0.0, 1.0]                // INFO: from S to S (no transition back to P)
        ]
    }

    /// Get the emission matrix for the HMM
    ///
    /// # Arguments
    ///
    /// * `weight` - Probability of emitting 'A' in state P
    ///
    /// # Returns
    ///
    /// * `Array2<f64>` - Emission matrix
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let weight = 0.99;
    /// let emission = HMM::__get_emission(weight);
    ///
    /// assert_eq!(emission, array![[0.99, 0.0033333333333333335,
    /// 0.0033333333333333335, 0.0033333333333333335], [0.25; 4]]);
    /// ```
    fn __get_emission(weight: f64) -> Array2<f64> {
        let uniform = 0.25;

        let emit_a = weight;
        let emit_not_a = (1.0 - emit_a) / 3.0;

        array![
            [emit_a, emit_not_a, emit_not_a, emit_not_a], // INFO: P state emissions
            [uniform; 4]                                  // INFO: S state emissions
        ]
    }

    /// Run the Viterbi algorithm on the given sequence
    ///
    /// # Arguments
    ///
    /// * `sequence` - A slice of usize representing the sequence
    ///
    /// # Returns
    ///
    /// * `Vec<State>` - The most likely state sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let hmm = HMM::init(0.9, 0.99);
    /// let sequence = vec![0, 1, 2, 3, 4];
    /// let path = hmm.__viterbi(&sequence);
    ///
    /// assert_eq!(path, vec![State(0), State(1), State(1), State(1), State(1)]);
    /// ```
    fn __viterbi(&self, sequence: &[usize]) -> Vec<State> {
        let (path, _) = viterbi(&self.model, sequence);
        path
    }

    /// Segment the sequence using the Viterbi algorithm
    ///
    /// # Arguments
    ///
    /// * `sequence` - A slice of usize representing the sequence
    ///
    /// # Returns
    ///
    /// * `usize` - The length of the tail
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let hmm = HMM::init(0.9, 0.99);
    /// let sequence = vec![0, 1, 2, 3, 4];
    /// let tail_len = hmm.segment(&sequence);
    ///
    /// assert_eq!(tail_len, 5);
    /// ```
    fn segment(&self, sequence: &[usize]) -> usize {
        let path = self.__viterbi(sequence);

        // INFO: iterate over path to find end of tail
        // INFO: state 0 corresponds to P
        let stop = path
            .iter()
            .rposition(|&state| state.0 == 0)
            .map_or(0, |pos| pos + 1);

        let tail = &sequence[..stop];
        tail.len()
    }
}

/// Representation of a BAM record in isotools
///
/// This struct is used to store the information of a BAM record
/// in a handy way for isotools.
///
/// # Example
///
/// ```rust, no_run
/// use iso::Read;
///
/// let read = Read::new();
/// read.set_name("R1");
/// read.set_strand(Strand::Forward);
/// read.set_five_clip(5);
/// read.set_three_clip(3);
///
///
/// assert_eq!(read.name, "R1");
/// assert_eq!(read.strand, Strand::Forward);
/// assert_eq!(read.five_clip, 5);
/// assert_eq!(read.three_clip, 3);
/// ```
#[derive(Debug, PartialEq, Clone)]
struct Read {
    name: String,
    strand: Strand,
    five_clip: usize,
    three_clip: usize,
    has_hard_clip_three: bool,
    sequence: Sequence,
    identity: f32,
    alignment_length: usize,
    matches: usize,
    end_site: usize,
    polya_len: usize,
    polya_read_len: usize,
}

#[allow(dead_code)]
impl Read {
    /// Creates a new `Read` instance
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// assert_eq!(read.name, "");
    /// ```
    fn new() -> Self {
        Self {
            name: String::new(),
            strand: Strand::Forward,
            five_clip: 0,
            three_clip: 0,
            has_hard_clip_three: false,
            sequence: Sequence::new(b""),
            identity: 0.0,
            alignment_length: 0,
            matches: 0,
            end_site: 0,
            polya_len: 0,
            polya_read_len: 0,
        }
    }

    /// Creates a new `Read` instance from a BAM record
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let record = Record::new();
    /// let read = Read::from(&record);
    ///
    /// assert_eq!(read.name, record.name().unwrap());
    /// ```
    fn from(record: &Record) -> Self {
        let mut read = Read::new();

        read.set_name(record);
        read.get_strand_from_record(record);
        read.set_mapping_features(record);
        read.set_sequence(record);

        read
    }

    /// Creates a Read from a mapping record (without sequence)
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    /// use noodles_sam::alignment::Record;
    ///
    /// let record = Record::new();
    /// let read = Read::from_mapping_record(&record);
    /// ```
    fn from_mapping_record(record: &Record) -> Self {
        let mut read = Read::new();
        read.get_strand_from_record(record);
        read.set_mapping_features(record);
        read
    }

    /// Get the strand from the BAM record
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let record = Record::new();
    /// let mut read = Read::from(&record);
    ///
    /// read.get_strand_from_record(&record);
    /// assert_eq!(read.strand, Strand::Forward);
    /// ```
    fn get_strand_from_record(&mut self, record: &Record) {
        let flag = record.flags();

        let strand = match flag.is_reverse_complemented() {
            true => Strand::Reverse,
            false => Strand::Forward,
        };

        self.set_strand(strand);
    }

    /// Set the strand of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_strand(Strand::Forward);
    ///
    /// assert_eq!(read.strand, Strand::Forward);
    /// ```
    fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
    }

    /// Set the polya length of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_polya_len(10);
    ///
    /// assert_eq!(read.polya_len, 10);
    /// ```
    fn set_polya_len(&mut self, polya: usize) {
        self.polya_len = polya;
    }

    /// Set the polya read length of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_polya_read_len(10);
    ///
    /// assert_eq!(read.polya_read_len, 10);
    /// ```
    fn set_polya_read_len(&mut self, polya: usize) {
        self.polya_read_len = polya;
    }

    /// Set the main mapping features of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_mapping_features(&record);
    ///
    /// assert_eq!(read.alignment_length, 100);
    /// assert_eq!(read.matches, 90);
    /// assert_eq!(read.identity, 90.0);
    /// assert_eq!(read.end_site, 200);
    /// ```
    fn set_mapping_features(&mut self, record: &Record) {
        let mut alignment_length = 0;
        let mut matches = 0;
        let mut expansion = 0;

        let mut five_clip = 0;
        let mut three_clip = 0;
        let mut hard_clip_five = false;
        let mut hard_clip_three = false;

        let mut end_site = record
            .alignment_start()
            .and_then(|start| start.ok())
            .unwrap_or_else(|| panic!("ERROR: could not get end site from read!"))
            .get();

        record.cigar().iter().for_each(|op| {
            let op = op.unwrap();

            match op.kind() {
                Kind::SequenceMatch => {
                    alignment_length += op.len();
                    matches += op.len();
                    expansion += op.len();
                }
                Kind::SequenceMismatch | Kind::Insertion | Kind::Deletion => {
                    alignment_length += op.len();
                    expansion += op.len();
                }
                Kind::HardClip | Kind::SoftClip => {
                    if alignment_length == 0 {
                        five_clip += op.len();

                        if op.kind() == Kind::HardClip {
                            hard_clip_five = true;
                        }
                    } else {
                        three_clip += op.len();

                        if op.kind() == Kind::HardClip {
                            hard_clip_three = true;
                        }
                    }
                }
                _ => {}
            }
        });

        match self.strand {
            Strand::Forward => {
                end_site += expansion - 1;

                self.five_clip = five_clip;
                self.three_clip = three_clip;
                self.has_hard_clip_three = hard_clip_three;
            }
            Strand::Reverse => {
                self.five_clip = three_clip;
                self.three_clip = five_clip;
                self.has_hard_clip_three = hard_clip_five;
            }
        }

        let identity = (matches as f32 / alignment_length as f32 * 1000.0).round() / 10.0;

        self.set_matches(matches);
        self.set_alignment_length(alignment_length);
        self.set_end_site(end_site);
        self.set_identity(identity);
    }

    /// Set matches of the read based on the CIGAR
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_matches(10);
    ///
    /// assert_eq!(read.matches, 10);
    /// ```
    fn set_matches(&mut self, matches: usize) {
        self.matches = matches;
    }

    /// Set the alignment length of the read based on the CIGAR
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_alignment_length(100);
    ///
    /// assert_eq!(read.alignment_length, 100);
    /// ```
    fn set_alignment_length(&mut self, alignment_length: usize) {
        self.alignment_length = alignment_length;
    }

    /// Set the end site of the read based on the CIGAR
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_end_site(200);
    ///
    /// assert_eq!(read.end_site, 200);
    /// ```
    fn set_end_site(&mut self, end_site: usize) {
        self.end_site = end_site;
    }

    /// Set the identity of the read based on the CIGAR
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_identity(90.0);
    ///
    /// assert_eq!(read.identity, 90.0);
    /// ```
    fn set_identity(&mut self, identity: f32) {
        self.identity = identity;
    }

    /// Set the sequence of the read from the BAM record
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let record = Record::new();
    /// let mut read = Read::from(&record);
    ///
    /// read.set_sequence(&record);
    /// assert_eq!(read.sequence, record.sequence());
    /// ```
    fn set_sequence(&mut self, record: &Record) {
        let sequence = Sequence::decode(record.sequence().as_ref());

        match self.strand {
            Strand::Forward => self.sequence = sequence,
            Strand::Reverse => self.sequence = sequence.reverse_complement(),
        }
    }

    /// Set the three clip of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_three_clip(10);
    ///
    /// assert_eq!(read.three_clip, 10);
    /// ```
    fn set_three_clip(&mut self, three_clip: usize) {
        self.three_clip = three_clip;
    }

    /// Set the name of read based on the BAM record name
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let record = Record::new();
    /// let mut read = Read::from(&record);
    ///
    /// read.set_name(&record);
    /// assert_eq!(read.name, record.name().unwrap());
    /// ```
    fn set_name(&mut self, record: &Record) {
        let name = match record.name() {
            Some(name) => name,
            None => {
                panic!("ERROR: could not get name for record (None returned)");
            }
        };
        self.name = name.to_string();
    }

    /// Get the sized reverse clipped sequence of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let rev_clipped_seq = read.get_rev_clipped_seq();
    ///
    /// assert_eq!(rev_clipped_seq, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_rev_clipped_seq_sized(&self) -> Vec<usize> {
        let seq_len = self.sequence.len();
        let three_clip = self.three_clip;

        // let start = if three_clip > seq_len {
        //     0
        // } else {
        //     seq_len - three_clip
        // };
        let start = seq_len.saturating_sub(three_clip);

        self.sequence.reverse_encode(start, seq_len)
    }

    /// Get the reverse clipped sequence of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let rev_clipped_seq = read.get_rev_clipped_seq();
    ///
    /// assert_eq!(rev_clipped_seq, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_rev_clipped_seq(&self) -> Vec<u8> {
        let seq_len = self.sequence.len();
        let three_clip = self.three_clip;

        // let start = if three_clip > seq_len {
        //     0
        // } else {
        //     seq_len - three_clip
        // };
        let start = seq_len.saturating_sub(three_clip);

        self.sequence.reverse_encode_u8(start, seq_len)
    }

    /// Get the clipped sequence of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let clipped_seq = read.get_clipped_seq();
    ///
    /// assert_eq!(clipped_seq, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_clipped_seq(&self) -> &[u8] {
        let seq_len = self.sequence.len();
        let three_clip = self.three_clip;

        assert!(
            three_clip < seq_len,
            "ERROR: suffix length is greater than sequence length -> {} : {}!",
            three_clip,
            seq_len
        );

        self.sequence.slice_as_bytes(seq_len - three_clip, seq_len)
    }

    /// Get the reverse suffix of the read sized
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let rev_suffix = read.get_rev_seq_suffix_sized(5);
    ///
    /// assert_eq!(rev_suffix, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_rev_seq_suffix_sized(&self, suffix: usize) -> Vec<usize> {
        let seq_len = self.sequence.len();

        // let start = if suffix > seq_len {
        //     0
        // } else {
        //     seq_len - suffix
        // };
        let start = seq_len.saturating_sub(suffix);
        self.sequence.reverse_encode(start, seq_len)
    }

    /// Tag read following the format:
    ///     > R[number]_[chr]::FC5:TC24:PA45:PR65:IY98
    /// where:
    ///     R[number] = read number
    ///     [chr] = chromosome name
    ///     FC = five_clip
    ///     TC = three_clip
    ///     PA = polya_len
    ///     PR = polya_read_len
    ///     IY = identity
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.tag_read();
    ///
    /// assert_eq!(read.name, "R1_chr1::FC5:TC24:PA45:PR65:IY98");
    /// ```
    fn tag_read(
        &self,
        index: u64,
        chr: &str,
        batch: &str,
        singleton: bool,
        fragmented: bool,
    ) -> String {
        let batch = if !batch.is_empty() {
            format!("@{batch}")
        } else {
            String::new()
        };

        let mut name = format!(
            "R{}{}_{chr}{BIG_SEP}FC{}{SEP}TC{}{SEP}PA{}{SEP}PR{}{SEP}IY{}",
            index,
            batch,
            self.five_clip,
            self.three_clip,
            self.polya_len,
            self.polya_read_len,
            (self.identity * 10.0) as u64
        );

        if singleton {
            name = format!("{name}{SEP}SG");
        }

        if fragmented {
            name = format!("{name}{SEP}FG");
        }

        name
    }
}

/// Strand type
///
/// This enum is used to store the strand of a sequence.
///
/// # Example
///
/// ```rust, no_run
/// use iso::Strand;
///
/// let forward = Strand::Forward;
/// let reverse = Strand::Reverse;
/// ```
#[derive(Debug, PartialEq, Clone, Hash, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Get the opposite strand
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Strand;
    ///
    /// let forward = Strand::Forward;
    /// let reverse = forward.opposite();
    ///
    /// assert_eq!(reverse, Strand::Reverse);
    /// ```
    pub fn opposite(&self) -> Self {
        match self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
        }
    }

    /// Get the string representation of the strand
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Strand;
    /// let forward = Strand::Forward;
    /// assert_eq!(forward.as_str(), "+");
    ///
    /// let reverse = Strand::Reverse;
    /// assert_eq!(reverse.as_str(), "-");
    /// ```
    pub fn as_str(&self) -> &str {
        match self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
        }
    }
}

impl std::str::FromStr for Strand {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => Err("ERROR: Cannot parse strand!".into()),
        }
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

/// Sequence struct
///
/// This struct is used to store a sequence of nucleotides.
///
/// # Example
/// ```rust, no_run
/// use iso::Sequence;
///
/// let seq = Sequence::new(b"ATCG");
/// assert_eq!(seq.len(), 4);
/// assert_eq!(seq.is_empty(), false);
/// assert_eq!(seq.as_bytes(), b"ATCG");
/// assert_eq!(seq.as_str(), "ATCG");
/// assert_eq!(seq.to_string(), "ATCG");
/// assert_eq!(seq.to_uppercase(), "ATCG");
/// assert_eq!(seq.to_lowercase(), "atcg");
/// assert_eq!(seq.reverse_complement().to_string(), "CGAT");
/// assert_eq!(seq.slice(0, 2), "AT");
/// assert_eq!(seq.slice_as_seq(0, 2).to_string(), "AT");
/// assert_eq!(seq.slice_as_bytes(0, 2), b"AT");
/// assert_eq!(seq.at_as_bytes(0), 65);
/// assert_eq!(seq.fill(2), "AATCG");
/// assert_eq!(seq.skip(1, 3).to_string(), "ACG");
/// ```
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Sequence {
    pub seq: String,
}

impl Sequence {
    /// Create a new sequence
    ///
    /// # Example
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.len(), 4);
    /// ```
    pub fn new(seq: &[u8]) -> Self {
        Self {
            seq: unsafe { std::str::from_utf8_unchecked(seq).to_string() },
        }
    }

    /// Decode a sequence from bytes
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.decode(b"ACGT"), "ACGT");
    /// ```
    pub fn decode(seq: &[u8]) -> Self {
        let base_count = seq.len() * 2;
        let mut capacity = String::with_capacity(base_count);

        for i in 0..base_count {
            let byte = seq[i / 2];
            let base_code = if i % 2 == 0 { byte >> 4 } else { byte & 0x0F };

            let base = Sequence::__decode_base(base_code);
            if base != b'=' {
                capacity.push(char::from(base));
            }
        }

        Self { seq: capacity }
    }

    /// Get the length of the sequence
    ///
    /// # Example
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.len(), 4);
    /// ```
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Check if the sequence is empty
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.is_empty(), false);
    /// ```
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Get the sequence as bytes
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.as_bytes(), b"ATCG");
    /// ```
    pub fn as_bytes(&self) -> &[u8] {
        self.seq.as_bytes()
    }

    /// Get the sequence as a string
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.as_str(), "ATCG");
    /// ```
    pub fn as_str(&self) -> &str {
        self.seq.as_str()
    }

    /// Get the sequence as a string
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.to_string(), String::from("ATCG"));
    /// ```
    #[allow(clippy::inherent_to_string_shadow_display)]
    pub fn to_string(&self) -> String {
        self.seq.clone()
    }

    /// Get the sequence as uppercase
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"atcg");
    /// assert_eq!(seq.to_uppercase(), "ATCG");
    /// ```
    pub fn to_uppercase(&self) -> String {
        self.seq.to_uppercase()
    }

    /// Get the sequence as lowercase
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.to_lowercase(), "atcg");
    /// ```
    pub fn to_lowercase(&self) -> String {
        self.seq.to_lowercase()
    }

    /// Get the complement of the sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.complement().to_string(), "GCTA");
    /// ```
    pub fn complement(&self) -> Self {
        let comp = self
            .seq
            .chars()
            .map(|c| COMPLEMENT[c as usize] as char)
            .collect::<String>();

        Self { seq: comp }
    }

    /// Get the reverse complement of the sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.reverse_complement().to_string(), "CGAT");
    /// ```
    pub fn reverse_complement(&self) -> Self {
        let mut rev = self.seq.chars().rev().collect::<String>();
        rev.make_ascii_uppercase();
        rev = rev
            .chars()
            .map(|c| COMPLEMENT[c as usize] as char)
            .collect::<String>();

        Self { seq: rev }
    }

    /// Get a slice of the sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.slice(0, 2), "AT");
    /// ```
    pub fn slice(&self, start: usize, end: usize) -> String {
        self.seq[start..end].to_string()
    }

    /// Get a slice of the sequence as a Sequence struct
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.slice_as_seq(0, 2).to_string(), "AT");
    /// ```
    pub fn slice_as_seq(&self, start: usize, end: usize) -> Self {
        Self {
            seq: self.seq[start..end].to_string(),
        }
    }

    /// Get a slice of the sequence as bytes
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.slice_as_bytes(0, 2), b"AT");
    /// ```
    pub fn slice_as_bytes(&self, start: usize, end: usize) -> &[u8] {
        &self.seq.as_bytes()[start..end]
    }

    /// Get the ASCII value of a nucleotide at a given index
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.at_as_bytes(0), 65);
    /// ```
    pub fn at_as_bytes(&self, idx: usize) -> usize {
        self.seq.as_bytes()[idx] as usize
    }

    /// Skip a given range of the sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.skip(1, 3).to_string(), "ACG");
    /// ```
    pub fn skip(&self, from: usize, to: usize) -> Sequence {
        Sequence {
            seq: self.seq[..from].to_string() + &self.seq[to..],
        }
    }

    /// Decode a base to a u8
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let base = Sequence::__decode_base(1);
    /// assert_eq!(base, b'A');
    /// ```
    pub fn __decode_base(nt: u8) -> u8 {
        match nt & 0x0f {
            0 => b'=',
            1 => b'A',
            2 => b'C',
            3 => b'M',
            4 => b'G',
            5 => b'R',
            6 => b'S',
            7 => b'V',
            8 => b'T',
            9 => b'W',
            10 => b'Y',
            11 => b'H',
            12 => b'K',
            13 => b'D',
            14 => b'B',
            15 => b'N',
            _ => panic!(
                "{}",
                format!("ERROR: invalid character in sequence: {}", nt)
            ),
        }
    }

    /// Encode a base to a u8
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let base = Sequence::__encode_base(b'A');
    /// assert_eq!(base, 1);
    /// ```
    pub fn __encode_base(nt: u8) -> u8 {
        match nt {
            b'=' => 0,
            b'A' => 1,
            b'C' => 2,
            b'M' => 3,
            b'G' => 4,
            b'R' => 5,
            b'S' => 6,
            b'V' => 7,
            b'T' => 8,
            b'W' => 9,
            b'Y' => 10,
            b'H' => 11,
            b'K' => 12,
            b'D' => 13,
            b'B' => 14,
            _ => panic!(
                "{}",
                format!("ERROR: invalid character in sequence: {}", nt)
            ),
        }
    }

    /// Encode a cannonical base to a u8
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let base = Sequence::__encode_base_2(b'A');
    /// assert_eq!(base, 0);
    /// ```
    pub fn __encode_base_2(nt: u8) -> Option<u8> {
        match nt {
            b'A' => Some(0),
            b'C' => Some(1),
            b'T' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    }

    /// Encode the reverse of a sequence to a usize
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.reverse_encode(0, 4), vec![3, 2, 1, 0]);
    /// ```
    pub fn reverse_encode(&self, start: usize, end: usize) -> Vec<usize> {
        self.slice_as_bytes(start, end)
            .iter()
            .rev()
            .filter_map(|b| Self::__encode_base_2(*b))
            .map(|nt| nt as usize)
            .collect::<Vec<usize>>()
    }

    /// Encode the reverse of a sequence to a u8
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.reverse_encode_u8(0, 4), vec![3, 2, 1, 0]);
    /// ```
    pub fn reverse_encode_u8(&self, start: usize, end: usize) -> Vec<u8> {
        self.slice_as_bytes(start, end)
            .iter()
            .rev()
            .filter_map(|b| Self::__encode_base_2(*b))
            .collect::<Vec<u8>>()
    }
}

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.seq)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_output_path_ignores_chromosome_when_split_disabled() {
        let output = build_output_path(
            &PathBuf::from("out"),
            &PathBuf::from("sample"),
            false,
            "@",
            Some("chr1"),
            PASS_PREFIX,
            "bam",
        );

        assert_eq!(output, PathBuf::from("out/sample.hq.bam"));
    }

    #[test]
    fn build_output_path_prefixes_chromosome_when_split_enabled() {
        let output = build_output_path(
            &PathBuf::from("out"),
            &PathBuf::from("sample"),
            true,
            "@",
            Some("chr10"),
            FAIL_PREFIX,
            "bed",
        );

        assert_eq!(output, PathBuf::from("out/chr10@sample.lq.bed"));
    }

    #[test]
    fn build_output_path_honors_custom_delimiter_when_split_enabled() {
        let output = build_output_path(
            &PathBuf::from("out"),
            &PathBuf::from("sample"),
            true,
            ":",
            Some("chrM"),
            PASS_PREFIX,
            "bam",
        );

        assert_eq!(output, PathBuf::from("out/chrM:sample.hq.bam"));
    }
}
