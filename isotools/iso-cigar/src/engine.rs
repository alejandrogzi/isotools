// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core modul for rescuing missed 3' splice junctions by CIGAR matching
//! Alejandro Gonzales-Irribarren, 2026
//!
//! This module contains the main function for rescuing missed 3' splice
//! junctions by CIGAR matching. It is a wrapper around the `engine` module
//! that handles the command-line interface and logging. In short, it reads
//! a BAM file and a GTF file, and writes a new BAM file with rescued junctions.
//!
//! A minimap2-aligned transcript can end exactly at an exon boundary
//! but carry a small 3' soft-clipped sequence that was not placed anywhere.
//! The hypothesis is: that clipped sequence is the start of the next downstream
//! exon, but the aligner missed the intron.
//!
//! In short, iso-cigar identifies candidate reads: those ending ±wiggle bp from a
//! known internal exon boundary, with soft-clip ≥ cutoff; retrieves the sequence of
//! the immediately downstream exon from the reference genome; checks whether the
//! soft-clipped bases match the beginning of that downstream exon; if yes: rewrites
//! the CIGAR to add an N intron gap and converts the soft-clip into a = match

use std::collections::{hash_map::Entry, BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;
use std::thread;

use anyhow::{anyhow, Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use log::{debug, info, warn};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam as sam;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::cigar::{op::Kind, Op};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::header::record::value::{
    map::{
        self, header::sort_order, header::tag as header_tag, program::tag as program_tag, Program,
    },
    Map,
};

use crate::annotation::{AnnotationIndex, Junction, Strand};
use crate::cli::Cli;
use crate::genome::{GenomeSession, GenomeSource};
use crate::logging::elapsed;

const PROGRAM_ID_PREFIX: &str = "iso-cigar";
const TAG_CORRECTED: Tag = Tag::new(b'x', b'c');
const TAG_JUNCTION: Tag = Tag::new(b'x', b'j');
const TAG_DELTA: Tag = Tag::new(b'x', b'd');
const TAG_RESCUED: Tag = Tag::new(b'x', b'r');
const TAG_INTRON: Tag = Tag::new(b'x', b'i');
const TAG_CS_LOWER: Tag = Tag::new(b'c', b's');

/// Runs the iso-cigar correction pipeline.
pub fn run(cli: Cli) -> Result<()> {
    info!(
        "[{}] start threads={} split_bam={} keep_additional={}",
        elapsed(),
        cli.threads,
        cli.split_bam,
        cli.keep_additional_corrections
    );

    fs::create_dir_all(&cli.outdir)
        .with_context(|| format!("failed to create output directory {}", cli.outdir.display()))?;

    info!(
        "[{}] load annotation {}",
        elapsed(),
        cli.annotation.display()
    );
    let annotation = Arc::new(AnnotationIndex::load(&cli.annotation)?);
    info!(
        "[{}] annotation ready contigs={} junctions={}",
        elapsed(),
        annotation.contig_count(),
        annotation.junction_count()
    );

    info!("[{}] open genome {}", elapsed(), cli.sequence.display());
    let genome_source = GenomeSource::open(&cli.sequence)?;
    info!(
        "[{}] genome ready format={}",
        elapsed(),
        genome_source.describe()
    );

    let bam_file = File::open(&cli.bam)
        .with_context(|| format!("failed to open BAM {}", cli.bam.display()))?;
    let mut reader = bam::io::Reader::new(bam_file);
    let header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header from {}", cli.bam.display()))?;

    let output_base = output_basename(&cli.bam)?;
    let extended_path = cli.outdir.join(format!("{output_base}.extended.bam"));
    let aligned_path = cli.outdir.join(format!("{output_base}.aligned.bam"));

    let mut output_header = header.clone();
    let program_id = install_program_line(&mut output_header, &cli)?;
    mark_header_unsorted(&mut output_header);

    info!("[{}] write extended {}", elapsed(), extended_path.display());
    let mut extended_writer = bam::io::Writer::new(
        File::create(&extended_path)
            .with_context(|| format!("failed to create {}", extended_path.display()))?,
    );
    extended_writer
        .write_header(&output_header)
        .with_context(|| format!("failed to write header to {}", extended_path.display()))?;

    let mut aligned_writer = if cli.split_bam {
        info!("[{}] write aligned {}", elapsed(), aligned_path.display());
        let mut writer = bam::io::Writer::new(
            File::create(&aligned_path)
                .with_context(|| format!("failed to create {}", aligned_path.display()))?,
        );
        writer
            .write_header(&output_header)
            .with_context(|| format!("failed to write header to {}", aligned_path.display()))?;
        Some(writer)
    } else {
        None
    };

    let header_reference_names: Vec<Vec<u8>> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| name.iter().copied().collect())
        .collect();
    let reference_to_annotation: Arc<Vec<Option<u32>>> = Arc::new(
        header_reference_names
            .iter()
            .map(|name| annotation.contig_id(name))
            .collect(),
    );

    let config = Config::from(cli);
    let worker_count = config.threads.max(1);
    let queue_capacity = (worker_count * 8).max(16);
    let (work_tx, work_rx) = bounded(queue_capacity);
    let (result_tx, result_rx) = bounded(queue_capacity);

    let reader_header = header.clone();
    let reader_handle = thread::Builder::new()
        .name("iso-cigar-reader".to_string())
        .spawn(move || read_records(reader, reader_header, work_tx))
        .context("failed to spawn reader thread")?;

    let mut worker_handles = Vec::with_capacity(worker_count);
    for worker_id in 0..worker_count {
        let worker_rx = work_rx.clone();
        let worker_tx = result_tx.clone();
        let worker_annotation = Arc::clone(&annotation);
        let worker_refmap = Arc::clone(&reference_to_annotation);
        let worker_genome = genome_source.clone();
        let worker_program_id = program_id.clone();
        let worker_config = config.clone();

        let handle = thread::Builder::new()
            .name(format!("iso-cigar-worker-{worker_id}"))
            .spawn(move || {
                worker_main(
                    worker_id,
                    worker_config,
                    worker_annotation,
                    worker_genome,
                    worker_program_id,
                    worker_refmap,
                    worker_rx,
                    worker_tx,
                )
            })
            .with_context(|| format!("failed to spawn worker thread {worker_id}"))?;
        worker_handles.push(handle);
    }
    drop(result_tx);
    drop(work_rx);

    info!("[{}] workers started count={}", elapsed(), worker_count);
    let stats = write_results(
        &output_header,
        &mut extended_writer,
        &mut aligned_writer,
        config.split_bam,
        result_rx,
    )?;

    join_thread(reader_handle, "reader")??;
    for (worker_id, handle) in worker_handles.into_iter().enumerate() {
        join_thread(handle, &format!("worker-{worker_id}"))??;
    }

    extended_writer.try_finish()?;
    if let Some(writer) = aligned_writer.as_mut() {
        writer.try_finish()?;
    }

    info!(
        "[{}] done input={} corrected={} unchanged={} additional={}",
        elapsed(),
        stats.total_input,
        stats.corrected_input,
        stats.unchanged_input,
        stats.additional_output
    );

    if stats.corrected_input == 0 {
        warn!("[{}] no primary alignments were corrected", elapsed());
    }

    Ok(())
}

/// Reads BAM records and sends them to the worker queue.
fn read_records<R>(
    mut reader: bam::io::Reader<R>,
    header: sam::Header,
    work_tx: Sender<WorkItem>,
) -> Result<()>
where
    R: Read,
{
    let mut record = RecordBuf::default();
    let mut index = 0_u64;

    loop {
        let bytes = reader.read_record_buf(&header, &mut record)?;
        if bytes == 0 {
            break;
        }

        work_tx
            .send(WorkItem {
                index,
                record: record.clone(),
            })
            .map_err(|_| anyhow!("worker queue closed while reading BAM"))?;
        index += 1;
    }

    info!("[{}] reader finished records={index}", elapsed());
    Ok(())
}

/// Worker thread main function for processing BAM records.
fn worker_main(
    worker_id: usize,
    config: Config,
    annotation: Arc<AnnotationIndex>,
    genome_source: GenomeSource,
    program_id: String,
    reference_to_annotation: Arc<Vec<Option<u32>>>,
    work_rx: Receiver<WorkItem>,
    result_tx: Sender<ProcessedItem>,
) -> Result<()> {
    debug!("[{}] worker-{worker_id} start", elapsed());
    let mut engine = Engine::new(config, annotation, genome_source, program_id)?;
    let mut processed = 0_u64;

    while let Ok(item) = work_rx.recv() {
        let outcome = engine.correct_record(item.record, reference_to_annotation.as_slice())?;
        result_tx
            .send(ProcessedItem {
                index: item.index,
                outcome,
            })
            .map_err(|_| anyhow!("writer queue closed while processing BAM"))?;
        processed += 1;
    }

    debug!(
        "[{}] worker-{worker_id} stop processed={processed}",
        elapsed()
    );
    Ok(())
}

/// Writes processed results to output BAM files.
fn write_results<W>(
    header: &sam::Header,
    extended_writer: &mut bam::io::Writer<W>,
    aligned_writer: &mut Option<bam::io::Writer<W>>,
    split_bam: bool,
    result_rx: Receiver<ProcessedItem>,
) -> Result<RunStats>
where
    W: Write,
{
    let mut pending = BTreeMap::new();
    let mut next_index = 0_u64;
    let mut stats = RunStats::default();

    while let Ok(item) = result_rx.recv() {
        pending.insert(item.index, item.outcome);

        while let Some(outcome) = pending.remove(&next_index) {
            stats.total_input += 1;
            match outcome {
                CorrectionOutcome::Unchanged(record) => {
                    stats.unchanged_input += 1;
                    if split_bam {
                        if let Some(writer) = aligned_writer.as_mut() {
                            writer.write_alignment_record(header, &record)?;
                        }
                    } else {
                        extended_writer.write_alignment_record(header, &record)?;
                    }
                }
                CorrectionOutcome::Corrected {
                    primary,
                    additional,
                } => {
                    stats.corrected_input += 1;
                    stats.additional_output += additional.len() as u64;
                    extended_writer.write_alignment_record(header, &primary)?;
                    for extra in additional {
                        extended_writer.write_alignment_record(header, &extra)?;
                    }
                }
            }
            next_index += 1;
        }
    }

    if !pending.is_empty() {
        return Err(anyhow!(
            "writer finished with {} out-of-order records buffered",
            pending.len()
        ));
    }

    Ok(stats)
}

/// Joins a thread and propagates panics.
fn join_thread<T>(handle: thread::JoinHandle<Result<T>>, label: &str) -> Result<Result<T>> {
    handle
        .join()
        .map_err(|_| anyhow!("{label} thread panicked"))
}

/// Configuration for the correction engine.
#[derive(Clone, Debug)]
struct Config {
    threads: usize,
    split_bam: bool,
    wiggle: u32,
    clip_cutoff: u32,
    keep_additional_corrections: bool,
}

impl From<Cli> for Config {
    /// Converts CLI args to engine configuration.
    fn from(value: Cli) -> Self {
        Self {
            threads: value.threads.max(1),
            split_bam: value.split_bam,
            wiggle: value.wiggle,
            clip_cutoff: value.clip_cutoff,
            keep_additional_corrections: value.keep_additional_corrections,
        }
    }
}

/// The correction engine for processing BAM records.
struct Engine {
    config: Config,
    annotation: Arc<AnnotationIndex>,
    genome: GenomeSession,
    program_id: String,
    exon_cache: HashMap<(usize, ExonSide), Vec<u8>>,
}

impl Engine {
    /// Creates a new correction engine.
    fn new(
        config: Config,
        annotation: Arc<AnnotationIndex>,
        genome_source: GenomeSource,
        program_id: String,
    ) -> Result<Self> {
        Ok(Self {
            config,
            annotation,
            genome: genome_source.session()?,
            program_id,
            exon_cache: HashMap::new(),
        })
    }

    /// Corrects a BAM record if possible.
    fn correct_record(
        &mut self,
        record: RecordBuf,
        reference_to_annotation: &[Option<u32>],
    ) -> Result<CorrectionOutcome> {
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            return Ok(CorrectionOutcome::Unchanged(record));
        }

        let reference_id = match record.reference_sequence_id() {
            Some(id) => id,
            None => return Ok(CorrectionOutcome::Unchanged(record)),
        };

        let contig_id = match reference_to_annotation.get(reference_id).and_then(|id| *id) {
            Some(id) => id,
            None => return Ok(CorrectionOutcome::Unchanged(record)),
        };

        let context =
            match RecordContext::from_record(&record, self.config.wiggle, self.config.clip_cutoff)?
            {
                Some(context) => context,
                None => return Ok(CorrectionOutcome::Unchanged(record)),
            };

        let candidate_junctions: Vec<_> = self
            .annotation
            .junctions_near(
                contig_id,
                context.strand,
                context.aligned_three_prime_boundary,
                self.config.wiggle,
            )
            .into_iter()
            .cloned()
            .collect();

        if candidate_junctions.is_empty() {
            return Ok(CorrectionOutcome::Unchanged(record));
        }

        let mut best_by_geometry: HashMap<CorrectionKey, CandidateCorrection> = HashMap::new();
        for junction in &candidate_junctions {
            if let Some(candidate) = self.try_candidate(&context, junction)? {
                let key = candidate.key();
                match best_by_geometry.entry(key) {
                    Entry::Vacant(slot) => {
                        slot.insert(candidate);
                    }
                    Entry::Occupied(mut slot) => {
                        if candidate.better_than(slot.get()) {
                            slot.insert(candidate);
                        }
                    }
                }
            }
        }

        if best_by_geometry.is_empty() {
            return Ok(CorrectionOutcome::Unchanged(record));
        }

        let mut candidates: Vec<CandidateCorrection> = best_by_geometry.into_values().collect();
        candidates.sort_by(|left, right| right.score_tuple().cmp(&left.score_tuple()));

        let best = candidates.remove(0);
        let primary = self.materialize_candidate(&record, &best, false)?;
        let additional = if self.config.keep_additional_corrections {
            let mut extra = Vec::with_capacity(candidates.len());
            for candidate in candidates {
                extra.push(self.materialize_candidate(&record, &candidate, true)?);
            }
            extra
        } else {
            Vec::new()
        };

        Ok(CorrectionOutcome::Corrected {
            primary,
            additional,
        })
    }

    /// Tries a candidate correction for a record and junction.
    fn try_candidate(
        &mut self,
        context: &RecordContext,
        junction: &Junction,
    ) -> Result<Option<CandidateCorrection>> {
        let delta = strand_normalized_delta(
            context.strand,
            context.aligned_three_prime_boundary,
            junction.boundary,
        );

        if delta.unsigned_abs() > u64::from(self.config.wiggle) {
            return Ok(None);
        }

        if delta > context.context_len as i64 {
            return Ok(None);
        }

        if delta < 0 && delta.unsigned_abs() as usize > context.softclip_len {
            return Ok(None);
        }

        let window_start = context.softclip_query_start - context.context_len;
        let window_end = context.softclip_query_start + context.softclip_len;
        let window = &context.transcript_query[window_start..window_end];

        if delta < 0 {
            let upstream_span = delta.unsigned_abs() as usize;
            let upstream_seq = self.exon_sequence(junction, ExonSide::Upstream)?;
            if upstream_span > upstream_seq.len() {
                return Ok(None);
            }

            let expected = &upstream_seq[upstream_seq.len() - upstream_span..];
            let observed = &window[context.context_len..context.context_len + upstream_span];
            if observed != expected {
                return Ok(None);
            }
        }

        let downstream_start = (context.context_len as i64 - delta) as usize;
        if downstream_start > window.len() {
            return Ok(None);
        }

        let downstream_seq = self.exon_sequence(junction, ExonSide::Downstream)?;
        let rescue_len = longest_common_prefix(&window[downstream_start..], downstream_seq);
        if rescue_len == 0 {
            return Ok(None);
        }

        let rescued_from_softclip = rescue_len.saturating_sub(delta.max(0) as usize);
        if rescued_from_softclip < self.config.clip_cutoff as usize {
            return Ok(None);
        }

        let adjusted_terminal_len = context.terminal_match_len as i64 - delta;
        if adjusted_terminal_len < 0 {
            return Ok(None);
        }

        let new_softclip_len = context.softclip_len as i64 - rescue_len as i64 + delta;
        if new_softclip_len < 0 {
            return Ok(None);
        }

        let mut transcript_ops = context.prefix_ops[..context.terminal_match_index].to_vec();
        if adjusted_terminal_len > 0 {
            transcript_ops.push(Op::new(Kind::Match, adjusted_terminal_len as usize));
        } else if transcript_ops.is_empty() {
            return Ok(None);
        }
        transcript_ops.push(Op::new(Kind::Skip, junction.intron_len as usize));
        transcript_ops.push(Op::new(Kind::SequenceMatch, rescue_len));
        if new_softclip_len > 0 {
            transcript_ops.push(Op::new(Kind::SoftClip, new_softclip_len as usize));
        }
        transcript_ops.extend_from_slice(&context.tail_hard_clips);
        normalize_cigar_ops(&mut transcript_ops);

        if transcript_ops.is_empty() {
            return Ok(None);
        }

        let genomic_ops = if context.strand == Strand::Reverse {
            let mut ops = transcript_ops.clone();
            ops.reverse();
            ops
        } else {
            transcript_ops
        };

        let new_start0 = if context.strand == Strand::Forward {
            context.start0
        } else {
            reverse_corrected_start(context.start0, junction.intron_len, rescue_len, delta)?
        };

        Ok(Some(CandidateCorrection {
            junction_id: junction.id,
            delta,
            intron_len: junction.intron_len,
            rescued_from_softclip,
            new_start0,
            genomic_ops,
        }))
    }

    /// Materializes a corrected BAM record from a candidate.
    fn materialize_candidate(
        &self,
        record: &RecordBuf,
        candidate: &CandidateCorrection,
        secondary: bool,
    ) -> Result<RecordBuf> {
        let mut corrected = record.clone();
        *corrected.alignment_start_mut() = Some(position_from_start0(candidate.new_start0)?);
        *corrected.cigar_mut() = candidate.genomic_ops.clone().into_iter().collect();

        let mut flags = corrected.flags();
        if secondary {
            flags.remove(Flags::SUPPLEMENTARY);
            flags.insert(Flags::SECONDARY);
        }
        *corrected.flags_mut() = flags;

        strip_stale_tags(corrected.data_mut());
        corrected
            .data_mut()
            .insert(Tag::PROGRAM, Value::from(self.program_id.clone()));
        corrected
            .data_mut()
            .insert(TAG_CORRECTED, Value::from(1_u8));
        corrected
            .data_mut()
            .insert(TAG_JUNCTION, Value::from(candidate.junction_id as u32));
        corrected
            .data_mut()
            .insert(TAG_DELTA, Value::from(candidate.delta as i32));
        corrected.data_mut().insert(
            TAG_RESCUED,
            Value::from(candidate.rescued_from_softclip as u32),
        );
        corrected
            .data_mut()
            .insert(TAG_INTRON, Value::from(candidate.intron_len as u32));

        Ok(corrected)
    }

    /// Gets the exon sequence for a junction side.
    fn exon_sequence(&mut self, junction: &Junction, side: ExonSide) -> Result<&Vec<u8>> {
        let key = (junction.id, side);
        match self.exon_cache.entry(key) {
            Entry::Occupied(entry) => Ok(entry.into_mut()),
            Entry::Vacant(entry) => {
                let exon = match side {
                    ExonSide::Upstream => junction.upstream,
                    ExonSide::Downstream => junction.downstream,
                };
                let contig = self.annotation.contig_name(junction.contig_id);
                let mut seq = self
                    .genome
                    .fetch(contig, exon.start, exon.end)
                    .with_context(|| {
                        format!(
                            "failed to fetch exon {}:{}-{}",
                            String::from_utf8_lossy(contig),
                            exon.start,
                            exon.end
                        )
                    })?;
                if junction.strand == Strand::Reverse {
                    seq = reverse_complement(&seq);
                }
                Ok(entry.insert(seq))
            }
        }
    }
}

/// Which side of a junction (upstream or downstream exon).
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum ExonSide {
    Upstream,
    Downstream,
}

/// Outcome of correcting a record.
enum CorrectionOutcome {
    Unchanged(RecordBuf),
    Corrected {
        primary: RecordBuf,
        additional: Vec<RecordBuf>,
    },
}

struct WorkItem {
    index: u64,
    record: RecordBuf,
}

struct ProcessedItem {
    index: u64,
    outcome: CorrectionOutcome,
}

#[derive(Default)]
struct RunStats {
    total_input: u64,
    corrected_input: u64,
    unchanged_input: u64,
    additional_output: u64,
}

/// Context extracted from a BAM record for correction.
#[derive(Clone, Debug)]
struct RecordContext {
    strand: Strand,
    start0: u64,
    aligned_three_prime_boundary: u64,
    transcript_query: Vec<u8>,
    prefix_ops: Vec<Op>,
    tail_hard_clips: Vec<Op>,
    terminal_match_index: usize,
    terminal_match_len: usize,
    softclip_len: usize,
    softclip_query_start: usize,
    context_len: usize,
}

impl RecordContext {
    /// Extracts record context from a BAM record.
    fn from_record(record: &RecordBuf, wiggle: u32, clip_cutoff: u32) -> Result<Option<Self>> {
        let start0 = match record.alignment_start() {
            Some(pos) => pos.get() as u64 - 1,
            None => return Ok(None),
        };

        let strand = if record.flags().is_reverse_complemented() {
            Strand::Reverse
        } else {
            Strand::Forward
        };

        let cigar_ops = record.cigar().as_ref().to_vec();
        if cigar_ops.is_empty() {
            return Ok(None);
        }

        let alignment_span = record.cigar().alignment_span() as u64;
        if alignment_span == 0 {
            return Ok(None);
        }

        let aligned_three_prime_boundary = if strand == Strand::Forward {
            start0 + alignment_span
        } else {
            start0
        };

        let mut transcript_ops = cigar_ops;
        if strand == Strand::Reverse {
            transcript_ops.reverse();
        }

        let mut hard_clip_count = 0usize;
        while hard_clip_count < transcript_ops.len()
            && transcript_ops[transcript_ops.len() - 1 - hard_clip_count].kind() == Kind::HardClip
        {
            hard_clip_count += 1;
        }

        let softclip_index = match transcript_ops.len().checked_sub(hard_clip_count + 1) {
            Some(index) => index,
            None => return Ok(None),
        };
        let softclip_op = transcript_ops[softclip_index];
        if softclip_op.kind() != Kind::SoftClip || softclip_op.len() < clip_cutoff as usize {
            return Ok(None);
        }
        if softclip_index == 0 {
            return Ok(None);
        }

        let prefix_ops = transcript_ops[..softclip_index].to_vec();
        let tail_hard_clips = transcript_ops[softclip_index + 1..].to_vec();
        let terminal_match_index = prefix_ops.len() - 1;
        let terminal_op = prefix_ops[terminal_match_index];
        if !is_match_like(terminal_op.kind()) {
            return Ok(None);
        }

        let mut query_consumed = 0usize;
        for op in &transcript_ops[..softclip_index] {
            if op.kind().consumes_read() {
                query_consumed += op.len();
            }
        }
        let softclip_query_start = query_consumed;

        let mut clean_query_bases = 0usize;
        for op in prefix_ops.iter().rev() {
            if !is_match_like(op.kind()) {
                if clean_query_bases < wiggle as usize {
                    return Ok(None);
                }
                break;
            }

            clean_query_bases += op.len();
            if clean_query_bases >= wiggle as usize {
                break;
            }
        }

        let context_len = clean_query_bases
            .min(wiggle as usize)
            .min(softclip_query_start);
        if context_len == 0 {
            return Ok(None);
        }

        let sequence = record.sequence().as_ref();
        let transcript_query = if strand == Strand::Forward {
            let mut query = sequence.to_vec();
            query.make_ascii_uppercase();
            query
        } else {
            reverse_complement(sequence)
        };

        Ok(Some(Self {
            strand,
            start0,
            aligned_three_prime_boundary,
            transcript_query,
            prefix_ops,
            tail_hard_clips,
            terminal_match_index,
            terminal_match_len: terminal_op.len(),
            softclip_len: softclip_op.len(),
            softclip_query_start,
            context_len,
        }))
    }
}

/// A candidate correction for a BAM record.
#[derive(Clone, Debug)]
struct CandidateCorrection {
    junction_id: usize,
    delta: i64,
    intron_len: u64,
    rescued_from_softclip: usize,
    new_start0: u64,
    genomic_ops: Vec<Op>,
}

impl CandidateCorrection {
    /// Returns a key for deduplication.
    fn key(&self) -> CorrectionKey {
        CorrectionKey {
            start0: self.new_start0,
            cigar: cigar_key(&self.genomic_ops),
        }
    }

    /// Compares if this candidate is better than another.
    fn better_than(&self, other: &Self) -> bool {
        self.score_tuple() > other.score_tuple()
    }

    /// Returns scoring tuple for comparison.
    fn score_tuple(&self) -> (usize, i64, i64) {
        (
            self.rescued_from_softclip,
            -(self.delta.unsigned_abs() as i64),
            -(self.junction_id as i64),
        )
    }
}

/// Key for deduplicating corrections (start + CIGAR).
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct CorrectionKey {
    start0: u64,
    cigar: Vec<(u8, u32)>,
}

/// Removes stale tags from corrected records.
fn strip_stale_tags(data: &mut sam::alignment::record_buf::Data) {
    data.remove(&Tag::EDIT_DISTANCE);
    data.remove(&Tag::MISMATCHED_POSITIONS);
    data.remove(&Tag::ALIGNMENT_SCORE);
    data.remove(&Tag::OTHER_ALIGNMENTS);
    data.remove(&Tag::CIGAR);
    data.remove(&Tag::ORIGINAL_ALIGNMENT);
    data.remove(&Tag::ORIGINAL_CIGAR);
    data.remove(&Tag::ORIGINAL_POSITION);
    data.remove(&TAG_CS_LOWER);
}

/// Extracts basename from a BAM path.
fn output_basename(path: &Path) -> Result<String> {
    path.file_stem()
        .and_then(|stem| stem.to_str())
        .map(|stem| stem.to_owned())
        .ok_or_else(|| anyhow!("cannot determine BAM basename for {}", path.display()))
}

/// Installs program line in BAM header.
fn install_program_line(header: &mut sam::Header, cli: &Cli) -> Result<String> {
    let mut id = PROGRAM_ID_PREFIX.to_string();
    let programs = header.programs_mut().as_mut();
    let mut suffix = 1usize;

    while programs.contains_key(id.as_bytes()) {
        suffix += 1;
        id = format!("{PROGRAM_ID_PREFIX}-{suffix}");
    }

    let cmdline = format!(
        "{} --bam {} --annotation {} --sequence {} --outdir {}{}{} --threads {} --wiggle {} --clip-cutoff {} --level {}",
        PROGRAM_ID_PREFIX,
        cli.bam.display(),
        cli.annotation.display(),
        cli.sequence.display(),
        cli.outdir.display(),
        if cli.split_bam { " --split-bam" } else { "" },
        if cli.keep_additional_corrections {
            " --keep-additional-corrections"
        } else {
            ""
        },
        cli.threads,
        cli.wiggle,
        cli.clip_cutoff,
        match cli.level {
            crate::cli::LogLevel::Trace => "trace",
            crate::cli::LogLevel::Debug => "debug",
            crate::cli::LogLevel::Info => "info",
            crate::cli::LogLevel::Warn => "warn",
            crate::cli::LogLevel::Error => "error",
            crate::cli::LogLevel::Off => "off",
        }
    );

    let program = Map::<Program>::builder()
        .insert(program_tag::NAME, PROGRAM_ID_PREFIX)
        .insert(program_tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(program_tag::COMMAND_LINE, cmdline)
        .build()?;
    programs.insert(id.clone().into(), program);
    Ok(id)
}

/// Marks BAM header as unsorted.
fn mark_header_unsorted(header: &mut sam::Header) {
    let hd = header
        .header_mut()
        .get_or_insert_with(Map::<map::Header>::default);
    hd.other_fields_mut()
        .insert(header_tag::SORT_ORDER, sort_order::UNSORTED.into());
}

/// Calculates corrected start position for reverse strand.
fn reverse_corrected_start(
    start0: u64,
    intron_len: u64,
    rescue_len: usize,
    delta: i64,
) -> Result<u64> {
    let shifted = start0
        .checked_sub(intron_len)
        .and_then(|value| value.checked_sub(rescue_len as u64))
        .ok_or_else(|| anyhow!("corrected reverse-strand position underflow"))?;

    if delta >= 0 {
        shifted
            .checked_add(delta as u64)
            .ok_or_else(|| anyhow!("corrected reverse-strand position overflow"))
    } else {
        shifted
            .checked_sub(delta.unsigned_abs())
            .ok_or_else(|| anyhow!("corrected reverse-strand position underflow"))
    }
}

/// Returns strand-normalized delta between aligned and junction boundary.
fn strand_normalized_delta(strand: Strand, aligned_boundary: u64, junction_boundary: u64) -> i64 {
    match strand {
        Strand::Forward => aligned_boundary as i64 - junction_boundary as i64,
        Strand::Reverse => junction_boundary as i64 - aligned_boundary as i64,
    }
}

/// Returns the length of the longest common prefix.
fn longest_common_prefix(left: &[u8], right: &[u8]) -> usize {
    left.iter()
        .zip(right.iter())
        .take_while(|(a, b)| a == b)
        .count()
}

/// Returns true if the CIGAR operation is match-like.
fn is_match_like(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

/// Normalizes CIGAR operations by merging adjacent same-kind ops.
fn normalize_cigar_ops(ops: &mut Vec<Op>) {
    ops.retain(|op| !op.is_empty());
    if ops.is_empty() {
        return;
    }

    let mut normalized: Vec<Op> = Vec::with_capacity(ops.len());
    for op in ops.drain(..) {
        if let Some(last) = normalized.last_mut() {
            if last.kind() == op.kind() {
                *last = Op::new(op.kind(), last.len() + op.len());
                continue;
            }
        }
        normalized.push(op);
    }

    *ops = normalized;
}

/// Converts 0-based start to SAM Position.
fn position_from_start0(start0: u64) -> Result<Position> {
    Position::try_from(start0 as usize + 1)
        .map_err(|_| anyhow!("alignment start {} is out of SAM range", start0))
}

/// Returns the reverse complement of a sequence.
fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    let mut rev = Vec::with_capacity(sequence.len());
    for base in sequence.iter().rev() {
        rev.push(complement(*base));
    }
    rev
}

/// Returns the complement of a base.
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        b'U' | b'u' => b'A',
        b'R' | b'r' => b'Y',
        b'Y' | b'y' => b'R',
        b'S' | b's' => b'S',
        b'W' | b'w' => b'W',
        b'K' | b'k' => b'M',
        b'M' | b'm' => b'K',
        b'B' | b'b' => b'V',
        b'D' | b'd' => b'H',
        b'H' | b'h' => b'D',
        b'V' | b'v' => b'B',
        _ => b'N',
    }
}

/// Returns a key for deduplication from CIGAR ops.
fn cigar_key(ops: &[Op]) -> Vec<(u8, u32)> {
    ops.iter()
        .map(|op| (kind_code(op.kind()), op.len() as u32))
        .collect()
}

/// Converts CIGAR kind to SAM character code.
fn kind_code(kind: Kind) -> u8 {
    match kind {
        Kind::Match => b'M',
        Kind::Insertion => b'I',
        Kind::Deletion => b'D',
        Kind::Skip => b'N',
        Kind::SoftClip => b'S',
        Kind::HardClip => b'H',
        Kind::Pad => b'P',
        Kind::SequenceMatch => b'=',
        Kind::SequenceMismatch => b'X',
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::Path;
    use std::sync::Arc;

    use super::*;
    use noodles_sam::alignment::record_buf::{Cigar, Sequence};
    use tempfile::tempdir;

    fn make_record(flags: Flags, start1: usize, cigar: &[Op], seq: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start1).unwrap())
            .set_cigar(Cigar::from(cigar.to_vec()))
            .set_sequence(Sequence::from(seq.to_vec()))
            .build()
    }

    fn set_bases(target: &mut [u8], start: usize, bases: &[u8]) {
        target[start..start + bases.len()].copy_from_slice(bases);
    }

    fn write_fasta(path: &Path, contig: &str, sequence: &[u8]) {
        let mut fasta = Vec::new();
        fasta.extend_from_slice(b">");
        fasta.extend_from_slice(contig.as_bytes());
        fasta.extend_from_slice(b"\n");
        fasta.extend_from_slice(sequence);
        fasta.extend_from_slice(b"\n");
        fs::write(path, fasta).unwrap();
    }

    fn cigar_signature(record: &RecordBuf) -> Vec<(u8, usize)> {
        record
            .cigar()
            .as_ref()
            .iter()
            .map(|op| (kind_code(op.kind()), op.len()))
            .collect()
    }

    fn test_engine(annotation_path: &Path, genome_path: &Path) -> Engine {
        let annotation = Arc::new(AnnotationIndex::load(annotation_path).unwrap());
        let genome = GenomeSource::open(genome_path).unwrap();
        Engine::new(
            Config {
                threads: 1,
                split_bam: false,
                wiggle: 2,
                clip_cutoff: 5,
                keep_additional_corrections: false,
            },
            annotation,
            genome,
            "iso-cigar".to_string(),
        )
        .unwrap()
    }

    #[test]
    fn record_context_detects_forward_softclip_tail() {
        let record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 5)],
            b"ACGTACGTTTTTT",
        );

        let context = RecordContext::from_record(&record, 4, 5).unwrap().unwrap();
        assert_eq!(context.strand, Strand::Forward);
        assert_eq!(context.softclip_len, 5);
        assert_eq!(context.context_len, 4);
        assert_eq!(context.aligned_three_prime_boundary, 108);
    }

    #[test]
    fn record_context_detects_reverse_softclip_tail() {
        let record = make_record(
            Flags::REVERSE_COMPLEMENTED,
            201,
            &[Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 8)],
            b"AAAAAACGTACGT",
        );

        let context = RecordContext::from_record(&record, 4, 5).unwrap().unwrap();
        assert_eq!(context.strand, Strand::Reverse);
        assert_eq!(context.softclip_len, 5);
        assert_eq!(context.context_len, 4);
        assert_eq!(context.aligned_three_prime_boundary, 200);
    }

    #[test]
    fn normalize_merges_adjacent_ops() {
        let mut ops = vec![
            Op::new(Kind::Match, 4),
            Op::new(Kind::Match, 3),
            Op::new(Kind::SoftClip, 0),
        ];
        normalize_cigar_ops(&mut ops);
        assert_eq!(ops, vec![Op::new(Kind::Match, 7)]);
    }

    #[test]
    fn reverse_complement_uppercases_unknowns() {
        assert_eq!(reverse_complement(b"acgtn"), b"NACGT");
    }

    #[test]
    fn strand_delta_is_normalized() {
        assert_eq!(strand_normalized_delta(Strand::Forward, 120, 125), -5);
        assert_eq!(strand_normalized_delta(Strand::Reverse, 120, 125), 5);
    }

    #[test]
    fn candidate_scoring_prefers_more_rescue_then_smaller_delta() {
        let a = CandidateCorrection {
            junction_id: 3,
            delta: 0,
            intron_len: 10,
            rescued_from_softclip: 6,
            new_start0: 10,
            genomic_ops: vec![Op::new(Kind::Match, 10)],
        };
        let b = CandidateCorrection {
            junction_id: 2,
            delta: 2,
            intron_len: 10,
            rescued_from_softclip: 5,
            new_start0: 10,
            genomic_ops: vec![Op::new(Kind::Match, 10)],
        };
        assert!(a.better_than(&b));
    }

    #[test]
    fn correct_record_rescues_forward_exact_boundary() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 200];
        set_bases(&mut genome, 100, b"ACGTACGT");
        set_bases(&mut genome, 118, b"TTGCAAAAAA");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            b"chr1\t100\t128\ttx1\t0\t+\t100\t128\t0,0,0\t2\t8,10\t0,18\n",
        )
        .unwrap();

        let mut engine = test_engine(&annotation_path, &genome_path);
        let record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 5)],
            b"ACGTACGTTTGCA",
        );

        let outcome = engine.correct_record(record, &[Some(0)]).unwrap();
        let corrected = match outcome {
            CorrectionOutcome::Corrected { primary, .. } => primary,
            CorrectionOutcome::Unchanged(_) => panic!("expected corrected record"),
        };

        assert_eq!(corrected.alignment_start().unwrap().get(), 101);
        assert_eq!(
            cigar_signature(&corrected),
            vec![(b'M', 8), (b'N', 10), (b'=', 5)]
        );
        assert_eq!(
            corrected.data().get(&TAG_RESCUED),
            Some(&Value::from(5_u32))
        );
    }

    #[test]
    fn correct_record_rescues_reverse_exact_boundary() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 200];
        set_bases(&mut genome, 100, b"TGCAT");
        set_bases(&mut genome, 120, b"ACGTACGT");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            b"chr1\t100\t128\ttx1\t0\t-\t100\t128\t0,0,0\t2\t5,8\t0,20\n",
        )
        .unwrap();

        let mut engine = test_engine(&annotation_path, &genome_path);
        let record = make_record(
            Flags::REVERSE_COMPLEMENTED,
            121,
            &[Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 8)],
            b"TGCATACGTACGT",
        );

        let outcome = engine.correct_record(record, &[Some(0)]).unwrap();
        let corrected = match outcome {
            CorrectionOutcome::Corrected { primary, .. } => primary,
            CorrectionOutcome::Unchanged(_) => panic!("expected corrected record"),
        };

        assert_eq!(corrected.alignment_start().unwrap().get(), 101);
        assert_eq!(
            cigar_signature(&corrected),
            vec![(b'=', 5), (b'N', 15), (b'M', 8)]
        );
    }

    #[test]
    fn correct_record_validates_upstream_extension_and_strips_stale_tags() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 200];
        set_bases(&mut genome, 100, b"ACGTACGT");
        set_bases(&mut genome, 118, b"TTGCAAAAAA");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            b"chr1\t100\t128\ttx1\t0\t+\t100\t128\t0,0,0\t2\t8,10\t0,18\n",
        )
        .unwrap();

        let mut engine = test_engine(&annotation_path, &genome_path);
        let mut record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 6), Op::new(Kind::SoftClip, 7)],
            b"ACGTACGTTTGCA",
        );
        record
            .data_mut()
            .insert(Tag::EDIT_DISTANCE, Value::from(3_u8));
        record
            .data_mut()
            .insert(Tag::MISMATCHED_POSITIONS, Value::from("6A1"));

        let outcome = engine.correct_record(record, &[Some(0)]).unwrap();
        let corrected = match outcome {
            CorrectionOutcome::Corrected { primary, .. } => primary,
            CorrectionOutcome::Unchanged(_) => panic!("expected corrected record"),
        };

        assert_eq!(
            cigar_signature(&corrected),
            vec![(b'M', 8), (b'N', 10), (b'=', 5)]
        );
        assert_eq!(corrected.data().get(&TAG_DELTA), Some(&Value::from(-2_i32)));
        assert!(corrected.data().get(&Tag::EDIT_DISTANCE).is_none());
        assert!(corrected.data().get(&Tag::MISMATCHED_POSITIONS).is_none());
    }
}
