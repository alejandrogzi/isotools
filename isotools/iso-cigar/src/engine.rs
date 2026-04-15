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
//! exon (or any other downstream), but the aligner missed the intron.
//!
//! In short, iso-cigar identifies candidate reads: those ending ±wiggle bp from a
//! known internal exon boundary, with soft-clip ≥ cutoff; retrieves the sequence of
//! the immediately downstream exon from the reference genome; checks whether the
//! soft-clipped bases match the beginning of that downstream exon; if yes: rewrites
//! the CIGAR to add an N intron gap and converts the soft-clip into a = match

use std::cmp::Reverse;
use std::collections::{hash_map::Entry, BTreeMap, HashMap};
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;

use anyhow::{anyhow, Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use log::{debug, info, warn};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_csi::binning_index::index::reference_sequence::{bin::Chunk, index::LinearIndex, Bin};
use noodles_csi::binning_index::index::ReferenceSequence as BinningReferenceSequence;
use noodles_sam as sam;
use noodles_sam::alignment::record::cigar::{op::Kind, Op};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::{io::Write as _, Record as _};
use noodles_sam::header::record::value::{
    map::{
        self, header::sort_order, header::tag as header_tag, program::tag as program_tag, Program,
    },
    Map,
};

use crate::annotation::{
    intron_len_between, AnnotationIndex, Junction, Strand, TranscriptJunctionRef, TranscriptModel,
};
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
        "[{}] start threads={} split_bam={} keep_additional={} extend_alignment={}",
        elapsed(),
        cli.threads,
        cli.split_bam,
        cli.keep_additional_corrections,
        cli.extend_alignment
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

    let config = Config::from(cli.clone());
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

    write_bam_index(&extended_path)?;
    if cli.split_bam {
        write_bam_index(&aligned_path)?;
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

        // INFO: remove extended if --split-bam is set
        if cli.split_bam {
            info!(
                "[{}] removing extended BAM {} of size 0",
                elapsed(),
                extended_path.display()
            );

            remove_bam_and_index(&extended_path)?;
        }
    }

    info!("[{}] done!", elapsed());

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
    extend_alignment: bool,
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
            extend_alignment: value.extend_alignment,
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
    transcript_exon_cache: HashMap<(usize, usize), Vec<u8>>,
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
            transcript_exon_cache: HashMap::new(),
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
            for candidate in self.try_candidates(&context, junction)? {
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

    /// Tries candidate corrections for a record and starting junction.
    fn try_candidates(
        &mut self,
        context: &RecordContext,
        junction: &Junction,
    ) -> Result<Vec<CandidateCorrection>> {
        let delta = strand_normalized_delta(
            context.strand,
            context.aligned_three_prime_boundary,
            junction.boundary,
        );

        if delta.unsigned_abs() > u64::from(self.config.wiggle) {
            return Ok(Vec::new());
        }

        if delta > context.context_len as i64 {
            return Ok(Vec::new());
        }

        if delta < 0 && delta.unsigned_abs() as usize > context.softclip_len {
            return Ok(Vec::new());
        }

        let window_start = context.softclip_query_start - context.context_len;
        let window_end = context.softclip_query_start + context.softclip_len;
        let window = &context.transcript_query[window_start..window_end];

        if delta < 0 {
            let upstream_span = delta.unsigned_abs() as usize;
            let upstream_seq = self.exon_sequence(junction, ExonSide::Upstream)?;
            if upstream_span > upstream_seq.len() {
                return Ok(Vec::new());
            }

            let expected = &upstream_seq[upstream_seq.len() - upstream_span..];
            let observed = &window[context.context_len..context.context_len + upstream_span];
            if observed != expected {
                return Ok(Vec::new());
            }
        }

        let downstream_start = (context.context_len as i64 - delta) as usize;
        if downstream_start > window.len() {
            return Ok(Vec::new());
        }

        let (first_rescue_len, first_exon_len) = {
            let downstream_seq = self.exon_sequence(junction, ExonSide::Downstream)?;
            (
                longest_common_prefix(&window[downstream_start..], downstream_seq),
                downstream_seq.len(),
            )
        };
        if first_rescue_len == 0 {
            return Ok(Vec::new());
        }

        let base_path = RescuePath::single(junction.intron_len, first_rescue_len);
        if !self.config.extend_alignment
            || first_rescue_len < first_exon_len
            || downstream_start + first_rescue_len >= window.len()
        {
            return Ok(self
                .build_candidate(context, junction.id, junction.id, delta, &base_path)?
                .into_iter()
                .collect());
        }

        let transcript_refs: Vec<_> = self.annotation.transcript_refs(junction.id).to_vec();
        if transcript_refs.is_empty() {
            return Ok(self
                .build_candidate(context, junction.id, junction.id, delta, &base_path)?
                .into_iter()
                .collect());
        }

        let mut candidates = Vec::with_capacity(transcript_refs.len());
        for transcript_ref in transcript_refs {
            let Some(path) = self.walk_transcript_rescue(
                window,
                downstream_start,
                junction,
                transcript_ref,
                first_rescue_len,
                first_exon_len,
            )?
            else {
                continue;
            };

            if let Some(candidate) = self.build_candidate(
                context,
                junction.id,
                transcript_ref.transcript_id,
                delta,
                &path,
            )? {
                candidates.push(candidate);
            }
        }

        if candidates.is_empty() {
            Ok(self
                .build_candidate(context, junction.id, junction.id, delta, &base_path)?
                .into_iter()
                .collect())
        } else {
            Ok(candidates)
        }
    }

    /// Walks downstream exons for a transcript-aware rescue path.
    fn walk_transcript_rescue(
        &mut self,
        window: &[u8],
        downstream_start: usize,
        junction: &Junction,
        transcript_ref: TranscriptJunctionRef,
        first_rescue_len: usize,
        first_exon_len: usize,
    ) -> Result<Option<RescuePath>> {
        let Some(transcript) = self
            .annotation
            .transcript(transcript_ref.transcript_id)
            .cloned()
        else {
            return Ok(None);
        };

        if transcript.contig_id != junction.contig_id || transcript.strand != junction.strand {
            return Ok(None);
        }

        let Some(&upstream_exon) = transcript.exons.get(transcript_ref.upstream_exon_index) else {
            return Ok(None);
        };
        let Some(&first_downstream_exon) =
            transcript.exons.get(transcript_ref.upstream_exon_index + 1)
        else {
            return Ok(None);
        };

        if upstream_exon != junction.upstream || first_downstream_exon != junction.downstream {
            return Ok(None);
        }

        let mut path = RescuePath::single(junction.intron_len, first_rescue_len);
        if first_rescue_len < first_exon_len {
            return Ok(Some(path));
        }

        let mut query_pos = downstream_start + first_rescue_len;
        let mut exon_index = transcript_ref.upstream_exon_index + 2;
        while exon_index < transcript.exons.len() && query_pos < window.len() {
            let upstream = transcript.exons[exon_index - 1];
            let downstream = transcript.exons[exon_index];
            let Some(intron_len) = intron_len_between(transcript.strand, upstream, downstream)
            else {
                break;
            };

            let exon_seq = self.transcript_exon_sequence(&transcript, exon_index)?;
            let rescue_len = longest_common_prefix(&window[query_pos..], exon_seq);
            if rescue_len == 0 {
                break;
            }

            path.push(intron_len, rescue_len)?;
            query_pos += rescue_len;

            if rescue_len < exon_seq.len() {
                break;
            }

            exon_index += 1;
        }

        Ok(Some(path))
    }

    /// Builds a concrete corrected candidate from a rescue path.
    fn build_candidate(
        &self,
        context: &RecordContext,
        junction_id: usize,
        tie_breaker: usize,
        delta: i64,
        path: &RescuePath,
    ) -> Result<Option<CandidateCorrection>> {
        let rescued_from_softclip = path.rescued_total.saturating_sub(delta.max(0) as usize);
        if rescued_from_softclip < self.config.clip_cutoff as usize {
            return Ok(None);
        }

        let adjusted_terminal_len = context.terminal_match_len as i64 - delta;
        if adjusted_terminal_len < 0 {
            return Ok(None);
        }

        let new_softclip_len = context.softclip_len as i64 - path.rescued_total as i64 + delta;
        if new_softclip_len < 0 {
            return Ok(None);
        }

        let mut transcript_ops = context.prefix_ops[..context.terminal_match_index].to_vec();
        if adjusted_terminal_len > 0 {
            transcript_ops.push(Op::new(Kind::Match, adjusted_terminal_len as usize));
        } else if transcript_ops.is_empty() {
            return Ok(None);
        }

        for segment in &path.segments {
            transcript_ops.push(Op::new(Kind::Skip, segment.intron_len as usize));
            transcript_ops.push(Op::new(Kind::SequenceMatch, segment.rescue_len));
        }

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

        let total_intron_len = path.total_intron_len()?;
        let total_reference_extension = total_intron_len
            .checked_add(path.rescued_total as u64)
            .ok_or_else(|| anyhow!("rescued reference span overflow"))?;
        let new_start0 = if context.strand == Strand::Forward {
            context.start0
        } else {
            reverse_corrected_start(context.start0, total_reference_extension, delta)?
        };

        Ok(Some(CandidateCorrection {
            junction_id,
            delta,
            intron_len: total_intron_len,
            rescued_from_softclip,
            rescued_exons: path.rescued_exons,
            tie_breaker,
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
        if !self.exon_cache.contains_key(&key) {
            let exon = match side {
                ExonSide::Upstream => junction.upstream,
                ExonSide::Downstream => junction.downstream,
            };
            let seq = self.fetch_exon_sequence(
                junction.contig_id,
                junction.strand,
                exon,
                "junction exon",
            )?;
            self.exon_cache.insert(key, seq);
        }
        Ok(self.exon_cache.get(&key).unwrap())
    }

    /// Gets a transcript exon sequence in transcript orientation.
    fn transcript_exon_sequence(
        &mut self,
        transcript: &TranscriptModel,
        exon_index: usize,
    ) -> Result<&Vec<u8>> {
        let key = (transcript.id, exon_index);
        if !self.transcript_exon_cache.contains_key(&key) {
            let exon = transcript.exons[exon_index];
            let seq = self.fetch_exon_sequence(
                transcript.contig_id,
                transcript.strand,
                exon,
                "transcript exon",
            )?;
            self.transcript_exon_cache.insert(key, seq);
        }
        Ok(self.transcript_exon_cache.get(&key).unwrap())
    }

    fn fetch_exon_sequence(
        &mut self,
        contig_id: u32,
        strand: Strand,
        exon: crate::annotation::Exon,
        label: &str,
    ) -> Result<Vec<u8>> {
        let contig = self.annotation.contig_name(contig_id);
        let mut seq = self
            .genome
            .fetch(contig, exon.start, exon.end)
            .with_context(|| {
                format!(
                    "failed to fetch {label} {}:{}-{}",
                    String::from_utf8_lossy(contig),
                    exon.start,
                    exon.end
                )
            })?;
        if strand == Strand::Reverse {
            seq = reverse_complement(&seq);
        }
        Ok(seq)
    }
}

/// Which side of a junction (upstream or downstream exon).
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum ExonSide {
    Upstream,
    Downstream,
}

#[derive(Clone, Debug)]
struct RescueSegment {
    intron_len: u64,
    rescue_len: usize,
}

#[derive(Clone, Debug)]
struct RescuePath {
    segments: Vec<RescueSegment>,
    rescued_total: usize,
    rescued_exons: usize,
}

impl RescuePath {
    fn single(intron_len: u64, rescue_len: usize) -> Self {
        Self {
            segments: vec![RescueSegment {
                intron_len,
                rescue_len,
            }],
            rescued_total: rescue_len,
            rescued_exons: 1,
        }
    }

    fn push(&mut self, intron_len: u64, rescue_len: usize) -> Result<()> {
        self.rescued_total = self
            .rescued_total
            .checked_add(rescue_len)
            .ok_or_else(|| anyhow!("rescued length overflow"))?;
        self.rescued_exons = self
            .rescued_exons
            .checked_add(1)
            .ok_or_else(|| anyhow!("rescued exon count overflow"))?;
        self.segments.push(RescueSegment {
            intron_len,
            rescue_len,
        });
        Ok(())
    }

    fn total_intron_len(&self) -> Result<u64> {
        self.segments.iter().try_fold(0_u64, |sum, segment| {
            sum.checked_add(segment.intron_len)
                .ok_or_else(|| anyhow!("rescued intron length overflow"))
        })
    }
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
    rescued_exons: usize,
    tie_breaker: usize,
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
    fn score_tuple(&self) -> (usize, usize, Reverse<u64>, Reverse<usize>) {
        (
            self.rescued_from_softclip,
            self.rescued_exons,
            Reverse(self.delta.unsigned_abs()),
            Reverse(self.tie_breaker),
        )
    }
}

/// Key for deduplicating corrections (start + CIGAR).
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct CorrectionKey {
    start0: u64,
    cigar: Vec<(u8, u32)>,
}

const BAI_MIN_SHIFT: u8 = 14;
const BAI_DEPTH: u8 = 5;

#[derive(Clone, Default)]
struct ReferenceSequenceIndexBuilder {
    bins: BTreeMap<usize, Vec<Chunk>>,
}

impl ReferenceSequenceIndexBuilder {
    fn add_record(&mut self, start: Position, end: Position, chunk: Chunk) {
        let bin_id = reg2bin(start, end, BAI_MIN_SHIFT, BAI_DEPTH);
        push_chunk(self.bins.entry(bin_id).or_default(), chunk);
    }

    fn build(self) -> BinningReferenceSequence<LinearIndex> {
        let bins = self
            .bins
            .into_iter()
            .map(|(id, chunks)| (id, Bin::new(chunks)))
            .collect();
        BinningReferenceSequence::new(bins, LinearIndex::default(), None)
    }
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

/// Returns the associated BAM index path (`<bam>.bai`).
fn bam_index_path(path: &Path) -> PathBuf {
    let mut value = OsString::from(path.as_os_str());
    value.push(".bai");
    PathBuf::from(value)
}

/// Builds and writes a BAI for a BAM file.
fn write_bam_index(path: &Path) -> Result<()> {
    info!(
        "[{}] write index {}",
        elapsed(),
        bam_index_path(path).display()
    );
    let index = build_bam_index(path)?;
    bam::bai::fs::write(bam_index_path(path), &index)
        .with_context(|| format!("failed to write BAM index for {}", path.display()))?;
    Ok(())
}

/// Removes a BAM and its associated index if present.
fn remove_bam_and_index(path: &Path) -> Result<()> {
    fs::remove_file(path).with_context(|| format!("failed to remove BAM {}", path.display()))?;

    let index_path = bam_index_path(path);
    if index_path.exists() {
        fs::remove_file(&index_path)
            .with_context(|| format!("failed to remove BAM index {}", index_path.display()))?;
    }

    Ok(())
}

/// Builds a BAM index from an output BAM without assuming coordinate sort order.
fn build_bam_index(path: &Path) -> Result<bam::bai::Index> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .with_context(|| format!("failed to open BAM for indexing {}", path.display()))?;
    let header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header from {}", path.display()))?;

    let mut reference_sequences =
        vec![ReferenceSequenceIndexBuilder::default(); header.reference_sequences().len()];
    let mut unplaced_unmapped_record_count = 0_u64;
    let mut record = bam::Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        match record_alignment_context(&record)? {
            Some((reference_sequence_id, start, end, _)) => {
                if let Some(reference_sequence) = reference_sequences.get_mut(reference_sequence_id)
                {
                    reference_sequence.add_record(start, end, chunk);
                } else {
                    return Err(anyhow!(
                        "record reference sequence {reference_sequence_id} exceeds header references"
                    ));
                }
            }
            None => {
                unplaced_unmapped_record_count += 1;
            }
        }

        start_position = end_position;
    }

    let reference_sequences = reference_sequences
        .into_iter()
        .map(ReferenceSequenceIndexBuilder::build)
        .collect();

    Ok(bam::bai::Index::builder()
        .set_reference_sequences(reference_sequences)
        .set_unplaced_unmapped_record_count(unplaced_unmapped_record_count)
        .build())
}

/// Returns the alignment context used for BAM indexing.
fn record_alignment_context(
    record: &bam::Record,
) -> Result<Option<(usize, Position, Position, bool)>> {
    Ok(
        match (
            record.reference_sequence_id().transpose()?,
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) {
            (Some(reference_sequence_id), Some(start), Some(end)) => Some((
                reference_sequence_id,
                start,
                end,
                !record.flags().is_unmapped(),
            )),
            _ => None,
        },
    )
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
        "{} --bam {} --annotation {} --sequence {} --outdir {}{}{}{} --threads {} --wiggle {} --clip-cutoff {} --level {}",
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
        if cli.extend_alignment {
            " --extend-alignment"
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
fn reverse_corrected_start(start0: u64, reference_extension: u64, delta: i64) -> Result<u64> {
    let shifted = start0
        .checked_sub(reference_extension)
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

/// Adds a chunk to a bin, merging directly adjacent overlaps in file order.
fn push_chunk(chunks: &mut Vec<Chunk>, chunk: Chunk) {
    if let Some(last_chunk) = chunks.last_mut() {
        if chunk.start() <= last_chunk.end() {
            *last_chunk = Chunk::new(last_chunk.start(), chunk.end());
            return;
        }
    }

    chunks.push(chunk);
}

/// Maps an alignment interval to its BAI bin.
fn reg2bin(start: Position, end: Position, min_shift: u8, depth: u8) -> usize {
    let beg = usize::from(start) - 1;
    let end = usize::from(end) - 1;

    let mut level = depth;
    let mut shift = min_shift;
    let mut offset = ((1 << (depth * 3)) - 1) / 7;

    while level > 0 {
        if beg >> shift == end >> shift {
            return offset + (beg >> shift);
        }

        level -= 1;
        shift += 3;
        offset -= 1 << (level * 3);
    }

    0
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::num::NonZeroUsize;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;

    use super::*;
    use noodles_sam::alignment::record_buf::{Cigar, Sequence};
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
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

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .set_reference_sequences(
                [(
                    "chr1".into(),
                    Map::<ReferenceSequence>::new(NonZeroUsize::try_from(512).unwrap()),
                )]
                .into_iter()
                .collect(),
            )
            .build()
    }

    fn test_engine_with_options(
        annotation_path: &Path,
        genome_path: &Path,
        extend_alignment: bool,
        keep_additional_corrections: bool,
    ) -> Engine {
        let annotation = Arc::new(AnnotationIndex::load(annotation_path).unwrap());
        let genome = GenomeSource::open(genome_path).unwrap();
        Engine::new(
            Config {
                threads: 1,
                split_bam: false,
                wiggle: 2,
                clip_cutoff: 5,
                keep_additional_corrections,
                extend_alignment,
            },
            annotation,
            genome,
            "iso-cigar".to_string(),
        )
        .unwrap()
    }

    fn test_engine(annotation_path: &Path, genome_path: &Path) -> Engine {
        test_engine_with_options(annotation_path, genome_path, false, false)
    }

    fn expect_corrected(outcome: CorrectionOutcome) -> (RecordBuf, Vec<RecordBuf>) {
        match outcome {
            CorrectionOutcome::Corrected {
                primary,
                additional,
            } => (primary, additional),
            CorrectionOutcome::Unchanged(_) => panic!("expected corrected record"),
        }
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
    fn bam_index_path_appends_bai_suffix() {
        assert_eq!(
            bam_index_path(Path::new("reads.extended.bam")),
            PathBuf::from("reads.extended.bam.bai")
        );
    }

    #[test]
    fn write_bam_index_creates_associated_bai() {
        let temp = tempdir().unwrap();
        let bam_path = temp.path().join("out.bam");
        let header = test_header();

        let mut writer = bam::io::Writer::new(File::create(&bam_path).unwrap());
        writer.write_header(&header).unwrap();
        writer
            .write_alignment_record(
                &header,
                &make_record(Flags::empty(), 101, &[Op::new(Kind::Match, 8)], b"ACGTACGT"),
            )
            .unwrap();
        writer
            .write_alignment_record(
                &header,
                &make_record(Flags::empty(), 21, &[Op::new(Kind::Match, 8)], b"TTTTGGGG"),
            )
            .unwrap();
        writer.try_finish().unwrap();

        write_bam_index(&bam_path).unwrap();

        let index_path = bam_index_path(&bam_path);
        assert!(index_path.exists());
        let index = bam::bai::fs::read(&index_path).unwrap();
        assert_eq!(index.reference_sequences().len(), 1);
    }

    #[test]
    fn remove_bam_and_index_removes_both_files() {
        let temp = tempdir().unwrap();
        let bam_path = temp.path().join("out.bam");
        let bai_path = bam_index_path(&bam_path);

        fs::write(&bam_path, b"bam").unwrap();
        fs::write(&bai_path, b"bai").unwrap();

        remove_bam_and_index(&bam_path).unwrap();

        assert!(!bam_path.exists());
        assert!(!bai_path.exists());
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
    fn candidate_scoring_prefers_more_exons_then_smaller_delta() {
        let a = CandidateCorrection {
            junction_id: 3,
            delta: 0,
            intron_len: 10,
            rescued_from_softclip: 6,
            rescued_exons: 2,
            tie_breaker: 3,
            new_start0: 10,
            genomic_ops: vec![Op::new(Kind::Match, 10)],
        };
        let b = CandidateCorrection {
            junction_id: 2,
            delta: 2,
            intron_len: 10,
            rescued_from_softclip: 6,
            rescued_exons: 1,
            tie_breaker: 2,
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
    fn correct_record_requires_extend_alignment_for_forward_multi_exon_rescue() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 220];
        set_bases(&mut genome, 100, b"ACGTACGT");
        set_bases(&mut genome, 118, b"TTGC");
        set_bases(&mut genome, 130, b"GGAAC");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            b"chr1\t100\t135\ttx1\t0\t+\t100\t135\t0,0,0\t3\t8,4,5\t0,18,30\n",
        )
        .unwrap();

        let record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 9)],
            b"ACGTACGTTTGCGGAAC",
        );

        let mut engine = test_engine(&annotation_path, &genome_path);
        let outcome = engine.correct_record(record.clone(), &[Some(0)]).unwrap();
        assert!(matches!(outcome, CorrectionOutcome::Unchanged(_)));

        let mut extend_engine =
            test_engine_with_options(&annotation_path, &genome_path, true, false);
        let (corrected, additional) =
            expect_corrected(extend_engine.correct_record(record, &[Some(0)]).unwrap());

        assert!(additional.is_empty());
        assert_eq!(corrected.alignment_start().unwrap().get(), 101);
        assert_eq!(
            cigar_signature(&corrected),
            vec![(b'M', 8), (b'N', 10), (b'=', 4), (b'N', 8), (b'=', 5)]
        );
        assert_eq!(
            corrected.data().get(&TAG_RESCUED),
            Some(&Value::from(9_u32))
        );
        assert_eq!(
            corrected.data().get(&TAG_INTRON),
            Some(&Value::from(18_u32))
        );
    }

    #[test]
    fn correct_record_keeps_one_exon_rescue_with_extend_alignment_enabled() {
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

        let record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 5)],
            b"ACGTACGTTTGCA",
        );

        let mut engine = test_engine(&annotation_path, &genome_path);
        let (baseline, _) =
            expect_corrected(engine.correct_record(record.clone(), &[Some(0)]).unwrap());

        let mut extend_engine =
            test_engine_with_options(&annotation_path, &genome_path, true, false);
        let (corrected, additional) =
            expect_corrected(extend_engine.correct_record(record, &[Some(0)]).unwrap());

        assert!(additional.is_empty());
        assert_eq!(
            corrected.alignment_start().unwrap().get(),
            baseline.alignment_start().unwrap().get()
        );
        assert_eq!(cigar_signature(&corrected), cigar_signature(&baseline));
        assert_eq!(
            corrected.data().get(&TAG_RESCUED),
            baseline.data().get(&TAG_RESCUED)
        );
    }

    #[test]
    fn correct_record_rescues_reverse_two_downstream_exons() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 240];
        set_bases(&mut genome, 100, b"GTTCC");
        set_bases(&mut genome, 113, b"GCAA");
        set_bases(&mut genome, 125, b"ACGTACGT");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            b"chr1\t100\t133\ttx1\t0\t-\t100\t133\t0,0,0\t3\t5,4,8\t0,13,25\n",
        )
        .unwrap();

        let mut engine = test_engine_with_options(&annotation_path, &genome_path, true, false);
        let record = make_record(
            Flags::REVERSE_COMPLEMENTED,
            126,
            &[Op::new(Kind::SoftClip, 9), Op::new(Kind::Match, 8)],
            b"GTTCCGCAAACGTACGT",
        );

        let (corrected, additional) =
            expect_corrected(engine.correct_record(record, &[Some(0)]).unwrap());

        assert!(additional.is_empty());
        assert_eq!(corrected.alignment_start().unwrap().get(), 101);
        assert_eq!(
            cigar_signature(&corrected),
            vec![(b'=', 5), (b'N', 8), (b'=', 4), (b'N', 8), (b'M', 8)]
        );
        assert_eq!(
            corrected.data().get(&TAG_RESCUED),
            Some(&Value::from(9_u32))
        );
    }

    #[test]
    fn correct_record_stops_after_partial_second_downstream_exon() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 220];
        set_bases(&mut genome, 100, b"ACGTACGT");
        set_bases(&mut genome, 118, b"TTGC");
        set_bases(&mut genome, 130, b"GGAAC");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            b"chr1\t100\t135\ttx1\t0\t+\t100\t135\t0,0,0\t3\t8,4,5\t0,18,30\n",
        )
        .unwrap();

        let mut engine = test_engine_with_options(&annotation_path, &genome_path, true, false);
        let mut record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 8)],
            b"ACGTACGTTTGCGGTT",
        );
        record
            .data_mut()
            .insert(Tag::EDIT_DISTANCE, Value::from(3_u8));
        record
            .data_mut()
            .insert(Tag::MISMATCHED_POSITIONS, Value::from("9A1"));

        let (corrected, additional) =
            expect_corrected(engine.correct_record(record, &[Some(0)]).unwrap());

        assert!(additional.is_empty());
        assert_eq!(
            cigar_signature(&corrected),
            vec![
                (b'M', 8),
                (b'N', 10),
                (b'=', 4),
                (b'N', 8),
                (b'=', 2),
                (b'S', 2)
            ]
        );
        assert_eq!(
            corrected.data().get(&TAG_RESCUED),
            Some(&Value::from(6_u32))
        );
        assert!(corrected.data().get(&Tag::EDIT_DISTANCE).is_none());
        assert!(corrected.data().get(&Tag::MISMATCHED_POSITIONS).is_none());
    }

    #[test]
    fn correct_record_rejects_multi_exon_when_second_exon_mismatches_immediately() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 220];
        set_bases(&mut genome, 100, b"ACGTACGT");
        set_bases(&mut genome, 118, b"TTGC");
        set_bases(&mut genome, 130, b"GGAAC");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            b"chr1\t100\t135\ttx1\t0\t+\t100\t135\t0,0,0\t3\t8,4,5\t0,18,30\n",
        )
        .unwrap();

        let mut engine = test_engine_with_options(&annotation_path, &genome_path, true, false);
        let record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 6)],
            b"ACGTACGTTTGCTT",
        );

        let outcome = engine.correct_record(record, &[Some(0)]).unwrap();
        assert!(matches!(outcome, CorrectionOutcome::Unchanged(_)));
    }

    #[test]
    fn correct_record_prefers_best_transcript_when_paths_diverge() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 260];
        set_bases(&mut genome, 100, b"ACGTACGT");
        set_bases(&mut genome, 118, b"TTGCA");
        set_bases(&mut genome, 131, b"GGAA");
        set_bases(&mut genome, 141, b"CCCC");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            concat!(
                "chr1\t100\t135\ttxA\t0\t+\t100\t135\t0,0,0\t3\t8,5,4\t0,18,31\n",
                "chr1\t100\t145\ttxB\t0\t+\t100\t145\t0,0,0\t3\t8,5,4\t0,18,41\n"
            ),
        )
        .unwrap();

        let mut engine = test_engine_with_options(&annotation_path, &genome_path, true, false);
        let record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 9)],
            b"ACGTACGTTTGCAGGAA",
        );

        let (corrected, additional) =
            expect_corrected(engine.correct_record(record, &[Some(0)]).unwrap());

        assert!(additional.is_empty());
        assert_eq!(
            cigar_signature(&corrected),
            vec![(b'M', 8), (b'N', 10), (b'=', 5), (b'N', 8), (b'=', 4)]
        );
        assert_eq!(
            corrected.data().get(&TAG_RESCUED),
            Some(&Value::from(9_u32))
        );
    }

    #[test]
    fn correct_record_suppresses_duplicate_transcript_outcomes() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 220];
        set_bases(&mut genome, 100, b"ACGTACGT");
        set_bases(&mut genome, 118, b"TTGC");
        set_bases(&mut genome, 130, b"GGAAC");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            concat!(
                "chr1\t100\t135\ttx1\t0\t+\t100\t135\t0,0,0\t3\t8,4,5\t0,18,30\n",
                "chr1\t100\t135\ttx2\t0\t+\t100\t135\t0,0,0\t3\t8,4,5\t0,18,30\n"
            ),
        )
        .unwrap();

        let mut engine = test_engine_with_options(&annotation_path, &genome_path, true, true);
        let record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 9)],
            b"ACGTACGTTTGCGGAAC",
        );

        let (corrected, additional) =
            expect_corrected(engine.correct_record(record, &[Some(0)]).unwrap());

        assert_eq!(
            cigar_signature(&corrected),
            vec![(b'M', 8), (b'N', 10), (b'=', 4), (b'N', 8), (b'=', 5)]
        );
        assert!(additional.is_empty());
    }

    #[test]
    fn correct_record_keeps_distinct_additional_corrections() {
        let temp = tempdir().unwrap();
        let genome_path = temp.path().join("genome.fa");
        let annotation_path = temp.path().join("annotation.bed");

        let mut genome = vec![b'A'; 260];
        set_bases(&mut genome, 100, b"ACGTACGT");
        set_bases(&mut genome, 118, b"TTGCA");
        set_bases(&mut genome, 131, b"GGAA");
        set_bases(&mut genome, 141, b"CCCC");
        write_fasta(&genome_path, "chr1", &genome);
        fs::write(
            &annotation_path,
            concat!(
                "chr1\t100\t135\ttxA\t0\t+\t100\t135\t0,0,0\t3\t8,5,4\t0,18,31\n",
                "chr1\t100\t145\ttxB\t0\t+\t100\t145\t0,0,0\t3\t8,5,4\t0,18,41\n"
            ),
        )
        .unwrap();

        let mut engine = test_engine_with_options(&annotation_path, &genome_path, true, true);
        let record = make_record(
            Flags::empty(),
            101,
            &[Op::new(Kind::Match, 8), Op::new(Kind::SoftClip, 9)],
            b"ACGTACGTTTGCAGGAA",
        );

        let (corrected, additional) =
            expect_corrected(engine.correct_record(record, &[Some(0)]).unwrap());

        assert_eq!(
            cigar_signature(&corrected),
            vec![(b'M', 8), (b'N', 10), (b'=', 5), (b'N', 8), (b'=', 4)]
        );
        assert_eq!(additional.len(), 1);
        assert!(additional[0].flags().is_secondary());
        assert_eq!(
            cigar_signature(&additional[0]),
            vec![(b'M', 8), (b'N', 10), (b'=', 5), (b'S', 4)]
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
