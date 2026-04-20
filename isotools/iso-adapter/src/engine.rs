// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Reader / worker / writer orchestration for iso-adapter.
//!
//! A single reader thread pulls `RecordBuf`s out of the input BAM and pushes
//! them onto a bounded work queue. `n = --threads` workers consume the queue,
//! run detection, optionally trim, and push the outcome onto an output
//! queue. A single writer consumes the output queue and writes records to
//! a new BAM when any trimming flag is set; otherwise the queue is drained
//! in place.
//!
//! Per-end processing runs in two passes so that a clip containing, say,
//! a cDNA primer *outward* of a polyA tail can have both handled:
//!
//! 1. exact + fuzzy match → if `Adapter` and `-R` is set, trim adapter
//!    + everything outward of it; stamp `ad:i:N`. Re-extract the clip.
//! 2. run the match again on the (now shorter) clip → if `Homopolymer` and
//!    `-P` is set, apply the polyA policy and stamp `pa:i:N` (3') or
//!    `pt:i:N` (5').

use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;
use std::thread;

use anyhow::{anyhow, Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use log::{debug, info, trace, warn};
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::RecordBuf;

use crate::cli::Cli;
use crate::detector::{
    extract_clips, AdapterDb, AdapterKind, AdapterMatch, ClipEnd, SharedAdapterDb,
    SharedUnknownClips, UnknownClipStore,
};
use crate::logging::elapsed;
use crate::stats::ClipStats;
use crate::writer::{trim_clip, TrimOutcome};

/// `ad:i:<N>` — number of bases trimmed by `--remove-adapters` at that end.
const TAG_ADAPTER: Tag = Tag::new(b'a', b'd');
/// `pa:i:<N>` — number of bases trimmed outward of a 3' polyA tail.
const TAG_POLYA: Tag = Tag::new(b'p', b'a');
/// `pt:i:<N>` — number of bases trimmed at a 5' polyT run (run + outward).
const TAG_POLYT: Tag = Tag::new(b'p', b't');

/// Runs the iso-adapter pipeline.
pub fn run(cli: Cli) -> Result<()> {
    let will_trim = cli.remove_adapters || cli.trim_polya;

    info!(
        "[{}] start threads={} remove_adapters={} trim_polya={} max_edit={} min_clip={} freq_threshold={}",
        elapsed(),
        cli.threads,
        cli.remove_adapters,
        cli.trim_polya,
        cli.max_edit_dist,
        cli.min_clip_len,
        cli.freq_threshold
    );

    let input_path = cli.input.clone();
    let bam_file = File::open(&input_path)
        .with_context(|| format!("failed to open BAM {}", input_path.display()))?;
    let mut reader = bam::io::Reader::new(bam_file);
    let header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header from {}", input_path.display()))?;

    let db: SharedAdapterDb = Arc::new(
        AdapterDb::new(cli.min_clip_len, cli.max_edit_dist)
            .context("failed to build adapter matcher")?,
    );
    info!(
        "[{}] adapter matcher ready patterns={} (includes reverse complements)",
        elapsed(),
        db.pattern_count()
    );

    let stats = Arc::new(ClipStats::default());
    let unknown: SharedUnknownClips = Arc::new(UnknownClipStore::default());

    let worker_count = cli.threads.max(1);
    let queue_capacity = (worker_count * 8).max(16);
    let (work_tx, work_rx) = bounded::<WorkItem>(queue_capacity);
    let (result_tx, result_rx) = bounded::<ProcessedItem>(queue_capacity);

    let reader_header = header.clone();
    let reader_handle = thread::Builder::new()
        .name("iso-adapter-reader".to_string())
        .spawn(move || read_records(reader, reader_header, work_tx))
        .context("failed to spawn reader thread")?;

    let config = WorkerConfig {
        remove_adapters: cli.remove_adapters,
        trim_polya: cli.trim_polya,
        min_clip_len: cli.min_clip_len,
    };

    let mut worker_handles = Vec::with_capacity(worker_count);
    for worker_id in 0..worker_count {
        let worker_rx = work_rx.clone();
        let worker_tx = result_tx.clone();
        let worker_db = Arc::clone(&db);
        let worker_stats = Arc::clone(&stats);
        let worker_unknown = Arc::clone(&unknown);
        let worker_config = config.clone();

        let handle = thread::Builder::new()
            .name(format!("iso-adapter-worker-{worker_id}"))
            .spawn(move || {
                worker_main(
                    worker_id,
                    worker_config,
                    worker_db,
                    worker_stats,
                    worker_unknown,
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

    let output_path = if will_trim {
        Some(cli.resolved_output_bam())
    } else {
        None
    };

    if let Some(out) = output_path.as_ref() {
        info!("[{}] writing trimmed BAM {}", elapsed(), out.display());
        write_trimmed_bam(out, &header, result_rx)?;
    } else {
        drain_results(result_rx)?;
    }

    join_thread(reader_handle, "reader")??;
    for (worker_id, handle) in worker_handles.into_iter().enumerate() {
        join_thread(handle, &format!("worker-{worker_id}"))??;
    }

    let novel = unknown.novel_candidates(cli.freq_threshold);
    if !novel.is_empty() {
        warn!(
            "[{}] {} novel-clip candidate(s) exceeded freq_threshold={}",
            elapsed(),
            novel.len(),
            cli.freq_threshold
        );
    }

    stats.report(&novel);
    info!("[{}] done!", elapsed());

    Ok(())
}

/// Per-worker configuration.
#[derive(Clone, Debug)]
struct WorkerConfig {
    remove_adapters: bool,
    trim_polya: bool,
    min_clip_len: usize,
}

/// Reads BAM records and ships them to the worker queue.
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
        let bytes = reader
            .read_record_buf(&header, &mut record)
            .with_context(|| format!("failed to read record #{index}"))?;
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

        if index % 1_000_000 == 0 {
            info!("[{}] read {index} records", elapsed());
        }
    }

    info!("[{}] reader finished records={index}", elapsed());
    Ok(())
}

/// Worker thread main function.
fn worker_main(
    worker_id: usize,
    config: WorkerConfig,
    db: SharedAdapterDb,
    stats: Arc<ClipStats>,
    unknown: SharedUnknownClips,
    work_rx: Receiver<WorkItem>,
    result_tx: Sender<ProcessedItem>,
) -> Result<()> {
    debug!("[{}] worker-{worker_id} start", elapsed());
    let mut processed = 0_u64;

    while let Ok(mut item) = work_rx.recv() {
        process_record(&mut item.record, &config, &db, &stats, &unknown);
        result_tx
            .send(ProcessedItem {
                index: item.index,
                record: item.record,
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

/// Runs detection (and optionally trimming) on a single record.
fn process_record(
    record: &mut RecordBuf,
    config: &WorkerConfig,
    db: &AdapterDb,
    stats: &ClipStats,
    unknown: &UnknownClipStore,
) {
    stats.record_seen();

    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        stats.record_skipped();
        return;
    }

    let mut had_adapter = false;
    let mut had_unknown = false;

    for end in [ClipEnd::FivePrime, ClipEnd::ThreePrime] {
        let result = process_end(record, end, config, db, stats, unknown);
        if result.had_adapter {
            had_adapter = true;
        }
        if result.had_unknown {
            had_unknown = true;
        }
    }

    if had_adapter {
        stats.record_with_adapter();
    }
    if had_unknown && !had_adapter {
        stats.record_with_unknown_clip();
    }
}

#[derive(Default)]
struct EndOutcome {
    had_adapter: bool,
    had_unknown: bool,
}

/// Handles detection + up-to-two-pass trimming for a single end.
fn process_end(
    record: &mut RecordBuf,
    end: ClipEnd,
    config: &WorkerConfig,
    db: &AdapterDb,
    stats: &ClipStats,
    unknown: &UnknownClipStore,
) -> EndOutcome {
    let mut out = EndOutcome::default();
    // First pass: adapter or homopolymer, whichever the matcher finds.
    let Some(first_clip) = clip_for_end(record, end) else {
        return out;
    };
    stats.observe_clip(end, first_clip.len());

    let first_match = db.match_adapter(&first_clip);
    match first_match {
        Some(m) => {
            trace!(
                "[{}] match end={:?} label={} edit={} range={:?} kind={:?}",
                elapsed(),
                end,
                m.label,
                m.edit_distance,
                m.clip_range,
                m.kind
            );
            out.had_adapter = true;
            stats.observe_match(&m);
            apply_policy(record, end, &first_clip, &m, config, stats);
        }
        None => {
            if first_clip.len() >= config.min_clip_len {
                unknown.record(&first_clip);
                out.had_unknown = true;
            }
            return out;
        }
    }

    // Second pass: after trimming an outward adapter, the clip may now end
    // in a polyA tail (or start with a polyT) that the adapter trim exposed.
    if !(config.remove_adapters && config.trim_polya) {
        return out;
    }
    let Some(second_clip) = clip_for_end(record, end) else {
        return out;
    };
    if second_clip.len() < config.min_clip_len {
        return out;
    }
    if let Some(m) = db.match_adapter(&second_clip) {
        if m.kind == AdapterKind::Homopolymer {
            trace!(
                "[{}] second-pass polyA match end={:?} label={} range={:?}",
                elapsed(),
                end,
                m.label,
                m.clip_range
            );
            stats.observe_match(&m);
            apply_policy(record, end, &second_clip, &m, config, stats);
        }
    }
    out
}

/// Dispatches the trim policy that applies to `(kind, end)` and stamps the
/// corresponding BAM tag on success.
fn apply_policy(
    record: &mut RecordBuf,
    end: ClipEnd,
    clip: &[u8],
    m: &AdapterMatch,
    config: &WorkerConfig,
    stats: &ClipStats,
) {
    let clip_len = clip.len();
    // For homopolymer hits the leftmost-longest automaton will anchor at the
    // first T/A in the clip, which may be one base earlier than the actual
    // polyA/polyT run if the preceding noise happens to end with the same
    // nucleotide. Widen the range to cover the full contiguous run so the
    // trim excises every A (or T) on that end.
    let effective_range = match m.kind {
        AdapterKind::Homopolymer => extend_homopolymer(clip, &m.clip_range),
        AdapterKind::Adapter => m.clip_range.clone(),
    };

    let Some(plan) = trim_plan(end, clip_len, m.kind, &effective_range, config) else {
        return;
    };

    match trim_clip(record, end, clip_len, plan.cut_start, plan.cut_end) {
        TrimOutcome::Trimmed { bases_removed } => {
            trace!(
                "[{}] trimmed end={:?} tag={:?} bases_removed={bases_removed}",
                elapsed(),
                end,
                plan.tag
            );
            stats.record_trimmed();
            stamp_tag(record, plan.tag, bases_removed);
        }
        TrimOutcome::Unchanged => {
            trace!(
                "[{}] trim skipped (degenerate) end={:?} range={:?}",
                elapsed(),
                end,
                m.clip_range
            );
        }
    }
}

/// Translates `(kind, end, flags)` into a concrete cut range within the clip.
///
/// Returns `None` when the user flag for this kind is not set, when the
/// match range is degenerate, or when there is nothing outward to trim.
fn trim_plan(
    end: ClipEnd,
    clip_len: usize,
    kind: AdapterKind,
    range: &std::ops::Range<usize>,
    config: &WorkerConfig,
) -> Option<TrimPlan> {
    if range.start >= range.end || range.end > clip_len {
        return None;
    }

    match (kind, end) {
        // Adapter at the 5' end: strip the adapter + any overhang outward
        // of it. Keep bases inward of the adapter (range.end .. clip_len).
        (AdapterKind::Adapter, ClipEnd::FivePrime) if config.remove_adapters => Some(TrimPlan {
            cut_start: 0,
            cut_end: range.end,
            tag: TAG_ADAPTER,
        }),
        // Adapter at the 3' end: strip the adapter + everything outward.
        (AdapterKind::Adapter, ClipEnd::ThreePrime) if config.remove_adapters => Some(TrimPlan {
            cut_start: range.start,
            cut_end: clip_len,
            tag: TAG_ADAPTER,
        }),
        // 3' homopolymer = polyA tail: keep the run, strip anything that
        // came after it (outward). Nothing to do if the run already reaches
        // the very end of the clip.
        (AdapterKind::Homopolymer, ClipEnd::ThreePrime) if config.trim_polya => {
            if range.end >= clip_len {
                None
            } else {
                Some(TrimPlan {
                    cut_start: range.end,
                    cut_end: clip_len,
                    tag: TAG_POLYA,
                })
            }
        }
        // 5' homopolymer = polyT with no biological meaning: hard-trim the
        // run itself *and* everything outward of it.
        (AdapterKind::Homopolymer, ClipEnd::FivePrime) if config.trim_polya => Some(TrimPlan {
            cut_start: 0,
            cut_end: range.end,
            tag: TAG_POLYT,
        }),
        _ => None,
    }
}

/// Widens `range` to cover the full contiguous run of the same base on both
/// sides inside `clip`. This compensates for leftmost-longest Aho–Corasick
/// picking a slightly shifted anchor when the surrounding noise includes
/// bases that match the homopolymer character.
fn extend_homopolymer(clip: &[u8], range: &std::ops::Range<usize>) -> std::ops::Range<usize> {
    if range.start >= clip.len() || range.start >= range.end {
        return range.clone();
    }
    let ch = clip[range.start];
    let mut start = range.start;
    while start > 0 && clip[start - 1] == ch {
        start -= 1;
    }
    let mut end = range.end;
    while end < clip.len() && clip[end] == ch {
        end += 1;
    }
    start..end
}

struct TrimPlan {
    cut_start: usize,
    cut_end: usize,
    tag: Tag,
}

/// Stamps a `*:i:N` BAM tag on a record. Overwrites any existing value at
/// that tag so that a second-pass trim adds cleanly on top of the first.
fn stamp_tag(record: &mut RecordBuf, tag: Tag, value: usize) {
    let data = record.data_mut();
    let existing = data.get(&tag).and_then(|v| match v {
        Value::UInt32(n) => Some(*n as usize),
        Value::Int32(n) => Some((*n).max(0) as usize),
        _ => None,
    });
    let total = existing.unwrap_or(0) + value;
    data.insert(tag, Value::from(total as u32));
}

/// Returns an owned copy of the specified end's soft-clip, if any.
fn clip_for_end(record: &RecordBuf, end: ClipEnd) -> Option<Vec<u8>> {
    let (five, three) = extract_clips(record);
    match end {
        ClipEnd::FivePrime => five.map(|s| s.to_vec()),
        ClipEnd::ThreePrime => three.map(|s| s.to_vec()),
    }
}

/// Drains the result queue without writing a BAM (scan-only mode).
fn drain_results(result_rx: Receiver<ProcessedItem>) -> Result<()> {
    while result_rx.recv().is_ok() {}
    Ok(())
}

/// Writes trimmed records to a new BAM, preserving input order.
fn write_trimmed_bam(
    output: &Path,
    header: &sam::Header,
    result_rx: Receiver<ProcessedItem>,
) -> Result<()> {
    use std::collections::BTreeMap;

    if let Some(parent) = output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).with_context(|| {
                format!("failed to create output directory {}", parent.display())
            })?;
        }
    }

    let file = File::create(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    let mut writer = bam::io::Writer::new(file);
    writer
        .write_header(header)
        .with_context(|| format!("failed to write header to {}", output.display()))?;

    let mut pending: BTreeMap<u64, RecordBuf> = BTreeMap::new();
    let mut next_index = 0_u64;

    for item in result_rx.iter() {
        let ProcessedItem { index, record } = item;
        pending.insert(index, record);
        while let Some(record) = pending.remove(&next_index) {
            writer
                .write_alignment_record(header, &record)
                .with_context(|| format!("failed to write record #{next_index}"))?;
            next_index += 1;
        }
    }

    if !pending.is_empty() {
        return Err(anyhow!(
            "writer finished with {} out-of-order records buffered",
            pending.len()
        ));
    }

    writer.try_finish()?;
    Ok(())
}

/// Joins a worker / reader thread, propagating panics.
fn join_thread<T>(handle: thread::JoinHandle<Result<T>>, label: &str) -> Result<Result<T>> {
    handle
        .join()
        .map_err(|_| anyhow!("{label} thread panicked"))
}

struct WorkItem {
    index: u64,
    record: RecordBuf,
}

struct ProcessedItem {
    index: u64,
    record: RecordBuf,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detector::{AdapterDb, UnknownClipStore};
    use noodles_sam::alignment::record::cigar::{op::Kind, Op};
    use noodles_sam::alignment::record_buf::{Cigar as CigarBuf, Sequence};

    fn make_record(cigar: Vec<Op>, seq: &[u8]) -> RecordBuf {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = sam::alignment::record::Flags::empty();
        *rec.cigar_mut() = CigarBuf::from(cigar);
        *rec.sequence_mut() = Sequence::from(seq.to_vec());
        rec
    }

    #[test]
    fn process_record_unmapped_is_skipped() {
        let db = AdapterDb::new(10, 3).unwrap();
        let stats = ClipStats::default();
        let unknown = UnknownClipStore::default();

        // RecordBuf::default() starts out with Flags::UNMAPPED.
        let mut rec = RecordBuf::default();
        *rec.cigar_mut() = CigarBuf::default();

        process_record(
            &mut rec,
            &WorkerConfig {
                remove_adapters: false,
                trim_polya: false,
                min_clip_len: 10,
            },
            &db,
            &stats,
            &unknown,
        );
        assert_eq!(stats.total(), 1);
    }

    #[test]
    fn adapter_5p_trim_removes_outward_and_adapter_keeps_inward() {
        let db = AdapterDb::new(10, 3).unwrap();
        let stats = ClipStats::default();
        let unknown = UnknownClipStore::default();

        // Layout (5'→3'): [noise 3 nt][primer 25 nt][inner overhang 2 nt] | [aligned 8 nt]
        //                                  5' soft-clip (30 nt)
        // Adapter-outward trim cuts [0 .. end_of_adapter) — that's noise + primer.
        // The inner overhang stays as a shortened soft-clip.
        let noise = b"GGG";
        let primer = b"AAGCAGTGGTATCAACGCAGAGTAC";
        let inner = b"AC";
        let aligned = b"ACGTACGT";
        let mut seq = Vec::new();
        seq.extend_from_slice(noise);
        seq.extend_from_slice(primer);
        seq.extend_from_slice(inner);
        seq.extend_from_slice(aligned);

        let clip_len = noise.len() + primer.len() + inner.len();
        let mut rec = make_record(
            vec![
                Op::new(Kind::SoftClip, clip_len),
                Op::new(Kind::Match, aligned.len()),
            ],
            &seq,
        );

        process_record(
            &mut rec,
            &WorkerConfig {
                remove_adapters: true,
                trim_polya: false,
                min_clip_len: 10,
            },
            &db,
            &stats,
            &unknown,
        );

        let mut expected = Vec::new();
        expected.extend_from_slice(inner);
        expected.extend_from_slice(aligned);
        assert_eq!(rec.sequence().as_ref(), expected.as_slice());

        let bases_trimmed = noise.len() + primer.len();
        let tag_value = rec.data().get(&TAG_ADAPTER).cloned();
        assert_eq!(tag_value, Some(Value::from(bases_trimmed as u32)));
    }

    #[test]
    fn polya_3p_trim_keeps_run_and_tags() {
        let db = AdapterDb::new(10, 3).unwrap();
        let stats = ClipStats::default();
        let unknown = UnknownClipStore::default();

        // 3' clip: 16 A's (polyA) + 6 arbitrary outward bases.
        let aligned = b"ACGTACGT";
        let polya = b"AAAAAAAAAAAAAAAA";
        let outward = b"GGTTCC";
        let mut seq = Vec::new();
        seq.extend_from_slice(aligned);
        seq.extend_from_slice(polya);
        seq.extend_from_slice(outward);

        let clip_len = polya.len() + outward.len();
        let mut rec = make_record(
            vec![
                Op::new(Kind::Match, aligned.len()),
                Op::new(Kind::SoftClip, clip_len),
            ],
            &seq,
        );

        process_record(
            &mut rec,
            &WorkerConfig {
                remove_adapters: false,
                trim_polya: true,
                min_clip_len: 10,
            },
            &db,
            &stats,
            &unknown,
        );

        // Aligned + polyA should be preserved; outward bases removed.
        let mut expected = Vec::new();
        expected.extend_from_slice(aligned);
        expected.extend_from_slice(polya);
        assert_eq!(rec.sequence().as_ref(), expected.as_slice());

        let pa_tag = rec.data().get(&TAG_POLYA).cloned();
        assert_eq!(pa_tag, Some(Value::from(outward.len() as u32)));
    }

    #[test]
    fn polyt_5p_trim_removes_run_and_outward() {
        let db = AdapterDb::new(10, 3).unwrap();
        let stats = ClipStats::default();
        let unknown = UnknownClipStore::default();

        // 5' clip: 4 noise bases + 16 T's (polyT). Both should disappear.
        let noise = b"GGAT";
        let polyt = b"TTTTTTTTTTTTTTTT";
        let aligned = b"ACGTACGT";
        let mut seq = Vec::new();
        seq.extend_from_slice(noise);
        seq.extend_from_slice(polyt);
        seq.extend_from_slice(aligned);

        let clip_len = noise.len() + polyt.len();
        let mut rec = make_record(
            vec![
                Op::new(Kind::SoftClip, clip_len),
                Op::new(Kind::Match, aligned.len()),
            ],
            &seq,
        );

        process_record(
            &mut rec,
            &WorkerConfig {
                remove_adapters: false,
                trim_polya: true,
                min_clip_len: 10,
            },
            &db,
            &stats,
            &unknown,
        );

        assert_eq!(rec.sequence().as_ref(), aligned);
        let pt_tag = rec.data().get(&TAG_POLYT).cloned();
        assert_eq!(pt_tag, Some(Value::from(clip_len as u32)));
    }

    #[test]
    fn remove_adapters_alone_does_not_touch_polya() {
        let db = AdapterDb::new(10, 3).unwrap();
        let stats = ClipStats::default();
        let unknown = UnknownClipStore::default();

        let aligned = b"ACGTACGT";
        let polya = b"AAAAAAAAAAAAAAAA";
        let mut seq = Vec::new();
        seq.extend_from_slice(aligned);
        seq.extend_from_slice(polya);

        let mut rec = make_record(
            vec![
                Op::new(Kind::Match, aligned.len()),
                Op::new(Kind::SoftClip, polya.len()),
            ],
            &seq,
        );

        process_record(
            &mut rec,
            &WorkerConfig {
                remove_adapters: true,
                trim_polya: false,
                min_clip_len: 10,
            },
            &db,
            &stats,
            &unknown,
        );

        // polyA untouched, no adapter tag stamped.
        assert_eq!(rec.sequence().as_ref().len(), aligned.len() + polya.len());
        assert!(rec.data().get(&TAG_ADAPTER).is_none());
        assert!(rec.data().get(&TAG_POLYA).is_none());
    }
}
