// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Reader / worker / writer orchestration for iso-adapter.
//!
//! A single reader thread pulls `RecordBuf`s out of the input BAM and pushes
//! them onto a bounded work queue. `n = --threads` workers consume the queue,
//! run detection, optionally trim, and push the outcome onto an output
//! queue. When `--remove-adapters` is set, a single writer consumes the
//! output queue and writes trimmed records to a new BAM; otherwise detection
//! runs in "scan" mode and the writer queue is drained in-place.

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
use noodles_sam::alignment::record_buf::RecordBuf;

use crate::cli::Cli;
use crate::detector::{
    extract_clips, AdapterDb, AdapterMatch, ClipEnd, SharedAdapterDb, SharedUnknownClips,
    UnknownClipStore,
};
use crate::logging::elapsed;
use crate::stats::ClipStats;
use crate::writer::{trim_record, TrimOutcome};

/// Runs the iso-adapter pipeline.
pub fn run(cli: Cli) -> Result<()> {
    info!(
        "[{}] start threads={} remove={} max_edit={} min_clip={} freq_threshold={}",
        elapsed(),
        cli.threads,
        cli.remove_adapters,
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

    let output_path = if cli.remove_adapters {
        Some(cli.resolved_output_bam())
    } else {
        None
    };

    if let Some(out) = output_path.as_ref() {
        info!("[{}] writing trimmed BAM {}", elapsed(), out.display());
        let writer_header = header.clone();
        write_trimmed_bam(
            out,
            &writer_header,
            result_rx,
            Arc::clone(&stats),
        )?;
    } else {
        drain_results(result_rx, Arc::clone(&stats))?;
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

/// Worker-side configuration passed by value to avoid sharing the whole CLI.
#[derive(Clone, Debug)]
struct WorkerConfig {
    remove_adapters: bool,
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
        let outcome = process_record(&mut item.record, &config, &db, &stats, &unknown);
        result_tx
            .send(ProcessedItem {
                index: item.index,
                record: item.record,
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

/// Runs detection (and optionally trimming) on a single record.
fn process_record(
    record: &mut RecordBuf,
    config: &WorkerConfig,
    db: &AdapterDb,
    stats: &ClipStats,
    unknown: &UnknownClipStore,
) -> RecordOutcome {
    stats.record_seen();

    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        stats.record_skipped();
        return RecordOutcome::Unchanged;
    }

    // We need the clip byte slices and their lengths before calling `trim_record`
    // (which borrows `record` mutably), so copy them to owned `Vec<u8>` first.
    let (five_prime_clip, three_prime_clip) = {
        let (five, three) = extract_clips(record);
        (five.map(|s| s.to_vec()), three.map(|s| s.to_vec()))
    };

    let mut had_adapter = false;
    let mut had_unknown = false;

    // Process 5' end.
    if let Some(clip) = five_prime_clip.as_ref() {
        stats.observe_clip(ClipEnd::FivePrime, clip.len());
        match db.match_adapter(clip) {
            Some(m) => {
                trace!(
                    "[{}] 5p hit label={} edit={} range={:?}",
                    elapsed(),
                    m.label,
                    m.edit_distance,
                    m.clip_range
                );
                had_adapter = true;
                stats.observe_match(&m);
                if config.remove_adapters {
                    let clip_len = clip.len();
                    apply_trim(record, ClipEnd::FivePrime, clip_len, &m, stats);
                }
            }
            None => {
                if clip.len() >= config.min_clip_len {
                    unknown.record(clip);
                    had_unknown = true;
                }
            }
        }
    }

    // Process 3' end. NB: after 5' trim the sequence may have changed, so
    // recompute the three-prime clip length if we trimmed at the 5' end.
    let recomputed_three = if config.remove_adapters && five_prime_clip.is_some() {
        extract_clips(record).1.map(|s| s.to_vec())
    } else {
        three_prime_clip
    };

    if let Some(clip) = recomputed_three.as_ref() {
        stats.observe_clip(ClipEnd::ThreePrime, clip.len());
        match db.match_adapter(clip) {
            Some(m) => {
                trace!(
                    "[{}] 3p hit label={} edit={} range={:?}",
                    elapsed(),
                    m.label,
                    m.edit_distance,
                    m.clip_range
                );
                had_adapter = true;
                stats.observe_match(&m);
                if config.remove_adapters {
                    let clip_len = clip.len();
                    apply_trim(record, ClipEnd::ThreePrime, clip_len, &m, stats);
                }
            }
            None => {
                if clip.len() >= config.min_clip_len {
                    unknown.record(clip);
                    had_unknown = true;
                }
            }
        }
    }

    if had_adapter {
        stats.record_with_adapter();
    }
    if had_unknown && !had_adapter {
        stats.record_with_unknown_clip();
    }

    RecordOutcome::Processed { had_adapter }
}

/// Applies a trim and updates stats.
fn apply_trim(
    record: &mut RecordBuf,
    end: ClipEnd,
    clip_len: usize,
    m: &AdapterMatch,
    stats: &ClipStats,
) {
    match trim_record(record, end, clip_len, m) {
        TrimOutcome::Trimmed { bases_removed } => {
            trace!(
                "[{}] trimmed {end:?} bases_removed={bases_removed}",
                elapsed()
            );
            stats.record_trimmed();
        }
        TrimOutcome::Unchanged => {
            trace!("[{}] trim skipped (degenerate match) end={end:?}", elapsed());
        }
    }
}

/// Drains the result queue without writing a BAM (scan-only mode).
fn drain_results(result_rx: Receiver<ProcessedItem>, _stats: Arc<ClipStats>) -> Result<()> {
    while result_rx.recv().is_ok() {}
    Ok(())
}

/// Writes trimmed records to a new BAM, preserving input order.
fn write_trimmed_bam(
    output: &Path,
    header: &sam::Header,
    result_rx: Receiver<ProcessedItem>,
    _stats: Arc<ClipStats>,
) -> Result<()> {
    use std::collections::BTreeMap;

    if let Some(parent) = output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("failed to create output directory {}", parent.display()))?;
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
        let ProcessedItem { index, record, .. } = item;
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

/// Item flowing from reader to workers.
struct WorkItem {
    index: u64,
    record: RecordBuf,
}

/// Item flowing from workers to writer.
struct ProcessedItem {
    index: u64,
    record: RecordBuf,
    #[allow(dead_code)]
    outcome: RecordOutcome,
}

/// Per-record processing outcome (currently only used for logging / future hooks).
#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
enum RecordOutcome {
    Unchanged,
    Processed { had_adapter: bool },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detector::{AdapterDb, UnknownClipStore};
    use noodles_sam::alignment::record::cigar::{op::Kind, Op};
    use noodles_sam::alignment::record_buf::{Cigar as CigarBuf, Sequence};

    #[test]
    fn process_record_unmapped_is_skipped() {
        let db = AdapterDb::new(10, 3).unwrap();
        let stats = ClipStats::default();
        let unknown = UnknownClipStore::default();

        // RecordBuf::default() already has Flags::UNMAPPED set.
        let mut rec = RecordBuf::default();
        *rec.cigar_mut() = CigarBuf::default();

        let outcome = process_record(
            &mut rec,
            &WorkerConfig {
                remove_adapters: false,
                min_clip_len: 10,
            },
            &db,
            &stats,
            &unknown,
        );
        assert!(matches!(outcome, RecordOutcome::Unchanged));
        assert_eq!(stats.total(), 1);
    }

    #[test]
    fn process_record_finds_adapter_in_5p_clip() {
        let db = AdapterDb::new(10, 3).unwrap();
        let stats = ClipStats::default();
        let unknown = UnknownClipStore::default();

        let primer = b"AAGCAGTGGTATCAACGCAGAGTAC";
        let mut seq = primer.to_vec();
        seq.extend_from_slice(b"ACGTACGT");

        let mut rec = RecordBuf::default();
        // Clear UNMAPPED so the record is treated as mapped.
        *rec.flags_mut() = sam::alignment::record::Flags::empty();
        *rec.cigar_mut() = CigarBuf::from(vec![
            Op::new(Kind::SoftClip, primer.len()),
            Op::new(Kind::Match, 8),
        ]);
        *rec.sequence_mut() = Sequence::from(seq);

        let _ = process_record(
            &mut rec,
            &WorkerConfig {
                remove_adapters: false,
                min_clip_len: 10,
            },
            &db,
            &stats,
            &unknown,
        );
        assert!(stats.matched() >= 1);
    }
}
