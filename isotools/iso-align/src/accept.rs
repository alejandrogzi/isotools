// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Accept mode: decide, per realigned read, whether the realigned record
//! should replace the read's pass-1 records.
//!
//! A read is accepted when both hold:
//!
//! 1. the realigned primary scores higher than the pass-1 primary
//!    (`AS` realigned > `AS` pass-1). `AS` is read from the pass-1 BAM because
//!    iso-cigar drops `AS`/`NM`/`SA` on the records it edits;
//! 2. the realigned primary is high quality: iso-segment identity >=
//!    `min_identity` and both clips <= `max_clip`.
//!
//! Accepted reads are written as names; ARK then drops their pass-1 records
//! (primary and supplementary) and keeps only accepted reads in the realigned
//! BAM.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, Result};
use log::info;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::extract::{write_names, CandidateSet};
use crate::logging::elapsed;
use crate::types::AcceptConfig;

/// Primary-record summary used by the acceptance test.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrimaryScore {
    /// `AS` tag, if present.
    pub alignment_score: Option<i64>,
    /// iso-segment identity: matches / (M/=/X + I + D), in percent, rounded
    /// to one decimal; matches = M/=/X - max(0, NM - (I + D)). `None` when
    /// the record has no `NM` tag (identity unknown: never accepted).
    pub identity: Option<f32>,
    /// Leading S+H (reference orientation).
    pub left_clip: u32,
    /// Trailing S+H (reference orientation).
    pub right_clip: u32,
    pub is_reverse: bool,
}

/// Outcome of the acceptance test for one realigned read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Decision {
    Accepted,
    /// `AS` realigned <= `AS` pass-1.
    LowerScore,
    /// Realigned record fails the identity / clip test, or has no `NM`.
    LowQuality,
    /// No pass-1 primary for this read, or `AS` missing on either side.
    MissingScore,
}

impl Decision {
    pub const fn label(self) -> &'static str {
        match self {
            Decision::Accepted => "accepted",
            Decision::LowerScore => "rejected_score",
            Decision::LowQuality => "rejected_quality",
            Decision::MissingScore => "rejected_missing_score",
        }
    }
}

/// The acceptance test (pure; see the module docs).
pub fn decide(
    pass1: Option<&PrimaryScore>,
    realigned: &PrimaryScore,
    cfg: &AcceptConfig,
) -> Decision {
    let (Some(old), Some(new)) = (
        pass1.and_then(|p| p.alignment_score),
        realigned.alignment_score,
    ) else {
        return Decision::MissingScore;
    };
    if new <= old {
        return Decision::LowerScore;
    }
    if realigned.identity.map_or(true, |id| id < cfg.min_identity)
        || realigned.left_clip > cfg.max_clip
        || realigned.right_clip > cfg.max_clip
    {
        return Decision::LowQuality;
    }
    Decision::Accepted
}

/// Computes clips and the iso-segment identity from a CIGAR op stream and NM.
pub fn identity_and_clips<I>(ops: I, nm: u64) -> std::io::Result<(f32, u32, u32)>
where
    I: IntoIterator<Item = std::io::Result<noodles_sam::alignment::record::cigar::Op>>,
{
    let mut aligned: u64 = 0;
    let mut indels: u64 = 0;
    let mut left_clip: u32 = 0;
    let mut right_clip: u32 = 0;
    let mut seen_body = false;
    for op in ops {
        let op = op?;
        let len = op.len() as u64;
        match op.kind() {
            Kind::SoftClip | Kind::HardClip => {
                if seen_body {
                    right_clip = right_clip.saturating_add(len as u32);
                } else {
                    left_clip = left_clip.saturating_add(len as u32);
                }
            }
            kind => {
                seen_body = true;
                right_clip = 0;
                match kind {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => aligned += len,
                    Kind::Insertion | Kind::Deletion => indels += len,
                    _ => {}
                }
            }
        }
    }
    let alignment_length = aligned + indels;
    let mismatches = nm.saturating_sub(indels);
    let matches = aligned.saturating_sub(mismatches);
    let identity = if alignment_length == 0 {
        0.0
    } else {
        (matches as f32 / alignment_length as f32 * 1000.0).round() / 10.0
    };
    Ok((identity, left_clip, right_clip))
}

/// Scans the primary records of `path`. With `only`, records of other reads
/// are skipped before their CIGAR is read.
pub fn scan_primaries(
    path: &Path,
    bgzf_workers: NonZeroUsize,
    only: Option<&FxHashMap<Vec<u8>, PrimaryScore>>,
) -> Result<FxHashMap<Vec<u8>, PrimaryScore>> {
    let file =
        File::open(path).with_context(|| format!("failed to open BAM {}", path.display()))?;
    let mut reader = bam::io::Reader::from(bgzf::MultithreadedReader::with_worker_count(
        bgzf_workers,
        file,
    ));
    reader
        .read_header()
        .with_context(|| format!("failed to read BAM header from {}", path.display()))?;

    let mut out = FxHashMap::default();
    let mut record = bam::Record::default();
    let mut total: u64 = 0;
    loop {
        let bytes = reader
            .read_record(&mut record)
            .with_context(|| format!("failed to read BAM record #{total} of {}", path.display()))?;
        if bytes == 0 {
            break;
        }
        total += 1;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let Some(name) = record.name() else {
            continue;
        };
        let name: &[u8] = name.as_ref();
        if only.is_some_and(|m| !m.contains_key(name)) {
            continue;
        }

        let data = record.data();
        let int_tag = |tag: Tag| -> Result<Option<i64>> {
            match data.get(&tag) {
                Some(v) => Ok(v
                    .with_context(|| format!("invalid {tag:?} tag at record #{total}"))?
                    .as_int()),
                None => Ok(None),
            }
        };
        let alignment_score = int_tag(Tag::ALIGNMENT_SCORE)?;
        let nm = int_tag(Tag::EDIT_DISTANCE)?;
        let (identity, left_clip, right_clip) =
            identity_and_clips(record.cigar().iter(), nm.unwrap_or(0).max(0) as u64)
                .with_context(|| format!("invalid CIGAR at record #{total}"))?;
        let identity = nm.map(|_| identity);

        out.insert(
            name.to_vec(),
            PrimaryScore {
                alignment_score,
                identity,
                left_clip,
                right_clip,
                is_reverse: flags.is_reverse_complemented(),
            },
        );
    }
    info!(
        "[{}] scanned {} records, {} primaries kept from {}",
        elapsed(),
        total,
        out.len(),
        path.display()
    );
    Ok(out)
}

/// Runs accept mode: writes accepted read names to `output`, and optionally a
/// per-read report.
pub fn run_accept(
    pass1_bam: &Path,
    realigned_bam: &Path,
    output: &Path,
    report: Option<&Path>,
    cfg: &AcceptConfig,
    bgzf_workers: NonZeroUsize,
) -> Result<()> {
    let realigned = scan_primaries(realigned_bam, bgzf_workers, None)?;
    let pass1 = scan_primaries(pass1_bam, bgzf_workers, Some(&realigned))?;

    let mut names: Vec<&Vec<u8>> = realigned.keys().collect();
    names.sort_unstable();

    let mut accepted: CandidateSet = FxHashSet::default();
    let mut counts: FxHashMap<&'static str, u64> = FxHashMap::default();
    let mut rows = Vec::with_capacity(names.len());
    for name in names {
        let new = &realigned[name];
        let old = pass1.get(name);
        let decision = decide(old, new, cfg);
        *counts.entry(decision.label()).or_default() += 1;
        if decision == Decision::Accepted {
            accepted.insert(name.clone());
        }
        rows.push((name, old, new, decision));
    }

    write_names(output, &accepted)?;

    if let Some(path) = report {
        let file =
            File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut w = BufWriter::new(file);
        writeln!(
            w,
            "name\tas_pass1\tas_realigned\tidentity_pass1\tidentity_realigned\tfive_clip_realigned\tthree_clip_realigned\tdecision"
        )?;
        let fmt = |v: Option<i64>| v.map_or_else(|| "NA".to_string(), |x| x.to_string());
        let fmt_id = |v: Option<f32>| v.map_or_else(|| "NA".to_string(), |x| format!("{x:.1}"));
        for (name, old, new, decision) in rows {
            let (five, three) = if new.is_reverse {
                (new.right_clip, new.left_clip)
            } else {
                (new.left_clip, new.right_clip)
            };
            w.write_all(name)?;
            writeln!(
                w,
                "\t{}\t{}\t{}\t{}\t{five}\t{three}\t{}",
                fmt(old.and_then(|o| o.alignment_score)),
                fmt(new.alignment_score),
                fmt_id(old.and_then(|o| o.identity)),
                fmt_id(new.identity),
                decision.label()
            )?;
        }
        w.flush()?;
    }

    info!(
        "[{}] accept done realigned={} accepted={} rejected_score={} rejected_quality={} rejected_missing_score={}",
        elapsed(),
        realigned.len(),
        counts.get("accepted").copied().unwrap_or(0),
        counts.get("rejected_score").copied().unwrap_or(0),
        counts.get("rejected_quality").copied().unwrap_or(0),
        counts.get("rejected_missing_score").copied().unwrap_or(0),
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record::cigar::Op;

    fn ops(items: &[(Kind, usize)]) -> Vec<std::io::Result<Op>> {
        items.iter().map(|(k, l)| Ok(Op::new(*k, *l))).collect()
    }

    fn score(as_: Option<i64>, identity: f32, lc: u32, rc: u32) -> PrimaryScore {
        PrimaryScore {
            alignment_score: as_,
            identity: Some(identity),
            left_clip: lc,
            right_clip: rc,
            is_reverse: false,
        }
    }

    const CFG: AcceptConfig = AcceptConfig {
        min_identity: 98.0,
        max_clip: 20,
    };

    #[test]
    fn accepts_higher_score_and_hq() {
        let old = score(Some(1500), 96.0, 0, 557);
        let new = score(Some(2100), 99.1, 3, 20);
        assert_eq!(decide(Some(&old), &new, &CFG), Decision::Accepted);
    }

    #[test]
    fn rejects_equal_or_lower_score() {
        let old = score(Some(2100), 99.0, 0, 0);
        assert_eq!(
            decide(Some(&old), &score(Some(2100), 99.5, 0, 0), &CFG),
            Decision::LowerScore
        );
        assert_eq!(
            decide(Some(&old), &score(Some(2000), 99.5, 0, 0), &CFG),
            Decision::LowerScore
        );
    }

    #[test]
    fn rejects_low_quality_realignment() {
        let old = score(Some(1000), 99.0, 0, 0);
        assert_eq!(
            decide(Some(&old), &score(Some(1200), 97.9, 0, 0), &CFG),
            Decision::LowQuality
        );
        assert_eq!(
            decide(Some(&old), &score(Some(1200), 99.0, 21, 0), &CFG),
            Decision::LowQuality
        );
        assert_eq!(
            decide(Some(&old), &score(Some(1200), 99.0, 0, 21), &CFG),
            Decision::LowQuality
        );
        assert_eq!(
            decide(Some(&old), &score(Some(1200), 98.0, 20, 20), &CFG),
            Decision::Accepted
        );
    }

    #[test]
    fn rejects_missing_nm() {
        let old = score(Some(1000), 99.0, 0, 0);
        let mut new = score(Some(1200), 99.0, 0, 0);
        new.identity = None;
        assert_eq!(decide(Some(&old), &new, &CFG), Decision::LowQuality);
    }

    #[test]
    fn rejects_missing_scores() {
        let new = score(Some(1200), 99.0, 0, 0);
        assert_eq!(decide(None, &new, &CFG), Decision::MissingScore);
        assert_eq!(
            decide(Some(&score(None, 99.0, 0, 0)), &new, &CFG),
            Decision::MissingScore
        );
        assert_eq!(
            decide(
                Some(&score(Some(1), 99.0, 0, 0)),
                &score(None, 99.0, 0, 0),
                &CFG
            ),
            Decision::MissingScore
        );
    }

    #[test]
    fn identity_matches_iso_segment() {
        // 5S 90= 2I 5X 100N 3D 10= 7S: aligned 105, indels 5, NM = 5 + 5 = 10
        // mismatches = 10 - 5 = 5, matches = 100, length 110 -> 90.9
        let (id, l, r) = identity_and_clips(
            ops(&[
                (Kind::SoftClip, 5),
                (Kind::SequenceMatch, 90),
                (Kind::Insertion, 2),
                (Kind::SequenceMismatch, 5),
                (Kind::Skip, 100),
                (Kind::Deletion, 3),
                (Kind::SequenceMatch, 10),
                (Kind::SoftClip, 7),
            ]),
            10,
        )
        .unwrap();
        assert_eq!((id, l, r), (90.9, 5, 7));
    }

    #[test]
    fn identity_rounds_like_iso_segment() {
        // 979 matches of 1000 columns = 97.9; 9795/10000 rounds half away to 98.0
        let (id, _, _) = identity_and_clips(ops(&[(Kind::Match, 1000)]), 21).unwrap();
        assert_eq!(id, 97.9);
        let (id, _, _) = identity_and_clips(ops(&[(Kind::Match, 10_000)]), 205).unwrap();
        assert_eq!(id, 98.0);
    }
}
