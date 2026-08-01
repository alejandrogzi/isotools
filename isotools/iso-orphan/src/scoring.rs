// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Scoring module for evaluating transcript-level evidence.
//!
//! Provides reusable functions for:
//! - Junction-based scoring of multi-exon reads against reference transcripts
//! - Reciprocal overlap scoring for single-exon reads
//! - Coherent structural rescue against a single reference transcript
//! - Intron chain extraction and comparison for clustering
//! - Reference-overlap checking for guided evaluation
//! - Leave-one-out junction support counting for de novo evaluation
//!
//! Two principles run through the module. Evidence must be a *complete* feature — a
//! whole splice junction, a whole exon, both transcript ends — because a single shared
//! coordinate is far too easy to produce by coincidence. And evidence must come from a
//! single coherent source: one reference transcript, or reads other than the one being
//! judged.

use genepred::GenePred;

/// Parameters controlling scoring and classification thresholds.
#[derive(Debug, Clone)]
pub struct ScoringParams {
    /// Tolerance in bp for matching splice junctions (default: 0 = exact)
    pub junction_tolerance: u64,
    /// Minimum fraction of query junctions that must match a reference (default: 0.5)
    pub min_junction_frac: f64,
    /// Minimum reciprocal overlap for single-exon reads (default: 0.5)
    pub min_overlap_frac: f64,
    /// Minimum reads per de novo intron-chain cluster to keep (default: 5)
    ///
    /// This is a *structural* support threshold: it counts reads reproducing the same
    /// intron chain. It is never applied to the size of the packed component, which is a
    /// computational grouping property rather than evidence for a transcript.
    pub min_cluster_support: usize,
    /// Minimum reads per de novo single-exon overlap cluster to keep (default: 5)
    ///
    /// Single-exon reads carry no junction that could provide orthogonal evidence, so
    /// abundance stays a terminal criterion for them.
    pub min_single_exon_support: usize,
    /// Tolerance for terminal transcript-end matching in bp (default: 0 = exact)
    pub end_tolerance: u64,
    /// Minimum fraction of a read's introns that must be individually supported (default: 0.5)
    pub min_intron_support_frac: f64,
    /// Minimum fraction of the *other* reads of the group carrying an intron for it to
    /// count as "supported" (default: 0.5)
    pub intron_support_threshold: f64,
    /// Minimum splice-site score from bigWig to keep a read (default: 0.5)
    ///
    /// Compared against the median across the read's splice sites, which does not let a
    /// single uncovered site speak for the whole read.
    pub min_splice_score: f32,
}

impl Default for ScoringParams {
    fn default() -> Self {
        Self {
            junction_tolerance: 0,
            min_junction_frac: 0.5,
            min_overlap_frac: 0.5,
            min_cluster_support: 5,
            min_single_exon_support: 5,
            end_tolerance: 0,
            min_intron_support_frac: 0.5,
            intron_support_threshold: 0.5,
            min_splice_score: 0.5,
        }
    }
}

/// Result of scoring a multi-exon query against a reference transcript by junction comparison.
#[derive(Debug, Clone)]
pub struct JunctionScore {
    pub junction_matches: usize,
    pub query_junction_count: usize,
    pub ref_junction_count: usize,
    pub query_junction_frac: f64,
    pub ref_junction_frac: f64,
}

impl JunctionScore {
    /// Whether this score meets the minimum junction fraction threshold.
    pub fn passes(&self, min_frac: f64) -> bool {
        self.query_junction_frac >= min_frac
    }
}

/// Score a multi-exon query against a single reference transcript by intron comparison.
///
/// For each query intron (donor-acceptor pair), checks if a matching reference intron
/// exists within `tolerance` bp at both sites. Each reference intron is matched at most once.
///
/// # Examples
///
/// ```text
/// Query:     [100-200]---[300-400]---[500-600]  → introns: (200,300), (400,500)
/// Reference: [100-200]---[300-400]---[500-600]  → introns: (200,300), (400,500)
/// → junction_matches=2, query_junction_frac=1.0, ref_junction_frac=1.0
/// ```
pub fn score_junctions(query: &GenePred, reference: &GenePred, tolerance: u64) -> JunctionScore {
    align_junctions(query, reference, tolerance).score()
}

/// Pairing of query introns to reference introns produced by [`align_junctions`].
///
/// Keeping the matched index pairs, rather than only their count, is what allows a
/// caller to ask whether the shared junctions appear in the same order in both
/// transcripts, which distinguishes a partial-transcript relationship from a chimera
/// that borrows junctions from unrelated positions.
#[derive(Debug, Clone)]
pub struct JunctionAlignment {
    /// Matched `(query intron index, reference intron index)` pairs, in query order.
    pub matches: Vec<(usize, usize)>,
    /// Total number of query introns.
    pub query_junction_count: usize,
    /// Total number of reference introns.
    pub ref_junction_count: usize,
}

impl JunctionAlignment {
    /// Number of matched query junctions.
    pub fn matched(&self) -> usize {
        self.matches.len()
    }

    /// Whether the matched junctions appear in the same order in query and reference.
    ///
    /// An empty alignment is not order consistent: there is nothing to be consistent
    /// about, and callers use this as evidence of a shared structure.
    pub fn is_order_consistent(&self) -> bool {
        if self.matches.is_empty() {
            return false;
        }
        self.matches
            .windows(2)
            .all(|pair| pair[0].1 < pair[1].1 && pair[0].0 < pair[1].0)
    }

    /// Collapse into the plain feature values of a [`JunctionScore`].
    pub fn score(&self) -> JunctionScore {
        let matched = self.matched();

        let query_junction_frac = if self.query_junction_count == 0 {
            0.0
        } else {
            matched as f64 / self.query_junction_count as f64
        };
        let ref_junction_frac = if self.ref_junction_count == 0 {
            0.0
        } else {
            matched as f64 / self.ref_junction_count as f64
        };

        JunctionScore {
            junction_matches: matched,
            query_junction_count: self.query_junction_count,
            ref_junction_count: self.ref_junction_count,
            query_junction_frac,
            ref_junction_frac,
        }
    }
}

/// Pair each query intron with a matching reference intron within `tolerance`.
///
/// Each reference intron is matched at most once, and matching requires both the donor
/// and the acceptor site to agree, so a match is a complete splice junction rather than
/// a single coordinate.
pub fn align_junctions(
    query: &GenePred,
    reference: &GenePred,
    tolerance: u64,
) -> JunctionAlignment {
    let q_introns = query.introns();
    let r_introns = reference.introns();

    let mut alignment = JunctionAlignment {
        matches: Vec::new(),
        query_junction_count: q_introns.len(),
        ref_junction_count: r_introns.len(),
    };

    if q_introns.is_empty() || r_introns.is_empty() {
        return alignment;
    }

    let mut used = vec![false; r_introns.len()];

    for (i, (qs, qe)) in q_introns.iter().enumerate() {
        for (j, (rs, re)) in r_introns.iter().enumerate() {
            if !used[j] && qs.abs_diff(*rs) <= tolerance && qe.abs_diff(*re) <= tolerance {
                alignment.matches.push((i, j));
                used[j] = true;
                break;
            }
        }
    }

    alignment
}

/// Best junction comparison of a multi-exon query against a set of reference transcripts.
///
/// Retains both the winning feature values and the identity of the reference transcript
/// that produced them, so guided classification and reporting share a single comparison.
#[derive(Debug, Clone)]
pub struct BestJunctionMatch {
    /// Name of the reference transcript that produced `score`, when it has one.
    pub reference_id: Option<Vec<u8>>,
    /// Junction features of the winning comparison.
    pub score: JunctionScore,
}

/// Score a multi-exon query against multiple reference transcripts, returning the best match.
///
/// Only considers references with introns (multi-exon). Returns None if no multi-exon
/// references exist.
///
/// Ties are broken deterministically in favour of the first reference in input order:
/// equal `query_junction_frac` values never displace an already selected reference.
pub fn best_junction_match(
    query: &GenePred,
    references: &[GenePred],
    tolerance: u64,
) -> Option<BestJunctionMatch> {
    let mut best: Option<BestJunctionMatch> = None;

    for reference in references.iter().filter(|r| !r.introns().is_empty()) {
        let score = score_junctions(query, reference, tolerance);

        let is_better = match &best {
            None => true,
            Some(current) => score.query_junction_frac > current.score.query_junction_frac,
        };

        if is_better {
            best = Some(BestJunctionMatch {
                reference_id: reference.name().map(|name| name.to_vec()),
                score,
            });
        }
    }

    best
}

/// Score a multi-exon query against multiple reference transcripts, returning the best score.
///
/// Thin wrapper over [`best_junction_match`] for callers that do not need the reference
/// identity.
pub fn best_junction_score(
    query: &GenePred,
    references: &[GenePred],
    tolerance: u64,
) -> Option<JunctionScore> {
    best_junction_match(query, references, tolerance).map(|best| best.score)
}

/// Structural relationship that can rescue a query failing junction-fraction scoring.
///
/// Every variant describes a *complete* feature shared with **one** reference
/// transcript. A single matching coordinate is deliberately not representable: it can be
/// produced by two incompatible references, by an internal boundary unrelated to the
/// read's structure, or by coincidence, and is far weaker than the name "boundary match"
/// suggests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryEvidence {
    /// At least one complete splice junction is shared, and every shared junction
    /// appears in the same order in both transcripts.
    SpliceJunction,
    /// Both transcript ends agree with the same reference within `end_tolerance`.
    BothTranscriptEnds,
    /// A whole exon agrees with a reference exon at both of its boundaries.
    ExonPair,
}

impl BoundaryEvidence {
    /// Stable textual representation used in the report.
    pub fn as_str(&self) -> &'static str {
        match self {
            BoundaryEvidence::SpliceJunction => "SPLICE_JUNCTION",
            BoundaryEvidence::BothTranscriptEnds => "BOTH_TRANSCRIPT_ENDS",
            BoundaryEvidence::ExonPair => "EXON_PAIR",
        }
    }

    /// Relative strength, used to pick the best evidence across references.
    fn rank(&self) -> u8 {
        match self {
            BoundaryEvidence::SpliceJunction => 3,
            BoundaryEvidence::BothTranscriptEnds => 2,
            BoundaryEvidence::ExonPair => 1,
        }
    }
}

/// A coherent structural relationship between a query and a single reference transcript.
#[derive(Debug, Clone)]
pub struct BoundaryMatch {
    /// Name of the reference transcript that carries the whole relationship.
    pub reference_id: Option<Vec<u8>>,
    /// The strongest relationship found with that reference.
    pub evidence: BoundaryEvidence,
}

/// Strongest coherent relationship between `query` and any *single* multi-exon reference.
///
/// Each reference is evaluated on its own, so a rescue can never be assembled from
/// coordinates contributed by different, mutually incompatible transcripts. Splice
/// junctions and exon boundaries are compared within `junction_tolerance`; transcript
/// ends within `end_tolerance`.
///
/// Single-exon references are skipped: a spliced read shares no junction with an
/// unspliced transcript, so any agreement is a coordinate coincidence rather than
/// support for the read's structure. Such reads are meant to fall through to de novo
/// evaluation instead.
///
/// Ties are broken deterministically: stronger evidence wins, and on equal strength the
/// first reference in input order wins.
pub fn best_boundary_match(
    query: &GenePred,
    references: &[GenePred],
    junction_tolerance: u64,
    end_tolerance: u64,
) -> Option<BoundaryMatch> {
    let q_exons = query.exons();
    let mut best: Option<BoundaryMatch> = None;

    for reference in references.iter().filter(|r| !r.introns().is_empty()) {
        let evidence =
            if align_junctions(query, reference, junction_tolerance).is_order_consistent() {
                Some(BoundaryEvidence::SpliceJunction)
            } else if query.start().abs_diff(reference.start()) <= end_tolerance
                && query.end().abs_diff(reference.end()) <= end_tolerance
            {
                Some(BoundaryEvidence::BothTranscriptEnds)
            } else if shares_whole_exon(&q_exons, reference, junction_tolerance) {
                Some(BoundaryEvidence::ExonPair)
            } else {
                None
            };

        if let Some(evidence) = evidence {
            let is_better = match &best {
                None => true,
                Some(current) => evidence.rank() > current.evidence.rank(),
            };

            if is_better {
                best = Some(BoundaryMatch {
                    reference_id: reference.name().map(|name| name.to_vec()),
                    evidence,
                });
            }
        }
    }

    best
}

/// Whether any query exon agrees with a reference exon at both of its boundaries.
fn shares_whole_exon(q_exons: &[(u64, u64)], reference: &GenePred, tolerance: u64) -> bool {
    let r_exons = reference.exons();

    q_exons.iter().any(|(qs, qe)| {
        r_exons
            .iter()
            .any(|(rs, re)| qs.abs_diff(*rs) <= tolerance && qe.abs_diff(*re) <= tolerance)
    })
}

/// Compute reciprocal overlap between two intervals [a_start, a_end) and [b_start, b_end).
///
/// Returns min(overlap / a_length, overlap / b_length).
/// Returns 0.0 if intervals don't overlap or have zero length.
pub fn reciprocal_overlap(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> f64 {
    let overlap_start = a_start.max(b_start);
    let overlap_end = a_end.min(b_end);

    if overlap_start >= overlap_end {
        return 0.0;
    }

    let overlap = (overlap_end - overlap_start) as f64;
    let a_len = (a_end - a_start) as f64;
    let b_len = (b_end - b_start) as f64;

    if a_len == 0.0 || b_len == 0.0 {
        return 0.0;
    }

    (overlap / a_len).min(overlap / b_len)
}

/// Best single-exon overlap of a query against a set of reference transcripts.
///
/// Retains both the winning reciprocal overlap and the identity of the reference
/// transcript that produced it, so guided classification and reporting share a single
/// comparison.
#[derive(Debug, Clone)]
pub struct BestSingleExonMatch {
    /// Name of the reference transcript that produced `reciprocal_overlap`, when any
    /// reference exon overlaps the query at all.
    pub reference_id: Option<Vec<u8>>,
    /// Best reciprocal overlap found.
    pub reciprocal_overlap: f64,
}

/// Score a single-exon query against all reference exons by reciprocal overlap.
///
/// Returns the best reciprocal overlap of the query's single exon with any reference
/// exon, together with the reference transcript that produced it.
///
/// Ties are broken deterministically in favour of the first reference in input order,
/// and in favour of its first matching exon. A query that overlaps no reference exon
/// yields an overlap of `0.0` and no reference identity.
pub fn best_single_exon_match(query: &GenePred, references: &[GenePred]) -> BestSingleExonMatch {
    let mut best = BestSingleExonMatch {
        reference_id: None,
        reciprocal_overlap: 0.0,
    };

    let q_exons = query.exons();
    if q_exons.is_empty() {
        return best;
    }
    let (qs, qe) = q_exons[0];

    for reference in references {
        for (rs, re) in reference.exons() {
            let overlap = reciprocal_overlap(qs, qe, rs, re);
            if overlap > best.reciprocal_overlap {
                best.reciprocal_overlap = overlap;
                best.reference_id = reference.name().map(|name| name.to_vec());
            }
        }
    }

    best
}

/// Score a single-exon query against all reference exons by reciprocal overlap.
///
/// Thin wrapper over [`best_single_exon_match`] for callers that do not need the
/// reference identity.
pub fn best_single_exon_overlap(query: &GenePred, references: &[GenePred]) -> f64 {
    best_single_exon_match(query, references).reciprocal_overlap
}

/// Extract the intron chain from a record (list of intron intervals).
pub fn intron_chain(record: &GenePred) -> Vec<(u64, u64)> {
    record.introns()
}

/// Check if two intron chains match within a given tolerance at each junction.
pub fn intron_chains_match(a: &[(u64, u64)], b: &[(u64, u64)], tolerance: u64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter().zip(b.iter()).all(|((as_, ae), (bs, be))| {
        as_.abs_diff(*bs) <= tolerance && ae.abs_diff(*be) <= tolerance
    })
}

// ---------------------------------------------------------------------------
// Reference-overlap helpers
// ---------------------------------------------------------------------------

/// Check if any exon of the query overlaps any exon of any reference transcript.
///
/// This is the guided-vs-de-novo partition test: a read whose exons touch at least
/// one reference exon is scored in guided mode; a read with no exon overlap at all
/// is redirected to de novo evaluation.
pub fn overlaps_any_reference_exon(read: &GenePred, references: &[GenePred]) -> bool {
    let read_exons = read.exons();
    if read_exons.is_empty() {
        return false;
    }
    for (rs, re) in &read_exons {
        for reference in references {
            for (es, ee) in reference.exons() {
                // Two intervals overlap iff start_a < end_b AND start_b < end_a
                if *rs < ee && es < *re {
                    return true;
                }
            }
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Per-intron support helpers
// ---------------------------------------------------------------------------

/// Per-junction support counted across the multi-exon reads of one group.
///
/// Two properties matter for this to be evidence rather than an artefact:
///
/// * **Tolerance consistency.** Introns are grouped under a canonical representative
///   using the same `tolerance` that clusters intron chains, so two reads that cluster
///   together also support each other's junctions. Counting raw coordinates would make
///   the tolerance apply to one criterion and not the other.
/// * **Leave-one-out counting.** A read never counts towards the support of its own
///   junctions. Otherwise a read in a small group certifies itself: with two multi-exon
///   reads, a junction unique to one of them reaches a support fraction of 0.5 and
///   clears the default threshold on its own evidence alone.
#[derive(Debug, Clone)]
pub struct IntronSupport {
    /// Canonical intron representatives, sorted ascending.
    canonical: Vec<(u64, u64)>,
    /// Number of reads carrying each canonical intron, indexed as `canonical`.
    counts: Vec<usize>,
    /// Canonical intron indices carried by each read, indexed as the input reads.
    per_read: Vec<Vec<usize>>,
    /// Tolerance used to canonicalize junctions.
    tolerance: u64,
}

impl IntronSupport {
    /// Count junction support across `reads`, grouping junctions within `tolerance`.
    pub fn build(reads: &[&GenePred], tolerance: u64) -> Self {
        // INFO: canonical representatives are assigned in sorted coordinate order so the
        // grouping does not depend on the order reads arrive in
        let mut all: Vec<(u64, u64)> = reads.iter().flat_map(|read| read.introns()).collect();
        all.sort_unstable();
        all.dedup();

        let mut canonical: Vec<(u64, u64)> = Vec::new();
        for intron in all {
            let matches_existing = canonical
                .last()
                .is_some_and(|last| within_tolerance(*last, intron, tolerance));

            if !matches_existing {
                canonical.push(intron);
            }
        }

        let mut counts = vec![0usize; canonical.len()];
        let mut per_read: Vec<Vec<usize>> = Vec::with_capacity(reads.len());

        for read in reads {
            let mut carried: Vec<usize> = read
                .introns()
                .iter()
                .filter_map(|intron| Self::lookup(&canonical, *intron, tolerance))
                .collect();

            carried.sort_unstable();
            carried.dedup();

            for &idx in &carried {
                counts[idx] += 1;
            }

            per_read.push(carried);
        }

        Self {
            canonical,
            counts,
            per_read,
            tolerance,
        }
    }

    /// Index of the canonical representative of `intron`, if any.
    ///
    /// `canonical` is sorted, so only the window of candidates whose start lies within
    /// `tolerance` of the query needs to be examined.
    fn lookup(canonical: &[(u64, u64)], intron: (u64, u64), tolerance: u64) -> Option<usize> {
        let lower = intron.0.saturating_sub(tolerance);
        let upper = intron.0.saturating_add(tolerance);
        let from = canonical.partition_point(|candidate| candidate.0 < lower);

        canonical[from..]
            .iter()
            .take_while(|candidate| candidate.0 <= upper)
            .position(|candidate| within_tolerance(*candidate, intron, tolerance))
            .map(|offset| from + offset)
    }

    /// Number of reads carrying the canonical representative of `intron`.
    pub fn support_of(&self, intron: (u64, u64)) -> usize {
        Self::lookup(&self.canonical, intron, self.tolerance)
            .map(|idx| self.counts[idx])
            .unwrap_or(0)
    }

    /// Number of distinct canonical junctions in the group.
    pub fn junction_count(&self) -> usize {
        self.canonical.len()
    }

    /// Fraction of read `idx`'s junctions supported by the *other* reads of the group.
    ///
    /// A junction counts as supported when at least `threshold` of the other reads carry
    /// it. Returns `None` when the read has no junction, or when there is no other read
    /// to draw support from: absent evidence is not the same as evidence of absence, and
    /// the caller decides what to do with an unavailable value.
    pub fn support_fraction(&self, idx: usize, threshold: f64) -> Option<f64> {
        let carried = self.per_read.get(idx)?;
        if carried.is_empty() {
            return None;
        }

        let others = self.per_read.len().checked_sub(1)?;
        if others == 0 {
            return None;
        }

        let supported = carried
            .iter()
            .filter(|&&intron| {
                // INFO: leave-one-out — the read itself is excluded from its own support
                let by_others = self.counts[intron].saturating_sub(1);
                (by_others as f64 / others as f64) >= threshold
            })
            .count();

        Some(supported as f64 / carried.len() as f64)
    }
}

/// Whether two introns agree at both boundaries within `tolerance`.
fn within_tolerance(a: (u64, u64), b: (u64, u64), tolerance: u64) -> bool {
    a.0.abs_diff(b.0) <= tolerance && a.1.abs_diff(b.1) <= tolerance
}

/// Compute the median of a mutable f32 slice. Returns None for empty slices.
pub fn median_f32(values: &mut [f32]) -> Option<f32> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = values.len() / 2;
    if values.len().is_multiple_of(2) {
        Some((values[mid - 1] + values[mid]) / 2.0)
    } else {
        Some(values[mid])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn make_record(
        name: &[u8],
        start: u64,
        end: u64,
        block_starts: Vec<u64>,
        block_ends: Vec<u64>,
    ) -> GenePred {
        GenePred {
            chrom: b"chr1".to_vec(),
            start,
            end,
            name: Some(name.to_vec()),
            strand: None,
            thick_start: Some(start),
            thick_end: Some(end),
            block_count: Some(block_starts.len() as u32),
            block_starts: Some(block_starts),
            block_ends: Some(block_ends),
            extras: HashMap::new(),
        }
    }

    // -- reciprocal overlap tests --

    #[test]
    fn test_reciprocal_overlap_full() {
        assert_eq!(reciprocal_overlap(100, 200, 100, 200), 1.0);
    }

    #[test]
    fn test_reciprocal_overlap_half() {
        assert_eq!(reciprocal_overlap(100, 200, 150, 250), 0.5);
    }

    #[test]
    fn test_reciprocal_overlap_none() {
        assert_eq!(reciprocal_overlap(100, 200, 300, 400), 0.0);
    }

    #[test]
    fn test_reciprocal_overlap_contained() {
        let ov = reciprocal_overlap(100, 400, 150, 200);
        assert!((ov - 50.0 / 300.0).abs() < 1e-10);
    }

    // -- junction scoring tests --

    #[test]
    fn test_score_junctions_exact_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let score = score_junctions(&query, &reference, 0);
        assert_eq!(score.junction_matches, 2);
        assert_eq!(score.query_junction_frac, 1.0);
        assert_eq!(score.ref_junction_frac, 1.0);
        assert!(score.passes(0.5));
    }

    #[test]
    fn test_score_junctions_partial_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(
            b"r1",
            100,
            800,
            vec![100, 300, 500, 700],
            vec![200, 400, 600, 800],
        );

        let score = score_junctions(&query, &reference, 0);
        assert_eq!(score.junction_matches, 2);
        assert_eq!(score.query_junction_frac, 1.0);
        assert!((score.ref_junction_frac - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_score_junctions_no_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(b"r1", 100, 600, vec![100, 350, 550], vec![250, 450, 600]);

        let score = score_junctions(&query, &reference, 0);
        assert_eq!(score.junction_matches, 0);
        assert!(!score.passes(0.5));
    }

    #[test]
    fn test_score_junctions_with_tolerance() {
        let query = make_record(b"q1", 100, 600, vec![100, 302, 500], vec![200, 398, 600]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let score = score_junctions(&query, &reference, 3);
        assert_eq!(score.junction_matches, 2);
        assert!(score.passes(0.5));
    }

    #[test]
    fn test_best_junction_score_prefers_best_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let ref1 = make_record(b"r1", 100, 600, vec![100, 350, 550], vec![250, 450, 600]);
        let ref2 = make_record(b"r2", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let best = best_junction_score(&query, &[ref1, ref2], 0).unwrap();
        assert_eq!(best.junction_matches, 2);
        assert_eq!(best.query_junction_frac, 1.0);
    }

    #[test]
    fn test_best_junction_score_no_multi_exon_refs() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let se_ref = make_record(b"r1", 100, 600, vec![100], vec![600]);

        assert!(best_junction_score(&query, &[se_ref], 0).is_none());
    }

    #[test]
    fn test_best_single_exon_overlap_exact() {
        let query = make_record(b"q1", 100, 200, vec![100], vec![200]);
        let ref1 = make_record(b"r1", 100, 200, vec![100], vec![200]);
        let ref2 = make_record(b"r2", 300, 400, vec![300], vec![400]);

        assert_eq!(best_single_exon_overlap(&query, &[ref1, ref2]), 1.0);
    }

    // -- best-reference identity tests --

    #[test]
    fn test_best_junction_match_reports_reference_id() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let ref1 = make_record(b"r1", 100, 600, vec![100, 350, 550], vec![250, 450, 600]);
        let ref2 = make_record(b"r2", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let best = best_junction_match(&query, &[ref1, ref2], 0).unwrap();

        assert_eq!(best.reference_id.as_deref(), Some(b"r2".as_slice()));
        assert_eq!(best.score.junction_matches, 2);
        assert_eq!(best.score.query_junction_frac, 1.0);
    }

    #[test]
    fn test_best_junction_match_tie_break_is_first_reference() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        // both references reproduce the full query intron chain
        let ref1 = make_record(b"first", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let ref2 = make_record(
            b"second",
            100,
            800,
            vec![100, 300, 500, 700],
            vec![200, 400, 600, 800],
        );

        let best = best_junction_match(&query, &[ref1.clone(), ref2.clone()], 0).unwrap();
        assert_eq!(best.reference_id.as_deref(), Some(b"first".as_slice()));

        // stable under input order: the first reference in input order always wins
        let best = best_junction_match(&query, &[ref2, ref1], 0).unwrap();
        assert_eq!(best.reference_id.as_deref(), Some(b"second".as_slice()));
    }

    #[test]
    fn test_best_junction_match_none_without_multi_exon_refs() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let se_ref = make_record(b"r1", 100, 600, vec![100], vec![600]);

        assert!(best_junction_match(&query, &[se_ref], 0).is_none());
    }

    #[test]
    fn test_best_junction_score_agrees_with_best_junction_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let ref1 = make_record(b"r1", 100, 600, vec![100, 350, 550], vec![250, 450, 600]);
        let ref2 = make_record(b"r2", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let refs = vec![ref1, ref2];

        let score = best_junction_score(&query, &refs, 0).unwrap();
        let matched = best_junction_match(&query, &refs, 0).unwrap();

        assert_eq!(score.junction_matches, matched.score.junction_matches);
        assert_eq!(score.query_junction_frac, matched.score.query_junction_frac);
    }

    #[test]
    fn test_best_single_exon_match_reports_reference_id() {
        let query = make_record(b"q1", 100, 200, vec![100], vec![200]);
        let ref1 = make_record(b"r1", 300, 400, vec![300], vec![400]);
        let ref2 = make_record(b"r2", 100, 200, vec![100], vec![200]);

        let best = best_single_exon_match(&query, &[ref1, ref2]);

        assert_eq!(best.reciprocal_overlap, 1.0);
        assert_eq!(best.reference_id.as_deref(), Some(b"r2".as_slice()));
    }

    #[test]
    fn test_best_single_exon_match_without_overlap_has_no_reference() {
        let query = make_record(b"q1", 100, 200, vec![100], vec![200]);
        let ref1 = make_record(b"r1", 300, 400, vec![300], vec![400]);

        let best = best_single_exon_match(&query, &[ref1]);

        assert_eq!(best.reciprocal_overlap, 0.0);
        assert!(best.reference_id.is_none());
    }

    #[test]
    fn test_best_single_exon_match_tie_break_is_first_reference() {
        let query = make_record(b"q1", 100, 200, vec![100], vec![200]);
        let ref1 = make_record(b"first", 100, 200, vec![100], vec![200]);
        let ref2 = make_record(b"second", 100, 200, vec![100], vec![200]);

        let best = best_single_exon_match(&query, &[ref1, ref2]);

        assert_eq!(best.reciprocal_overlap, 1.0);
        assert_eq!(best.reference_id.as_deref(), Some(b"first".as_slice()));
    }

    // -- coherent boundary rescue tests --

    #[test]
    fn test_boundary_rescue_on_a_shared_junction() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(b"r1", 100, 400, vec![100, 300], vec![200, 400]);

        let rescue = best_boundary_match(&query, &[reference], 0, 0).unwrap();

        assert_eq!(rescue.evidence, BoundaryEvidence::SpliceJunction);
        assert_eq!(rescue.reference_id.as_deref(), Some(b"r1".as_slice()));
    }

    #[test]
    fn test_boundary_rescue_on_both_transcript_ends() {
        // no shared junction, but the same locus within tolerance
        let query = make_record(b"q1", 105, 595, vec![105, 350, 450], vec![250, 420, 595]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let refs = std::slice::from_ref(&reference);
        assert_eq!(
            best_boundary_match(&query, refs, 0, 10).unwrap().evidence,
            BoundaryEvidence::BothTranscriptEnds
        );
        // outside the end tolerance nothing rescues it
        assert!(best_boundary_match(&query, refs, 0, 0).is_none());
    }

    #[test]
    fn test_boundary_rescue_on_a_whole_shared_exon() {
        // reuses the whole internal exon [300,400) and shares nothing else
        let query = make_record(b"q1", 300, 5100, vec![300, 5000], vec![400, 5100]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let rescue = best_boundary_match(&query, &[reference], 0, 0).unwrap();
        assert_eq!(rescue.evidence, BoundaryEvidence::ExonPair);
    }

    #[test]
    fn test_boundary_rescue_rejects_a_single_shared_coordinate() {
        // 300 coincides with a reference exon start; no complete feature is shared
        let query = make_record(b"q1", 300, 5100, vec![300, 5000], vec![450, 5100]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        assert!(best_boundary_match(&query, &[reference], 0, 0).is_none());
    }

    #[test]
    fn test_boundary_rescue_rejects_a_single_shared_transcript_end() {
        // only the 5' end agrees, the 3' end is 4.5 kb away
        let query = make_record(b"q1", 100, 5100, vec![100, 5000], vec![150, 5100]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        assert!(best_boundary_match(&query, &[reference], 0, 10).is_none());
    }

    #[test]
    fn test_boundary_rescue_skips_single_exon_references() {
        // a spliced read shares both ends with an unspliced reference: not structural
        // support for its junctions
        let query = make_record(b"q1", 1000, 2000, vec![1000, 1500], vec![1200, 2000]);
        let reference = make_record(b"r1", 1000, 2000, vec![1000], vec![2000]);

        assert!(best_boundary_match(&query, &[reference], 0, 0).is_none());
    }

    #[test]
    fn test_boundary_rescue_never_pools_two_references() {
        // ref_a supplies the query's 5' end, ref_b its 3' end; neither supplies both and
        // no single reference shares a junction or a whole exon with the query
        let ref_a = make_record(b"ref_a", 100, 300, vec![100, 250], vec![150, 300]);
        let ref_b = make_record(b"ref_b", 4000, 5100, vec![4000, 5000], vec![4200, 5100]);
        let query = make_record(b"q1", 100, 5100, vec![100, 2000], vec![900, 5100]);

        assert!(best_boundary_match(&query, &[ref_a, ref_b], 0, 0).is_none());
    }

    #[test]
    fn test_boundary_rescue_prefers_the_strongest_evidence() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        // exon-pair only
        let weak = make_record(b"weak", 3000, 3400, vec![3000, 100], vec![3200, 200]);
        // shares a complete junction
        let strong = make_record(b"strong", 100, 400, vec![100, 300], vec![200, 400]);

        let rescue = best_boundary_match(&query, &[weak, strong], 0, 0).unwrap();
        assert_eq!(rescue.evidence, BoundaryEvidence::SpliceJunction);
        assert_eq!(rescue.reference_id.as_deref(), Some(b"strong".as_slice()));
    }

    #[test]
    fn test_junction_alignment_order_consistency() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let alignment = align_junctions(&query, &reference, 0);
        assert_eq!(alignment.matched(), 2);
        assert!(alignment.is_order_consistent());

        // no shared junction at all is not order consistent
        let unrelated = make_record(b"r2", 5000, 5600, vec![5000, 5300], vec![5200, 5600]);
        assert!(!align_junctions(&query, &unrelated, 0).is_order_consistent());
    }

    #[test]
    fn test_intron_chains_match_exact() {
        assert!(intron_chains_match(
            &[(200, 300), (400, 500)],
            &[(200, 300), (400, 500)],
            0
        ));
    }

    #[test]
    fn test_intron_chains_match_different_length() {
        assert!(!intron_chains_match(
            &[(200, 300)],
            &[(200, 300), (400, 500)],
            0
        ));
    }

    #[test]
    fn test_intron_chains_match_with_tolerance() {
        let a = vec![(200, 302), (398, 500)];
        let b = vec![(200, 300), (400, 500)];
        assert!(intron_chains_match(&a, &b, 3));
        assert!(!intron_chains_match(&a, &b, 1));
    }

    #[test]
    fn test_junction_score_beats_weak_overlap() {
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let good_query = make_record(b"q1", 110, 590, vec![110, 300, 500], vec![200, 400, 590]);
        let bad_query = make_record(b"q2", 100, 600, vec![100, 350, 550], vec![250, 450, 600]);

        let good_score = score_junctions(&good_query, &reference, 0);
        let bad_score = score_junctions(&bad_query, &reference, 0);

        assert!(good_score.query_junction_frac > bad_score.query_junction_frac);
        assert!(good_score.passes(0.5));
        assert!(!bad_score.passes(0.5));
    }

    #[test]
    fn test_multi_isoform_reference() {
        let ref1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let ref2 = make_record(b"r2", 100, 800, vec![100, 350, 700], vec![300, 600, 800]);
        let query = make_record(b"q1", 100, 800, vec![100, 350, 700], vec![300, 600, 800]);

        let best = best_junction_score(&query, &[ref1, ref2], 0).unwrap();
        assert_eq!(best.junction_matches, 2);
        assert_eq!(best.query_junction_frac, 1.0);
        assert!(best.passes(0.5));
    }

    #[test]
    fn test_single_exon_vs_multi_exon_ref() {
        let query = make_record(b"q1", 300, 400, vec![300], vec![400]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let overlap = best_single_exon_overlap(&query, &[reference]);
        assert_eq!(overlap, 1.0);
    }

    // -- reference overlap tests --

    #[test]
    fn test_overlaps_ref_exon_exact() {
        let read = make_record(b"r", 300, 400, vec![300], vec![400]);
        let refs = vec![make_record(
            b"ref",
            100,
            600,
            vec![100, 300, 500],
            vec![200, 400, 600],
        )];
        // Read exon [300,400) overlaps ref exon [300,400)
        assert!(overlaps_any_reference_exon(&read, &refs));
    }

    #[test]
    fn test_overlaps_ref_exon_partial() {
        let read = make_record(b"r", 150, 250, vec![150], vec![250]);
        let refs = vec![make_record(
            b"ref",
            100,
            600,
            vec![100, 300, 500],
            vec![200, 400, 600],
        )];
        // Read exon [150,250) overlaps ref exon [100,200)
        assert!(overlaps_any_reference_exon(&read, &refs));
    }

    #[test]
    fn test_no_overlap_intronic() {
        // Read sits entirely within a reference intron
        let read = make_record(b"r", 220, 280, vec![220], vec![280]);
        let refs = vec![make_record(
            b"ref",
            100,
            600,
            vec![100, 300, 500],
            vec![200, 400, 600],
        )];
        // Ref exons: [100,200), [300,400), [500,600). Read [220,280) is in intron [200,300).
        assert!(!overlaps_any_reference_exon(&read, &refs));
    }

    #[test]
    fn test_no_overlap_intergenic() {
        let read = make_record(b"r", 700, 800, vec![700], vec![800]);
        let refs = vec![make_record(
            b"ref",
            100,
            600,
            vec![100, 300, 500],
            vec![200, 400, 600],
        )];
        assert!(!overlaps_any_reference_exon(&read, &refs));
    }

    #[test]
    fn test_overlaps_multi_exon_read_partial() {
        // Multi-exon read: one exon overlaps reference, the other does not
        let read = make_record(b"r", 150, 700, vec![150, 650], vec![250, 700]);
        let refs = vec![make_record(
            b"ref",
            100,
            600,
            vec![100, 300, 500],
            vec![200, 400, 600],
        )];
        // Read exon [150,250) overlaps ref exon [100,200) → true
        assert!(overlaps_any_reference_exon(&read, &refs));
    }

    // -- intron support tests --

    #[test]
    fn test_intron_support_counts_each_junction() {
        // 3 reads sharing intron (200,300); 2 sharing (400,500); 1 unique (600,700)
        let r1 = make_record(
            b"r1",
            100,
            800,
            vec![100, 300, 500, 700],
            vec![200, 400, 600, 800],
        );
        let r2 = make_record(b"r2", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let r3 = make_record(b"r3", 100, 400, vec![100, 300], vec![200, 400]);

        let reads: Vec<&GenePred> = vec![&r1, &r2, &r3];
        let support = IntronSupport::build(&reads, 0);

        assert_eq!(support.support_of((200, 300)), 3);
        assert_eq!(support.support_of((400, 500)), 2);
        assert_eq!(support.support_of((600, 700)), 1);
        assert_eq!(support.support_of((999, 1999)), 0);
        assert_eq!(support.junction_count(), 3);
    }

    #[test]
    fn test_support_fraction_excludes_the_read_itself() {
        // two reads with completely different chains: neither corroborates the other
        let a = make_record(b"a", 100, 400, vec![100, 300], vec![200, 400]);
        let b = make_record(b"b", 100, 900, vec![100, 800], vec![700, 900]);

        let reads: Vec<&GenePred> = vec![&a, &b];
        let support = IntronSupport::build(&reads, 0);

        assert_eq!(
            support.support_fraction(0, 0.5),
            Some(0.0),
            "a read must not certify its own junctions"
        );
        assert_eq!(support.support_fraction(1, 0.5), Some(0.0));
    }

    #[test]
    fn test_support_fraction_is_unavailable_without_another_read() {
        let lone = make_record(b"lone", 100, 400, vec![100, 300], vec![200, 400]);
        let reads: Vec<&GenePred> = vec![&lone];
        let support = IntronSupport::build(&reads, 0);

        assert_eq!(
            support.support_fraction(0, 0.5),
            None,
            "no other read means support cannot be measured, not that it is zero"
        );
    }

    #[test]
    fn test_support_fraction_all_supported() {
        // five reads with an identical chain
        let reads_data: Vec<GenePred> = (0..5u64)
            .map(|i| {
                make_record(
                    format!("r{}", i).as_bytes(),
                    100 + i,
                    600,
                    vec![100 + i, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect();
        let reads: Vec<&GenePred> = reads_data.iter().collect();
        let support = IntronSupport::build(&reads, 0);

        assert_eq!(support.support_fraction(0, 0.5), Some(1.0));
    }

    #[test]
    fn test_support_fraction_partially_supported() {
        // six reads share (200,300) and (400,500), the variant shares only the first
        let mut reads_data: Vec<GenePred> = (0..6u64)
            .map(|i| {
                make_record(
                    format!("r{}", i).as_bytes(),
                    100 + i,
                    600,
                    vec![100 + i, 300, 500],
                    vec![200, 400, 600],
                )
            })
            .collect();
        reads_data.push(make_record(
            b"variant",
            100,
            650,
            vec![100, 300, 550],
            vec![200, 400, 650],
        ));

        let reads: Vec<&GenePred> = reads_data.iter().collect();
        let support = IntronSupport::build(&reads, 0);

        assert_eq!(support.support_fraction(6, 0.5), Some(0.5));
    }

    #[test]
    fn test_support_fraction_uses_the_tolerance() {
        // junctions differing by a few bases: grouped when the tolerance allows it
        let reads_data: Vec<GenePred> = (0..4u64)
            .map(|i| {
                make_record(
                    format!("r{}", i).as_bytes(),
                    100,
                    600,
                    vec![100, 300 + i],
                    vec![200 + i, 600],
                )
            })
            .collect();
        let reads: Vec<&GenePred> = reads_data.iter().collect();

        let exact = IntronSupport::build(&reads, 0);
        assert_eq!(exact.support_fraction(0, 0.5), Some(0.0));
        assert_eq!(exact.junction_count(), 4);

        let tolerant = IntronSupport::build(&reads, 5);
        assert_eq!(tolerant.support_fraction(0, 0.5), Some(1.0));
        assert_eq!(tolerant.junction_count(), 1);
    }

    #[test]
    fn test_support_fraction_threshold_is_applied() {
        // one read carries a junction that only one of four others shares
        let shared = make_record(b"shared", 100, 400, vec![100, 300], vec![200, 400]);
        let partner = make_record(b"partner", 100, 400, vec![100, 300], vec![200, 400]);
        let others: Vec<GenePred> = (0..3u64)
            .map(|i| {
                make_record(
                    format!("other{}", i).as_bytes(),
                    1000,
                    1400,
                    vec![1000, 1300],
                    vec![1200, 1400],
                )
            })
            .collect();

        let mut reads: Vec<&GenePred> = vec![&shared, &partner];
        reads.extend(others.iter());
        let support = IntronSupport::build(&reads, 0);

        // 1 of the 4 other reads carries it -> 0.25
        assert_eq!(support.support_fraction(0, 0.2), Some(1.0));
        assert_eq!(support.support_fraction(0, 0.5), Some(0.0));
    }

    #[test]
    fn test_support_is_independent_of_read_order() {
        let build = |order: [usize; 3]| -> Vec<f64> {
            let data: Vec<GenePred> = vec![
                make_record(b"a", 100, 600, vec![100, 302], vec![200, 600]),
                make_record(b"b", 100, 600, vec![100, 300], vec![202, 600]),
                make_record(b"c", 100, 600, vec![100, 301], vec![201, 600]),
            ];
            let reads: Vec<&GenePred> = order.iter().map(|&i| &data[i]).collect();
            let support = IntronSupport::build(&reads, 5);
            (0..3)
                .map(|i| support.support_fraction(i, 0.5).unwrap())
                .collect()
        };

        assert_eq!(build([0, 1, 2]), vec![1.0, 1.0, 1.0]);
        assert_eq!(build([2, 0, 1]), vec![1.0, 1.0, 1.0]);
    }

    // -- median tests --

    #[test]
    fn test_median_f32_odd() {
        assert_eq!(median_f32(&mut [3.0, 1.0, 2.0]), Some(2.0));
    }

    #[test]
    fn test_median_f32_even() {
        assert_eq!(median_f32(&mut [4.0, 1.0, 3.0, 2.0]), Some(2.5));
    }

    #[test]
    fn test_median_f32_single() {
        assert_eq!(median_f32(&mut [5.0]), Some(5.0));
    }

    #[test]
    fn test_median_f32_empty() {
        assert_eq!(median_f32(&mut []), None);
    }
}
