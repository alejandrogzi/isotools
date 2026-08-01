// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Splice-site score provider using SpliceAI BigWig files.
//!
//! Discovers four required BigWig files (donor/acceptor × plus/minus strand),
//! pre-loads significant scores, and provides per-read median splice-site scoring.

use bigtools::{utils::reopen::Reopen, BigWigRead};
use dashmap::DashMap;
use genepred::{GenePred, Strand as GpStrand};
use log::info;
use rayon::prelude::*;

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::{fs, io};

use crate::scoring::median_f32;

// ---------------------------------------------------------------------------
// SpliceScoreProvider
// ---------------------------------------------------------------------------

/// Splice-site evidence for one read, aggregated across all of its junctions.
///
/// `median` is the value classification compares against `min_splice_score`. A minimum
/// would be brittle: BigWig coverage is not uniform, and a single site that happens to be
/// uncovered scores 0.0 and would sink a read whose other sites are all strong. The
/// median summarises the read instead of letting its worst site speak for it.
///
/// `weakest` is reported alongside as a diagnostic. It gates nothing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpliceEvidence {
    /// Median score across all donor and acceptor sites.
    pub median: f32,
    /// Lowest score across all donor and acceptor sites. Diagnostic only.
    pub weakest: f32,
    /// Number of scored sites, two per junction.
    pub sites: usize,
}

impl SpliceEvidence {
    /// Whether the read's median splice-site score reaches `min_score`.
    pub fn passes(&self, min_score: f32) -> bool {
        self.median >= min_score
    }
}

/// Pre-loaded splice-site scores for querying at specific genomic positions.
///
/// Stores scores above a minimum threshold from four SpliceAI BigWig files.
/// Provides a `splice_evidence` method for computing per-read splice quality.
pub struct SpliceScoreProvider {
    donor_fwd: DashMap<String, HashMap<u64, f32>>,
    donor_rev: DashMap<String, HashMap<u64, f32>>,
    acceptor_fwd: DashMap<String, HashMap<u64, f32>>,
    acceptor_rev: DashMap<String, HashMap<u64, f32>>,
}

impl SpliceScoreProvider {
    /// Load splice scores from a directory containing four SpliceAI BigWig files.
    ///
    /// Only positions with scores >= `min_score` are stored in memory.
    pub fn from_dir(dir: &Path, min_score: f32) -> Result<Self, SpliceAiFileError> {
        let resolved = discover_spliceai_bigwigs(dir)?;
        info!("Loading BigWig splice scores from {:?}", dir);

        let provider = Self {
            donor_fwd: DashMap::new(),
            donor_rev: DashMap::new(),
            acceptor_fwd: DashMap::new(),
            acceptor_rev: DashMap::new(),
        };

        rayon::scope(|s| {
            s.spawn(|_| load_bigwig_scores(&resolved.donor_plus, &provider.donor_fwd, min_score));
            s.spawn(|_| {
                load_bigwig_scores(&resolved.acceptor_plus, &provider.acceptor_fwd, min_score)
            });
            s.spawn(|_| load_bigwig_scores(&resolved.donor_minus, &provider.donor_rev, min_score));
            s.spawn(|_| {
                load_bigwig_scores(&resolved.acceptor_minus, &provider.acceptor_rev, min_score)
            });
        });

        let total: usize = [
            &provider.donor_fwd,
            &provider.donor_rev,
            &provider.acceptor_fwd,
            &provider.acceptor_rev,
        ]
        .iter()
        .map(|m| {
            let local_size = m.iter().map(|e| e.value().len()).sum::<usize>();

            log::info!("INFO: Loaded {} splice scores", local_size,);
            local_size
        })
        .sum();
        info!("Loaded {} significant splice-site scores", total);

        Ok(provider)
    }

    /// Create from pre-built maps (for testing).
    #[cfg(test)]
    pub(crate) fn from_maps(
        donor_fwd: DashMap<String, HashMap<u64, f32>>,
        donor_rev: DashMap<String, HashMap<u64, f32>>,
        acceptor_fwd: DashMap<String, HashMap<u64, f32>>,
        acceptor_rev: DashMap<String, HashMap<u64, f32>>,
    ) -> Self {
        Self {
            donor_fwd,
            donor_rev,
            acceptor_fwd,
            acceptor_rev,
        }
    }

    /// Donor score at a single genomic base.
    fn donor_score_at(&self, chrom: &str, base: u64, is_forward: bool) -> f32 {
        let map = if is_forward {
            &self.donor_fwd
        } else {
            &self.donor_rev
        };

        map.get(chrom)
            .and_then(|m| m.get(&base).copied())
            .unwrap_or(0.0)
    }

    /// Acceptor score at a single genomic base.
    fn acceptor_score_at(&self, chrom: &str, base: u64, is_forward: bool) -> f32 {
        let map = if is_forward {
            &self.acceptor_fwd
        } else {
            &self.acceptor_rev
        };

        map.get(chrom)
            .and_then(|m| m.get(&base).copied())
            .unwrap_or(0.0)
    }

    /// Donor and acceptor score of one intron, in that order.
    ///
    /// BigWig scores are 0-based per-base values while a BED intron `(start, end)` has an
    /// exclusive end, so the two intronic bases carrying splice signal are `start` and
    /// `end - 1`. Which of them is the donor depends on the strand: transcription runs
    /// left to right on the plus strand, so `start` is the donor and `end - 1` the
    /// acceptor; on the minus strand the roles are reversed.
    ///
    /// Each site is read from the track that can actually carry its signal. Taking the
    /// maximum of both tracks at both ends, as this previously did, let a strong donor
    /// stand in for a missing acceptor and read the acceptor track at a base belonging to
    /// the preceding exon.
    fn junction_scores(&self, chrom: &str, intron: (u64, u64), is_forward: bool) -> (f32, f32) {
        let (start, last) = (intron.0, intron.1.saturating_sub(1));

        if is_forward {
            (
                self.donor_score_at(chrom, start, true),
                self.acceptor_score_at(chrom, last, true),
            )
        } else {
            (
                self.donor_score_at(chrom, last, false),
                self.acceptor_score_at(chrom, start, false),
            )
        }
    }

    /// Splice-site evidence across every junction of a record.
    ///
    /// Returns None if the record has no introns, no strand info, or unknown strand.
    pub fn splice_evidence(&self, record: &GenePred) -> Option<SpliceEvidence> {
        let introns = record.introns();
        if introns.is_empty() {
            return None;
        }

        let chrom = std::str::from_utf8(record.chrom()).ok()?;
        let is_forward = match record.strand()? {
            GpStrand::Forward => true,
            GpStrand::Reverse => false,
            GpStrand::Unknown => return None,
        };

        let mut scores: Vec<f32> = Vec::with_capacity(introns.len() * 2);
        for &intron in &introns {
            let (donor, acceptor) = self.junction_scores(chrom, intron, is_forward);
            scores.push(donor);
            scores.push(acceptor);
        }

        log::debug!("DEBUG: splice scores for {:?}: {:?}", record.name(), scores);

        let weakest = scores.iter().copied().fold(f32::INFINITY, f32::min);
        let median = median_f32(&mut scores)?;

        Some(SpliceEvidence {
            median,
            weakest,
            sites: scores.len(),
        })
    }

    /// Median splice-site score across all junctions of a record.
    ///
    /// Thin wrapper over [`SpliceScoreProvider::splice_evidence`] for callers that only
    /// need the median.
    pub fn median_splice_score(&self, record: &GenePred) -> Option<f32> {
        self.splice_evidence(record).map(|evidence| evidence.median)
    }
}

/// Load scores from a single BigWig file into a DashMap.
fn load_bigwig_scores(path: &Path, target: &DashMap<String, HashMap<u64, f32>>, min_score: f32) {
    let bwread = BigWigRead::open_file(path).unwrap_or_else(|e| {
        panic!("ERROR: Cannot open BigWig file {:?}: {:?}", path, e);
    });
    let chroms: Vec<_> = bwread.chroms().to_vec();

    chroms.into_par_iter().for_each(|chr| {
        let mut reader = BigWigRead::reopen(&bwread).unwrap_or_else(|e| {
            panic!("ERROR: Cannot reopen BigWig file {:?}: {:?}", path, e);
        });

        let values = reader.values(&chr.name, 0, chr.length).unwrap_or_else(|e| {
            panic!(
                "ERROR: Cannot read values from BigWig {:?} chr {}: {:?}",
                path, chr.name, e
            );
        });

        let scores: HashMap<u64, f32> = values
            .into_iter()
            .enumerate()
            .filter(|(_, v)| *v >= min_score)
            .map(|(i, v)| (i as u64, v))
            .collect();

        if !scores.is_empty() {
            target.insert(chr.name.clone(), scores);
        }
    });
}

// ---------------------------------------------------------------------------
// BigWig file discovery (adapted from existing infrastructure)
// ---------------------------------------------------------------------------

/// DNA strand orientation for splice sites.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

impl Strand {
    fn token(self) -> &'static str {
        match self {
            Strand::Forward => "plus",
            Strand::Reverse => "minus",
        }
    }
}

/// Splice site type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpliceSite {
    Donor,
    Acceptor,
}

impl std::fmt::Display for SpliceSite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpliceSite::Donor => write!(f, "D"),
            SpliceSite::Acceptor => write!(f, "A"),
        }
    }
}

impl SpliceSite {
    fn token(self) -> &'static str {
        match self {
            SpliceSite::Donor => "donor",
            SpliceSite::Acceptor => "acceptor",
        }
    }
}

/// Classification of a SpliceAI BigWig file by splice site type and strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct BigWigClass {
    splice_site: SpliceSite,
    strand: Strand,
}

impl BigWigClass {
    const fn new(splice_site: SpliceSite, strand: Strand) -> Self {
        Self {
            splice_site,
            strand,
        }
    }

    fn label(self) -> &'static str {
        match (self.splice_site, self.strand) {
            (SpliceSite::Donor, Strand::Forward) => "donor plus",
            (SpliceSite::Donor, Strand::Reverse) => "donor minus",
            (SpliceSite::Acceptor, Strand::Forward) => "acceptor plus",
            (SpliceSite::Acceptor, Strand::Reverse) => "acceptor minus",
        }
    }
}

const EXPECTED_BIGWIGS: [BigWigClass; 4] = [
    BigWigClass::new(SpliceSite::Donor, Strand::Forward),
    BigWigClass::new(SpliceSite::Acceptor, Strand::Forward),
    BigWigClass::new(SpliceSite::Donor, Strand::Reverse),
    BigWigClass::new(SpliceSite::Acceptor, Strand::Reverse),
];

/// Resolved paths to the four required SpliceAI BigWig files.
#[derive(Debug, Clone)]
struct ResolvedBigWigs {
    donor_plus: PathBuf,
    acceptor_plus: PathBuf,
    donor_minus: PathBuf,
    acceptor_minus: PathBuf,
}

/// Errors raised while discovering the four required SpliceAI BigWig files.
#[derive(Debug)]
pub enum SpliceAiFileError {
    InvalidDirectory {
        path: PathBuf,
    },
    ReadDirectory {
        path: PathBuf,
        source: io::Error,
    },
    ReadDirectoryEntry {
        path: PathBuf,
        source: io::Error,
    },
    InvalidFilename {
        path: PathBuf,
        reason: &'static str,
    },
    DuplicateClassification {
        classification: &'static str,
        paths: Vec<PathBuf>,
    },
    MissingClassifications {
        path: PathBuf,
        classifications: Vec<&'static str>,
    },
}

impl std::fmt::Display for SpliceAiFileError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidDirectory { path } => {
                write!(
                    f,
                    "SpliceAI BigWig path '{}' is not a directory",
                    path.display()
                )
            }
            Self::ReadDirectory { path, source } => write!(
                f,
                "failed to read SpliceAI BigWig directory '{}': {}",
                path.display(),
                source
            ),
            Self::ReadDirectoryEntry { path, source } => write!(
                f,
                "failed to read an entry from SpliceAI BigWig directory '{}': {}",
                path.display(),
                source
            ),
            Self::InvalidFilename { path, reason } => write!(
                f,
                "invalid SpliceAI BigWig filename '{}': {}",
                path.display(),
                reason
            ),
            Self::DuplicateClassification {
                classification,
                paths,
            } => write!(
                f,
                "multiple SpliceAI BigWig files matched {}: {}",
                classification,
                display_paths(paths)
            ),
            Self::MissingClassifications {
                path,
                classifications,
            } => write!(
                f,
                "missing required SpliceAI BigWig files in '{}': {}",
                path.display(),
                classifications.join(", ")
            ),
        }
    }
}

impl std::error::Error for SpliceAiFileError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ReadDirectory { source, .. } => Some(source),
            Self::ReadDirectoryEntry { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Scans a directory and resolves exactly four SpliceAI BigWig files.
///
/// Filenames are matched after removing non-alphanumeric characters and lowercasing the stem, so
/// patterns like `prefixDonorPlusSuffix.bw`, `ACCEPTOR_MINUS.bigWig`, and
/// `run42spliceaiacceptorplusv2.BW` are all supported.
fn discover_spliceai_bigwigs(dir: &Path) -> Result<ResolvedBigWigs, SpliceAiFileError> {
    if !dir.is_dir() {
        return Err(SpliceAiFileError::InvalidDirectory {
            path: dir.to_path_buf(),
        });
    }

    let mut matches: HashMap<BigWigClass, Vec<PathBuf>> = HashMap::new();

    for entry in fs::read_dir(dir).map_err(|source| SpliceAiFileError::ReadDirectory {
        path: dir.to_path_buf(),
        source,
    })? {
        let entry = entry.map_err(|source| SpliceAiFileError::ReadDirectoryEntry {
            path: dir.to_path_buf(),
            source,
        })?;
        let path = entry.path();

        if !path.is_file() || !is_bigwig_path(&path) {
            continue;
        }

        if let Some(classification) = classify_bigwig_path(&path)? {
            matches.entry(classification).or_default().push(path);
        }
    }

    for classification in EXPECTED_BIGWIGS {
        if let Some(paths) = matches.get_mut(&classification) {
            paths.sort();

            if paths.len() > 1 {
                return Err(SpliceAiFileError::DuplicateClassification {
                    classification: classification.label(),
                    paths: paths.clone(),
                });
            }
        }
    }

    let missing = EXPECTED_BIGWIGS
        .iter()
        .copied()
        .filter(|classification| !matches.contains_key(classification))
        .map(BigWigClass::label)
        .collect::<Vec<_>>();

    if !missing.is_empty() {
        return Err(SpliceAiFileError::MissingClassifications {
            path: dir.to_path_buf(),
            classifications: missing,
        });
    }

    Ok(ResolvedBigWigs {
        donor_plus: take_classified_path(
            &mut matches,
            BigWigClass::new(SpliceSite::Donor, Strand::Forward),
        ),
        acceptor_plus: take_classified_path(
            &mut matches,
            BigWigClass::new(SpliceSite::Acceptor, Strand::Forward),
        ),
        donor_minus: take_classified_path(
            &mut matches,
            BigWigClass::new(SpliceSite::Donor, Strand::Reverse),
        ),
        acceptor_minus: take_classified_path(
            &mut matches,
            BigWigClass::new(SpliceSite::Acceptor, Strand::Reverse),
        ),
    })
}

fn take_classified_path(
    matches: &mut HashMap<BigWigClass, Vec<PathBuf>>,
    classification: BigWigClass,
) -> PathBuf {
    matches
        .remove(&classification)
        .and_then(|mut paths| paths.pop())
        .unwrap_or_else(|| {
            panic!(
                "ERROR: missing validated SpliceAI BigWig classification {}",
                classification.label()
            )
        })
}

fn classify_bigwig_path(path: &Path) -> Result<Option<BigWigClass>, SpliceAiFileError> {
    let stem = path
        .file_stem()
        .and_then(|value| value.to_str())
        .ok_or_else(|| SpliceAiFileError::InvalidFilename {
            path: path.to_path_buf(),
            reason: "filename stem is not valid UTF-8",
        })?;

    classify_bigwig_stem(path, stem)
}

fn classify_bigwig_stem(path: &Path, stem: &str) -> Result<Option<BigWigClass>, SpliceAiFileError> {
    let normalized = normalize_bigwig_stem(stem);
    let has_donor = normalized.contains(SpliceSite::Donor.token());
    let has_acceptor = normalized.contains(SpliceSite::Acceptor.token());
    let has_plus = normalized.contains(Strand::Forward.token());
    let has_minus = normalized.contains(Strand::Reverse.token());

    let site_count = usize::from(has_donor) + usize::from(has_acceptor);
    let strand_count = usize::from(has_plus) + usize::from(has_minus);

    if site_count == 0 && strand_count == 0 {
        return Ok(None);
    }

    if site_count != 1 || strand_count != 1 {
        return Err(SpliceAiFileError::InvalidFilename {
            path: path.to_path_buf(),
            reason: classification_error_reason(site_count, strand_count),
        });
    }

    let splice_site = if has_donor {
        SpliceSite::Donor
    } else {
        SpliceSite::Acceptor
    };
    let strand = if has_plus {
        Strand::Forward
    } else {
        Strand::Reverse
    };

    Ok(Some(BigWigClass::new(splice_site, strand)))
}

fn classification_error_reason(site_count: usize, strand_count: usize) -> &'static str {
    match (site_count, strand_count) {
        (0, 1) | (0, 2) => "missing donor/acceptor token",
        (1, 0) | (2, 0) => "missing plus/minus token",
        (2, 1) => "contains both donor and acceptor tokens",
        (1, 2) => "contains both plus and minus tokens",
        (2, 2) => "contains both donor/acceptor and plus/minus token pairs",
        _ => "must contain exactly one donor/acceptor token and one plus/minus token",
    }
}

fn normalize_bigwig_stem(stem: &str) -> String {
    stem.chars()
        .filter(|character| character.is_ascii_alphanumeric())
        .map(|character| character.to_ascii_lowercase())
        .collect()
}

fn is_bigwig_path(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| {
            extension.eq_ignore_ascii_case("bw") || extension.eq_ignore_ascii_case("bigwig")
        })
}

fn display_paths(paths: &[PathBuf]) -> String {
    paths
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        env, fs,
        fs::File,
        time::{SystemTime, UNIX_EPOCH},
    };

    struct TempDir {
        path: PathBuf,
    }

    impl TempDir {
        fn new() -> Self {
            let unique = format!(
                "splicing-spliceai-test-{}-{}",
                std::process::id(),
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            );
            let path = env::temp_dir().join(unique);
            fs::create_dir_all(&path).unwrap();

            Self { path }
        }

        fn path(&self) -> &Path {
            &self.path
        }

        fn touch(&self, filename: &str) {
            File::create(self.path.join(filename)).unwrap();
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.path);
        }
    }

    // -- discovery tests --

    #[test]
    fn classify_bigwig_stem_accepts_affixes_casing_and_extensions() {
        let donor_plus = classify_bigwig_path(Path::new("run42spliceAiDONORPLUSv2.bigWig"))
            .unwrap()
            .unwrap();
        let acceptor_minus = classify_bigwig_path(Path::new("prefix_acceptor-minus_suffix.BW"))
            .unwrap()
            .unwrap();

        assert_eq!(
            donor_plus,
            BigWigClass::new(SpliceSite::Donor, Strand::Forward)
        );
        assert_eq!(
            acceptor_minus,
            BigWigClass::new(SpliceSite::Acceptor, Strand::Reverse)
        );
    }

    #[test]
    fn classify_bigwig_stem_rejects_partial_matches() {
        let error = classify_bigwig_path(Path::new("spliceai_donor_only.bw")).unwrap_err();
        match error {
            SpliceAiFileError::InvalidFilename { reason, .. } => {
                assert_eq!(reason, "missing plus/minus token");
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn discover_spliceai_bigwigs_resolves_mixed_extensions_and_affixes() {
        let dir = TempDir::new();
        dir.touch("run42DonorPlusv2.BW");
        dir.touch("prefix_acceptor-plus_suffix.bigWig");
        dir.touch("sampleDONORMINUSextra.bw");
        dir.touch("xxAcceptorMinusyy.bigwig");

        let resolved = discover_spliceai_bigwigs(dir.path()).unwrap();

        assert_eq!(
            resolved.donor_plus.file_name().unwrap(),
            "run42DonorPlusv2.BW"
        );
        assert_eq!(
            resolved.acceptor_plus.file_name().unwrap(),
            "prefix_acceptor-plus_suffix.bigWig"
        );
        assert_eq!(
            resolved.donor_minus.file_name().unwrap(),
            "sampleDONORMINUSextra.bw"
        );
        assert_eq!(
            resolved.acceptor_minus.file_name().unwrap(),
            "xxAcceptorMinusyy.bigwig"
        );
    }

    #[test]
    fn discover_spliceai_bigwigs_rejects_duplicate_classifications() {
        let dir = TempDir::new();
        dir.touch("sample_donor_plus_v1.bw");
        dir.touch("sample_donor_plus_v2.bigwig");
        dir.touch("sample_acceptor_plus.bw");
        dir.touch("sample_donor_minus.bw");
        dir.touch("sample_acceptor_minus.bw");

        let error = discover_spliceai_bigwigs(dir.path()).unwrap_err();
        match error {
            SpliceAiFileError::DuplicateClassification {
                classification,
                paths,
            } => {
                assert_eq!(classification, "donor plus");
                assert_eq!(paths.len(), 2);
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn discover_spliceai_bigwigs_rejects_missing_classifications() {
        let dir = TempDir::new();
        dir.touch("sample_donor_plus.bw");
        dir.touch("sample_acceptor_plus.bw");
        dir.touch("sample_acceptor_minus.bw");

        let error = discover_spliceai_bigwigs(dir.path()).unwrap_err();
        match error {
            SpliceAiFileError::MissingClassifications {
                classifications, ..
            } => {
                assert_eq!(classifications, vec!["donor minus"]);
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    // -- SpliceScoreProvider query tests --

    #[test]
    fn test_splice_score_at_queries_correct_strand() {
        let donor_fwd = DashMap::new();
        let mut chr1 = HashMap::new();
        chr1.insert(200u64, 0.9f32);
        donor_fwd.insert("chr1".to_string(), chr1);

        let acceptor_fwd = DashMap::new();
        let mut chr1_acc = HashMap::new();
        chr1_acc.insert(300u64, 0.8f32);
        acceptor_fwd.insert("chr1".to_string(), chr1_acc);

        let provider =
            SpliceScoreProvider::from_maps(donor_fwd, DashMap::new(), acceptor_fwd, DashMap::new());

        // donor of intron (200,301) is its first base, acceptor its last one
        assert_eq!(
            provider.junction_scores("chr1", (200, 301), true),
            (0.9, 0.8)
        );
        // on the minus strand the two roles swap ends
        assert_eq!(
            provider.junction_scores("chr1", (200, 301), false),
            (0.0, 0.0),
            "no reverse-strand scores were supplied"
        );
        // a junction at unscored positions yields zeros rather than a panic
        assert_eq!(
            provider.junction_scores("chr1", (999, 1999), true),
            (0.0, 0.0)
        );
    }

    #[test]
    fn test_junction_scores_use_the_site_specific_track() {
        // a strong donor must not stand in for a missing acceptor
        let donor_fwd = DashMap::new();
        let mut chr1 = HashMap::new();
        chr1.insert(200u64, 0.9f32);
        donor_fwd.insert("chr1".to_string(), chr1);

        let provider = SpliceScoreProvider::from_maps(
            donor_fwd,
            DashMap::new(),
            DashMap::new(),
            DashMap::new(),
        );

        assert_eq!(
            provider.junction_scores("chr1", (200, 300), true),
            (0.9, 0.0)
        );
    }

    #[test]
    fn test_junction_scores_swap_roles_on_the_reverse_strand() {
        // on the minus strand the intron's last base is the donor
        let donor_rev = DashMap::new();
        let acceptor_rev = DashMap::new();
        let mut d = HashMap::new();
        let mut a = HashMap::new();
        d.insert(299u64, 0.7f32);
        a.insert(200u64, 0.6f32);
        donor_rev.insert("chr1".to_string(), d);
        acceptor_rev.insert("chr1".to_string(), a);

        let provider =
            SpliceScoreProvider::from_maps(DashMap::new(), donor_rev, DashMap::new(), acceptor_rev);

        assert_eq!(
            provider.junction_scores("chr1", (200, 300), false),
            (0.7, 0.6)
        );
    }

    #[test]
    fn test_splice_evidence_reports_median_and_weakest() {
        let donor_fwd = DashMap::new();
        let acceptor_fwd = DashMap::new();
        let mut d = HashMap::new();
        let mut a = HashMap::new();
        d.insert(200u64, 0.9f32);
        d.insert(400u64, 0.8f32);
        a.insert(299u64, 0.7f32);
        // acceptor of intron (400,500) is absent
        donor_fwd.insert("chr1".to_string(), d);
        acceptor_fwd.insert("chr1".to_string(), a);

        let provider =
            SpliceScoreProvider::from_maps(donor_fwd, DashMap::new(), acceptor_fwd, DashMap::new());

        let record = GenePred {
            chrom: b"chr1".to_vec(),
            start: 100,
            end: 600,
            name: Some(b"test".to_vec()),
            strand: Some(GpStrand::Forward),
            thick_start: Some(100),
            thick_end: Some(600),
            block_count: Some(3),
            block_starts: Some(vec![100, 300, 500]),
            block_ends: Some(vec![200, 400, 600]),
            extras: HashMap::new(),
        };

        let evidence = provider.splice_evidence(&record).unwrap();

        assert_eq!(evidence.sites, 4);
        // scores [0.9, 0.7, 0.8, 0.0] -> median (0.7 + 0.8) / 2
        assert!((evidence.median - 0.75).abs() < 1e-6);
        assert_eq!(evidence.weakest, 0.0);
        assert!(
            evidence.passes(0.5),
            "one uncovered site does not sink an otherwise supported read"
        );
    }

    #[test]
    fn test_median_decides_not_the_weakest_site() {
        // three junctions scoring 0.9/0.8, 0.9/0.6 and 0.9/absent: the median is 0.850
        // and the read is supported, even though its weakest site is 0.0
        let donor_fwd = DashMap::new();
        let acceptor_fwd = DashMap::new();
        let mut d = HashMap::new();
        let mut a = HashMap::new();
        d.insert(200u64, 0.9f32);
        d.insert(400u64, 0.9f32);
        d.insert(600u64, 0.9f32);
        a.insert(299u64, 0.8f32);
        a.insert(499u64, 0.6f32);
        // acceptor of intron (600,700) is outside BigWig coverage
        donor_fwd.insert("chr1".to_string(), d);
        acceptor_fwd.insert("chr1".to_string(), a);

        let provider =
            SpliceScoreProvider::from_maps(donor_fwd, DashMap::new(), acceptor_fwd, DashMap::new());

        let record = GenePred {
            chrom: b"chr1".to_vec(),
            start: 100,
            end: 800,
            name: Some(b"test".to_vec()),
            strand: Some(GpStrand::Forward),
            thick_start: Some(100),
            thick_end: Some(800),
            block_count: Some(4),
            block_starts: Some(vec![100, 300, 500, 700]),
            block_ends: Some(vec![200, 400, 600, 800]),
            extras: HashMap::new(),
        };

        let evidence = provider.splice_evidence(&record).unwrap();

        assert_eq!(evidence.sites, 6);
        // sorted [0.0, 0.6, 0.8, 0.9, 0.9, 0.9] -> median (0.8 + 0.9) / 2
        assert!((evidence.median - 0.85).abs() < 1e-6);
        assert_eq!(evidence.weakest, 0.0);
        assert!(evidence.passes(0.5));
    }

    #[test]
    fn test_median_splice_score_forward_strand() {
        let donor_fwd = DashMap::new();
        let acceptor_fwd = DashMap::new();

        // BigWig scores are per-base: the donor sits on the first intronic base and the
        // acceptor on the last one, so an intron (start, end) is scored at start and
        // at end - 1
        let mut d = HashMap::new();
        d.insert(200u64, 0.9f32); // first base of intron (200,300) (donor)
        d.insert(400u64, 0.7f32); // first base of intron (400,500) (donor)
        donor_fwd.insert("chr1".to_string(), d);

        let mut a = HashMap::new();
        a.insert(299u64, 0.8f32); // last base of intron (200,300) (acceptor)
        a.insert(499u64, 0.6f32); // last base of intron (400,500) (acceptor)
        acceptor_fwd.insert("chr1".to_string(), a);

        let provider =
            SpliceScoreProvider::from_maps(donor_fwd, DashMap::new(), acceptor_fwd, DashMap::new());

        // Read: 3 exons [100-200, 300-400, 500-600], introns (200,300) (400,500), forward
        let record = GenePred {
            chrom: b"chr1".to_vec(),
            start: 100,
            end: 600,
            name: Some(b"test".to_vec()),
            strand: Some(GpStrand::Forward),
            thick_start: Some(100),
            thick_end: Some(600),
            block_count: Some(3),
            block_starts: Some(vec![100, 300, 500]),
            block_ends: Some(vec![200, 400, 600]),
            extras: HashMap::new(),
        };

        // Scores: donor@200=0.9, acceptor@299=0.8, donor@400=0.7, acceptor@499=0.6
        // Sorted: [0.6, 0.7, 0.8, 0.9], median = (0.7 + 0.8) / 2 = 0.75
        let median = provider.median_splice_score(&record).unwrap();
        assert!((median - 0.75).abs() < 1e-6);
    }

    #[test]
    fn test_median_splice_score_no_strand_returns_none() {
        let provider = SpliceScoreProvider::from_maps(
            DashMap::new(),
            DashMap::new(),
            DashMap::new(),
            DashMap::new(),
        );

        let record = GenePred {
            chrom: b"chr1".to_vec(),
            start: 100,
            end: 600,
            name: Some(b"test".to_vec()),
            strand: None,
            thick_start: Some(100),
            thick_end: Some(600),
            block_count: Some(3),
            block_starts: Some(vec![100, 300, 500]),
            block_ends: Some(vec![200, 400, 600]),
            extras: HashMap::new(),
        };

        assert!(provider.median_splice_score(&record).is_none());
    }

    #[test]
    fn test_median_splice_score_single_exon_returns_none() {
        let provider = SpliceScoreProvider::from_maps(
            DashMap::new(),
            DashMap::new(),
            DashMap::new(),
            DashMap::new(),
        );

        let record = GenePred {
            chrom: b"chr1".to_vec(),
            start: 100,
            end: 200,
            name: Some(b"test".to_vec()),
            strand: Some(GpStrand::Forward),
            thick_start: Some(100),
            thick_end: Some(200),
            block_count: Some(1),
            block_starts: Some(vec![100]),
            block_ends: Some(vec![200]),
            extras: HashMap::new(),
        };

        assert!(provider.median_splice_score(&record).is_none());
    }
}
