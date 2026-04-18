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

/// Pre-loaded splice-site scores for querying at specific genomic positions.
///
/// Stores scores above a minimum threshold from four SpliceAI BigWig files.
/// Provides a `median_splice_score` method for computing per-read splice quality.
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

    /// Query the best splice-site score at a given position and strand.
    ///
    /// Returns max(donor_score, acceptor_score) at that position, or 0.0 if absent.
    fn splice_score_at(&self, chrom: &str, pos: u64, is_forward: bool) -> f32 {
        let (donor_map, acceptor_map) = if is_forward {
            (&self.donor_fwd, &self.acceptor_fwd)
        } else {
            (&self.donor_rev, &self.acceptor_rev)
        };

        let d = donor_map
            .get(chrom)
            .and_then(|m| m.get(&pos).copied())
            .unwrap_or(0.0);
        let a = acceptor_map
            .get(chrom)
            .and_then(|m| m.get(&(pos - 1)).copied())
            .unwrap_or(0.0);
        d.max(a)
    }

    /// Compute median splice-site score across all intron boundaries of a record.
    ///
    /// For each intron (start, end), queries the splice score at both boundaries
    /// on the appropriate strand. Returns None if the record has no introns,
    /// no strand info, or unknown strand.
    pub fn median_splice_score(&self, record: &GenePred) -> Option<f32> {
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
        for &(start, end) in &introns {
            scores.push(self.splice_score_at(chrom, start, is_forward));
            scores.push(self.splice_score_at(chrom, end, is_forward));
        }

        log::debug!("DEBUG: splice scores for {:?}: {:?}", record.name(), scores);

        median_f32(&mut scores)
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

        assert_eq!(provider.splice_score_at("chr1", 200, true), 0.9);
        assert_eq!(provider.splice_score_at("chr1", 300, true), 0.8);
        assert_eq!(provider.splice_score_at("chr1", 200, false), 0.0); // no rev scores
        assert_eq!(provider.splice_score_at("chr1", 999, true), 0.0); // missing position
    }

    #[test]
    fn test_median_splice_score_forward_strand() {
        let donor_fwd = DashMap::new();
        let acceptor_fwd = DashMap::new();

        let mut d = HashMap::new();
        d.insert(200u64, 0.9f32); // intron start (donor)
        d.insert(400u64, 0.7f32); // intron start (donor)
        donor_fwd.insert("chr1".to_string(), d);

        let mut a = HashMap::new();
        a.insert(300u64, 0.8f32); // intron end (acceptor)
        a.insert(500u64, 0.6f32); // intron end (acceptor)
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

        // Scores: donor@200=0.9, acceptor@300=0.8, donor@400=0.7, acceptor@500=0.6
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
