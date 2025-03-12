use dashmap::DashMap;
use hashbrown::HashSet;
use serde::{Deserialize, Serialize};

use std::borrow::Borrow;
use std::env;
use std::str::from_utf8_unchecked;

mod fns;
pub use fns::*;

mod mods;
pub use mods::*;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// numeric values
pub const SCALE: u64 = 100000000000; // 100Gb
pub const MIN_THREADS: usize = 1;
pub const MIN_BED_FIELDS: usize = 12;
pub const MIN_BED4_FIELDS: usize = 4;

// truncation numeric values
pub const TRUNCATION_THRESHOLD: f32 = 0.5;
pub const TRUNCATION_RECOVERY_THRESHOLD: f32 = 0.5;

// intron-retention numeric values || intron-classification numeric values
pub const RETENTION_RATIO_THRESHOLD: f32 = 0.001; // WARN: allowing everthing to enter recover step
pub const INTRON_RETENTION_RECOVERY_THRESHOLD: f32 = 0.5;
pub const INTRON_FREQUENCY_RECOVERY_THRESHOLD: f64 = 0.5;
pub const SPLICE_AI_SCORE_RECOVERY_THRESHOLD: f32 = 0.01; // INFO: if both splice sites are above, is a true intron
pub const MAX_ENT_SCORE_RECOVERY_THRESHOLD: f32 = 1.5; // INFO: if both splice sites are above, is a true intron

// fusion numeric values
pub const FUSION_RATIO_THRESHOLD: f32 = 0.5;

// polya numeric values
pub const INTRAPRIMING_RATIO_THRESHOLD: f32 = 0.5;

// file names
pub const INTRON_RETENTIONS: &str = "intron.retentions.bed";
pub const INTRON_RETENTION_FREE: &str = "intron.retentions.free.bed";
pub const TRUNCATIONS: &str = "truncations.bed";
pub const TRUNCATION_FREE: &str = "truncations.free.bed";
pub const BED3: &str = "ir.bed";
pub const INTERGENIC_REGIONS: &str = "intergenic.bed";
pub const FUSIONS: &str = "fusions.bed";
pub const FUSION_FREE: &str = "fusions.free.bed";
pub const FUSION_REVIEW: &str = "fusions.review.bed";
pub const FUSION_FAKES: &str = "fusions.fakes.bed";
pub const INTRON_CLASSIFICATION: &str = "reference_introns.tsv";
pub const MAXENTSCAN_ACCEPTOR_DB: &str = "db.tsv";
pub const MAXENTSCAN_DONOR_DB: &str = "donor.tsv";
pub const ORF_ASSIGNED_READS: &str = "orf_reads.bed";
pub const INTRAPRIMING_REVIEW: &str = "intrapriming_review.bed";
pub const POLYA_PASS: &str = "polya_pass.bed";
pub const POLYA_INTRAPRIMING: &str = "polya_intrapriming.bed";

// spliceai-related names
pub const ACCEPTOR_MINUS: &str = "spliceAiAcceptorMinus.bw";
pub const ACCEPTOR_PLUS: &str = "spliceAiAcceptorPlus.bw";
pub const DONOR_MINUS: &str = "spliceAiDonorMinus.bw";
pub const DONOR_PLUS: &str = "spliceAiDonorPlus.bw";

// flags
pub const COLORIZE: bool = false;

// collections
pub const COMPLEMENT: [u8; 128] = {
    let mut nt = [0; 128];
    nt[b'A' as usize] = b'T';
    nt[b'T' as usize] = b'A';
    nt[b'C' as usize] = b'G';
    nt[b'G' as usize] = b'C';
    nt[b'a' as usize] = b't';
    nt[b't' as usize] = b'a';
    nt[b'c' as usize] = b'g';
    nt[b'g' as usize] = b'c';
    nt[b'N' as usize] = b'N';
    nt[b'n' as usize] = b'n';
    nt
};

pub const BGD: [f64; 128] = {
    let mut bgd = [0.0; 128];
    bgd[b'A' as usize] = 0.27;
    bgd[b'T' as usize] = 0.27;
    bgd[b'C' as usize] = 0.23;
    bgd[b'G' as usize] = 0.23;
    bgd
};

pub const CONS1: [f64; 128] = {
    let mut bgd = [0.0; 128];
    bgd[b'A' as usize] = 0.9903;
    bgd[b'C' as usize] = 0.0032;
    bgd[b'G' as usize] = 0.0034;
    bgd[b'T' as usize] = 0.0030;
    bgd
};

pub const CONS2: [f64; 128] = {
    let mut bgd = [0.0; 128];
    bgd[b'A' as usize] = 0.0027;
    bgd[b'C' as usize] = 0.0037;
    bgd[b'G' as usize] = 0.9905;
    bgd[b'T' as usize] = 0.0030;
    bgd
};

// dirnames
pub const CLASSIFY_ASSETS: &str = "assets";

// types
pub type SpliceMap = (StrandSpliceMap, StrandSpliceMap);
pub type StrandSpliceMap = DashMap<String, DashMap<usize, f32>>;
pub type SharedSpliceMap = (Option<DashMap<usize, f32>>, Option<DashMap<usize, f32>>);
pub type SpliceScores = (Vec<StrandSpliceMap>, Vec<StrandSpliceMap>);

// traits
pub trait BedParser: Send + Sync + Sized {
    fn parse(
        line: &str,
        overlap: OverlapType,
        is_ref: bool,
    ) -> Result<Self, Box<dyn std::error::Error>>
    where
        Self: Sized;
    fn chrom(&self) -> &str;
    fn coord(&self) -> (u64, u64);
    fn intronic_coords(&self) -> HashSet<(u64, u64)>; // WARN: will not work for Bed[4,6,8]
    fn exonic_coords(&self) -> HashSet<(u64, u64)>; // WARN: will not work for Bed[4,6,8]
    fn name(&self) -> &str; // WARN: will not work with Bed[3]
    fn strand(&self) -> Strand;
    fn score(&self) -> f32;
    fn start(&self) -> u64;
    fn end(&self) -> u64;
    fn cds_start(&self) -> u64;
    fn cds_end(&self) -> u64;
    fn block_sizes(&self) -> Vec<u64>;
    fn block_starts(&self) -> Vec<u64>;
    fn block_count(&self) -> u64;
    fn rgb(&self) -> &str;
}

// public enums
pub enum SpliceSite {
    Donor,
    Acceptor,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum OverlapType {
    CDS,      // CDS-overlap
    Exon,     // exon-overlap
    Boundary, // boundary-overlap
    CDSBound, // CDS-overlap-UTR-boundary
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum SupportType {
    Splicing,
    RT,
    Unclear,
}

impl std::fmt::Display for SupportType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SupportType::Splicing => write!(f, "SPLICED"),
            SupportType::RT => write!(f, "RT"),
            SupportType::Unclear => write!(f, "UNCLEAR"),
        }
    }
}

impl std::str::FromStr for SupportType {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "SPLICED" => Ok(SupportType::Splicing),
            "RT" => Ok(SupportType::RT),
            "UNCLEAR" => Ok(SupportType::Unclear),
            _ => Err("ERROR: Cannot parse support type!".into()),
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum CoordType {
    Bounds,
    Intronic,
    Exonic,
}

#[derive(Debug, PartialEq, Clone, Hash, Eq, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
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

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum MatchType {
    Intron,
    SpliceSite,
}

impl Default for MatchType {
    fn default() -> Self {
        MatchType::SpliceSite
    }
}

impl From<bool> for MatchType {
    fn from(value: bool) -> Self {
        match value {
            true => MatchType::Intron,
            false => MatchType::SpliceSite,
        }
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum BedColumn {
    Chrom,
    Start,
    End,
    Name,
    Score,
    Strand,
    ThickStart,
    ThickEnd,
    ItemRgb,
    BlockCount,
    BlockSizes,
    BlockStarts,
}

#[derive(Debug, PartialEq, Clone)]
pub enum BedColumnValue {
    Chrom(String),
    Start(u64),
    End(u64),
    Name(String),
    Score(Vec<f32>), // WARN: trick for duplicated rows!
    Strand(Strand),
    ThickStart(u64),
    ThickEnd(u64),
    ItemRgb(String),
    BlockCount(u64),
    BlockSizes(Vec<u64>),
    BlockStarts(Vec<u64>),
}

impl BedColumnValue {
    /// Get the max score (if applicable)
    ///
    /// # Example
    /// ```rust, no_run
    /// use iso::BedColumnValue;
    ///
    /// let score = BedColumnValue::Score(vec![1.0, 2.0, 3.0]);
    /// assert_eq!(score.max_score(), Some(3.0));
    ///
    /// let score = BedColumnValue::Score(vec![1.0]);
    /// assert_eq!(score.max_score(), Some(1.0));
    /// ```
    pub fn max_score(&self) -> Option<f32> {
        match self {
            BedColumnValue::Score(scores) => scores
                .iter()
                .cloned()
                .max_by(|a, b| a.partial_cmp(b).unwrap()),
            _ => None,
        }
    }

    /// Get the average score (if applicable)
    ///
    /// # Example
    /// ```rust, no_run
    /// use iso::BedColumnValue;
    /// let score = BedColumnValue::Score(vec![1.0, 2.0, 3.0]);
    ///
    /// assert_eq!(score.avg_score(), Some(2.0));
    /// ```
    pub fn avg_score(&self) -> Option<f32> {
        match self {
            BedColumnValue::Score(scores) => {
                if scores.is_empty() {
                    None
                } else {
                    Some(scores.iter().sum::<f32>() / scores.len() as f32)
                }
            }
            _ => None,
        }
    }
}

// public structs
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Sequence {
    pub seq: String,
}

impl Sequence {
    pub fn new(seq: &[u8]) -> Self {
        Self {
            seq: unsafe { from_utf8_unchecked(seq).to_string() },
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn as_bytes(&self) -> &[u8] {
        self.seq.as_bytes()
    }

    pub fn as_str(&self) -> &str {
        self.seq.as_str()
    }

    pub fn to_string(&self) -> String {
        self.seq.clone()
    }

    pub fn to_uppercase(&self) -> String {
        self.seq.to_uppercase()
    }

    pub fn to_lowercase(&self) -> String {
        self.seq.to_lowercase()
    }

    pub fn reverse_complement(&self) -> Self {
        let mut rev = self.seq.chars().rev().collect::<String>();
        rev.make_ascii_uppercase();
        rev = rev
            .chars()
            .map(|c| COMPLEMENT[c as usize] as char)
            .collect::<String>();

        Self { seq: rev }
    }

    pub fn slice(&self, start: usize, end: usize) -> String {
        self.seq[start..end].to_string()
    }

    pub fn slice_as_seq(&self, start: usize, end: usize) -> Self {
        Self {
            seq: self.seq[start..end].to_string(),
        }
    }

    pub fn slice_as_bytes(&self, start: usize, end: usize) -> &[u8] {
        self.seq[start..end].as_bytes()
    }

    pub fn at_as_bytes(&self, idx: usize) -> usize {
        self.seq.as_bytes()[idx] as usize
    }

    pub fn fill(&self, kmer: usize) -> String {
        let mut seq = "A".repeat(kmer);
        seq.push_str(self.seq.as_str());

        seq
    }

    pub fn skip(&self, from: usize, to: usize) -> Sequence {
        Sequence {
            seq: self.seq[..from].to_string() + &self.seq[to..],
        }
    }
}

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.seq)
    }
}

impl Borrow<String> for Sequence {
    fn borrow(&self) -> &String {
        &self.seq
    }
}
