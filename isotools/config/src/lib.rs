use dashmap::DashMap;
use hashbrown::HashSet;
use serde::{Deserialize, Serialize};

use std::borrow::Borrow;
use std::env;
use std::str::from_utf8_unchecked;
use std::str::FromStr;

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
/// Trait for parsing bed files
///
/// # Example
///
/// ```rust, no_run
/// use iso::BedParser;
///
/// #[derive(Debug)]
/// struct Bed6 {
///    chrom: String,
///    start: u64,
///    end: u64,
///    name: String,
///    score: f32,
///    strand: Strand,
/// }
///
/// impl BedParser for Bed6 {
///   fn parse(
///     line: &str,
///     overlap: OverlapType,
///     is_ref: bool,
/// ) -> Result<Self, Box<dyn std::error::Error>> {
///   let fields: Vec<&str> = line.split('\t').collect();
///
///    Ok(Self {
///     chrom: fields[0].to_string(),
///     start: fields[1].parse::<u64>()?,
///     end: fields[2].parse::<u64>()?,
///     name: fields[3].to_string(),
///     score: fields[4].parse::<f32>()?,
///     strand: fields[5].parse::<Strand>()?,
///     })
/// }
///
/// fn chrom(&self) -> &str {
///     &self.chrom
/// }
/// ```
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

/// Tab delimited parser
///
/// # Example
/// ```rust, no_run
/// use iso::TsvParser;
///
/// #[derive(Debug)]
/// struct IsoformParser {
///    fields: Vec<String>,
/// }
///
/// impl TsvParser for IsoformParser {
///    fn parse(line: &str) -> Result<Self, anyhow::Error> {
///       Ok(Self {
///          fields: line.split('\t').map(|s| s.to_string()).collect(),
///      })
///   }
///
///  fn key(&self, index: usize) -> &str {
///     &self.fields[index]
/// }
///
/// fn value<V: FromStr>(&self, index: usize) -> Result<V, V::Err> {
///    self.fields[index].parse::<V>()
/// }
/// }
/// ```
pub trait TsvParser {
    fn parse(line: &str) -> Result<Self, anyhow::Error>
    where
        Self: Sized;
    fn key(&self, index: usize) -> &str;
    fn value<V: FromStr>(&self, index: usize) -> Result<V, V::Err>; // INFO:value column can be any type that implements FromStr
}

// public enums
/// Splice site type
///
/// This enum is used to store the type of splice site.
///
/// # Example
///
/// ```rust, no_run
/// use iso::SpliceSite;
///
/// let donor = SpliceSite::Donor;
/// let acceptor = SpliceSite::Acceptor;
/// ```
pub enum SpliceSite {
    Donor,
    Acceptor,
}

/// Overlap type
///
/// This enum is used to store the type of overlap.
///
/// # Example
///
/// ```rust, no_run
/// use iso::OverlapType;
///
/// let cds = OverlapType::CDS;
/// let exon = OverlapType::Exon;
/// let boundary = OverlapType::Boundary;
/// let cds_bound = OverlapType::CDSBound;
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum OverlapType {
    CDS,      // CDS-overlap
    Exon,     // exon-overlap
    Boundary, // boundary-overlap
    CDSBound, // CDS-overlap-UTR-boundary
}

impl From<&str> for OverlapType {
    fn from(value: &str) -> Self {
        match value {
            "cds" => OverlapType::CDS,
            "exon" => OverlapType::Exon,
            "bounds" => OverlapType::Boundary,
            "cds-bounded" => OverlapType::CDSBound,
            _ => panic!("ERROR: Cannot parse overlap type!"),
        }
    }
}

/// Support type
///
/// This enum is used to store the type of support.
///
/// # Example
///
/// ```rust, no_run
/// use iso::SupportType;
///
/// let spliced = SupportType::Splicing;
/// let rt = SupportType::RT;
/// let unclear = SupportType::Unclear;
/// ```
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

/// Coordinate type
///
/// This enum is used to store the type of coordinate.
///
/// # Example
///
/// ```rust, no_run
/// use iso::CoordType;
///
/// let bounds = CoordType::Bounds;
/// let intronic = CoordType::Intronic;
/// let exonic = CoordType::Exonic;
/// ```
#[derive(Debug, PartialEq, Clone)]
pub enum CoordType {
    Bounds,
    Intronic,
    Exonic,
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

/// Match type
///
/// This enum is used to store the type of match.
///
/// # Example
///
/// ```rust, no_run
/// use iso::MatchType;
///
/// let intron = MatchType::Intron;
/// let splice_site = MatchType::SpliceSite;
/// ```
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

/// Bed column names
///
/// This enum is used to store the names of a bed file.
///
/// # Example
///
/// ```rust, no_run
/// use iso::BedColumn;
///
/// let chrom = BedColumn::Chrom;
/// let start = BedColumn::Start;
/// ```
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

/// Bed column values
///
/// This enum is used to store the values of a bed file.
///
/// # Example
///
/// ```rust, no_run
/// use iso::BedColumnValue;
///
/// let chrom = BedColumnValue::Chrom(String::from("chr1"));
/// let start = BedColumnValue::Start(1);
/// ```
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
            seq: unsafe { from_utf8_unchecked(seq).to_string() },
        }
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
        self.seq[start..end].as_bytes()
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

    /// Fill the sequence with a given kmer
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"TCG");
    /// assert_eq!(seq.fill(2), "AATCG");
    /// ```
    pub fn fill(&self, kmer: usize) -> String {
        let mut seq = "A".repeat(kmer);
        seq.push_str(self.seq.as_str());

        seq
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
