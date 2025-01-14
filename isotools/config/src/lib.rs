use dashmap::{DashMap, DashSet};
use hashbrown::HashSet;
use indicatif::{ProgressBar, ProgressStyle};
use num_traits::{Num, NumCast};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use thiserror::Error;

use std::any::Any;
use std::borrow::Borrow;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::str::from_utf8_unchecked;
use std::time::Duration;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// numeric values
pub const SCALE: u64 = 100000000000; // 100Gb
pub const MIN_THREADS: usize = 1;
pub const MIN_BED_FIELDS: usize = 12;
pub const MIN_BED4_FIELDS: usize = 4;

// truncation numeric values
pub const TRUNCATION_THRESHOLD: f32 = 0.5;
pub const TRUNCATION_RECOVERY_THRESHOLD: f32 = 0.5;

// intron-retention numeric values
pub const RETENTION_RATIO_THRESHOLD: f32 = 0.001; // WARN: allowing everthing to enter recover step
pub const INTRON_RETENTION_RECOVERY_THRESHOLD: f32 = 0.5;
pub const SPLICE_AI_SCORE_RECOVERY_THRESHOLD: f32 = 0.01; // INFO: if both splice sites are above, is a true intron

// fusion numeric values
pub const FUSION_RATIO_THRESHOLD: f32 = 0.5;

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
pub const INTRON_CLASSIFICATION: &str = "reference_introns.bed";
pub const MAXENTSCAN_ACCEPTOR_DB: &str = "db.tsv";
pub const MAXENTSCAN_DONOR_DB: &str = "donor.tsv";

// spliceai-related names
pub const ACCEPTOR_MINUS: &str = "spliceAiAcceptorMinus.bw";
pub const ACCEPTOR_PLUS: &str = "spliceAiAcceptorPlus.bw";
pub const DONOR_MINUS: &str = "spliceAiDonorMinus.bw";
pub const DONOR_PLUS: &str = "spliceAiDonorPlus.bw";

// flags
pub const COLORIZE: bool = false;
pub const OVERLAP_CDS: bool = false;
pub const OVERLAP_EXON: bool = true;

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
pub trait BedRecord: Send + Sync {
    fn parse(line: String) -> Result<Self, Box<dyn std::error::Error>>
    where
        Self: Sized;
    fn chrom(&self) -> &str;
    fn coord(&self) -> (u64, u64);
    fn intronic_coords(&self) -> HashSet<(u64, u64)>; // WARN: will not work for Bed[4,6,8]
    fn exonic_coords(&self) -> HashSet<(u64, u64)>; // WARN: will not work for Bed[4,6,8]
}

// public enums
pub enum SpliceSite {
    Donor,
    Acceptor,
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

// os
#[cfg(not(windows))]
const TICK_SETTINGS: (&str, u64) = ("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏ ", 80);
#[cfg(windows)]
const TICK_SETTINGS: (&str, u64) = (r"+-x| ", 200);

/// return a pre-configured progress bar
pub fn get_progress_bar(length: u64, msg: &str) -> ProgressBar {
    let progressbar_style = ProgressStyle::default_spinner()
        .tick_chars(TICK_SETTINGS.0)
        .template(" {spinner} {msg:<30} {wide_bar} ETA {eta_precise} ")
        .expect("no template error");

    let progress_bar = ProgressBar::new(length);

    progress_bar.set_style(progressbar_style);
    progress_bar.enable_steady_tick(Duration::from_millis(TICK_SETTINGS.1));
    progress_bar.set_message(msg.to_owned());

    progress_bar
}

/// write a DashSet to a file
pub fn write_objs<T>(data: &DashSet<T>, fname: &str)
where
    T: AsRef<str> + Sync + Send + Eq + std::hash::Hash,
{
    log::info!("Reads in {}: {:?}. Writing...", fname, data.len());
    let f = match File::create(fname) {
        Ok(f) => f,
        Err(e) => panic!("Error creating file: {}", e),
    };
    let mut writer = BufWriter::new(f);

    for line in data.iter() {
        writeln!(writer, "{}", line.as_ref()).unwrap_or_else(|e| {
            panic!("Error writing to file: {}", e);
        });
    }
}

/// write any collection to a file
pub fn write_collection(data: &Vec<String>, fname: &str) {
    log::info!("Reads in {}: {:?}. Writing...", fname, data.len());
    let f = match File::create(fname) {
        Ok(f) => f,
        Err(e) => panic!("Error creating file: {}", e),
    };
    let mut writer = BufWriter::new(f);

    for line in data.iter() {
        writeln!(writer, "{}", line).unwrap_or_else(|e| {
            panic!("Error writing to file: {}", e);
        });
    }
}

/// argument checker for all subcommands
pub trait ArgCheck {
    fn check(&self) -> Result<(), CliError> {
        self.validate_args()
    }

    fn validate_args(&self) -> Result<(), CliError> {
        self.check_dbs()?;

        if !self.get_blacklist().is_empty() {
            self.check_blacklist()?;
        } else {
            log::warn!("No blacklist provided. Skipping...");
        };

        Ok(())
    }

    fn check_dbs(&self) -> Result<(), CliError> {
        if self.get_ref().is_empty() {
            let err = "No reference files provided".to_string();
            return Err(CliError::InvalidInput(err));
        }
        for db in self.get_ref() {
            validate(db)?;
        }

        if self.get_query().is_empty() {
            let err = "No query file provided".to_string();
            return Err(CliError::InvalidInput(err));
        }
        for query in self.get_query() {
            validate(query)?;
        }

        Ok(())
    }

    fn check_blacklist(&self) -> Result<(), CliError> {
        for bl in self.get_blacklist() {
            validate(bl)?;
        }
        Ok(())
    }

    fn get_blacklist(&self) -> &Vec<PathBuf>;
    fn get_ref(&self) -> &Vec<PathBuf>;
    fn get_query(&self) -> &Vec<PathBuf>;
}

/// error handling for CLI
#[derive(Debug, Error)]
pub enum CliError {
    #[error("Invalid input: {0}")]
    InvalidInput(String),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

/// argument validation
pub fn validate(arg: &PathBuf) -> Result<(), CliError> {
    if !arg.exists() {
        return Err(CliError::InvalidInput(format!(
            "ERROR: {:?} does not exist",
            arg
        )));
    }

    if !arg.is_file() {
        return Err(CliError::InvalidInput(format!(
            "ERROR: {:?} is not a file",
            arg
        )));
    }

    match arg.extension() {
        Some(ext) if ext == "bed" => (),
        _ => {
            return Err(CliError::InvalidInput(format!(
                "ERROR: file {:?} is not a BED file",
                arg
            )))
        }
    }

    match std::fs::metadata(arg) {
        Ok(metadata) if metadata.len() == 0 => Err(CliError::InvalidInput(format!(
            "ERROR: file {:?} is empty",
            arg
        ))),
        Ok(_) => Ok(()),
        Err(e) => Err(CliError::IoError(e)),
    }
}

// module descriptors
#[derive(Debug)]
pub enum ModuleType {
    IntronRetention,
    StartTruncation,
    FusionDetection,
}

pub trait ModuleMap: Any {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value>;
    fn set_value(&mut self, key: Box<dyn Any>, value: serde_json::Value) -> Result<(), String>;
    fn as_any(&self) -> &dyn Any;
}

macro_rules! downcast_dbg {
    ($formatter:expr, $module:expr, $($type:ty),+) => {
        {
            let mut result: Option<std::fmt::Result> = None;
            $(
                if result.is_none() {
                    if let Some(debuggable) = $module.as_any().downcast_ref::<$type>() {
                        result = Some(write!($formatter, "{:?}", debuggable));
                    }
                }
            )+
            result.unwrap_or_else(|| write!($formatter, "Unknown ModuleMap implementation"))
        }
    };
}

impl std::fmt::Debug for dyn ModuleMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        downcast_dbg!(
            f,
            self,
            IntronRetentionDescriptor,
            StartTruncationDescriptor,
            FusionDetectionDescriptor
        )
    }
}

#[allow(dead_code)]
#[derive(Debug)]
pub struct ModuleDescriptor {
    module: ModuleType,
}

impl ModuleDescriptor {
    pub fn with_schema(module: ModuleType) -> Box<dyn ModuleMap> {
        match module {
            ModuleType::IntronRetention => IntronRetentionDescriptor::new(),
            ModuleType::StartTruncation => StartTruncationDescriptor::new(),
            ModuleType::FusionDetection => FusionDetectionDescriptor::new(),
        }
    }
}

pub struct IntronRetentionDescriptor {
    pub intron_retention: Value,
    pub is_retention_supported: Value,
    pub is_retention_supported_map: Value,
    pub component_size: Value,
    pub ref_component_size: Value,
    pub query_component_size: Value,
    pub component_retention_ratio: Value,
    pub is_dirty_component: Value,
    pub intron_support_ratio: Value,
    pub exon_support_ratio: Value,
    pub number_of_retentions: Value,
    pub number_of_true_retentions: Value,
    pub number_of_partial_retentions: Value,
    pub number_of_false_retentions: Value,
    pub number_of_recovers: Value,
    pub number_of_unrecovers: Value,
    pub number_of_true_retentions_supported: Value,
    pub number_of_partial_retentions_supported: Value,
    pub number_of_false_retentions_supported: Value,
    pub location_of_retention: Value,
    pub retention_acceptor_score: Value,
    pub retention_donor_score: Value,
    pub retention_in_cds: Value,
    pub retention_in_utr: Value,
    pub is_intron_retained_in_frame: Value,
    pub is_toga_intron: Value,
}

impl IntronRetentionDescriptor {
    pub fn new() -> Box<Self> {
        Box::new(Self {
            intron_retention: Value::Bool(false),
            is_retention_supported: Value::Null,
            is_retention_supported_map: Value::Null,
            component_size: Value::Null,
            ref_component_size: Value::Null,
            query_component_size: Value::Null,
            component_retention_ratio: Value::Null,
            is_dirty_component: Value::Bool(false),
            intron_support_ratio: Value::Null,
            exon_support_ratio: Value::Null,
            number_of_retentions: Value::Number(0.into()),
            number_of_true_retentions: Value::Null,
            number_of_partial_retentions: Value::Null,
            number_of_false_retentions: Value::Null,
            number_of_recovers: Value::Number(0.into()),
            number_of_unrecovers: Value::Number(0.into()),
            number_of_true_retentions_supported: Value::Null,
            number_of_partial_retentions_supported: Value::Null,
            number_of_false_retentions_supported: Value::Null,
            location_of_retention: Value::Null,
            retention_acceptor_score: Value::Null,
            retention_donor_score: Value::Null,
            retention_in_cds: Value::Null,
            retention_in_utr: Value::Null,
            is_intron_retained_in_frame: Value::Null,
            is_toga_intron: Value::Null,
        })
    }
}

#[derive(Debug, Clone)]
pub enum IntronRetentionValue {
    IsIntronRetention,
    IsRetentionSupported,
    IsRetentionSupportedMap,
    ComponentSize,
    RefComponentSize,
    QueryComponentSize,
    ComponentRetentionRatio,
    IsDirtyComponent,
    IntronSupportRatio,
    ExonSupportRatio,
    NumberOfRetentions,
    NumberOfTrueRetentions,
    NumberOfPartialRetentions,
    NumberOfFalseRetentions,
    NumberOfRecovers,
    NumberOfUnrecovers,
    NumberOfTrueRententionsSupported,
    NumberOfPartialRententionsSupported,
    NumberOfFalseRententionsSupported,
    RetentionLocation,
    RetentionAcceptorScore,
    RetentionDonorScore,
    IsRetentionInCds,
    IsRetentionInUtr,
    IsIntronRetainedInFrame,
    IsTogaIntron,
}

impl ModuleMap for IntronRetentionDescriptor {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<IntronRetentionValue>() {
            match *key {
                IntronRetentionValue::IsIntronRetention => Some(self.intron_retention.clone()),
                IntronRetentionValue::IsRetentionSupported => {
                    Some(self.is_retention_supported.clone())
                }
                IntronRetentionValue::IsRetentionSupportedMap => {
                    Some(self.is_retention_supported_map.clone())
                }
                IntronRetentionValue::ComponentSize => Some(self.component_size.clone()),
                IntronRetentionValue::RefComponentSize => Some(self.ref_component_size.clone()),
                IntronRetentionValue::QueryComponentSize => Some(self.query_component_size.clone()),
                IntronRetentionValue::ComponentRetentionRatio => {
                    Some(self.intron_support_ratio.clone())
                }
                IntronRetentionValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
                IntronRetentionValue::IntronSupportRatio => Some(self.intron_support_ratio.clone()),
                IntronRetentionValue::ExonSupportRatio => Some(self.exon_support_ratio.clone()),
                IntronRetentionValue::NumberOfRetentions => Some(self.number_of_retentions.clone()),
                IntronRetentionValue::NumberOfTrueRetentions => {
                    Some(self.number_of_true_retentions.clone())
                }
                IntronRetentionValue::NumberOfPartialRetentions => {
                    Some(self.number_of_partial_retentions.clone())
                }
                IntronRetentionValue::NumberOfFalseRetentions => {
                    Some(self.number_of_false_retentions.clone())
                }
                IntronRetentionValue::NumberOfRecovers => Some(self.number_of_recovers.clone()),
                IntronRetentionValue::NumberOfUnrecovers => Some(self.number_of_unrecovers.clone()),
                IntronRetentionValue::NumberOfTrueRententionsSupported => {
                    Some(self.number_of_true_retentions_supported.clone())
                }
                IntronRetentionValue::NumberOfPartialRententionsSupported => {
                    Some(self.number_of_partial_retentions_supported.clone())
                }
                IntronRetentionValue::NumberOfFalseRententionsSupported => {
                    Some(self.number_of_false_retentions_supported.clone())
                }
                IntronRetentionValue::RetentionLocation => Some(self.location_of_retention.clone()),
                IntronRetentionValue::RetentionAcceptorScore => {
                    Some(self.retention_acceptor_score.clone())
                }
                IntronRetentionValue::RetentionDonorScore => {
                    Some(self.retention_donor_score.clone())
                }
                IntronRetentionValue::IsRetentionInCds => Some(self.retention_in_cds.clone()),
                IntronRetentionValue::IsRetentionInUtr => Some(self.retention_in_utr.clone()),
                IntronRetentionValue::IsIntronRetainedInFrame => {
                    Some(self.retention_in_utr.clone())
                }
                IntronRetentionValue::IsTogaIntron => Some(self.is_toga_intron.clone()),
            }
        } else {
            None
        }
    }

    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<IntronRetentionValue>() {
            match *key {
                IntronRetentionValue::IsIntronRetention => {
                    self.intron_retention = value;
                    Ok(())
                }
                IntronRetentionValue::IsRetentionSupported => {
                    self.is_retention_supported = value;
                    Ok(())
                }
                IntronRetentionValue::IsRetentionSupportedMap => {
                    self.is_retention_supported_map = value;
                    Ok(())
                }
                IntronRetentionValue::ComponentSize => {
                    self.component_size = value;
                    Ok(())
                }
                IntronRetentionValue::RefComponentSize => {
                    self.ref_component_size = value;
                    Ok(())
                }
                IntronRetentionValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                IntronRetentionValue::ComponentRetentionRatio => {
                    self.component_retention_ratio = value;
                    Ok(())
                }
                IntronRetentionValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
                IntronRetentionValue::IntronSupportRatio => {
                    self.intron_support_ratio = value;
                    Ok(())
                }
                IntronRetentionValue::ExonSupportRatio => {
                    self.exon_support_ratio = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfRetentions => {
                    self.number_of_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfTrueRetentions => {
                    self.number_of_true_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfPartialRetentions => {
                    self.number_of_partial_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfFalseRetentions => {
                    self.number_of_false_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfRecovers => {
                    self.number_of_recovers = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfUnrecovers => {
                    self.number_of_unrecovers = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfTrueRententionsSupported => {
                    self.number_of_true_retentions_supported = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfPartialRententionsSupported => {
                    self.number_of_partial_retentions_supported = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfFalseRententionsSupported => {
                    self.number_of_false_retentions_supported = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionLocation => {
                    self.location_of_retention = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionAcceptorScore => {
                    self.retention_acceptor_score = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionDonorScore => {
                    self.retention_donor_score = value;
                    Ok(())
                }
                IntronRetentionValue::IsRetentionInCds => {
                    self.retention_in_cds = value;
                    Ok(())
                }
                IntronRetentionValue::IsRetentionInUtr => {
                    self.retention_in_utr = value;
                    Ok(())
                }
                IntronRetentionValue::IsIntronRetainedInFrame => {
                    self.is_intron_retained_in_frame = value;
                    Ok(())
                }
                IntronRetentionValue::IsTogaIntron => {
                    self.is_toga_intron = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl std::fmt::Debug for IntronRetentionDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_intron_retention: {:?},
            is_retention_supported: {:?},
            is_retention_supported_map: {:?},
            component_size: {:?},
            ref_component_size: {:?},
            query_component_size: {:?},
            component_retention_ratio: {:?},
            is_dirty_component: {:?},
            intron_support_ratio: {:?},
            exon_support_ratio: {:?},
            number_of_retentions: {:?},
            number_of_true_retentions: {:?},
            number_of_partial_retentions: {:?},
            number_of_false_retentions: {:?},
            number_of_recovers: {:?},
            number_of_unrecovers: {:?},
            number_of_true_retentions_supported: {:?},
            number_of_partial_retentions_supported: {:?},
            number_of_false_retentions_supported: {:?},
            retention_location: {:?},
            retention_acceptor_score: {:?},
            retention_donor_score: {:?},
            is_retention_in_cds: {:?},
            is_retention_in_utr: {:?},
            is_intron_retained_in_frame: {:?},
            is_toga_intron: {:?}
            }}",
            self.intron_retention,
            self.is_retention_supported,
            self.is_retention_supported_map,
            self.component_size,
            self.ref_component_size,
            self.query_component_size,
            self.component_retention_ratio,
            self.is_dirty_component,
            self.intron_support_ratio,
            self.exon_support_ratio,
            self.number_of_retentions,
            self.number_of_true_retentions,
            self.number_of_partial_retentions,
            self.number_of_false_retentions,
            self.number_of_recovers,
            self.number_of_unrecovers,
            self.number_of_true_retentions_supported,
            self.number_of_partial_retentions_supported,
            self.number_of_false_retentions_supported,
            self.location_of_retention,
            self.retention_acceptor_score,
            self.retention_donor_score,
            self.retention_in_cds,
            self.retention_in_utr,
            self.is_intron_retained_in_frame,
            self.is_toga_intron
        )
    }
}

pub struct StartTruncationDescriptor {
    pub is_read_truncated: Value,
    pub is_novel_start: Value,
    pub is_dirty_component: Value,
    pub component_size: Value,
    pub ref_component_size: Value,
    pub query_component_size: Value,
    pub truncation_support_ratio: Value,
    pub is_truncation_supported: Value,
    pub component_truncation_ratio: Value,
}

impl StartTruncationDescriptor {
    pub fn new() -> Box<Self> {
        Box::new(Self {
            is_read_truncated: Value::Null,
            is_novel_start: Value::Null,
            is_dirty_component: Value::Null,
            component_size: Value::Null,
            ref_component_size: Value::Null,
            query_component_size: Value::Null,
            truncation_support_ratio: Value::Null,
            is_truncation_supported: Value::Null,
            component_truncation_ratio: Value::Null,
        })
    }
}

#[derive(Debug, Clone)]
pub enum StartTruncationValue {
    IsReadTruncated,
    IsNovelStart,
    TruncationSupportRatio,
    IsTruncationSupported,
    ComponentSize,
    RefComponentSize,
    QueryComponentSize,
    ComponentTruncationRatio,
    IsDirtyComponent,
}

impl ModuleMap for StartTruncationDescriptor {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<StartTruncationValue>() {
            match *key {
                StartTruncationValue::IsReadTruncated => Some(self.is_read_truncated.clone()),
                StartTruncationValue::IsNovelStart => Some(self.is_novel_start.clone()),
                StartTruncationValue::TruncationSupportRatio => {
                    Some(self.truncation_support_ratio.clone())
                }
                StartTruncationValue::IsTruncationSupported => {
                    Some(self.is_truncation_supported.clone())
                }
                StartTruncationValue::ComponentSize => Some(self.component_size.clone()),
                StartTruncationValue::RefComponentSize => Some(self.ref_component_size.clone()),
                StartTruncationValue::QueryComponentSize => Some(self.query_component_size.clone()),
                StartTruncationValue::ComponentTruncationRatio => {
                    Some(self.component_truncation_ratio.clone())
                }
                StartTruncationValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
            }
        } else {
            None
        }
    }

    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<StartTruncationValue>() {
            match *key {
                StartTruncationValue::IsReadTruncated => {
                    self.is_read_truncated = value;
                    Ok(())
                }
                StartTruncationValue::IsNovelStart => {
                    self.is_novel_start = value;
                    Ok(())
                }
                StartTruncationValue::TruncationSupportRatio => {
                    self.truncation_support_ratio = value;
                    Ok(())
                }
                StartTruncationValue::IsTruncationSupported => {
                    self.is_truncation_supported = value;
                    Ok(())
                }
                StartTruncationValue::ComponentSize => {
                    self.component_size = value;
                    Ok(())
                }
                StartTruncationValue::RefComponentSize => {
                    self.ref_component_size = value;
                    Ok(())
                }
                StartTruncationValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                StartTruncationValue::ComponentTruncationRatio => {
                    self.component_truncation_ratio = value;
                    Ok(())
                }
                StartTruncationValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl std::fmt::Debug for StartTruncationDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_read_truncated: {:?},
            is_novel_start: {:?},
            is_dirty_component: {:?},
            component_size: {:?},
            ref_component_size: {:?},
            query_component_size: {:?},
            truncation_support_ratio: {:?},
            is_truncation_supported: {:?},
            component_truncation_ratio: {:?}
            }}",
            self.is_read_truncated,
            self.is_novel_start,
            self.is_dirty_component,
            self.component_size,
            self.ref_component_size,
            self.query_component_size,
            self.truncation_support_ratio,
            self.is_truncation_supported,
            self.component_truncation_ratio
        )
    }
}

pub struct FusionDetectionDescriptor {
    is_fused_read: Value,
    is_fusion_supported: Value,
    component_size: Value,
    ref_component_size: Value,
    query_component_size: Value,
    component_fusion_ratio: Value,
    is_dirty_component: Value,
    location_of_fusion: Value,
    fusion_in_frame: Value,
}

impl FusionDetectionDescriptor {
    pub fn new() -> Box<Self> {
        Box::new(Self {
            is_fused_read: Value::Null,
            is_fusion_supported: Value::Null,
            component_size: Value::Null,
            ref_component_size: Value::Null,
            query_component_size: Value::Null,
            component_fusion_ratio: Value::Null,
            is_dirty_component: Value::Null,
            location_of_fusion: Value::Null,
            fusion_in_frame: Value::Null,
        })
    }
}

pub enum FusionDetectionValue {
    IsFusedRead,
    IsFusionSupported,
    ComponentSize,
    RefComponentSize,
    QueryComponentSize,
    ComponentFusionRatio,
    IsDirtyComponent,
    LocationOfFusion,
    FusionInFrame,
}

impl ModuleMap for FusionDetectionDescriptor {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<FusionDetectionValue>() {
            match *key {
                FusionDetectionValue::IsFusedRead => Some(self.is_fused_read.clone()),
                FusionDetectionValue::IsFusionSupported => Some(self.is_fusion_supported.clone()),
                FusionDetectionValue::ComponentSize => Some(self.component_size.clone()),
                FusionDetectionValue::RefComponentSize => Some(self.ref_component_size.clone()),
                FusionDetectionValue::QueryComponentSize => Some(self.query_component_size.clone()),
                FusionDetectionValue::ComponentFusionRatio => {
                    Some(self.component_fusion_ratio.clone())
                }
                FusionDetectionValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
                FusionDetectionValue::LocationOfFusion => Some(self.location_of_fusion.clone()),
                FusionDetectionValue::FusionInFrame => Some(self.fusion_in_frame.clone()),
            }
        } else {
            None
        }
    }

    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<FusionDetectionValue>() {
            match *key {
                FusionDetectionValue::IsFusedRead => {
                    self.is_fused_read = value;
                    Ok(())
                }
                FusionDetectionValue::IsFusionSupported => {
                    self.is_fusion_supported = value;
                    Ok(())
                }
                FusionDetectionValue::ComponentSize => {
                    self.component_size = value;
                    Ok(())
                }
                FusionDetectionValue::RefComponentSize => {
                    self.ref_component_size = value;
                    Ok(())
                }
                FusionDetectionValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                FusionDetectionValue::ComponentFusionRatio => {
                    self.component_fusion_ratio = value;
                    Ok(())
                }
                FusionDetectionValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
                FusionDetectionValue::LocationOfFusion => {
                    self.location_of_fusion = value;
                    Ok(())
                }
                FusionDetectionValue::FusionInFrame => {
                    self.fusion_in_frame = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl std::fmt::Debug for FusionDetectionDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_fused_read: {:?},
            is_fusion_supported: {:?},
            component_size: {:?},
            ref_component_size: {:?},
            query_component_size: {:?}
            component_fusion_ratio: {:?},
            is_dirty_component: {:?},
            location_of_fusion: {:?},
            fusion_in_frame: {:?}
            }}",
            self.is_fused_read,
            self.is_fusion_supported,
            self.component_size,
            self.ref_component_size,
            self.query_component_size,
            self.component_fusion_ratio,
            self.is_dirty_component,
            self.location_of_fusion,
            self.fusion_in_frame,
        )
    }
}

// quality of life improvement fns
#[inline(always)]
pub fn exonic_overlap<N, I>(exons_a: &I, exons_b: &I) -> bool
where
    N: Num + NumCast + Copy + PartialOrd,
    I: IntoIterator,
    I::Item: Borrow<(N, N)>,
    for<'a> &'a I: IntoIterator<Item = &'a I::Item>,
{
    let mut iter_a = exons_a.into_iter();
    let mut iter_b = exons_b.into_iter();

    let mut exon_a = iter_a.next();
    let mut exon_b = iter_b.next();

    loop {
        match (exon_a, exon_b) {
            (Some(start_end_a), Some(start_end_b)) => {
                let (start_a, end_a) = start_end_a.borrow();
                let (start_b, end_b) = start_end_b.borrow();

                if *start_a < *end_b && *start_b < *end_a {
                    return true;
                }

                if *end_a < *end_b {
                    exon_a = iter_a.next();
                } else {
                    exon_b = iter_b.next();
                }
            }
            _ => break,
        }
    }

    false
}
