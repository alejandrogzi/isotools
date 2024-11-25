use dashmap::DashSet;
use indicatif::{ProgressBar, ProgressStyle};
use serde_json::Value;
use std::any::Any;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Duration;
use thiserror::Error;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// numeric values
pub const SCALE: u64 = 100000000000; // 100Gb
pub const MIN_THREADS: usize = 1;
pub const MIN_BED_FIELDS: usize = 12;
pub const MIN_BED4_FIELDS: usize = 4;
pub const TRUNCATION_THRESHOLD: f32 = 0.5;
pub const TRUNCATION_RECOVERY_THRESHOLD: f32 = 0.5;
pub const RETENTION_RATIO_THRESHOLD: f32 = 0.5;
pub const EXON_RETENTION_RECOVERY_THRESHOLD: f32 = 0.5;
pub const INTRON_RETENTION_RECOVERY_THRESHOLD: f32 = 0.5;

// file names
pub const HIT: &str = "hits.bed";
pub const PASS: &str = "pass.bed";
pub const BED3: &str = "ir.bed";
pub const F5: &str = "5ends.txt";
pub const INTERGENIC_REGIONS: &str = "intergenic.bed";
pub const CHIMERAS: &str = "chimeras.bed";
pub const CHIMERIC_FREE: &str = "chimeric.free.bed";

// flags
pub const COLORIZE: bool = false;
pub const OVERLAP_CDS: bool = false;
pub const OVERLAP_EXON: bool = true;

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
        return Err(CliError::InvalidInput(format!("{:?} does not exist", arg)));
    }

    if !arg.is_file() {
        return Err(CliError::InvalidInput(format!("{:?} is not a file", arg)));
    }

    match arg.extension() {
        Some(ext) if ext == "bed" => (),
        _ => {
            return Err(CliError::InvalidInput(format!(
                "file {:?} is not a BED file",
                arg
            )))
        }
    }

    match std::fs::metadata(arg) {
        Ok(metadata) if metadata.len() == 0 => {
            Err(CliError::InvalidInput(format!("file {:?} is empty", arg)))
        }
        Ok(_) => Ok(()),
        Err(e) => Err(CliError::IoError(e)),
    }
}

// module descriptors
#[derive(Debug)]
pub enum ModuleType {
    IntronRetention,
    StartTruncation,
    FusionRead,
}

pub trait ModuleMap: Any {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value>;
    fn set_value(&mut self, key: Box<dyn Any>, value: serde_json::Value) -> Result<(), String>;
    fn as_any(&self) -> &dyn Any;
}

impl std::fmt::Debug for dyn ModuleMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(debuggable) = self.as_any().downcast_ref::<IntronRetentionDescriptor>() {
            write!(f, "{:?}", debuggable)
        } else {
            write!(f, "Unknown ModuleMap implementation")
        }
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
            ModuleType::StartTruncation => unimplemented!(),
            ModuleType::FusionRead => unimplemented!(),
        }
    }
}

pub struct IntronRetentionDescriptor {
    pub intron_retention: Value,
    pub is_retention_supported: Value,
    pub is_retention_supported_map: Value,
    pub retention_support_ratio: Value,
    pub exon_support_ratio: Value,
    pub number_of_retentions: Value,
    pub number_of_recovers: Value,
    pub number_of_unrecovers: Value,
    pub location_of_retention: Value,
    pub retention_in_cds: Value,
    pub retention_in_utr: Value,
}

impl IntronRetentionDescriptor {
    pub fn new() -> Box<Self> {
        Box::new(Self {
            intron_retention: Value::Bool(false),
            is_retention_supported: Value::Bool(false),
            is_retention_supported_map: Value::Array(vec![]),
            retention_support_ratio: Value::Array(vec![]),
            exon_support_ratio: Value::Array(vec![]),
            number_of_retentions: Value::Number(0.into()),
            number_of_recovers: Value::Number(0.into()),
            number_of_unrecovers: Value::Number(0.into()),
            location_of_retention: Value::Array(vec![]),
            retention_in_cds: Value::Array(vec![]),
            retention_in_utr: Value::Array(vec![]),
        })
    }
}

#[derive(Debug, Clone)]
pub enum IntronRetentionValue {
    IsIntronRetention,
    IsRetentionSupported,
    IsRetentionSupportedMap,
    RetentionSupportRatio,
    ExonSupportRatio,
    NumberOfRetentions,
    NumberOfRecovers,
    NumberOfUnrecovers,
    RetentionLocation,
    RetentionInCds,
    RetentionInUtr,
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
                IntronRetentionValue::RetentionSupportRatio => {
                    Some(self.retention_support_ratio.clone())
                }
                IntronRetentionValue::ExonSupportRatio => Some(self.exon_support_ratio.clone()),
                IntronRetentionValue::NumberOfRetentions => Some(self.number_of_retentions.clone()),
                IntronRetentionValue::NumberOfRecovers => Some(self.number_of_recovers.clone()),
                IntronRetentionValue::NumberOfUnrecovers => Some(self.number_of_unrecovers.clone()),
                IntronRetentionValue::RetentionLocation => Some(self.location_of_retention.clone()),
                IntronRetentionValue::RetentionInCds => Some(self.retention_in_cds.clone()),
                IntronRetentionValue::RetentionInUtr => Some(self.retention_in_utr.clone()),
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
                IntronRetentionValue::RetentionSupportRatio => {
                    self.retention_support_ratio = value;
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
                IntronRetentionValue::NumberOfRecovers => {
                    self.number_of_recovers = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfUnrecovers => {
                    self.number_of_unrecovers = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionLocation => {
                    self.location_of_retention = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionInCds => {
                    self.retention_in_cds = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionInUtr => {
                    self.retention_in_utr = value;
                    Ok(())
                }
            }
        } else {
            Err("Invalid key".to_string())
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
            intron_retention: {:?},
            is_retention_supported: {:?},
            is_retention_supported_map: {:?},
            retention_support_ratio: {:?},
            exon_support_ratio: {:?},
            number_of_retentions: {:?},
            number_of_recovers: {:?},
            number_of_unrecovers: {:?},
            location_of_retention: {:?},
            retention_in_cds: {:?},
            retention_in_utr: {:?}
            }}",
            self.intron_retention,
            self.is_retention_supported,
            self.is_retention_supported_map,
            self.retention_support_ratio,
            self.exon_support_ratio,
            self.number_of_retentions,
            self.number_of_recovers,
            self.number_of_unrecovers,
            self.location_of_retention,
            self.retention_in_cds,
            self.retention_in_utr
        )
    }
}
