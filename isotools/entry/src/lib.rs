use config::{
    bed_to_struct_collection, merge, write_descriptor, IntronRetentionValue, ModuleMap,
    PolyAPredictionValue, StartTruncationValue,
};
use dashmap::DashSet;
use hashbrown::HashMap;
use packbed::GenePred;
use rayon::prelude::*;

use iso_classify::lib_iso_classify;
use iso_intron::lib_iso_intron;
use iso_polya::lib_iso_polya;
use iso_utr::lib_iso_utr;

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Arc,
};

const KEYS: [&str; 6] = [
    "--query",
    "--toga",
    "--aparent",
    "--bigwig",
    "--twobit",
    "--outdir",
];
const GLOBAL_DESCRIPTOR: &str = "global_descriptor.tsv";

pub fn lib(mut args: Vec<String>) {
    __check_args(&args);

    let introns = lib_iso_classify(args.clone())
        .expect("ERROR: Failed to classify introns")
        .display()
        .to_string();

    // WARN: will expect to always have outdir as last argument [ last 2 ]
    // INFO: fmt -> '[--arg1, <VALUE>, --arg2, <VALUE>]'
    let outdir = args.pop().expect(&format!(
        "ERROR: Missing output directory argument, you had: {:?}",
        args
    ));
    let mut args = args[..args.len() - 1].to_vec(); // INFO: dropping --outdir

    args.extend(vec![
        "--introns".to_string(),
        introns,
        "--in-memory".to_string(),
    ]);

    let args = Arc::new(args);

    let retentions = lib_iso_intron(args.clone());
    let truncations = lib_iso_utr(args.clone());
    let intraprimings = lib_iso_polya(args.clone());

    let global_descriptor = merge(vec![
        retentions,
        truncations,
        dashmap::DashMap::new(),
        intraprimings,
    ]);

    write_descriptor(
        &global_descriptor,
        format!("{}/{}", outdir, GLOBAL_DESCRIPTOR).as_str(),
    );

    // descriptor_to_bed(global_descriptor, outdir, args[1].clone());
}

/// Check if all required arguments are present
///
/// # Arguments
///
/// * `args` - A vector of strings representing the command line arguments
///
/// # Returns
///
/// None
///
/// # Example
///
/// ```rust, no_run
/// use iso_classify::lib;
///
/// let args = vec![
///    "--query".to_string(),
///   "--toga".to_string(),
///   "--aparent".to_string(),
///   "--bigwig".to_string(),
///   "--twobit".to_string(),
/// ];
///
/// lib::lib(args);
/// ```
fn __check_args(args: &Vec<String>) {
    for key in KEYS.iter() {
        if !args.contains(&key.to_string()) {
            log::error!("Missing required argument: {}", key);
            std::process::exit(1);
        }
    }
}

pub fn descriptor_to_bed(
    descriptor: dashmap::DashMap<String, Box<dyn ModuleMap>>,
    outdir: String,
    bed: String,
) {
    let outdir = PathBuf::from(outdir);

    let clean = DashSet::new();
    let intrapriming = DashSet::new();
    let rts = DashSet::new();
    let retentions = DashSet::new();
    let trash = DashSet::new();

    descriptor.par_iter().for_each(|record| {
        let key = record.key().clone();
        let value = record.value();

        let mut flaw_count = 0;

        let is_intrapriming = value
            .get_value(Box::new(PolyAPredictionValue::IsIntrapriming))
            .unwrap_or_else(|| {
                log::warn!("Missing value for IsIntrapriming in key: {}", key);
                std::process::exit(1);
            });
        let is_truncated = value
            .get_value(Box::new(StartTruncationValue::IsReadTruncated))
            .unwrap_or_else(|| {
                log::warn!("Missing value for IsReadTruncated in key: {}", key);
                std::process::exit(1);
            });
        let exonic_status = value
            .get_value(Box::new(IntronRetentionValue::ExonicStatus))
            .unwrap_or_else(|| {
                log::warn!("Missing value for ExonicStatus in key: {}", key);
                std::process::exit(1);
            });
        let intronic_status = value
            .get_value(Box::new(IntronRetentionValue::IntronicStatus))
            .unwrap_or_else(|| {
                log::warn!("Missing value for IntronicStatus in key: {}", key);
                std::process::exit(1);
            });
        let has_rt_intron = value
            .get_value(Box::new(IntronRetentionValue::HasRTIntron))
            .unwrap_or_else(|| {
                log::warn!("Missing value for HasRTIntron in key: {}", key);
                std::process::exit(1);
            });

        if is_intrapriming == serde_json::Value::Bool(true) {
            flaw_count += 1;
        }

        if is_truncated == serde_json::Value::Bool(true) {
            flaw_count += 1;
        }

        if has_rt_intron == serde_json::Value::Bool(true) {
            flaw_count += 1;
        }

        if exonic_status == serde_json::Value::String("DISCARD".to_string())
            && intronic_status == serde_json::Value::String("DISCARD".to_string())
        {
            flaw_count += 1;
        }

        if flaw_count >= 2 {
            trash.insert(key);
        } else if has_rt_intron == serde_json::Value::Bool(true) {
            rts.insert(key);
        } else if is_intrapriming == serde_json::Value::Bool(true) {
            intrapriming.insert(key);
        } else if exonic_status == serde_json::Value::String("DISCARD".to_string())
            && intronic_status == serde_json::Value::String("DISCARD".to_string())
        {
            retentions.insert(key);
        } else {
            if flaw_count == 0 {
                clean.insert(key);
            }
        }
    });

    // INFO: read the bed file to be able to mutate GenePred
    let records = bed_to_struct_collection::<GenePred>(Arc::new(bed), config::BedColumn::Name)
        .unwrap_or_else(|e| {
            log::error!("Failed to read BED file: {}", e);
            std::process::exit(1);
        });

    for (collection, category) in std::iter::zip(
        vec![clean, intrapriming, rts, retentions, trash],
        vec![
            ReadCategory::Clean,
            ReadCategory::Intrapriming,
            ReadCategory::RTintron,
            ReadCategory::Retention,
            ReadCategory::Trash,
        ],
    ) {
        extract_from_bed(collection, outdir.clone(), &records, &descriptor, category);
    }
}

pub enum ReadCategory {
    Clean,
    Intrapriming,
    RTintron,
    Retention,
    Trash,
}

fn extract_from_bed(
    collection: DashSet<String>, // INFO: key is read name!
    outdir: PathBuf,
    records: &dashmap::DashMap<String, HashMap<String, GenePred>>,
    _descriptor: &dashmap::DashMap<String, Box<dyn ModuleMap>>,
    category: ReadCategory,
) {
    let file = match category {
        ReadCategory::Clean => outdir.join("clean.bed"),
        ReadCategory::Intrapriming => outdir.join("intrapriming.bed"),
        ReadCategory::RTintron => outdir.join("rt.bed"),
        ReadCategory::Retention => outdir.join("retention.bed"),
        ReadCategory::Trash => outdir.join("trash.bed"),
    };

    let mut writer = BufWriter::new(File::open(file).unwrap());

    // INFO: getting chr from read name R{}_{chr}__
    for read in collection {
        let chr = read
            .split("__")
            .next()
            .unwrap_or_else(|| panic!("ERROR: could not split by BIG_SEP -> {read:?}"))
            .split("_")
            .last()
            .unwrap_or_else(|| panic!("ERROR: could not split by SMALL_SEP -> {read:?}"));

        let mut rc = records
            .get_mut(chr)
            .unwrap_or_else(|| panic!("ERROR: {chr:?} not in records!"));
        let gp = rc
            .get_mut(&read)
            .unwrap_or_else(|| panic!("ERROR: {read:?} not found in records!"));

        let _ = writeln!(writer, "{}", gp.line);
    }

    // INFO: descriptor is passed to build bigBed
}
