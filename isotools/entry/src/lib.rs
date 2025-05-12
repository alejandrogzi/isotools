use config::{merge, write_descriptor};

use iso_classify::lib_iso_classify;
use iso_intron::lib_iso_intron;
use iso_polya::lib_iso_polya;
use iso_utr::lib_iso_utr;

use std::sync::Arc;

const KEYS: [&str; 5] = ["--query", "--toga", "--aparent", "--bigwig", "--twobit"];
const GLOBAL_DESCRIPTOR: &str = "global_descriptor.tsv";

pub fn lib(mut args: Vec<String>) {
    __check_args(&args);

    // WARN: will expect to always have outdir as last argument
    let outdir = args.pop().expect(&format!(
        "ERROR: Missing output directory argument, you had: {:?}",
        args
    ));

    let introns = lib_iso_classify(args.clone())
        .expect("ERROR: Failed to classify introns")
        .display()
        .to_string();

    args.extend(vec![
        "--introns".to_string(),
        introns,
        "--in-memory".to_string(),
    ]);

    let args = Arc::new(args);

    // TODO: -F/--fusions -> will include lib_iso_fusion(args.clone());
    let retentions = lib_iso_intron(args.clone());
    let truncations = lib_iso_utr(args.clone());
    let intraprimings = lib_iso_polya(args);

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
