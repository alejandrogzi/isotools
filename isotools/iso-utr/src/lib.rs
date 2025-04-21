pub mod cli;
pub mod core;
pub mod utils;

use config::ModuleMap;
use dashmap::DashMap;
use std::sync::Arc;

pub fn lib_iso_utr(args: Arc<Vec<String>>) -> DashMap<String, Box<dyn ModuleMap>> {
    let args = cli::Args::from(args);
    let descriptor =
        crate::core::detect_truncations(args).expect("ERROR: Failed to detect truncations");

    return descriptor;
}
