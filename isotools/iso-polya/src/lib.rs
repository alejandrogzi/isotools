pub mod cli;
pub mod core;
pub mod utils;

use config::ModuleMap;
use dashmap::DashMap;
use std::sync::Arc;

pub fn lib_iso_polya(args: Arc<Vec<String>>) -> DashMap<String, Box<dyn ModuleMap>> {
    let args = cli::CallerArgs::from(args);
    let descriptor =
        crate::core::pas::pas_caller(args).expect("ERROR: Failed to detect intron retentions");

    return descriptor;
}
