use dashmap::DashSet;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Duration;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const INTRON: &str = "intron";
pub const FIVEND: &str = "fivend";
pub const THREEND: &str = "threend";

// numeric values
pub const SCALE: u64 = 100000000000; // 100Gb
pub const MIN_THREADS: usize = 1;
pub const MIN_BED_FIELDS: usize = 12;
pub const MIN_BED4_FIELDS: usize = 4;

// file names
pub const HIT: &str = "hits.bed";
pub const PASS: &str = "pass.bed";
pub const BED3: &str = "ir.bed";
pub const F5: &str = "5ends.txt";

// flags
pub const COLORIZE: bool = false;
pub const OVERLAP: bool = true;

#[cfg(not(windows))]
const TICK_SETTINGS: (&str, u64) = ("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏ ", 80);

#[cfg(windows)]
const TICK_SETTINGS: (&str, u64) = (r"+-x| ", 200);

/// Return a pre-configured progress bar
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

pub fn write_objs<T>(data: &DashSet<T>, fname: &str)
where
    T: AsRef<str> + Sync + Send + Eq + std::hash::Hash,
{
    log::info!("Reads in {}: {:?}", fname, data.len());
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
