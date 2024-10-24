use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;

pub const MIN_THREADS: usize = 1;
pub const MIN_BED_FIELDS: usize = 12;
pub const MIN_BED4_FIELDS: usize = 4;
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const INTRON: &str = "intron";
pub const FIVEND: &str = "fivend";
pub const THREEND: &str = "threend";
pub const SCALE: u64 = 100000000000; // 100Gb
pub const HIT: &str = "hits.bed";
pub const PASS: &str = "pass.bed";
pub const BED3: &str = "ir.bed";
pub const F5: &str = "5ends.txt";

pub const COLORIZE: bool = true;
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
