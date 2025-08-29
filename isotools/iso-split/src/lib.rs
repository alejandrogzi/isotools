//! Core module for splitting a .fa/.fq file into chunks
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for splitting .fa/.fq files
//! based on custom requirements in parallel.
//!
//! In short, the module accepts any type of .fa or .fq file
//! and process the reads or sequences inside them in parallel
//! when is possible. Compressed files are also accepted. The
//! user has the ability to specify is the splitting process should
//! be done based on specific chunk sizes or number of files, and
//! the amount of parallelization that should be used in the process.

use std::{
    fs::{create_dir_all, File},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    os::unix::fs::symlink,
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::{Ok, Result};
use config::{ChunkRegion, SplitMode};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use memchr::memchr_iter;
use memmap2::Mmap;
use rayon::prelude::*;

const FA_NEEDLE: u8 = b'>';

pub mod cli;
use cli::Args;

/// Dispatches file processing based on its suffix.
///
/// This macro inspects the file name of the given `Path` and attempts to match
/// its suffix against a predefined set of patterns. For each matched suffix,
/// it executes the corresponding action. If no suffix matches, it logs an error
/// and bails out, indicating an unrecognized file format.
///
/// This provides a convenient way to route different file types to specific
/// processing functions, often based on common bioinformatics file extensions
/// like `.fa.gz` or `.fastq`.
///
/// # Arguments
///
/// * `$file` - An expression that evaluates to a reference to a `Path` or `PathBuf`,
///             representing the input file.
/// * `$suffixes_and_actions` - A block containing `suffix => $action` pairs,
///                             where `suffix` is a literal string to match
///                             against the end of the file name, and `$action`
///                             is a block of code to execute if the suffix matches.
///
/// # Panics
///
/// This macro will `unwrap_or_default()` on `file_name().to_str()`, which
/// means it expects valid UTF-8 file names for matching. Non-UTF-8 file names
/// will result in an empty string for `f`, which might not match any suffix.
///
/// # Errors
///
/// Returns an `anyhow::Error` if no suffix matches the input file's name.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use anyhow::Result;
///
/// // Assume these functions exist for the example
/// fn process_fasta(path: &PathBuf) -> Result<()> { /* ... */ Ok(()) }
/// fn process_fastq(path: &PathBuf) -> Result<()> { /* ... */ Ok(()) }
///
/// let my_file = PathBuf::from("data/sequences.fa.gz");
///
/// dispatch!(&my_file, {
///     "fa.gz" => process_fasta(&my_file)?,
///     "fq.gz" => process_fastq(&my_file)?,
///     "fa" => process_fasta(&my_file)?,
/// });
/// ```
#[macro_export]
macro_rules! dispatch {
    ($file:expr, { $($suffix:literal => $action:expr),* $(,)?}) => {{
        let f = $file.file_name().and_then(|f| f.to_str()).unwrap_or_default();
        let mut matched = false;
        $(
            if f.ends_with($suffix) {
                matched = true;
                $action
            }
        )*
        if !matched {
            dbg!(f);
            anyhow::bail!("ERROR: unrecognized file format: {}", $file.display());
        }
    }};
}

/// Splits large sequencing files (FASTA, FASTQ, gzipped versions) into smaller chunks or files.
///
/// This is the main entry point for the `iso-split` functionality. It parses
/// command-line arguments, determines the input file type based on its extension,
/// and dispatches to the appropriate splitting function (`split_fa`, `split_fa_gz`, `split_fq`).
///
/// The splitting logic depends on the `SplitMode` specified in the `Args`
/// (either by a fixed number of records per chunk/file or by a target number of output files).
///
/// # Arguments
///
/// * `args` - A `Vec<String>` representing the command-line arguments passed to the program.
///
/// # Returns
///
/// * `Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                 argument parsing fails, the file format is unrecognized, or
///                 any splitting operation encounters an error.
///
/// # Errors
///
/// * Returns an error if argument parsing (`cli::Args::from`) fails.
/// * Returns an error if the input file's suffix is not recognized by the `dispatch!` macro.
/// * Returns any error propagated from the called splitting functions (`split_fa`, etc.).
///
/// # Example
///
/// ```rust, no_run
/// fn main() -> Result<()> {
///     lib_iso_split(vec!["--file".to_string(), "input.fasta.gz".to_string(), "--chunk-size".to_string(), "1000".to_string()])
///     Ok(())
/// }
/// ```
pub fn lib_iso_split(args: Vec<String>) -> Result<()> {
    let args = cli::Args::from(args);

    let _ = dispatch!(&args.file, {
        "fa.gz" => split_fa_gz(&args)?,
        "fasta.gz" =>  split_fa_gz(&args)?,
        "fq.gz" =>  split_fq(&args)?,
        "fastq.gz" =>  split_fq(&args)?,
        "fa" =>  split_fa(&args)?,
        "fasta" =>  split_fa(&args)?,
    });

    Ok(())
}

/// Splits a non-gzipped FASTA file into multiple smaller FASTA files.
///
/// This function reads a FASTA file, identifies the start positions of all records
/// using `memchr_iter` to find `FA_NEEDLE` (typically `>`). It then divides the
/// file's content (memory-mapped for efficiency) into chunks based on the
/// specified `SplitMode` (either `ChunkSize` or `NumFiles`). Each chunk is then
/// written to a new output file within the designated output directory, utilizing
/// a Rayon thread pool for parallel processing.
///
/// # Arguments
///
/// * `args` - A reference to an `Args` struct containing the input file path,
///            output directory, number of threads, splitting mode, and an optional suffix.
///
/// # Returns
///
/// * `Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                 any operation (file opening, memory mapping, directory creation,
///                 writing to files, or thread pool building) fails.
///
/// # Errors
///
/// * Returns an error if the input FASTA file cannot be opened or memory-mapped.
/// * Returns an error if no FASTA records are found in the input file.
/// * Returns an error if the output directory cannot be created.
/// * Returns an error if `SplitMode::NumFiles` is 0.
/// * Returns any `std::io::Error` during file writing.
///
/// # Parallelism
///
/// This function uses `rayon` for parallel processing of chunks, improving
/// performance for large files. The number of threads is configured via `args.threads`.
///
/// # Example
///
/// ```rust, no_run
/// use anyhow::Result;
/// use std::path::PathBuf;
/// // Assuming Args and SplitMode are defined as in lib_iso_split example
///
/// fn main() -> Result<()> {
///     let args = cli::Args {
///         file: PathBuf::from("input.fa"),
///         outdir: PathBuf::from("fa_chunks"),
///         threads: 4,
///         suffix: Some("part".to_string()),
///         mode_chunk_size: Some(100), // Split into chunks of 100 records
///         mode_num_files: None,
///     };
///     // split_fa(&args)?;
///     println!("Successfully split FASTA file.");
///     Ok(())
/// }
/// ```
pub fn split_fa(args: &Args) -> Result<()> {
    log::info!("INFO: running in FASTA mode with args: {:?}", &args);

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;

    let file = File::open(&args.file)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data = Arc::new(mmap);

    let header: Vec<usize> = memchr_iter(FA_NEEDLE, data.as_ref()).collect();
    if header.is_empty() {
        anyhow::bail!("ERROR: No FASTA records found");
    }

    let mode = args.mode()?;
    let suffix = &args.suffix.clone().unwrap_or(String::new());

    match mode {
        SplitMode::ChunkSize(nchunks) => {
            let chunk_size = (header.len() + nchunks - 1) / nchunks;
            let mut chunks = Vec::new();
            for i in 0..nchunks {
                let start = *header.get(i * chunk_size).unwrap_or(&data.len());
                let end = *header.get((i + 1) * chunk_size).unwrap_or(&data.len());
                chunks.push(ChunkRegion { start, end });
            }

            create_dir_all(&args.outdir)?;

            chunks
                .into_par_iter()
                .enumerate()
                .try_for_each(|(i, chunk)| -> Result<()> {
                    let output_file =
                        PathBuf::from(&args.outdir).join(format!("tmp_chunk_{:03}_{suffix}.fa", i));
                    let mut writer = BufWriter::new(File::create(output_file)?);
                    writer.write_all(&data[chunk.start..chunk.end])?;
                    writer.flush()?;
                    Ok(())
                })?;

            Ok(())
        }
        SplitMode::NumFiles(files) => {
            if files == 0 {
                anyhow::bail!("ERROR: --files must be greater than 0");
            }

            let records_per_file = (header.len() + files - 1) / files;
            create_dir_all(&args.outdir)?;

            (0..files)
                .map(|i| {
                    let start_idx = i * records_per_file;
                    let end_idx = ((i + 1) * records_per_file).min(header.len());

                    let start = *header.get(start_idx).unwrap_or(&data.len());
                    let end = *header.get(end_idx).unwrap_or(&data.len());

                    ChunkRegion { start, end }
                })
                .enumerate()
                .try_for_each(|(i, chunk)| -> Result<()> {
                    let output_file =
                        PathBuf::from(&args.outdir).join(format!("tmp_part_{:03}_{suffix}.fa", i));
                    let mut writer = BufWriter::new(File::create(output_file)?);
                    writer.write_all(&data[chunk.start..chunk.end])?;
                    writer.flush()?;
                    Ok(())
                })?;

            Ok(())
        }
    }
}

/// Splits a gzipped FASTA file (`.fa.gz` or `.fasta.gz`) into multiple smaller gzipped FASTA files.
///
/// This function first decompresses the entire gzipped input file into memory.
/// It then identifies the start positions of all FASTA records using `memchr_iter`.
/// The decompressed data is divided into chunks based on the `SplitMode` (either
/// `ChunkSize` for records per file or `NumFiles` for a target number of files).
/// Each chunk is then individually compressed and written to a new gzipped output file
/// within the specified output directory, leveraging a global Rayon thread pool for parallel execution.
///
/// A special case is handled: if the total number of records is less than or equal to
/// the `records_per_file` when in `ChunkSize` mode, it creates a symlink to the original file
/// instead of splitting, to avoid unnecessary work.
///
/// # Arguments
///
/// * `args` - A reference to an `Args` struct containing the input file path,
///            output directory, number of threads, splitting mode, and an optional suffix.
///
/// # Returns
///
/// * `Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                 any operation (file opening, decompression, directory creation,
///                 writing to files, or thread pool building) fails.
///
/// # Errors
///
/// * Returns an error if the input gzipped FASTA file cannot be opened or decompressed.
/// * Returns an error if no FASTA records are found in the decompressed data.
/// * Returns an error if the output directory cannot be created.
/// * Returns an error if `SplitMode::NumFiles` is 0.
/// * Returns any `std::io::Error` during file writing or compression.
///
/// # Parallelism
///
/// This function uses `rayon` for parallel processing of chunks, improving
/// performance for large files. It builds a global thread pool with `args.threads`.
///
/// # Example
///
/// ```rust, no_run
/// use anyhow::Result;
/// use std::path::PathBuf;
/// // Assuming Args and SplitMode are defined as in lib_iso_split example
///
/// fn main() -> Result<()> {
///     let args = cli::Args {
///         file: PathBuf::from("input.fa.gz"),
///         outdir: PathBuf::from("fa_gz_chunks"),
///         threads: 4,
///         suffix: Some("part".to_string()),
///         mode_num_files: Some(10), // Split into 10 output files
///         mode_chunk_size: None,
///     };
///     // split_fa_gz(&args)?;
///     println!("Successfully split gzipped FASTA file.");
///     Ok(())
/// }
/// ```
pub fn split_fa_gz(args: &Args) -> Result<()> {
    log::info!("INFO: running in FASTA.GZ mode with args: {:?}", &args);

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?; // ensure global thread pool

    let mut gz = MultiGzDecoder::new(File::open(&args.file)?);
    let mut decompressed = Vec::new();
    gz.read_to_end(&mut decompressed)?;
    let data = Arc::new(decompressed);

    let header: Vec<usize> = memchr_iter(FA_NEEDLE, data.as_ref()).collect();
    log::info!(
        "INFO: Found {} records in {}",
        header.len(),
        &args.file.display()
    );

    if header.is_empty() {
        log::error!("ERROR: No FASTA records found in decompressed data");
        anyhow::bail!("ERROR: No FASTA records found in decompressed data");
    }

    let mode = args.mode()?;
    let suffix = args.suffix.clone().unwrap_or_default();
    create_dir_all(&args.outdir)?;

    let chunks: Vec<_> = match mode {
        SplitMode::ChunkSize(records_per_file) => {
            if header.len() <= records_per_file {
                // INFO: fewer records than chunk size -> just symlink
                let outpath =
                    PathBuf::from(&args.outdir).join(format!("tmp_chunk_000_{}.fa.gz", suffix));
                std::fs::create_dir_all(&args.outdir)?;
                log::warn!(
                            "Only {} records found, less than chunk size {}, creating symlink to original file...",
                            header.len(),
                            records_per_file
                        );
                symlink(&args.file, &outpath)?;
                return Ok(());
            }

            let nchunks = (header.len() + records_per_file - 1) / records_per_file;

            (0..nchunks)
                .map(|i| {
                    let start = *header.get(i * records_per_file).unwrap_or(&data.len());
                    let end = *header
                        .get((i + 1) * records_per_file)
                        .unwrap_or(&data.len());
                    ChunkRegion { start, end }
                })
                .collect()
        }

        SplitMode::NumFiles(files) => {
            if files == 0 {
                anyhow::bail!("ERROR: --files must be greater than 0");
            }

            let records_per_file = (header.len() + files - 1) / files;
            (0..files)
                .map(|i| {
                    let start_idx = i * records_per_file;
                    let end_idx = ((i + 1) * records_per_file).min(header.len());
                    let start = *header.get(start_idx).unwrap_or(&data.len());
                    let end = *header.get(end_idx).unwrap_or(&data.len());
                    ChunkRegion { start, end }
                })
                .collect()
        }
    };

    chunks
        .into_par_iter()
        .enumerate()
        .try_for_each(|(i, chunk)| -> Result<()> {
            let output_file =
                PathBuf::from(&args.outdir).join(format!("tmp_chunk_{:03}_{suffix}.fa.gz", i));
            let file = File::create(output_file)?;
            let writer = BufWriter::new(file);

            // You can customize compression level, time, filename, etc.
            let mut encoder = flate2::write::GzEncoder::new(writer, flate2::Compression::fast());

            encoder.write_all(&data[chunk.start..chunk.end])?;
            encoder.finish()?; // ensures footer is written
            Ok(())
        })?;

    Ok(())
}

/// Splits a gzipped FASTQ file (`.fq.gz` or `.fastq.gz`) into multiple smaller gzipped FASTQ files.
///
/// This function acts as a dispatcher for FASTQ splitting, determining the
/// records-per-file based on the `SplitMode` (`ChunkSize` or `NumFiles`).
/// If splitting by `NumFiles`, it first counts all records in the input file
/// to calculate the appropriate `records_per_file` for even distribution.
/// It then delegates the actual splitting and writing to `split_by_chunk_size`.
///
/// # Arguments
///
/// * `args` - A reference to an `Args` struct containing the input file path,
///            output directory, number of threads (though not directly used here,
///            passed to underlying functions), splitting mode, and an optional suffix.
///
/// # Returns
///
/// * `anyhow::Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                         the splitting mode is invalid, record counting fails, or
///                         `split_by_chunk_size` encounters an error.
///
/// # Errors
///
/// * Returns an error if `args.mode()` returns an error (e.g., no split mode specified).
/// * Returns an error if `count_fastq_records` fails.
/// * Propagates any error from `split_by_chunk_size`.
///
/// # Example
///
/// ```rust, no_run
/// use anyhow::Result;
/// use std::path::PathBuf;
/// // Assuming Args and SplitMode are defined as in lib_iso_split example
///
/// fn main() -> Result<()> {
///     let args = cli::Args {
///         file: PathBuf::from("input.fastq.gz"),
///         outdir: PathBuf::from("fq_chunks"),
///         threads: 4, // Not directly used in split_fq, but passed to helpers
///         suffix: Some("batch".to_string()),
///         mode_chunk_size: None,
///         mode_num_files: Some(5), // Split into 5 output files
///     };
///     // split_fq(&args)?;
///     println!("Successfully split FASTQ file.");
///     Ok(())
/// }
/// ```
pub fn split_fq(args: &Args) -> anyhow::Result<()> {
    log::info!("INFO: running in FASTQ mode with args: {:?}", &args);

    let mode = args.mode()?;
    let suffix = &args.suffix.clone().unwrap_or(String::new());

    match mode {
        SplitMode::ChunkSize(records_per_file) => {
            log::info!("INFO: splitting by chunk size!");
            split_by_chunk_size(&args.file, &args.outdir, records_per_file, &suffix)
        }
        SplitMode::NumFiles(num_files) => {
            log::info!("INFO: splitting by number of files!");
            let total_records = count_fastq_records(&args.file)?;
            let records_per_file = (total_records + num_files - 1) / num_files;
            split_by_chunk_size(&args.file, &args.outdir, records_per_file, &suffix)
        }
    }
}

/// Reads a gzipped FASTQ file and splits it into multiple gzipped FASTQ files
/// based on a specified number of records per output file.
///
/// This function iteratively reads four lines at a time (representing one FASTQ record)
/// from the input gzipped file. It writes these records to an output file.
/// Once the `records_per_file` limit is reached for the current output file,
/// a new gzipped output file is created, and writing continues to the new file.
/// Output files are named `tmp_chunk_{index:03}_{suffix}.fastq.gz`.
///
/// # Arguments
///
/// * `input` - A reference to the path of the input gzipped FASTQ file.
/// * `out_dir` - A reference to the path of the directory where output files will be saved.
/// * `records_per_file` - The maximum number of FASTQ records to write into each output file.
/// * `suffix` - An optional string suffix to include in the names of the output files.
///
/// # Returns
///
/// * `anyhow::Result<()>` - An `Ok(())` on successful completion, or an `anyhow::Error` if
///                         any I/O operation (file opening, reading, writing, directory creation) fails,
///                         or if the FASTQ file format is malformed (e.g., unexpected EOF).
///
/// # Errors
///
/// * Returns an error if the input file cannot be opened or read.
/// * Returns an error if the output directory cannot be created.
/// * Returns an `anyhow::Error` if an unexpected end-of-file is encountered, indicating a malformed FASTQ record.
/// * Returns any `std::io::Error` during file writing or compression.
///
/// # Example
///
/// ```rust, no_run
/// use anyhow::Result;
/// use std::path::PathBuf;
///
/// fn main() -> Result<()> {
///     // Assuming an existing 'input.fastq.gz' and an 'output' directory
///     let input_file = PathBuf::from("input.fastq.gz");
///     let output_directory = PathBuf::from("output_chunks");
///     let records_per_chunk = 1000;
///     let file_suffix = "batch1".to_string();
///
///     // Make sure you create a dummy input.fastq.gz for a runnable example
///     // split_by_chunk_size(&input_file, &output_directory, records_per_chunk, &file_suffix)?;
///     println!("FASTQ file split by chunk size successfully.");
///     Ok(())
/// }
/// ```
fn split_by_chunk_size<P: AsRef<Path>>(
    input: P,
    out_dir: P,
    records_per_file: usize,
    suffix: &String,
) -> anyhow::Result<()> {
    let file = File::open(&input)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let mut lines = reader.lines();

    create_dir_all(&out_dir)?;

    let mut file_index = 0;
    let mut record_index = 0;
    let mut writer = new_writer(&out_dir, file_index, &suffix)?;

    while let Some(h) = lines.next() {
        let header = h?;
        let seq = lines
            .next()
            .ok_or_else(|| anyhow::anyhow!("ERROR: unexpected EOF"))??;
        let plus = lines
            .next()
            .ok_or_else(|| anyhow::anyhow!("ERROR: unexpected EOF"))??;
        let qual = lines
            .next()
            .ok_or_else(|| anyhow::anyhow!("ERROR: unexpected EOF"))??;

        if record_index == records_per_file {
            file_index += 1;
            writer = new_writer(&out_dir, file_index, &suffix)?;
            record_index = 0;
        }

        writeln!(writer, "{header}\n{seq}\n{plus}\n{qual}")?;
        record_index += 1;
    }

    Ok(())
}

/// Counts the total number of FASTQ records in a gzipped FASTQ file.
///
/// This function reads through the entire gzipped FASTQ file line by line.
/// Since each FASTQ record consists of exactly four lines (header, sequence,
/// plus line, quality scores), the total number of records is determined by
/// dividing the total line count by four.
///
/// # Arguments
///
/// * `input` - A reference to the path of the input gzipped FASTQ file.
///
/// # Returns
///
/// * `anyhow::Result<usize>` - An `Ok(count)` containing the total number of
///                            FASTQ records, or an `anyhow::Error` if the
///                            file cannot be opened or read.
///
/// # Errors
///
/// * Returns an error if the input file cannot be opened or if there are any
///   issues during decompression and reading of lines.
///
/// # Example
///
/// ```rust, no_run
/// use anyhow::Result;
/// use std::path::PathBuf;
///
/// fn main() -> Result<()> {
///     // Assuming an existing 'input.fastq.gz'
///     let input_file = PathBuf::from("input.fastq.gz");
///     // let num_records = count_fastq_records(&input_file)?;
///     // println!("Total FASTQ records: {}", num_records);
///     Ok(())
/// }
/// ```
fn count_fastq_records<P: AsRef<Path>>(input: P) -> anyhow::Result<usize> {
    let file = File::open(&input)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let lines = reader.lines();

    let count = lines.count() / 4;
    Ok(count)
}

/// Creates and returns a new `BufWriter` wrapped around a `GzEncoder` for a new gzipped FASTQ file.
///
/// This helper function constructs a `PathBuf` for the new output file, using
/// the provided output directory, an index for numbering, and a suffix.
/// It then creates the file, initializes a `GzEncoder` with fast compression,
/// and wraps it in a `BufWriter` for efficient writing.
///
/// # Arguments
///
/// * `out_dir` - A reference to the path of the output directory.
/// * `index` - The zero-based index for the output file, used for sequential naming
///             (e.g., `tmp_chunk_000`, `tmp_chunk_001`).
/// * `suffix` - A reference to a `String` containing an optional suffix to append
///              to the filename before the `.fastq.gz` extension.
///
/// # Returns
///
/// * `anyhow::Result<BufWriter<GzEncoder<File>>>` - An `Ok` containing the configured
///                                                 writer, or an `anyhow::Error` if
///                                                 file creation or encoder initialization fails.
///
/// # Errors
///
/// * Returns an error if the output file cannot be created at the specified path.
///
/// # Example
///
/// ```rust, no_run
/// use anyhow::Result;
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::BufWriter;
/// use flate2::write::GzEncoder;
/// use flate2::Compression;
///
/// fn main() -> Result<()> {
///     let out_dir = PathBuf::from("my_output_dir");
///     std::fs::create_dir_all(&out_dir)?; // Ensure directory exists
///
///     let file_index = 0;
///     let file_suffix = "sample".to_string();
///
///     let mut writer = new_writer(&out_dir, file_index, &file_suffix)?;
///     writer.write_all(b"@header\nATGC\n+\n!!!!\n")?;
///     writer.flush()?;
///     writer.finish() // would be called by the GzEncoder on drop or explicit finish()
///     println!("New gzipped FASTQ writer created and used.");
///     Ok(())
/// }
/// ```
fn new_writer<P: AsRef<Path>>(
    out_dir: P,
    index: usize,
    suffix: &String,
) -> anyhow::Result<BufWriter<GzEncoder<File>>> {
    let path = out_dir
        .as_ref()
        .join(format!("tmp_chunk_{index:03}_{suffix}.fastq.gz"));
    let file = File::create(path)?;
    let encoder = GzEncoder::new(file, Compression::fast());
    Ok(BufWriter::new(encoder))
}
