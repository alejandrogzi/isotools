use std::{
    fs::{create_dir_all, File},
    io::{BufRead, BufReader, BufWriter, Write},
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

#[macro_export]
macro_rules! dispatch {
    ($file:expr, { $($suffix:literal => $action:expr),* $(,)? }) => {{
        let f = $file.file_name().and_then(|f| f.to_str()).unwrap_or_default();
        $(
            if f.ends_with($suffix) {
                $action
            } else
        )* {
            anyhow::bail!("ERROR: unrecognized file format: {}", $file.display());
        }
    }};
}

pub fn lib_iso_split(args: Vec<String>) -> Result<()> {
    let args = cli::Args::from(args);

    let _ = dispatch!(args.file,
        {
            ".fa" => split_fa(args),
            ".fasta" => split_fa(args),
            ".fa.gz" => split_fa(args),
            ".fasta.gz" => split_fa(args),
            ".fq.gz" => split_fq(args),
            ".fastq.gz" => split_fq(args),
        }
    );

    Ok(())
}

pub fn split_fa(args: Args) -> Result<()> {
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
    let suffix = args.suffix.unwrap_or(String::new());

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
                        PathBuf::from(&args.outdir).join(format!("chunk_{:03}_{suffix}.fa", i));
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
                        PathBuf::from(&args.outdir).join(format!("part_{:03}_{suffix}.fa", i));
                    let mut writer = BufWriter::new(File::create(output_file)?);
                    writer.write_all(&data[chunk.start..chunk.end])?;
                    writer.flush()?;
                    Ok(())
                })?;

            Ok(())
        }
    }
}

pub fn split_fq(args: Args) -> anyhow::Result<()> {
    log::info!("INFO: running in FASTQ mode with args: {:?}", &args);

    let mode = args.mode()?;
    let suffix = args.suffix.unwrap_or(String::new());

    match mode {
        SplitMode::ChunkSize(records_per_file) => {
            split_by_chunk_size(args.file, args.outdir, records_per_file, suffix)
        }
        SplitMode::NumFiles(num_files) => {
            let total_records = count_fastq_records(&args.file)?;
            let records_per_file = (total_records + num_files - 1) / num_files;
            split_by_chunk_size(args.file, args.outdir, records_per_file, suffix)
        }
    }
}

fn split_by_chunk_size<P: AsRef<Path>>(
    input: P,
    out_dir: P,
    records_per_file: usize,
    suffix: String,
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

fn count_fastq_records<P: AsRef<Path>>(input: P) -> anyhow::Result<usize> {
    let file = File::open(&input)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let lines = reader.lines();

    let count = lines.count() / 4;
    Ok(count)
}

fn new_writer<P: AsRef<Path>>(
    out_dir: P,
    index: usize,
    suffix: &String,
) -> anyhow::Result<BufWriter<GzEncoder<File>>> {
    let path = out_dir
        .as_ref()
        .join(format!("chunk_{index:03}_{suffix}.fastq.gz"));
    let file = File::create(path)?;
    let encoder = GzEncoder::new(file, Compression::default());
    Ok(BufWriter::new(encoder))
}
