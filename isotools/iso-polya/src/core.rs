use config::{get_progress_bar, Sequence, Strand};
use dashmap::DashSet;
use iso_classify::{core::Genome, utils::get_sequences};
use packbed::{par_reader, unpack};
use rayon::prelude::*;

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};

use crate::{cli::Args, utils::PolyAPred};

pub const CHUNK_SIZE: usize = 500;
pub const PARA: &str = "para";
pub const APPARENT_PY: &str = "apparent.py";
pub const ISO_POLYA: &str = "iso-polya";
pub const ASSETS: &str = "assets";
pub const RAM_PER_SITE: f32 = 0.025;
pub const JOBLIST: &str = "joblist";

pub fn calculate_polya(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    let isoseqs = unpack::<PolyAPred, _>(args.bed, config::OverlapType::Exon, true)
        .expect("ERROR: Could not unpack bed file!");
    let genome = get_sequences(args.twobit).expect("ERROR: Could not get read .2bit file!");

    let pb = get_progress_bar(isoseqs.len() as u64, "Processing reads...");
    let accumulator = ParallelAccumulator::default();

    // INFO: a bucket represents Chromosome [String] -> Reads [Vec<PolyAPred>]
    isoseqs.into_par_iter().for_each(|bucket| {
        let chr = bucket.0;
        let reads = bucket.1;

        distribute(reads, &genome, chr, &accumulator);

        pb.inc(1);
    });

    pb.finish_and_clear();

    let paths = chunk_writer(&accumulator);

    // INFO: push each path to cluster and run APPARENT
    // INFO: in the future the arg should be parsed as Executor::Something
    if args.para {
        submit_jobs(paths, args.use_max_peak);

        // INFO: wait until all jobs are don and then merge the results
        merge_results();
    }

    Ok(())
}

#[inline(always)]
fn distribute(
    reads: Vec<PolyAPred>,
    genome: &Genome,
    chr: String,
    accumulator: &ParallelAccumulator,
) {
    reads.into_par_iter().for_each(|read| {
        // INFO: getting the sequence from the genome for each read
        let seq = match read.strand {
            Strand::Forward => Sequence::new(
                genome
                    .get(&chr)
                    .expect("ERROR: Chromosome not in assembly, check your .2bit!")
                    [read.start as usize..read.end as usize]
                    .as_ref(),
            ),
            Strand::Reverse => Sequence::new(
                genome
                    .get(&chr)
                    .expect("ERROR: Chromosome not in assembly, check your .2bit!")
                    [read.start as usize..read.end as usize]
                    .as_ref(),
            )
            .reverse_complement(),
        };

        let row = format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            read.chrom, read.start, read.end, read.name, read.strand, seq
        );
        accumulator.lines.insert(row);
    });
}

struct ParallelAccumulator {
    lines: DashSet<String>,
    paths: DashSet<String>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            lines: DashSet::new(),
            paths: DashSet::new(),
        }
    }
}

fn chunk_writer(accumulator: &ParallelAccumulator) -> &ParallelAccumulator {
    let counter = AtomicUsize::new(0);

    accumulator
        .lines
        .par_iter()
        .collect::<Vec<_>>()
        .chunks(CHUNK_SIZE)
        .for_each(|chunk| {
            let index = counter.fetch_add(1, Ordering::Relaxed);
            let filename = PathBuf::from(format!("chunk_{}", index));
            let dest = get_assets_dir().join(filename);

            if let Ok(file) = File::create(&dest) {
                let mut writer = BufWriter::new(file);

                for line in chunk {
                    let _ = writer.write_all(line.as_bytes());
                    let _ = writer.write_all(b"\n");
                }

                let _ = writer.flush();
                accumulator.paths.insert(dest.to_string_lossy().to_string());
            }
        });

    accumulator
}

fn submit_jobs(accumulator: &ParallelAccumulator, max_peak: bool) {
    let mem = CHUNK_SIZE as f32 * RAM_PER_SITE * 1024.0;
    let joblist = create_joblist(accumulator, max_peak);

    let code = std::process::Command::new(PARA)
        .arg("make")
        .arg("apparent")
        .arg(joblist)
        .arg("-q")
        .arg("shortmed")
        .arg("-memoryMb")
        .arg(mem.to_string())
        .output();

    if let Ok(_) = code {
        log::info!("SUCCESS: Job submitted and finished succesfully!");
    } else {
        log::error!("ERROR: Job failed to submit!");
    }

    // INFO: removing chunks!
    for entry in std::fs::read_dir(get_assets_dir())
        .expect("Failed to read assets directory")
        .flatten()
    {
        let path = entry.path();
        if path
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .starts_with("chunk_")
        {
            let _ = std::fs::remove_file(path);
        }
    }
}

fn create_joblist(accumulator: &ParallelAccumulator, max_peak: bool) -> PathBuf {
    let joblist = PathBuf::from(JOBLIST);
    let file = File::create(joblist.clone()).expect("Failed to create joblist");
    let mut writer = BufWriter::new(file);
    let executable = get_assets_dir().join(APPARENT_PY);

    for path in accumulator.paths.iter() {
        let mut cmd = format!(
            "python3 {} -p {}",
            executable.to_string_lossy(),
            path.as_str()
        );

        if max_peak {
            cmd.push_str(" -mp");
        }

        let _ = writer.write_all(cmd.as_bytes());
        let _ = writer.write_all(b"\n");
    }

    let _ = writer.flush();

    return joblist;
}

fn get_assets_dir() -> PathBuf {
    let mut assets = std::env::current_dir().expect("Failed to get executable path");

    if !assets.ends_with(ISO_POLYA) {
        let rest = PathBuf::from(ISO_POLYA).join(ASSETS);
        assets.push(rest);

        return assets;
    } else {
        return assets.join(ASSETS);
    }
}

fn merge_results() {
    let assets = get_assets_dir();

    let mut beds = Vec::new();
    let mut bgs = Vec::new();

    for entry in std::fs::read_dir(assets.clone())
        .expect("ERROR: Failed to read assets directory")
        .flatten()
    {
        let path = entry.path();
        if let Some(ext) = path.extension() {
            match ext.to_str() {
                Some("bed") => beds.push(path),
                Some("bedgraph") => bgs.push(path),
                _ => {}
            }
        }
    }

    let bed = par_reader(beds).expect("ERROR: Failed to merge bed files");
    let bg = par_reader(bgs).expect("ERROR: Failed to merge bedGraph files");

    // INFO: write bed and bedgraph files and remove the chunks
    let bed_dest = assets.join("iso_polya_aparent.bed");
    let bg_dest = assets.join("iso_polya_aparent.bedGraph");

    vec![(bed_dest, bed), (bg_dest, bg)]
        .into_par_iter()
        .for_each(|(dest, data)| {
            let file = File::create(&dest).expect("ERROR :Failed to create file");
            let mut writer = BufWriter::new(file);

            writer
                .write_all(data.as_bytes())
                .expect("ERROR: Failed to write file");
        });

    // INFO: deleting APPARENT chunked output
    for entry in std::fs::read_dir(assets).expect("ERROR: Failed to read assets directory") {
        if let Ok(entry) = entry {
            let path = entry.path();
            if path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .starts_with("polya")
            {
                let _ = std::fs::remove_file(path);
            }
        }
    }
}
