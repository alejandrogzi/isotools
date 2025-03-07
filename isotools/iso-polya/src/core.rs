use bigtools::BigWigWrite;
use config::{get_progress_bar, Sequence, Strand};
use dashmap::DashSet;
use iso_classify::core::Genome;
use packbed::{par_reader, unpack};
use rayon::prelude::*;

use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};

use crate::{
    cli::{AparentArgs, FilterArgs, CHUNK_SIZE},
    utils::{bg_par_reader, get_sequences, PolyAPred},
};

pub const PARA: &str = "para";
pub const APPARENT_PY: &str = "run_aparent.py";
pub const ISO_POLYA: &str = "iso-polya";
pub const ASSETS: &str = "assets";
pub const RAM_PER_SITE: f32 = 0.025;
pub const JOBLIST: &str = "joblist";
pub const FILTER_MINIMAP: &str = "filterMinimapQuality.perl";
const TOKIO_RUNTIME_THREADS: usize = 8;

pub fn calculate_polya(args: AparentArgs) -> Result<(), Box<dyn std::error::Error>> {
    let isoseqs = unpack::<PolyAPred, _>(args.bed, config::OverlapType::Exon, true)
        .expect("ERROR: Could not unpack bed file!");
    let (genome, chrom_sizes) =
        get_sequences(args.twobit).expect("ERROR: Could not get read .2bit file!");

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
        // INFO: check if a .para dir is present
        check_para_dir();
        submit_jobs(paths, args.use_max_peak);

        // INFO: wait until all jobs are don and then merge the results
        merge_results(chrom_sizes);
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
    let joblist = create_joblist_aparent(accumulator, max_peak);

    let code = std::process::Command::new(PARA)
        .arg("make")
        .arg("aparent")
        .arg(joblist)
        .arg("-q")
        .arg("shortmed")
        .arg("-memoryMb")
        .arg(mem.to_string())
        .output()
        .expect("ERROR: Failed to submit job");

    if !code.status.success() {
        let err = String::from_utf8_lossy(&code.stderr);
        log::error!("ERROR: Job failed to submit! {}", err);
        std::process::exit(1);
    } else {
        log::info!("SUCCESS: Jobs submitted!");
    }

    log::info!("INFO: Jobs finished successfully!");

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

fn create_joblist_aparent(accumulator: &ParallelAccumulator, max_peak: bool) -> PathBuf {
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

fn merge_results(chrom_sizes: HashMap<String, u32>) {
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
                Some("bedGraph") => bgs.push(path),
                _ => {}
            }
        }
    }

    let bed = par_reader(beds).expect("ERROR: Failed to merge bed files");
    let (bg_plus, bg_minus) = bg_par_reader(bgs).expect("ERROR: Failed to merge bedGraph files");

    // INFO: write bed and bedgraph files and remove the chunks
    let bed_dest = assets.join("iso_polya_aparent.bed");
    let bg_dest_plus = assets.join("iso_polya_aparent_plus.bedGraph");
    let bg_dest_minus = assets.join("iso_polya_aparent_minus.bedGraph");
    let bw_dest_plus = assets.join("iso_polya_aparent_plus.bw");
    let bw_dest_minus = assets.join("iso_polya_aparent_minus.bw");

    vec![
        (bed_dest, bed),
        (bg_dest_plus.clone(), bg_plus),
        (bg_dest_minus.clone(), bg_minus),
    ]
    .into_par_iter()
    .for_each(|(dest, data)| {
        let file = File::create(&dest).expect("ERROR :Failed to create file");
        let mut writer = BufWriter::new(file);

        writer
            .write_all(data.as_bytes())
            .expect("ERROR: Failed to write file");
    });

    log::info!("INFO: Merged chunks and cleaning...");

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

    log::info!("INFO: Writing BigWig file from bedGraph fragments...");

    vec![(bg_dest_plus, bw_dest_plus), (bg_dest_minus, bw_dest_minus)]
        .into_par_iter()
        .for_each(|(bg_dest, bw_dest)| {
            let runtime = tokio::runtime::Builder::new_multi_thread()
                .worker_threads(TOKIO_RUNTIME_THREADS)
                .build()
                .expect("Unable to create runtime.");

            write_bigwig(bg_dest, bw_dest, chrom_sizes.clone(), runtime);
        });

    log::info!("SUCCESS: APPARENT finished successfully!");
}

fn write_bigwig(
    bg_dest: PathBuf,
    bw_dest: PathBuf,
    chrom_sizes: HashMap<String, u32>,
    handle: tokio::runtime::Runtime,
) {
    let bedgraph = File::open(bg_dest.clone()).expect("ERROR: Failed to open bedgraph file");
    let bw = BigWigWrite::create_file(bw_dest.to_string_lossy().as_ref(), chrom_sizes)
        .expect("ERROR: Cannot create BigWig file");

    let data = bigtools::beddata::BedParserStreamingIterator::from_bedgraph_file(bedgraph, true);
    bw.write(data, handle)
        .expect("ERROR: Failed to write BigWig file");

    std::fs::remove_file(bg_dest).expect("ERROR: Failed to remove bedGraph file");
}

fn check_para_dir() {
    let para = std::env::current_dir()
        .expect("ERROR: Failed to get current directory")
        .join(format!(".{}", PARA));

    if para.exists() {
        std::fs::remove_dir_all(&para).expect("ERROR: Failed to remove .para directory");
    }
}

pub fn filter_minimap(args: FilterArgs) -> Result<(), Box<dyn std::error::Error>> {
    let joblist = create_job_filter_minimap(args.clone());

    if args.para {
        let mut additional_args = vec![];
        if let Some(queue) = args.queue {
            additional_args.push(format!("-q {}", queue));
        }

        if let Some(mem) = args.mem {
            additional_args.push(format!("-memoryMb {}", mem));
        }

        let code = std::process::Command::new(PARA)
            .arg("make")
            .arg("filter_minimap")
            .arg(joblist)
            .args(additional_args)
            .output()
            .expect("ERROR: Failed to submit job");

        if !code.status.success() {
            let err = String::from_utf8_lossy(&code.stderr);
            log::error!("ERROR: Job failed to submit! {}", err);
            std::process::exit(1);
        } else {
            log::info!("SUCCESS: Jobs submitted!");
        }

        log::info!("INFO: Jobs finished successfully!");
    } else {
        log::info!(
            "{}",
            format!(
                "INFO: Only writing joblist to: {}",
                joblist.to_string_lossy()
            )
        );
    }

    Ok(())
}

fn create_job_filter_minimap(args: FilterArgs) -> PathBuf {
    let joblist = PathBuf::from(JOBLIST);
    let file = File::create(joblist.clone()).expect("Failed to create joblist");
    let mut writer = BufWriter::new(file);
    let executable = get_assets_dir().join(FILTER_MINIMAP);

    let mut cmd = format!(
        "perl {} {} -perID {} -clip3 {} -clip5 {} -P2P {} -emitA {}",
        executable.to_string_lossy(),
        args.sam
            .get(0)
            .expect("ERROR: No sam file provided")
            .to_string_lossy(),
        args.per_id,
        args.clip3,
        args.clip5,
        args.p2p,
        args.emit_a
    );

    if args.stat {
        cmd.push_str(" -statFile");
    }

    if args.keep {
        cmd.push_str(" -keepBad5Prime");
    }

    if let Some(suffix) = args.suffix {
        let suffix = format!(" -polyAReadSuffix {}", suffix);
        cmd.push_str(&suffix);
    }

    let _ = writer.write_all(cmd.as_bytes());
    let _ = writer.write_all(b"\n");

    let _ = writer.flush();

    return joblist;
}

pub fn pas_caller() {
    unimplemented!();
}
