use bigtools::BigWigWrite;
use config::{get_progress_bar, Sequence, Strand};
use dashmap::DashSet;
use iso_classify::core::Genome;
use packbed::{par_reader, unpack};
use rand::Rng;
use rayon::prelude::*;

use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};

use crate::{
    cli::{AparentArgs, CHUNK_SIZE},
    utils::{bg_par_reader, get_assets_dir, get_sequences, MiniPolyAPred},
};

const PARA: &str = "para";
const APPARENT_PY: &str = "run_aparent.py";
const RAM_PER_SITE: f32 = 0.025;
const JOBLIST: &str = "joblist";
const TOKIO_RUNTIME_THREADS: usize = 8;

pub fn calculate_polya(args: AparentArgs) -> Result<(), Box<dyn std::error::Error>> {
    let isoseqs = unpack::<MiniPolyAPred, _>(args.bed, config::OverlapType::Exon, true)
        .expect("ERROR: Could not unpack bed file!");
    let (genome, chrom_sizes) = get_sequences(
        args.twobit
            .expect("ERROR: Provide a .2bit file using --twobit!"),
    )
    .expect("ERROR: Could not get read .2bit file!");

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
    reads: Vec<MiniPolyAPred>,
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

/// ParallelAccumulator struct
///
/// # Fields
///
/// * `lines` - DashSet with the lines to be written
/// * `paths` - DashSet with the paths to the chunks
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.lines.len(), 0);
/// ```
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

/// Write lines in ParallelAccumulator to chunks of a given size
///
/// # Arguments
///
/// * `accumulator` - ParallelAccumulator struct
///
/// # Example
///
/// ```rust, no_run
/// chunk_writer(accumulator);
///
/// assert_eq!(std::fs::metadata("chunk_0").is_ok(), true);
/// ```
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

/// Submit the APPARENT jobs to the cluster
///
/// # Arguments
///
/// * `accumulator` - ParallelAccumulator struct
/// * `max_peak` - bool to use the max peak
///
/// # Example
///
/// ```rust, no_run
/// submit_jobs(accumulator, max_peak);
///
/// assert_eq!(std::fs::metadata("joblist").is_ok(), true);
/// ```
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

/// Create a joblist file for APPARENT
///
/// # Arguments
///
/// * `accumulator` - ParallelAccumulator struct
/// * `max_peak` - bool to use the max peak
///
/// # Example
///
/// ```rust, no_run
/// let joblist = create_joblist_aparent(accumulator, max_peak);
///
/// assert_eq!(std::fs::metadata("joblist").is_ok(), true);
/// ```
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

/// Write a bed file from APARENT predictions
///
/// # Arguments
///
/// * `bed_dest` - PathBuf to the bed file
/// * `bed` - String with the bed data
///
/// # Example
///
/// ```rust, no_run
/// bed_write(bed_dest, bed);
///
/// assert_eq!(std::fs::metadata("iso_polya_aparent.bed").is_ok(), true);
/// ```
fn bed_write(bed_dest: PathBuf, bed: String) {
    let file = File::create(&bed_dest).expect("ERROR: Failed to create bed file");
    let mut writer = BufWriter::new(file);
    writer
        .write_all(bed.as_bytes())
        .expect("ERROR: Failed to write bed file");
}

/// Write bedGraph files from APARENT predictions
///
/// # Arguments
///
/// * `bg_dest_plus` - PathBuf to the bedGraph file for the plus strand
/// * `bg_plus` - String with the bedGraph data for the plus strand
/// * `bg_dest_minus` - PathBuf to the bedGraph file for the minus strand
/// * `bg_minus` - String with the bedGraph data for the minus strand
///
/// # Example
///
/// ```rust, no_run
/// bg_write(bg_dest_plus, bg_plus, bg_dest_minus, bg_minus);
///
/// assert_eq!(std::fs::metadata("iso_polya_aparent_plus.bedGraph").is_ok(), true);
/// ```
fn bg_write(bg_dest_plus: PathBuf, bg_plus: String, bg_dest_minus: PathBuf, bg_minus: String) {
    vec![(bg_dest_plus, bg_plus), (bg_dest_minus, bg_minus)]
        .into_par_iter()
        .for_each(|(dest, data)| {
            let file = File::create(&dest).expect("ERROR: Failed to create bedGraph file");
            let mut writer = BufWriter::new(file);
            writer
                .write_all(data.as_bytes())
                .expect("ERROR: Failed to write bedGraph file");
        });
}

/// Write BigWig files from APARENT predictions
///
/// # Arguments
///
/// * `bg_dest_plus` - PathBuf to the bedGraph file for the plus strand
/// * `bw_dest_plus` - PathBuf to the BigWig file for the plus strand
/// * `bg_dest_minus` - PathBuf to the bedGraph file for the minus strand
/// * `bw_dest_minus` - PathBuf to the BigWig file for the minus strand
///
/// # Example
///
/// ```rust, no_run
/// bw_write(bg_dest_plus, bw_dest_plus, bg_dest_minus, bw_dest_minus, chrom_sizes);
///
/// assert_eq!(std::fs::metadata("iso_polya_aparent_plus.bw").is_ok(), true);
/// ```
fn bw_write(
    bg_dest_plus: PathBuf,
    bw_dest_plus: PathBuf,
    bg_dest_minus: PathBuf,
    bw_dest_minus: PathBuf,
    chrom_sizes: HashMap<String, u32>,
) {
    vec![(bg_dest_plus, bw_dest_plus), (bg_dest_minus, bw_dest_minus)]
        .into_par_iter()
        .for_each(|(bg_dest, bw_dest)| {
            let runtime = tokio::runtime::Builder::new_multi_thread()
                .worker_threads(TOKIO_RUNTIME_THREADS)
                .build()
                .expect("Unable to create runtime.");
            write_bigwig(bg_dest, bw_dest, chrom_sizes.clone(), runtime);
        });
}

/// Merge the results from the APPARENT chunks
///
/// # Arguments
///
/// * `chrom_sizes` - HashMap with the chromosome names and sizes
///
/// # Example
///
/// ```rust, no_run
/// merge_results(chrom_sizes);
///
/// assert_eq!(std::fs::metadata("iso_polya_aparent.bed").is_ok(), true);
/// ```
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

    let bed_dest = assets.join("iso_polya_aparent.bed");
    let bg_dest_plus = assets.join("iso_polya_aparent_plus.bedGraph");
    let bg_dest_minus = assets.join("iso_polya_aparent_minus.bedGraph");
    let bw_dest_plus = assets.join("iso_polya_aparent_plus.bw");
    let bw_dest_minus = assets.join("iso_polya_aparent_minus.bw");

    bed_write(bed_dest, bed);
    bg_write(
        bg_dest_plus.clone(),
        bg_plus,
        bg_dest_minus.clone(),
        bg_minus,
    );

    log::info!("INFO: Merged chunks and cleaning...");
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
    bw_write(
        bg_dest_plus,
        bw_dest_plus,
        bg_dest_minus,
        bw_dest_minus,
        chrom_sizes,
    );
    log::info!("SUCCESS: APPARENT finished successfully!");
}

/// Write a BigWig file from a bedGraph file
///
/// # Arguments
///
/// * `bg_dest` - PathBuf to the bedGraph file
/// * `bw_dest` - PathBuf to the BigWig file
/// * `chrom_sizes` - HashMap with the chromosome names and sizes
/// * `handle` - tokio runtime handle
///
/// # Example
///
/// ```rust, no_run
/// write_bigwig(bg_dest, bw_dest, chrom_sizes, handle);
///
/// assert_eq!(std::fs::metadata("iso_polya_aparent_plus.bw").is_ok(), true);
/// ```
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

/// Check if a .para directory is present and remove it
///
/// # Example
///
/// ```rust, no_run
/// check_para_dir();
///
/// assert_eq!(std::fs::metadata(".para").is_err(), true);
/// ```
fn check_para_dir() {
    let para = std::env::current_dir()
        .expect("ERROR: Failed to get current directory")
        .join(format!(".{}", PARA));

    if para.exists() {
        std::fs::remove_dir_all(&para).expect("ERROR: Failed to remove .para directory");
    }
}

/// Simulate reads from/to a given bp space
/// where the last k bp are polyadenylated on purpose
///
/// # Arguments
///
/// * `args` - AparentArgs struct
///
/// # Example
///
/// ```rust, no_run
/// simulate_polya_reads(args);
/// ```
pub fn simulate_polya_reads(args: AparentArgs) -> Result<(), Box<dyn std::error::Error>> {
    let accumulator = ParallelAccumulator::default();

    make_reads(args.clone(), &accumulator);
    let paths = chunk_writer(&accumulator);

    // INFO: push each path to cluster and run APPARENT
    // INFO: in the future the arg should be parsed as Executor::Something
    if args.para {
        // INFO: check if a .para dir is present
        check_para_dir();
        submit_jobs(paths, args.use_max_peak);

        // INFO: wait until all jobs are don and then merge the results
        let chrom_sizes = make_chrom_sizes();
        merge_results(chrom_sizes);
    }

    Ok(())
}

/// Simulate reads from/to a given bp space
///
/// # Arguments
///
/// * `args` - AparentArgs struct
/// * 'accumulator' - ParallelAccumulator struct
///
/// # Example
///
/// ```rust, no_run
/// make_reads(args, accumulator);
/// ```
fn make_reads(args: AparentArgs, accumulator: &ParallelAccumulator) {
    log::info!("INFO: Simulating reads...");

    let mut rng = rand::thread_rng();
    let chroms = make_chrom_sizes();

    for counter in 0..args.number_of_reads {
        let chrom = chroms
            .keys()
            .nth(rng.gen_range(0..chroms.len()))
            .expect("ERROR: Chromosome not found!");
        let size = chroms.get(chrom).expect("ERROR: Chromosome not found!");

        if *size < args.read_length as u32 {
            continue;
        }

        let start = rng.gen_range(0..size.saturating_sub(args.read_length as u32));
        let end = start + args.read_length as u32;

        let strand = if args.stranded {
            if rng.gen_bool(0.5) {
                Strand::Reverse
            } else {
                Strand::Forward
            }
        } else {
            Strand::Forward
        };

        let polya_length = rng.gen_range(0..=args.polya_range);
        let seq = Sequence::random(args.read_length - polya_length).fill_back(polya_length);

        let name = format!("read_{}_{}", counter, polya_length);

        let row = format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            chrom, start, end, name, strand, seq
        );
        accumulator.lines.insert(row);
    }
}

/// Create a HashMap with the chromosome names and sizes
/// for the simulation
///
/// # Example
///
/// ```rust, no_run
/// let chrom_sizes = make_chrom_sizes();
///
/// assert_eq!(chrom_sizes.get("1"), Some(&100000));
/// ```
fn make_chrom_sizes() -> HashMap<String, u32> {
    let chroms = vec!["chr1", "chr2", "chr3"];
    let sizes = vec![100000, 200000, 300000];

    let mut chrom_sizes = HashMap::new();
    chroms.iter().zip(sizes.iter()).for_each(|(chrom, size)| {
        chrom_sizes.insert(chrom.to_string(), *size);
    });

    chrom_sizes
}
