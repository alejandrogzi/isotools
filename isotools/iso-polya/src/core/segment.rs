use std::{
    fmt::Write as FmtWrite,
    fs::File,
    hash::{Hash, Hasher},
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Arc,
};

use crate::cli::SegmentArgs;

use bio::stats::hmm::{discrete_emission::Model, viterbi, State};
use config::*;
use dashmap::DashSet;
use ndarray::{array, Array2};
use noodles_bam::{io::writer::Builder, Record};
use noodles_core::{region::Interval, Position};
use noodles_sam::{
    alignment::{
        io::Write as BamWrite,
        record::{cigar::op::Kind, Cigar},
        RecordBuf,
    },
    Header,
};
use rayon::prelude::*;

pub fn segment(args: SegmentArgs) -> Result<(), String> {
    let (chroms, refs, header) = read_bam(&args.bam);
    let accumulator = ParallelAccumulator::default();

    let pb = get_progress_bar(chroms.len() as u64, "Processing reads...");

    // INFO: process each reference region in parallel
    chroms.par_iter().for_each(|chr| {
        // INFO: thread-local copy of the reference sequences
        let references = Arc::clone(&refs);

        let mut reader = noodles_bam::io::indexed_reader::Builder::default()
            .build_from_path(&args.bam)
            .expect(&format!(
                "ERROR: could not open file: {}",
                args.bam.display()
            ));

        // INFO: must re-read header per thread since IndexedReader owns its header
        let header = reader.read_header().expect(&format!(
            "ERROR: could not read header for file: {}",
            args.bam.display()
        ));

        let mut track = 0;

        let end = Position::new(references[chr.as_bytes()].length().get()).unwrap();
        let bounds = Interval::from(Position::MIN..=end);
        let region = noodles_core::Region::new(chr.as_str(), bounds);

        reader
            .query(&header, &region)
            .expect(&format!("ERROR: could not query region: {}", region))
            .for_each(|result| match result {
                Ok(record) => {
                    track += 1;
                    process_record(record, &header, track, chr, &args, &accumulator);
                }
                Err(_) => {}
            });

        pb.inc(1);
    });

    pb.finish_and_clear();
    write_output(accumulator, &args, &header, refs);

    Ok(())
}

/// Writes processed reads to either BAM or BED format output files based on configuration
///
/// This function routes the accumulated reads to the appropriate output writer
/// based on the `bed` flag in the arguments. It handles both accepted and rejected
/// reads according to the output paths specified in the arguments.
///
/// # Arguments
///
/// * `accumulator` - Parallel accumulator containing processed reads (both accepted and rejected)
/// * `args` - Configuration parameters containing output paths and format selection
/// * `header` - SAM/BAM header reference (only used for BAM output)
///
/// # Behavior
///
/// - When `args.bed` is true:
///   * Writes output in BED format using `write_beds`
///   * Header is not used
/// - When `args.bed` is false:
///   * Writes output in BAM format using `write_bams`
///   * Uses the provided header for BAM writing
///
/// # Example
///
/// ```rust, no_run
/// use noodles_sam::header::Header;
///
/// let accumulator = ParallelAccumulator::new();
/// let args = SegmentArgs {
///     bed: true,  // Output as BED format
///     ..Default::default()
/// };
/// let header = Header::default();
///
/// write_output(accumulator, &args, &header);
/// ```
fn write_output(
    accumulator: ParallelAccumulator,
    args: &SegmentArgs,
    header: &Header,
    refs: Arc<noodles_sam::header::ReferenceSequences>,
) {
    match args.bed {
        true => write_beds(accumulator, args, refs),
        false => write_bams(accumulator, args, header),
    }
}

/// Writes a BAM file from a ParallelAccumulator
///
/// # Arguments
///
/// * `accumulator` - ParallelAccumulator collection
/// *  args - Module arguments
/// *  header - BAM file header
///
/// # Example
///
/// ```rust, no_run
/// write_bams(accumulator, args, header);
/// ```
fn write_bams(accumulator: ParallelAccumulator, args: &SegmentArgs, header: &Header) {
    log::info!("INFO: Writing BAM files from filtered records!");

    let header = Arc::new(header.clone());
    let prefix = &args.prefix;

    let accept = args.outdir.join(format!("{}.good.bam", prefix.display()));
    let reject = args.outdir.join(format!("{}.bad.bam", prefix.display()));

    let spawn_writer =
        |output_path: PathBuf, records: DashSet<HashedRecord>, header: Arc<Header>| {
            std::thread::spawn(move || {
                let mut writer = Builder::default()
                    .build_from_path(&output_path)
                    .unwrap_or_else(|_| {
                        panic!("ERROR: could not create file: {}", output_path.display())
                    });

                writer
                    .write_header(&header)
                    .expect("ERROR: failed to write header");

                for record in records {
                    let record = Arc::try_unwrap(record.0).unwrap_or_else(|arc| (*arc).clone());
                    writer
                        .write_alignment_record(&header, &record)
                        .expect("ERROR: failed to write record");
                }

                writer.finish(&header)
            })
        };

    let _ = spawn_writer(accept, accumulator.accept, Arc::clone(&header))
        .join()
        .expect("ERROR: could not join acceptance writer!");
    let _ = spawn_writer(reject, accumulator.reject, Arc::clone(&header))
        .join()
        .expect("ERROR: could not join rejection writer");
}

/// Writes a BED file from a ParallelAccumulator of BAM/SAM records
///
/// # Arguments
///
/// * `accumulator` - ParallelAccumulator collection
/// *  args - Module arguments
/// *  refs - BAM reference sequences index
///
/// # Example
///
/// ```rust, no_run
/// write_bams(accumulator, args, refs);
/// ```
fn write_beds(
    accumulator: ParallelAccumulator,
    args: &SegmentArgs,
    refs: Arc<noodles_sam::header::ReferenceSequences>,
) {
    let prefix = &args.prefix;

    let accept = args.outdir.join(format!("{}.good.bed", prefix.display()));
    let reject = args.outdir.join(format!("{}.bad.bed", prefix.display()));

    log::info!("INFO: Converting BAM records into BED12 lines!");
    let accept_lines = convert(&accumulator.accept, Arc::clone(&refs), RGB_ACCEPT);
    let reject_lines = convert(&accumulator.reject, Arc::clone(&refs), RGB_REJECT);

    fn spawn_writer(output_path: PathBuf, lines: Vec<String>) -> std::thread::JoinHandle<()> {
        std::thread::spawn(move || {
            let mut writer = BufWriter::new(File::create(&output_path).expect(&format!(
                "ERROR: could not create file: {}",
                output_path.display()
            )));

            let joined = lines.join("\n");
            writer.write_all(joined.as_bytes()).unwrap();
        })
    }

    log::info!("INFO: Writing BED files from filtered records!");

    let _ = spawn_writer(accept, accept_lines)
        .join()
        .expect("ERROR: failed to write good.bed");
    let _ = spawn_writer(reject, reject_lines)
        .join()
        .expect("ERROR: failed to write bad.bed");
}

/// Reads a BAM file and returns the chromosome names,
/// reference sequences, and header
///
/// # Arguments
///
/// * `bam` - Path to the BAM file
///
/// # Returns
///
/// * A tuple containing:
///     * `Vec<String>` - Vector of chromosome names
///     * `Arc<noodles_sam::header::ReferenceSequences>` - Reference sequences
///     * `noodles_sam::header::Header` - Header of the BAM file
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
///
/// let bam_path = PathBuf::from("path/to/bam/file.bam");
/// let (chroms, refs, header) = read_bam(&bam_path);
///
/// assert_eq!(chroms, vec!["chr1", "chr2"]);
/// ```
fn read_bam(
    bam: &PathBuf,
) -> (
    Vec<String>,
    Arc<noodles_sam::header::ReferenceSequences>,
    Header,
) {
    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .build_from_path(&bam)
        .expect(&format!("ERROR: could not open file: {}", bam.display()));

    let header = reader.read_header().unwrap();
    let refs = Arc::new(header.reference_sequences().clone()); // clone for sharing

    let chroms: Vec<String> = refs.iter().map(|(name, _)| name.to_string()).collect();

    return (chroms, refs, header);
}

/// Processes a single alignment record to analyze poly-A tails and apply quality filters
///
/// This function performs several operations:
/// 1. Filters out low-quality reads based on minimum identity threshold
/// 2. Predicts poly-A tail length using HMM models
/// 3. Optionally tags reads with additional information
/// 4. Distributes reads to accept/reject accumulators based on quality criteria
///
/// # Arguments
///
/// * `record` - The alignment record to process
/// * `header` - SAM/BAM header reference
/// * `track` - Track identifier used for read tagging
/// * `chr` - Chromosome name string reference
/// * `args` - Configuration parameters for processing
/// * `accumulator` - Parallel accumulator for collecting accepted/rejected reads
///
/// # Behavior
///
/// - Reads failing the minimum identity threshold (`args.min_identity`) are immediately discarded
/// - For reads with 3' end clipping:
///   - Predicts poly-A tail length using HMM if not hard-clipped
///   - Optionally adjusts prediction using suffix analysis if configured
/// - Updates the record's clipping information based on poly-A prediction
/// - Tags the read name if configured (`args.tag`)
/// - Distributes to accept/reject accumulators based on:
///   - Identity threshold (`args.identity`)
///   - Maximum 5' clip length (`args.max_clip_five`)
///   - Maximum 3' clip length (`args.max_clip_three`)
///
/// # Example
///
/// ```rust, no_run
/// use std::sync::Arc;
/// use noodles_sam::header::Header;
///
/// let record = get_alignment_record(); // Assume this gets a record
/// let header = Header::default();
/// let args = SegmentArgs::default();
/// let accumulator = ParallelAccumulator::new();
///
/// process_record(
///     record,
///     &header,
///     1,
///     &"chr1".to_string(),
///     &args,
///     &accumulator
/// );
/// ```
fn process_record(
    record: Record,
    header: &Header,
    track: u64,
    chr: &String,
    args: &SegmentArgs,
    accumulator: &ParallelAccumulator,
) {
    let mut read = Read::from(&record);

    // INFO: enforces minimum identity to filter
    // out low quality reads from the analysis
    if read.identity < args.min_identity {
        return;
    }

    let hmm = HMM::init(args.p2p, args.emit_a);

    if read.three_clip > 0 {
        if !read.has_hard_clip_three {
            let clipped_seq = read.get_rev_clipped_seq_sized();
            read.set_polya_len(predict_tail(clipped_seq, &hmm));
        }
    }

    if args.tail_suffix > 0 {
        let suffix_len = read.three_clip + args.tail_suffix;
        let mut polya_suffix_len =
            predict_tail_with_suffix(suffix_len, &read, &hmm, args.suffix_step_size);

        if read.polya_len > polya_suffix_len {
            polya_suffix_len = read.polya_len;
        }

        read.set_polya_read_len(polya_suffix_len);
    }

    read.set_three_clip(read.three_clip - read.polya_len);

    let mut record = RecordBuf::try_from_alignment_record(&header, &record).expect(&format!(
        "ERROR: could not convert record to RecordBuf: {:?}",
        record
    ));

    if args.tag {
        *record.name_mut() = Some(read.tag_read(track, chr).into());
    }

    if read.identity >= args.identity
        && read.five_clip <= args.max_clip_five
        && read.three_clip <= args.max_clip_three
    {
        accumulator.accept(Arc::from(record));
    } else {
        accumulator.reject(Arc::from(record));
    }
}

/// Predict the length of the polyA tail using a suffix
///
/// # Arguments
///
/// * `suffix_len` - Length of the suffix
/// * `read` - A reference to the Read object
/// * `hmm` - A reference to the HMM object
/// * `step_size` - Step size for the suffix length
///
/// # Returns
///
/// * `usize` - The length of the tail
///
/// # Example
///
/// ```rust, no_run
/// use iso::HMM;
///
/// let suffix_len = 10;
/// let read = Read::new();
/// let hmm = HMM::init(0.9, 0.99);
/// let step_size = 5;
/// let polya_len = predict_tail_with_suffix(suffix_len, &read, &hmm, step_size);
///
/// assert_eq!(polya_len, 10);
/// ```
fn predict_tail_with_suffix(
    mut suffix_len: usize,
    read: &Read,
    hmm: &HMM,
    step_size: usize,
) -> usize {
    loop {
        let suffix = read.get_rev_seq_suffix_sized(suffix_len);
        let polya_len = predict_tail(suffix, &hmm);

        if polya_len == suffix_len {
            suffix_len += step_size;
        } else {
            return polya_len;
        }
    }
}

/// Parallel accumulator for BAM records
///
/// # Fields
///
/// * `pass` - DashSet to store the pass reads
/// * `intrapriming` - DashSet to store the intrapriming reads
/// * `review` - DashSet to store the review reads
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.pass.len(), 0);
/// ```
#[derive(Debug)]
struct ParallelAccumulator {
    accept: DashSet<HashedRecord>,
    reject: DashSet<HashedRecord>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            accept: DashSet::new(),
            reject: DashSet::new(),
        }
    }
}

impl ParallelAccumulator {
    pub fn accept(&self, record: Arc<RecordBuf>) {
        self.accept.insert(HashedRecord(record));
    }

    pub fn reject(&self, record: Arc<RecordBuf>) {
        self.reject.insert(HashedRecord(record));
    }
}

/// A wrapper around Arc<Record> that implements Hash and Eq
#[derive(Debug)]
struct HashedRecord(Arc<RecordBuf>);

impl PartialEq for HashedRecord {
    fn eq(&self, other: &Self) -> bool {
        self.0.name() == other.0.name()
            && self.0.alignment_start().unwrap() == other.0.alignment_start().unwrap()
    }
}
impl Eq for HashedRecord {}

impl Hash for HashedRecord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.name().hash(state);
        self.0.alignment_start().unwrap().hash(state);
    }
}

impl HashedRecord {
    /// Convert a BAM record into a BED12 line
    ///
    /// # Fields
    ///
    /// * `chr` - Chromosome where the record belongs to
    /// * `rgb` - RGB color of the record in the bed file
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let chr = String::from("chr1");
    /// let rgb = Arc::from(String::from("255,0,0"));
    /// let record: HashedRecord = record;
    /// record.to_bed(chr, rgb);
    /// ```
    #[inline(always)]
    fn to_bed(&self, chr: &str, rgb: &str) -> Option<String> {
        let record = &self.0;

        let start = record.alignment_start()?.get() - 1;
        let score = record
            .mapping_quality()
            .map(|q| q.get())
            .unwrap_or(0)
            .to_string();

        let name = record.name().unwrap();
        let strand = if record.flags().is_reverse_complemented() {
            '-'
        } else {
            '+'
        };

        let mut blocks = Vec::new();
        let mut ref_pos = start;
        let mut block_start = ref_pos;
        let mut block_len = 0;

        let cigar = record.cigar();
        for c in cigar.iter() {
            let c = c.unwrap();
            let len = c.len();
            match c.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                    block_len += len;
                    ref_pos += len;
                }
                Kind::Skip => {
                    // INFO: end current block
                    if block_len > 0 {
                        blocks.push((block_start, block_len));
                    }

                    // INFO: skip the intron
                    ref_pos += len;

                    // INFO: start new block
                    block_start = ref_pos;
                    block_len = 0;
                }
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => continue,
            }
        }

        if block_len > 0 {
            blocks.push((block_start, block_len));
        }

        if blocks.is_empty() {
            return None;
        }

        let chrom_end = blocks.last().map(|(s, l)| s + l).unwrap_or(start);
        let thick_start = start;
        let thick_end = chrom_end;

        let mut block_sizes = String::with_capacity(16 * blocks.len());
        let mut block_starts = String::with_capacity(16 * blocks.len());

        for (s, l) in &blocks {
            let _ = write!(block_sizes, "{},", l);
            let _ = write!(block_starts, "{},", s - start);
        }
        block_sizes.pop();
        block_starts.pop();

        let mut line = String::with_capacity(256);
        write!(
            &mut line,
            "{chr}\t{start}\t{chrom_end}\t{name}\t{score}\t{strand}\t\
             {thick_start}\t{thick_end}\t{rgb}\t{}\t{block_sizes}\t{block_starts}",
            blocks.len()
        )
        .unwrap();

        Some(line)
    }
}

/// Convert a set of BAM records into a
/// set of BED records in parallel
///
/// # Arguments
///
/// * `records` - Set of BAM records as HashedRecords
/// * `refs` - Reference sequences from the BAM file
/// * `rgb` - HTML color code for the group
///
/// # Returns
///
/// * `Vec<String>` - BED records
///
/// # Example
///
/// ```rust, no_run
/// let bed_records = convert(records, refs, rgb);
/// ```
#[inline(always)]
fn convert(
    records: &DashSet<HashedRecord>,
    refs: Arc<noodles_sam::header::ReferenceSequences>,
    rgb: &str,
) -> Vec<String> {
    let pb = get_progress_bar(records.len() as u64, "Converting...");

    // WARN: filtering errors directly!
    records
        .par_iter()
        .filter_map(|record| {
            let chr = refs
                .get_index(record.0.reference_sequence_id()?)
                .map(|(name, _)| name.to_string())?;

            let line = record.to_bed(&chr, rgb);
            pb.inc(1);

            line
        })
        .collect()
}

/// Predict the length of the polyA tail
/// using a two-state Hidden Markov Model (HMM)
///
/// # Arguments
///
/// * `sequence` - A vector of usize representing the sequence
///
/// # Returns
///
/// * `usize` - The length of the tail
///
/// # Example
///
/// ```rust, no_run
/// use iso::HMM;
///
/// let sequence = vec![0, 1, 2, 3, 4];
/// let polya_len = predict_tail(sequence);
///
/// assert_eq!(polya_len, 5);
/// ```
fn predict_tail(sequence: Vec<usize>, hmm: &HMM) -> usize {
    let tail = hmm.segment(&sequence);

    return tail;
}

/// A two-state Hidden Markov Model (HMM)
/// for polyA tail detection
///
/// # Example
///
/// ```rust, no_run
/// use iso::HMM;
///
/// let hmm = HMM::init(0.9, 0.99);
/// let sequence = vec![0, 1, 2, 3, 4];
/// let polya_len = hmm.segment(&sequence);
///
/// assert_eq!(polya_len, 5);
/// ```
struct HMM {
    model: Model,
}

impl HMM {
    /// Initialize a two-state HMM for polyA tail detection
    ///
    /// # Arguments
    ///
    /// * `forward` - Probability of staying in state P
    /// * `weight` - Probability of emitting 'A' in state P
    ///
    /// # Returns
    ///
    /// * `HMM` - A new instance of the HMM
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let hmm = HMM::init(0.9, 0.99);
    /// let sequence = vec![0, 1, 2, 3, 4];
    /// let polya_len = hmm.segment(&sequence);
    ///
    /// assert_eq!(polya_len, 5);
    /// ```
    fn init(forward: f64, weight: f64) -> Self {
        let transition = HMM::__get_transition(forward);
        let emission = HMM::__get_emission(weight);
        let initial = array![0.5, 0.5];

        let model = Model::with_float(&transition, &emission, &initial).expect(&format!(
            "Failed to create HMM with transition: {:?}, emission: {:?}, initial: {:?}",
            transition, emission, initial
        ));

        Self { model }
    }

    /// Get the transition matrix for the HMM
    ///
    /// # Arguments
    ///
    /// * `forward` - Probability of staying in state P
    ///
    /// # Returns
    ///
    /// * `Array2<f64>` - Transition matrix
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let forward = 0.9;
    /// let transition = HMM::__get_transition(forward);
    ///
    /// assert_eq!(transition, array![[0.9, 0.1], [0.0, 1.0]]);
    /// ```
    fn __get_transition(forward: f64) -> Array2<f64> {
        array![
            [forward, 1.0 - forward], // INFO: from P to P and P to S
            [0.0, 1.0]                // INFO: from S to S (no transition back to P)
        ]
    }

    /// Get the emission matrix for the HMM
    ///
    /// # Arguments
    ///
    /// * `weight` - Probability of emitting 'A' in state P
    ///
    /// # Returns
    ///
    /// * `Array2<f64>` - Emission matrix
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let weight = 0.99;
    /// let emission = HMM::__get_emission(weight);
    ///
    /// assert_eq!(emission, array![[0.99, 0.0033333333333333335,
    /// 0.0033333333333333335, 0.0033333333333333335], [0.25; 4]]);
    /// ```
    fn __get_emission(weight: f64) -> Array2<f64> {
        let uniform = 0.25;

        let emit_a = weight;
        let emit_not_a = (1.0 - emit_a) / 3.0;

        array![
            [emit_a, emit_not_a, emit_not_a, emit_not_a], // INFO: P state emissions
            [uniform; 4]                                  // INFO: S state emissions
        ]
    }

    /// Run the Viterbi algorithm on the given sequence
    ///
    /// # Arguments
    ///
    /// * `sequence` - A slice of usize representing the sequence
    ///
    /// # Returns
    ///
    /// * `Vec<State>` - The most likely state sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let hmm = HMM::init(0.9, 0.99);
    /// let sequence = vec![0, 1, 2, 3, 4];
    /// let path = hmm.__viterbi(&sequence);
    ///
    /// assert_eq!(path, vec![State(0), State(1), State(1), State(1), State(1)]);
    /// ```
    fn __viterbi(&self, sequence: &[usize]) -> Vec<State> {
        let (path, _) = viterbi(&self.model, sequence);
        return path;
    }

    /// Segment the sequence using the Viterbi algorithm
    ///
    /// # Arguments
    ///
    /// * `sequence` - A slice of usize representing the sequence
    ///
    /// # Returns
    ///
    /// * `usize` - The length of the tail
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::HMM;
    ///
    /// let hmm = HMM::init(0.9, 0.99);
    /// let sequence = vec![0, 1, 2, 3, 4];
    /// let tail_len = hmm.segment(&sequence);
    ///
    /// assert_eq!(tail_len, 5);
    /// ```
    fn segment(&self, sequence: &[usize]) -> usize {
        let path = self.__viterbi(sequence);

        // INFO: iterate over path to find end of tail
        // INFO: state 0 corresponds to P
        let stop = path
            .iter()
            .rposition(|&state| state.0 == 0)
            .map_or(0, |pos| pos + 1);

        let tail = &sequence[..stop];
        return tail.len();
    }
}

/// Representation of a BAM record in isotools
///
/// This struct is used to store the information of a BAM record
/// in a handy way for isotools.
///
/// # Example
///
/// ```rust, no_run
/// use iso::Read;
///
/// let read = Read::new();
/// read.set_name("R1");
/// read.set_strand(Strand::Forward);
/// read.set_five_clip(5);
/// read.set_three_clip(3);
///
///
/// assert_eq!(read.name, "R1");
/// assert_eq!(read.strand, Strand::Forward);
/// assert_eq!(read.five_clip, 5);
/// assert_eq!(read.three_clip, 3);
/// ```
#[derive(Debug, PartialEq, Clone)]
struct Read {
    name: String,
    strand: Strand,
    five_clip: usize,
    three_clip: usize,
    has_hard_clip_three: bool,
    sequence: Sequence,
    identity: f32,
    alignment_length: usize,
    matches: usize,
    end_site: usize,
    polya_len: usize,
    polya_read_len: usize,
}

#[allow(dead_code)]
impl Read {
    /// Creates a new `Read` instance
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// assert_eq!(read.name, "");
    /// ```
    fn new() -> Self {
        Self {
            name: String::new(),
            strand: Strand::Forward,
            five_clip: 0,
            three_clip: 0,
            has_hard_clip_three: false,
            sequence: Sequence::new(b""),
            identity: 0.0,
            alignment_length: 0,
            matches: 0,
            end_site: 0,
            polya_len: 0,
            polya_read_len: 0,
        }
    }

    /// Creates a new `Read` instance from a BAM record
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let record = Record::new();
    /// let read = Read::from(&record);
    ///
    /// assert_eq!(read.name, record.name().unwrap());
    /// ```
    fn from(record: &Record) -> Self {
        let mut read = Read::new();

        read.set_name(record);
        read.get_strand_from_record(record);
        read.set_mapping_features(record);
        read.set_sequence(record);

        read
    }

    /// Get the strand from the BAM record
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let record = Record::new();
    /// let mut read = Read::from(&record);
    ///
    /// read.get_strand_from_record(&record);
    /// assert_eq!(read.strand, Strand::Forward);
    /// ```
    fn get_strand_from_record(&mut self, record: &Record) {
        let flag = record.flags();

        let strand = match flag.is_reverse_complemented() {
            true => Strand::Reverse,
            false => Strand::Forward,
        };

        self.set_strand(strand);
    }

    /// Set the strand of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_strand(Strand::Forward);
    ///
    /// assert_eq!(read.strand, Strand::Forward);
    /// ```
    fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
    }

    /// Set the polya length of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_polya_len(10);
    ///
    /// assert_eq!(read.polya_len, 10);
    /// ```
    fn set_polya_len(&mut self, polya: usize) {
        self.polya_len = polya;
    }

    /// Set the polya read length of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_polya_read_len(10);
    ///
    /// assert_eq!(read.polya_read_len, 10);
    /// ```
    fn set_polya_read_len(&mut self, polya: usize) {
        self.polya_read_len = polya;
    }

    /// Set the main mapping features of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_mapping_features(&record);
    ///
    /// assert_eq!(read.alignment_length, 100);
    /// assert_eq!(read.matches, 90);
    /// assert_eq!(read.identity, 90.0);
    /// assert_eq!(read.end_site, 200);
    /// ```
    fn set_mapping_features(&mut self, record: &Record) {
        let mut alignment_length = 0;
        let mut matches = 0;
        let mut expansion = 0;

        let mut five_clip = 0;
        let mut three_clip = 0;
        let mut hard_clip_five = false;
        let mut hard_clip_three = false;

        let mut end_site = record
            .alignment_start()
            .and_then(|start| start.ok())
            .expect(&format!(
                "ERROR: could not get alignment start for record: {:?}",
                record.name()
            ))
            .get();

        record.cigar().iter().for_each(|op| {
            let op = op.unwrap();

            match op.kind() {
                Kind::SequenceMatch => {
                    alignment_length += op.len();
                    matches += op.len();
                    expansion += op.len();
                }
                Kind::SequenceMismatch | Kind::Insertion | Kind::Deletion => {
                    alignment_length += op.len();
                    expansion += op.len();
                }
                Kind::HardClip | Kind::SoftClip => {
                    if alignment_length == 0 {
                        five_clip += op.len();

                        if op.kind() == Kind::HardClip {
                            hard_clip_five = true;
                        }
                    } else {
                        three_clip += op.len();

                        if op.kind() == Kind::HardClip {
                            hard_clip_three = true;
                        }
                    }
                }
                _ => {}
            }
        });

        match self.strand {
            Strand::Forward => {
                end_site += expansion - 1;

                self.five_clip = five_clip;
                self.three_clip = three_clip;
                self.has_hard_clip_three = hard_clip_three;
            }
            Strand::Reverse => {
                self.five_clip = three_clip;
                self.three_clip = five_clip;
                self.has_hard_clip_three = hard_clip_five;
            }
        }

        let identity = (matches as f32 / alignment_length as f32 * 1000.0).round() / 10.0;

        self.set_matches(matches);
        self.set_alignment_length(alignment_length);
        self.set_end_site(end_site);
        self.set_identity(identity);
    }

    /// Set matches of the read based on the CIGAR
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_matches(10);
    ///
    /// assert_eq!(read.matches, 10);
    /// ```
    fn set_matches(&mut self, matches: usize) {
        self.matches = matches;
    }

    /// Set the alignment length of the read based on the CIGAR
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_alignment_length(100);
    ///
    /// assert_eq!(read.alignment_length, 100);
    /// ```
    fn set_alignment_length(&mut self, alignment_length: usize) {
        self.alignment_length = alignment_length;
    }

    /// Set the end site of the read based on the CIGAR
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_end_site(200);
    ///
    /// assert_eq!(read.end_site, 200);
    /// ```
    fn set_end_site(&mut self, end_site: usize) {
        self.end_site = end_site;
    }

    /// Set the identity of the read based on the CIGAR
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_identity(90.0);
    ///
    /// assert_eq!(read.identity, 90.0);
    /// ```
    fn set_identity(&mut self, identity: f32) {
        self.identity = identity;
    }

    /// Set the sequence of the read from the BAM record
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let record = Record::new();
    /// let mut read = Read::from(&record);
    ///
    /// read.set_sequence(&record);
    /// assert_eq!(read.sequence, record.sequence());
    /// ```
    fn set_sequence(&mut self, record: &Record) {
        let sequence = Sequence::decode(record.sequence().as_ref());

        match self.strand {
            Strand::Forward => self.sequence = sequence,
            Strand::Reverse => self.sequence = sequence.reverse_complement(),
        }
    }

    /// Set the three clip of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let mut read = Read::new();
    /// read.set_three_clip(10);
    ///
    /// assert_eq!(read.three_clip, 10);
    /// ```
    fn set_three_clip(&mut self, three_clip: usize) {
        self.three_clip = three_clip;
    }

    /// Set the name of read based on the BAM record name
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let record = Record::new();
    /// let mut read = Read::from(&record);
    ///
    /// read.set_name(&record);
    /// assert_eq!(read.name, record.name().unwrap());
    /// ```
    fn set_name(&mut self, record: &Record) {
        self.name = record
            .name()
            .expect(&format!(
                "ERROR: could not get name for record: {:?}",
                record.name()
            ))
            .to_string();
    }

    /// Get the sized reverse clipped sequence of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let rev_clipped_seq = read.get_rev_clipped_seq();
    ///
    /// assert_eq!(rev_clipped_seq, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_rev_clipped_seq_sized(&self) -> Vec<usize> {
        let seq_len = self.sequence.len();
        let three_clip = self.three_clip;

        let start = if three_clip > seq_len {
            0
        } else {
            seq_len - three_clip
        };

        self.sequence.reverse_encode(start, seq_len)
    }

    /// Get the reverse clipped sequence of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let rev_clipped_seq = read.get_rev_clipped_seq();
    ///
    /// assert_eq!(rev_clipped_seq, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_rev_clipped_seq(&self) -> Vec<u8> {
        let seq_len = self.sequence.len();
        let three_clip = self.three_clip;

        let start = if three_clip > seq_len {
            0
        } else {
            seq_len - three_clip
        };

        self.sequence.reverse_encode_u8(start, seq_len)
    }

    /// Get the clipped sequence of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let clipped_seq = read.get_clipped_seq();
    ///
    /// assert_eq!(clipped_seq, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_clipped_seq(&self) -> &[u8] {
        let seq_len = self.sequence.len();
        let three_clip = self.three_clip;

        assert!(
            three_clip < seq_len,
            "ERROR: suffix length is greater than sequence length -> {} : {}!",
            three_clip,
            seq_len
        );

        self.sequence.slice_as_bytes(seq_len - three_clip, seq_len)
    }

    /// Get the suffix of the read
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let suffix = read.get_seq_suffix(5);
    ///
    /// assert_eq!(suffix, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_seq_suffix(&self, suffix: usize) -> &[u8] {
        let seq_len = self.sequence.len();

        assert!(
            suffix < seq_len,
            "ERROR: suffix length is greater than sequence length -> {} : {}!",
            suffix,
            seq_len
        );

        return self.sequence.slice_as_bytes(seq_len - suffix, seq_len);
    }

    /// Get the reverse suffix of the read sized
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.set_three_clip(10);
    /// let rev_suffix = read.get_rev_seq_suffix_sized(5);
    ///
    /// assert_eq!(rev_suffix, vec![1, 2, 3, 4, 5]);
    /// ```
    fn get_rev_seq_suffix_sized(&self, suffix: usize) -> Vec<usize> {
        let seq_len = self.sequence.len();

        let start = if suffix > seq_len {
            0
        } else {
            seq_len - suffix
        };

        return self.sequence.reverse_encode(start, seq_len);
    }

    /// Tag read following the format:
    ///     > R[number]_[chr]::FC5:TC24:PA45:PR65:IY98
    /// where:
    ///     R[number] = read number
    ///     [chr] = chromosome name
    ///     FC = five_clip
    ///     TC = three_clip
    ///     PA = polya_len
    ///     PR = polya_read_len
    ///     IY = identity
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Read;
    ///
    /// let read = Read::new();
    /// read.tag_read();
    ///
    /// assert_eq!(read.name, "R1_chr1::FC5:TC24:PA45:PR65:IY98");
    /// ```
    fn tag_read(&self, index: u64, chr: &String) -> String {
        format!(
            "R{}_{chr}::FC{}:TC{}:PA{}:PR{}:IY{}",
            index,
            self.five_clip,
            self.three_clip,
            self.polya_len,
            self.polya_read_len,
            (self.identity * 10.0) as u64
        )
    }
}
