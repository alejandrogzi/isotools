use anyhow::Result;
use dashmap::{DashMap, DashSet};
use hashbrown::{HashMap, HashSet};
use packbed::{packbed, record::IntronPredStats, BedPackage, IntronBucket, PackMode};
use rayon::prelude::*;

use config::{
    get_progress_bar, write_objs, Sequence, SharedSpliceMap, Strand, SupportType,
    INTRON_CLASSIFICATION, INTRON_FREQUENCY_RECOVERY_THRESHOLD, MAX_ENT_SCORE_RECOVERY_THRESHOLD,
    SCALE, SPLICE_AI_SCORE_RECOVERY_THRESHOLD,
};

use crate::cli::IntronArgs as Args;
use crate::utils::{
    calculate_acceptor_score, create_splice_map, get_sequences, get_splice_scores,
    load_scan_scores, unpack_blacklist, SpliceScoreMap,
};

const WINDOW_SIZE: usize = 8;
const RT_REPEAT: usize = 12;
const MISMATCHES: u32 = 1;
const MAX_ENT_DONOR_MIN_SIZE: usize = 9;
const NAG_PATTERNS: [&str; 3] = ["CAG", "TAG", "AAG"];
const WIGGLE_SWITCH: [usize; 2] = [2, 4];

type ScanScores = Option<(SpliceScoreMap, SpliceScoreMap)>;
pub type Genome = DashMap<String, Vec<u8>>;

pub fn classify_introns(args: Args) -> Result<()> {
    let isoseqs = packbed(
        args.iso,
        args.toga,
        config::OverlapType::Exon,
        PackMode::Intron,
    )?;

    let (splice_plus, splice_minus) = get_splice_scores(args.spliceai);
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let genome = if let Some(twobit) = args.twobit {
        get_sequences(twobit)
    } else {
        None
    };

    // if args.scan is true, we need to load both databases
    let scan_scores = if args.scan { load_scan_scores() } else { None };

    let pb = get_progress_bar(isoseqs.len() as u64, "Processing...");

    let accumulator = ParallelAccumulator::default();

    isoseqs.into_par_iter().for_each(|bucket| {
        let chr = bucket.0;
        let components = bucket.1;

        let binding = HashSet::new();
        let banned = blacklist.get(&chr).unwrap_or(&binding);
        let splice_map = create_splice_map(&chr, &splice_plus, &splice_minus);

        distribute(
            components,
            banned,
            &splice_map,
            &scan_scores,
            &genome,
            &accumulator,
            args.nag,
        );

        pb.inc(1);
    });

    pb.finish_and_clear();
    write_objs(&accumulator.introns, INTRON_CLASSIFICATION);

    Ok(())
}

struct ParallelAccumulator {
    introns: DashSet<String>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            introns: DashSet::new(),
        }
    }
}

#[inline(always)]
fn distribute(
    components: Vec<Box<dyn BedPackage>>,
    banned: &HashSet<(u64, u64)>,
    splice_map: &(SharedSpliceMap, SharedSpliceMap),
    scan_scores: &ScanScores,
    genome: &Option<Genome>,
    accumulator: &ParallelAccumulator,
    nag: bool,
) {
    components.into_par_iter().for_each(|mut comp| {
        let comp = comp
            .as_any_mut()
            .downcast_mut::<IntronBucket>()
            .expect("ERROR: Could not downcast to IntronPred!");

        let info = process_component(comp, banned, splice_map, scan_scores, genome, nag);

        info.into_iter().for_each(|(_, intron_descriptor)| {
            if !intron_descriptor.is_empty() {
                accumulator.introns.insert(intron_descriptor);
            }
        });
    });
}

#[inline(always)]
fn process_component(
    component: &mut IntronBucket,
    banned: &HashSet<(u64, u64)>,
    splice_map: &(SharedSpliceMap, SharedSpliceMap),
    scan_scores: &ScanScores,
    genome: &Option<Genome>,
    nag: bool,
) -> HashMap<(u64, u64), String> {
    let chr = component.chrom.clone();
    let strand = component.strand.clone();

    let mut acc = if nag {
        HashMap::new()
    } else {
        HashMap::with_capacity(component.introns.len())
    };

    // INFO: getting stranded spliceAi scores
    let splice_scores = match component.strand {
        Strand::Forward => Some(&splice_map.0),
        Strand::Reverse => Some(&splice_map.1),
    };

    for (intron, descriptor) in component.introns.iter_mut() {
        let intron_start = intron.0 as u64;
        let intron_end = intron.1 as u64;

        if banned.contains(&(intron_start, intron_end)) {
            continue;
        }

        get_sj_context(intron, descriptor, &strand, &chr, genome, scan_scores);
        get_sj_ai_scores(intron, descriptor, splice_scores, &strand);

        if nag && descriptor.acceptor_sequence == "AG" && descriptor.is_toga_supported {
            scan_nag_repeats(
                intron,
                descriptor,
                &strand,
                &chr,
                genome,
                scan_scores,
                splice_scores,
                &mut acc,
            );
        }

        // WARN: breaking point to convert coords!
        let (intron_start, intron_end) = match component.strand {
            Strand::Forward => (intron_start, intron_end),
            Strand::Reverse => ((SCALE - intron_end), (SCALE - intron_start)),
        };

        // INFO: not including NAG-derived introns bc they are already Splicing
        // WARN: strange TOGA supported RT introns -> s14	9965456	9968743
        if descriptor.is_rt_intron {
            if descriptor.is_toga_supported {
                descriptor.support = SupportType::Splicing;
            } else {
                descriptor.support = SupportType::RT;
            }
        } else {
            if descriptor.is_toga_supported
                || (descriptor.splice_ai_donor >= SPLICE_AI_SCORE_RECOVERY_THRESHOLD
                    && descriptor.splice_ai_acceptor >= SPLICE_AI_SCORE_RECOVERY_THRESHOLD)
            {
                descriptor.support = SupportType::Splicing;
            } else if (descriptor.seen as f64 / descriptor.spanned as f64) // WARN: how do we deal with 1/2 (50%) cases?
                >= INTRON_FREQUENCY_RECOVERY_THRESHOLD
            {
                descriptor.support = SupportType::Splicing;
            } else if (descriptor.max_ent_donor >= MAX_ENT_SCORE_RECOVERY_THRESHOLD
                && descriptor.max_ent_acceptor >= MAX_ENT_SCORE_RECOVERY_THRESHOLD)
                && (descriptor.splice_ai_donor > 0.0 && descriptor.splice_ai_acceptor > 0.0)
            {
                // INFO: new branch for MaxEnt only -> not trusting it alone
                // INFO: here we test if maxEnt is significant + if there is
                // INFO: spliceAi signal [> 0.0]
                descriptor.support = SupportType::Splicing;
            } else {
                descriptor.support = SupportType::Unclear;
            }
        }

        acc.insert(
            (intron_start, intron_end),
            descriptor.fmt(&chr, &strand, intron_start, intron_end),
        );
    }

    acc
}

fn get_sj_context(
    intron: &(u64, u64),
    descriptor: &mut IntronPredStats,
    strand: &Strand,
    chr: &String,
    genome: &Option<Genome>,
    scan_scores: &ScanScores,
) {
    let intron_start = intron.0 as u64;
    let intron_end = intron.1 as u64;

    if let Some(genome) = genome {
        match strand {
            Strand::Forward => {
                // INFO: For TOGA nag we need +3 upstream and +3 downstream
                // INFO: For RT repeats we need 12-mers (11 exon + 1 intron) at 5'
                // INFO: and 11 intron + 1 exon at 3'

                let donor_context = Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read donor context!")
                        [intron_start as usize - 12..intron_start as usize + 5]
                        .as_ref(),
                );
                let donor_seq = donor_context.slice(11, 13);
                let donor_rt_context = donor_context.slice(0, 12);

                let acceptor_context = Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read acceptor context!")
                        [intron_end as usize - 19..intron_end as usize + 4]
                        .as_ref(),
                );
                let acceptor_seq = acceptor_context.slice(18, 20);
                let acceptor_rt_context = acceptor_context.slice(9, 21);

                descriptor.donor_sequence = donor_seq;
                descriptor.donor_context = donor_context;

                descriptor.acceptor_sequence = acceptor_seq;
                descriptor.acceptor_context = acceptor_context;

                descriptor.donor_rt_context = donor_rt_context;
                descriptor.acceptor_rt_context = acceptor_rt_context;
            }
            Strand::Reverse => {
                let donor_context = Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read donor context!")
                        [(SCALE - intron_start) as usize - 5..(SCALE - intron_start) as usize + 12]
                        .as_ref(),
                )
                .reverse_complement();
                let donor_seq = donor_context.slice(11, 13);
                let donor_rt_context = donor_context.slice(0, 12);

                let acceptor_context = Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read acceptor context!")
                        [(SCALE - intron_end) as usize - 4..(SCALE - intron_end) as usize + 19]
                        .as_ref(),
                )
                .reverse_complement();
                let acceptor_seq = acceptor_context.slice(18, 20);
                let acceptor_rt_context = acceptor_context.slice(9, 21);

                descriptor.donor_sequence = donor_seq;
                descriptor.donor_context = donor_context;

                descriptor.acceptor_sequence = acceptor_seq;
                descriptor.acceptor_context = acceptor_context;

                descriptor.donor_rt_context = donor_rt_context;
                descriptor.acceptor_rt_context = acceptor_rt_context;
            }
        }

        get_sj_max_entropy(descriptor, scan_scores);
        scan_rt_repeats(descriptor);
    }
}

fn get_sj_max_entropy(descriptor: &mut IntronPredStats, scan_scores: &ScanScores) {
    if let Some(scan_scores) = scan_scores {
        let (donor_score_map, acceptor_score_map): &(SpliceScoreMap, SpliceScoreMap) = scan_scores;

        let donor_max_ent_context = &descriptor.donor_context.slice(8, 17);

        if donor_max_ent_context.len() != MAX_ENT_DONOR_MIN_SIZE {
            eprintln!("ERROR: Donor context is not 9 bases long!");
        }

        let donor_score = donor_score_map
            .get(donor_max_ent_context)
            .and_then(|r| r.get(0))
            .unwrap_or(&0.0);

        descriptor.max_ent_donor = *donor_score as f32;

        let acceptor_score =
            calculate_acceptor_score(&descriptor.acceptor_context, acceptor_score_map);

        descriptor.max_ent_acceptor = acceptor_score as f32;
    }
}

fn get_sj_ai_scores(
    intron: &(u64, u64),
    descriptor: &mut IntronPredStats,
    splice_scores: Option<&SharedSpliceMap>,
    strand: &Strand,
) {
    let intron_start = intron.0 as u64;
    let intron_end = intron.1 as u64;

    if let Some(splice_scores) = splice_scores {
        if let Some(donor_score_map) = splice_scores.0.as_ref() {
            let acceptor_score_map = splice_scores
                .1
                .as_ref()
                .expect("ERROR: Acceptor score map is None, this is a bug!");

            let (intron_donor, intron_acceptor) = match strand {
                // donor [-1 to match bigtools coords]
                Strand::Forward => (intron_start as usize - 1, intron_end as usize),
                // acceptor [-1 to match bigtools coords]
                Strand::Reverse => (
                    (SCALE - intron_start) as usize,
                    (SCALE - intron_end) as usize - 1,
                ),
            };

            let (donor_score, acceptor_score) = (
                donor_score_map
                    .get(&intron_donor)
                    .map(|r| *r)
                    .unwrap_or(0.0),
                acceptor_score_map
                    .get(&intron_acceptor)
                    .map(|r| *r)
                    .unwrap_or(0.0),
            );

            descriptor.splice_ai_donor = donor_score;
            descriptor.splice_ai_acceptor = acceptor_score;
        }
    }
}

fn process_nag_pattern(
    base_intron: &(u64, u64),
    offset: i64,
    strand: &Strand,
    chr: &String,
    genome: &Option<Genome>,
    scan_scores: &ScanScores,
    splice_scores: Option<&SharedSpliceMap>,
) -> String {
    let intron = (base_intron.0 as u64, (base_intron.1 as i64 + offset) as u64);

    let mut new_descriptor = IntronPredStats::new();
    get_sj_context(
        &intron,
        &mut new_descriptor,
        strand,
        chr,
        genome,
        scan_scores,
    );
    get_sj_ai_scores(&intron, &mut new_descriptor, splice_scores, strand);
    new_descriptor.is_nag_intron = true;

    let scaled_intron = match strand {
        Strand::Forward => intron,
        Strand::Reverse => ((SCALE - intron.1), (SCALE - intron.0)),
    };

    // WARN: NAG-derived introns will be always SupportType::Splicing
    new_descriptor.support = SupportType::Splicing;
    new_descriptor.fmt(chr, strand, scaled_intron.0, scaled_intron.1)
}

fn scan_nag_repeats(
    intron: &(u64, u64),
    descriptor: &mut IntronPredStats,
    strand: &Strand,
    chr: &String,
    genome: &Option<Genome>,
    scan_scores: &ScanScores,
    splice_scores: Option<&SharedSpliceMap>,
    // acc: &mut Vec<String>,
    acc: &mut HashMap<(u64, u64), String>,
) {
    let pre_acceptor = &descriptor.acceptor_context.slice(14, 17);
    let post_acceptor = &descriptor.acceptor_context.slice(20, 23);

    if NAG_PATTERNS.contains(&pre_acceptor.as_str()) {
        let rs = process_nag_pattern(intron, -3, strand, chr, genome, scan_scores, splice_scores);
        acc.insert((intron.0, intron.1 - 3), rs);
    }

    if NAG_PATTERNS.contains(&post_acceptor.as_str()) {
        let rs = process_nag_pattern(intron, 3, strand, chr, genome, scan_scores, splice_scores);
        acc.insert((intron.0, intron.1 + 3), rs);
    }
}

fn scan_rt_repeats(descriptor: &mut IntronPredStats) {
    unsafe {
        scan_sequence(descriptor);
    }
}

#[inline(always)]
unsafe fn scan_sequence(descriptor: &mut IntronPredStats) {
    let donor = &descriptor.donor_rt_context.as_bytes();
    let acceptor = &descriptor.acceptor_rt_context.as_bytes();

    if donor.len() != RT_REPEAT || acceptor.len() != RT_REPEAT {
        log::error!("ERROR: RT context is not 12 bases long, this is a bug!");
        std::process::exit(1);
    }

    #[inline(always)]
    fn base_to_bits(base: u8) -> u8 {
        match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 0,
            _ => panic!("Invalid base"),
        }
    }

    // INFO: convert 8-mer to u16 (2 bits per base = 16 bits total)
    #[inline(always)]
    fn window_to_int(window: &[u8]) -> u16 {
        let mut value: u16 = 0;
        for &base in window {
            value = (value << 2) | (base_to_bits(base) as u16);
        }
        value
    }

    #[inline(always)]
    fn hamming_distance(a: u16, b: u16) -> u32 {
        (a ^ b).count_ones()
    }

    if donor.len() < WINDOW_SIZE || acceptor.len() < WINDOW_SIZE {
        return;
    }

    // INFO: build a set of all 8-mers in seq1
    let mut kmers = HashSet::with_capacity(donor.len() - WINDOW_SIZE + 1);
    for idx in 0..=donor.len() - WINDOW_SIZE {
        let kmer = window_to_int(&donor[idx..idx + WINDOW_SIZE]);
        kmers.insert((kmer, idx));
    }

    // INFO: compare k-mers from seq2 against the set
    for idx in 0..=acceptor.len() - WINDOW_SIZE {
        let kmer2 = window_to_int(&acceptor[idx..idx + WINDOW_SIZE]);

        for &(kmer1, _) in &kmers {
            if hamming_distance(kmer1, kmer2) <= MISMATCHES {
                // INFO: returning bool and skipping early
                // repeats.push((i, j));
                // dbg!("Match found: seq1[{}..{}], seq2[{}..{}]", _, _ + 8, idx, idx + 8);
                descriptor.is_rt_intron = true;
                return;
            }
        }
    }

    descriptor.is_rt_intron = false;
}

#[allow(unused)]
fn wiggle_splice_sites(
    acc: &mut HashMap<(u64, u64), String>,
    intron: &(u64, u64),
    ref_introns: &HashSet<(u64, u64)>,
    strand: &Strand,
) {
    fn process_intron(
        acc: &mut HashMap<(u64, u64), String>,
        ref_introns: &HashSet<(u64, u64)>,
        intron_start: u64,
        intron_end: u64,
        strand: &Strand,
    ) {
        // INFO: look if swithed intron is inside ref_introns
        if ref_introns.contains(&(intron_start, intron_end)) {
            let (intron_start, intron_end) = match strand {
                Strand::Forward => (intron_start, intron_end),
                Strand::Reverse => ((SCALE - intron_end), (SCALE - intron_start)),
            };

            // INFO: if its already classified, change support and toga_support
            if acc.contains_key(&(intron_start, intron_end)) {
                let mut fields = acc
                    .get(&(intron_start, intron_end))
                    .unwrap()
                    .split('\t')
                    .collect::<Vec<_>>();

                fields[15] = "true";
                fields[21] = "SPLICED";

                acc.insert((intron_start, intron_end), fields.join("\t"));
            } else {
                // WARN: if its not classified, we need to insert it, leaving the rest of the fields empty
                // WARN: the String inserted here should be unreacheable!
                acc.insert((intron_start, intron_end), String::new());
            }
        }
    }

    for wiggle in WIGGLE_SWITCH {
        // INFO: backward
        process_intron(
            acc,
            ref_introns,
            intron.0 as u64 - wiggle as u64,
            intron.1 as u64 - wiggle as u64,
            strand,
        );

        // INFO: forward
        process_intron(
            acc,
            ref_introns,
            intron.0 as u64 + wiggle as u64,
            intron.1 as u64 + wiggle as u64,
            strand,
        );
    }
}
