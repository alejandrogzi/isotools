use anyhow::Result;
use dashmap::DashMap;
use dashmap::DashSet;
use hashbrown::HashSet;
use packbed::{packbed, BedPackage, IntronPred, PackMode};
use rayon::prelude::*;

use config::{
    get_progress_bar, write_objs, Sequence, SharedSpliceMap, Strand, INTRON_CLASSIFICATION,
    OVERLAP_CDS, OVERLAP_EXON, SCALE,
};

use crate::cli::IntronArgs as Args;
use crate::utils::{
    calculate_acceptor_score, create_splice_map, get_sequences, get_splice_scores,
    load_scan_scores, unpack_blacklist, SpliceScoreMap,
};

type ScanScores = Option<(SpliceScoreMap, SpliceScoreMap)>;
type Genome = DashMap<String, Vec<u8>>;

pub fn classify_introns(args: Args) -> Result<()> {
    let isoseqs = packbed(
        args.iso,
        args.toga,
        OVERLAP_CDS,
        OVERLAP_EXON,
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
            .downcast_mut::<IntronPred>()
            .expect("ERROR: Could not downcast to IntronPred!");

        let info = process_component(comp, banned, splice_map, scan_scores, genome, nag);

        info.into_iter().for_each(|intron| {
            accumulator.introns.insert(intron);
        });
    });
}

// working with a single component (IntronPred)
// WARN: This fn should be granularized!
#[inline(always)]
fn process_component(
    component: &mut IntronPred,
    banned: &HashSet<(u64, u64)>,
    splice_map: &(SharedSpliceMap, SharedSpliceMap),
    scan_scores: &ScanScores,
    genome: &Option<Genome>,
    nag: bool,
) -> Vec<String> {
    let mut acc = Vec::with_capacity(component.introns.len());

    // get SpliceAi scores
    let splice_scores = match component.strand {
        Strand::Forward => Some(&splice_map.0),
        Strand::Reverse => Some(&splice_map.1),
    };

    // QUESTION: should we parallelize this?
    for (intron, descriptor) in component.introns.iter_mut() {
        let intron_start = intron.0 as u64;
        let intron_end = intron.1 as u64;

        if banned.contains(&(intron_start, intron_end)) {
            continue;
        }

        // get sequences
        if let Some(genome) = genome {
            // depends on the strand
            match component.strand {
                Strand::Forward => {
                    let donor_context = Sequence::new(
                        genome
                            .get(component.chrom.as_str())
                            .expect("ERROR: Could not read donor context!")
                            [intron_start as usize - 4..intron_start as usize + 5]
                            .as_ref(),
                    );
                    let donor_seq = donor_context.slice(3, 5);

                    let acceptor_context = Sequence::new(
                        genome
                            .get(component.chrom.as_str())
                            .expect("ERROR: Could not read acceptor context!")
                            [intron_end as usize - 19..intron_end as usize + 4]
                            .as_ref(),
                    );
                    let acceptor_seq = acceptor_context.slice(18, 20);

                    descriptor.donor_sequence = donor_seq;
                    descriptor.acceptor_sequence = acceptor_seq;
                    descriptor.donor_context = donor_context;
                    descriptor.acceptor_context = acceptor_context;
                }
                Strand::Reverse => {
                    let donor_context = Sequence::new(
                        genome
                            .get(component.chrom.as_str())
                            .expect("ERROR: Could not read donor context!")
                            [(SCALE - intron_start) as usize - 5
                                ..(SCALE - intron_start) as usize + 4]
                            .as_ref(),
                    )
                    .reverse_complement();
                    let donor_seq = donor_context.slice(3, 5);

                    let acceptor_context = Sequence::new(
                        genome
                            .get(component.chrom.as_str())
                            .expect("ERROR: Could not read acceptor context!")
                            [(SCALE - intron_end) as usize - 2..(SCALE - intron_end) as usize + 21]
                            .as_ref(),
                    )
                    .reverse_complement();
                    let acceptor_seq = acceptor_context.slice(18, 20);

                    descriptor.donor_sequence = donor_seq;
                    descriptor.acceptor_sequence = acceptor_seq;
                    descriptor.donor_context = donor_context;
                    descriptor.acceptor_context = acceptor_context;
                }
            }

            // get MaxEnt scores
            if let Some(scan_scores) = scan_scores {
                let (donor_score_map, acceptor_score_map): &(SpliceScoreMap, SpliceScoreMap) =
                    scan_scores;

                let donor_score = donor_score_map
                    .get(&descriptor.donor_context)
                    .and_then(|r| r.get(0))
                    .unwrap_or(&0.0);

                descriptor.max_ent_donor = *donor_score as f32;

                let acceptor_score =
                    calculate_acceptor_score(&descriptor.acceptor_context, acceptor_score_map);

                descriptor.max_ent_acceptor = acceptor_score as f32;
            }
        }

        // get spliceAi scores
        if let Some(splice_scores) = splice_scores {
            if let Some(donor_score_map) = splice_scores.0.as_ref() {
                let acceptor_score_map = splice_scores
                    .1
                    .as_ref()
                    .expect("ERROR: Acceptor score map is None, this is a bug!");

                let (intron_donor, intron_acceptor) = match component.strand {
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

        let (intron_start, intron_end) = match component.strand {
            Strand::Forward => (intron_start, intron_end),
            Strand::Reverse => ((SCALE - intron_end), (SCALE - intron_start)),
        };

        acc.push(descriptor.fmt(
            component.chrom.clone(),
            component.strand.clone(),
            intron_start,
            intron_end,
        ));
    }

    acc
}
