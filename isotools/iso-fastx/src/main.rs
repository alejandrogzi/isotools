// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! iso-fastx: classify (`inspect`) and strand-orient (`orient`) long-read FASTA/FASTQ files
//! before alignment.
//!
//! `inspect` samples the first records of a file and decides what the pipeline must do with
//! it (lima, orient, align as is, or stop); `orient` flips the reads that carry a 5' polyT
//! instead of a 3' polyA. Primers come from iso-adapter's database.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use flate2::{write::GzEncoder, Compression};
use iso_adapter::adapters::ADAPTER_DB;
use iso_adapter::cli::{DEFAULT_MAX_EDIT_DIST, DEFAULT_MIN_CLIP_LEN};
use iso_adapter::detector::AdapterDb;
use needletail::errors::ParseErrorKind;
use needletail::parser::{write_fasta, write_fastq, LineEnding, SequenceRecord};
use needletail::{parse_fastx_file, Sequence};

/// Nucleotides at each read end scanned for a polyA / polyT tail.
const TAIL_WINDOW: usize = 60;
/// Consecutive A / T that make a tail.
const TAIL_RUN: usize = 15;
/// Nucleotides at each read end scanned for primers.
const PRIMER_WINDOW: usize = 150;
/// Library kits, weakest evidence first: a read takes the strongest kit among its hits.
const KITS: [&str; 3] = ["express", "isoseqx", "custom"];

#[derive(Parser)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = "Inspect and orient long-read FASTA/FASTQ files"
)]
enum Cli {
    /// Classify a file (state, kit, tails, length, quality, name style) into <prefix>.inspect.tsv
    Inspect {
        /// FASTA/FASTQ file, optionally compressed
        #[arg(long, value_name = "PATH")]
        fastx: PathBuf,
        /// Prefix of the output report
        #[arg(long)]
        prefix: String,
        /// Records to sample from the start of the file
        #[arg(long, default_value_t = 5000)]
        reads: usize,
        /// Extra primers (FASTA), matched as `user:<id>`
        #[arg(long, value_name = "PATH")]
        primers: Option<PathBuf>,
    },
    /// Reverse-complement the polyT reads of a file and report in <prefix>.orient.tsv
    Orient {
        /// FASTA/FASTQ file, optionally compressed
        #[arg(long, value_name = "PATH")]
        fastx: PathBuf,
        /// Output file, same format as the input; gzipped when it ends in .gz
        #[arg(long, value_name = "PATH")]
        output: PathBuf,
        /// Prefix of the output report
        #[arg(long)]
        prefix: String,
    },
}

fn main() {
    let result = match Cli::parse() {
        Cli::Inspect {
            fastx,
            prefix,
            reads,
            primers,
        } => inspect(&fastx, &prefix, reads, primers.as_deref()),
        Cli::Orient {
            fastx,
            output,
            prefix,
        } => orient(&fastx, &output, &prefix),
    };

    if let Err(err) = result {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}

/// Streams the records of `path` (compression detected) into `f`, which returns `false` to
/// stop. A file without records, a zero-byte one included, is not an error.
fn each(path: &Path, mut f: impl FnMut(&SequenceRecord) -> Result<bool>) -> Result<()> {
    let mut reader = match parse_fastx_file(path) {
        Ok(reader) => reader,
        Err(e) if e.kind == ParseErrorKind::EmptyFile => return Ok(()),
        Err(e) => return Err(e).with_context(|| format!("failed to open {}", path.display())),
    };

    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("failed to parse {}", path.display()))?;
        if !f(&record)? {
            break;
        }
    }

    Ok(())
}

/// `(polyA in the last 60 nt, polyT in the first 60 nt)`, case-insensitive.
fn tails(seq: &[u8]) -> (bool, bool) {
    let run = |window: &[u8], base: u8| {
        window
            .split(|b| !b.eq_ignore_ascii_case(&base))
            .any(|r| r.len() >= TAIL_RUN)
    };

    (
        run(&seq[seq.len().saturating_sub(TAIL_WINDOW)..], b'A'),
        run(&seq[..seq.len().min(TAIL_WINDOW)], b'T'),
    )
}

/// Matcher over the Iso-Seq / Clontech primers of iso-adapter's database plus the user's.
///
/// Homopolymers, SMRTbell, ONT and Illumina entries are left out on purpose: their exact and
/// fuzzy hits (a polyA tail sits right before the 3' primer) would shadow the primer hit.
fn primer_db(user: Option<&Path>) -> Result<AdapterDb> {
    let mut entries: Vec<(Vec<u8>, &'static str)> = ADAPTER_DB
        .iter()
        .filter(|(_, label)| label.starts_with("pacbio:isoseq") || label.starts_with("clontech:"))
        .map(|&(seq, label)| (seq.to_vec(), label))
        .collect();

    if let Some(path) = user {
        each(path, |rec| {
            let id = rec.id().split(|b| b.is_ascii_whitespace()).next();
            let label = format!("user:{}", String::from_utf8_lossy(id.unwrap_or_default()));
            // ponytail: AdapterDb labels are &'static str, so the few user labels are leaked
            entries.push((
                rec.seq().to_ascii_uppercase(),
                Box::leak(label.into_boxed_str()),
            ));
            Ok(true)
        })?;
    }

    AdapterDb::from_entries(
        entries.iter().map(|(seq, label)| (seq.as_slice(), *label)),
        DEFAULT_MIN_CLIP_LEN,
        DEFAULT_MAX_EDIT_DIST,
    )
    .context("failed to build the primer matcher")
}

/// Strongest kit (index into `KITS`) among the primer hits in the first and last 150 nt.
fn primer_kit(db: &AdapterDb, seq: &[u8]) -> Option<usize> {
    [
        &seq[..seq.len().min(PRIMER_WINDOW)],
        &seq[seq.len().saturating_sub(PRIMER_WINDOW)..],
    ]
    .iter()
    .filter_map(|window| db.match_adapter(&window.to_ascii_uppercase()))
    .map(|hit| match hit.label {
        l if l.starts_with("user:") => 2,
        l if l.starts_with("pacbio:isoseqx") => 1,
        _ => 0,
    })
    .max()
}

/// How a read header looks: where the file probably comes from.
fn name_style(header: &[u8]) -> &'static str {
    // the remainder of `s` after at least one leading digit
    fn num(s: &str) -> Option<&str> {
        let rest = s.trim_start_matches(|c: char| c.is_ascii_digit());
        (rest.len() < s.len()).then_some(rest)
    }

    let header = String::from_utf8_lossy(header);
    let name = header.split_whitespace().next().unwrap_or_default();
    // <prefix><digits>.<digits>
    let dotted = |prefix: &str| {
        header
            .strip_prefix(prefix)
            .and_then(num)
            .and_then(|rest| rest.strip_prefix('.'))
            .and_then(num)
            .is_some()
    };
    // .../<zmw>/<start>_<end>
    let subread = {
        let mut parts = name.rsplitn(3, '/');
        match (parts.next(), parts.next(), parts.next()) {
            (Some(range), Some(zmw), Some(_)) => {
                num(zmw) == Some("")
                    && range
                        .split_once('_')
                        .is_some_and(|(a, b)| num(a) == Some("") && num(b) == Some(""))
            }
            _ => false,
        }
    };

    if header.strip_prefix("transcript/").and_then(num).is_some()
        || dotted("PB.")
        || header.contains("full_length_coverage=")
    {
        "clustered"
    } else if subread && !header.contains("/ccs") {
        "pacbio_subread"
    } else if header.contains("/ccs") {
        "pacbio"
    } else if ["SRR", "ERR", "DRR"].iter().any(|prefix| dotted(prefix)) {
        "sra"
    } else {
        "other"
    }
}

/// Samples the first `reads` records and writes `<prefix>.inspect.tsv`.
fn inspect(fastx: &Path, prefix: &str, reads: usize, primers: Option<&Path>) -> Result<()> {
    let db = primer_db(primers)?;

    let (mut n, mut primed, mut polya, mut polyt, mut len_sum) = (0, 0, 0, 0, 0);
    let (mut qv_sum, mut bases, mut qv_first, mut qv_constant) = (0u64, 0usize, None, true);
    let (mut kits, mut styles) = ([0usize; 3], BTreeMap::new());

    each(fastx, |rec| {
        let seq = rec.seq();
        let (a, t) = tails(&seq);
        n += 1;
        len_sum += seq.len();
        polya += usize::from(a);
        polyt += usize::from(t);

        if let Some(kit) = primer_kit(&db, &seq) {
            primed += 1;
            kits[kit] += 1;
        }

        *styles.entry(name_style(rec.id())).or_insert(0usize) += 1;

        if let Some(qual) = rec.qual() {
            bases += qual.len();
            qv_sum += qual
                .iter()
                .map(|&c| u64::from(c.saturating_sub(33)))
                .sum::<u64>();
            qv_constant = qv_constant && qual.iter().all(|&c| c == *qv_first.get_or_insert(c));
        }

        Ok(n < reads)
    })?;

    let rate = |x: usize| if n == 0 { 0.0 } else { x as f64 / n as f64 };
    let (primer_rate, polya3, polyt5, mean_len) =
        (rate(primed), rate(polya), rate(polyt), rate(len_sum));
    // only FASTQ has quality bases
    let mean_qv = (bases > 0).then(|| qv_sum as f64 / bases as f64);
    let style = styles
        .iter()
        .max_by_key(|&(_, &count)| count)
        .map_or("other", |(style, _)| style);
    let kit = kits
        .iter()
        .enumerate()
        .max_by_key(|&(_, &count)| count)
        .filter(|&(_, &count)| count > 0)
        .map_or("none", |(i, _)| KITS[i]);

    // constant QVs alone are not a subreads signal: SRA Lite sets every base to Q30
    let state = if n == 0 {
        "empty"
    } else if style == "clustered" {
        "clustered"
    } else if style == "pacbio_subread" || mean_qv.is_some_and(|qv| qv < 15.0) {
        "subreads"
    } else if primer_rate >= 0.30 {
        "ccs"
    } else if primer_rate >= 0.02 || mean_len > 8000.0 {
        "ambiguous"
    } else if polyt5 > 0.05 {
        "mixed"
    } else if polya3 >= 0.10 {
        "fl"
    } else {
        "flnc"
    };

    let note = match state {
        "ambiguous" => [
            (primer_rate >= 0.02).then_some("primer_rate>=0.02_partial_primer_signal"),
            (mean_len > 8000.0).then_some("mean_len>8000_possible_unsegmented_arrays_need_skera"),
        ]
        .into_iter()
        .flatten()
        .collect::<Vec<_>>()
        .join(";"),
        "mixed" if polyt5 >= 0.3 => "balanced_tails_unknown_primers?".to_string(),
        _ => "-".to_string(),
    };

    let na = || "NA".to_string();
    let qv = mean_qv.map_or_else(na, |qv| format!("{qv:.2}"));
    let constant = mean_qv.map_or_else(na, |_| qv_constant.to_string());
    let file = file_name(fastx);
    let path = format!("{prefix}.inspect.tsv");

    std::fs::write(
        &path,
        format!(
            "file\treads\tstate\tkit\tprimer_rate\tpolya3\tpolyt5\tmean_len\tmean_qv\tqv_constant\tname_style\tnote\n\
             {file}\t{n}\t{state}\t{kit}\t{primer_rate:.4}\t{polya3:.4}\t{polyt5:.4}\t{mean_len:.1}\t{qv}\t{constant}\t{style}\t{note}\n"
        ),
    )
    .with_context(|| format!("failed to write {path}"))
}

/// Writes `fastx` to `output` with every polyT-only read reverse-complemented and
/// `<prefix>.orient.tsv`. Reads are never dropped or renamed.
fn orient(fastx: &Path, output: &Path, prefix: &str) -> Result<()> {
    let file =
        File::create(output).with_context(|| format!("failed to create {}", output.display()))?;
    let mut file = BufWriter::new(file);

    let counts = if output.extension().is_some_and(|ext| ext == "gz") {
        let mut gz = GzEncoder::new(file, Compression::default());
        let counts = flip(fastx, &mut gz)?;
        file = gz.finish()?;
        counts
    } else {
        flip(fastx, &mut file)?
    };
    file.flush()?;

    let [kept, flipped, ambiguous, no_signal] = counts;
    let path = format!("{prefix}.orient.tsv");

    std::fs::write(
        &path,
        format!(
            "file\treads\tkept\tflipped\tambiguous\tno_signal\n{}\t{}\t{kept}\t{flipped}\t{ambiguous}\t{no_signal}\n",
            file_name(fastx),
            counts.iter().sum::<usize>()
        ),
    )
    .with_context(|| format!("failed to write {path}"))
}

/// Streams `fastx` into `out`; returns the reads `[kept, flipped, ambiguous, no_signal]`.
///
/// polyA only keeps the read, polyT only flips it, and both or neither keep it unchanged.
fn flip(fastx: &Path, out: &mut impl Write) -> Result<[usize; 4]> {
    let mut counts = [0; 4];

    each(fastx, |rec| {
        let mut seq = rec.seq().into_owned();
        let mut qual = rec.qual().map(<[u8]>::to_vec);

        let class = match tails(&seq) {
            (true, false) => 0,
            (false, true) => 1,
            (true, true) => 2,
            (false, false) => 3,
        };
        counts[class] += 1;

        if class == 1 {
            seq = seq.reverse_complement();
            if let Some(qual) = &mut qual {
                qual.reverse();
            }
        }

        match qual {
            Some(qual) => write_fastq(rec.id(), &seq, Some(&qual), out, LineEnding::Unix)?,
            None => write_fasta(rec.id(), &seq, out, LineEnding::Unix)?,
        }

        Ok(true)
    })?;

    Ok(counts)
}

/// File name of `path`, for the `file` column of the reports.
fn file_name(path: &Path) -> String {
    path.file_name()
        .unwrap_or(path.as_os_str())
        .to_string_lossy()
        .into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// Deterministic ACGT body: no tails, no primers.
    fn body(n: usize, seed: usize) -> String {
        let mut x = seed as u32;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1664525).wrapping_add(1013904223);
                b"ACGT"[(x >> 24) as usize & 3] as char
            })
            .collect()
    }

    fn rc(seq: &str) -> String {
        let comp = |c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        };
        seq.chars().rev().map(comp).collect()
    }

    /// Read generator: a 300 nt body, wrapped in `head` and `tail` for the first `k` reads.
    fn spiked(k: usize, head: &str, tail: &str) -> impl Fn(usize) -> String {
        let (head, tail) = (head.to_string(), tail.to_string());
        move |i| match i < k {
            true => format!("{head}{}{tail}", body(300, i)),
            false => body(300, i),
        }
    }

    /// 10 records named by `name`: FASTQ with every base at `qual` if given, else FASTA.
    fn records(
        qual: Option<char>,
        name: impl Fn(usize) -> String,
        seq: impl Fn(usize) -> String,
    ) -> String {
        let record = |i| match qual {
            Some(q) => format!(
                "@{}\n{}\n+\n{}\n",
                name(i),
                seq(i),
                q.to_string().repeat(seq(i).len())
            ),
            None => format!(">{}\n{}\n", name(i), seq(i)),
        };
        (0..10).map(record).collect()
    }

    #[test]
    fn inspect_reaches_every_state() {
        let dir = tempfile::tempdir().unwrap();
        let dir = dir.path();
        let (neb, x, user) = (
            "GCAATGAAGTCGCAGGGTTGGG",
            "CTACACGACGCTCTTCCGATCT",
            "ACGTTGCATGCCGATTACGGATCC",
        );
        let user_fa = dir.join("user.fa");
        std::fs::write(
            &user_fa,
            format!(">mine my primer\n{}\n", user.to_lowercase()),
        )
        .unwrap();

        let polya_primer = format!("{}{}", "A".repeat(20), rc("AAGCAGTGGTATCAACGCAGAGTAC"));
        let read = |i| format!("read{i}");
        let long = |i| body(9000, i);
        // exact prefix at the head, or reverse-complemented with one mismatch at the tail
        let isoseqx = |i| match i % 2 {
            0 => format!("{x}{}", body(300, i)),
            _ => format!(
                "{}{}",
                body(300, i),
                rc(&format!("{}G{}", &x[..5], &x[6..]))
            ),
        };
        let custom = |i| {
            let (head, tail) = (if i < 8 { neb } else { "" }, if i < 6 { user } else { "" });
            format!("{head}{}{tail}", body(200, i))
        };
        let mixed = |i| match i {
            0..=2 => format!("{}{}", "T".repeat(20), body(300, i)),
            3 | 4 => format!("{}{}", body(300, i), "A".repeat(20)),
            _ => body(300, i),
        };

        #[rustfmt::skip]
        let cases = [
            ("empty.fa", String::new(), "reads=0 state=empty kit=none"),
            ("clustered.fa", records(None, |i| format!("transcript/{i}"), spiked(0, "", "")), "state=clustered"),
            ("subreads.fa", records(None, |i| format!("m1/{i}/0_300"), spiked(0, "", "")), "state=subreads"),
            ("lowqv.fq", records(Some('+'), read, spiked(0, "", "")), "state=subreads mean_qv=10.00"),
            // a constant Q40 is no subreads signal
            ("ccs.fq", records(Some('I'), |i| format!("m1/{i}/ccs"), spiked(4, neb, "")), "state=ccs kit=express primer_rate=0.4000 mean_qv=40.00 qv_constant=true name_style=pacbio"),
            ("isoseqx.fa", records(None, read, isoseqx), "state=ccs kit=isoseqx primer_rate=1.0000 mean_qv=NA qv_constant=NA"),
            // the user primer outranks the NEB one within a read
            ("custom.fa", records(None, read, custom), "state=ccs kit=custom primer_rate=0.8000"),
            // the only primer sits behind a polyA run, which must not shadow it
            ("partial.fa", records(None, read, spiked(1, "", &polya_primer)), "state=ambiguous kit=express note=primer_rate>=0.02_partial_primer_signal"),
            ("long.fa", records(None, read, long), "state=ambiguous kit=none note=mean_len>8000_possible_unsegmented_arrays_need_skera"),
            // polyT outranks polyA
            ("mixed.fa", records(None, read, mixed), "state=mixed polyt5=0.3000 polya3=0.2000 note=balanced_tails_unknown_primers?"),
            ("fl.fa", records(None, |i| format!("SRR1.{i}"), spiked(2, "", &"A".repeat(20))), "state=fl name_style=sra note=-"),
            ("flnc.fa", records(None, read, spiked(0, "", "")), "state=flnc name_style=other"),
        ];

        for (file, text, want) in cases {
            let fastx = dir.join(file);
            std::fs::write(&fastx, text).unwrap();
            let prefix = fastx.display().to_string();
            let primers = (file == "custom.fa").then_some(user_fa.as_path());
            inspect(&fastx, &prefix, 5000, primers).unwrap();

            let tsv = std::fs::read_to_string(format!("{prefix}.inspect.tsv")).unwrap();
            let mut lines = tsv.lines().map(|line| line.split('\t'));
            let row: Vec<_> = lines.next().unwrap().zip(lines.next().unwrap()).collect();
            assert_eq!(row[0], ("file", file));
            for pair in want.split(' ') {
                let pair_kv = pair.split_once('=').unwrap();
                assert!(row.contains(&pair_kv), "{file}: {pair} not in {row:?}");
            }
        }

        for (header, style) in [
            ("transcript/5", "clustered"),
            ("PB.3.2 x", "clustered"),
            ("a full_length_coverage=2", "clustered"),
            ("m1/7/0_900", "pacbio_subread"),
            ("m1/7/0_900 RQ=0.9", "pacbio_subread"),
            ("m1/7/ccs", "pacbio"),
            ("SRR1.2 x", "sra"),
            ("DRR12.3", "sra"),
            ("SRR1", "other"),
            ("PB.1", "other"),
            ("transcript/x", "other"),
            ("m1/a/0_900", "other"),
        ] {
            assert_eq!(name_style(header.as_bytes()), style, "{header}");
        }
    }

    #[test]
    fn orient_flips_only_polyt_reads() {
        let dir = tempfile::tempdir().unwrap();
        let dir = dir.path();
        let (a, t) = ("A".repeat(20), "T".repeat(20));
        let reads = [
            ("fwd", body(80, 1) + &a),
            ("rev", t.clone() + &body(80, 2)),
            ("both", t + &body(80, 3) + &a),
            ("none", body(100, 4)),
        ];
        // distinct qualities along the read make their reversal visible
        let qual = |s: &str| -> String {
            (0..s.len())
                .map(|i| (b'!' + i as u8 % 40) as char)
                .collect()
        };
        // the file as written (`flipped` false) or as orient must write it
        let text = |fastq: bool, flipped: bool| -> String {
            let record = |(name, s): &(&str, String)| {
                let (s, q) = match name == &"rev" && flipped {
                    true => (rc(s), qual(s).chars().rev().collect()),
                    false => (s.clone(), qual(s)),
                };
                match fastq {
                    true => format!("@{name}\n{s}\n+\n{q}\n"),
                    false => format!(">{name}\n{s}\n"),
                }
            };
            reads.iter().map(record).collect()
        };

        for (input, output, fastq) in [
            ("in.fastq", "out.fastq.gz", true),
            ("in.fa", "out.fa", false),
        ] {
            let (input, output) = (dir.join(input), dir.join(output));
            std::fs::write(&input, text(fastq, false)).unwrap();
            let prefix = input.display().to_string();
            orient(&input, &output, &prefix).unwrap();

            let bytes = std::fs::read(&output).unwrap();
            let mut written = String::new();
            match output.extension().unwrap() == "gz" {
                true => flate2::read::GzDecoder::new(&bytes[..]).read_to_string(&mut written),
                false => bytes.as_slice().read_to_string(&mut written),
            }
            .unwrap();
            assert_eq!(written, text(fastq, true));

            let tsv = std::fs::read_to_string(format!("{prefix}.orient.tsv")).unwrap();
            let name = input.file_name().unwrap().to_string_lossy();
            let header = "file\treads\tkept\tflipped\tambiguous\tno_signal";
            assert_eq!(tsv, format!("{header}\n{name}\t4\t1\t1\t1\t1\n"));
        }
    }
}
