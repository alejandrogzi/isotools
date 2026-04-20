// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Known adapter / primer sequences found at the ends of long-read
//! transcriptomics reads.
//!
//! Each entry is `(sequence_bytes, human_readable_label)`. Labels follow the
//! convention `VENDOR:KIT:ROLE[:ORIENT]` so that downstream tooling can
//! group hits. `sequence_bytes` is always stored in its 5'->3' form; both the
//! sequence itself and its reverse complement are enrolled into the matcher
//! automatically by `detector::AdapterDb`.
//!
//! Source references (where public) are noted inline.

/// Minimum length of a sequence we accept into the database. Anything
/// shorter would produce an unacceptable false-positive rate.
pub const MIN_ADAPTER_LEN: usize = 10;

/// `(sequence, label)` pairs of known long-read adapters / primers.
///
/// All sequences use uppercase IUPAC `ACGT`; any ambiguity base is dropped
/// from the stored form so that the exact-match automaton stays pure-`ACGT`.
/// The detector enrolls each entry **plus its reverse complement** so that
/// both 5' and 3' clip orientations are covered.
pub static ADAPTER_DB: &[(&[u8], &str)] = &[
    //
    // ── PacBio Iso-Seq / SMRTbell ─────────────────────────────────────────
    //

    // SMRTbell hairpin adapter (linear form of the blunt-ended hairpin).
    // Source: PacBio Technical Note "Preparing Iso-Seq Libraries Using
    // SMRTbell Express Template Prep Kit 2.0". The canonical linear
    // representation is the concatenation of the two strands around the
    // loop; reads will only ever see the 5' half because the loop is cut
    // during CCS. We keep the half that actually appears in subreads.
    (
        b"ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT",
        "pacbio:smrtbell:hairpin",
    ),
    (
        b"ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT",
        "pacbio:smrtbell:hairpin_rc",
    ),

    // Iso-Seq Clontech / NEB TSO (template switching oligo) — 5' primer.
    // Source: PacBio "Procedure & Checklist - Iso-Seq Express".
    (
        b"AAGCAGTGGTATCAACGCAGAGTACATGGGG",
        "pacbio:isoseq:primer_5p_tso",
    ),
    // Same TSO without the terminal GGG (as sometimes seen after CCS).
    (
        b"AAGCAGTGGTATCAACGCAGAGTACATGGG",
        "pacbio:isoseq:primer_5p_tso_short",
    ),
    // Iso-Seq 3' cDNA primer (without polyT tail); followed in reads by
    // a run of Ts corresponding to the polyA captured by oligo(dT).
    (
        b"AAGCAGTGGTATCAACGCAGAGTAC",
        "pacbio:isoseq:primer_3p_cdna",
    ),
    // NEBNext Single-Cell / Low-Input cDNA primer.
    // Source: NEB E6421 manual.
    (b"AAGCAGTGGTATCAACGCAGAGT", "pacbio:isoseq:nebnext_primer"),

    // SMARTer-Seq v4 (Clontech/Takara) TSO — commonly used upstream of
    // PacBio Iso-Seq library prep.
    (
        b"AAGCAGTGGTATCAACGCAGAGTACATGGG",
        "clontech:smarter_v4:tso",
    ),
    (
        b"AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTT",
        "clontech:smarter_v4:cds_primer",
    ),

    //
    // ── PolyA / PolyT tails ───────────────────────────────────────────────
    //
    // Long homopolymers captured by oligo(dT) priming. Stored once at a
    // representative length; banded edit-distance matching tolerates
    // length variation on either side.
    (b"AAAAAAAAAAAAAAAA", "homopolymer:polyA"),
    (b"TTTTTTTTTTTTTTTT", "homopolymer:polyT"),

    //
    // ── Oxford Nanopore ───────────────────────────────────────────────────
    //

    // SQK-LSK (ligation sequencing kit) Y-adapter top strand ("LA").
    // Source: ONT community — "Sequencing adapter sequences".
    (
        b"AATGTACTTCGTTCAGTTACGTATTGCT",
        "ont:sqk_lsk:y_adapter_top",
    ),
    // Reverse-complement side that actually appears at the 3' end of reads
    // after ligation (the "bottom" strand of the Y).
    (
        b"GCAATACGTAACTGAACGAAGTACATT",
        "ont:sqk_lsk:y_adapter_bottom",
    ),

    // SQK-PCS / PCB (cDNA-PCR) Strand-Switching Primer (SSP).
    // Source: ONT PCS109 / PCS110 protocol.
    (b"TTTCTGTTGGTGCTGATATTGCTGGG", "ont:sqk_pcs:ssp"),
    // VN Primer (VNP) — 3' oligo(dT) primer used to capture polyA.
    // The 20×T run plus VN degenerate bases is represented by its constant
    // stem only; the polyT is captured by the `homopolymer:polyT` entry.
    (b"ACTTGCCTGTCGCTCTATCTTC", "ont:sqk_pcs:vnp"),

    // cDNA RT Adapter (CRTA) / cDNA adapter used in PCS109+.
    (b"GCAATACGTAACTGAACGAAGT", "ont:sqk_pcs:crta"),

    // RTA (Reverse Transcription Adapter) — SQK-RNA002/SQK-RNA004 direct-RNA
    // RT adapter top and bottom strands. Source: ONT SQK-RNA002 / SQK-RNA004
    // protocols (the RTA stem is shared between the two chemistries).
    (
        b"GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC",
        "ont:sqk_rna:rta_top",
    ),
    (
        b"GAACCTACTACTATGAGAAGAACAAGAAGCC",
        "ont:sqk_rna:rta_bottom",
    ),

    // RMX — SQK-RNA Motor Protein Adapter (partial, conserved stem).
    (b"CCTGTACTTCGTTCAGTTACGTATTGCT", "ont:sqk_rna:rmx"),

    //
    // ── Common ligation / barcoding adapters ──────────────────────────────
    //

    // Native Barcoding Kit (SQK-NBD) conserved flank around the variable
    // barcode. Source: ONT Native Barcoding Kit (SQK-NBD112/114) manual.
    (b"AAGGTTAACACAAAGACACCGACAACTTTCTT", "ont:sqk_nbd:flank_5p"),
    (b"CAGCACCTCAAGTCTTGATTGCC", "ont:sqk_nbd:flank_3p"),

    // ONT rapid adapter (SQK-RAD / SQK-RBK) — tagmentation-based prep.
    // Source: ONT SQK-RAD004 protocol.
    (
        b"GCTTGGGTGTTTAACCTTTTTGTCAGT",
        "ont:sqk_rad:rapid_adapter",
    ),

    // Illumina TruSeq adapter — occasionally seen when short-read libraries
    // are re-sequenced on long-read platforms, and a common contaminant.
    // Source: Illumina adapter sequences document, 2022.
    (
        b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "illumina:truseq:adapter",
    ),
    (
        b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "illumina:truseq:adapter_r2",
    ),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adapter_db_not_empty() {
        assert!(!ADAPTER_DB.is_empty());
    }

    #[test]
    fn adapter_db_sequences_are_acgt_and_long_enough() {
        for (seq, label) in ADAPTER_DB {
            assert!(
                seq.len() >= MIN_ADAPTER_LEN,
                "{label} sequence shorter than MIN_ADAPTER_LEN",
            );
            for &base in *seq {
                assert!(
                    matches!(base, b'A' | b'C' | b'G' | b'T'),
                    "{label} contains non-ACGT base {:?}",
                    base as char,
                );
            }
        }
    }

    #[test]
    fn adapter_db_labels_unique() {
        use std::collections::HashSet;
        let mut seen: HashSet<&str> = HashSet::new();
        for (_, label) in ADAPTER_DB {
            assert!(seen.insert(label), "duplicate adapter label: {label}");
        }
    }
}
