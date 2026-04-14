// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core modul for rescuing missed 3' splice junctions by CIGAR matching
//! Alejandro Gonzales-Irribarren, 2026
//!
//! This module contains the main function for rescuing missed 3' splice
//! junctions by CIGAR matching. It is a wrapper around the `engine` module
//! that handles the command-line interface and logging. In short, it reads
//! a BAM file and a GTF file, and writes a new BAM file with rescued junctions.
//!
//! A minimap2-aligned transcript can end exactly at an exon boundary
//! but carry a small 3' soft-clipped sequence that was not placed anywhere.
//! The hypothesis is: that clipped sequence is the start of the next downstream
//! exon, but the aligner missed the intron.
//!
//! In short, iso-cigar identifies candidate reads: those ending ±wiggle bp from a
//! known internal exon boundary, with soft-clip ≥ cutoff; retrieves the sequence of
//! the immediately downstream exon from the reference genome; checks whether the
//! soft-clipped bases match the beginning of that downstream exon; if yes: rewrites
//! the CIGAR to add an N intron gap and converts the soft-clip into a = match

fn main() {
    if let Err(err) = iso_cigar::run() {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}
