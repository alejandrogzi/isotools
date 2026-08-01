// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for detecting orphans in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the core functions for detecting orphans in a query set of reads
//! and processing the components in a parallel fashion.
//!
//! In short, each component of reads is subjected to any of the two modes: guided or
//! self-guided. Guided mode is the default mode and is used when the user provides a
//! reference file as input. Self-guided mode is used when the user does not.
//! Both, guided and self-guided, cover an extensive amount of curated oprhan cases under
//! the assumption that they do not represent a valid source of evidence for transcription.
//! The process is heavily parallellized to offer fast performance on large datasets.
//!
//! # Evidence
//!
//! Every query read is judged against two independent sources of evidence, applied in
//! sequence: agreement with a reference transcript, and support from the other query
//! reads of its component. Reference overlap is evidence, not a routing gate, so a read
//! that fails the reference comparison still gets to be judged on its own merits.
//!
//! Two rules make each source trustworthy. Reference evidence must be a *complete*
//! structural feature shared with a *single* reference transcript, never a lone
//! coordinate pooled across incompatible transcripts. Query evidence is counted
//! leave-one-out, so a read can never certify its own junctions. Component size, being a
//! property of how reads were grouped rather than of the transcripts themselves, is
//! reported as context but never decides an outcome for a multi-exon read.
//!
//! # Outputs
//!
//! Every run writes three files to `<outdir>/orphans/`:
//!
//! * `<prefix>.hq.bed` — records classified as `PASS`
//! * `<prefix>.scraps.bed` — records classified as `SCRAP`
//! * `<prefix>.report.tsv` — one row per processed query record, explaining its
//!   classification
//!
//! The report is always generated and needs no flag. BED selection and report rows come
//! from the same [`report::ClassifiedRead`] values, so they cannot disagree. See
//! [`report`] and the crate README for the column definitions and reason vocabulary.
//!
//! # Usage
//!
//! The program can be run using the following command:
//!
//! ```bash
//! iso-orphan <ARGUMENTS>
//! ```
//!
//! ## Examples
//!
//! ```bash
//! iso-orphan --bed tests/data/mm39.bed --toga tests/data/mm39.toga.bed --all --outdir tests/data --name test --threads 1
//! ```
//!
//! ```bash
//! iso-orphan --bed tests/data/mm39.bed --all --outdir tests/data --name test --threads 1
//! ```
//!
//! ```bash
//!  iso-orphan -b tests/data/mm39.bed -t tests/data/mm39.toga.bed -O test -T 1
//! ```

pub mod cli;
pub mod core;
pub mod denovo;
pub mod report;
pub mod scoring;
pub mod splice;
pub mod utils;
