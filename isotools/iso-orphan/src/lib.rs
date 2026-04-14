// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for detecting orphans in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the core functions for detecting orphans in a query set of reads
//! and processing the components in a parallel fashion.
//!
//! In short, each component of reads is subjected to any of the two modes: guided or
//! self-guided. Guided mode is the default mode and is used when the user provides a TOGA file
//! as input. Self-guided mode is used when the user does not provide a TOGA file as input.
//! Both, guided and self-guided, cover an extensive amount of curated oprhan cases under
//! the assumption that they do not represent a valid source of evidence for transcription.
//! The process is heavily parallellized to offer fast performance on large datasets.
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
pub mod utils;
