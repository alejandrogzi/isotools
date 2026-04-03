//! Core module for detecting intron retentions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2026
//!
//! This module contains the main algorithm for detecting truncations
//! in a query set of reads.
//!
//! In short it takes a set of query reads and a set of reference reads and
//! detects truncations in the query reads. It does this by checking if the
//! query reads overlap with any middle exon from the reference set of reads.
//! A recovery step can be performed by evaluating the support of the middle
//! exon in the reference set of reads.

pub mod cli;
pub mod core;
pub mod utils;
