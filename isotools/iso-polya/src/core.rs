//! Core module for scoring, segmenting and clustering reads
//! based on their polyA features
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for grouping reads
//! and processing components based on polyA features in parallel.
//!
//! In short, this modules provides three subtools, namely: aparent,
//! caller and segment. Each one with a specific goal. The first one
//! runs APARENT, a machine-learning model, to score each read's end.
//! The segment module filters reads based on alignment quality and
//! predicts the polyA tail using a two-state HMM model. Finally, the
//! caller module groups all the previous information and tries to
//! determine the intraprimming potential for each read,

pub mod apa;
pub mod pas;
pub mod segment;
