//! Core types for the allwave library

use std::fmt;

/// A sequence with an identifier and data
#[derive(Debug, Clone)]
pub struct Sequence {
    pub id: String,
    pub seq: Vec<u8>,
}

/// Result of a pairwise alignment
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Indices into the original sequence array
    pub query_idx: usize,
    pub target_idx: usize,
    
    /// Alignment coordinates (0-based)
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    
    /// Orientation - true if query was reverse-complemented
    pub is_reverse: bool,
    
    /// Alignment details
    pub cigar_bytes: Vec<u8>,       // Raw CIGAR from WFA2
    pub score: i32,                  // WFA2 score (lower is better)
    pub num_matches: usize,          // Number of matching bases
    pub alignment_length: usize,     // Total alignment length
}

/// Parameters for sequence alignment
#[derive(Debug, Clone)]
pub struct AlignmentParams {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub gap2_open: Option<i32>,
    pub gap2_extend: Option<i32>,
    pub max_divergence: Option<f64>,
}

impl Default for AlignmentParams {
    fn default() -> Self {
        Self {
            match_score: 0,
            mismatch_penalty: 5,
            gap_open: 8,
            gap_extend: 2,
            gap2_open: Some(24),
            gap2_extend: Some(1),
            max_divergence: None,
        }
    }
}

impl AlignmentParams {
    /// Create parameters for edit distance alignment
    pub fn edit_distance() -> Self {
        Self {
            match_score: 0,
            mismatch_penalty: 1,
            gap_open: 1,
            gap_extend: 1,
            gap2_open: None,
            gap2_extend: None,
            max_divergence: None,
        }
    }
}

/// Sparsification strategy for all-vs-all alignments
#[derive(Debug, Clone)]
pub enum SparsificationStrategy {
    /// Align all pairs
    None,
    /// Random sampling with given probability
    Random(f64),
    /// Automatic sparsification based on sequence count
    Auto,
}

/// Alignment mode (edit distance, single-piece affine, two-piece affine)
#[derive(Debug, Clone, Copy)]
pub enum AlignmentMode {
    EditDistance,
    SinglePieceAffine,
    TwoPieceAffine,
}

impl AlignmentMode {
    /// Determine alignment mode from parameters
    pub fn from_params(params: &AlignmentParams) -> Self {
        if params.gap2_open.is_some() && params.gap2_extend.is_some() {
            AlignmentMode::TwoPieceAffine
        } else if params.gap_open == params.gap_extend && params.gap_open == params.mismatch_penalty {
            AlignmentMode::EditDistance
        } else {
            AlignmentMode::SinglePieceAffine
        }
    }
}

/// Error type for alignment operations
#[derive(Debug)]
pub struct AlignmentError {
    pub message: String,
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for AlignmentError {}