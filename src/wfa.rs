use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus};
use std::error::Error;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub enum AlignmentMode {
    EditDistance,
    SinglePieceAffine,
    TwoPieceAffine,
}

#[derive(Debug)]
pub struct AlignmentError {
    message: String,
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl Error for AlignmentError {}

pub struct Penalties {
    pub mismatch: i32,
    pub gap_opening1: i32,
    pub gap_extension1: i32,
    pub gap_opening2: i32,
    pub gap_extension2: i32,
}

pub struct AlignmentResult {
    #[allow(dead_code)]
    pub score: i32,
    pub cigar: String,
    pub matches: usize,
    #[allow(dead_code)]
    pub mismatches: usize,
    #[allow(dead_code)]
    pub insertions: usize,
    #[allow(dead_code)]
    pub deletions: usize,
    pub alignment_length: usize,
}

fn cigar_bytes_to_string(cigar_bytes: &[u8]) -> String {
    let mut cigar_str = String::new();
    let mut i = 0;
    
    while i < cigar_bytes.len() {
        let op = cigar_bytes[i];
        let mut count = 1;
        let mut j = i + 1;
        
        // Count consecutive same operations
        while j < cigar_bytes.len() && cigar_bytes[j] == op {
            count += 1;
            j += 1;
        }
        
        // Convert to CIGAR format
        let op_char = match op {
            b'M' => '=',    // Match
            b'X' => 'X',    // Mismatch
            b'I' => 'I',    // Insertion
            b'D' => 'D',    // Deletion
            _ => '?',
        };
        
        cigar_str.push_str(&format!("{count}{op_char}"));
        i = j;
    }
    
    cigar_str
}

fn count_cigar_operations(cigar_bytes: &[u8]) -> (usize, usize, usize, usize, usize) {
    let mut matches = 0;
    let mut mismatches = 0;
    let mut insertions = 0;
    let mut deletions = 0;
    
    for &op in cigar_bytes {
        match op {
            b'M' => matches += 1,
            b'X' => mismatches += 1,
            b'I' => insertions += 1,
            b'D' => deletions += 1,
            _ => {}
        }
    }
    
    let alignment_length = matches + mismatches;
    (matches, mismatches, insertions, deletions, alignment_length)
}

pub fn align_sequences(
    pattern: &[u8],
    text: &[u8],
    penalties: &Penalties,
    mode: AlignmentMode,
) -> Result<AlignmentResult, AlignmentError> {
    // Create wavefront aligner based on mode
    let wf = match mode {
        AlignmentMode::EditDistance => {
            // For edit distance, gap_open = gap_ext = mismatch
            AffineWavefronts::with_penalties(
                0,  // match
                penalties.mismatch,
                penalties.mismatch, // gap_opening = mismatch for edit distance
                penalties.mismatch, // gap_extension = mismatch for edit distance
            )
        }
        AlignmentMode::SinglePieceAffine => {
            // Single-piece affine gap penalties
            AffineWavefronts::with_penalties(
                0,  // match
                penalties.mismatch,
                penalties.gap_opening1,
                penalties.gap_extension1,
            )
        }
        AlignmentMode::TwoPieceAffine => {
            // Two-piece affine gap penalties
            AffineWavefronts::with_penalties_affine2p(
                0,  // match
                penalties.mismatch,
                penalties.gap_opening1,
                penalties.gap_extension1,
                penalties.gap_opening2,
                penalties.gap_extension2,
            )
        }
    };
    
    // Perform alignment
    let status = wf.align(pattern, text);
    
    match status {
        AlignmentStatus::Completed => {
            let score = wf.score();
            let cigar_bytes = wf.cigar();
            let cigar = cigar_bytes_to_string(cigar_bytes);
            let (matches, mismatches, insertions, deletions, alignment_length) = 
                count_cigar_operations(cigar_bytes);
            
            Ok(AlignmentResult {
                score,
                cigar,
                matches,
                mismatches,
                insertions,
                deletions,
                alignment_length,
            })
        }
        _ => Err(AlignmentError {
            message: format!("Alignment failed with status: {status:?}"),
        }),
    }
}