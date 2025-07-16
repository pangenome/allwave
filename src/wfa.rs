use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, HeuristicStrategy, MemoryMode,
};
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
        // IMPORTANT: WFA2 uses opposite convention for I/D operations!
        // WFA2: I = consume reference, D = consume query
        // Standard: I = consume query, D = consume reference
        // So we swap them here
        let op_char = match op {
            b'M' => '=', // Match (WFA2 uses M for exact match)
            b'X' => 'X', // Mismatch
            b'I' => 'D', // WFA2 'I' means standard 'D'
            b'D' => 'I', // WFA2 'D' means standard 'I'
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
            // WFA2 uses opposite convention - swap I and D
            b'I' => deletions += 1,  // WFA2 'I' means standard 'D'
            b'D' => insertions += 1, // WFA2 'D' means standard 'I'
            _ => {}
        }
    }

    let alignment_length = matches + mismatches;
    (matches, mismatches, insertions, deletions, alignment_length)
}

pub fn validate_cigar_alignment(
    cigar: &[u8],
    query: &[u8],
    reference: &[u8],
) -> Result<(), String> {
    let mut q_pos = 0;
    let mut r_pos = 0;

    // IMPORTANT: WFA2 CIGAR describes transforming reference into query
    // This is opposite of standard CIGAR convention!
    // In WFA2: I = insert from query, D = delete from reference
    for &op in cigar {
        match op {
            b'M' | b'=' | b'X' => {
                if q_pos >= query.len() || r_pos >= reference.len() {
                    return Err(format!(
                        "CIGAR extends beyond sequences at M/=/X op: q_pos={}, r_pos={}, query_len={}, ref_len={}",
                        q_pos, r_pos, query.len(), reference.len()
                    ));
                }
                q_pos += 1;
                r_pos += 1;
            }
            b'I' => {
                // In WFA2: I means consume reference (opposite of standard)
                if r_pos >= reference.len() {
                    return Err(format!(
                        "CIGAR extends beyond reference at I op: r_pos={}, ref_len={}",
                        r_pos,
                        reference.len()
                    ));
                }
                r_pos += 1;
            }
            b'D' => {
                // In WFA2: D means consume query (opposite of standard)
                if q_pos >= query.len() {
                    return Err(format!(
                        "CIGAR extends beyond query at D op: q_pos={}, query_len={}",
                        q_pos,
                        query.len()
                    ));
                }
                q_pos += 1;
            }
            _ => {
                return Err(format!(
                    "Invalid CIGAR operation: {} (0x{:02x})",
                    op as char, op
                ));
            }
        }
    }

    if q_pos != query.len() {
        return Err(format!(
            "CIGAR doesn't cover full query: {} vs {}",
            q_pos,
            query.len()
        ));
    }

    if r_pos != reference.len() {
        return Err(format!(
            "CIGAR doesn't cover full reference: {} vs {}",
            r_pos,
            reference.len()
        ));
    }

    Ok(())
}

pub fn align_sequences(
    pattern: &[u8],
    text: &[u8],
    penalties: &Penalties,
    mode: AlignmentMode,
) -> Result<AlignmentResult, AlignmentError> {
    // Create wavefront aligner based on mode
    let mut wf = match mode {
        AlignmentMode::EditDistance => {
            // For edit distance, gap_open = gap_ext = mismatch
            AffineWavefronts::with_penalties_and_memory_mode(
                0, // match
                penalties.mismatch,
                penalties.mismatch, // gap_opening = mismatch for edit distance
                penalties.mismatch, // gap_extension = mismatch for edit distance
                MemoryMode::Ultralow,
            )
        }
        AlignmentMode::SinglePieceAffine => {
            // Single-piece affine gap penalties
            AffineWavefronts::with_penalties_and_memory_mode(
                0, // match
                penalties.mismatch,
                penalties.gap_opening1,
                penalties.gap_extension1,
                MemoryMode::Ultralow,
            )
        }
        AlignmentMode::TwoPieceAffine => {
            // Two-piece affine gap penalties with ultralow memory mode
            AffineWavefronts::with_penalties_affine2p_and_memory_mode(
                0, // match
                penalties.mismatch,
                penalties.gap_opening1,
                penalties.gap_extension1,
                penalties.gap_opening2,
                penalties.gap_extension2,
                MemoryMode::Ultralow,
            )
        }
    };

    // Critical configuration settings for reliable global alignment
    wf.set_alignment_scope(AlignmentScope::Alignment);
    wf.set_alignment_span(AlignmentSpan::End2End);
    wf.set_heuristic(&HeuristicStrategy::None);

    // Perform alignment
    let status = wf.align(pattern, text);

    match status {
        AlignmentStatus::Completed => {
            let score = wf.score();
            let cigar_bytes = wf.cigar();

            // Validate CIGAR before using it
            if let Err(e) = validate_cigar_alignment(cigar_bytes, pattern, text) {
                return Err(AlignmentError {
                    message: format!("CIGAR validation failed: {e}"),
                });
            }

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
