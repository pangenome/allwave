//! Core alignment functionality

use crate::types::{AlignmentError, AlignmentMode, AlignmentParams, AlignmentResult, Sequence};
use edlib_rs::edlibrs::*;
use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, HeuristicStrategy, MemoryMode,
};
use std::cell::RefCell;

thread_local! {
    static WFA_ALIGNER: RefCell<Option<AffineWavefronts>> = const { RefCell::new(None) };
    static WFA_ORIENTATION: RefCell<Option<AffineWavefronts>> = const { RefCell::new(None) };
}

/// Perform alignment between two sequences
pub fn align_pair(
    query: &Sequence,
    target: &Sequence,
    query_idx: usize,
    target_idx: usize,
    params: &AlignmentParams,
    _orientation_params: &AlignmentParams,
    use_mash_orientation: bool,
) -> AlignmentResult {
    // Determine best orientation using selected method
    let (query_seq, is_reverse) = if use_mash_orientation {
        determine_orientation_mash(&query.seq, &target.seq)
    } else {
        determine_orientation_edlib(&query.seq, &target.seq)
    };

    // Perform the actual alignment with WFA2
    match perform_wfa_alignment(&query_seq, &target.seq, params) {
        Ok(mut result) => {
            result.query_idx = query_idx;
            result.target_idx = target_idx;
            result.is_reverse = is_reverse;
            result
        }
        Err(_) => {
            // Return empty alignment on failure
            AlignmentResult {
                query_idx,
                target_idx,
                query_start: 0,
                query_end: 0,
                target_start: 0,
                target_end: 0,
                is_reverse,
                cigar_bytes: vec![],
                score: i32::MAX,
                num_matches: 0,
                alignment_length: 0,
            }
        }
    }
}

/// Determine best orientation for alignment using strand-specific mash distance
fn determine_orientation_mash(query: &[u8], target: &[u8]) -> (Vec<u8>, bool) {
    const ORIENTATION_KMER_SIZE: usize = 15;

    // Use fixed sketch size for consistent performance
    // 1000 is a good balance between accuracy and speed
    // For very low ANI (< 70%), consider increasing to 2000
    let sketch_size = 1000;

    // Create strand-specific sketches (without canonicalization)
    let target_sketch = sketch_sequence_stranded(target, ORIENTATION_KMER_SIZE, sketch_size);
    let fwd_sketch = sketch_sequence_stranded(query, ORIENTATION_KMER_SIZE, sketch_size);

    let rev_seq = reverse_complement(query);
    let rev_sketch = sketch_sequence_stranded(&rev_seq, ORIENTATION_KMER_SIZE, sketch_size);

    // Compute Jaccard similarities
    let fwd_jaccard = jaccard_similarity(&fwd_sketch, &target_sketch);
    let rev_jaccard = jaccard_similarity(&rev_sketch, &target_sketch);

    // Choose the better orientation (higher Jaccard = more similar)
    if fwd_jaccard >= rev_jaccard {
        (query.to_vec(), false)
    } else {
        (rev_seq, true)
    }
}

/// Create strand-specific MinHash sketch (no canonicalization)
fn sketch_sequence_stranded(sequence: &[u8], k: usize, sketch_size: usize) -> Vec<u64> {
    if sequence.len() < k {
        return Vec::new();
    }

    let mut hashes = Vec::new();

    // Extract all k-mers and hash them (no canonicalization)
    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];

        // Skip k-mers containing non-ACGT characters
        if kmer.iter().any(|&b| !is_dna_base(b)) {
            continue;
        }

        // Hash the k-mer directly (no reverse complement comparison)
        let hash = hash_kmer(kmer);
        hashes.push(hash);
    }

    // Sort and take the smallest hashes (MinHash)
    hashes.sort_unstable();
    hashes.truncate(sketch_size);
    hashes
}

/// Compute Jaccard similarity between two sketches
fn jaccard_similarity(sketch1: &[u64], sketch2: &[u64]) -> f64 {
    use std::collections::HashSet;

    let set1: HashSet<_> = sketch1.iter().collect();
    let set2: HashSet<_> = sketch2.iter().collect();

    let intersection_size = set1.intersection(&set2).count();
    let union_size = set1.union(&set2).count();

    if union_size == 0 {
        0.0
    } else {
        intersection_size as f64 / union_size as f64
    }
}

/// Hash a k-mer using a simple hash function
fn hash_kmer(kmer: &[u8]) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut hasher = DefaultHasher::new();
    kmer.hash(&mut hasher);
    hasher.finish()
}

/// Check if byte represents a DNA base
fn is_dna_base(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Determine best orientation for alignment using fast edit distance (edlib)
fn determine_orientation_edlib(query: &[u8], target: &[u8]) -> (Vec<u8>, bool) {
    let config = EdlibAlignConfigRs {
        k: -1,                                       // No limit on edit distance
        mode: EdlibAlignModeRs::EDLIB_MODE_NW,       // Global alignment
        task: EdlibAlignTaskRs::EDLIB_TASK_DISTANCE, // Just get distance
        additionalequalities: &[],
    };

    // Get edit distance for forward orientation
    let fwd_result = edlibAlignRs(query, target, &config);
    let fwd_distance = fwd_result.editDistance;

    // Get edit distance for reverse complement orientation
    let rev_seq = reverse_complement(query);
    let rev_result = edlibAlignRs(&rev_seq, target, &config);
    let rev_distance = rev_result.editDistance;

    // Choose the better orientation
    if fwd_distance <= rev_distance {
        (query.to_vec(), false)
    } else {
        (rev_seq, true)
    }
}

/// Reverse complement a DNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => b'N',
        })
        .collect()
}

/// Perform WFA2 alignment
fn perform_wfa_alignment(
    query: &[u8],
    target: &[u8],
    params: &AlignmentParams,
) -> Result<AlignmentResult, AlignmentError> {
    let mode = AlignmentMode::from_params(params);

    // Create or reuse WFA2 aligner
    WFA_ALIGNER.with(|aligner_cell| {
        let mut aligner_opt = aligner_cell.borrow_mut();

        let wf = match &mut *aligner_opt {
            Some(wf) => wf,
            None => {
                // Create new aligner based on mode
                let new_wf = match mode {
                    AlignmentMode::EditDistance => {
                        AffineWavefronts::with_penalties_and_memory_mode(
                            params.match_score,
                            params.mismatch_penalty,
                            params.mismatch_penalty,
                            params.mismatch_penalty,
                            MemoryMode::Ultralow,
                        )
                    }
                    AlignmentMode::SinglePieceAffine => {
                        AffineWavefronts::with_penalties_and_memory_mode(
                            params.match_score,
                            params.mismatch_penalty,
                            params.gap_open,
                            params.gap_extend,
                            MemoryMode::Ultralow,
                        )
                    }
                    AlignmentMode::TwoPieceAffine => {
                        AffineWavefronts::with_penalties_affine2p_and_memory_mode(
                            params.match_score,
                            params.mismatch_penalty,
                            params.gap_open,
                            params.gap_extend,
                            params.gap2_open.unwrap_or(params.gap_open),
                            params.gap2_extend.unwrap_or(params.gap_extend),
                            MemoryMode::Ultralow,
                        )
                    }
                };
                *aligner_opt = Some(new_wf);
                aligner_opt.as_mut().unwrap()
            }
        };

        // Configure aligner
        wf.set_alignment_scope(AlignmentScope::Alignment);
        wf.set_alignment_span(AlignmentSpan::End2End);
        wf.set_heuristic(&HeuristicStrategy::None);

        // Perform alignment
        let status = wf.align(query, target);

        match status {
            AlignmentStatus::Completed => {
                let score = wf.score();
                let cigar_bytes = wf.cigar().to_vec();

                // Validate and process CIGAR
                let (num_matches, alignment_length) = count_cigar_operations(&cigar_bytes);
                let (query_end, target_end) = parse_cigar_lengths(&cigar_bytes);

                Ok(AlignmentResult {
                    query_idx: 0,  // Will be set by caller
                    target_idx: 0, // Will be set by caller
                    query_start: 0,
                    query_end,
                    target_start: 0,
                    target_end,
                    is_reverse: false, // Will be set by caller
                    cigar_bytes,
                    score,
                    num_matches,
                    alignment_length,
                })
            }
            _ => Err(AlignmentError {
                message: format!("Alignment failed with status: {status:?}"),
            }),
        }
    })
}

/// Count operations in CIGAR bytes
fn count_cigar_operations(cigar_bytes: &[u8]) -> (usize, usize) {
    let mut matches = 0;
    let mut alignment_length = 0;

    for &op in cigar_bytes {
        match op {
            b'M' => {
                matches += 1;
                alignment_length += 1;
            }
            b'X' => {
                alignment_length += 1;
            }
            _ => {}
        }
    }

    (matches, alignment_length)
}

/// Parse CIGAR to get alignment lengths
fn parse_cigar_lengths(cigar_bytes: &[u8]) -> (usize, usize) {
    let mut query_len = 0;
    let mut target_len = 0;

    // WFA2 uses opposite convention for I/D operations
    for &op in cigar_bytes {
        match op {
            b'M' | b'X' => {
                query_len += 1;
                target_len += 1;
            }
            b'I' => {
                // WFA2 'I' means consume reference (standard 'D')
                target_len += 1;
            }
            b'D' => {
                // WFA2 'D' means consume query (standard 'I')
                query_len += 1;
            }
            _ => {}
        }
    }

    (query_len, target_len)
}

/// Convert CIGAR bytes to standard string representation
pub fn cigar_bytes_to_string(cigar_bytes: &[u8]) -> String {
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

        // Convert to CIGAR format with I/D swap for WFA2
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
