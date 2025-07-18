//! Allwave - High-performance pairwise sequence aligner
//!
//! This library provides efficient all-vs-all sequence alignment using
//! bidirectional wavefront alignment (biWFA).

pub mod alignment;
pub mod iterator;
pub mod knn_graph;
pub mod mash;
pub mod neighbor_joining;
pub mod types;

// Keep existing modules for compatibility
pub mod test_framework;
pub mod validation;
pub mod validation_correct;
pub mod validation_simple;
pub mod wfa;

// Re-export main types and functions
pub use alignment::{cigar_bytes_to_string, reverse_complement};
pub use iterator::{AllPairIterator, AllPairParallelIterator};
pub use types::{
    AlignmentError, AlignmentMode, AlignmentParams, AlignmentResult, Sequence,
    SparsificationStrategy,
};

/// Process all-vs-all alignments with a callback for streaming results
///
/// # Arguments
/// * `sequences` - Vector of sequences to align
/// * `params` - Alignment parameters
/// * `sparsification` - Sparsification strategy
/// * `callback` - Function called for each alignment result
///
/// # Example
/// ```
/// use allwave::{process_alignments_with_callback, Sequence, AlignmentParams, SparsificationStrategy};
///
/// let sequences = vec![
///     Sequence { id: "seq1".to_string(), seq: b"ATCG".to_vec() },
///     Sequence { id: "seq2".to_string(), seq: b"CGTA".to_vec() },
/// ];
///
/// let params = AlignmentParams::edit_distance();
///
/// process_alignments_with_callback(
///     &sequences,
///     params,
///     SparsificationStrategy::None,
///     |alignment| {
///         println!("Alignment: {} vs {}", alignment.query_idx, alignment.target_idx);
///         Ok(())
///     }
/// ).unwrap();
/// ```
pub fn process_alignments_with_callback<F>(
    sequences: &[Sequence],
    params: AlignmentParams,
    sparsification: SparsificationStrategy,
    callback: F,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
where
    F: Fn(AlignmentResult) -> Result<(), Box<dyn std::error::Error + Send + Sync>> + Send + Sync,
{
    let aligner = AllPairIterator::with_options(sequences, params, true, true, sparsification);
    aligner.for_each_with_callback(callback)
}

/// Format alignment result as PAF record
pub fn alignment_to_paf(result: &AlignmentResult, sequences: &[Sequence]) -> String {
    let query = &sequences[result.query_idx];
    let target = &sequences[result.target_idx];

    let query_len = query.seq.len();
    let target_len = target.seq.len();

    let query_aligned_len = result.query_end - result.query_start;
    let target_aligned_len = result.target_end - result.target_start;
    let block_len = target_aligned_len.max(query_aligned_len);

    // Calculate identity
    let identity = if result.alignment_length > 0 {
        result.num_matches as f64 / result.alignment_length as f64
    } else {
        0.0
    };

    // Convert CIGAR bytes to string
    let cigar = cigar_bytes_to_string(&result.cigar_bytes);

    // Determine strand
    let strand = if result.is_reverse { '-' } else { '+' };

    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgi:f:{:.6}\tcg:Z:{}",
        query.id,
        query_len,
        result.query_start,
        result.query_end,
        strand,
        target.id,
        target_len,
        result.target_start,
        result.target_end,
        result.num_matches,
        block_len,
        60, // mapping quality placeholder
        identity,
        cigar
    )
}

/// Parse alignment scores from string format
/// Format: "match,mismatch,gap_open,gap_ext[,gap_open2,gap_ext2]"
pub fn parse_scores(scores_str: &str) -> Result<AlignmentParams, String> {
    let scores: Vec<i32> = scores_str
        .split(',')
        .map(|s| s.trim().parse::<i32>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("Failed to parse scores: {e}"))?;

    match scores.len() {
        4 => {
            // Edit distance or single-piece affine
            Ok(AlignmentParams {
                match_score: scores[0],
                mismatch_penalty: scores[1],
                gap_open: scores[2],
                gap_extend: scores[3],
                gap2_open: None,
                gap2_extend: None,
                max_divergence: None,
            })
        }
        6 => {
            // Two-piece affine
            Ok(AlignmentParams {
                match_score: scores[0],
                mismatch_penalty: scores[1],
                gap_open: scores[2],
                gap_extend: scores[3],
                gap2_open: Some(scores[4]),
                gap2_extend: Some(scores[5]),
                max_divergence: None,
            })
        }
        _ => Err(format!(
            "Invalid number of scores: {}. Expected 4 or 6 values.",
            scores.len()
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_scores() {
        // Test edit distance
        let params = parse_scores("0,1,1,1").unwrap();
        assert_eq!(params.match_score, 0);
        assert_eq!(params.mismatch_penalty, 1);
        assert_eq!(params.gap_open, 1);
        assert_eq!(params.gap_extend, 1);
        assert!(params.gap2_open.is_none());

        // Test two-piece affine
        let params = parse_scores("0,5,8,2,24,1").unwrap();
        assert_eq!(params.match_score, 0);
        assert_eq!(params.mismatch_penalty, 5);
        assert_eq!(params.gap_open, 8);
        assert_eq!(params.gap_extend, 2);
        assert_eq!(params.gap2_open, Some(24));
        assert_eq!(params.gap2_extend, Some(1));
    }

    #[test]
    fn test_alignment_mode_detection() {
        let edit_params = AlignmentParams::edit_distance();
        assert!(matches!(
            AlignmentMode::from_params(&edit_params),
            AlignmentMode::EditDistance
        ));

        let two_piece = parse_scores("0,5,8,2,24,1").unwrap();
        assert!(matches!(
            AlignmentMode::from_params(&two_piece),
            AlignmentMode::TwoPieceAffine
        ));
    }
}
