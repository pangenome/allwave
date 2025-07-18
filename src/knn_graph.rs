//! K-nearest neighbor graph construction for intelligent sequence alignment sampling
//!
//! This module implements k-nearest neighbor graphs based on mash distances between sequences.
//! Each sequence is connected to its k most similar sequences, creating a structured graph
//! that maintains good connectivity while being much sparser than all-vs-all comparisons.

use crate::types::Sequence;

/// Extract sequence pairs using k-nearest neighbor graph with optional stranger-joining
pub fn extract_tree_pairs(
    sequences: &[Sequence],
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
    kmer_size: usize,
) -> Vec<(usize, usize)> {
    if sequences.len() < 2 {
        return Vec::new();
    }

    // Compute distance matrix
    let distance_matrix =
        crate::mash::compute_distance_matrix_with_params(sequences, kmer_size, 1000);

    // Build k-nearest neighbor graph
    let mut all_pairs = Vec::new();

    if k_nearest > 0 {
        let nearest_pairs = build_knn_graph(&distance_matrix, k_nearest, false);
        all_pairs.extend(nearest_pairs);
    }

    // Build k-farthest neighbor graph (stranger-joining)
    if k_farthest > 0 {
        let farthest_pairs = build_knn_graph(&distance_matrix, k_farthest, true);
        all_pairs.extend(farthest_pairs);
    }

    // Add random pairs if requested
    if random_fraction > 0.0 {
        let random_pairs = generate_random_pairs(sequences.len(), random_fraction, sequences);
        all_pairs.extend(random_pairs);
    }

    // Remove duplicates
    all_pairs.sort_unstable();
    all_pairs.dedup();

    all_pairs
}

/// Extract sequence pairs using k-nearest neighbor graph (backward compatibility)
pub fn extract_knn_pairs(
    sequences: &[Sequence],
    k_neighbors: usize,
    random_fraction: f64,
    kmer_size: usize,
) -> Vec<(usize, usize)> {
    extract_tree_pairs(sequences, k_neighbors, 0, random_fraction, kmer_size)
}

/// Build k-nearest or k-farthest neighbor graph from distance matrix
fn build_knn_graph(
    distance_matrix: &[Vec<f64>],
    k_neighbors: usize,
    farthest: bool,
) -> Vec<(usize, usize)> {
    let n = distance_matrix.len();
    let mut pairs = Vec::new();

    // For each sequence, find its k nearest or farthest neighbors
    for (i, row) in distance_matrix.iter().enumerate() {
        // Create list of (distance, index) pairs for all other sequences
        let mut neighbors: Vec<(f64, usize)> =
            (0..n).filter(|&j| i != j).map(|j| (row[j], j)).collect();

        // Sort by distance
        if farthest {
            // For farthest neighbors, sort descending (largest distances first)
            neighbors.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        } else {
            // For nearest neighbors, sort ascending (smallest distances first)
            neighbors.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        }

        // Take the k closest/farthest neighbors
        let k_actual = k_neighbors.min(neighbors.len());
        for &(_, j) in neighbors.iter().take(k_actual) {
            pairs.push((i, j));
        }
    }

    pairs
}

/// Generate random pairs with deterministic hashing
fn generate_random_pairs(n: usize, fraction: f64, sequences: &[Sequence]) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();

    for i in 0..n {
        for j in 0..n {
            if i != j && should_include_pair(i, j, fraction, sequences) {
                pairs.push((i, j));
            }
        }
    }

    pairs
}

/// Deterministic pair inclusion based on sequence IDs
fn should_include_pair(i: usize, j: usize, fraction: f64, sequences: &[Sequence]) -> bool {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut hasher = DefaultHasher::new();
    let seq_i = &sequences[i].id;
    let seq_j = &sequences[j].id;
    let combined = format!("{seq_i}:{seq_j}");
    combined.hash(&mut hasher);

    let hash = hasher.finish();
    let normalized = (hash as f64) / (u64::MAX as f64);
    normalized < fraction
}

/// Estimate the number of pairs in a tree-based sampling strategy
pub fn estimate_tree_pair_count(
    n: usize,
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
) -> usize {
    let nearest_pairs = n * k_nearest.min(n.saturating_sub(1));
    let farthest_pairs = n * k_farthest.min(n.saturating_sub(1));
    let total_possible = n * (n - 1);
    let random_pairs = (total_possible as f64 * random_fraction).round() as usize;
    (nearest_pairs + farthest_pairs + random_pairs).min(total_possible)
}

/// Estimate the number of pairs in a k-nearest neighbor graph (backward compatibility)
pub fn estimate_knn_pair_count(n: usize, k_neighbors: usize, random_fraction: f64) -> usize {
    estimate_tree_pair_count(n, k_neighbors, 0, random_fraction)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Sequence;

    #[test]
    fn test_knn_graph_basic() {
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq2".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            }, // identical to seq1
            Sequence {
                id: "seq3".to_string(),
                seq: b"GGGGGGGGGGGGGGGG".to_vec(),
            }, // different
        ];

        let pairs = extract_knn_pairs(&sequences, 1, 0.0, 15);

        // Each sequence should connect to its 1 nearest neighbor
        // seq1 and seq2 are identical, so they should be each other's nearest neighbors
        // seq3 should connect to one of them
        assert!(pairs.len() >= 2); // At least 2 directed edges
        assert!(pairs.len() <= 3); // At most 3 directed edges
    }

    #[test]
    fn test_build_knn_graph() {
        // Simple 3x3 distance matrix
        let distances = vec![
            vec![0.0, 0.1, 0.9], // seq0 is close to seq1, far from seq2
            vec![0.1, 0.0, 0.8], // seq1 is close to seq0, far from seq2
            vec![0.9, 0.8, 0.0], // seq2 is far from both
        ];

        let pairs = build_knn_graph(&distances, 1, false);

        // Each sequence should have exactly 1 outgoing edge to its nearest neighbor
        assert_eq!(pairs.len(), 3);

        // Check that seq0 connects to seq1 (distance 0.1)
        assert!(pairs.contains(&(0, 1)));
        // Check that seq1 connects to seq0 (distance 0.1)
        assert!(pairs.contains(&(1, 0)));
        // seq2 should connect to either seq0 or seq1 (both closer than the other)
        assert!(pairs.contains(&(2, 0)) || pairs.contains(&(2, 1)));
    }

    #[test]
    fn test_knn_with_k_equals_2() {
        let distances = vec![
            vec![0.0, 0.1, 0.5, 0.9],
            vec![0.1, 0.0, 0.6, 0.8],
            vec![0.5, 0.6, 0.0, 0.2],
            vec![0.9, 0.8, 0.2, 0.0],
        ];

        let pairs = build_knn_graph(&distances, 2, false);

        // Each of 4 sequences should have 2 outgoing edges
        assert_eq!(pairs.len(), 8);
    }

    #[test]
    fn test_estimate_pair_count() {
        assert_eq!(estimate_knn_pair_count(4, 1, 0.0), 4); // 4 sequences, 1 neighbor each
        assert_eq!(estimate_knn_pair_count(4, 2, 0.0), 8); // 4 sequences, 2 neighbors each

        // With random sampling
        let count = estimate_knn_pair_count(4, 1, 0.1);
        assert!(count >= 4); // At least the knn pairs
        assert!(count <= 12); // At most all possible pairs (4*3)
    }

    #[test]
    fn test_empty_sequences() {
        let sequences = vec![];
        let pairs = extract_knn_pairs(&sequences, 1, 0.0, 15);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_single_sequence() {
        let sequences = vec![Sequence {
            id: "seq1".to_string(),
            seq: b"ATCG".to_vec(),
        }];
        let pairs = extract_knn_pairs(&sequences, 1, 0.0, 15);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_tree_sampling_with_strangers() {
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq2".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            }, // identical to seq1
            Sequence {
                id: "seq3".to_string(),
                seq: b"GGGGGGGGGGGGGGGG".to_vec(),
            }, // different from seq1/seq2
        ];

        // Test with 1 nearest + 1 farthest
        let pairs = extract_tree_pairs(&sequences, 1, 1, 0.0, 15);

        // Should have 2 pairs per sequence (1 nearest + 1 farthest)
        // But duplicates are removed, so we might have slightly fewer
        // For 3 sequences: each gets 1 nearest + 1 farthest = 6 total, but some overlap
        assert!(pairs.len() >= 4); // At least some pairs
        assert!(pairs.len() <= 6); // At most all combinations
    }

    #[test]
    fn test_build_knn_graph_farthest() {
        // Simple 3x3 distance matrix
        let distances = vec![
            vec![0.0, 0.1, 0.9], // seq0 is close to seq1, far from seq2
            vec![0.1, 0.0, 0.8], // seq1 is close to seq0, far from seq2
            vec![0.9, 0.8, 0.0], // seq2 is far from both
        ];

        let pairs = build_knn_graph(&distances, 1, true); // Get 1 farthest neighbor

        // Each sequence should have exactly 1 outgoing edge to its farthest neighbor
        assert_eq!(pairs.len(), 3);

        // Check that seq0 connects to seq2 (distance 0.9, the farthest)
        assert!(pairs.contains(&(0, 2)));
        // Check that seq1 connects to seq2 (distance 0.8, the farthest)
        assert!(pairs.contains(&(1, 2)));
        // seq2 should connect to either seq0 or seq1 (both are equidistant in terms of being farthest)
        assert!(pairs.contains(&(2, 0)) || pairs.contains(&(2, 1)));
    }

    #[test]
    fn test_estimate_tree_pair_count() {
        // Test with both nearest and farthest neighbors
        assert_eq!(estimate_tree_pair_count(4, 1, 1, 0.0), 8); // 4 sequences, 1 nearest + 1 farthest each
        assert_eq!(estimate_tree_pair_count(4, 2, 1, 0.0), 12); // 4 sequences, 2 nearest + 1 farthest each

        // With random sampling
        let count = estimate_tree_pair_count(4, 1, 1, 0.1);
        assert!(count >= 8); // At least the tree pairs
        assert!(count <= 12); // At most all possible pairs (4*3)
    }

    #[test]
    fn test_only_nearest_neighbors() {
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq2".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            },
        ];

        // Only nearest neighbors, no strangers
        let pairs = extract_tree_pairs(&sequences, 1, 0, 0.0, 15);
        assert_eq!(pairs.len(), 2); // Each sequence connects to the other as nearest neighbor
    }

    #[test]
    fn test_only_strangers() {
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq2".to_string(),
                seq: b"GGGGGGGGGGGGGGGG".to_vec(),
            },
        ];

        // Only strangers, no nearest neighbors
        let pairs = extract_tree_pairs(&sequences, 0, 1, 0.0, 15);
        assert_eq!(pairs.len(), 2); // Each sequence connects to the other as farthest neighbor
    }
}
