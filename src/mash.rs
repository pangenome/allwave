//! Mash distance computation using MinHash k-mer sketching
//!
//! This module provides functionality to compute pairwise mash distances between sequences
//! using k-mer sketching with MinHash. The implementation uses k=15 by default for good
//! resolution while maintaining reasonable performance.

use crate::types::Sequence;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};

/// Default k-mer size for mash distance computation
pub const DEFAULT_KMER_SIZE: usize = 15;

/// Default sketch size (number of minimizers to keep)
pub const DEFAULT_SKETCH_SIZE: usize = 1000;

/// A k-mer sketch using MinHash
#[derive(Debug, Clone)]
pub struct KmerSketch {
    /// Minimizers (smallest hash values) for this sequence
    pub minimizers: Vec<u64>,
    /// K-mer size used
    pub k: usize,
    /// Original sequence length
    pub length: usize,
}

impl KmerSketch {
    /// Create a new sketch from a sequence
    pub fn from_sequence(sequence: &[u8], k: usize, sketch_size: usize) -> Self {
        let minimizers = sketch_sequence(sequence, k, sketch_size);
        Self {
            minimizers,
            k,
            length: sequence.len(),
        }
    }

    /// Compute Jaccard index with another sketch
    pub fn jaccard(&self, other: &KmerSketch) -> f64 {
        if self.k != other.k {
            return 0.0;
        }

        let set1: HashSet<_> = self.minimizers.iter().collect();
        let set2: HashSet<_> = other.minimizers.iter().collect();

        let intersection_size = set1.intersection(&set2).count();
        let union_size = set1.union(&set2).count();

        if union_size == 0 {
            0.0
        } else {
            intersection_size as f64 / union_size as f64
        }
    }

    /// Compute mash distance from Jaccard index
    pub fn mash_distance(&self, other: &KmerSketch) -> f64 {
        let jaccard = self.jaccard(other);
        if jaccard <= 0.0 {
            1.0 // Maximum distance
        } else {
            // Mash distance: -1/k * ln(2*J/(1+J))
            // where J is Jaccard index
            let k = self.k as f64;
            let ratio = (2.0 * jaccard) / (1.0 + jaccard);
            if ratio <= 0.0 {
                1.0
            } else {
                (-1.0 / k) * ratio.ln()
            }
        }
    }
}

/// Create MinHash sketch from sequence
fn sketch_sequence(sequence: &[u8], k: usize, sketch_size: usize) -> Vec<u64> {
    if sequence.len() < k {
        return Vec::new();
    }

    let mut hashes = Vec::new();

    // Extract all k-mers and hash them
    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];

        // Skip k-mers containing non-ACGT characters
        if kmer.iter().any(|&b| !is_dna_base(b)) {
            continue;
        }

        // Hash the k-mer (considering both orientations)
        let hash_fwd = hash_kmer(kmer);
        let hash_rev = hash_kmer(&reverse_complement_kmer(kmer));

        // Use canonical k-mer (lexicographically smaller)
        let canonical_hash = hash_fwd.min(hash_rev);
        hashes.push(canonical_hash);
    }

    // Sort and take the smallest hashes (MinHash)
    hashes.sort_unstable();
    hashes.truncate(sketch_size);
    hashes
}

/// Hash a k-mer using a simple polynomial rolling hash
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    kmer.hash(&mut hasher);
    hasher.finish()
}

/// Check if byte represents a DNA base
fn is_dna_base(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Compute reverse complement of a k-mer
fn reverse_complement_kmer(kmer: &[u8]) -> Vec<u8> {
    kmer.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b, // Keep non-DNA characters as-is
        })
        .collect()
}

/// Compute all-vs-all mash distances between sequences
pub fn compute_distance_matrix(sequences: &[Sequence]) -> Vec<Vec<f64>> {
    compute_distance_matrix_with_params(sequences, DEFAULT_KMER_SIZE, DEFAULT_SKETCH_SIZE)
}

/// Compute all-vs-all mash distances with custom parameters
pub fn compute_distance_matrix_with_params(
    sequences: &[Sequence],
    k: usize,
    sketch_size: usize,
) -> Vec<Vec<f64>> {
    let n = sequences.len();
    let mut matrix = vec![vec![0.0; n]; n];

    // Create sketches for all sequences
    let sketches: Vec<_> = sequences
        .iter()
        .map(|seq| KmerSketch::from_sequence(&seq.seq, k, sketch_size))
        .collect();

    // Compute pairwise distances
    for i in 0..n {
        for j in i + 1..n {
            let distance = sketches[i].mash_distance(&sketches[j]);
            matrix[i][j] = distance;
            matrix[j][i] = distance; // Symmetric
        }
    }

    matrix
}

/// Print distance matrix in a tab-separated format
pub fn print_distance_matrix(sequences: &[Sequence], matrix: &[Vec<f64>]) {
    // Print header row
    print!("sequence");
    for seq in sequences {
        print!("\t{}", seq.id);
    }
    println!();

    // Print matrix rows
    for (i, row) in matrix.iter().enumerate() {
        print!("{}", sequences[i].id);
        for &distance in row {
            print!("\t{distance:.6}");
        }
        println!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_sketch() {
        let seq = b"ATCGATCGATCG";
        let sketch = KmerSketch::from_sequence(seq, 4, 10);
        assert!(!sketch.minimizers.is_empty());
        assert_eq!(sketch.k, 4);
        assert_eq!(sketch.length, seq.len());
    }

    #[test]
    fn test_jaccard_identical() {
        let seq = b"ATCGATCGATCG";
        let sketch1 = KmerSketch::from_sequence(seq, 4, 10);
        let sketch2 = KmerSketch::from_sequence(seq, 4, 10);

        let jaccard = sketch1.jaccard(&sketch2);
        assert!((jaccard - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mash_distance_identical() {
        let seq = b"ATCGATCGATCG";
        let sketch1 = KmerSketch::from_sequence(seq, 4, 10);
        let sketch2 = KmerSketch::from_sequence(seq, 4, 10);

        let distance = sketch1.mash_distance(&sketch2);
        assert!(distance < 1e-10); // Should be very close to 0
    }

    #[test]
    fn test_reverse_complement() {
        let kmer = b"ATCG";
        let rc = reverse_complement_kmer(kmer);
        assert_eq!(rc, b"CGAT");
    }

    #[test]
    fn test_distance_matrix() {
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq2".to_string(),
                seq: b"ATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq3".to_string(),
                seq: b"GGGGGGGGGGGGGGGG".to_vec(),
            },
        ];

        let matrix = compute_distance_matrix(&sequences);
        assert_eq!(matrix.len(), 3);
        assert_eq!(matrix[0].len(), 3);

        // Diagonal should be 0
        assert!(matrix[0][0] < 1e-6);
        assert!(matrix[1][1] < 1e-6);
        assert!(matrix[2][2] < 1e-6);

        // Identical sequences should have distance ~0 (allow some tolerance for floating point)
        assert!(matrix[0][1] < 1e-6);
        assert!(matrix[1][0] < 1e-6);

        // Different sequences should have non-zero distance
        assert!(matrix[0][2] > 0.0);
        assert!(matrix[2][0] > 0.0);
    }
}
