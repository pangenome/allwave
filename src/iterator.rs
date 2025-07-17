//! Iterator for all-vs-all sequence alignments

use crate::alignment::align_pair;
use crate::knn_graph::{estimate_knn_pair_count, extract_knn_pairs};
use crate::mash::DEFAULT_KMER_SIZE;
use crate::types::{AlignmentParams, AlignmentResult, Sequence, SparsificationStrategy};
use rayon::prelude::*;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// Iterator for all-pair alignments
pub struct AllPairIterator<'a> {
    sequences: &'a [Sequence],
    params: AlignmentParams,
    orientation_params: AlignmentParams,
    exclude_self: bool,
    _sparsification: SparsificationStrategy,
    pair_iter: Box<dyn Iterator<Item = (usize, usize)> + Send + 'a>,
}

impl<'a> AllPairIterator<'a> {
    /// Create iterator for all-vs-all alignments with no sparsification
    pub fn new(sequences: &'a [Sequence], params: AlignmentParams) -> Self {
        Self::with_options(sequences, params, true, SparsificationStrategy::None)
    }

    /// Create iterator with custom options
    pub fn with_options(
        sequences: &'a [Sequence],
        params: AlignmentParams,
        exclude_self: bool,
        sparsification: SparsificationStrategy,
    ) -> Self {
        let n = sequences.len();

        // Generate all pairs
        let mut pairs: Vec<(usize, usize)> = if exclude_self {
            (0..n)
                .flat_map(|i| (0..n).filter(move |&j| i != j).map(move |j| (i, j)))
                .collect()
        } else {
            (0..n).flat_map(|i| (0..n).map(move |j| (i, j))).collect()
        };

        // Apply sparsification
        match &sparsification {
            SparsificationStrategy::None => {}
            SparsificationStrategy::Random(keep_fraction) => {
                pairs = apply_random_sparsification(pairs, *keep_fraction, sequences);
            }
            SparsificationStrategy::Auto => {
                // Use giant component model with 0.95 probability for auto mode
                let keep_fraction = compute_connectivity_probability(n, 0.95);
                pairs = apply_random_sparsification(pairs, keep_fraction, sequences);
            }
            SparsificationStrategy::Connectivity(connectivity_prob) => {
                let keep_fraction = compute_connectivity_probability(n, *connectivity_prob);
                pairs = apply_random_sparsification(pairs, keep_fraction, sequences);
            }
            SparsificationStrategy::NeighborJoining(k_neighbors, random_fraction, kmer_size) => {
                pairs = extract_knn_pairs(
                    sequences,
                    *k_neighbors,
                    *random_fraction,
                    kmer_size.unwrap_or(DEFAULT_KMER_SIZE),
                );
            }
        }

        Self {
            sequences,
            params,
            orientation_params: AlignmentParams::edit_distance(),
            exclude_self,
            _sparsification: sparsification,
            pair_iter: Box::new(pairs.into_iter()),
        }
    }

    /// Set custom orientation detection parameters
    pub fn with_orientation_params(mut self, params: AlignmentParams) -> Self {
        self.orientation_params = params;
        self
    }

    /// Set sparsification strategy
    pub fn with_sparsification(self, strategy: SparsificationStrategy) -> Self {
        // Need to regenerate pairs with new sparsification
        Self::with_options(self.sequences, self.params, self.exclude_self, strategy)
    }

    /// Convert to parallel iterator
    pub fn into_par_iter(self) -> AllPairParallelIterator<'a> {
        // Collect pairs for parallel processing
        let pairs: Vec<_> = self.pair_iter.collect();

        AllPairParallelIterator {
            sequences: self.sequences,
            params: self.params,
            orientation_params: self.orientation_params,
            pairs,
        }
    }

    /// Process alignments with a callback function for streaming results
    pub fn for_each_with_callback<F>(
        self,
        callback: F,
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
    where
        F: Fn(AlignmentResult) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
            + Send
            + Sync,
    {
        self.into_par_iter().for_each_with_callback(callback)
    }

    /// Get the actual number of pairs that will be processed
    pub fn pair_count(&self) -> usize {
        let n = self.sequences.len();
        let base_pairs = if self.exclude_self {
            n * (n - 1)
        } else {
            n * n
        };

        match &self._sparsification {
            SparsificationStrategy::None => base_pairs,
            SparsificationStrategy::Random(keep_fraction) => {
                (base_pairs as f64 * keep_fraction).round() as usize
            }
            SparsificationStrategy::Auto => {
                // Use connectivity model with 0.95 probability for auto mode
                let keep_fraction = compute_connectivity_probability(n, 0.95);
                (base_pairs as f64 * keep_fraction).round() as usize
            }
            SparsificationStrategy::Connectivity(connectivity_prob) => {
                let keep_fraction = compute_connectivity_probability(n, *connectivity_prob);
                (base_pairs as f64 * keep_fraction).round() as usize
            }
            SparsificationStrategy::NeighborJoining(k_neighbors, random_fraction, _kmer_size) => {
                estimate_knn_pair_count(n, *k_neighbors, *random_fraction)
            }
        }
    }
}

impl<'a> Iterator for AllPairIterator<'a> {
    type Item = AlignmentResult;

    fn next(&mut self) -> Option<Self::Item> {
        // Get next pair
        let (i, j) = self.pair_iter.next()?;

        // Perform alignment
        let result = align_pair(
            &self.sequences[i],
            &self.sequences[j],
            i,
            j,
            &self.params,
            &self.orientation_params,
        );

        Some(result)
    }
}

/// Parallel iterator for all-pair alignments
pub struct AllPairParallelIterator<'a> {
    sequences: &'a [Sequence],
    params: AlignmentParams,
    orientation_params: AlignmentParams,
    pairs: Vec<(usize, usize)>,
}

impl<'a> ParallelIterator for AllPairParallelIterator<'a> {
    type Item = AlignmentResult;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: rayon::iter::plumbing::UnindexedConsumer<Self::Item>,
    {
        self.pairs
            .into_par_iter()
            .map(|(i, j)| {
                align_pair(
                    &self.sequences[i],
                    &self.sequences[j],
                    i,
                    j,
                    &self.params,
                    &self.orientation_params,
                )
            })
            .drive_unindexed(consumer)
    }
}

impl<'a> AllPairParallelIterator<'a> {
    /// Process alignments with a callback function for streaming results
    pub fn for_each_with_callback<F>(
        self,
        callback: F,
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
    where
        F: Fn(AlignmentResult) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
            + Send
            + Sync,
    {
        use rayon::prelude::*;
        use std::sync::Mutex;

        let error_holder = Mutex::new(None);

        self.pairs
            .into_par_iter()
            .try_for_each(|(i, j)| {
                let alignment = align_pair(
                    &self.sequences[i],
                    &self.sequences[j],
                    i,
                    j,
                    &self.params,
                    &self.orientation_params,
                );

                match callback(alignment) {
                    Ok(()) => Ok(()),
                    Err(e) => {
                        let mut error_guard = error_holder.lock().unwrap();
                        if error_guard.is_none() {
                            *error_guard = Some(e);
                        }
                        Err(())
                    }
                }
            })
            .map_err(|_| {
                error_holder
                    .into_inner()
                    .unwrap()
                    .unwrap_or_else(|| "Unknown error during parallel processing".into())
            })
    }
}

/// Apply random sparsification to pairs using deterministic hashing of sequence names
fn apply_random_sparsification(
    mut pairs: Vec<(usize, usize)>,
    keep_fraction: f64,
    sequences: &[Sequence],
) -> Vec<(usize, usize)> {
    pairs.retain(|(i, j)| {
        let mut hasher = DefaultHasher::new();

        // Hash the sequence names (IDs) instead of indices
        // This ensures consistent results regardless of sequence order
        let seq_i = &sequences[*i].id;
        let seq_j = &sequences[*j].id;

        // Create directed pair hash: hash(A,B) != hash(B,A)
        // This is important because we're dealing with directed alignments
        // Use simple concatenation to avoid potential hasher bias
        let combined = format!("{seq_i}:{seq_j}");

        combined.hash(&mut hasher);

        let hash = hasher.finish();

        // Convert hash to a value between 0 and 1
        let normalized = (hash as f64) / (u64::MAX as f64);
        normalized < keep_fraction
    });

    pairs
}

/// Compute edge probability for Erdős-Rényi random graph giant component
///
/// For a random graph G(n,p) to maintain a single giant component with probability
/// `connectivity_prob`, we need: p = (log n - log(-log(connectivity_prob)))/n
///
/// This is based on the sharp threshold result:
/// If p = (log n + c)/n, then P_connected → e^(-e^(-c)) as n → ∞
///
/// # Arguments
/// * `n` - Number of nodes (sequences)
/// * `connectivity_prob` - Desired probability that graph has giant component (0 < x < 1)
///
/// # Returns
/// Edge probability p for the random graph
fn compute_connectivity_probability(n: usize, connectivity_prob: f64) -> f64 {
    if n <= 1 {
        return 1.0;
    }

    // Clamp connectivity probability to reasonable range
    let x = connectivity_prob.clamp(0.001, 0.999);

    // For very small n, use simpler heuristics
    if n <= 10 {
        return match n {
            2 => 1.0, // 2 nodes: need the single edge
            3 => 0.8, // 3 nodes: need high probability
            4 => 0.7, // 4 nodes
            5 => 0.6, // 5 nodes
            _ => 0.5, // 6-10 nodes
        };
    }

    let n_f = n as f64;
    let log_n = n_f.ln();

    // c = -log(-log(x))
    let c = -(-x.ln()).ln();

    // p = (log n + c) / n
    let p = (log_n + c) / n_f;

    // For directed graphs (all-vs-all), we need to account for the fact
    // that we're considering n*(n-1) directed edges, not n*(n-1)/2 undirected ones
    // So we can use the same formula directly

    // Ensure p is reasonable (not too small or too large)
    p.clamp(0.001, 1.0)
}
