//! Iterator for all-vs-all sequence alignments

use crate::alignment::align_pair;
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
    /// Create iterator for all-vs-all alignments
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
            (0..n).flat_map(|i| (0..n).filter(move |&j| i != j).map(move |j| (i, j))).collect()
        } else {
            (0..n).flat_map(|i| (0..n).map(move |j| (i, j))).collect()
        };
        
        // Apply sparsification
        match &sparsification {
            SparsificationStrategy::None => {},
            SparsificationStrategy::Random(keep_fraction) => {
                pairs = apply_random_sparsification(pairs, *keep_fraction);
            },
            SparsificationStrategy::Auto => {
                let keep_fraction = compute_auto_sparsification(n);
                pairs = apply_random_sparsification(pairs, keep_fraction);
            },
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
                let keep_fraction = compute_auto_sparsification(n);
                (base_pairs as f64 * keep_fraction).round() as usize
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

/// Apply random sparsification to pairs using deterministic hashing
fn apply_random_sparsification(
    mut pairs: Vec<(usize, usize)>, 
    keep_fraction: f64
) -> Vec<(usize, usize)> {
    pairs.retain(|(i, j)| {
        let mut hasher = DefaultHasher::new();
        i.hash(&mut hasher);
        j.hash(&mut hasher);
        let hash = hasher.finish();
        
        // Convert hash to a value between 0 and 1
        let normalized = (hash as f64) / (u64::MAX as f64);
        normalized < keep_fraction
    });
    
    pairs
}

/// Compute automatic sparsification factor based on sequence count
fn compute_auto_sparsification(n: usize) -> f64 {
    match n {
        0..=50 => 1.0,          // No sparsification for small datasets
        51..=200 => 0.8,        // Keep 80% of pairs
        201..=500 => 0.5,       // Keep 50% of pairs
        501..=1000 => 0.3,      // Keep 30% of pairs
        1001..=5000 => 0.1,     // Keep 10% of pairs
        _ => 0.05,              // Keep 5% for very large datasets
    }
}