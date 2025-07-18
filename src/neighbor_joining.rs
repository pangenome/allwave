//! Neighbor-joining phylogenetic tree construction
//!
//! This module implements the neighbor-joining algorithm for constructing phylogenetic trees
//! from distance matrices. The algorithm is used to create trees that guide intelligent
//! sequence alignment sampling.

use std::collections::HashMap;

/// A node in the neighbor-joining tree
#[derive(Debug, Clone)]
pub struct TreeNode {
    /// Node identifier
    pub id: usize,
    /// If this is a leaf, the original sequence index
    pub sequence_index: Option<usize>,
    /// Left child (None for leaves)
    pub left: Option<Box<TreeNode>>,
    /// Right child (None for leaves)
    pub right: Option<Box<TreeNode>>,
    /// Branch length to parent
    pub branch_length: f64,
}

impl TreeNode {
    /// Create a new leaf node
    pub fn leaf(id: usize, sequence_index: usize) -> Self {
        Self {
            id,
            sequence_index: Some(sequence_index),
            left: None,
            right: None,
            branch_length: 0.0,
        }
    }

    /// Create a new internal node
    pub fn internal(id: usize, left: TreeNode, right: TreeNode) -> Self {
        Self {
            id,
            sequence_index: None,
            left: Some(Box::new(left)),
            right: Some(Box::new(right)),
            branch_length: 0.0,
        }
    }

    /// Check if this node is a leaf
    pub fn is_leaf(&self) -> bool {
        self.left.is_none() && self.right.is_none()
    }

    /// Get all leaf nodes (sequences) in this subtree
    pub fn get_leaves(&self) -> Vec<usize> {
        if let Some(seq_idx) = self.sequence_index {
            vec![seq_idx]
        } else {
            let mut leaves = Vec::new();
            if let Some(left) = &self.left {
                leaves.extend(left.get_leaves());
            }
            if let Some(right) = &self.right {
                leaves.extend(right.get_leaves());
            }
            leaves
        }
    }

    /// Get all edges in the tree as (parent_sequences, child_sequences) pairs
    pub fn get_edges(&self) -> Vec<(Vec<usize>, Vec<usize>)> {
        let mut edges = Vec::new();

        if let (Some(left), Some(right)) = (&self.left, &self.right) {
            let left_leaves = left.get_leaves();
            let right_leaves = right.get_leaves();
            let parent_leaves = self.get_leaves();

            // Add edge from this node to left child
            edges.push((parent_leaves.clone(), left_leaves.clone()));
            // Add edge from this node to right child
            edges.push((parent_leaves, right_leaves));

            // Recursively get edges from children
            edges.extend(left.get_edges());
            edges.extend(right.get_edges());
        }

        edges
    }
}

/// Construct neighbor-joining tree from distance matrix
pub fn neighbor_joining(distance_matrix: &[Vec<f64>]) -> Option<TreeNode> {
    let n = distance_matrix.len();
    if n < 2 {
        return None;
    }

    if n == 2 {
        // Special case: only two sequences
        let mut left = TreeNode::leaf(0, 0);
        let mut right = TreeNode::leaf(1, 1);
        left.branch_length = distance_matrix[0][1] / 2.0;
        right.branch_length = distance_matrix[0][1] / 2.0;
        return Some(TreeNode::internal(2, left, right));
    }

    // Initialize active nodes (all sequences start as leaves)
    let mut active_nodes: HashMap<usize, TreeNode> =
        (0..n).map(|i| (i, TreeNode::leaf(i, i))).collect();

    // Working distance matrix (will be modified)
    #[allow(clippy::map_clone)]
    let mut distances = distance_matrix
        .iter()
        .map(|row| row.clone())
        .collect::<Vec<_>>();

    let mut next_node_id = n;

    // Main neighbor-joining loop
    while active_nodes.len() > 2 {
        let active_indices: Vec<usize> = active_nodes.keys().cloned().collect();
        let m = active_indices.len();

        // Compute Q matrix for neighbor-joining criterion
        let mut q_matrix = vec![vec![0.0; m]; m];
        let mut row_sums = vec![0.0; m];

        // Calculate row sums
        for i in 0..m {
            for j in 0..m {
                if i != j {
                    row_sums[i] += distances[active_indices[i]][active_indices[j]];
                }
            }
        }

        // Calculate Q matrix
        for i in 0..m {
            for j in i + 1..m {
                let d_ij = distances[active_indices[i]][active_indices[j]];
                q_matrix[i][j] = (m as f64 - 2.0) * d_ij - row_sums[i] - row_sums[j];
                q_matrix[j][i] = q_matrix[i][j];
            }
        }

        // Find minimum Q value
        let mut min_q = f64::INFINITY;
        let mut min_i = 0;
        let mut min_j = 1;

        #[allow(clippy::needless_range_loop)]
        for i in 0..m {
            for j in i + 1..m {
                if q_matrix[i][j] < min_q {
                    min_q = q_matrix[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        let node_i_idx = active_indices[min_i];
        let node_j_idx = active_indices[min_j];

        // Calculate branch lengths
        let d_ij = distances[node_i_idx][node_j_idx];
        let branch_i = d_ij / 2.0 + (row_sums[min_i] - row_sums[min_j]) / (2.0 * (m as f64 - 2.0));
        let branch_j = d_ij - branch_i;

        // Create new internal node
        let mut node_i = active_nodes.remove(&node_i_idx).unwrap();
        let mut node_j = active_nodes.remove(&node_j_idx).unwrap();

        node_i.branch_length = branch_i.max(0.0);
        node_j.branch_length = branch_j.max(0.0);

        let new_node = TreeNode::internal(next_node_id, node_i, node_j);

        // Update distance matrix
        // Calculate distances from new node to all remaining nodes
        let mut new_distances = vec![0.0; distances.len().max(next_node_id + 1)];
        for &k_idx in &active_indices {
            if k_idx != node_i_idx && k_idx != node_j_idx {
                let d_ik = distances[node_i_idx][k_idx];
                let d_jk = distances[node_j_idx][k_idx];
                new_distances[k_idx] = (d_ik + d_jk - d_ij) / 2.0;
            }
        }

        // Extend the matrix if needed
        let new_size = distances.len().max(next_node_id + 1);
        for row in &mut distances {
            row.resize(new_size, 0.0);
        }
        distances.resize(new_size, vec![0.0; new_size]);

        // Update the distances matrix
        for &k_idx in &active_indices {
            if k_idx != node_i_idx && k_idx != node_j_idx {
                distances[next_node_id][k_idx] = new_distances[k_idx];
                distances[k_idx][next_node_id] = new_distances[k_idx];
            }
        }

        active_nodes.insert(next_node_id, new_node);
        next_node_id += 1;
    }

    // Handle the final two nodes
    if active_nodes.len() == 2 {
        let indices: Vec<usize> = active_nodes.keys().cloned().collect();
        let mut node1 = active_nodes.remove(&indices[0]).unwrap();
        let mut node2 = active_nodes.remove(&indices[1]).unwrap();

        let final_distance = if indices[0] < distances.len() && indices[1] < distances[0].len() {
            distances[indices[0]][indices[1]]
        } else {
            1.0 // Default distance
        };

        node1.branch_length = final_distance / 2.0;
        node2.branch_length = final_distance / 2.0;

        Some(TreeNode::internal(next_node_id, node1, node2))
    } else {
        None
    }
}

/// Extract all sequence pairs that should be aligned based on the tree structure
pub fn extract_tree_pairs(tree: &TreeNode, random_fraction: f64) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();

    // Get all edges in the tree
    let edges = tree.get_edges();

    // For each edge, sample pairs between the two groups
    for (group1, group2) in edges {
        for &i in &group1 {
            for &j in &group2 {
                if i != j {
                    // Include this pair with the specified probability
                    if sample_with_probability(i, j, random_fraction) {
                        pairs.push((i, j));
                    }
                }
            }
        }
    }

    // Remove duplicates
    pairs.sort_unstable();
    pairs.dedup();

    pairs
}

/// Sample a pair with given probability using deterministic hashing
fn sample_with_probability(i: usize, j: usize, probability: f64) -> bool {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut hasher = DefaultHasher::new();
    (i, j).hash(&mut hasher);
    let hash = hasher.finish();
    let normalized = (hash as f64) / (u64::MAX as f64);
    normalized < probability
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leaf_node() {
        let node = TreeNode::leaf(0, 5);
        assert!(node.is_leaf());
        assert_eq!(node.sequence_index, Some(5));
        assert_eq!(node.get_leaves(), vec![5]);
    }

    #[test]
    fn test_internal_node() {
        let left = TreeNode::leaf(0, 1);
        let right = TreeNode::leaf(1, 2);
        let internal = TreeNode::internal(2, left, right);

        assert!(!internal.is_leaf());
        assert_eq!(internal.sequence_index, None);

        let mut leaves = internal.get_leaves();
        leaves.sort();
        assert_eq!(leaves, vec![1, 2]);
    }

    #[test]
    fn test_simple_neighbor_joining() {
        // Simple 3x3 distance matrix
        let distances = vec![
            vec![0.0, 1.0, 3.0],
            vec![1.0, 0.0, 2.0],
            vec![3.0, 2.0, 0.0],
        ];

        let tree = neighbor_joining(&distances);
        assert!(tree.is_some());

        let tree = tree.unwrap();
        assert!(!tree.is_leaf());

        let leaves = tree.get_leaves();
        assert_eq!(leaves.len(), 3);
    }

    #[test]
    fn test_two_sequence_case() {
        let distances = vec![vec![0.0, 1.5], vec![1.5, 0.0]];

        let tree = neighbor_joining(&distances);
        assert!(tree.is_some());

        let tree = tree.unwrap();
        let leaves = tree.get_leaves();
        assert_eq!(leaves.len(), 2);
    }

    #[test]
    fn test_extract_tree_pairs() {
        let left = TreeNode::leaf(0, 0);
        let right = TreeNode::leaf(1, 1);
        let tree = TreeNode::internal(2, left, right);

        let pairs = extract_tree_pairs(&tree, 1.0); // Sample all pairs
        assert!(!pairs.is_empty());
    }
}
