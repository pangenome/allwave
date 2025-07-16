# Sparsification Strategies in Allwave

Allwave supports several sparsification strategies to reduce the number of pairwise alignments computed in all-vs-all scenarios. This document explains the mathematical foundations and practical usage of these strategies.

## Available Strategies

### 1. None (`none`)
Compute all pairwise alignments. For n sequences, this results in n² alignments (or n(n-1) if excluding self-alignments).

### 2. Random (`random:<fraction>`)
Keep a random fraction of all pairs. For example, `random:0.5` keeps 50% of all pairs using deterministic hashing based on sequence names.

### 3. Auto (`auto`)
Automatically selects the giant component strategy with 95% probability. This is equivalent to `giant:0.95`.

### 4. Giant Component (`giant:<probability>`)
Uses Erdős-Rényi random graph theory to maintain a single giant (connected) component with the specified probability.

## Giant Component Strategy: Mathematical Foundation

The giant component strategy is based on random graph theory and provides a principled way to sparsify alignments while maintaining graph connectivity.

### Key Concepts

**Giant Component vs. Full Connectivity:**
- **Giant Component**: A single large connected component containing most nodes
- **Full Connectivity**: All nodes are connected in a single component

For a random graph G(n,p) to be **fully connected** (single component containing all n nodes) with probability x, we need a much higher edge probability than for just a giant component.

### Connectivity Threshold

The connectivity threshold is **p = log(n)/n**.

More precisely, for the probability P_connected(n,p) that the graph is connected:

**1. Sharp threshold result:**
- If p = (log n + c)/n, then as n → ∞:
  ```
  P_connected → e^(-e^(-c))
  ```

**2. To achieve connectivity with probability x:**
Set e^(-e^(-c)) = x and solve for c:
- c = -log(-log(x))
- Therefore: **p = (log n - log(-log(x)))/n**

### Examples

For different connectivity probabilities:
- For x = 0.5: p ≈ (log n + 0.3665)/n
- For x = 0.9: p ≈ (log n + 2.2504)/n  
- For x = 0.99: p ≈ (log n + 4.6001)/n
- For x = 0.999: p ≈ (log n + 6.9068)/n

### Key Insights

1. **Scaling**: The edge probability must scale as Θ(log n/n) for connectivity
2. **Threshold**: This is much higher than the Θ(1/n) threshold for a giant component
3. **Precision**: The multiplier of log n/n determines the exact probability through the double exponential function
4. **High probability**: To guarantee connectivity with high probability (say 0.99), you need approximately p ≈ 5 log(n)/n

### Implementation

Allwave implements this using the function `compute_connectivity_probability(n, connectivity_prob)`:

```rust
/// Compute edge probability for Erdős-Rényi random graph giant component
/// 
/// For a random graph G(n,p) to maintain a single giant component with probability 
/// `connectivity_prob`, we need: p = (log n - log(-log(connectivity_prob)))/n
fn compute_connectivity_probability(n: usize, connectivity_prob: f64) -> f64 {
    // For large n, use the theoretical formula
    let log_n = (n as f64).ln();
    let c = -(-connectivity_prob.ln()).ln();
    let p = (log_n + c) / (n as f64);
    
    // Clamp to reasonable bounds
    p.max(0.0).min(1.0)
}
```

### Finite Graph Corrections

For finite graphs, a more accurate approximation includes:
```
P_connected ≈ exp(-n·p^(n-1)·e^(-p·n(n-1)/2))
```
when p is near the threshold.

## Usage Examples

### Command Line
```bash
# Use giant component with 99% probability (default)
allwave -i sequences.fa -p giant:0.99

# Use giant component with 95% probability
allwave -i sequences.fa -p giant:0.95

# Use random 50% sampling
allwave -i sequences.fa -p random:0.5

# Use automatic strategy (equivalent to giant:0.95)
allwave -i sequences.fa -p auto

# Compute all pairs
allwave -i sequences.fa -p none
```

### API Usage
```rust
use allwave::{AllPairIterator, SparsificationStrategy};

// Giant component with 99% probability
let strategy = SparsificationStrategy::Connectivity(0.99);
let aligner = AllPairIterator::with_options(&sequences, params, true, strategy);

// Random sampling
let strategy = SparsificationStrategy::Random(0.5);
let aligner = AllPairIterator::with_options(&sequences, params, true, strategy);
```

## Practical Considerations

1. **Default Strategy**: `giant:0.99` provides a good balance between connectivity and computational efficiency
2. **Large Datasets**: For very large datasets (n > 10,000), consider lower probabilities like `giant:0.90`
3. **Small Datasets**: For small datasets (n < 100), `none` or `random:0.8` may be more appropriate
4. **Deterministic**: All strategies use deterministic hashing based on sequence names, ensuring reproducible results

## References

- Erdős, P. and Rényi, A. (1960). "On the evolution of random graphs"
- Bollobás, B. (2001). "Random Graphs" (2nd ed.)
- Frieze, A. and Karoński, M. (2016). "Introduction to Random Graphs"