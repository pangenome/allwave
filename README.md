# AllWave

A high-performance pairwise sequence aligner using bidirectional wavefront alignment (biWFA) with advanced sparsification strategies and intelligent orientation detection. Available as both a command-line tool and a Rust library.

## Features

- **Bidirectional Wavefront Alignment**: Uses the WFA2-lib library for efficient sequence alignment
- **Advanced Sparsification**: Multiple strategies including neighbor-joining trees, giant component connectivity, and stranger-joining
- **Intelligent Orientation Detection**: Mash distance-based orientation detection with edlib fallback option
- **Two-piece Affine Gap Penalties**: Supports complex gap penalty models with dual gap extension penalties
- **PAF Output Format**: Generates standard PAF (Pairwise Alignment Format) output compatible with downstream tools
- **Multi-sequence Support**: Aligns all pairs of sequences in a FASTA file with configurable sparsification
- **Gzipped Input Support**: Automatically detects and handles gzipped FASTA files (.fa.gz)
- **Multithreading Support**: Parallel alignment processing with configurable thread count
- **Library API**: Use AllWave's alignment capabilities in your own Rust projects

## Installation

### Prerequisites

- Rust 1.70 or later
- Git (for fetching dependencies)
- Linux operating system (macOS support is planned for future releases)

#### System Dependencies

AllWave requires the following system libraries to be installed:

**Ubuntu/Debian:**
```bash
sudo apt-get update
# Remove conflicting package if it exists
sudo apt-get remove -y libcurl4-gnutls-dev || true
sudo apt-get install -y libhts-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev
```

**CentOS/RHEL/Fedora:**
```bash
# CentOS/RHEL 8+
sudo dnf install -y htslib-devel zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel

# CentOS/RHEL 7
sudo yum install -y htslib-devel zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel
```

**macOS:**
```bash
# Using Homebrew
brew install htslib zlib bzip2 xz curl openssl
```

These dependencies are required for the faigz-rs library which provides gzipped FASTA file support.

### Building from Source

```bash
git clone https://github.com/pangenome/allwave.git
cd allwave
cargo build --release
```

The binary will be available at `target/release/allwave`.

**Note**: Currently, AllWave only supports Linux due to OpenMP dependencies in the WFA2 library. macOS and Windows support is planned for future releases.

## Usage

### Command Line Tool

```bash
allwave --input <FASTA_FILE> [OPTIONS]
```

#### Core Options

- `-i, --input <FILE>`: Input FASTA file containing sequences to align (required)
- `-o, --output <FILE>`: Output PAF file (default: stdout)
- `-s, --scores <SCORES>`: Alignment scoring parameters (default: "0,5,8,2,24,1")
- `-t, --threads <N>`: Number of threads to use for parallel processing (default: 1)
- `-p, --sparsification <STRATEGY>`: Sparsification strategy (default: "giant:0.99")

#### Specialized Options

- `--no-progress`: Disable progress bar output
- `--edlib-orientation`: Use edlib edit distance for orientation detection instead of mash
- `--mash-matrix`: Output mash distance matrix and exit (for debugging/comparison with mash)

## Sparsification Strategies

Sparsification reduces the number of pairwise alignments computed while maintaining good coverage and connectivity. AllWave supports several strategies:

### Basic Strategies

- **`none`**: Compute all pairwise alignments (n² complexity)
- **`auto`**: Equivalent to `giant:0.95` - automatic giant component connectivity
- **`random:<fraction>`**: Keep random fraction of pairs (0.0 to 1.0)
- **`giant:<probability>`**: Erdős-Rényi random graph model maintaining giant component connectivity

### Advanced Tree-Based Strategies

- **`tree:<near>:<far>:<random>[:<kmer>]`**: Intelligent sampling using k-mer similarity trees
  - `<near>`: Number of k-nearest neighbors (most similar sequences)
  - `<far>`: Number of k-farthest neighbors (most distant sequences - "stranger-joining")
  - `<random>`: Fraction of additional random pairs (0.0 to 1.0)
  - `<kmer>`: K-mer size for similarity computation (default: 15)

### Sparsification Examples

```bash
# All pairwise alignments
allwave -i sequences.fa -p none

# Giant component connectivity (99% probability)
allwave -i sequences.fa -p giant:0.99

# Random 50% of pairs
allwave -i sequences.fa -p random:0.5

# Tree-based: 2 nearest + 1 farthest + 10% random
allwave -i sequences.fa -p tree:2:1:0.1

# Tree-based with custom k-mer size
allwave -i sequences.fa -p tree:3:0:0.2:21

# Stranger-joining only (most distant sequences)
allwave -i sequences.fa -p tree:0:5:0.0
```

### When to Use Each Strategy

- **`none`**: Small datasets (< 100 sequences) where all pairs are needed
- **`giant:0.99`**: Default for most cases - maintains connectivity with high probability
- **`random:0.1-0.5`**: Quick sampling for large datasets
- **`tree:k:0:0.0`**: Neighbor-joining tree (k=1 for pure tree, k>1 for k-nearest neighbor net)
- **`tree:k:j:f`**: Comprehensive sampling combining similarity-based and diversity-based selection

## Orientation Detection

AllWave automatically detects the best orientation (forward/reverse complement) for each sequence pair:

- **Default**: Mash distance-based using strand-specific k-mer sketching (fast, accurate)
- **Fallback**: Use `--edlib-orientation` for edit distance-based detection (slower, sometimes more accurate)

The mash-based method uses k=15 sketches and is significantly faster for longer sequences while maintaining high accuracy.

## Scoring Parameters

The `--scores` parameter controls alignment scoring with a compact interface:

### Parameter Format

**4 parameters**: `match,mismatch,gap_open,gap_ext`
- Used for edit distance or single-piece affine gap penalties
- Example: `0,1,1,1` (edit distance), `0,5,8,2` (affine)

**6 parameters**: `match,mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2`
- Two-piece affine gap model (default)
- Allows different penalties for short vs. long gaps
- Example: `0,5,8,2,24,1` (default genomic parameters)

### Parameter Meanings

- **match**: Score for matching bases (usually 0)
- **mismatch**: Penalty for mismatched bases
- **gap_open**: Penalty for opening a gap (higher = fewer gaps)
- **gap_extend**: Penalty for extending a gap (higher = shorter gaps)
- **gap_open2**: Penalty for opening a second gap type (two-piece model)
- **gap_extend2**: Penalty for extending the second gap type

### Common Scoring Schemes

| Use Case | Parameters | Description |
|----------|------------|-------------|
| **Edit Distance** | `0,1,1,1` | All operations cost 1 |
| **Affine Gaps** | `0,5,8,2` | Simple gap penalty model |
| **Genomic (default)** | `0,5,8,2,24,1` | Optimized for genomic sequences |
| **High Sensitivity** | `0,1,2,1` | Sensitive alignment for similar sequences |
| **Penalize Long Gaps** | `0,5,10,2,50,1` | Heavily penalize long indels |

### Examples by Use Case

```bash
# High sensitivity for similar sequences
allwave -i seqs.fa -s "0,1,2,1"

# Standard genomic alignment (default)
allwave -i seqs.fa -s "0,5,8,2,24,1"

# Penalize long gaps more heavily
allwave -i seqs.fa -s "0,5,10,2,50,1"

# Edit distance for simple comparisons
allwave -i seqs.fa -s "0,1,1,1"
```

## Examples

### Basic Usage

```bash
# Align all sequence pairs with default parameters
allwave -i sequences.fa -o alignments.paf

# Use edit distance scoring
allwave -i sequences.fa -s "0,1,1,1"

# Process gzipped FASTA with multiple threads
allwave -i sequences.fa.gz -o alignments.paf -t 8

# Use sparsification for large datasets
allwave -i large_dataset.fa -p giant:0.99 -t 16
```

### Advanced Sparsification

```bash
# Neighbor-joining tree approach
allwave -i sequences.fa -p tree:1:0:0.0 -o tree_alignments.paf

# Combined similarity and diversity sampling
allwave -i sequences.fa -p tree:3:2:0.1 -o comprehensive_alignments.paf

# Debug mash distances
allwave -i sequences.fa --mash-matrix > distances.tsv
```

### Orientation Detection

```bash
# Use default mash-based orientation (fast)
allwave -i sequences.fa -o alignments.paf

# Use edlib-based orientation (accurate)
allwave -i sequences.fa -o alignments.paf --edlib-orientation
```

## Library API

AllWave can be used as a library in your Rust projects. Add it to your `Cargo.toml`:

```toml
[dependencies]
allwave = { git = "https://github.com/pangenome/allwave.git" }
```

### Basic Usage

```rust
use allwave::{AllPairIterator, Sequence, AlignmentParams, SparsificationStrategy, alignment_to_paf};

// Load your sequences
let sequences = vec![
    Sequence {
        id: "seq1".to_string(),
        seq: b"ACGTACGTACGT".to_vec(),
    },
    Sequence {
        id: "seq2".to_string(), 
        seq: b"ACGTACGTTCGT".to_vec(),
    },
];

// Set alignment parameters (match, mismatch, gap_open, gap_extend, gap2_open, gap2_extend)
let params = AlignmentParams {
    match_score: 0,
    mismatch_penalty: 5,
    gap_open: 8,
    gap_extend: 2,
    gap2_open: Some(24),
    gap2_extend: Some(1),
    max_divergence: None,
};

// Create iterator with sparsification
let aligner = AllPairIterator::with_options(
    &sequences, 
    params,
    true,  // exclude_self
    true,  // use_mash_orientation
    SparsificationStrategy::Connectivity(0.99)  // giant component
);

// Process alignments
for alignment in aligner {
    let paf_line = alignment_to_paf(&alignment, &sequences);
    println!("{}", paf_line);
}
```

### Advanced Sparsification API

```rust
use allwave::SparsificationStrategy;

// Tree-based sampling: 2 nearest + 1 farthest + 10% random
let aligner = AllPairIterator::with_options(
    &sequences,
    params,
    true,
    true,
    SparsificationStrategy::TreeSampling(2, 1, 0.1, Some(15))
);

// Parallel processing with callback
aligner.for_each_with_callback(|alignment| {
    let paf_line = alignment_to_paf(&alignment, &sequences);
    println!("{}", paf_line);
    Ok(())
})?;
```

## Output Format

AllWave produces PAF (Pairwise Alignment Format) output with the following fields:

1. Query sequence name
2. Query sequence length
3. Query start (0-based)
4. Query end
5. Strand (+ or -)
6. Target sequence name
7. Target sequence length
8. Target start (0-based)
9. Target end
10. Number of matches
11. Alignment block length
12. Mapping quality (currently fixed at 60)
13. Tags:
    - `gi:f`: Identity fraction
    - `cg:Z`: CIGAR string

## Performance

AllWave is optimized for speed using:
- Efficient WFA2 algorithm implementation
- Intelligent sparsification strategies
- Mash distance-based orientation detection
- Rust's zero-cost abstractions
- Parallel processing with rayon

Performance scales from small datasets (< 100 sequences) to large pangenome datasets (> 10,000 sequences) through adaptive sparsification.

## Testing

Run the test suite:

```bash
cargo test
```

The test suite includes comprehensive tests for:
- Basic mutations (SNPs and small indels)
- Microsatellite expansions and contractions
- Copy Number Variations (CNVs)
- Complex mutations (all types combined)
- High divergence sequences
- Orientation detection comparison (mash vs edlib)
- Sparsification strategies
- Tree-based sampling methods

## Dependencies

- [bio](https://github.com/rust-bio/rust-bio) - Bioinformatics library for Rust
- [lib_wfa2](https://github.com/ekg/lib_wfa2) - Rust bindings for WFA2-lib
- [edlib_rs](https://github.com/jean-pierreBoth/edlib-rs) - Rust bindings for edlib (fast edit distance)
- [clap](https://github.com/clap-rs/clap) - Command line argument parsing
- [rayon](https://github.com/rayon-rs/rayon) - Data parallelism library
- [faigz-rs](https://github.com/pangenome/faigz-rs) - Fast gzipped FASTA indexing

## License

This project is licensed under the MIT License - see the LICENSE file for details.

Copyright (c) 2025 Erik Garrison

## Acknowledgments

- [WFA2-lib](https://github.com/smarco/WFA2-lib) by Santiago Marco-Sola for the wavefront alignment algorithm
- [edlib](https://github.com/Martinsos/edlib) by Martin Šošić for fast edit distance computation
- [mash](https://github.com/marbl/Mash) by Brian Ondov for the MinHash distance concept
- The Rust bioinformatics community for excellent libraries and tools