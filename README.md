# AllWave

A high-performance pairwise sequence aligner using bidirectional wavefront alignment (biWFA) with two-piece affine gap penalties.

## Features

- **Bidirectional Wavefront Alignment**: Uses the WFA2-lib library for efficient sequence alignment
- **Two-piece Affine Gap Penalties**: Supports complex gap penalty models with two different gap extension penalties
- **Automatic Orientation Detection**: Uses edit distance to determine the best orientation (forward/reverse complement) for alignment
- **PAF Output Format**: Generates standard PAF (Pairwise Alignment Format) output compatible with downstream tools
- **Multi-sequence Support**: Aligns all pairs of sequences in a FASTA file

## Installation

### Prerequisites

- Rust 1.70 or later
- Git (for fetching dependencies)

### Building from Source

```bash
git clone https://github.com/yourusername/allwave.git
cd allwave
cargo build --release
```

The binary will be available at `target/release/allwave`.

## Usage

```bash
allwave --input <FASTA_FILE> [--output <PAF_FILE>] [--scores <SCORES>]
```

### Options

- `-i, --input <FILE>`: Input FASTA file containing sequences to align (required)
- `-o, --output <FILE>`: Output PAF file (default: stdout)
- `-s, --scores <SCORES>`: Alignment scoring parameters (default: "0,5,8,2,24,1")

### Examples

```bash
# Align all sequence pairs in a FASTA file with default parameters
allwave --input sequences.fa --output alignments.paf

# Use edit distance scoring
allwave --input sequences.fa --scores "0,1,1,1"

# Use single-piece affine gap penalties
allwave --input sequences.fa --scores "0,5,8,1"

# Use custom two-piece affine gap penalties
allwave --input sequences.fa --scores "0,5,10,2,30,1"
```

## Alignment Parameters

The `--scores` parameter accepts comma-separated values to control alignment scoring.

### Common Scoring Schemes

| Mode | Values | Description |
|------|--------|-------------|
| **Edit Distance** | `0,1,1,1` | All operations cost 1 |
| **Single-Piece Affine** | `0,5,8,1` | Simple gap penalty model |
| **Two-Piece Affine** (default) | `0,5,8,2,24,1` | Complex gap penalty for genomic sequences |

### Parameter Format

**4 parameters**: `match,mismatch,gap_open,gap_ext`
- Used for edit distance (when all non-match operations have equal cost)
- Used for single-piece affine gap penalties

**6 parameters**: `match,mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2`
- Two-piece affine gap model
- Better for sequences with both short and long indels

### Examples by Use Case

```bash
# High sensitivity for similar sequences
allwave --input seqs.fa --scores "0,1,2,1"

# Standard genomic alignment (default)
allwave --input seqs.fa --scores "0,5,8,2,24,1"

# Penalize long gaps more heavily
allwave --input seqs.fa --scores "0,5,10,2,50,1"

# Edit distance for simple comparisons
allwave --input seqs.fa --scores "0,1,1,1"
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

## Testing

Run the test suite:

```bash
cargo test
```

The test suite includes:
- Basic mutations (SNPs and small indels)
- Microsatellite expansions and contractions
- Copy Number Variations (CNVs)
- Complex mutations (all types combined)
- High divergence sequences
- Identical sequences

## Performance

AllWave is optimized for speed using:
- Efficient WFA2 algorithm implementation
- Early termination in orientation detection
- Rust's zero-cost abstractions

## Dependencies

- [bio](https://github.com/rust-bio/rust-bio) - Bioinformatics library for Rust
- [lib_wfa2](https://github.com/AndreaGuarracino/lib_wfa2) - Rust bindings for WFA2-lib
- [clap](https://github.com/clap-rs/clap) - Command line argument parsing

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [WFA2-lib](https://github.com/smarco/WFA2-lib) by Santiago Marco-Sola for the wavefront alignment algorithm
- The Rust bioinformatics community for excellent libraries and tools