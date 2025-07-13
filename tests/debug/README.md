# Debug and Testing Tools

This directory contains debugging and testing utilities that were created during the development of allwave to diagnose and fix various issues with WFA2 integration.

## Building Debug Tools

These tools are not built by default to speed up compilation. To build a specific tool:

```bash
cargo build --release --features debug-bins --bin <tool_name>
```

To run a tool:

```bash
cargo run --release --features debug-bins --bin <tool_name> [args]
```

## Available Tools

### check_wfa_ops
Tests WFA2 CIGAR operation codes output. Verifies that WFA2 correctly distinguishes between matches (M) and mismatches (X).

### debug_align
Simple sequence comparison tool for checking character-by-character differences between two sequences.

### debug_cigar
Tests CIGAR interpretation, specifically our fix for WFA2's opposite I/D convention.

### debug_problem
Investigates specific problematic alignments from x.fa that were failing pafcheck validation.

### debug_wfa
General-purpose WFA2 debugging tool that can analyze alignments from any FASTA file.

### test_allwave
Comprehensive integration test suite with synthetic sequences and various mutation types.

### test_cigar_interpretation
Specifically tests the I/D operation swap fix.

### test_wfa_order
Verifies correct parameter order for WFA2 alignment functions.

### verify_memory_mode
Tests that memory mode settings are working correctly with lib_wfa2.

## Historical Context

These tools were instrumental in fixing several issues:
1. WFA2's opposite convention for I/D operations in CIGAR strings
2. Memory mode settings not being applied correctly
3. CIGAR strings extending beyond sequence boundaries
4. Understanding WFA2's parameter order and alignment behavior