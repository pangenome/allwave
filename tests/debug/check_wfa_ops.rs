//! Check WFA2 CIGAR operation codes
//!
//! Purpose: Tests what operation codes (M, X, I, D) WFA2 outputs for different
//! alignment scenarios. Used to verify that WFA2 correctly distinguishes between
//! matches (M) and mismatches (X) in its CIGAR output.
//!
//! Usage: cargo run --bin check_wfa_ops

use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, HeuristicStrategy,
};

fn main() {
    // Test sequences with known differences
    let query = b"ACGTACGTACGT";
    let reference = b"ACGTACGTTCGT"; // One mismatch at position 8 (A->T)

    println!("Query:     {}", std::str::from_utf8(query).unwrap());
    println!("Reference: {}", std::str::from_utf8(reference).unwrap());
    println!("           --------^--- (mismatch here)");

    let mut aligner = AffineWavefronts::default();
    aligner.set_alignment_scope(AlignmentScope::Alignment);
    aligner.set_alignment_span(AlignmentSpan::End2End);
    aligner.set_heuristic(&HeuristicStrategy::None);

    let status = aligner.align(query, reference);
    if matches!(status, AlignmentStatus::Completed) {
        let cigar = aligner.cigar();
        print!("\nCIGAR bytes: ");
        for &op in cigar {
            print!("{}", op as char);
        }
        println!();

        print!("As hex: ");
        for &op in cigar {
            print!("{:02x} ", op);
        }
        println!();

        // Count operation types
        let mut op_counts = std::collections::HashMap::new();
        for &op in cigar {
            *op_counts.entry(op).or_insert(0) += 1;
        }

        println!("\nOperation counts:");
        for (op, count) in &op_counts {
            println!("  {} (0x{:02x}): {}", *op as char, op, count);
        }
    }
}
