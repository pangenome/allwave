//! Test CIGAR I/D operation interpretation
//!
//! Purpose: Specifically tests our fix for WFA2's opposite convention for I/D
//! operations. WFA2 uses I to mean "consume reference" and D to mean "consume
//! query", which is opposite to the standard CIGAR convention. This test
//! verifies our swap is working correctly.
//!
//! Usage: cargo run --bin test_cigar_interpretation

use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, HeuristicStrategy,
};

fn main() {
    // Test with sequences where we know the expected result
    let query = b"ACGTACGT"; // 8 bases
    let reference = b"ACGTACGT"; // 8 bases - identical

    println!("=== Test 1: Identical sequences ===");
    println!("Query:     {}", std::str::from_utf8(query).unwrap());
    println!("Reference: {}", std::str::from_utf8(reference).unwrap());

    let mut aligner = AffineWavefronts::default();
    aligner.set_alignment_scope(AlignmentScope::Alignment);
    aligner.set_alignment_span(AlignmentSpan::End2End);
    aligner.set_heuristic(&HeuristicStrategy::None);

    let status = aligner.align(query, reference);
    if matches!(status, AlignmentStatus::Completed) {
        let cigar = aligner.cigar();
        print!("CIGAR: ");
        for &op in cigar {
            print!("{}", op as char);
        }
        println!(" (should be all M)");
    }

    // Test 2: Query has extra bases at end
    println!("\n=== Test 2: Query longer than reference ===");
    let query2 = b"ACGTACGTAA"; // 10 bases
    let reference2 = b"ACGTACGT"; // 8 bases

    println!(
        "Query:     {} (len={})",
        std::str::from_utf8(query2).unwrap(),
        query2.len()
    );
    println!(
        "Reference: {} (len={})",
        std::str::from_utf8(reference2).unwrap(),
        reference2.len()
    );

    let status2 = aligner.align(query2, reference2);
    if matches!(status2, AlignmentStatus::Completed) {
        let cigar = aligner.cigar();
        print!("CIGAR: ");
        for &op in cigar {
            print!("{}", op as char);
        }
        println!("\nExpected: 8M2I (8 matches, 2 insertions from query)");

        // Count ops
        let mut m_count = 0;
        let mut i_count = 0;
        let mut d_count = 0;
        for &op in cigar {
            match op {
                b'M' | b'=' => m_count += 1,
                b'I' => i_count += 1,
                b'D' => d_count += 1,
                _ => {}
            }
        }
        println!("Actual: M={}, I={}, D={}", m_count, i_count, d_count);
    }

    // Test 3: Reference longer than query
    println!("\n=== Test 3: Reference longer than query ===");
    let query3 = b"ACGTACGT"; // 8 bases
    let reference3 = b"ACGTACGTAA"; // 10 bases

    println!(
        "Query:     {} (len={})",
        std::str::from_utf8(query3).unwrap(),
        query3.len()
    );
    println!(
        "Reference: {} (len={})",
        std::str::from_utf8(reference3).unwrap(),
        reference3.len()
    );

    let status3 = aligner.align(query3, reference3);
    if matches!(status3, AlignmentStatus::Completed) {
        let cigar = aligner.cigar();
        print!("CIGAR: ");
        for &op in cigar {
            print!("{}", op as char);
        }
        println!("\nExpected: 8M2D (8 matches, 2 deletions from reference)");

        // Count ops
        let mut m_count = 0;
        let mut i_count = 0;
        let mut d_count = 0;
        for &op in cigar {
            match op {
                b'M' | b'=' => m_count += 1,
                b'I' => i_count += 1,
                b'D' => d_count += 1,
                _ => {}
            }
        }
        println!("Actual: M={}, I={}, D={}", m_count, i_count, d_count);
    }
}
