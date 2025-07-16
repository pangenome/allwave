//! Test WFA2 parameter order
//!
//! Purpose: Verifies the correct parameter order for WFA2 alignment functions.
//! The Rust wrapper uses (query, reference) order, which maps to (pattern, text)
//! in WFA2 terminology. This test was created to ensure we weren't swapping
//! parameters incorrectly.
//!
//! Usage: cargo run --bin test_wfa_order

use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode};

fn main() {
    let seq1 = b"ACGTACGTACGT"; // 12 bases - will be pattern/query
    let seq2 = b"ACGTACGTAC"; // 10 bases - will be text/reference

    println!("Testing WFA2 alignment order\n");
    println!(
        "seq1 (pattern): {} ({} bases)",
        std::str::from_utf8(seq1).unwrap(),
        seq1.len()
    );
    println!(
        "seq2 (text):    {} ({} bases)",
        std::str::from_utf8(seq2).unwrap(),
        seq2.len()
    );

    let wf = AffineWavefronts::with_penalties_and_memory_mode(0, 5, 8, 2, MemoryMode::Ultralow);

    println!("\n=== Align seq1 vs seq2 ===");
    let status1 = wf.align(seq1, seq2);
    println!("Status: {:?}", status1);

    let cigar1 = wf.cigar();
    print!("CIGAR: ");
    for &op in cigar1 {
        print!("{}", op as char);
    }
    println!();

    // Count operations
    let mut ops = [0; 256];
    for &op in cigar1 {
        ops[op as usize] += 1;
    }

    println!("Operations:");
    for (i, &count) in ops.iter().enumerate() {
        if count > 0 {
            println!("  {} ({}): {}", i as u8 as char, i, count);
        }
    }

    println!("\n=== Align seq2 vs seq1 (swapped) ===");
    let status2 = wf.align(seq2, seq1);
    println!("Status: {:?}", status2);

    let cigar2 = wf.cigar();
    print!("CIGAR: ");
    for &op in cigar2 {
        print!("{}", op as char);
    }
    println!();
}
