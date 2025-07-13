//! Debug CIGAR interpretation for I/D operations
//! 
//! Purpose: Tests our understanding of WFA2's CIGAR convention where I and D
//! operations are swapped compared to standard CIGAR format. This was crucial
//! for fixing the pafcheck validation failures.
//! 
//! Usage: cargo run --bin debug_cigar

use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, HeuristicStrategy, MemoryMode,
};

fn main() {
    // Create two sequences of different lengths
    let seq1 = b"ACGTACGTACGT"; // 12 bases
    let seq2 = b"ACGTACGTAC";   // 10 bases

    println!("Seq1: {} (len={})", std::str::from_utf8(seq1).unwrap(), seq1.len());
    println!("Seq2: {} (len={})", std::str::from_utf8(seq2).unwrap(), seq2.len());

    // Create aligner
    let mut aligner = AffineWavefronts::with_penalties_affine2p_and_memory_mode(
        0, 5, 8, 2, 24, 1, MemoryMode::Ultralow
    );
    aligner.set_alignment_scope(AlignmentScope::Alignment);
    aligner.set_alignment_span(AlignmentSpan::End2End);
    aligner.set_heuristic(&HeuristicStrategy::None);

    // Align seq1 (query) vs seq2 (reference)
    println!("\n=== Aligning seq1 (query) vs seq2 (reference) ===");
    let status = aligner.align(seq1, seq2);
    println!("Status: {:?}", status);

    if matches!(status, AlignmentStatus::Completed) {
        let cigar = aligner.cigar();
        print!("CIGAR: ");
        for &op in cigar {
            print!("{}", op as char);
        }
        println!(" (len={})", cigar.len());

        // Count operations
        let mut q_pos = 0;
        let mut r_pos = 0;
        for &op in cigar {
            match op {
                b'M' | b'=' | b'X' => {
                    q_pos += 1;
                    r_pos += 1;
                }
                b'I' => q_pos += 1,
                b'D' => r_pos += 1,
                _ => {}
            }
        }
        println!("Query consumed: {}/{}", q_pos, seq1.len());
        println!("Reference consumed: {}/{}", r_pos, seq2.len());
    }

    // Try the opposite
    println!("\n=== Aligning seq2 (query) vs seq1 (reference) ===");
    let status2 = aligner.align(seq2, seq1);
    println!("Status: {:?}", status2);

    if matches!(status2, AlignmentStatus::Completed) {
        let cigar = aligner.cigar();
        print!("CIGAR: ");
        for &op in cigar {
            print!("{}", op as char);
        }
        println!(" (len={})", cigar.len());

        // Count operations
        let mut q_pos = 0;
        let mut r_pos = 0;
        for &op in cigar {
            match op {
                b'M' | b'=' | b'X' => {
                    q_pos += 1;
                    r_pos += 1;
                }
                b'I' => q_pos += 1,
                b'D' => r_pos += 1,
                _ => {}
            }
        }
        println!("Query consumed: {}/{}", q_pos, seq2.len());
        println!("Reference consumed: {}/{}", r_pos, seq1.len());
    }
}