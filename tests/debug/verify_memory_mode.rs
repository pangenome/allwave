//! Verify memory mode configuration
//! 
//! Purpose: Tests that the memory mode setting is working correctly with the
//! updated lib_wfa2. This was crucial for verifying that the fix in ekg's fork
//! actually resolved the memory mode issue where settings weren't being applied.
//! 
//! Usage: cargo run --bin verify_memory_mode

use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, HeuristicStrategy, MemoryMode,
};

fn main() {
    println!("=== Testing Memory Mode with Updated lib_wfa2 ===");
    
    // Test 1: Create with default (should be High)
    let aligner1 = AffineWavefronts::with_penalties_affine2p(0, 5, 8, 2, 24, 1);
    println!("Default constructor memory mode: {:?}", aligner1.get_memory_mode());
    
    // Test 2: Create with explicit Ultralow
    let aligner2 = AffineWavefronts::with_penalties_affine2p_and_memory_mode(
        0, 5, 8, 2, 24, 1, MemoryMode::Ultralow
    );
    println!("Ultralow constructor memory mode: {:?}", aligner2.get_memory_mode());
    
    // Test 3: Try each memory mode
    for mode in [MemoryMode::High, MemoryMode::Medium, MemoryMode::Low, MemoryMode::Ultralow] {
        let aligner = AffineWavefronts::with_penalties_affine2p_and_memory_mode(
            0, 5, 8, 2, 24, 1, mode.clone()
        );
        println!("Created with {:?}, got: {:?}", mode, aligner.get_memory_mode());
    }
    
    // Test alignment to make sure it works
    let seq1 = b"ACGTACGTACGT";
    let seq2 = b"ACGTACGTACGT";
    
    let mut aligner = AffineWavefronts::with_penalties_affine2p_and_memory_mode(
        0, 5, 8, 2, 24, 1, MemoryMode::Ultralow
    );
    aligner.set_alignment_scope(AlignmentScope::Alignment);
    aligner.set_alignment_span(AlignmentSpan::End2End);
    aligner.set_heuristic(&HeuristicStrategy::None);
    
    let status = aligner.align(seq1, seq2);
    println!("\nTest alignment with Ultralow mode: {:?}", status);
    if matches!(status, AlignmentStatus::Completed) {
        println!("Score: {}", aligner.score());
    }
}