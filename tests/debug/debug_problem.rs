//! Debug specific problematic alignments
//! 
//! Purpose: Investigates alignment issues with specific sequences from x.fa that
//! were failing pafcheck validation. This tool helped identify that WFA2 was
//! producing CIGARs that extended beyond sequence boundaries.
//! 
//! Usage: cargo run --bin debug_problem
//! Expects: x.fa file in current directory

use bio::io::fasta;
use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, HeuristicStrategy, MemoryMode,
};
use std::fs::File;
use std::io::BufReader;

fn main() {
    let file = File::open("x.fa").expect("Cannot open file");
    let reader = BufReader::new(file);
    let fasta_reader = fasta::Reader::new(reader);
    let mut sequences = Vec::new();

    for result in fasta_reader.records() {
        let record = result.expect("Error reading record");
        sequences.push((record.id().to_string(), record.seq().to_vec()));
    }

    // Get the two problematic sequences
    let (name1, seq1) = &sequences[0];
    let (name2, seq2) = &sequences[1];

    println!("Seq1: {} (len={})", name1, seq1.len());
    println!("Seq2: {} (len={})", name2, seq2.len());

    // Create aligner
    let mut aligner = AffineWavefronts::with_penalties_affine2p_and_memory_mode(
        0, 5, 8, 2, 24, 1, MemoryMode::Ultralow
    );
    aligner.set_alignment_scope(AlignmentScope::Alignment);
    aligner.set_alignment_span(AlignmentSpan::End2End);
    aligner.set_heuristic(&HeuristicStrategy::None);

    // Align
    let status = aligner.align(seq1, seq2);
    println!("\nAlignment status: {:?}", status);

    if matches!(status, AlignmentStatus::Completed) {
        let cigar = aligner.cigar();
        println!("CIGAR length: {}", cigar.len());

        // Count operation types
        let mut counts = std::collections::HashMap::new();
        for &op in cigar {
            *counts.entry(op).or_insert(0) += 1;
        }

        println!("\nOperation counts:");
        for (op, count) in &counts {
            println!("  {} ({}): {}", *op as char, op, count);
        }

        // Find where the issue occurs
        let mut q_pos = 0;
        let mut r_pos = 0;
        let mut _last_good_pos = 0;

        for (i, &op) in cigar.iter().enumerate() {
            match op {
                b'M' | b'=' | b'X' => {
                    if q_pos < seq1.len() && r_pos < seq2.len() {
                        _last_good_pos = i;
                    }
                    q_pos += 1;
                    r_pos += 1;
                }
                b'I' => q_pos += 1,
                b'D' => r_pos += 1,
                _ => {}
            }

            if q_pos > seq1.len() || r_pos > seq2.len() {
                println!("\nIssue at CIGAR position {}: op={} ({})", i, op as char, op);
                println!("Query position: {} / {}", q_pos, seq1.len());
                println!("Reference position: {} / {}", r_pos, seq2.len());
                
                // Show context
                if i >= 5 {
                    print!("Previous 5 ops: ");
                    for j in (i-5)..i {
                        print!("{}", cigar[j] as char);
                    }
                    println!();
                }
                
                break;
            }
        }

        println!("\nFinal positions:");
        println!("Query: {} / {}", q_pos, seq1.len());
        println!("Reference: {} / {}", r_pos, seq2.len());
    }
}