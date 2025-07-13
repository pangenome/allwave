//! General WFA2 debugging tool
//! 
//! Purpose: A flexible debugging tool that can analyze alignments from any FASTA
//! file. Tests alignments in both orientations and provides detailed CIGAR analysis.
//! This was one of the main tools used to understand WFA2's behavior.
//! 
//! Usage: cargo run --bin debug_wfa <fasta_file>

use bio::io::fasta;
use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentStatus, HeuristicStrategy, MemoryMode,
};
use std::env;
use std::fs::File;
use std::io::BufReader;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <fasta_file>", args[0]);
        std::process::exit(1);
    }

    let file = File::open(&args[1]).expect("Cannot open file");
    let reader = BufReader::new(file);
    let fasta_reader = fasta::Reader::new(reader);
    let mut sequences = Vec::new();

    for result in fasta_reader.records() {
        let record = result.expect("Error reading record");
        sequences.push((record.id().to_string(), record.seq().to_vec()));
    }

    if sequences.len() < 2 {
        eprintln!("Need at least 2 sequences");
        return;
    }

    // Test with the problematic sequences
    for i in 0..sequences.len().min(2) {
        for j in i + 1..sequences.len().min(2) {
            let (name1, seq1) = &sequences[i];
            let (name2, seq2) = &sequences[j];

            println!("\n=== Aligning {} vs {} ===", name1, name2);
            println!("Seq1 length: {}", seq1.len());
            println!("Seq2 length: {}", seq2.len());

            // Create aligner with two-piece affine and ultralow memory
            let mut aligner = AffineWavefronts::with_penalties_affine2p_and_memory_mode(
                0, 5, 8, 2, 24, 1, MemoryMode::Ultralow
            );

            // Set critical configuration
            aligner.set_alignment_scope(AlignmentScope::Alignment);
            aligner.set_heuristic(&HeuristicStrategy::None);

            // Try alignment
            let status = aligner.align(seq1, seq2);
            println!("Alignment status: {:?}", status);

            if matches!(status, AlignmentStatus::Completed) {
                let score = aligner.score();
                let cigar_bytes = aligner.cigar();
                println!("Score: {}", score);
                println!("CIGAR length: {}", cigar_bytes.len());

                // Count CIGAR operations
                let mut query_pos = 0;
                let mut ref_pos = 0;
                for &op in cigar_bytes {
                    match op {
                        b'M' | b'=' | b'X' => {
                            query_pos += 1;
                            ref_pos += 1;
                        }
                        b'I' => query_pos += 1,
                        b'D' => ref_pos += 1,
                        _ => println!("Unknown op: {} (0x{:02x})", op as char, op),
                    }
                }

                println!("CIGAR covers query: {} / {}", query_pos, seq1.len());
                println!("CIGAR covers reference: {} / {}", ref_pos, seq2.len());

                // Show first and last 10 operations
                if cigar_bytes.len() > 20 {
                    print!("First 10 ops: ");
                    for &op in &cigar_bytes[..10] {
                        print!("{}", op as char);
                    }
                    println!();

                    print!("Last 10 ops: ");
                    for &op in &cigar_bytes[cigar_bytes.len() - 10..] {
                        print!("{}", op as char);
                    }
                    println!();
                }

                // Try swapping sequences
                println!("\n--- Trying swapped order ---");
                let status2 = aligner.align(seq2, seq1);
                println!("Swapped status: {:?}", status2);

                if matches!(status2, AlignmentStatus::Completed) {
                    let cigar_bytes2 = aligner.cigar();
                    let mut query_pos2 = 0;
                    let mut ref_pos2 = 0;
                    for &op in cigar_bytes2 {
                        match op {
                            b'M' | b'=' | b'X' => {
                                query_pos2 += 1;
                                ref_pos2 += 1;
                            }
                            b'I' => query_pos2 += 1,
                            b'D' => ref_pos2 += 1,
                            _ => {}
                        }
                    }
                    println!("Swapped CIGAR covers query: {} / {}", query_pos2, seq2.len());
                    println!("Swapped CIGAR covers reference: {} / {}", ref_pos2, seq1.len());
                }
            }
        }
    }
}