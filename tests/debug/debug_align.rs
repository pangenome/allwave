//! Debug sequence alignment differences
//! 
//! Purpose: Simple debugging tool to check character-by-character differences
//! between two sequences. Used during development to understand why certain
//! alignments were failing or producing unexpected results.
//! 
//! Usage: cargo run --bin debug_align

fn main() {
    // Read sequences from files (created by debug scripts)
    let seq1 = std::fs::read_to_string("seq1.txt").unwrap_or_default().trim().to_string();
    let seq2 = std::fs::read_to_string("seq2.txt").unwrap_or_default().trim().to_string();
    
    println!("Seq1 len: {}", seq1.len());
    println!("Seq2 len: {}", seq2.len());
    
    // Check difference
    for (i, (a, b)) in seq1.chars().zip(seq2.chars()).enumerate() {
        if a != b {
            println!("Diff at {}: {} vs {}", i, a, b);
            if i > 10 {
                println!("Context: ...{}[{}]{}", 
                    &seq1[i.saturating_sub(10)..i],
                    a,
                    &seq1[i+1..std::cmp::min(i+11, seq1.len())]);
                println!("         ...{}[{}]{}", 
                    &seq2[i.saturating_sub(10)..i],
                    b,
                    &seq2[i+1..std::cmp::min(i+11, seq2.len())]);
            }
            break;
        }
    }
}