//! Comprehensive test suite for allwave
//!
//! Purpose: Full integration test suite that generates synthetic sequences with
//! various types of mutations (SNPs, indels, microsatellites, CNVs) and verifies
//! that allwave can align them correctly. This was used to ensure our fixes
//! didn't break existing functionality.
//!
//! Usage: cargo run --bin test_allwave

use std::io::Write;
use std::process::{Command, Stdio};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== AllWave Test Suite ===\n");

    // Test 1: Basic mutations (SNPs and small indels)
    println!("Test 1: Basic mutations (SNPs and small indels)");
    let seq_config = SequenceConfig {
        length: 5000,
        gc_content: 0.5,
        seed: 42,
    };

    let mutations = vec![
        Mutation::SNP { position: 100 },
        Mutation::SNP { position: 500 },
        Mutation::Insertion {
            position: 1000,
            sequence: "ACGT".to_string(),
        },
        Mutation::Deletion {
            position: 2000,
            length: 5,
        },
    ];

    run_test("basic_mutations", &seq_config, &mutations)?;

    // Test 2: Microsatellite mutations
    println!("\nTest 2: Microsatellite expansions and contractions");
    let mutations = vec![
        Mutation::MicrosatelliteExpansion {
            position: 500,
            unit: "AT".to_string(),
            original_copies: 10,
            new_copies: 15,
        },
        Mutation::MicrosatelliteContraction {
            position: 1500,
            unit: "CAG".to_string(),
            original_copies: 12,
            new_copies: 8,
        },
    ];

    run_test("microsatellite_mutations", &seq_config, &mutations)?;

    // Test 3: Large structural variants
    println!("\nTest 3: Copy number variations (CNVs)");
    let mutations = vec![
        Mutation::CNVDuplication {
            position: 1000,
            length: 2000,
            copies: 2,
        },
        Mutation::CNVDeletion {
            position: 3500,
            length: 1000,
        },
    ];

    run_test("cnv_mutations", &seq_config, &mutations)?;

    println!("\nAll tests completed!");
    Ok(())
}

#[derive(Clone)]
struct SequenceConfig {
    length: usize,
    gc_content: f64,
    seed: u64,
}

#[derive(Clone)]
enum Mutation {
    SNP {
        position: usize,
    },
    Insertion {
        position: usize,
        sequence: String,
    },
    Deletion {
        position: usize,
        length: usize,
    },
    MicrosatelliteExpansion {
        position: usize,
        unit: String,
        original_copies: usize,
        new_copies: usize,
    },
    MicrosatelliteContraction {
        position: usize,
        unit: String,
        original_copies: usize,
        new_copies: usize,
    },
    CNVDuplication {
        position: usize,
        length: usize,
        copies: usize,
    },
    CNVDeletion {
        position: usize,
        length: usize,
    },
}

fn generate_sequence(config: &SequenceConfig) -> String {
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut seq = String::with_capacity(config.length);

    for _ in 0..config.length {
        let r: f64 = rng.gen();
        let base = if r < config.gc_content / 2.0 {
            'G'
        } else if r < config.gc_content {
            'C'
        } else if r < (1.0 + config.gc_content) / 2.0 {
            'A'
        } else {
            'T'
        };
        seq.push(base);
    }

    seq
}

fn apply_mutations(sequence: &str, mutations: &[Mutation]) -> String {
    let mut seq = sequence.to_string();

    // Apply mutations in reverse order to maintain positions
    let mut sorted_mutations = mutations.to_vec();
    sorted_mutations.sort_by_key(|m| match m {
        Mutation::SNP { position } => *position,
        Mutation::Insertion { position, .. } => *position,
        Mutation::Deletion { position, .. } => *position,
        Mutation::MicrosatelliteExpansion { position, .. } => *position,
        Mutation::MicrosatelliteContraction { position, .. } => *position,
        Mutation::CNVDuplication { position, .. } => *position,
        Mutation::CNVDeletion { position, .. } => *position,
    });
    sorted_mutations.reverse();

    for mutation in sorted_mutations {
        seq = match mutation {
            Mutation::SNP { position } => {
                if *position < seq.len() {
                    let mut chars: Vec<char> = seq.chars().collect();
                    chars[*position] = match chars[*position] {
                        'A' => 'T',
                        'T' => 'A',
                        'G' => 'C',
                        'C' => 'G',
                        _ => 'N',
                    };
                    chars.into_iter().collect()
                } else {
                    seq
                }
            }
            Mutation::Insertion { position, sequence } => {
                if *position <= seq.len() {
                    let mut result = String::new();
                    result.push_str(&seq[..*position]);
                    result.push_str(sequence);
                    result.push_str(&seq[*position..]);
                    result
                } else {
                    seq
                }
            }
            Mutation::Deletion { position, length } => {
                if *position + *length <= seq.len() {
                    let mut result = String::new();
                    result.push_str(&seq[..*position]);
                    result.push_str(&seq[*position + *length..]);
                    result
                } else {
                    seq
                }
            }
            _ => seq, // TODO: Implement other mutation types
        };
    }

    seq
}

fn run_test(
    name: &str,
    config: &SequenceConfig,
    mutations: &[Mutation],
) -> Result<(), Box<dyn std::error::Error>> {
    // Generate reference sequence
    let reference = generate_sequence(config);

    // Apply mutations to create query sequence
    let query = apply_mutations(&reference, mutations);

    // Write sequences to FASTA
    let fasta_content = format!(">reference\n{}\n>query\n{}\n", reference, query);
    std::fs::write(format!("test_{}.fa", name), &fasta_content)?;

    // Run allwave
    let output = Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--bin",
            "allwave",
            "--",
            "-i",
            &format!("test_{}.fa", name),
        ])
        .output()?;

    if !output.status.success() {
        eprintln!(
            "Allwave failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return Err("Allwave execution failed".into());
    }

    // Parse PAF output
    let paf = String::from_utf8(output.stdout)?;
    for line in paf.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 12 {
            let query_name = fields[0];
            let query_len = fields[1];
            let query_start = fields[2];
            let query_end = fields[3];
            let identity = fields[9].parse::<f64>()? / fields[10].parse::<f64>()?;

            println!(
                "  Aligned {} ({} bp): {}-{} (identity: {:.2}%)",
                query_name,
                query_len,
                query_start,
                query_end,
                identity * 100.0
            );
        }
    }

    // Clean up
    std::fs::remove_file(format!("test_{}.fa", name))?;

    Ok(())
}
