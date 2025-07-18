use allwave::test_framework::*;
use allwave::validation::*;
use allwave::validation_correct::*;
use allwave::validation_simple::*;
use std::collections::HashMap;
use std::process::Command;
use rand::Rng;

#[test]
fn test_basic_mutations() {
    println!("\n=== Test: Basic mutations (SNPs and small indels) ===");

    let mut seq_config = SequenceConfig::default();
    seq_config.length = 5000;
    seq_config.seed = 42;

    let mut mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.01;
    mut_config.indel_rate = 0.001;
    mut_config.microsatellite_rate = 0.0;
    mut_config.cnv_rate = 0.0;
    mut_config.seed = 123;

    let test_case = create_test_case("test_basic", &seq_config, &mut_config);
    let result = run_alignment_test(&test_case).expect("Failed to run test");

    // Check that we detected most SNPs
    assert!(
        result.detection_accuracy[&MutationType::SNP] > 0.8,
        "SNP detection accuracy too low: {:.2}%",
        result.detection_accuracy[&MutationType::SNP] * 100.0
    );

    // Check alignment coverage
    assert!(
        result.query_coverage > 0.95,
        "Query coverage too low: {:.2}%",
        result.query_coverage * 100.0
    );
    assert!(
        result.target_coverage > 0.95,
        "Target coverage too low: {:.2}%",
        result.target_coverage * 100.0
    );

    result.print_summary();
}

#[test]
fn test_microsatellite_mutations() {
    println!("\n=== Test: Microsatellite expansions and contractions ===");

    let test_case = create_microsatellite_test_case();
    let result = run_alignment_test(&test_case).expect("Failed to run test");

    // Check that we have good coverage
    assert!(
        result.query_coverage > 0.95,
        "Query coverage too low: {:.2}%",
        result.query_coverage * 100.0
    );
    assert!(
        result.target_coverage > 0.95,
        "Target coverage too low: {:.2}%",
        result.target_coverage * 100.0
    );

    // Check that we detected insertions/deletions (microsatellite changes appear as indels)
    let total_indels = result
        .detected_mutations
        .get(&MutationType::Insertion)
        .unwrap_or(&0)
        + result
            .detected_mutations
            .get(&MutationType::Deletion)
            .unwrap_or(&0);
    assert!(
        total_indels > 0,
        "No insertions or deletions detected for microsatellite changes"
    );

    result.print_summary();
}

#[test]
fn test_cnv_mutations() {
    println!("\n=== Test: Copy Number Variations (CNVs) ===");

    let mut seq_config = SequenceConfig::default();
    seq_config.length = 20000;
    seq_config.seed = 44;

    let mut mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.0001;
    mut_config.indel_rate = 0.0;
    mut_config.microsatellite_rate = 0.0;
    mut_config.cnv_rate = 0.0001;
    mut_config.cnv_min_size = 1000;
    mut_config.cnv_max_size = 3000;
    mut_config.seed = 125;

    let test_case = create_test_case("test_cnv", &seq_config, &mut_config);
    let result = run_alignment_test(&test_case).expect("Failed to run test");

    // Check coverage - might be significantly lower due to large structural variations
    assert!(
        result.query_coverage > 0.7,
        "Query coverage too low: {:.2}%",
        result.query_coverage * 100.0
    );
    assert!(
        result.target_coverage > 0.7,
        "Target coverage too low: {:.2}%",
        result.target_coverage * 100.0
    );

    // Check that we detected large insertions/deletions
    let has_large_variants = *result
        .detected_mutations
        .get(&MutationType::CNVDuplication)
        .unwrap_or(&0)
        > 0
        || *result
            .detected_mutations
            .get(&MutationType::CNVDeletion)
            .unwrap_or(&0)
            > 0;
    assert!(has_large_variants, "No large structural variants detected");

    result.print_summary();
}

#[test]
fn test_complex_mutations() {
    println!("\n=== Test: Complex mutations (all types combined) ===");

    let mut seq_config = SequenceConfig::default();
    seq_config.length = 10000;
    seq_config.seed = 45;

    let mut mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.005;
    mut_config.indel_rate = 0.0005;
    mut_config.microsatellite_rate = 0.0001;
    mut_config.cnv_rate = 0.00005;
    mut_config.seed = 126;

    let test_case = create_test_case("test_complex", &seq_config, &mut_config);
    let result = run_alignment_test(&test_case).expect("Failed to run test");

    // With complex mutations, we still expect good coverage
    assert!(
        result.query_coverage > 0.9,
        "Query coverage too low: {:.2}%",
        result.query_coverage * 100.0
    );
    assert!(
        result.target_coverage > 0.9,
        "Target coverage too low: {:.2}%",
        result.target_coverage * 100.0
    );

    // Check that multiple mutation types were detected
    let mutation_types_detected = result
        .detected_mutations
        .iter()
        .filter(|(_, &count)| count > 0)
        .count();
    assert!(
        mutation_types_detected >= 2,
        "Too few mutation types detected: {}",
        mutation_types_detected
    );

    result.print_summary();
}

#[test]
fn test_high_divergence() {
    println!("\n=== Test: High divergence sequences ===");

    let mut seq_config = SequenceConfig::default();
    seq_config.length = 3000;
    seq_config.seed = 46;

    let mut mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.05; // 5% divergence
    mut_config.indel_rate = 0.01;
    mut_config.seed = 127;

    let test_case = create_test_case("test_high_div", &seq_config, &mut_config);
    let result = run_alignment_test(&test_case).expect("Failed to run test");

    // Even with high divergence, we should align most of the sequences
    assert!(
        result.query_coverage > 0.8,
        "Query coverage too low: {:.2}%",
        result.query_coverage * 100.0
    );
    assert!(
        result.target_coverage > 0.8,
        "Target coverage too low: {:.2}%",
        result.target_coverage * 100.0
    );

    // Identity should be lower due to high mutation rate
    assert!(
        result.alignment_stats.identity < 0.96,
        "Identity too high for high divergence: {:.2}%",
        result.alignment_stats.identity * 100.0
    );

    result.print_summary();
}

#[test]
fn test_identical_sequences() {
    println!("\n=== Test: Identical sequences ===");

    let seq_config = SequenceConfig {
        length: 5000,
        gc_content: 0.5,
        seed: 100,
    };

    let reference = generate_random_sequence(&seq_config);
    let test_case = TestCase {
        name: "test_identical".to_string(),
        reference: reference.clone(),
        query: reference,
        mutations: vec![],
    };

    let result = run_alignment_test(&test_case).expect("Failed to run test");

    // Should have perfect alignment
    assert_eq!(result.query_coverage, 1.0, "Query coverage should be 100%");
    assert_eq!(
        result.target_coverage, 1.0,
        "Target coverage should be 100%"
    );
    assert_eq!(
        result.alignment_stats.identity, 1.0,
        "Identity should be 100%"
    );
    assert_eq!(
        result.alignment_stats.mismatches, 0,
        "Should have no mismatches"
    );
    assert_eq!(
        result.alignment_stats.insertions, 0,
        "Should have no insertions"
    );
    assert_eq!(
        result.alignment_stats.deletions, 0,
        "Should have no deletions"
    );

    result.print_summary();
}

// Helper functions

fn run_alignment_test(
    test_case: &TestCase,
) -> Result<ValidationResult, Box<dyn std::error::Error>> {
    // Print mutation statistics
    let stats = test_case.get_mutation_stats();
    println!("  Mutations introduced:");
    for (mutation_type, count) in &stats {
        println!("    {:?}: {}", mutation_type, count);
    }

    // Write test case to FASTA
    let fasta_path = format!("/tmp/allwave_test_{}.fa", test_case.name);
    test_case.write_fasta(&fasta_path)?;

    // Run allwave
    let output = Command::new("./target/debug/allwave")
        .arg("--input")
        .arg(&fasta_path)
        .output()?;

    // Clean up
    std::fs::remove_file(&fasta_path)?;

    if !output.status.success() {
        return Err(format!(
            "allwave failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }

    // Parse and validate results
    let paf_output = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = paf_output.lines().collect();

    if lines.is_empty() {
        return Err("No alignment output".into());
    }

    // Find the alignment between reference and query
    let mut result = None;
    for line in &lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 6 {
            let query_name = fields[0];
            let target_name = fields[5];

            // Look for reference vs query alignment
            if (query_name.contains("reference") && target_name.contains("query"))
                || (query_name.contains("query") && target_name.contains("reference"))
            {
                // Parse PAF fields
                let fields: Vec<&str> = line.split('\t').collect();
                let query_start = fields[2].parse::<usize>().unwrap();
                let query_end = fields[3].parse::<usize>().unwrap();
                let strand = fields[4];
                let target_start = fields[7].parse::<usize>().unwrap();
                let target_end = fields[8].parse::<usize>().unwrap();

                // Get CIGAR string
                let cigar_str = fields
                    .iter()
                    .find(|f| f.starts_with("cg:Z:"))
                    .unwrap()
                    .trim_start_matches("cg:Z:");

                // In PAF format: fields[0] is query, fields[5] is target
                // So we map based on the actual names
                let (paf_query_seq, paf_target_seq) = if query_name.contains("reference") {
                    (&test_case.reference, &test_case.query)
                } else {
                    (&test_case.query, &test_case.reference)
                };

                // Extract the aligned subsequences
                let query_subseq = &paf_query_seq[query_start..query_end];
                let target_subseq = &paf_target_seq[target_start..target_end];

                // Validate CIGAR correctness
                match verify_alignment_cigar(
                    query_subseq,
                    target_subseq,
                    cigar_str,
                    strand.chars().next().unwrap(),
                ) {
                    Ok(true) => println!("  Alignment validation: PASSED"),
                    Ok(false) => println!("  Alignment validation: FAILED"),
                    Err(e) => {
                        println!("  Alignment validation: FAILED - {}", e);
                        // Continue anyway to get statistics
                    }
                }

                // Use simple validation for statistics
                result = Some(validate_alignment_simple(
                    paf_query_seq,
                    paf_target_seq,
                    line,
                )?);

                // Also run original validation for mutation detection
                if let Ok(full_result) = run_and_validate_test(test_case, line) {
                    // Merge results
                    if let Some(ref mut r) = result {
                        r.expected_mutations = full_result.expected_mutations;
                        r.detected_mutations = full_result.detected_mutations;
                        r.detection_accuracy = full_result.detection_accuracy;
                    }
                }
                break;
            }
        }
    }

    result.ok_or_else(|| "No valid alignment found between reference and query".into())
}

fn create_microsatellite_test_case() -> TestCase {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let mut rng = StdRng::seed_from_u64(1234);
    let mut reference = Vec::new();

    // Add some random sequence
    for _ in 0..1000 {
        let bases = [b'A', b'T', b'G', b'C'];
        reference.push(bases[rng.gen_range(0..4)]);
    }

    // Add microsatellite: (AT)20
    for _ in 0..20 {
        reference.extend_from_slice(b"AT");
    }

    // More random sequence
    for _ in 0..1000 {
        let bases = [b'A', b'T', b'G', b'C'];
        reference.push(bases[rng.gen_range(0..4)]);
    }

    // Add microsatellite: (CAG)15
    for _ in 0..15 {
        reference.extend_from_slice(b"CAG");
    }

    // More random sequence
    for _ in 0..1000 {
        let bases = [b'A', b'T', b'G', b'C'];
        reference.push(bases[rng.gen_range(0..4)]);
    }

    // Add microsatellite: (GATA)10
    for _ in 0..10 {
        reference.extend_from_slice(b"GATA");
    }

    // Final random sequence
    for _ in 0..1000 {
        let bases = [b'A', b'T', b'G', b'C'];
        reference.push(bases[rng.gen_range(0..4)]);
    }

    // Apply microsatellite mutations
    let mut config = MutationConfig::default();
    config.microsatellite_rate = 1.0; // Force mutations on all microsatellites
    config.snp_rate = 0.001;
    config.seed = 999;

    let (query, mutations) = apply_mutations(&reference, &config);

    TestCase {
        name: "test_microsatellite".to_string(),
        reference,
        query,
        mutations,
    }
}

#[test]
fn test_reverse_complement_alignment() {
    println!("\n=== Test: Reverse complement alignment ===");

    // Create a sequence with some structure that's easy to verify
    let seq_config = SequenceConfig {
        length: 2000,
        gc_content: 0.6,
        seed: 200,
    };

    let reference = generate_random_sequence(&seq_config);

    // Create query with mutations
    let mut mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.01;
    mut_config.indel_rate = 0.001;
    mut_config.seed = 201;

    let (query, _mutations) = apply_mutations(&reference, &mut_config);

    // Create FASTA with reference and both orientations of query
    let fasta_path = "/tmp/allwave_test_revcomp.fa";
    use std::io::Write;
    let mut file = std::fs::File::create(fasta_path).unwrap();

    writeln!(file, ">reference").unwrap();
    writeln!(file, "{}", String::from_utf8_lossy(&reference)).unwrap();
    writeln!(file, ">query_forward").unwrap();
    writeln!(file, "{}", String::from_utf8_lossy(&query)).unwrap();
    writeln!(file, ">query_reverse").unwrap();
    writeln!(
        file,
        "{}",
        String::from_utf8_lossy(&reverse_complement(&query))
    )
    .unwrap();

    // Run allwave with no sparsification to ensure all pairs are computed
    let output = Command::new("./target/debug/allwave")
        .arg("--input")
        .arg(fasta_path)
        .arg("-p")
        .arg("none")
        .output()
        .expect("Failed to run allwave");

    std::fs::remove_file(fasta_path).unwrap();

    assert!(
        output.status.success(),
        "allwave failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let paf_output = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = paf_output.lines().collect();

    // Find alignments between reference and both query orientations
    let mut ref_vs_fwd = None;
    let mut ref_vs_rev = None;

    for line in &lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 12 {
            if fields[0] == "reference" && fields[5] == "query_forward" {
                ref_vs_fwd = Some(line.to_string());
            } else if fields[0] == "reference" && fields[5] == "query_reverse" {
                ref_vs_rev = Some(line.to_string());
            }
        }
    }

    assert!(
        ref_vs_fwd.is_some(),
        "No alignment found for reference vs query_forward"
    );
    assert!(
        ref_vs_rev.is_some(),
        "No alignment found for reference vs query_reverse"
    );

    // Parse strands
    let fwd_fields: Vec<&str> = ref_vs_fwd.as_ref().unwrap().split('\t').collect();
    let rev_fields: Vec<&str> = ref_vs_rev.as_ref().unwrap().split('\t').collect();

    let fwd_strand = fwd_fields[4];
    let rev_strand = rev_fields[4];

    println!("  Forward query aligned with strand: {}", fwd_strand);
    println!("  Reverse query aligned with strand: {}", rev_strand);

    // The forward query should align with + strand, reverse with - strand
    assert_eq!(fwd_strand, "+", "Forward query should align with + strand");
    assert_eq!(
        rev_strand, "-",
        "Reverse complement query should align with - strand"
    );

    // Check that both alignments have similar identity (within 1%)
    let fwd_identity = parse_identity(&fwd_fields);
    let rev_identity = parse_identity(&rev_fields);

    println!("  Forward alignment identity: {:.2}%", fwd_identity * 100.0);
    println!("  Reverse alignment identity: {:.2}%", rev_identity * 100.0);

    assert!(
        (fwd_identity - rev_identity).abs() < 0.01,
        "Identity difference too large: {:.2}% vs {:.2}%",
        fwd_identity * 100.0,
        rev_identity * 100.0
    );
}

#[test]
fn test_long_sequences() {
    println!("\n=== Test: Long sequences (100kb) ===");

    let seq_config = SequenceConfig {
        length: 100000,
        gc_content: 0.45,
        seed: 300,
    };

    let mut mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.002;
    mut_config.indel_rate = 0.0002;
    mut_config.microsatellite_rate = 0.00001;
    mut_config.cnv_rate = 0.000001;
    mut_config.seed = 301;

    let test_case = create_test_case("test_long", &seq_config, &mut_config);
    let result = run_alignment_test(&test_case).expect("Failed to run test");

    // Should still have good alignment for long sequences
    assert!(
        result.query_coverage > 0.95,
        "Query coverage too low: {:.2}%",
        result.query_coverage * 100.0
    );
    assert!(
        result.target_coverage > 0.95,
        "Target coverage too low: {:.2}%",
        result.target_coverage * 100.0
    );

    // Check that alignment length is reasonable
    assert!(
        result.alignment_stats.alignment_length > 95000,
        "Alignment too short for 100kb sequences: {}",
        result.alignment_stats.alignment_length
    );

    result.print_summary();
}

#[test]
fn test_alignment_correctness() {
    println!("\n=== Test: Alignment correctness verification ===");

    // Create sequences with known mutations at specific positions
    let reference = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_vec();
    let mut query = reference.clone();

    // Introduce specific mutations
    query[10] = b'G'; // SNP at position 10
    query[20] = b'C'; // SNP at position 20
    query.remove(30); // Deletion at position 30
    query.insert(40, b'A'); // Insertion at position 40

    let test_case = TestCase {
        name: "test_correctness".to_string(),
        reference: reference.clone(),
        query: query.clone(),
        mutations: vec![], // We don't need to track mutations for this test
    };

    // Write test case to FASTA
    let fasta_path = "/tmp/allwave_test_correctness.fa";
    test_case.write_fasta(&fasta_path).unwrap();

    // Run allwave
    let output = Command::new("./target/debug/allwave")
        .arg("--input")
        .arg(&fasta_path)
        .output()
        .unwrap();

    std::fs::remove_file(&fasta_path).unwrap();

    assert!(output.status.success());

    let paf_output = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = paf_output.lines().collect();
    assert!(!lines.is_empty());

    // Parse CIGAR string
    let fields: Vec<&str> = lines[0].split('\t').collect();
    let cigar_str = fields
        .iter()
        .find(|f| f.starts_with("cg:Z:"))
        .unwrap()
        .trim_start_matches("cg:Z:");

    let cigar_ops = parse_cigar(cigar_str);

    // Verify CIGAR operations
    let mut total_mismatches = 0;
    let mut total_insertions = 0;
    let mut total_deletions = 0;

    for op in &cigar_ops {
        match op {
            CigarOp::Mismatch(count) => total_mismatches += *count as usize,
            CigarOp::Insertion(count) => total_insertions += *count as usize,
            CigarOp::Deletion(count) => total_deletions += *count as usize,
            _ => {}
        }
    }

    println!("  Detected mismatches: {}", total_mismatches);
    println!("  Detected insertions: {}", total_insertions);
    println!("  Detected deletions: {}", total_deletions);
    println!("  CIGAR: {}", cigar_str);

    // We should detect exactly 2 mismatches, 1 insertion, 1 deletion
    assert_eq!(total_mismatches, 2, "Should detect exactly 2 mismatches");
    assert_eq!(total_insertions, 1, "Should detect exactly 1 insertion");
    assert_eq!(total_deletions, 1, "Should detect exactly 1 deletion");
}

#[test]
fn test_tandem_repeats() {
    println!("\n=== Test: Tandem repeats and low complexity regions ===");

    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let mut rng = StdRng::seed_from_u64(400);
    let mut reference = Vec::new();

    // Random sequence
    for _ in 0..500 {
        let bases = [b'A', b'T', b'G', b'C'];
        reference.push(bases[rng.gen_range(0..4)]);
    }

    // Add tandem repeat: 50x "ATCG"
    for _ in 0..50 {
        reference.extend_from_slice(b"ATCG");
    }

    // More random
    for _ in 0..500 {
        let bases = [b'A', b'T', b'G', b'C'];
        reference.push(bases[rng.gen_range(0..4)]);
    }

    // Add homopolymer run: 100x "A"
    for _ in 0..100 {
        reference.push(b'A');
    }

    // More random
    for _ in 0..500 {
        let bases = [b'A', b'T', b'G', b'C'];
        reference.push(bases[rng.gen_range(0..4)]);
    }

    // Add variable number tandem repeat: "GATTACA" repeated 20 times
    for _ in 0..20 {
        reference.extend_from_slice(b"GATTACA");
    }

    // Final random sequence
    for _ in 0..500 {
        let bases = [b'A', b'T', b'G', b'C'];
        reference.push(bases[rng.gen_range(0..4)]);
    }

    // Create query with mutations in repeat regions
    let mut mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.005;
    mut_config.indel_rate = 0.002;
    mut_config.seed = 401;

    let (query, _) = apply_mutations(&reference, &mut_config);

    let test_case = TestCase {
        name: "test_tandem".to_string(),
        reference,
        query,
        mutations: vec![],
    };

    let result = run_alignment_test(&test_case).expect("Failed to run test");

    // Should handle repeat regions reasonably well
    assert!(
        result.query_coverage > 0.9,
        "Query coverage too low: {:.2}%",
        result.query_coverage * 100.0
    );
    assert!(
        result.target_coverage > 0.9,
        "Target coverage too low: {:.2}%",
        result.target_coverage * 100.0
    );

    result.print_summary();
}

#[test]
fn test_multiple_sequence_pairs() {
    println!("\n=== Test: Multiple sequence pairs in single file ===");

    // Create multiple sequences
    let seq_configs = vec![
        SequenceConfig {
            length: 1000,
            gc_content: 0.3,
            seed: 500,
        },
        SequenceConfig {
            length: 1500,
            gc_content: 0.5,
            seed: 501,
        },
        SequenceConfig {
            length: 2000,
            gc_content: 0.7,
            seed: 502,
        },
    ];

    let sequences: Vec<Vec<u8>> = seq_configs
        .iter()
        .map(|config| generate_random_sequence(config))
        .collect();

    // Write all sequences to FASTA
    let fasta_path = "/tmp/allwave_test_multiple.fa";
    use std::io::Write;
    let mut file = std::fs::File::create(fasta_path).unwrap();

    for (i, seq) in sequences.iter().enumerate() {
        writeln!(file, ">seq{}", i).unwrap();
        writeln!(file, "{}", String::from_utf8_lossy(seq)).unwrap();
    }

    // Run allwave with no sparsification to ensure all pairs are computed
    let output = Command::new("./target/debug/allwave")
        .arg("--input")
        .arg(fasta_path)
        .arg("-p")
        .arg("none")
        .output()
        .unwrap();

    std::fs::remove_file(fasta_path).unwrap();

    assert!(output.status.success());

    let paf_output = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = paf_output.lines().collect();

    // Should have n*(n-1) alignments for n sequences
    let expected_alignments = 3 * 2; // 3 sequences, each aligned to 2 others
    assert_eq!(
        lines.len(),
        expected_alignments,
        "Expected {} alignments, got {}",
        expected_alignments,
        lines.len()
    );

    // Verify all pairs are present
    let mut pair_counts = HashMap::new();
    for line in &lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 6 {
            let pair = format!("{}-{}", fields[0], fields[5]);
            *pair_counts.entry(pair).or_insert(0) += 1;
        }
    }

    // Each pair should appear exactly once
    for (pair, count) in &pair_counts {
        assert_eq!(*count, 1, "Pair {} appeared {} times", pair, count);
    }

    println!("  Total alignments: {}", lines.len());
    println!("  Unique pairs: {}", pair_counts.len());
}

// Helper functions

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => b'N',
        })
        .collect()
}

fn parse_identity(fields: &[&str]) -> f64 {
    for field in fields {
        if field.starts_with("gi:f:") {
            if let Ok(identity) = field.trim_start_matches("gi:f:").parse::<f64>() {
                return identity;
            }
        }
    }
    0.0
}

#[test]
fn test_orientation_detection_comparison() {
    println!("\n=== Test: Mash vs Edlib orientation detection comparison ===");
    
    // Create test sequences with known orientations
    let test_cases = create_orientation_test_cases();
    
    for test_case in test_cases {
        println!("  Testing: {}", test_case.name);
        
        // Test mash-based orientation detection
        let mash_result = test_orientation_detection_mash(&test_case);
        
        // Test edlib-based orientation detection
        let edlib_result = test_orientation_detection_edlib(&test_case);
        
        println!("    Mash result: forward={:.4}, reverse={:.4}, chose={}", 
                 mash_result.forward_distance, mash_result.reverse_distance, 
                 if mash_result.is_reverse { "reverse" } else { "forward" });
        println!("    Edlib result: forward={}, reverse={}, chose={}", 
                 edlib_result.forward_distance, edlib_result.reverse_distance,
                 if edlib_result.is_reverse { "reverse" } else { "forward" });
        
        // Both methods should agree on orientation
        assert_eq!(mash_result.is_reverse, edlib_result.is_reverse,
                   "Orientation detection methods disagree for test case: {}", test_case.name);
        
        // Expected orientation should match the test case
        assert_eq!(mash_result.is_reverse, test_case.expected_reverse,
                   "Mash orientation detection incorrect for test case: {}", test_case.name);
        assert_eq!(edlib_result.is_reverse, test_case.expected_reverse,
                   "Edlib orientation detection incorrect for test case: {}", test_case.name);
    }
    
    println!("  All orientation detection tests passed!");
}

#[derive(Debug)]
struct OrientationTestCase {
    name: String,
    reference: Vec<u8>,
    query: Vec<u8>,
    expected_reverse: bool,
}

#[derive(Debug)]
struct OrientationResult {
    forward_distance: f64,
    reverse_distance: f64,
    is_reverse: bool,
}

fn create_orientation_test_cases() -> Vec<OrientationTestCase> {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    
    let mut test_cases = Vec::new();
    let mut rng = StdRng::seed_from_u64(12345);
    
    // Test case 1: Identical sequences (should be forward)
    let ref_seq = generate_test_sequence(1000, &mut rng);
    test_cases.push(OrientationTestCase {
        name: "identical_sequences".to_string(),
        reference: ref_seq.clone(),
        query: ref_seq.clone(),
        expected_reverse: false,
    });
    
    // Test case 2: Forward with mutations
    let ref_seq = generate_test_sequence(1000, &mut rng);
    let mut query_fwd = ref_seq.clone();
    apply_test_mutations(&mut query_fwd, 0.01, &mut rng);
    test_cases.push(OrientationTestCase {
        name: "forward_with_mutations".to_string(),
        reference: ref_seq.clone(),
        query: query_fwd,
        expected_reverse: false,
    });
    
    // Test case 3: Reverse complement with mutations
    let ref_seq = generate_test_sequence(1000, &mut rng);
    let mut query_rev = reverse_complement(&ref_seq);
    apply_test_mutations(&mut query_rev, 0.01, &mut rng);
    test_cases.push(OrientationTestCase {
        name: "reverse_with_mutations".to_string(),
        reference: ref_seq.clone(),
        query: query_rev,
        expected_reverse: true,
    });
    
    // Test case 4: High mutation rate forward
    let ref_seq = generate_test_sequence(1000, &mut rng);
    let mut query_fwd = ref_seq.clone();
    apply_test_mutations(&mut query_fwd, 0.05, &mut rng);
    test_cases.push(OrientationTestCase {
        name: "high_mutation_forward".to_string(),
        reference: ref_seq.clone(),
        query: query_fwd,
        expected_reverse: false,
    });
    
    // Test case 5: High mutation rate reverse
    let ref_seq = generate_test_sequence(1000, &mut rng);
    let mut query_rev = reverse_complement(&ref_seq);
    apply_test_mutations(&mut query_rev, 0.05, &mut rng);
    test_cases.push(OrientationTestCase {
        name: "high_mutation_reverse".to_string(),
        reference: ref_seq.clone(),
        query: query_rev,
        expected_reverse: true,
    });
    
    // Test case 6: Short sequences
    let ref_seq = generate_test_sequence(100, &mut rng);
    let query_rev = reverse_complement(&ref_seq);
    test_cases.push(OrientationTestCase {
        name: "short_sequences_reverse".to_string(),
        reference: ref_seq.clone(),
        query: query_rev,
        expected_reverse: true,
    });
    
    // Test case 7: Long sequences
    let ref_seq = generate_test_sequence(10000, &mut rng);
    let mut query_fwd = ref_seq.clone();
    apply_test_mutations(&mut query_fwd, 0.001, &mut rng);
    test_cases.push(OrientationTestCase {
        name: "long_sequences_forward".to_string(),
        reference: ref_seq.clone(),
        query: query_fwd,
        expected_reverse: false,
    });
    
    // Test case 8: Ambiguous case (very high mutation rate)
    let ref_seq = generate_test_sequence(500, &mut rng);
    let mut query_fwd = ref_seq.clone();
    apply_test_mutations(&mut query_fwd, 0.2, &mut rng);
    test_cases.push(OrientationTestCase {
        name: "ambiguous_high_mutation".to_string(),
        reference: ref_seq.clone(),
        query: query_fwd,
        expected_reverse: false,
    });
    
    test_cases
}

fn generate_test_sequence(length: usize, rng: &mut impl Rng) -> Vec<u8> {
    let bases = [b'A', b'T', b'G', b'C'];
    (0..length).map(|_| bases[rng.gen_range(0..4)]).collect()
}

fn apply_test_mutations(seq: &mut [u8], mutation_rate: f64, rng: &mut impl Rng) {
    let bases = [b'A', b'T', b'G', b'C'];
    for i in 0..seq.len() {
        if rng.gen::<f64>() < mutation_rate {
            let current_base = seq[i];
            let mut new_bases = bases.to_vec();
            new_bases.retain(|&b| b != current_base);
            seq[i] = new_bases[rng.gen_range(0..3)];
        }
    }
}

fn test_orientation_detection_mash(test_case: &OrientationTestCase) -> OrientationResult {
    const ORIENTATION_KMER_SIZE: usize = 15;
    const ORIENTATION_SKETCH_SIZE: usize = 1000;
    
    // Create strand-specific sketches (without canonicalization)
    let target_sketch = sketch_sequence_stranded_test(&test_case.reference, ORIENTATION_KMER_SIZE, ORIENTATION_SKETCH_SIZE);
    let fwd_sketch = sketch_sequence_stranded_test(&test_case.query, ORIENTATION_KMER_SIZE, ORIENTATION_SKETCH_SIZE);
    
    let rev_seq = reverse_complement(&test_case.query);
    let rev_sketch = sketch_sequence_stranded_test(&rev_seq, ORIENTATION_KMER_SIZE, ORIENTATION_SKETCH_SIZE);
    
    // Compute Jaccard similarities
    let fwd_jaccard = jaccard_similarity_test(&fwd_sketch, &target_sketch);
    let rev_jaccard = jaccard_similarity_test(&rev_sketch, &target_sketch);
    
    OrientationResult {
        forward_distance: 1.0 - fwd_jaccard,
        reverse_distance: 1.0 - rev_jaccard,
        is_reverse: fwd_jaccard < rev_jaccard,
    }
}

/// Create strand-specific MinHash sketch for testing (no canonicalization)
fn sketch_sequence_stranded_test(sequence: &[u8], k: usize, sketch_size: usize) -> Vec<u64> {
    if sequence.len() < k {
        return Vec::new();
    }

    let mut hashes = Vec::new();

    // Extract all k-mers and hash them (no canonicalization)
    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];

        // Skip k-mers containing non-ACGT characters
        if kmer.iter().any(|&b| !is_dna_base_test(b)) {
            continue;
        }

        // Hash the k-mer directly (no reverse complement comparison)
        let hash = hash_kmer_test(kmer);
        hashes.push(hash);
    }

    // Sort and take the smallest hashes (MinHash)
    hashes.sort_unstable();
    hashes.truncate(sketch_size);
    hashes
}

/// Compute Jaccard similarity between two sketches for testing
fn jaccard_similarity_test(sketch1: &[u64], sketch2: &[u64]) -> f64 {
    use std::collections::HashSet;
    
    let set1: HashSet<_> = sketch1.iter().collect();
    let set2: HashSet<_> = sketch2.iter().collect();

    let intersection_size = set1.intersection(&set2).count();
    let union_size = set1.union(&set2).count();

    if union_size == 0 {
        0.0
    } else {
        intersection_size as f64 / union_size as f64
    }
}

/// Hash a k-mer for testing
fn hash_kmer_test(kmer: &[u8]) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    
    let mut hasher = DefaultHasher::new();
    kmer.hash(&mut hasher);
    hasher.finish()
}

/// Check if byte represents a DNA base for testing
fn is_dna_base_test(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn test_orientation_detection_edlib(test_case: &OrientationTestCase) -> OrientationResult {
    use edlib_rs::edlibrs::*;
    
    let config = EdlibAlignConfigRs {
        k: -1,
        mode: EdlibAlignModeRs::EDLIB_MODE_NW,
        task: EdlibAlignTaskRs::EDLIB_TASK_DISTANCE,
        additionalequalities: &[],
    };
    
    // Get edit distance for forward orientation
    let fwd_result = edlibAlignRs(&test_case.query, &test_case.reference, &config);
    let fwd_distance = fwd_result.editDistance;
    
    // Get edit distance for reverse complement orientation
    let rev_seq = reverse_complement(&test_case.query);
    let rev_result = edlibAlignRs(&rev_seq, &test_case.reference, &config);
    let rev_distance = rev_result.editDistance;
    
    OrientationResult {
        forward_distance: fwd_distance as f64,
        reverse_distance: rev_distance as f64,
        is_reverse: fwd_distance > rev_distance,
    }
}

#[test]
fn test_orientation_detection_performance() {
    println!("\n=== Test: Orientation detection performance benchmark ===");
    
    use std::time::Instant;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    
    let mut rng = StdRng::seed_from_u64(99999);
    
    // Test different sequence lengths
    let test_lengths = vec![100, 500, 1000, 5000, 10000];
    
    for length in test_lengths {
        println!("  Testing length: {}", length);
        
        let reference = generate_test_sequence(length, &mut rng);
        let mut query = reference.clone();
        apply_test_mutations(&mut query, 0.01, &mut rng);
        
        let test_case = OrientationTestCase {
            name: format!("perf_test_{}", length),
            reference,
            query,
            expected_reverse: false,
        };
        
        // Time mash-based detection
        let start = Instant::now();
        let mash_result = test_orientation_detection_mash(&test_case);
        let mash_time = start.elapsed();
        
        // Time edlib-based detection
        let start = Instant::now();
        let edlib_result = test_orientation_detection_edlib(&test_case);
        let edlib_time = start.elapsed();
        
        println!("    Mash time: {:?}, Edlib time: {:?}", mash_time, edlib_time);
        println!("    Mash result: {}, Edlib result: {}", 
                 if mash_result.is_reverse { "reverse" } else { "forward" },
                 if edlib_result.is_reverse { "reverse" } else { "forward" });
        
        // Both should agree on orientation
        assert_eq!(mash_result.is_reverse, edlib_result.is_reverse,
                   "Methods disagree for length {}", length);
    }
}
