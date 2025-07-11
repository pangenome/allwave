use allwave::test_framework::*;
use allwave::validation::*;
use std::process::Command;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== AllWave Test Suite ===\n");
    
    // Test 1: Basic mutations (SNPs and small indels)
    println!("Test 1: Basic mutations (SNPs and small indels)");
    let mut seq_config = SequenceConfig::default();
    seq_config.length = 5000;
    seq_config.seed = 42;
    
    let mut mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.01;
    mut_config.indel_rate = 0.001;
    mut_config.microsatellite_rate = 0.0;
    mut_config.cnv_rate = 0.0;
    mut_config.seed = 123;
    
    let test_case = create_test_case("test1_basic", &seq_config, &mut_config);
    run_test(&test_case)?;
    
    // Test 2: Microsatellite mutations
    println!("\nTest 2: Microsatellite expansions and contractions");
    seq_config.seed = 43;
    mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.001;
    mut_config.microsatellite_rate = 0.01;
    mut_config.cnv_rate = 0.0;
    mut_config.seed = 124;
    
    // Add some microsatellites to the sequence
    let test_case = create_microsatellite_test();
    run_test(&test_case)?;
    
    // Test 3: CNV mutations
    println!("\nTest 3: Copy Number Variations (CNVs)");
    seq_config.length = 20000;
    seq_config.seed = 44;
    mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.0001;
    mut_config.cnv_rate = 0.0001;
    mut_config.cnv_min_size = 1000;
    mut_config.cnv_max_size = 3000;
    mut_config.seed = 125;
    
    let test_case = create_test_case("test3_cnv", &seq_config, &mut_config);
    run_test(&test_case)?;
    
    // Test 4: Complex mutations (all types)
    println!("\nTest 4: Complex mutations (all types combined)");
    seq_config.length = 10000;
    seq_config.seed = 45;
    mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.005;
    mut_config.indel_rate = 0.0005;
    mut_config.microsatellite_rate = 0.0001;
    mut_config.cnv_rate = 0.00005;
    mut_config.seed = 126;
    
    let test_case = create_test_case("test4_complex", &seq_config, &mut_config);
    run_test(&test_case)?;
    
    // Test 5: Edge cases
    println!("\nTest 5: Edge cases (high mutation rate)");
    seq_config.length = 3000;
    seq_config.seed = 46;
    mut_config = MutationConfig::default();
    mut_config.snp_rate = 0.05;
    mut_config.indel_rate = 0.01;
    mut_config.seed = 127;
    
    let test_case = create_test_case("test5_edge", &seq_config, &mut_config);
    run_test(&test_case)?;
    
    println!("\n=== All tests completed ===");
    Ok(())
}

fn run_test(test_case: &TestCase) -> Result<(), Box<dyn std::error::Error>> {
    // Print mutation statistics
    let stats = test_case.get_mutation_stats();
    println!("  Mutations introduced:");
    for (mutation_type, count) in &stats {
        println!("    {mutation_type:?}: {count}");
    }
    
    // Write test case to FASTA
    let fasta_path = format!("/tmp/{}.fa", test_case.name);
    test_case.write_fasta(&fasta_path)?;
    
    // Run allwave
    let output = Command::new("./target/debug/allwave")
        .arg("--input")
        .arg(&fasta_path)
        .output()?;
    
    if !output.status.success() {
        eprintln!("  ERROR: allwave failed to run");
        eprintln!("  stderr: {}", String::from_utf8_lossy(&output.stderr));
        return Ok(());
    }
    
    // Parse and validate results
    let paf_output = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = paf_output.lines().collect();
    
    if lines.is_empty() {
        eprintln!("  ERROR: No alignment output");
        return Ok(());
    }
    
    // Validate the first alignment (reference vs query)
    let result = run_and_validate_test(test_case, lines[0])?;
    result.print_summary();
    
    // Clean up
    std::fs::remove_file(&fasta_path)?;
    
    Ok(())
}

fn create_microsatellite_test() -> TestCase {
    use rand::{Rng, SeedableRng};
    use rand::rngs::StdRng;
    
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
        name: "test2_microsatellite".to_string(),
        reference,
        query,
        mutations,
    }
}