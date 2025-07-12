use crate::test_framework::{MutationType, TestCase};
use std::collections::HashMap;

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

/// Represents a CIGAR operation
#[derive(Debug, Clone, PartialEq)]
pub enum CigarOp {
    Match(u32),
    Mismatch(u32),
    Insertion(u32),
    Deletion(u32),
}

/// Parse CIGAR string into operations
pub fn parse_cigar(cigar_str: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num_str = String::new();

    for ch in cigar_str.chars() {
        if ch.is_numeric() {
            num_str.push(ch);
        } else if let Ok(count) = num_str.parse::<u32>() {
            let op = match ch {
                '=' | 'M' => CigarOp::Match(count), // M can be match or mismatch
                'X' => CigarOp::Mismatch(count),
                'I' => CigarOp::Insertion(count),
                'D' => CigarOp::Deletion(count),
                _ => continue,
            };
            ops.push(op);
            num_str.clear();
        }
    }

    ops
}

/// Calculate alignment statistics from CIGAR
pub fn calculate_alignment_stats(cigar_ops: &[CigarOp]) -> AlignmentStats {
    let mut stats = AlignmentStats::default();

    for op in cigar_ops {
        match op {
            CigarOp::Match(len) => {
                stats.matches += *len as usize;
                stats.alignment_length += *len as usize;
            }
            CigarOp::Mismatch(len) => {
                stats.mismatches += *len as usize;
                stats.alignment_length += *len as usize;
            }
            CigarOp::Insertion(len) => {
                stats.insertions += *len as usize;
                stats.gap_opens += 1;
            }
            CigarOp::Deletion(len) => {
                stats.deletions += *len as usize;
                stats.gap_opens += 1;
            }
        }
    }

    stats.identity = if stats.alignment_length > 0 {
        stats.matches as f64 / stats.alignment_length as f64
    } else {
        0.0
    };

    stats
}

#[derive(Debug, Default)]
pub struct AlignmentStats {
    pub matches: usize,
    pub mismatches: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub gap_opens: usize,
    pub alignment_length: usize,
    pub identity: f64,
}

/// Apply CIGAR operations to verify alignment correctness
pub fn verify_cigar_alignment(
    query_seq: &[u8],
    target_seq: &[u8],
    cigar_ops: &[CigarOp],
    query_start: usize,
    target_start: usize,
) -> Result<bool, String> {
    let mut query_pos = query_start;
    let mut target_pos = target_start;

    for (i, op) in cigar_ops.iter().enumerate() {
        match op {
            CigarOp::Match(len) => {
                let len = *len as usize;
                if query_pos + len > query_seq.len() || target_pos + len > target_seq.len() {
                    return Err(format!(
                        "CIGAR operation {i} extends beyond sequence bounds"
                    ));
                }
                // Note: 'M' operations can contain both matches and mismatches
                // so we don't validate exact matches here
                query_pos += len;
                target_pos += len;
            }
            CigarOp::Mismatch(len) => {
                let len = *len as usize;
                if query_pos + len > query_seq.len() || target_pos + len > target_seq.len() {
                    return Err(format!(
                        "CIGAR operation {i} extends beyond sequence bounds"
                    ));
                }
                for j in 0..len {
                    if query_seq[query_pos + j] == target_seq[target_pos + j] {
                        return Err(format!(
                            "Match in 'X' operation at query pos {} ({}), target pos {} ({})",
                            query_pos + j,
                            query_seq[query_pos + j] as char,
                            target_pos + j,
                            target_seq[target_pos + j] as char
                        ));
                    }
                }
                query_pos += len;
                target_pos += len;
            }
            CigarOp::Insertion(len) => {
                let len = *len as usize;
                if query_pos + len > query_seq.len() {
                    return Err("Insertion extends beyond query sequence bounds".to_string());
                }
                query_pos += len;
            }
            CigarOp::Deletion(len) => {
                let len = *len as usize;
                if target_pos + len > target_seq.len() {
                    return Err("Deletion extends beyond target sequence bounds".to_string());
                }
                target_pos += len;
            }
        }
    }

    Ok(true)
}

/// Validate alignment against known mutations
pub fn validate_alignment(
    test_case: &TestCase,
    cigar_str: &str,
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
    strand: char,
) -> ValidationResult {
    let cigar_ops = parse_cigar(cigar_str);
    let stats = calculate_alignment_stats(&cigar_ops);

    // Get the appropriate query sequence based on strand
    let query_seq = if strand == '-' {
        reverse_complement(&test_case.query)
    } else {
        test_case.query.clone()
    };

    // Verify CIGAR alignment correctness
    let alignment_correct = if query_start < query_seq.len()
        && query_end <= query_seq.len()
        && target_start < test_case.reference.len()
        && target_end <= test_case.reference.len()
    {
        match verify_cigar_alignment(
            &query_seq[query_start..query_end],
            &test_case.reference[target_start..target_end],
            &cigar_ops,
            0,
            0,
        ) {
            Ok(_) => true,
            Err(e) => {
                eprintln!("Alignment verification failed: {e}");
                false
            }
        }
    } else {
        eprintln!(
            "Alignment coordinates out of bounds: query {}..{} (len {}), target {}..{} (len {})",
            query_start,
            query_end,
            query_seq.len(),
            target_start,
            target_end,
            test_case.reference.len()
        );
        false
    };

    // Check if alignment covers the expected regions
    let query_coverage = (query_end - query_start) as f64 / test_case.query.len() as f64;
    let target_coverage = (target_end - target_start) as f64 / test_case.reference.len() as f64;

    // Count mutations by type
    let mut mutation_counts = HashMap::new();
    for mutation in &test_case.mutations {
        *mutation_counts
            .entry(mutation.mutation_type.clone())
            .or_insert(0) += 1;
    }

    // Detect mutations from CIGAR
    let detected_mutations = detect_mutations_from_cigar(&cigar_ops);

    // Compare expected vs detected mutations
    let mut detection_accuracy = HashMap::new();
    for (mutation_type, expected_count) in &mutation_counts {
        let detected_count = detected_mutations.get(mutation_type).unwrap_or(&0);
        let accuracy = if *expected_count > 0 {
            *detected_count as f64 / *expected_count as f64
        } else {
            1.0
        };
        detection_accuracy.insert(mutation_type.clone(), accuracy);
    }

    ValidationResult {
        alignment_stats: stats,
        query_coverage,
        target_coverage,
        expected_mutations: mutation_counts,
        detected_mutations,
        detection_accuracy,
        is_valid: query_coverage > 0.95 && target_coverage > 0.95 && alignment_correct,
        alignment_correct,
    }
}

/// Detect mutation types from CIGAR operations
fn detect_mutations_from_cigar(cigar_ops: &[CigarOp]) -> HashMap<MutationType, usize> {
    let mut mutations = HashMap::new();
    mutations.insert(MutationType::SNP, 0);
    mutations.insert(MutationType::Insertion, 0);
    mutations.insert(MutationType::Deletion, 0);

    for op in cigar_ops {
        match op {
            CigarOp::Mismatch(count) => {
                *mutations.get_mut(&MutationType::SNP).unwrap() += *count as usize;
            }
            CigarOp::Insertion(count) => {
                *mutations.get_mut(&MutationType::Insertion).unwrap() += 1;
                // Large insertions might be CNV duplications
                if *count >= 1000 {
                    *mutations.entry(MutationType::CNVDuplication).or_insert(0) += 1;
                }
            }
            CigarOp::Deletion(count) => {
                *mutations.get_mut(&MutationType::Deletion).unwrap() += 1;
                // Large deletions might be CNV deletions
                if *count >= 1000 {
                    *mutations.entry(MutationType::CNVDeletion).or_insert(0) += 1;
                }
            }
            _ => {}
        }
    }

    mutations
}

#[derive(Debug)]
pub struct ValidationResult {
    pub alignment_stats: AlignmentStats,
    pub query_coverage: f64,
    pub target_coverage: f64,
    pub expected_mutations: HashMap<MutationType, usize>,
    pub detected_mutations: HashMap<MutationType, usize>,
    pub detection_accuracy: HashMap<MutationType, f64>,
    pub is_valid: bool,
    pub alignment_correct: bool,
}

impl ValidationResult {
    pub fn print_summary(&self) {
        println!("=== Alignment Validation Summary ===");
        println!("Alignment Stats:");
        println!("  Identity: {:.2}%", self.alignment_stats.identity * 100.0);
        println!("  Matches: {}", self.alignment_stats.matches);
        println!("  Mismatches: {}", self.alignment_stats.mismatches);
        println!("  Insertions: {}", self.alignment_stats.insertions);
        println!("  Deletions: {}", self.alignment_stats.deletions);
        println!("  Gap opens: {}", self.alignment_stats.gap_opens);
        println!("\nCoverage:");
        println!("  Query coverage: {:.2}%", self.query_coverage * 100.0);
        println!("  Target coverage: {:.2}%", self.target_coverage * 100.0);
        println!("\nMutation Detection:");

        for mutation_type in [
            MutationType::SNP,
            MutationType::Insertion,
            MutationType::Deletion,
            MutationType::CNVDuplication,
            MutationType::CNVDeletion,
        ] {
            let expected = self.expected_mutations.get(&mutation_type).unwrap_or(&0);
            let detected = self.detected_mutations.get(&mutation_type).unwrap_or(&0);
            let accuracy = self.detection_accuracy.get(&mutation_type).unwrap_or(&0.0);

            println!(
                "  {:?}: expected={}, detected={}, accuracy={:.2}%",
                mutation_type,
                expected,
                detected,
                accuracy * 100.0
            );
        }

        println!(
            "\nAlignment correctness: {}",
            if self.alignment_correct {
                "CORRECT"
            } else {
                "INCORRECT"
            }
        );
        println!(
            "Validation: {}",
            if self.is_valid { "PASSED" } else { "FAILED" }
        );
    }
}

/// Run alignment and validate results
pub fn run_and_validate_test(
    test_case: &TestCase,
    paf_line: &str,
) -> Result<ValidationResult, Box<dyn std::error::Error>> {
    // Parse PAF line
    let fields: Vec<&str> = paf_line.split('\t').collect();
    if fields.len() < 12 {
        return Err("Invalid PAF format".into());
    }

    let query_name = fields[0];
    let _target_name = fields[5];
    let query_start = fields[2].parse::<usize>()?;
    let query_end = fields[3].parse::<usize>()?;
    let strand = fields[4].chars().next().ok_or("No strand found")?;
    let target_start = fields[7].parse::<usize>()?;
    let target_end = fields[8].parse::<usize>()?;

    // Determine which sequence is which based on names
    // In our test cases, we write "test_name_reference" and "test_name_query"
    let (_actual_query_seq, _actual_target_seq) = if query_name.ends_with("_reference") {
        (&test_case.reference, &test_case.query)
    } else if query_name.ends_with("_query") {
        (&test_case.query, &test_case.reference)
    } else {
        // Default assumption: first sequence in PAF is reference
        (&test_case.reference, &test_case.query)
    };

    // Find CIGAR string
    let cigar_str = fields
        .iter()
        .find(|f| f.starts_with("cg:Z:"))
        .ok_or("No CIGAR string found")?
        .trim_start_matches("cg:Z:");

    Ok(validate_alignment(
        test_case,
        cigar_str,
        query_start,
        query_end,
        target_start,
        target_end,
        strand,
    ))
}
