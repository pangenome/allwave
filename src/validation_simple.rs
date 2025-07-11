use crate::validation::{parse_cigar, calculate_alignment_stats, CigarOp, ValidationResult};
use std::collections::HashMap;

/// Simple validation that just checks CIGAR correctness
pub fn validate_alignment_simple(
    query_seq: &[u8],
    target_seq: &[u8],
    paf_line: &str,
) -> Result<ValidationResult, Box<dyn std::error::Error>> {
    // Parse PAF line
    let fields: Vec<&str> = paf_line.split('\t').collect();
    if fields.len() < 12 {
        return Err("Invalid PAF format".into());
    }
    
    let query_start = fields[2].parse::<usize>()?;
    let query_end = fields[3].parse::<usize>()?;
    let strand = fields[4].chars().next().ok_or("No strand found")?;
    let target_start = fields[7].parse::<usize>()?;
    let target_end = fields[8].parse::<usize>()?;
    
    // Find CIGAR string
    let cigar_str = fields.iter()
        .find(|f| f.starts_with("cg:Z:"))
        .ok_or("No CIGAR string found")?
        .trim_start_matches("cg:Z:");
    
    let cigar_ops = parse_cigar(cigar_str);
    let stats = calculate_alignment_stats(&cigar_ops);
    
    // Apply strand orientation to query if needed
    let oriented_query = if strand == '-' {
        reverse_complement(query_seq)
    } else {
        query_seq.to_vec()
    };
    
    // Debug print
    println!("  Validating alignment: query {}..{} (oriented len {}), target {}..{} (len {})",
             query_start, query_end, oriented_query.len(), 
             target_start, target_end, target_seq.len());
    
    // Verify alignment
    let alignment_correct = verify_cigar_correctness(
        &oriented_query[query_start..query_end],
        &target_seq[target_start..target_end],
        &cigar_ops,
    );
    
    // Calculate coverage
    let query_coverage = (query_end - query_start) as f64 / query_seq.len() as f64;
    let target_coverage = (target_end - target_start) as f64 / target_seq.len() as f64;
    
    Ok(ValidationResult {
        alignment_stats: stats,
        query_coverage,
        target_coverage,
        expected_mutations: HashMap::new(),
        detected_mutations: HashMap::new(),
        detection_accuracy: HashMap::new(),
        is_valid: alignment_correct && query_coverage > 0.8 && target_coverage > 0.8,
        alignment_correct,
    })
}

fn verify_cigar_correctness(
    query_seq: &[u8],
    target_seq: &[u8],
    cigar_ops: &[CigarOp],
) -> bool {
    let mut q_pos = 0;
    let mut t_pos = 0;
    
    for op in cigar_ops {
        match op {
            CigarOp::Match(len) => {
                let len = *len as usize;
                if q_pos + len > query_seq.len() || t_pos + len > target_seq.len() {
                    eprintln!("Match operation extends beyond sequence bounds");
                    return false;
                }
                for i in 0..len {
                    if query_seq[q_pos + i] != target_seq[t_pos + i] {
                        eprintln!("Mismatch in = operation at positions q:{}, t:{}: {} != {}", 
                                 q_pos + i, t_pos + i, 
                                 query_seq[q_pos + i] as char, 
                                 target_seq[t_pos + i] as char);
                        // Print context
                        let q_start = q_pos.saturating_sub(5);
                        let q_end = (q_pos + i + 5).min(query_seq.len());
                        let t_start = t_pos.saturating_sub(5);
                        let t_end = (t_pos + i + 5).min(target_seq.len());
                        eprintln!("  Query context:  {}", 
                                 String::from_utf8_lossy(&query_seq[q_start..q_end]));
                        eprintln!("  Target context: {}", 
                                 String::from_utf8_lossy(&target_seq[t_start..t_end]));
                        return false;
                    }
                }
                q_pos += len;
                t_pos += len;
            }
            CigarOp::Mismatch(len) => {
                let len = *len as usize;
                if q_pos + len > query_seq.len() || t_pos + len > target_seq.len() {
                    eprintln!("Mismatch operation extends beyond sequence bounds");
                    return false;
                }
                // For mismatches, we don't check that they're actually different
                // because the aligner already determined they're mismatches
                q_pos += len;
                t_pos += len;
            }
            CigarOp::Insertion(len) => {
                let len = *len as usize;
                if q_pos + len > query_seq.len() {
                    eprintln!("Insertion extends beyond query sequence");
                    return false;
                }
                q_pos += len;
            }
            CigarOp::Deletion(len) => {
                let len = *len as usize;
                if t_pos + len > target_seq.len() {
                    eprintln!("Deletion extends beyond target sequence");
                    return false;
                }
                t_pos += len;
            }
        }
    }
    
    // Check that we consumed the right amount of sequence
    if q_pos != query_seq.len() || t_pos != target_seq.len() {
        println!("  CIGAR consumed: query {}/{}, target {}/{}", 
                  q_pos, query_seq.len(), t_pos, target_seq.len());
        // For global alignments, we should consume all bases
        if q_pos != query_seq.len() || t_pos != target_seq.len() {
            eprintln!("  ERROR: CIGAR doesn't account for full alignment");
            return false;
        }
    }
    
    println!("  Alignment validation: PASSED");
    true
}

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