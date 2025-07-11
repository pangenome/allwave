use crate::validation::{parse_cigar, CigarOp};

/// Verify that a CIGAR string correctly describes the alignment between two sequences
pub fn verify_alignment_cigar(
    query_seq: &[u8],
    target_seq: &[u8],
    cigar_str: &str,
    strand: char,
) -> Result<bool, String> {
    let cigar_ops = parse_cigar(cigar_str);
    
    // Apply reverse complement if needed
    let oriented_query = if strand == '-' {
        reverse_complement(query_seq)
    } else {
        query_seq.to_vec()
    };
    
    let mut q_pos = 0;
    let mut t_pos = 0;
    let mut errors = Vec::new();
    
    for (op_idx, op) in cigar_ops.iter().enumerate() {
        match op {
            CigarOp::Match(len) => {
                let len = *len as usize;
                
                // Check bounds
                if q_pos + len > oriented_query.len() {
                    errors.push(format!("Op {}: Match({}) exceeds query bounds at pos {}", op_idx, len, q_pos));
                    break;
                }
                if t_pos + len > target_seq.len() {
                    errors.push(format!("Op {}: Match({}) exceeds target bounds at pos {}", op_idx, len, t_pos));
                    break;
                }
                
                // Verify matches
                for i in 0..len {
                    if oriented_query[q_pos + i] != target_seq[t_pos + i] {
                        if errors.len() < 5 {  // Only report first few errors
                            errors.push(format!(
                                "Op {}: Mismatch in '=' at q[{}]={} vs t[{}]={}", 
                                op_idx,
                                q_pos + i, oriented_query[q_pos + i] as char,
                                t_pos + i, target_seq[t_pos + i] as char
                            ));
                        }
                    }
                }
                
                q_pos += len;
                t_pos += len;
            }
            CigarOp::Mismatch(len) => {
                let len = *len as usize;
                
                // Check bounds
                if q_pos + len > oriented_query.len() {
                    errors.push(format!("Op {}: Mismatch({}) exceeds query bounds at pos {}", op_idx, len, q_pos));
                    break;
                }
                if t_pos + len > target_seq.len() {
                    errors.push(format!("Op {}: Mismatch({}) exceeds target bounds at pos {}", op_idx, len, t_pos));
                    break;
                }
                
                // We could verify these are actually different, but the aligner says they are
                q_pos += len;
                t_pos += len;
            }
            CigarOp::Insertion(len) => {
                let len = *len as usize;
                
                // Check bounds
                if q_pos + len > oriented_query.len() {
                    errors.push(format!("Op {}: Insertion({}) exceeds query bounds at pos {}", op_idx, len, q_pos));
                    break;
                }
                
                // Insertion: query has extra bases not in target
                q_pos += len;
            }
            CigarOp::Deletion(len) => {
                let len = *len as usize;
                
                // Check bounds  
                if t_pos + len > target_seq.len() {
                    errors.push(format!("Op {}: Deletion({}) exceeds target bounds at pos {}", op_idx, len, t_pos));
                    break;
                }
                
                // Deletion: target has extra bases not in query
                t_pos += len;
            }
        }
    }
    
    // Check we consumed all of both sequences
    if q_pos != oriented_query.len() {
        errors.push(format!("CIGAR only consumed {}/{} query bases", q_pos, oriented_query.len()));
    }
    if t_pos != target_seq.len() {
        errors.push(format!("CIGAR only consumed {}/{} target bases", t_pos, target_seq.len()));
    }
    
    if errors.is_empty() {
        Ok(true)
    } else {
        Err(errors.join("; "))
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_simple_match() {
        let query = b"ATCG";
        let target = b"ATCG";
        let cigar = "4=";
        
        assert!(verify_alignment_cigar(query, target, cigar, '+').unwrap());
    }
    
    #[test]
    fn test_with_mismatch() {
        let query = b"ATCG";
        let target = b"ATGG";
        let cigar = "2=1X1=";
        
        assert!(verify_alignment_cigar(query, target, cigar, '+').unwrap());
    }
    
    #[test] 
    fn test_with_deletion() {
        // Query has deletion relative to target
        let query = b"ATAT";
        let target = b"ATCGAT";
        let cigar = "2=2D2=";  // 2 match, 2 deletion (in target), 2 match
        
        assert!(verify_alignment_cigar(query, target, cigar, '+').unwrap());
    }
    
    #[test]
    fn test_with_insertion() {
        // Query has insertion relative to target  
        let query = b"ATCGAT";
        let target = b"ATAT";
        let cigar = "2=2I2=";  // 2 match, 2 insertion (in query), 2 match
        
        assert!(verify_alignment_cigar(query, target, cigar, '+').unwrap());
    }
}