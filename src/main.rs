use bio::io::fasta;
use clap::Parser;
use edlib_rs::edlibrs::*;
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

mod wfa;
use wfa::{align_sequences, AlignmentMode, AlignmentResult, Penalties};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input FASTA file
    #[arg(short, long)]
    input: PathBuf,

    /// Output PAF file (default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Alignment scores: match,mismatch,gap_open,gap_ext[,gap_open2,gap_ext2]
    #[arg(short, long, default_value = "0,5,8,2,24,1")]
    scores: String,
    
    /// Number of threads to use for parallel processing
    #[arg(short, long, default_value = "1")]
    threads: usize,
}

fn parse_scores(
    scores_str: &str,
) -> Result<(AlignmentMode, Penalties), Box<dyn std::error::Error>> {
    let scores: Vec<i32> = scores_str
        .split(',')
        .map(|s| s.trim().parse::<i32>())
        .collect::<Result<Vec<_>, _>>()?;

    match scores.len() {
        4 => {
            // Edit distance or single-piece affine
            let penalties = Penalties {
                mismatch: scores[1],
                gap_opening1: scores[2],
                gap_extension1: scores[3],
                gap_opening2: scores[2], // Same as gap_opening1 for single-piece
                gap_extension2: scores[3], // Same as gap_extension1 for single-piece
            };

            // If gap_open == gap_ext == mismatch, it's edit distance
            if scores[2] == scores[3] && scores[2] == scores[1] {
                Ok((AlignmentMode::EditDistance, penalties))
            } else {
                Ok((AlignmentMode::SinglePieceAffine, penalties))
            }
        }
        6 => {
            // Two-piece affine
            let penalties = Penalties {
                mismatch: scores[1],
                gap_opening1: scores[2],
                gap_extension1: scores[3],
                gap_opening2: scores[4],
                gap_extension2: scores[5],
            };
            Ok((AlignmentMode::TwoPieceAffine, penalties))
        }
        _ => Err(format!(
            "Invalid number of scores: {}. Expected 4 or 6 values.",
            scores.len()
        )
        .into()),
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


struct PafRecord<'a> {
    query_name: &'a str,
    query_len: usize,
    query_start: usize,
    query_end: usize,
    query_strand: char,
    target_name: &'a str,
    target_len: usize,
    target_start: usize,
    target_end: usize,
    alignment: &'a AlignmentResult,
}

fn parse_cigar_lengths(cigar: &str) -> (usize, usize) {
    let mut query_len = 0;
    let mut target_len = 0;
    let mut num_str = String::new();
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            if let Ok(count) = num_str.parse::<usize>() {
                match ch {
                    '=' | 'X' | 'M' => {
                        query_len += count;
                        target_len += count;
                    }
                    'I' => {
                        query_len += count;
                    }
                    'D' => {
                        target_len += count;
                    }
                    _ => {}
                }
            }
            num_str.clear();
        }
    }
    
    (query_len, target_len)
}

#[allow(dead_code)]
fn truncate_cigar_to_bounds(cigar: &str, max_query: usize, max_target: usize) -> String {
    let mut result = String::new();
    let mut query_pos = 0;
    let mut target_pos = 0;
    let mut num_str = String::new();
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            if let Ok(count) = num_str.parse::<usize>() {
                let mut actual_count = count;
                
                match ch {
                    '=' | 'X' => {
                        // Check if this operation would exceed bounds
                        let query_space = max_query.saturating_sub(query_pos);
                        let target_space = max_target.saturating_sub(target_pos);
                        actual_count = actual_count.min(query_space).min(target_space);
                        
                        if actual_count > 0 {
                            result.push_str(&format!("{}{}", actual_count, ch));
                            query_pos += actual_count;
                            target_pos += actual_count;
                        }
                        
                        // If we've reached the limit, stop processing
                        if query_pos >= max_query || target_pos >= max_target {
                            break;
                        }
                    }
                    'I' => {
                        let query_space = max_query.saturating_sub(query_pos);
                        actual_count = actual_count.min(query_space);
                        
                        if actual_count > 0 {
                            result.push_str(&format!("{}{}", actual_count, ch));
                            query_pos += actual_count;
                        }
                        
                        if query_pos >= max_query {
                            break;
                        }
                    }
                    'D' => {
                        let target_space = max_target.saturating_sub(target_pos);
                        actual_count = actual_count.min(target_space);
                        
                        if actual_count > 0 {
                            result.push_str(&format!("{}{}", actual_count, ch));
                            target_pos += actual_count;
                        }
                        
                        if target_pos >= max_target {
                            break;
                        }
                    }
                    _ => {
                        // Unknown operation, include as-is
                        result.push_str(&format!("{}{}", count, ch));
                    }
                }
            }
            num_str.clear();
        }
    }
    
    result
}

fn write_paf_record(out: &mut dyn Write, record: &PafRecord) -> io::Result<()> {
    let query_aligned_len = record.query_end - record.query_start;
    let target_aligned_len = record.target_end - record.target_start;
    let block_len = target_aligned_len.max(query_aligned_len);

    // Calculate identity
    let identity = record.alignment.matches as f64 / record.alignment.alignment_length as f64;

    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgi:f:{:.6}\tcg:Z:{}",
        record.query_name,
        record.query_len,
        record.query_start,
        record.query_end,
        record.query_strand,
        record.target_name,
        record.target_len,
        record.target_start,
        record.target_end,
        record.alignment.matches,
        block_len,
        60, // mapping quality placeholder
        identity,
        record.alignment.cigar
    )?;

    writeln!(out)?;
    Ok(())
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    
    // Set the number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    // Read FASTA file (handle both plain and gzipped)
    let mut sequences = Vec::new();
    
    // Check if file is gzipped by extension and read accordingly
    if args.input.extension().and_then(|s| s.to_str()) == Some("gz") {
        // Handle gzipped file
        let file = File::open(&args.input)?;
        let decoder = GzDecoder::new(file);
        let fasta_reader = fasta::Reader::new(decoder);
        
        for result in fasta_reader.records() {
            let record = result?;
            sequences.push((record.id().to_string(), record.seq().to_vec()));
        }
    } else {
        // Handle plain file
        let file = File::open(&args.input)?;
        let fasta_reader = fasta::Reader::new(file);
        
        for result in fasta_reader.records() {
            let record = result?;
            sequences.push((record.id().to_string(), record.seq().to_vec()));
        }
    }

    // Parse alignment parameters
    let (alignment_mode, penalties) = parse_scores(&args.scores)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;

    // Create all alignment pairs
    let mut pairs = Vec::new();
    for i in 0..sequences.len() {
        for j in 0..sequences.len() {
            if i != j {
                pairs.push((i, j));
            }
        }
    }

    // Process alignments in parallel
    let results: Vec<_> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let (query_name, query_seq) = &sequences[i];
            let (target_name, target_seq) = &sequences[j];

            // Check both orientations using edlib for fast edit distance
            let config = EdlibAlignConfigRs {
                k: -1,  // No limit on edit distance
                mode: EdlibAlignModeRs::EDLIB_MODE_NW,  // Global alignment
                task: EdlibAlignTaskRs::EDLIB_TASK_DISTANCE,  // Just get distance, no alignment path
                additionalequalities: &[],  // No additional equalities
            };
            
            // Get edit distance for forward orientation
            let fwd_result = edlibAlignRs(query_seq, target_seq, &config);
            let fwd_distance = fwd_result.editDistance;
            
            // Get edit distance for reverse complement orientation
            let rev_seq = reverse_complement(query_seq);
            let rev_result = edlibAlignRs(&rev_seq, target_seq, &config);
            let rev_distance = rev_result.editDistance;
            
            // Choose the better orientation based on edit distance
            let (seq_to_align, strand) = if fwd_distance <= rev_distance {
                (query_seq.to_vec(), '+')
            } else {
                (rev_seq, '-')
            };

            // Now perform the actual alignment with user-specified parameters
            match align_sequences(&seq_to_align, target_seq, &penalties, alignment_mode) {
                Ok(alignment) => {
                    Some((i, j, strand, alignment))
                }
                Err(e) => {
                    eprintln!("Warning: Failed to align {} vs {}: {}", query_name, target_name, e);
                    None
                }
            }
        })
        .collect();

    // Write results to output
    let output = Arc::new(Mutex::new(if let Some(output_path) = args.output {
        Box::new(File::create(output_path)?) as Box<dyn Write + Send>
    } else {
        Box::new(io::stdout()) as Box<dyn Write + Send>
    }));

    for result in results {
        if let Some((i, j, strand, alignment)) = result {
            let (query_name, query_seq) = &sequences[i];
            let (target_name, target_seq) = &sequences[j];
            
            let (query_end, target_end) = parse_cigar_lengths(&alignment.cigar);
            
            let record = PafRecord {
                query_name,
                query_len: query_seq.len(),
                query_start: 0,
                query_end,
                query_strand: strand,
                target_name,
                target_len: target_seq.len(),
                target_start: 0,
                target_end,
                alignment: &alignment,
            };
            
            let mut output = output.lock().unwrap();
            write_paf_record(&mut *output, &record)?;
        }
    }

    Ok(())
}
