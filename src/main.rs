use bio::io::fasta;
use clap::Parser;
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

fn edit_distance(seq1: &[u8], seq2: &[u8], max_dist: usize) -> usize {
    let len1 = seq1.len().min(1000); // Sample first 1000bp for speed
    let len2 = seq2.len().min(1000);

    let mut current = vec![0; len2 + 1];
    let mut previous = vec![0; len2 + 1];

    // Initialize first row
    for (j, prev) in previous.iter_mut().enumerate().take(len2 + 1) {
        *prev = j;
    }

    for i in 1..=len1 {
        current[0] = i;
        let mut min_in_row = i;

        for j in 1..=len2 {
            let cost = if seq1[i - 1] == seq2[j - 1] { 0 } else { 1 };
            current[j] = std::cmp::min(
                std::cmp::min(previous[j] + 1, current[j - 1] + 1),
                previous[j - 1] + cost,
            );
            min_in_row = min_in_row.min(current[j]);
        }

        // Early termination if all values exceed max_dist
        if min_in_row > max_dist {
            return max_dist + 1;
        }

        std::mem::swap(&mut current, &mut previous);
    }

    previous[len2]
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
                    '=' | 'X' => {
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

            // Check both orientations using edit distance
            let fwd_dist = edit_distance(query_seq, target_seq, 100);
            let rev_seq = reverse_complement(query_seq);
            let rev_dist = edit_distance(&rev_seq, target_seq, 100);

            // Choose the better orientation
            let (seq_to_align, strand) = if fwd_dist <= rev_dist {
                (query_seq.to_vec(), '+')
            } else {
                (rev_seq, '-')
            };

            // Perform alignment
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
            
            let (cigar_query_len, cigar_target_len) = parse_cigar_lengths(&alignment.cigar);
            
            // For reverse strand alignments, we need to adjust the query coordinates
            let (query_start, query_end) = if strand == '-' {
                (query_seq.len().saturating_sub(cigar_query_len), query_seq.len())
            } else {
                (0, cigar_query_len.min(query_seq.len()))
            };
            
            let record = PafRecord {
                query_name,
                query_len: query_seq.len(),
                query_start,
                query_end,
                query_strand: strand,
                target_name,
                target_len: target_seq.len(),
                target_start: 0,
                target_end: cigar_target_len.min(target_seq.len()),
                alignment: &alignment,
            };
            
            let mut output = output.lock().unwrap();
            write_paf_record(&mut *output, &record)?;
        }
    }

    Ok(())
}
