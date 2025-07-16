use allwave::{AllPairIterator, Sequence, alignment_to_paf, parse_scores, SparsificationStrategy};
use bio::io::fasta;
use clap::Parser;
use faigz_rs::{FastaFormat, FastaIndex, FastaReader};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

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
    
    /// Sparsification strategy: none, random:<fraction>, or auto
    #[arg(short = 'p', long, default_value = "none")]
    sparsification: String,
    
    /// Disable progress bar output
    #[arg(long)]
    no_progress: bool,
}





fn main() -> io::Result<()> {
    let args = Args::parse();
    
    // Set the number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .map_err(io::Error::other)?;

    // Parse sparsification strategy
    let sparsification = match args.sparsification.as_str() {
        "none" => SparsificationStrategy::None,
        "auto" => SparsificationStrategy::Auto,
        s if s.starts_with("random:") => {
            let fraction_str = &s[7..];
            let fraction: f64 = fraction_str.parse()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid random fraction"))?;
            SparsificationStrategy::Random(fraction)
        }
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid sparsification strategy")),
    };
    
    // Read FASTA file (handle both plain and gzipped)
    let mut sequences = Vec::new();
    
    // Check if file is gzipped by extension and read accordingly
    if args.input.extension().and_then(|s| s.to_str()) == Some("gz") {
        // Use faigz-rs for gzipped files
        let index = FastaIndex::new(args.input.to_str().unwrap(), FastaFormat::Fasta)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to load faigz index: {}", e)))?;
        
        let reader = FastaReader::new(&index)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to create faigz reader: {}", e)))?;
        
        for i in 0..index.num_sequences() {
            if let Some(name) = index.sequence_name(i) {
                match reader.fetch_seq_all(&name) {
                    Ok(seq_str) => {
                        sequences.push(Sequence {
                            id: name.clone(),
                            seq: seq_str.into_bytes(),
                        });
                    }
                    Err(e) => {
                        return Err(io::Error::new(io::ErrorKind::InvalidData, 
                            format!("Failed to read sequence {}: {}", name, e)));
                    }
                }
            }
        }
    } else {
        // Handle plain file using bio crate
        let file = File::open(&args.input)?;
        let fasta_reader = fasta::Reader::new(file);
        
        for result in fasta_reader.records() {
            let record = result?;
            sequences.push(Sequence {
                id: record.id().to_string(),
                seq: record.seq().to_vec(),
            });
        }
    }

    // Parse alignment parameters
    let params = parse_scores(&args.scores)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // Create the alignment iterator with sparsification
    let aligner = AllPairIterator::with_options(&sequences, params, true, sparsification);
    
    // Get actual number of pairs to process
    let total_pairs = aligner.pair_count();

    // Create progress bar (if enabled)
    let progress = if args.no_progress {
        ProgressBar::hidden()
    } else {
        let pb = ProgressBar::new(total_pairs as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({percent}%) {per_sec} ETA: {eta}")
                .unwrap()
                .progress_chars("=>-")
        );
        pb.set_message("Aligning sequences");
        pb
    };

    // Create atomic counter for progress updates
    let completed = Arc::new(AtomicUsize::new(0));
    let progress_clone = progress.clone();
    
    // Convert to parallel iterator and collect results
    let paf_records: Vec<String> = aligner
        .into_par_iter()
        .map(|alignment| {
            let result = alignment_to_paf(&alignment, &sequences);
            
            // Update progress atomically (only if progress bar is enabled)
            if !args.no_progress {
                let prev = completed.fetch_add(1, Ordering::Relaxed);
                if prev % 100 == 0 || prev == total_pairs - 1 {
                    progress_clone.set_position((prev + 1) as u64);
                }
            }
            
            result
        })
        .collect();

    // Ensure progress bar shows completion
    if !args.no_progress {
        progress.finish_with_message("Alignment complete!");
    }

    // Write results to output
    let mut output: Box<dyn Write> = if let Some(output_path) = args.output {
        Box::new(File::create(output_path)?)
    } else {
        Box::new(io::stdout())
    };

    for record in paf_records {
        writeln!(output, "{record}")?;
    }

    Ok(())
}
