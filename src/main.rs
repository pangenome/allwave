use allwave::{AllPairIterator, Sequence, alignment_to_paf, parse_scores};
use bio::io::fasta;
use clap::Parser;
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;

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





fn main() -> io::Result<()> {
    let args = Args::parse();
    
    // Set the number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .map_err(io::Error::other)?;

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
            sequences.push(Sequence {
                id: record.id().to_string(),
                seq: record.seq().to_vec(),
            });
        }
    } else {
        // Handle plain file
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

    // Create the alignment iterator
    let aligner = AllPairIterator::new(&sequences, params);
    
    // Convert to parallel iterator and collect results
    let paf_records: Vec<String> = aligner
        .into_par_iter()
        .map(|alignment| alignment_to_paf(&alignment, &sequences))
        .collect();

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
