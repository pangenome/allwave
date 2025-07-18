use allwave::{alignment_to_paf, parse_scores, AllPairIterator, Sequence, SparsificationStrategy};
use bio::io::fasta;
use clap::Parser;
use faigz_rs::{FastaFormat, FastaIndex, FastaReader};
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::sync::Arc;
use std::time::Instant;

/// Format duration in seconds to human-readable string
fn format_duration(seconds: f64) -> String {
    if seconds < 60.0 {
        format!("{seconds:.1}s")
    } else if seconds < 3600.0 {
        let minutes = (seconds / 60.0).floor();
        let remaining_seconds = seconds % 60.0;
        format!("{minutes:.0}m{remaining_seconds:.0}s")
    } else {
        let hours = (seconds / 3600.0).floor();
        let remaining_minutes = ((seconds % 3600.0) / 60.0).floor();
        format!("{hours:.0}h{remaining_minutes:.0}m")
    }
}

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

    /// Sparsification strategy:
    /// none | auto | random:<frac> | giant:<prob> | tree:<near>:<far>:<random>[:<kmer>]
    #[arg(short = 'p', long, default_value = "giant:0.99")]
    sparsification: String,

    /// Disable progress bar output
    #[arg(long)]
    no_progress: bool,

    /// Output mash distance matrix and exit (for debugging/comparison with mash)
    #[arg(long)]
    mash_matrix: bool,

    /// Use edlib edit distance for orientation detection instead of mash
    #[arg(long)]
    edlib_orientation: bool,
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
            if fraction <= 0.0 || fraction > 1.0 {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "Random fraction must be between 0 and 1"));
            }
            SparsificationStrategy::Random(fraction)
        }
        s if s.starts_with("giant:") => {
            let prob_str = &s[6..];
            let prob: f64 = prob_str.parse()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid giant component probability"))?;
            if prob <= 0.0 || prob >= 1.0 {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "Giant component probability must be between 0 and 1"));
            }
            SparsificationStrategy::Connectivity(prob)
        }
        // Legacy support for 'connectivity' - can be removed later
        s if s.starts_with("connectivity:") => {
            let prob_str = &s[13..];
            let prob: f64 = prob_str.parse()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid connectivity probability"))?;
            if prob <= 0.0 || prob >= 1.0 {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "Connectivity probability must be between 0 and 1"));
            }
            SparsificationStrategy::Connectivity(prob)
        }
        s if s.starts_with("tree:") => {
            let parts: Vec<&str> = s[5..].split(':').collect();
            if parts.len() < 3 || parts.len() > 4 {
                return Err(io::Error::new(io::ErrorKind::InvalidInput,
                    "Invalid tree format. Use: tree:<k_nearest>:<k_farthest>:<random_fraction>[:<kmer_size>]"));
            }

            let k_nearest: usize = parts[0].parse()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid k nearest count"))?;
            let k_farthest: usize = parts[1].parse()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid k farthest count"))?;
            let random_frac: f64 = parts[2].parse()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid random fraction"))?;

            if k_nearest == 0 && k_farthest == 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "At least one of k_nearest or k_farthest must be greater than 0"));
            }
            if !(0.0..=1.0).contains(&random_frac) {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "Random fraction must be between 0 and 1"));
            }

            let kmer_size = if parts.len() == 4 {
                let k: usize = parts[3].parse()
                    .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid k-mer size"))?;
                if !(3..=31).contains(&k) {
                    return Err(io::Error::new(io::ErrorKind::InvalidInput, "K-mer size must be between 3 and 31"));
                }
                Some(k)
            } else {
                None
            };

            SparsificationStrategy::TreeSampling(k_nearest, k_farthest, random_frac, kmer_size)
        }
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput,
            "Invalid sparsification strategy. Use: none, auto, giant:<probability>, random:<fraction>, or tree:<near>:<far>:<random>[:<kmer>]")),
    };

    // Read FASTA file (handle both plain and gzipped)
    let mut sequences = Vec::new();

    // Check if file is gzipped by extension and read accordingly
    if args.input.extension().and_then(|s| s.to_str()) == Some("gz") {
        // Use faigz-rs for gzipped files
        let index =
            FastaIndex::new(args.input.to_str().unwrap(), FastaFormat::Fasta).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Failed to load faigz index: {e}"),
                )
            })?;

        let reader = FastaReader::new(&index).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Failed to create faigz reader: {e}"),
            )
        })?;

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
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Failed to read sequence {name}: {e}"),
                        ));
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

    // Handle mash matrix output if requested
    if args.mash_matrix {
        // Extract k-mer size from sparsification parameter if using tree sampling
        let kmer_size = match &sparsification {
            SparsificationStrategy::TreeSampling(_, _, _, Some(k)) => *k,
            _ => 15, // Default k-mer size
        };

        // Compute and print mash distance matrix
        let matrix = allwave::mash::compute_distance_matrix_with_params(&sequences, kmer_size, 1000);
        allwave::mash::print_distance_matrix(&sequences, &matrix);
        return Ok(());
    }

    // Parse alignment parameters
    let params =
        parse_scores(&args.scores).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // Create the alignment iterator with sparsification
    let aligner = AllPairIterator::with_options(&sequences, params, true, !args.edlib_orientation, sparsification);

    // Get actual number of pairs to process
    let total_pairs = aligner.pair_count();

    // Determine if we're in interactive mode (output to terminal, not file)
    let is_interactive = args.output.is_none() && atty::is(atty::Stream::Stderr);

    // Create progress tracking
    let completed = Arc::new(AtomicUsize::new(0));
    let start_time = Instant::now();
    let last_log_time = Arc::new(std::sync::Mutex::new(Instant::now()));

    // Create progress bar for interactive mode only
    let progress = if args.no_progress {
        ProgressBar::hidden()
    } else if is_interactive {
        let pb = ProgressBar::new(total_pairs as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{elapsed_precise} {pos}/{len} ({percent}%) ETA: {eta}")
                .unwrap(),
        );
        pb
    } else {
        ProgressBar::hidden()
    };

    let progress_clone = progress.clone();
    let last_log_time_clone = last_log_time.clone();

    // Create channel for streaming PAF records
    let (tx, rx) = mpsc::channel::<String>();

    // Spawn writer thread
    let writer_handle = if let Some(output_path) = args.output {
        let output_path = output_path.clone();
        std::thread::spawn(move || -> io::Result<()> {
            let mut output = File::create(output_path)?;
            for paf_record in rx {
                writeln!(output, "{paf_record}")?;
            }
            Ok(())
        })
    } else {
        std::thread::spawn(move || -> io::Result<()> {
            let mut output = io::stdout();
            for paf_record in rx {
                writeln!(output, "{paf_record}")?;
            }
            Ok(())
        })
    };

    // Process alignments with streaming callback
    let result = aligner.for_each_with_callback(|alignment| {
        let paf_record = alignment_to_paf(&alignment, &sequences);

        // Send PAF record through channel
        tx.send(paf_record).map_err(|e| {
            Box::new(std::io::Error::new(
                std::io::ErrorKind::BrokenPipe,
                e.to_string(),
            )) as Box<dyn std::error::Error + Send + Sync>
        })?;

        // Update progress
        if !args.no_progress {
            let current = completed.fetch_add(1, Ordering::Relaxed) + 1;

            if is_interactive {
                // Interactive mode: update progress bar every alignment
                progress_clone.set_position(current as u64);
            } else {
                // File mode: log at most once per second
                let should_log = if current == total_pairs {
                    true // Always log completion
                } else {
                    // Check if at least 1 second has passed since last log
                    if let Ok(mut last_log) = last_log_time_clone.try_lock() {
                        let now = Instant::now();
                        if now.duration_since(*last_log).as_secs() >= 1 {
                            *last_log = now;
                            true
                        } else {
                            false
                        }
                    } else {
                        false
                    }
                };

                if should_log {
                    let elapsed = start_time.elapsed();
                    let percentage = (current as f64 / total_pairs as f64) * 100.0;
                    let rate = current as f64 / elapsed.as_secs_f64();

                    if current == total_pairs {
                        eprintln!(
                            "[{:.1}s] {}/{} ({:.1}%) {:.1} alignments/sec",
                            elapsed.as_secs_f64(),
                            current,
                            total_pairs,
                            percentage,
                            rate
                        );
                    } else {
                        // Calculate ETA
                        let remaining = total_pairs - current;
                        let eta_seconds = remaining as f64 / rate;
                        let eta_str = format_duration(eta_seconds);

                        eprintln!(
                            "[{:.1}s] {}/{} ({:.1}%) {:.1} alignments/sec ETA: {}",
                            elapsed.as_secs_f64(),
                            current,
                            total_pairs,
                            percentage,
                            rate,
                            eta_str
                        );
                    }
                }
            }
        }

        Ok(())
    });

    // Close the channel to signal writer thread to finish
    drop(tx);

    // Handle any errors from the streaming process
    result.map_err(io::Error::other)?;

    // Wait for writer thread to finish
    writer_handle
        .join()
        .map_err(|e| io::Error::other(format!("Writer thread panicked: {e:?}")))??;

    // Show completion message
    if !args.no_progress {
        let elapsed = start_time.elapsed();
        if is_interactive {
            progress.finish_with_message(format!(
                "Completed {} alignments in {:.1}s",
                total_pairs,
                elapsed.as_secs_f64()
            ));
        } else {
            let rate = total_pairs as f64 / elapsed.as_secs_f64();
            eprintln!(
                "[{:.1}s] {}/{} (100.0%) {:.1} alignments/sec - Complete!",
                elapsed.as_secs_f64(),
                total_pairs,
                total_pairs,
                rate
            );
        }
    }

    Ok(())
}
