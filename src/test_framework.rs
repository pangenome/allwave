use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;

/// Represents a mutation that was applied to a sequence
#[derive(Debug, Clone)]
pub struct Mutation {
    pub mutation_type: MutationType,
    pub position: usize,
    pub original: String,
    pub mutated: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MutationType {
    SNP,
    Insertion,
    Deletion,
    MicrosatelliteExpansion,
    MicrosatelliteContraction,
    CNVDuplication,
    CNVDeletion,
}

/// Configuration for sequence generation
pub struct SequenceConfig {
    pub length: usize,
    pub gc_content: f64,
    pub seed: u64,
}

impl Default for SequenceConfig {
    fn default() -> Self {
        SequenceConfig {
            length: 10000,
            gc_content: 0.5,
            seed: 42,
        }
    }
}

/// Configuration for mutations
pub struct MutationConfig {
    pub snp_rate: f64,
    pub indel_rate: f64,
    pub microsatellite_rate: f64,
    pub cnv_rate: f64,
    pub cnv_min_size: usize,
    pub cnv_max_size: usize,
    pub microsatellite_min_unit: usize,
    pub microsatellite_max_unit: usize,
    pub microsatellite_min_repeats: usize,
    pub microsatellite_max_repeats: usize,
    pub indel_max_size: usize,
    pub seed: u64,
}

impl Default for MutationConfig {
    fn default() -> Self {
        MutationConfig {
            snp_rate: 0.001,
            indel_rate: 0.0001,
            microsatellite_rate: 0.00005,
            cnv_rate: 0.00001,
            cnv_min_size: 1000,
            cnv_max_size: 5000,
            microsatellite_min_unit: 1,
            microsatellite_max_unit: 6,
            microsatellite_min_repeats: 5,
            microsatellite_max_repeats: 20,
            indel_max_size: 10,
            seed: 42,
        }
    }
}

/// Generate a random DNA sequence
pub fn generate_random_sequence(config: &SequenceConfig) -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut sequence = Vec::with_capacity(config.length);

    let gc_bases = ['G', 'C'];
    let at_bases = ['A', 'T'];

    for _ in 0..config.length {
        let base = if rng.gen::<f64>() < config.gc_content {
            gc_bases[rng.gen_range(0..2)]
        } else {
            at_bases[rng.gen_range(0..2)]
        };
        sequence.push(base as u8);
    }

    sequence
}

/// Find microsatellite regions in a sequence
fn find_microsatellites(
    sequence: &[u8],
    min_unit: usize,
    max_unit: usize,
    min_repeats: usize,
) -> Vec<(usize, usize, Vec<u8>)> {
    let mut microsatellites = Vec::new();

    for unit_size in min_unit..=max_unit {
        let mut i = 0;
        while i + unit_size * min_repeats <= sequence.len() {
            let unit = &sequence[i..i + unit_size];
            let mut repeats = 1;
            let mut j = i + unit_size;

            while j + unit_size <= sequence.len() && &sequence[j..j + unit_size] == unit {
                repeats += 1;
                j += unit_size;
            }

            if repeats >= min_repeats {
                microsatellites.push((i, repeats * unit_size, unit.to_vec()));
                i = j;
            } else {
                i += 1;
            }
        }
    }

    microsatellites
}

/// Apply mutations to a sequence
pub fn apply_mutations(sequence: &[u8], config: &MutationConfig) -> (Vec<u8>, Vec<Mutation>) {
    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut mutated = sequence.to_vec();
    let mut mutations = Vec::new();
    let mut offset: i32 = 0; // Track position changes due to indels

    // Apply SNPs
    for (i, &base) in sequence.iter().enumerate() {
        if rng.gen::<f64>() < config.snp_rate {
            let original_base = base as char;
            let bases = ['A', 'T', 'G', 'C'];
            let mut new_base = bases[rng.gen_range(0..4)];
            while new_base == original_base {
                new_base = bases[rng.gen_range(0..4)];
            }

            let adjusted_pos = (i as i32 + offset) as usize;
            if adjusted_pos < mutated.len() {
                mutated[adjusted_pos] = new_base as u8;
                mutations.push(Mutation {
                    mutation_type: MutationType::SNP,
                    position: i,
                    original: original_base.to_string(),
                    mutated: new_base.to_string(),
                });
            }
        }
    }

    // Apply indels
    let mut i = 0;
    while i < sequence.len() {
        if rng.gen::<f64>() < config.indel_rate {
            let is_insertion = rng.gen::<bool>();
            let size = rng.gen_range(1..=config.indel_max_size);
            let adjusted_pos = (i as i32 + offset) as usize;

            if is_insertion {
                // Generate random insertion
                let mut insertion = Vec::with_capacity(size);
                for _ in 0..size {
                    let bases = ['A', 'T', 'G', 'C'];
                    insertion.push(bases[rng.gen_range(0..4)] as u8);
                }

                if adjusted_pos <= mutated.len() {
                    mutated.splice(adjusted_pos..adjusted_pos, insertion.clone());
                    mutations.push(Mutation {
                        mutation_type: MutationType::Insertion,
                        position: i,
                        original: String::new(),
                        mutated: String::from_utf8(insertion).unwrap(),
                    });
                    offset += size as i32;
                }
            } else {
                // Deletion
                let end = std::cmp::min(adjusted_pos + size, mutated.len());
                if adjusted_pos < mutated.len() && end > adjusted_pos {
                    let deleted = mutated[adjusted_pos..end].to_vec();
                    mutated.splice(adjusted_pos..end, std::iter::empty());
                    mutations.push(Mutation {
                        mutation_type: MutationType::Deletion,
                        position: i,
                        original: String::from_utf8(deleted).unwrap(),
                        mutated: String::new(),
                    });
                    offset -= (end - adjusted_pos) as i32;
                }
            }
            i += size; // Skip ahead to avoid overlapping mutations
        }
        i += 1;
    }

    // Apply microsatellite mutations
    let microsatellites = find_microsatellites(
        sequence,
        config.microsatellite_min_unit,
        config.microsatellite_max_unit,
        config.microsatellite_min_repeats,
    );

    for (start, length, unit) in microsatellites {
        if rng.gen::<f64>() < config.microsatellite_rate {
            let current_repeats = length / unit.len();
            let is_expansion = rng.gen::<bool>();

            let new_repeats = if is_expansion {
                rng.gen_range(current_repeats + 1..=config.microsatellite_max_repeats)
            } else {
                rng.gen_range(2..current_repeats)
            };

            let adjusted_start = (start as i32 + offset) as usize;
            let adjusted_end = (start as i32 + length as i32 + offset) as usize;

            if adjusted_start < mutated.len() && adjusted_end <= mutated.len() {
                let mut new_microsatellite = Vec::new();
                for _ in 0..new_repeats {
                    new_microsatellite.extend_from_slice(&unit);
                }

                let original_str =
                    String::from_utf8(mutated[adjusted_start..adjusted_end].to_vec()).unwrap();
                let mutated_str = String::from_utf8(new_microsatellite.clone()).unwrap();

                mutated.splice(adjusted_start..adjusted_end, new_microsatellite);

                mutations.push(Mutation {
                    mutation_type: if is_expansion {
                        MutationType::MicrosatelliteExpansion
                    } else {
                        MutationType::MicrosatelliteContraction
                    },
                    position: start,
                    original: original_str,
                    mutated: mutated_str,
                });

                offset += (new_repeats * unit.len()) as i32 - length as i32;
            }
        }
    }

    // Apply CNVs
    let mut i = 0;
    while i + config.cnv_min_size < sequence.len() {
        if rng.gen::<f64>() < config.cnv_rate {
            let cnv_size = rng.gen_range(config.cnv_min_size..=config.cnv_max_size);
            let cnv_size = std::cmp::min(cnv_size, sequence.len() - i);
            let is_duplication = rng.gen::<bool>();

            let adjusted_pos = (i as i32 + offset) as usize;

            if is_duplication {
                // CNV duplication - duplicate 2-5 times
                let copies = rng.gen_range(2..=5);
                if adjusted_pos + cnv_size <= mutated.len() {
                    let segment = mutated[adjusted_pos..adjusted_pos + cnv_size].to_vec();
                    let original_str = String::from_utf8_lossy(&segment).to_string();

                    let mut duplicated = Vec::new();
                    for _ in 0..copies {
                        duplicated.extend_from_slice(&segment);
                    }
                    let mutated_str = String::from_utf8_lossy(&duplicated).to_string();

                    // Insert additional copies after the original
                    for _ in 1..copies {
                        mutated.splice(
                            adjusted_pos + cnv_size..adjusted_pos + cnv_size,
                            segment.clone(),
                        );
                        offset += cnv_size as i32;
                    }

                    mutations.push(Mutation {
                        mutation_type: MutationType::CNVDuplication,
                        position: i,
                        original: original_str,
                        mutated: mutated_str,
                    });
                }
            } else {
                // CNV deletion
                let end = std::cmp::min(adjusted_pos + cnv_size, mutated.len());
                if adjusted_pos < mutated.len() && end > adjusted_pos {
                    let deleted = mutated[adjusted_pos..end].to_vec();
                    mutated.splice(adjusted_pos..end, std::iter::empty());

                    mutations.push(Mutation {
                        mutation_type: MutationType::CNVDeletion,
                        position: i,
                        original: String::from_utf8_lossy(&deleted).to_string(),
                        mutated: String::new(),
                    });
                    offset -= (end - adjusted_pos) as i32;
                }
            }
            i += cnv_size; // Skip ahead to avoid overlapping CNVs
        }
        i += 1;
    }

    (mutated, mutations)
}

/// Create a test case with known mutations
pub fn create_test_case(
    name: &str,
    seq_config: &SequenceConfig,
    mut_config: &MutationConfig,
) -> TestCase {
    let reference = generate_random_sequence(seq_config);
    let (query, mutations) = apply_mutations(&reference, mut_config);

    TestCase {
        name: name.to_string(),
        reference,
        query,
        mutations,
    }
}

#[derive(Debug)]
pub struct TestCase {
    pub name: String,
    pub reference: Vec<u8>,
    pub query: Vec<u8>,
    pub mutations: Vec<Mutation>,
}

impl TestCase {
    /// Get mutation statistics
    pub fn get_mutation_stats(&self) -> HashMap<MutationType, usize> {
        let mut stats = HashMap::new();
        for mutation in &self.mutations {
            *stats.entry(mutation.mutation_type.clone()).or_insert(0) += 1;
        }
        stats
    }

    /// Write test case to FASTA file
    pub fn write_fasta(&self, path: &str) -> std::io::Result<()> {
        use std::io::Write;
        let mut file = std::fs::File::create(path)?;

        writeln!(file, ">{}_reference", self.name)?;
        writeln!(file, "{}", String::from_utf8_lossy(&self.reference))?;
        writeln!(file, ">{}_query", self.name)?;
        writeln!(file, "{}", String::from_utf8_lossy(&self.query))?;

        Ok(())
    }
}
