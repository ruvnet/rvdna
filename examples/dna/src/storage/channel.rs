//! Synthetic sequencing-error channel.
//!
//! Models the noise of a real sequencing run: per-base substitutions,
//! insertions and deletions, whole-strand dropout, and configurable read
//! coverage. All randomness is driven by a single seeded RNG so a given
//! `(strands, model, seed)` always produces identical output.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

/// Per-base error rates plus whole-strand dropout and sequencing coverage.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorModel {
    /// Per-base substitution probability.
    pub p_sub: f64,
    /// Per-base insertion probability.
    pub p_ins: f64,
    /// Per-base deletion probability.
    pub p_del: f64,
    /// Probability that an entire strand is lost (never sequenced).
    pub p_drop: f64,
    /// Reads produced per surviving strand (sequencing depth).
    pub coverage: usize,
}

impl Default for ErrorModel {
    fn default() -> Self {
        Self {
            p_sub: 0.01,
            p_ins: 0.002,
            p_del: 0.002,
            p_drop: 0.05,
            coverage: 1,
        }
    }
}

const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// Pick a uniformly random base.
fn random_base(rng: &mut StdRng) -> char {
    BASES[rng.gen_range(0..4)]
}

/// Pick a uniformly random base different from `exclude`.
fn random_base_excluding(rng: &mut StdRng, exclude: char) -> char {
    loop {
        let b = random_base(rng);
        if b != exclude {
            return b;
        }
    }
}

/// Push every strand through the channel, returning the read pool.
///
/// `seed` makes the corruption deterministic. Dropped strands contribute no
/// reads; surviving strands contribute up to `coverage` (independently noisy)
/// reads each.
pub fn apply(strands: &[String], model: &ErrorModel, seed: u64) -> Vec<String> {
    let mut rng = StdRng::seed_from_u64(seed);
    let coverage = if model.coverage == 0 {
        1
    } else {
        model.coverage
    };

    let mut reads: Vec<String> = Vec::new();

    for strand in strands {
        // Whole-strand dropout.
        if rng.gen::<f64>() < model.p_drop {
            continue;
        }

        for _ in 0..coverage {
            let mut read = String::with_capacity(strand.len());
            for base in strand.chars() {
                // Deletion: skip this base entirely.
                if rng.gen::<f64>() < model.p_del {
                    continue;
                }
                // Insertion: emit a random base BEFORE copying.
                if rng.gen::<f64>() < model.p_ins {
                    read.push(random_base(&mut rng));
                }
                // Substitution vs. faithful copy.
                if rng.gen::<f64>() < model.p_sub {
                    read.push(random_base_excluding(&mut rng, base));
                } else {
                    read.push(base);
                }
            }
            reads.push(read);
        }
    }

    reads
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_strand(rng: &mut StdRng, len: usize) -> String {
        (0..len).map(|_| random_base(rng)).collect()
    }

    #[test]
    fn substitution_rate_in_tolerance() {
        let model = ErrorModel {
            p_sub: 0.05,
            p_ins: 0.0,
            p_del: 0.0,
            p_drop: 0.0,
            coverage: 1,
        };

        // Many long strands so reads stay aligned (no ins/del) and we can
        // measure per-base substitution rate directly.
        let mut gen = StdRng::seed_from_u64(42);
        let strands: Vec<String> = (0..200).map(|_| random_strand(&mut gen, 200)).collect();

        let reads = apply(&strands, &model, 7);
        assert_eq!(reads.len(), strands.len());

        let mut total = 0usize;
        let mut subs = 0usize;
        for (orig, read) in strands.iter().zip(reads.iter()) {
            assert_eq!(orig.len(), read.len(), "no indels expected");
            for (a, b) in orig.chars().zip(read.chars()) {
                total += 1;
                if a != b {
                    subs += 1;
                }
            }
        }
        let rate = subs as f64 / total as f64;
        assert!(
            (rate - 0.05).abs() < 0.01,
            "measured sub rate {rate} not near 0.05"
        );
    }

    #[test]
    fn dropout_and_coverage() {
        let mut gen = StdRng::seed_from_u64(1);
        let strands: Vec<String> = (0..100).map(|_| random_strand(&mut gen, 50)).collect();

        // Near-certain dropout => (near) zero reads.
        let drop_model = ErrorModel {
            p_sub: 0.0,
            p_ins: 0.0,
            p_del: 0.0,
            p_drop: 1.0,
            coverage: 3,
        };
        let dropped = apply(&strands, &drop_model, 99);
        assert!(dropped.is_empty(), "p_drop=1.0 should yield no reads");

        // No dropout, coverage 3 => exactly 3x reads.
        let cov_model = ErrorModel {
            p_sub: 0.0,
            p_ins: 0.0,
            p_del: 0.0,
            p_drop: 0.0,
            coverage: 3,
        };
        let covered = apply(&strands, &cov_model, 99);
        assert_eq!(covered.len(), strands.len() * 3);
    }

    #[test]
    fn determinism() {
        let mut gen = StdRng::seed_from_u64(5);
        let strands: Vec<String> = (0..50).map(|_| random_strand(&mut gen, 80)).collect();
        let model = ErrorModel::default();

        let a = apply(&strands, &model, 12345);
        let b = apply(&strands, &model, 12345);
        assert_eq!(a, b, "same seed must give identical output");
    }

    #[test]
    fn only_valid_bases_emitted() {
        let mut gen = StdRng::seed_from_u64(3);
        let strands: Vec<String> = (0..20).map(|_| random_strand(&mut gen, 60)).collect();
        let model = ErrorModel {
            p_sub: 0.1,
            p_ins: 0.1,
            p_del: 0.1,
            p_drop: 0.0,
            coverage: 2,
        };
        let reads = apply(&strands, &model, 8);
        for read in &reads {
            assert!(read.chars().all(|c| BASES.contains(&c)));
        }
    }
}
