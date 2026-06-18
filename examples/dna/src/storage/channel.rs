//! Synthetic sequencing-error channel.
//!
//! STUB — implemented by a swarm agent. See the contract below.

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

/// Push every strand through the channel, returning the shuffled read pool.
///
/// `seed` makes the corruption deterministic. Dropped strands contribute no
/// reads; surviving strands contribute up to `coverage` (independently noisy)
/// reads each.
pub fn apply(_strands: &[String], _model: &ErrorModel, _seed: u64) -> Vec<String> {
    unimplemented!("channel::apply")
}
