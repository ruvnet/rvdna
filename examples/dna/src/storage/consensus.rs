//! Read clustering + majority-vote consensus.
//!
//! STUB — implemented by a swarm agent. See the contract below.

/// Greedily cluster reads by sequence similarity (normalised edit/Hamming
/// distance ≤ `max_dist`). Returns groups of read indices into `reads`.
pub fn cluster(_reads: &[String], _max_dist: f64) -> Vec<Vec<usize>> {
    unimplemented!("consensus::cluster")
}

/// Position-wise majority-vote consensus over a group of reads. Handles small
/// length differences (caused by indels) by voting up to the modal length.
pub fn consensus(_reads: &[&str]) -> String {
    unimplemented!("consensus::consensus")
}
