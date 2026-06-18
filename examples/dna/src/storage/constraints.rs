//! Reversible bytes ⇄ DNA mapping under biological constraints.
//!
//! STUB — implemented by a swarm agent. See the function contracts below.

/// Encode bytes into a homopolymer-free, GC-balanced `ACGT` string.
pub fn encode_bytes(_data: &[u8]) -> String {
    unimplemented!("constraints::encode_bytes")
}

/// Inverse of [`encode_bytes`]. Returns the original bytes.
pub fn decode_bytes(_seq: &str) -> Result<Vec<u8>, String> {
    unimplemented!("constraints::decode_bytes")
}

/// GC content fraction (0.0–1.0) of a sequence.
pub fn gc_content(_seq: &str) -> f64 {
    unimplemented!("constraints::gc_content")
}

/// Longest run of an identical nucleotide.
pub fn max_homopolymer_run(_seq: &str) -> usize {
    unimplemented!("constraints::max_homopolymer_run")
}
