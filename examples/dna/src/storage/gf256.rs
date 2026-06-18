//! Reed–Solomon inner code over GF(256).
//!
//! STUB — implemented by a swarm agent. See the contract below.

/// Systematic Reed–Solomon codec over GF(256) with `nsym` parity symbols.
#[derive(Debug, Clone)]
pub struct ReedSolomon {
    /// Number of parity symbols appended by [`ReedSolomon::encode`].
    pub nsym: usize,
}

impl ReedSolomon {
    /// New codec appending `nsym` parity bytes (corrects up to `nsym / 2` errors).
    pub fn new(nsym: usize) -> Self {
        Self { nsym }
    }

    /// Append `nsym` parity bytes; returns `data.len() + nsym` bytes.
    /// `data.len() + nsym` must be ≤ 255.
    pub fn encode(&self, _data: &[u8]) -> Vec<u8> {
        unimplemented!("ReedSolomon::encode")
    }

    /// Correct up to `nsym / 2` substitution errors and return the message
    /// (parity stripped), or `Err` if the codeword is unrecoverable.
    pub fn decode(&self, _received: &[u8]) -> Result<Vec<u8>, String> {
        unimplemented!("ReedSolomon::decode")
    }
}
