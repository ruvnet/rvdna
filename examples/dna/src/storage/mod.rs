//! # DNA Storage Codec Simulator
//!
//! Turn an arbitrary file into a pool of synthetic DNA strands ("oligos"),
//! push them through a noisy sequencing channel (substitutions, insertions,
//! deletions, and whole-strand dropout), then reconstruct the original bytes —
//! demonstrating end-to-end error correction with no wet lab required.
//!
//! ## Pipeline
//!
//! ```text
//!  file bytes
//!     │  CRC32 + manifest
//!     ▼
//!  source blocks ──► Fountain (LT) outer code ──► droplets  (rateless erasure code)
//!     │                                              │
//!     │                                  per-strand header [idx|seed|degree]
//!     │                                              ▼
//!     │                              Reed–Solomon inner code (GF(256) parity)
//!     │                                              ▼
//!     │                       constraint-aware DNA mapping (homopolymer-free,
//!     │                                GC-balanced base-3 transform)
//!     ▼                                              ▼
//!  DnaArchive  ◄──────────────────────────────  Vec<Strand> = ACGT strings
//!
//!  ── noisy channel (sub / ins / del / dropout, configurable coverage) ──►
//!
//!  read pool ──► cluster + consensus ──► RS decode ──► Fountain peel ──►
//!            reassemble ──► CRC32 verify ──► original file bytes
//! ```
//!
//! The heavy lifting lives in focused submodules:
//! - [`constraints`] — reversible bytes⇄DNA mapping with GC / homopolymer constraints.
//! - [`gf256`] — Reed–Solomon (GF(256)) inner error-correction code.
//! - [`fountain`] — Luby-Transform rateless erasure (outer) code.
//! - [`channel`] — synthetic sequencing-error injection.
//! - [`consensus`] — read clustering + majority-vote consensus.

use crate::error::{DnaError, Result};
use serde::{Deserialize, Serialize};

pub mod channel;
pub mod constraints;
pub mod consensus;
pub mod fountain;
pub mod gf256;

// ============================================================================
// Shared configuration & data types (the swarm contract)
// ============================================================================

/// Tunable parameters for the storage codec.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EncodeParams {
    /// Number of payload bytes per source block (and per fountain droplet).
    pub block_size: usize,
    /// Reed–Solomon parity symbols appended to every strand (corrects `nsym/2`
    /// substitution errors per strand).
    pub rs_parity: usize,
    /// Redundancy factor: emit `ceil(overhead * num_blocks)` droplets/strands.
    /// Values > 1.0 over-provision so the fountain decoder tolerates lost strands.
    pub overhead: f64,
    /// Maximum allowed homopolymer run length (informational; the base-3
    /// transform guarantees runs of length 1).
    pub max_homopolymer: usize,
    /// Seed for deterministic droplet generation.
    pub seed: u64,
}

impl Default for EncodeParams {
    fn default() -> Self {
        Self {
            block_size: 32,
            rs_parity: 8,
            overhead: 1.8,
            max_homopolymer: 1,
            seed: 0xC0FFEE,
        }
    }
}

/// A single synthetic DNA strand (oligo).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Strand {
    /// Logical strand index (recoverable from the strand's own header).
    pub index: u32,
    /// The nucleotide sequence (`ACGT`).
    pub sequence: String,
}

/// Self-describing container: everything a decoder needs plus the strand pool.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaArchive {
    /// Original file name (best-effort, for the `decode` convenience path).
    pub filename: String,
    /// Original payload length in bytes.
    pub byte_len: usize,
    /// CRC32 of the original payload (integrity check after reassembly).
    pub crc32: u32,
    /// Number of source blocks the payload was split into.
    pub num_blocks: usize,
    /// Codec parameters used to produce this archive.
    pub params: EncodeParams,
    /// The synthetic DNA strands.
    pub strands: Vec<Strand>,
}

impl DnaArchive {
    /// Total nucleotides across all strands.
    pub fn total_bases(&self) -> usize {
        self.strands.iter().map(|s| s.sequence.len()).sum()
    }

    /// Information density in bits stored per nucleotide (payload bits / bases).
    pub fn bits_per_base(&self) -> f64 {
        let bases = self.total_bases();
        if bases == 0 {
            return 0.0;
        }
        (self.byte_len as f64 * 8.0) / bases as f64
    }

    /// Mean GC content across strands (0.0–1.0).
    pub fn mean_gc(&self) -> f64 {
        if self.strands.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .strands
            .iter()
            .map(|s| constraints::gc_content(&s.sequence))
            .sum();
        sum / self.strands.len() as f64
    }
}

/// Outcome of decoding a (possibly corrupted) read pool.
#[derive(Debug, Clone)]
pub struct DecodeReport {
    /// Recovered payload bytes (empty if reconstruction failed).
    pub bytes: Vec<u8>,
    /// Whether the recovered payload matched the archive CRC32.
    pub crc_ok: bool,
    /// How many distinct source blocks the fountain decoder recovered.
    pub blocks_recovered: usize,
    /// How many reads survived RS decoding into valid strands.
    pub strands_recovered: usize,
    /// Number of reads fed into the decoder.
    pub reads_in: usize,
}

// ============================================================================
// CRC32 (IEEE) — small dependency-free integrity check
// ============================================================================

/// Compute the IEEE CRC32 of `data`.
pub fn crc32(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &b in data {
        crc ^= b as u32;
        for _ in 0..8 {
            let mask = (crc & 1).wrapping_neg();
            crc = (crc >> 1) ^ (0xEDB8_8320 & mask);
        }
    }
    !crc
}

// ============================================================================
// Orchestrator (filled in during integration)
// ============================================================================

/// End-to-end DNA storage codec.
#[derive(Debug, Clone)]
pub struct DnaStorageCodec {
    params: EncodeParams,
}

impl DnaStorageCodec {
    /// Create a codec with the given parameters.
    pub fn new(params: EncodeParams) -> Self {
        Self { params }
    }

    /// Encode raw bytes into a DNA archive.
    pub fn encode(&self, _filename: &str, _data: &[u8]) -> Result<DnaArchive> {
        Err(DnaError::Storage("orchestrator not yet wired".into()))
    }

    /// Decode a pool of (noisy) reads back into the original bytes.
    pub fn decode(&self, _archive: &DnaArchive, _reads: &[String]) -> Result<DecodeReport> {
        Err(DnaError::Storage("orchestrator not yet wired".into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc32_known_vector() {
        // CRC32("123456789") == 0xCBF43926
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
    }
}
