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
pub mod consensus;
pub mod constraints;
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
            overhead: 2.5,
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
// Orchestrator
// ============================================================================

use channel::ErrorModel;
use fountain::{Droplet, LtDecoder, LtEncoder};
use gf256::ReedSolomon;

/// Per-strand header: `index` (u32 BE) + fountain `seed` (u32 BE).
const HEADER_LEN: usize = 8;

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

    /// Codec parameters.
    pub fn params(&self) -> &EncodeParams {
        &self.params
    }

    /// Split `data` into `num_blocks` source blocks of `block_size` bytes,
    /// zero-padding the final block.
    fn split_blocks(&self, data: &[u8]) -> (usize, Vec<Vec<u8>>) {
        let bs = self.params.block_size.max(1);
        let num_blocks = (data.len() + bs - 1) / bs;
        let num_blocks = num_blocks.max(1);
        let mut blocks = Vec::with_capacity(num_blocks);
        for i in 0..num_blocks {
            let start = i * bs;
            let end = (start + bs).min(data.len());
            let mut block = vec![0u8; bs];
            if start < data.len() {
                block[..end - start].copy_from_slice(&data[start..end]);
            }
            blocks.push(block);
        }
        (num_blocks, blocks)
    }

    /// Encode raw bytes into a DNA archive.
    pub fn encode(&self, filename: &str, data: &[u8]) -> Result<DnaArchive> {
        let bs = self.params.block_size.max(1);
        let codeword_len = HEADER_LEN + bs + self.params.rs_parity;
        if codeword_len > 255 {
            return Err(DnaError::Storage(format!(
                "strand codeword {codeword_len} bytes exceeds GF(256) limit of 255; \
                 reduce block_size or rs_parity"
            )));
        }

        let (num_blocks, blocks) = self.split_blocks(data);
        let encoder = LtEncoder::new(num_blocks, bs);
        let rs = ReedSolomon::new(self.params.rs_parity);
        let base_seed = self.params.seed as u32;

        let num_strands = ((self.params.overhead.max(1.0)) * num_blocks as f64).ceil() as usize;
        // Add a small additive headroom so tiny files (small K), where the
        // robust-soliton peeling is least reliable, still over-provision enough
        // droplets to decode.
        let num_strands = num_strands.max(num_blocks + 8);

        let mut strands = Vec::with_capacity(num_strands);
        for i in 0..num_strands {
            let seed = base_seed.wrapping_add(i as u32);
            let droplet = encoder.encode(&blocks, seed);

            // Strand payload bytes: [index(4) | seed(4) | droplet.data]
            let mut strand_bytes = Vec::with_capacity(codeword_len - self.params.rs_parity);
            strand_bytes.extend_from_slice(&(i as u32).to_be_bytes());
            strand_bytes.extend_from_slice(&seed.to_be_bytes());
            strand_bytes.extend_from_slice(&droplet.data);

            let codeword = rs.encode(&strand_bytes);
            let sequence = constraints::encode_bytes(&codeword);
            strands.push(Strand {
                index: i as u32,
                sequence,
            });
        }

        Ok(DnaArchive {
            filename: filename.to_string(),
            byte_len: data.len(),
            crc32: crc32(data),
            num_blocks,
            params: self.params.clone(),
            strands,
        })
    }

    /// Decode a pool of (noisy) reads back into the original bytes.
    ///
    /// `archive` supplies the manifest (block count, sizes, CRC); only its
    /// `strands` are ignored — decoding works purely from `reads`.
    pub fn decode(&self, archive: &DnaArchive, reads: &[String]) -> Result<DecodeReport> {
        let bs = archive.params.block_size.max(1);
        let nsym = archive.params.rs_parity;
        let codeword_len = HEADER_LEN + bs + nsym;
        let expected_bases = codeword_len * 6;
        let rs = ReedSolomon::new(nsym);

        let mut droplets: Vec<Droplet> = Vec::new();
        let mut seen_seeds: std::collections::HashSet<u32> = std::collections::HashSet::new();

        // Demap a single read/consensus string into a strand codeword and
        // RS-decode it into a droplet, recording it (deduped by fountain seed).
        // Returns true if a *new* droplet was recovered.
        let mut try_recover = |seq: &str,
                               droplets: &mut Vec<Droplet>,
                               seen: &mut std::collections::HashSet<u32>|
         -> bool {
            let bases: Vec<char> = seq.chars().collect();
            // Normalise length to the expected codeword (indels shift length).
            let usable = if bases.len() >= expected_bases {
                expected_bases
            } else {
                (bases.len() / 6) * 6
            };
            if usable == 0 {
                return false;
            }
            let trimmed: String = bases[..usable].iter().collect();
            let bytes = match constraints::decode_bytes(&trimmed) {
                Ok(b) => b,
                Err(_) => return false,
            };
            if bytes.len() < codeword_len {
                return false;
            }
            let msg = match rs.decode(&bytes[..codeword_len]) {
                Ok(m) => m,
                Err(_) => return false,
            };
            if msg.len() < HEADER_LEN {
                return false;
            }
            let seed = u32::from_be_bytes([msg[4], msg[5], msg[6], msg[7]]);
            if !seen.insert(seed) {
                return false; // already have this strand
            }
            let mut payload = msg[HEADER_LEN..].to_vec();
            payload.resize(bs, 0);
            let (degree, _idx) = fountain::neighbours(seed, archive.num_blocks);
            droplets.push(Droplet {
                seed,
                degree,
                data: payload,
            });
            true
        };

        // Pass 1 — decode every read independently. RS recovers each strand's
        // exact [index|seed] header, so identity comes from the code, not from
        // fuzzy sequence matching. Reads RS can't fix are deferred to pass 2.
        let mut residual: Vec<usize> = Vec::new();
        for (i, read) in reads.iter().enumerate() {
            if !try_recover(read, &mut droplets, &mut seen_seeds) {
                residual.push(i);
            }
        }

        // Pass 2 (salvage) — when coverage > 1, the same strand was sequenced
        // several times. Cluster the *leftover* noisy reads and majority-vote a
        // consensus per cluster, which cancels random substitutions, then retry
        // RS. This only ever ADDS strands, so it never harms the clean path.
        if residual.len() > 1 {
            let residual_reads: Vec<String> = residual.iter().map(|&i| reads[i].clone()).collect();
            for group in consensus::cluster(&residual_reads, 0.3) {
                if group.len() < 2 {
                    continue;
                }
                let members: Vec<&str> =
                    group.iter().map(|&g| residual_reads[g].as_str()).collect();
                let cons = consensus::consensus(&members);
                try_recover(&cons, &mut droplets, &mut seen_seeds);
            }
        }

        let strands_recovered = droplets.len();

        // 3. Fountain peeling to recover the source blocks.
        let decoder = LtDecoder::new(archive.num_blocks, bs);
        let (bytes, blocks_recovered, decoded) = match decoder.decode(&droplets) {
            Some(blocks) => {
                let mut out = Vec::with_capacity(blocks.len() * bs);
                for b in &blocks {
                    out.extend_from_slice(b);
                }
                out.truncate(archive.byte_len);
                (out, archive.num_blocks, true)
            }
            None => (Vec::new(), 0, false),
        };

        // `decoded` (not `!bytes.is_empty()`) gates success so a zero-byte file
        // still verifies: an empty payload that peels cleanly must report crc_ok.
        let crc_ok = decoded && crc32(&bytes) == archive.crc32;

        Ok(DecodeReport {
            bytes,
            crc_ok,
            blocks_recovered,
            strands_recovered,
            reads_in: reads.len(),
        })
    }

    /// Convenience: encode `data`, push it through `model`, then decode —
    /// the full simulator round-trip used by the demo and CLI.
    pub fn simulate(
        &self,
        filename: &str,
        data: &[u8],
        model: &ErrorModel,
        channel_seed: u64,
    ) -> Result<(DnaArchive, DecodeReport)> {
        let archive = self.encode(filename, data)?;
        let sequences: Vec<String> = archive.strands.iter().map(|s| s.sequence.clone()).collect();
        let reads = channel::apply(&sequences, model, channel_seed);
        let report = self.decode(&archive, &reads)?;
        Ok((archive, report))
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

    #[test]
    fn roundtrip_no_errors() {
        let codec = DnaStorageCodec::new(EncodeParams::default());
        let data: Vec<u8> = (0..512).map(|i| (i * 31 + 7) as u8).collect();
        let archive = codec.encode("blob.bin", &data).unwrap();
        assert!(archive.mean_gc() > 0.3 && archive.mean_gc() < 0.7);
        // Feed the pristine strands straight back in.
        let reads: Vec<String> = archive.strands.iter().map(|s| s.sequence.clone()).collect();
        let report = codec.decode(&archive, &reads).unwrap();
        assert!(report.crc_ok, "clean round-trip must verify");
        assert_eq!(report.bytes, data);
    }

    #[test]
    fn roundtrip_with_substitutions_and_dropout() {
        let params = EncodeParams {
            block_size: 24,
            rs_parity: 12,
            overhead: 2.5,
            max_homopolymer: 1,
            seed: 42,
        };
        let codec = DnaStorageCodec::new(params);
        let data: Vec<u8> = (0..400).map(|i| (i * 17 + 3) as u8).collect();
        let model = ErrorModel {
            p_sub: 0.01,
            p_ins: 0.0,
            p_del: 0.0,
            p_drop: 0.15,
            coverage: 3,
        };
        let (_archive, report) = codec.simulate("photo.png", &data, &model, 7).unwrap();
        assert!(
            report.crc_ok,
            "recovery failed: blocks {}/{}, strands {}",
            report.blocks_recovered, _archive.num_blocks, report.strands_recovered
        );
        assert_eq!(report.bytes, data);
    }

    #[test]
    fn roundtrip_empty_and_tiny_files() {
        let codec = DnaStorageCodec::new(EncodeParams::default());
        for data in [vec![], vec![0x42u8], vec![1u8, 2, 3, 4, 5]] {
            let archive = codec.encode("tiny", &data).unwrap();
            let reads: Vec<String> = archive.strands.iter().map(|s| s.sequence.clone()).collect();
            let report = codec.decode(&archive, &reads).unwrap();
            assert!(report.crc_ok, "len {} must verify", data.len());
            assert_eq!(report.bytes, data, "len {} mismatch", data.len());
        }
    }

    #[test]
    fn recovers_under_indels_with_coverage() {
        // Insertions/deletions frame-shift a strand; high coverage + the
        // consensus salvage pass + fountain over-provisioning should still
        // reconstruct the payload.
        let params = EncodeParams {
            block_size: 24,
            rs_parity: 12,
            overhead: 3.0,
            max_homopolymer: 1,
            seed: 11,
        };
        let codec = DnaStorageCodec::new(params);
        let data: Vec<u8> = (0..300).map(|i| (i * 13 + 5) as u8).collect();
        let model = ErrorModel {
            p_sub: 0.005,
            p_ins: 0.002,
            p_del: 0.002,
            p_drop: 0.05,
            coverage: 6,
        };
        let (archive, report) = codec.simulate("indel.bin", &data, &model, 3).unwrap();
        assert!(
            report.crc_ok,
            "indel recovery failed: blocks {}/{}, strands {}",
            report.blocks_recovered, archive.num_blocks, report.strands_recovered
        );
        assert_eq!(report.bytes, data);
    }
}
