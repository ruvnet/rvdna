//! Luby-Transform (LT) rateless fountain code — the outer erasure code.
//!
//! STUB — implemented by a swarm agent. See the contract below.

/// One LT encoded symbol: the XOR of a seed-determined subset of source blocks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Droplet {
    /// PRNG seed identifying which source blocks were XORed (the decoder
    /// reproduces the same neighbour set from this seed).
    pub seed: u32,
    /// Number of source blocks combined into this droplet.
    pub degree: u32,
    /// XOR of the selected blocks (length == `block_size`).
    pub data: Vec<u8>,
}

/// LT encoder over `num_blocks` source blocks of `block_size` bytes each.
#[derive(Debug, Clone)]
pub struct LtEncoder {
    pub num_blocks: usize,
    pub block_size: usize,
}

impl LtEncoder {
    pub fn new(num_blocks: usize, block_size: usize) -> Self {
        Self {
            num_blocks,
            block_size,
        }
    }

    /// Deterministically produce the droplet for `seed`.
    pub fn encode(&self, _blocks: &[Vec<u8>], _seed: u32) -> Droplet {
        unimplemented!("LtEncoder::encode")
    }
}

/// LT belief-propagation ("peeling") decoder.
#[derive(Debug, Clone)]
pub struct LtDecoder {
    pub num_blocks: usize,
    pub block_size: usize,
}

impl LtDecoder {
    pub fn new(num_blocks: usize, block_size: usize) -> Self {
        Self {
            num_blocks,
            block_size,
        }
    }

    /// Peel the droplets; return all source blocks once fully recovered.
    pub fn decode(&self, _droplets: &[Droplet]) -> Option<Vec<Vec<u8>>> {
        unimplemented!("LtDecoder::decode")
    }
}

/// Shared neighbour selection: given a `seed` and `num_blocks`, return the
/// droplet degree and the sorted, de-duplicated source-block indices it XORs.
/// MUST be identical on the encode and decode sides.
pub fn neighbours(_seed: u32, _num_blocks: usize) -> (u32, Vec<usize>) {
    unimplemented!("fountain::neighbours")
}
