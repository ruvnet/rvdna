//! Luby-Transform (LT) rateless fountain code — the outer erasure code.
//!
//! The encoder produces an unbounded stream of [`Droplet`]s, each the XOR of a
//! seed-determined subset of source blocks. The [`LtDecoder`] reconstructs the
//! original blocks via belief-propagation ("peeling") once enough droplets are
//! collected. Encode and decode agree on the per-droplet neighbour set because
//! both derive it deterministically from the droplet's `seed` via
//! [`neighbours`], which uses an inline SplitMix64 PRNG (portable, no reliance
//! on `rand`'s generators).

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

/// Deterministic SplitMix64 PRNG. Implemented inline so neighbour selection is
/// fully portable and identical on encode and decode sides.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform f64 in [0, 1).
    fn next_f64(&mut self) -> f64 {
        // Use the top 53 bits for a uniform double.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform integer in `0..n` (n > 0).
    fn next_below(&mut self, n: u64) -> u64 {
        // Simple modulo; bias is negligible for the small ranges used here.
        self.next_u64() % n
    }
}

/// Build the normalized Robust Soliton PMF over degrees `1..=k`.
///
/// Returns a vector of length `k` where index `i` is the probability of degree
/// `i + 1`. Depends only on `k`.
fn robust_soliton_pmf(k: usize) -> Vec<f64> {
    if k == 0 {
        return Vec::new();
    }
    if k == 1 {
        return vec![1.0];
    }

    let c = 0.03_f64;
    let delta = 0.05_f64;
    let kf = k as f64;

    // Ideal soliton: rho(1) = 1/K, rho(d) = 1/(d*(d-1)) for d = 2..=K.
    let mut rho = vec![0.0_f64; k];
    rho[0] = 1.0 / kf;
    for d in 2..=k {
        rho[d - 1] = 1.0 / (d as f64 * (d as f64 - 1.0));
    }

    // Robust component tau, parameterized by R = c * ln(K/delta) * sqrt(K).
    let r = c * (kf / delta).ln() * kf.sqrt();
    // Spike location m = round(K/R), clamped to 1..=K.
    let mut m = if r > 0.0 {
        (kf / r).round() as usize
    } else {
        k
    };
    if m < 1 {
        m = 1;
    }
    if m > k {
        m = k;
    }

    let mut tau = vec![0.0_f64; k];
    // tau(d) = R/(d*K) for d = 1..=m-1
    for d in 1..m {
        tau[d - 1] = r / (d as f64 * kf);
    }
    // tau(m) = R * ln(R/delta) / K
    {
        let rln = if r > 0.0 { (r / delta).ln() } else { 0.0 };
        tau[m - 1] = r * rln / kf;
    }
    // tau(d) = 0 for d > m (already zero).

    // mu(d) = (rho + tau) / beta, beta = sum(rho + tau).
    let mut pmf = vec![0.0_f64; k];
    let mut beta = 0.0_f64;
    for d in 0..k {
        pmf[d] = rho[d] + tau[d];
        beta += pmf[d];
    }
    if beta <= 0.0 {
        // Degenerate fallback: uniform.
        return vec![1.0 / kf; k];
    }
    for p in pmf.iter_mut() {
        *p /= beta;
    }
    pmf
}

/// Sample a degree `d` in `1..=k` from the robust soliton PMF using `u` in
/// [0, 1) (inverse-CDF sampling).
fn sample_degree(pmf: &[f64], u: f64) -> usize {
    let mut acc = 0.0_f64;
    for (i, &p) in pmf.iter().enumerate() {
        acc += p;
        if u < acc {
            return i + 1;
        }
    }
    // Floating point slack: return the last degree.
    pmf.len()
}

/// Shared neighbour selection: given a `seed` and `num_blocks`, return the
/// droplet degree and the sorted, de-duplicated source-block indices it XORs.
/// MUST be identical on the encode and decode sides.
pub fn neighbours(seed: u32, num_blocks: usize) -> (u32, Vec<usize>) {
    let k = num_blocks;
    if k == 0 {
        return (0, Vec::new());
    }
    if k == 1 {
        return (1, vec![0]);
    }

    let mut rng = SplitMix64::new(seed as u64);
    let pmf = robust_soliton_pmf(k);

    // First draw selects the degree.
    let u = rng.next_f64();
    let mut d = sample_degree(&pmf, u);
    if d < 1 {
        d = 1;
    }
    if d > k {
        d = k;
    }

    // Sample d distinct indices in 0..k via partial Fisher–Yates over a
    // working permutation array.
    let mut perm: Vec<usize> = (0..k).collect();
    for i in 0..d {
        let j = i + (rng.next_below((k - i) as u64) as usize);
        perm.swap(i, j);
    }
    let mut idxs: Vec<usize> = perm[..d].to_vec();
    idxs.sort_unstable();

    (d as u32, idxs)
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
    pub fn encode(&self, blocks: &[Vec<u8>], seed: u32) -> Droplet {
        let (d, idxs) = neighbours(seed, self.num_blocks);
        let mut data = vec![0u8; self.block_size];
        for &i in &idxs {
            let block = &blocks[i];
            for (b, &src) in data.iter_mut().zip(block.iter()) {
                *b ^= src;
            }
        }
        Droplet {
            seed,
            degree: d,
            data,
        }
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
    pub fn decode(&self, droplets: &[Droplet]) -> Option<Vec<Vec<u8>>> {
        let k = self.num_blocks;
        if k == 0 {
            return Some(Vec::new());
        }

        // Working copies of each droplet's data and its set of still-unresolved
        // neighbour block indices.
        let mut work_data: Vec<Vec<u8>> = Vec::with_capacity(droplets.len());
        let mut unresolved: Vec<Vec<usize>> = Vec::with_capacity(droplets.len());
        for dr in droplets {
            let (_d, idxs) = neighbours(dr.seed, k);
            work_data.push(dr.data.clone());
            unresolved.push(idxs);
        }

        let mut recovered: Vec<Option<Vec<u8>>> = vec![None; k];
        let mut num_recovered = 0usize;

        loop {
            if num_recovered == k {
                break;
            }

            // Find droplets with exactly one unresolved neighbour (degree-1 in
            // the residual graph) and peel them.
            let mut progressed = false;
            for di in 0..droplets.len() {
                if unresolved[di].len() != 1 {
                    continue;
                }
                let block_idx = unresolved[di][0];
                if recovered[block_idx].is_some() {
                    // Stale; clear it out.
                    unresolved[di].clear();
                    continue;
                }

                // This droplet's current data IS the recovered block.
                let block_value = work_data[di].clone();
                recovered[block_idx] = Some(block_value.clone());
                num_recovered += 1;
                unresolved[di].clear();
                progressed = true;

                // Peel this block out of every other droplet that references it.
                for dj in 0..droplets.len() {
                    if dj == di {
                        continue;
                    }
                    if let Some(pos) = unresolved[dj].iter().position(|&x| x == block_idx) {
                        for (b, &src) in work_data[dj].iter_mut().zip(block_value.iter()) {
                            *b ^= src;
                        }
                        unresolved[dj].swap_remove(pos);
                    }
                }
            }

            if !progressed {
                // No ripple — stuck.
                break;
            }
        }

        if num_recovered == k {
            Some(recovered.into_iter().map(|b| b.unwrap()).collect())
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_blocks(k: usize, block_size: usize, seed: u32) -> Vec<Vec<u8>> {
        // Deterministic pseudo-random source blocks.
        let mut rng = SplitMix64::new(0xABCD_0000 ^ seed as u64);
        (0..k)
            .map(|_| {
                (0..block_size)
                    .map(|_| (rng.next_u64() & 0xFF) as u8)
                    .collect()
            })
            .collect()
    }

    #[test]
    fn neighbours_deterministic_and_in_range() {
        let k = 20;
        for seed in 0..200u32 {
            let (d, idxs) = neighbours(seed, k);
            let (d2, idxs2) = neighbours(seed, k);
            assert_eq!((d, &idxs), (d2, &idxs2));
            assert!(d >= 1 && (d as usize) <= k);
            assert_eq!(idxs.len(), d as usize);
            // distinct + sorted + in range
            for w in idxs.windows(2) {
                assert!(w[0] < w[1]);
            }
            assert!(idxs.iter().all(|&i| i < k));
        }
    }

    #[test]
    fn decode_all_droplets() {
        let k = 20;
        let block_size = 16;
        let blocks = make_blocks(k, block_size, 1);
        let enc = LtEncoder::new(k, block_size);
        let droplets: Vec<Droplet> = (0..(2 * k as u32))
            .map(|s| enc.encode(&blocks, s))
            .collect();
        let dec = LtDecoder::new(k, block_size);
        let out = dec.decode(&droplets).expect("should decode");
        assert_eq!(out, blocks);
    }

    #[test]
    fn decode_with_drops() {
        let k = 20;
        let block_size = 16;
        let blocks = make_blocks(k, block_size, 7);
        let enc = LtEncoder::new(k, block_size);

        // Generate 3*K droplets, then keep a pseudo-random ~2*K subset.
        let all: Vec<Droplet> = (0..(3 * k as u32))
            .map(|s| enc.encode(&blocks, s))
            .collect();
        let mut rng = SplitMix64::new(0xFEED_BEEF);
        let mut kept: Vec<Droplet> = all
            .into_iter()
            .filter(|_| rng.next_f64() < (2.0 / 3.0))
            .collect();
        // Ensure we kept at least 2*K; top up deterministically if short.
        let mut extra = 3 * k as u32;
        while kept.len() < 2 * k {
            kept.push(enc.encode(&blocks, extra));
            extra += 1;
        }

        let dec = LtDecoder::new(k, block_size);
        let out = dec.decode(&kept).expect("should decode from subset");
        assert_eq!(out, blocks);
    }

    #[test]
    fn single_block() {
        let k = 1;
        let block_size = 8;
        let blocks = make_blocks(k, block_size, 99);
        let enc = LtEncoder::new(k, block_size);
        let (d, idxs) = neighbours(12345, k);
        assert_eq!(d, 1);
        assert_eq!(idxs, vec![0]);
        let droplets: Vec<Droplet> = (0..3u32).map(|s| enc.encode(&blocks, s)).collect();
        let dec = LtDecoder::new(k, block_size);
        let out = dec.decode(&droplets).expect("should decode single");
        assert_eq!(out, blocks);
    }
}
