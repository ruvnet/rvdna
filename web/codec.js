// =============================================================================
// DNA Storage Codec Simulator — pure-JS codec (the core)
// =============================================================================
//
// This module is a faithful, dependency-free JavaScript mirror of the Rust
// crate at  examples/dna/src/storage/.  Each section below maps 1:1 to a Rust
// submodule so the browser demo is *honest* about what the real codec does:
//
//   constraints.rs  ->  encodeBytes / decodeBytes / gcContent / maxHomopolymerRun
//   gf256.rs        ->  ReedSolomon  (GF(256), prim poly 0x11D, generator 2)
//   channel.rs      ->  applyChannel (sub / ins / del / dropout / coverage)
//   consensus.rs    ->  consensus    (position-wise majority vote across reads)
//   mod.rs          ->  crc32, EncodeParams, Strand, DnaArchive, DecodeReport,
//                       DnaStorageCodec (encode / decode orchestrator)
//
//   fountain.rs     ->  SplitMix64 / robustSolitonPmf / neighbours / LtEncoder /
//                       LtDecoder  (Luby-Transform rateless erasure outer code)
//
// The outer erasure code is now a REAL Luby-Transform (LT) fountain code, ported
// bit-for-bit from fountain.rs: the same SplitMix64 PRNG, robust-soliton degree
// distribution, deterministic seed->neighbour selection, and belief-propagation
// peeling decoder. Encode emits seed-tagged droplets (XORs of source blocks);
// decode RS-corrects each read, dedupes droplets by seed, then peels them back
// into the source blocks — exactly the Rust pipeline. See README.md.
//
// Everything here is plain ES-module JavaScript: no build step, no network.
// =============================================================================

// -----------------------------------------------------------------------------
// 1) constraints.rs  —  reversible bytes <-> DNA mapping under bio constraints
// -----------------------------------------------------------------------------
//
// Homopolymer-free, GC-balanced base-3 transform. We never emit the same base
// twice in a row, so the maximum homopolymer run is always 1.
//
// alphabet: A=0, C=1, G=2, T=3.  A running `prev` base index starts at 0 (A).
// Each input byte (0..255) becomes 6 ternary digits, big-endian (3^6 = 729 > 255).
// For each trit t in {0,1,2}:  cur = (prev + 1 + t) % 4  -> guarantees cur != prev.

export const BASES = ['A', 'C', 'G', 'T'];
const BASE_INDEX = { A: 0, C: 1, G: 2, T: 3 };

/** Encode bytes into a homopolymer-free `ACGT` string. */
export function encodeBytes(data) {
  let prev = 0;          // running previous base index, init 0 (= 'A')
  let out = '';
  for (let k = 0; k < data.length; k++) {
    const b = data[k];
    // big-endian ternary digits: d[i] = floor(b / 3^(5-i)) % 3
    for (let i = 0; i < 6; i++) {
      const t = Math.floor(b / Math.pow(3, 5 - i)) % 3;
      const cur = (prev + 1 + t) % 4;   // cur != prev  => no homopolymer run
      out += BASES[cur];
      prev = cur;
    }
  }
  return out;
}

/**
 * Inverse of `encodeBytes`. Length must be a multiple of 6.
 *
 * `tolerant` (default false): a noisy read can contain a base transition that
 * the clean encoder never emits (prev == cur, i.e. trit value 3). In strict
 * mode this throws; in tolerant mode we clamp the illegal trit to 0 and keep
 * going, yielding a best-effort byte that the Reed–Solomon layer can then
 * correct. The demo's recovery path uses tolerant mode so framing survives.
 */
export function decodeBytes(seq, tolerant = false) {
  if (seq.length % 6 !== 0) {
    throw new Error(`decodeBytes: length ${seq.length} not a multiple of 6`);
  }
  const out = new Uint8Array(seq.length / 6);
  let prev = 0;
  let bi = 0;
  for (let k = 0; k < seq.length; k += 6) {
    let b = 0;
    for (let i = 0; i < 6; i++) {
      let cur = BASE_INDEX[seq[k + i]];
      if (cur === undefined) {
        if (!tolerant) throw new Error(`decodeBytes: bad base '${seq[k + i]}'`);
        cur = (prev + 1) % 4; // treat unknown char as trit 0
      }
      let t = (cur + 4 - prev - 1) % 4;   // inverse of the encode step
      if (t > 2) {
        if (!tolerant) throw new Error('decodeBytes: invalid trit (corrupt sequence)');
        t = 0; // clamp illegal transition; RS will repair the resulting byte
      }
      b = b * 3 + t;
      prev = cur;
    }
    out[bi++] = b;
  }
  return out;
}

/** GC content fraction (0.0–1.0) of a sequence. */
export function gcContent(seq) {
  if (seq.length === 0) return 0;
  let gc = 0;
  for (let i = 0; i < seq.length; i++) {
    const c = seq[i];
    if (c === 'G' || c === 'C') gc++;
  }
  return gc / seq.length;
}

/** Longest run of an identical nucleotide. */
export function maxHomopolymerRun(seq) {
  if (seq.length === 0) return 0;
  let best = 1, run = 1;
  for (let i = 1; i < seq.length; i++) {
    run = seq[i] === seq[i - 1] ? run + 1 : 1;
    if (run > best) best = run;
  }
  return best;
}

// -----------------------------------------------------------------------------
// 2) gf256.rs  —  Reed–Solomon inner code over GF(256)
// -----------------------------------------------------------------------------
//
// Field: GF(2^8) with primitive polynomial 0x11D, generator (primitive elt) 2.
// Systematic encoder appends `nsym` parity bytes; decoder uses the classic
// syndrome / Berlekamp–Massey / Chien-search / Forney pipeline and corrects up
// to `nsym/2` substitution errors per codeword.

// Log/antilog tables for GF(256). GF_EXP is doubled so we can index
// GF_EXP[GF_LOG[a] + GF_LOG[b]] without taking a modulo.
const GF_EXP = new Uint8Array(512);
const GF_LOG = new Uint8Array(256);
(function initGF() {
  let x = 1;
  for (let i = 0; i < 255; i++) {
    GF_EXP[i] = x;
    GF_LOG[x] = i;
    x <<= 1;                    // multiply by the generator (2)
    if (x & 0x100) x ^= 0x11D;  // reduce mod the primitive polynomial 0x11D
  }
  for (let i = 255; i < 512; i++) GF_EXP[i] = GF_EXP[i - 255];
})();

function gfMul(a, b) {
  if (a === 0 || b === 0) return 0;
  return GF_EXP[GF_LOG[a] + GF_LOG[b]];
}
function gfInv(a) {
  return GF_EXP[255 - GF_LOG[a]];
}
function gfPow(a, n) {
  // a^n in the field; n may be negative.
  let e = (GF_LOG[a] * n) % 255;
  if (e < 0) e += 255;
  return GF_EXP[e];
}

// --- Polynomials are stored coefficient-array, INDEX 0 = HIGHEST degree. ------
// (This is the convention used by the canonical Wikiversity RS reference, which
// this decoder follows closely so the GF(256) logic is easy to audit.)

function polyScale(p, x) {
  const r = new Uint8Array(p.length);
  for (let i = 0; i < p.length; i++) r[i] = gfMul(p[i], x);
  return r;
}
// Add (== subtract, XOR) two polynomials, aligned on the *high* degree (index 0).
function polyAdd(p, q) {
  const r = new Uint8Array(Math.max(p.length, q.length));
  for (let i = 0; i < p.length; i++) r[i + r.length - p.length] = p[i];
  for (let i = 0; i < q.length; i++) r[i + r.length - q.length] ^= q[i];
  return r;
}
function polyMul(p, q) {
  const r = new Uint8Array(p.length + q.length - 1);
  for (let j = 0; j < q.length; j++) {
    if (q[j] === 0) continue;
    for (let i = 0; i < p.length; i++) {
      if (p[i] !== 0) r[i + j] ^= gfMul(p[i], q[j]);
    }
  }
  return r;
}
// Evaluate p(x) via Horner (index 0 = highest degree).
function polyEval(p, x) {
  let y = p[0];
  for (let i = 1; i < p.length; i++) y = gfMul(y, x) ^ p[i];
  return y;
}

/** Systematic Reed–Solomon codec over GF(256) appending `nsym` parity bytes. */
export class ReedSolomon {
  constructor(nsym) {
    this.nsym = nsym;
    this.gen = ReedSolomon.generator(nsym);
  }

  /** Generator polynomial g(x) = prod_{i=0}^{nsym-1} (x - a^i). */
  static generator(nsym) {
    let g = Uint8Array.of(1);
    for (let i = 0; i < nsym; i++) {
      g = polyMul(g, Uint8Array.of(1, gfPow(2, i))); // (x + a^i)
    }
    return g;
  }

  /** Append `nsym` parity bytes; returns data.length + nsym bytes (<= 255). */
  encode(data) {
    // Remainder of (data * x^nsym) / gen via synthetic division.
    const out = new Uint8Array(data.length + this.nsym);
    out.set(data, 0);
    for (let i = 0; i < data.length; i++) {
      const coef = out[i];
      if (coef !== 0) {
        for (let j = 1; j < this.gen.length; j++) {
          out[i + j] ^= gfMul(this.gen[j], coef);
        }
      }
    }
    // The systematic message bytes were clobbered by the division; restore them.
    out.set(data, 0);
    return out;
  }

  /** Syndromes S_i = received(a^i), i = 0..nsym-1. */
  syndromes(received) {
    const s = new Uint8Array(this.nsym);
    let allZero = true;
    for (let i = 0; i < this.nsym; i++) {
      s[i] = polyEval(received, gfPow(2, i));
      if (s[i] !== 0) allZero = false;
    }
    return { s, allZero };
  }

  /**
   * Correct up to nsym/2 substitution errors and return the message
   * (parity stripped), or throw if the codeword is unrecoverable.
   */
  decode(received) {
    received = Uint8Array.from(received);
    const { s, allZero } = this.syndromes(received);
    if (allZero) return received.slice(0, received.length - this.nsym);

    // Syndrome polynomial, index 0 = highest degree: [S_{nsym-1} .. S_0].
    const synd = new Uint8Array(this.nsym);
    for (let i = 0; i < this.nsym; i++) synd[i] = s[this.nsym - 1 - i];

    // --- Berlekamp–Massey: error-locator sigma(x). ---
    let sigma = Uint8Array.of(1);
    let old = Uint8Array.of(1);
    for (let i = 0; i < this.nsym; i++) {
      old = polyAppendZero(old); // old(x) *= x
      // discrepancy delta = sum_j sigma[j] * S_{i-j}
      let delta = 0;
      for (let j = 0; j < sigma.length; j++) {
        delta ^= gfMul(sigma[sigma.length - 1 - j], s[i - j] || 0);
      }
      if (delta !== 0) {
        if (old.length > sigma.length) {
          const newOld = polyScale(sigma, gfInv(delta));
          sigma = polyAdd(sigma, polyScale(old, delta));
          old = newOld;
        } else {
          sigma = polyAdd(sigma, polyScale(old, delta));
        }
      }
    }
    sigma = polyTrim(sigma);
    const errCount = sigma.length - 1;
    if (errCount === 0 || errCount * 2 > this.nsym) {
      throw new Error('RS: too many errors to correct');
    }

    // --- Chien search: error positions are i where sigma(a^{-i}) == 0. ---
    const n = received.length;
    const positions = [];
    for (let i = 0; i < n; i++) {
      if (polyEval(sigma, gfPow(2, -i)) === 0) positions.push(i); // power of a^-i
    }
    if (positions.length !== errCount) {
      throw new Error('RS: could not locate all errors');
    }

    // --- Forney: error magnitudes. omega(x) = S(x)*sigma(x) mod x^nsym. ---
    let omega = polyMul(synd, sigma);
    omega = omega.slice(Math.max(0, omega.length - this.nsym)); // mod x^nsym
    omega = polyTrim(omega);
    const sigmaPrime = formalDerivative(sigma);

    const corrected = Uint8Array.from(received);
    for (const i of positions) {
      const xi = gfPow(2, i);          // a^i  (error locator value X_k)
      const xiInv = gfInv(xi);         // a^{-i}
      const num = polyEval(omega, xiInv);
      const den = polyEval(sigmaPrime, xiInv);
      if (den === 0) throw new Error('RS: Forney denominator zero');
      const mag = gfMul(xi, gfMul(num, gfInv(den)));
      // position i counts from the END (highest power) -> array index n-1-i
      corrected[n - 1 - i] ^= mag;
    }

    const check = this.syndromes(corrected);
    if (!check.allZero) throw new Error('RS: correction failed verification');
    return corrected.slice(0, corrected.length - this.nsym);
  }
}

// --- small poly helpers used only by the RS decoder --------------------------
function polyAppendZero(p) {
  // multiply by x: append a low-order zero coefficient (index grows at the tail)
  const r = new Uint8Array(p.length + 1);
  r.set(p, 0);
  return r;
}
function polyTrim(p) {
  let i = 0;
  while (i < p.length - 1 && p[i] === 0) i++;
  return p.slice(i);
}
function formalDerivative(p) {
  // index 0 = highest degree; degree of p[i] is (L-1-i).
  // In GF(2^k), d/dx keeps coefficients of odd degree, zeros even ones.
  const L = p.length;
  if (L <= 1) return Uint8Array.of(0);
  const r = new Uint8Array(L - 1);
  for (let i = 0; i < L - 1; i++) {
    const degree = L - 1 - i;
    r[i] = (degree & 1) ? p[i] : 0;
  }
  return r;
}

// -----------------------------------------------------------------------------
// 3) mod.rs  —  CRC32 (IEEE), reflected, poly 0xEDB88320
// -----------------------------------------------------------------------------
export function crc32(data) {
  let crc = 0xffffffff;
  for (let i = 0; i < data.length; i++) {
    crc ^= data[i];
    for (let k = 0; k < 8; k++) {
      const mask = -(crc & 1);
      crc = (crc >>> 1) ^ (0xedb88320 & mask);
    }
  }
  return (~crc) >>> 0;
}

// -----------------------------------------------------------------------------
// 4) channel.rs  —  synthetic sequencing-error channel
// -----------------------------------------------------------------------------
//
// Deterministic via a tiny mulberry32 PRNG seeded per run. Each surviving
// strand contributes up to `coverage` independently-noisy reads. Returns a
// shuffled read pool. We tag the returned reads with the *true* logical index
// so the demo can colour diffs and regroup — the index is also carried inside
// each strand's own header (see DnaStorageCodec), matching the Rust contract
// where reads are regrouped from their own data.

export class Rng {
  constructor(seed) { this.s = (seed >>> 0) || 1; }
  next() {
    // mulberry32
    let t = (this.s += 0x6d2b79f5) >>> 0;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }
  int(n) { return Math.floor(this.next() * n); }
}

export const DEFAULT_MODEL = {
  pSub: 0.01, pIns: 0.002, pDel: 0.002, pDrop: 0.05, coverage: 1,
};

/**
 * Push every strand through the channel.
 * `strands` is an array of { index, sequence }.
 * Returns an array of reads: { index, sequence } (shuffled).
 */
export function applyChannel(strands, model, seed) {
  const rng = new Rng(seed);
  const reads = [];
  for (const strand of strands) {
    if (rng.next() < model.pDrop) continue; // whole-strand dropout
    const copies = Math.max(1, model.coverage | 0);
    for (let c = 0; c < copies; c++) {
      reads.push({ index: strand.index, sequence: mutateSeq(strand.sequence, model, rng) });
    }
  }
  // Fisher–Yates shuffle for a realistic unordered read pool
  for (let i = reads.length - 1; i > 0; i--) {
    const j = rng.int(i + 1);
    [reads[i], reads[j]] = [reads[j], reads[i]];
  }
  return reads;
}

function mutateSeq(seq, model, rng) {
  let out = '';
  for (let i = 0; i < seq.length; i++) {
    if (rng.next() < model.pDel) continue;                 // deletion: drop base
    if (rng.next() < model.pIns) out += BASES[rng.int(4)]; // insertion before base
    let c = seq[i];
    if (rng.next() < model.pSub) {                          // substitution
      let nb = BASES[rng.int(4)];
      while (nb === c) nb = BASES[rng.int(4)];
      c = nb;
    }
    out += c;
  }
  return out;
}

// -----------------------------------------------------------------------------
// 5) consensus.rs  —  read clustering + majority-vote consensus
// -----------------------------------------------------------------------------
//
// We regroup reads by their strand header index (the decoder reads each strand's
// own header), then take a position-wise majority vote up to the modal read
// length. This cancels random substitutions whenever coverage > 1 and helps
// realign reads that suffered small indels.

/** Position-wise majority-vote consensus over a group of read sequences. */
export function consensus(reads) {
  if (reads.length === 0) return '';
  if (reads.length === 1) return reads[0];
  // modal length keeps indel-shortened/lengthened reads from skewing the result
  const counts = {};
  for (const r of reads) counts[r.length] = (counts[r.length] || 0) + 1;
  let modalLen = reads[0].length, best = 0;
  for (const [len, n] of Object.entries(counts)) {
    if (n > best) { best = n; modalLen = +len; }
  }
  let out = '';
  for (let pos = 0; pos < modalLen; pos++) {
    const tally = { A: 0, C: 0, G: 0, T: 0 };
    for (const r of reads) {
      if (pos < r.length && tally[r[pos]] !== undefined) tally[r[pos]]++;
    }
    let bestBase = 'A', bestN = -1;
    for (const b of BASES) if (tally[b] > bestN) { bestN = tally[b]; bestBase = b; }
    out += bestBase;
  }
  return out;
}

/**
 * Greedy single-linkage clustering of reads by normalized Hamming similarity,
 * mirroring `consensus::cluster` in the Rust crate. Two reads join the same
 * cluster when their prefix mismatch fraction is <= `maxDist`. Returns an array
 * of clusters, each an array of indices into `reads`.
 */
export function cluster(reads, maxDist = 0.3) {
  const clusters = [];
  for (let i = 0; i < reads.length; i++) {
    let placed = false;
    for (const c of clusters) {
      // compare against the cluster's first (representative) member
      const rep = reads[c[0]];
      const a = reads[i];
      const m = Math.min(rep.length, a.length);
      if (m === 0) continue;
      let diff = Math.abs(rep.length - a.length);
      for (let p = 0; p < m; p++) if (rep[p] !== a[p]) diff++;
      const denom = Math.max(rep.length, a.length);
      if (denom > 0 && diff / denom <= maxDist) { c.push(i); placed = true; break; }
    }
    if (!placed) clusters.push([i]);
  }
  return clusters;
}

// -----------------------------------------------------------------------------
// fountain.rs  —  Luby-Transform (LT) rateless erasure outer code
// -----------------------------------------------------------------------------
//
// Ported bit-for-bit from examples/dna/src/storage/fountain.rs. The encoder
// turns `numBlocks` source blocks into an unbounded stream of seed-tagged
// droplets (each the XOR of a seed-determined subset of blocks). The decoder
// peels droplets (belief propagation) back into the source blocks. Encode and
// decode agree because both derive the neighbour set deterministically from the
// droplet seed via `neighbours`, using an inline SplitMix64 PRNG.

const U64_MASK = (1n << 64n) - 1n;

/** Deterministic SplitMix64 PRNG (64-bit state via BigInt). Matches Rust. */
export class SplitMix64 {
  constructor(seed) {
    // accept BigInt or number; mask to 64 bits.
    this.state = BigInt.asUintN(64, typeof seed === 'bigint' ? seed : BigInt(seed >>> 0));
  }

  /** Rust `next_u64`: returns a BigInt in [0, 2^64). */
  nextU64() {
    this.state = (this.state + 0x9E3779B97F4A7C15n) & U64_MASK;
    let z = this.state;
    z = ((z ^ (z >> 30n)) * 0xBF58476D1CE4E5B9n) & U64_MASK;
    z = ((z ^ (z >> 27n)) * 0x94D049BB133111EBn) & U64_MASK;
    z = z ^ (z >> 31n);
    return z & U64_MASK;
  }

  /** Rust `next_f64`: uniform double in [0, 1) from the top 53 bits. */
  nextF64() {
    return Number(this.nextU64() >> 11n) / 9007199254740992; // 2^53
  }

  /** Rust `next_below`: uniform integer in 0..n (n > 0), returns a Number. */
  nextBelow(n) {
    return Number(this.nextU64() % BigInt(n));
  }
}

/**
 * Build the normalized Robust Soliton PMF over degrees 1..=k (index i is the
 * probability of degree i+1). Identical math to `robust_soliton_pmf` in Rust.
 */
export function robustSolitonPmf(k) {
  if (k === 0) return [];
  if (k === 1) return [1.0];

  const c = 0.03;
  const delta = 0.05;
  const kf = k;

  // Ideal soliton.
  const rho = new Array(k).fill(0);
  rho[0] = 1.0 / kf;
  for (let d = 2; d <= k; d++) rho[d - 1] = 1.0 / (d * (d - 1));

  const r = c * Math.log(kf / delta) * Math.sqrt(kf);
  let m = r > 0 ? Math.round(kf / r) : k;
  if (m < 1) m = 1;
  if (m > k) m = k;

  const tau = new Array(k).fill(0);
  for (let d = 1; d < m; d++) tau[d - 1] = r / (d * kf);
  {
    const rln = r > 0 ? Math.log(r / delta) : 0.0;
    tau[m - 1] = (r * rln) / kf;
  }

  const pmf = new Array(k).fill(0);
  let beta = 0.0;
  for (let d = 0; d < k; d++) { pmf[d] = rho[d] + tau[d]; beta += pmf[d]; }
  if (beta <= 0.0) return new Array(k).fill(1.0 / kf);
  for (let d = 0; d < k; d++) pmf[d] /= beta;
  return pmf;
}

/** Inverse-CDF degree sampling from `pmf` using u in [0,1). Matches Rust. */
export function sampleDegree(pmf, u) {
  let acc = 0.0;
  for (let i = 0; i < pmf.length; i++) {
    acc += pmf[i];
    if (u < acc) return i + 1;
  }
  return pmf.length;
}

/**
 * Shared neighbour selection: from `seed` and `numBlocks`, return the droplet
 * degree and the sorted, de-duplicated source-block indices it XORs. MUST be
 * identical on encode and decode. Matches `neighbours` in Rust.
 */
export function neighbours(seed, numBlocks) {
  const k = numBlocks;
  if (k === 0) return { degree: 0, idxs: [] };
  if (k === 1) return { degree: 1, idxs: [0] };

  const rng = new SplitMix64(BigInt(seed >>> 0));
  const pmf = robustSolitonPmf(k);

  const u = rng.nextF64();
  let d = sampleDegree(pmf, u);
  if (d < 1) d = 1;
  if (d > k) d = k;

  // Partial Fisher–Yates over a working permutation array.
  const perm = new Array(k);
  for (let i = 0; i < k; i++) perm[i] = i;
  for (let i = 0; i < d; i++) {
    const j = i + rng.nextBelow(k - i);
    const tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
  }
  const idxs = perm.slice(0, d);
  idxs.sort((a, b) => a - b);
  return { degree: d, idxs };
}

/** XOR `src` into `dst` (Uint8Array, in place) up to the shorter length. */
function xorInto(dst, src) {
  const n = Math.min(dst.length, src.length);
  for (let i = 0; i < n; i++) dst[i] ^= src[i];
}

/** LT encoder over `numBlocks` source blocks of `blockSize` bytes each. */
export class LtEncoder {
  constructor(numBlocks, blockSize) {
    this.numBlocks = numBlocks;
    this.blockSize = blockSize;
  }
  /** Deterministically produce the droplet { seed, degree, data } for `seed`. */
  encode(blocks, seed) {
    const { degree, idxs } = neighbours(seed, this.numBlocks);
    const data = new Uint8Array(this.blockSize);
    for (const i of idxs) xorInto(data, blocks[i]);
    return { seed: seed >>> 0, degree, data };
  }
}

/** LT belief-propagation ("peeling") decoder. */
export class LtDecoder {
  constructor(numBlocks, blockSize) {
    this.numBlocks = numBlocks;
    this.blockSize = blockSize;
  }
  /** Peel `droplets`; returns an array of source blocks, or null if stuck. */
  decode(droplets) {
    const k = this.numBlocks;
    if (k === 0) return [];

    const workData = [];
    const unresolved = [];
    for (const dr of droplets) {
      const { idxs } = neighbours(dr.seed, k);
      workData.push(Uint8Array.from(dr.data));
      unresolved.push(idxs.slice());
    }

    const recovered = new Array(k).fill(null);
    let numRecovered = 0;

    for (;;) {
      if (numRecovered === k) break;
      let progressed = false;
      for (let di = 0; di < droplets.length; di++) {
        if (unresolved[di].length !== 1) continue;
        const blockIdx = unresolved[di][0];
        if (recovered[blockIdx] !== null) { unresolved[di].length = 0; continue; }

        const blockValue = Uint8Array.from(workData[di]);
        recovered[blockIdx] = blockValue;
        numRecovered++;
        unresolved[di].length = 0;
        progressed = true;

        // Peel this block out of every other droplet referencing it.
        for (let dj = 0; dj < droplets.length; dj++) {
          if (dj === di) continue;
          const pos = unresolved[dj].indexOf(blockIdx);
          if (pos !== -1) {
            xorInto(workData[dj], blockValue);
            // swap_remove
            unresolved[dj][pos] = unresolved[dj][unresolved[dj].length - 1];
            unresolved[dj].pop();
          }
        }
      }
      if (!progressed) break;
    }

    if (numRecovered === k) return recovered;
    return null;
  }
}

// -----------------------------------------------------------------------------
// 6) mod.rs  —  orchestrator: DnaStorageCodec (encode / decode)
// -----------------------------------------------------------------------------
//
// Per-strand wire format (bytes, before DNA mapping), matching mod.rs exactly:
//
//   [ index(4) | seed(4) | droplet payload(block_size) | RS parity(rs_parity) ]
//      u32 BE strand id    u32 BE fountain seed   LT droplet      GF(256) parity
//
// The 8-byte header is itself protected by RS (it sits inside the codeword), so
// the decoder recovers the exact fountain `seed` after RS correction; droplets
// are deduped by seed and peeled by the LT decoder to recover the source blocks.

const HEADER_LEN = 8; // index(4) + seed(4)

export const DEFAULT_PARAMS = {
  blockSize: 32,    // payload bytes per strand
  rsParity: 8,      // RS parity bytes per strand (corrects up to 4 subs)
  // Redundancy factor: emit ceil(overhead * numBlocks) droplets. Now that the
  // outer code is a REAL LT fountain (not repetition), this matches the Rust
  // default of 2.5.
  overhead: 2.5,
  maxHomopolymer: 1,
  seed: 0xc0ffee,
};

function u32be(n) { n >>>= 0; return [(n >>> 24) & 0xff, (n >>> 16) & 0xff, (n >>> 8) & 0xff, n & 0xff]; }
function readU32be(a, o) { return ((a[o] << 24) | (a[o + 1] << 16) | (a[o + 2] << 8) | a[o + 3]) >>> 0; }

export class DnaStorageCodec {
  constructor(params = {}) {
    this.params = { ...DEFAULT_PARAMS, ...params };
    this.rs = new ReedSolomon(this.params.rsParity);
  }

  /**
   * Encode raw bytes (Uint8Array) into a DNA archive.
   * Returns { filename, byteLen, crc32, numBlocks, params, strands }
   * where strands = [{ index, sequence }] and `index` is the *physical*
   * strand id; the logical block id is index % numBlocks.
   */
  encode(filename, data) {
    const { blockSize, rsParity, overhead, seed } = this.params;
    const byteLen = data.length;
    const checksum = crc32(data);

    const codewordLen = HEADER_LEN + blockSize + rsParity;
    if (codewordLen > 255) {
      throw new Error(
        `strand codeword ${codewordLen} bytes exceeds GF(256) limit of 255; ` +
        `reduce blockSize or rsParity`);
    }

    // Split payload into fixed-size source blocks (zero-padded tail).
    const numBlocks = Math.max(1, Math.ceil(byteLen / blockSize));
    const blocks = [];
    for (let blk = 0; blk < numBlocks; blk++) {
      const block = new Uint8Array(blockSize);
      const start = blk * blockSize;
      for (let i = 0; i < blockSize; i++) {
        const p = start + i;
        block[i] = p < byteLen ? data[p] : 0;
      }
      blocks.push(block);
    }

    const encoder = new LtEncoder(numBlocks, blockSize);
    const baseSeed = seed >>> 0;

    // Mirror Rust headroom: ceil(overhead*K) droplets, with additive floor.
    let numStrands = Math.ceil(Math.max(1.0, overhead) * numBlocks);
    numStrands = Math.max(numStrands, numBlocks + 8);

    const strands = [];
    for (let i = 0; i < numStrands; i++) {
      const dropletSeed = (baseSeed + i) >>> 0; // (mod 2^32)
      const droplet = encoder.encode(blocks, dropletSeed);

      // frame = [index(4) | seed(4) | droplet.data]
      const frame = new Uint8Array(HEADER_LEN + blockSize);
      frame.set(u32be(i), 0);
      frame.set(u32be(dropletSeed), 4);
      frame.set(droplet.data, HEADER_LEN);

      const codeword = this.rs.encode(frame); // appends rsParity bytes
      strands.push({ index: i, sequence: encodeBytes(codeword) });
    }

    return {
      filename, byteLen, crc32: checksum, numBlocks,
      params: { ...this.params }, strands,
    };
  }

  /**
   * Decode a pool of (noisy) reads back into the original bytes.
   * `reads` = [{ index, sequence }]  (index is only used as a fallback hint;
   * primary regrouping uses the RS-protected header inside each read).
   * Returns a DecodeReport.
   */
  decode(archive, reads) {
    const { blockSize } = archive.params;
    const numBlocks = archive.numBlocks;
    const rsParity = archive.params.rsParity;
    const rs = new ReedSolomon(rsParity);
    const codewordLen = HEADER_LEN + blockSize + rsParity;
    const expectedBases = codewordLen * 6;

    // reads may be raw ACGT strings or { index, sequence } objects (the demo's
    // channel emits the latter). Normalise to plain sequence strings.
    const seqs = reads.map(r => (typeof r === 'string' ? r : r.sequence));

    const droplets = [];
    const seenSeeds = new Set();
    let errorsCorrected = 0;

    // Demap one sequence -> RS-decode -> droplet, deduped by fountain seed.
    // Returns true if a *new* droplet was recovered.
    const tryRecover = (seq) => {
      // Normalise length to the expected codeword (indels shift length).
      let usable;
      if (seq.length >= expectedBases) usable = expectedBases;
      else usable = Math.floor(seq.length / 6) * 6;
      if (usable === 0) return false;
      const trimmed = seq.slice(0, usable);

      let bytes;
      try { bytes = decodeBytes(trimmed, true); } catch { return false; }
      if (bytes.length < codewordLen) return false;

      let frame;
      try {
        const before = bytes.slice(0, codewordLen);
        frame = rs.decode(before);
        const fixed = rs.encode(frame);
        for (let i = 0; i < fixed.length; i++) if (fixed[i] !== before[i]) errorsCorrected++;
      } catch { return false; }

      if (frame.length < HEADER_LEN) return false;
      const seed = readU32be(frame, 4);
      if (seenSeeds.has(seed)) return false; // already have this strand
      seenSeeds.add(seed);

      const payload = new Uint8Array(blockSize);
      payload.set(frame.subarray(HEADER_LEN, Math.min(frame.length, HEADER_LEN + blockSize)));
      const { degree } = neighbours(seed, numBlocks);
      droplets.push({ seed, degree, data: payload });
      return true;
    };

    // Pass 1 — decode every read independently. RS recovers each strand's exact
    // [index|seed] header, so identity comes from the code, not fuzzy matching.
    // Reads RS can't fix are deferred to the salvage pass.
    const residual = [];
    for (let i = 0; i < seqs.length; i++) {
      if (!tryRecover(seqs[i])) residual.push(i);
    }

    // Pass 2 (salvage) — when coverage > 1 the same strand was sequenced several
    // times. Cluster the leftover noisy reads, majority-vote a consensus per
    // cluster to cancel random substitutions, then retry RS. Only ever ADDS
    // strands, so it never harms the clean path. Mirrors Rust pass-2.
    if (residual.length > 1) {
      const residualReads = residual.map(i => seqs[i]);
      for (const group of cluster(residualReads, 0.3)) {
        if (group.length < 2) continue;
        const members = group.map(g => residualReads[g]);
        const cons = consensus(members);
        tryRecover(cons);
      }
    }

    const strandsRecovered = droplets.length;

    // Fountain peeling to recover the source blocks.
    const decoder = new LtDecoder(numBlocks, blockSize);
    const recoveredBlocks = decoder.decode(droplets);

    let bytes, blocksRecovered;
    if (recoveredBlocks !== null) {
      const out = new Uint8Array(numBlocks * blockSize);
      for (let b = 0; b < numBlocks; b++) out.set(recoveredBlocks[b], b * blockSize);
      bytes = out.slice(0, archive.byteLen);
      blocksRecovered = numBlocks;
    } else {
      bytes = new Uint8Array(0);
      blocksRecovered = 0;
    }

    const crcOk = blocksRecovered === numBlocks && bytes.length > 0 &&
      crc32(bytes) === archive.crc32;

    return {
      bytes, crcOk,
      blocksRecovered,
      numBlocks,
      strandsRecovered,
      errorsCorrected,
      readsIn: reads.length,
    };
  }
}

// -----------------------------------------------------------------------------
// Archive-level derived stats (mirrors DnaArchive methods in mod.rs)
// -----------------------------------------------------------------------------
export function totalBases(archive) {
  return archive.strands.reduce((a, s) => a + s.sequence.length, 0);
}
export function bitsPerBase(archive) {
  const bases = totalBases(archive);
  return bases === 0 ? 0 : (archive.byteLen * 8) / bases;
}
export function meanGc(archive) {
  if (archive.strands.length === 0) return 0;
  let sum = 0;
  for (const s of archive.strands) sum += gcContent(s.sequence);
  return sum / archive.strands.length;
}
export function meanMaxHomopolymer(archive) {
  if (archive.strands.length === 0) return 0;
  let best = 0;
  for (const s of archive.strands) best = Math.max(best, maxHomopolymerRun(s.sequence));
  return best;
}

// -----------------------------------------------------------------------------
// Self-test (run in a console with:  node codec.js  ... or import & call).
// Kept here so the codec is auditable. Harmless in the browser (guarded).
// -----------------------------------------------------------------------------
export function selfTest() {
  const log = [];
  const assert = (c, m) => { if (!c) throw new Error('SELFTEST FAIL: ' + m); log.push('ok: ' + m); };

  // CRC known vector
  assert(crc32(new TextEncoder().encode('123456789')) === 0xcbf43926, 'crc32("123456789")==0xCBF43926');

  // constraints round-trip + homopolymer-free + reversibility for all bytes
  const all = new Uint8Array(256); for (let i = 0; i < 256; i++) all[i] = i;
  const dna = encodeBytes(all);
  assert(maxHomopolymerRun(dna) === 1, 'encoded stream is homopolymer-free (run==1)');
  const back = decodeBytes(dna);
  let rt = true; for (let i = 0; i < 256; i++) if (back[i] !== i) rt = false;
  assert(rt, 'encodeBytes/decodeBytes round-trips all 256 byte values');

  // GF(256) inverse sanity (a * a^-1 == 1 for all non-zero a)
  let invOk = true;
  for (let a = 1; a < 256; a++) if (gfMul(a, gfInv(a)) !== 1) invOk = false;
  assert(invOk, 'gfInv: a * a^-1 == 1 for all 255 non-zero field elements');

  // RS encode/decode with injected errors
  const rs = new ReedSolomon(8); // corrects 4 errors
  const msg = new Uint8Array(20); for (let i = 0; i < 20; i++) msg[i] = (i * 7 + 3) & 0xff;
  const cw = rs.encode(msg);
  cw[1] ^= 0x5a; cw[9] ^= 0xff; cw[15] ^= 0x01; cw[27] ^= 0x80; // 4 substitutions
  const dec = rs.decode(cw);
  let rsOk = true; for (let i = 0; i < 20; i++) if (dec[i] !== msg[i]) rsOk = false;
  assert(rsOk, 'ReedSolomon corrects 4 substitution errors (nsym=8)');

  // fountain (LT) determinism + peeling sanity, matching fountain.rs
  {
    // neighbours: deterministic, in-range, sorted, distinct.
    const k = 20;
    let okN = true;
    for (let s = 0; s < 200; s++) {
      const a = neighbours(s, k), b = neighbours(s, k);
      if (a.degree !== b.degree || a.idxs.length !== a.degree) okN = false;
      if (!(a.degree >= 1 && a.degree <= k)) okN = false;
      for (let i = 1; i < a.idxs.length; i++) if (a.idxs[i - 1] >= a.idxs[i]) okN = false;
      for (const idx of a.idxs) if (idx < 0 || idx >= k) okN = false;
      for (let i = 0; i < a.idxs.length; i++) if (a.idxs[i] !== b.idxs[i]) okN = false;
    }
    assert(okN, 'fountain neighbours are deterministic, sorted, distinct, in range');

    const one = neighbours(12345, 1);
    assert(one.degree === 1 && one.idxs.length === 1 && one.idxs[0] === 0,
      'fountain neighbours(k=1) == degree 1, [0]');

    // SplitMix64 known sequence: state seeded 0, first next_u64.
    const sm = new SplitMix64(0n);
    assert(sm.nextU64() === 16294208416658607535n, 'SplitMix64(0) first next_u64 matches Rust');

    // LT encode -> peel round-trips deterministic source blocks.
    const blockSize = 16;
    const blocks = [];
    const brng = new SplitMix64(0xABCD0000n);
    for (let i = 0; i < k; i++) {
      const blk = new Uint8Array(blockSize);
      for (let j = 0; j < blockSize; j++) blk[j] = Number(brng.nextU64() & 0xFFn);
      blocks.push(blk);
    }
    const enc = new LtEncoder(k, blockSize);
    const drs = [];
    for (let s = 0; s < 2 * k; s++) drs.push(enc.encode(blocks, s));
    const dec = new LtDecoder(k, blockSize);
    const outBlocks = dec.decode(drs);
    let peelOk = outBlocks !== null;
    if (peelOk) for (let i = 0; i < k; i++)
      for (let j = 0; j < blockSize; j++) if (outBlocks[i][j] !== blocks[i][j]) peelOk = false;
    assert(peelOk, 'LT encode+peel recovers all source blocks from 2K droplets');
  }

  // end-to-end with the channel (small payload)
  const codec = new DnaStorageCodec();
  const data = new TextEncoder().encode('I stored my GitHub repo in DNA. '.repeat(4));
  const arc = codec.encode('hello.txt', data);
  const reads = applyChannel(arc.strands, { pSub: 0.02, pIns: 0, pDel: 0, pDrop: 0.2, coverage: 3 }, 42);
  const rep = codec.decode(arc, reads);
  assert(rep.crcOk, 'end-to-end recovery succeeds with sub+dropout (crcOk)');

  // full multi-block file survives sub + strand dropout via the REAL fountain
  {
    const big = new Uint8Array(600);
    const rng2 = new Rng(0x1234abcd);
    for (let i = 0; i < big.length; i++) big[i] = rng2.int(256);
    const codec2 = new DnaStorageCodec();
    const arc2 = codec2.encode('blob.bin', big);
    const reads2 = applyChannel(arc2.strands,
      { pSub: 0.02, pIns: 0, pDel: 0, pDrop: 0.15, coverage: 3 }, 99);
    const rep2 = codec2.decode(arc2, reads2);
    assert(rep2.crcOk, 'full 600-byte file recovers byte-exact via LT fountain (sub+dropout)');
    let exact = rep2.bytes.length === big.length;
    if (exact) for (let i = 0; i < big.length; i++) if (rep2.bytes[i] !== big[i]) exact = false;
    assert(exact, 'recovered 600-byte payload is byte-for-byte identical');
  }

  return log;
}

// Node CLI entry point for auditing (ignored by browsers):
if (typeof process !== 'undefined' && process.argv && process.argv[1] &&
    process.argv[1].endsWith('codec.js')) {
  try {
    const out = selfTest();
    console.log(out.join('\n'));
    console.log('\nALL SELF-TESTS PASSED');
  } catch (e) {
    console.error(e.message);
    process.exit(1);
  }
}
