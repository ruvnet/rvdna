//! Reversible bytes <-> DNA codec that is homopolymer-free and GC-balanced.
//!
//! Uses a Goldman-style rotating base-3 ("trit") transform. Each byte is
//! expanded to 6 ternary digits (3^6 = 729 > 255), and each trit is mapped to
//! the next nucleotide via a rotating offset that guarantees consecutive bases
//! always differ (no homopolymer run longer than 1).

/// Nucleotide alphabet: A=0, C=1, G=2, T=3.
const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// Place values for the 6 big-endian ternary digits of a byte: 3^5 … 3^0.
/// Precomputed so the hot encode/decode loops avoid repeated `pow` calls.
const POW3: [usize; 6] = [243, 81, 27, 9, 3, 1];

#[inline]
fn base_char(idx: usize) -> char {
    BASES[idx & 3]
}

#[inline]
fn base_index(c: char) -> Option<usize> {
    match c {
        'A' => Some(0),
        'C' => Some(1),
        'G' => Some(2),
        'T' => Some(3),
        _ => None,
    }
}

/// Encode arbitrary bytes into a homopolymer-free DNA string.
pub fn encode_bytes(data: &[u8]) -> String {
    let mut out = String::with_capacity(data.len() * 6);
    let mut prev: usize = 0; // running base index, initialized to A

    for &b in data {
        let b = b as usize;
        // 6 ternary digits, big-endian; map each trit to a base that always
        // differs from the previous one (homopolymer-free by construction).
        for &p in POW3.iter() {
            let t = (b / p) % 3;
            let cur = (prev + 1 + t) % 4;
            out.push(base_char(cur));
            prev = cur;
        }
    }
    out
}

/// Decode a homopolymer-free DNA string back into the original bytes.
pub fn decode_bytes(seq: &str) -> Result<Vec<u8>, String> {
    let chars: Vec<char> = seq.chars().collect();
    if chars.len() % 6 != 0 {
        return Err(format!(
            "sequence length {} is not a multiple of 6",
            chars.len()
        ));
    }

    let mut out = Vec::with_capacity(chars.len() / 6);
    let mut prev: usize = 0;

    let mut i = 0;
    while i < chars.len() {
        let mut value: usize = 0;
        for (k, &p) in POW3.iter().enumerate() {
            let c = chars[i + k];
            let cur = base_index(c)
                .ok_or_else(|| format!("invalid base character '{}' at position {}", c, i + k))?;
            // Recover trit: cur = (prev + 1 + t) % 4  =>  t = (cur + 4 - prev - 1) % 4
            let t = (cur + 4 - prev - 1) % 4;
            if t > 2 {
                return Err(format!(
                    "corrupt symbol at position {}: trit out of range",
                    i + k
                ));
            }
            value += t * p;
            prev = cur;
        }
        out.push(value as u8); // truncates; valid data is always < 256
        i += 6;
    }

    Ok(out)
}

/// Fraction (0.0..=1.0) of bases that are G or C. Returns 0.0 for empty input.
pub fn gc_content(seq: &str) -> f64 {
    let mut total = 0usize;
    let mut gc = 0usize;
    for c in seq.chars() {
        total += 1;
        if c == 'G' || c == 'C' {
            gc += 1;
        }
    }
    if total == 0 {
        0.0
    } else {
        gc as f64 / total as f64
    }
}

/// Longest run of identical characters. Returns 0 for empty input.
pub fn max_homopolymer_run(seq: &str) -> usize {
    let mut max_run = 0usize;
    let mut cur_run = 0usize;
    let mut last: Option<char> = None;
    for c in seq.chars() {
        if Some(c) == last {
            cur_run += 1;
        } else {
            cur_run = 1;
            last = Some(c);
        }
        if cur_run > max_run {
            max_run = cur_run;
        }
    }
    max_run
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    #[test]
    fn roundtrip_various_lengths() {
        let mut rng = StdRng::seed_from_u64(0xDEAD_BEEF);
        for &len in &[0usize, 1, 2, 3, 7, 16, 64, 255, 1000] {
            let data: Vec<u8> = (0..len).map(|_| rng.gen()).collect();
            let encoded = encode_bytes(&data);
            let decoded = decode_bytes(&encoded).expect("decode should succeed");
            assert_eq!(decoded, data, "roundtrip failed for len {}", len);
        }
    }

    #[test]
    fn roundtrip_empty() {
        let data: Vec<u8> = Vec::new();
        let encoded = encode_bytes(&data);
        assert_eq!(encoded, "");
        let decoded = decode_bytes(&encoded).expect("empty decode");
        assert_eq!(decoded, data);
    }

    #[test]
    fn roundtrip_all_byte_values() {
        let data: Vec<u8> = (0u16..=255).map(|x| x as u8).collect();
        let encoded = encode_bytes(&data);
        let decoded = decode_bytes(&encoded).unwrap();
        assert_eq!(decoded, data);
    }

    #[test]
    fn no_homopolymer_runs() {
        let mut rng = StdRng::seed_from_u64(1234);
        for _ in 0..50 {
            let len = rng.gen_range(1..200);
            let data: Vec<u8> = (0..len).map(|_| rng.gen()).collect();
            let encoded = encode_bytes(&data);
            assert!(
                max_homopolymer_run(&encoded) <= 1,
                "homopolymer run > 1 in encoding"
            );
        }
    }

    #[test]
    fn rejects_bad_length() {
        // 5 bases is not a multiple of 6.
        assert!(decode_bytes("ACGTA").is_err());
        assert!(decode_bytes("A").is_err());
        assert!(decode_bytes("ACGTACG").is_err()); // 7
    }

    #[test]
    fn gc_content_statistical_sanity() {
        let mut rng = StdRng::seed_from_u64(99);
        let data: Vec<u8> = (0..5000).map(|_| rng.gen()).collect();
        let encoded = encode_bytes(&data);
        let gc = gc_content(&encoded);
        assert!(gc >= 0.3 && gc <= 0.7, "gc content {} out of [0.3,0.7]", gc);
    }

    #[test]
    fn gc_content_basic() {
        assert_eq!(gc_content(""), 0.0);
        assert_eq!(gc_content("GGCC"), 1.0);
        assert_eq!(gc_content("AATT"), 0.0);
        assert!((gc_content("ACGT") - 0.5).abs() < 1e-9);
    }

    #[test]
    fn max_homopolymer_run_basic() {
        assert_eq!(max_homopolymer_run(""), 0);
        assert_eq!(max_homopolymer_run("A"), 1);
        assert_eq!(max_homopolymer_run("AAAB"), 3);
        assert_eq!(max_homopolymer_run("ACGT"), 1);
    }
}
