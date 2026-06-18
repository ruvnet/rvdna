//! Read clustering and consensus calling for DNA storage decoding.

use std::collections::HashMap;

/// Normalized distance between two reads: (Hamming over shared min-length prefix
/// + absolute length difference) / max(len_a, len_b).
fn normalized_distance(a: &str, b: &str) -> f64 {
    let ab: Vec<char> = a.chars().collect();
    let bb: Vec<char> = b.chars().collect();
    let la = ab.len();
    let lb = bb.len();
    let max_len = la.max(lb);
    if max_len == 0 {
        return 0.0;
    }
    let min_len = la.min(lb);
    let mut hamming = 0usize;
    for i in 0..min_len {
        if ab[i] != bb[i] {
            hamming += 1;
        }
    }
    let len_diff = if la > lb { la - lb } else { lb - la };
    (hamming + len_diff) as f64 / max_len as f64
}

/// Greedy single-pass clustering. Each read is compared against the first
/// member (representative) of every existing cluster; it joins the first
/// cluster within `max_dist`, otherwise it forms a new cluster.
pub fn cluster(reads: &[String], max_dist: f64) -> Vec<Vec<usize>> {
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for (idx, read) in reads.iter().enumerate() {
        let mut placed = false;
        for c in clusters.iter_mut() {
            let rep = &reads[c[0]];
            if normalized_distance(read, rep) <= max_dist {
                c.push(idx);
                placed = true;
                break;
            }
        }
        if !placed {
            clusters.push(vec![idx]);
        }
    }
    clusters
}

/// Position-wise majority-vote consensus over a group of reads.
///
/// The consensus length L is the modal (most common) read length in the group.
/// For each position 0..L the majority base among reads having a character at
/// that position is chosen; ties are broken deterministically (first-seen, then
/// lexicographic). Empty group => "".
pub fn consensus(reads: &[&str]) -> String {
    if reads.is_empty() {
        return String::new();
    }

    // Determine modal length (ties broken by smaller length for determinism).
    let mut len_counts: HashMap<usize, usize> = HashMap::new();
    for r in reads {
        *len_counts.entry(r.chars().count()).or_insert(0) += 1;
    }
    let mut modal_len = 0usize;
    let mut best_count = 0usize;
    for (&len, &count) in len_counts.iter() {
        if count > best_count || (count == best_count && len < modal_len) {
            best_count = count;
            modal_len = len;
        }
    }

    if modal_len == 0 {
        return String::new();
    }

    // Pre-collect chars per read for indexed access.
    let read_chars: Vec<Vec<char>> = reads.iter().map(|r| r.chars().collect()).collect();

    let mut result = String::with_capacity(modal_len);
    for pos in 0..modal_len {
        // Count bases at this position, tracking first-seen order for tie-break.
        let mut counts: HashMap<char, usize> = HashMap::new();
        let mut order: Vec<char> = Vec::new();
        for rc in &read_chars {
            if pos < rc.len() {
                let ch = rc[pos];
                let e = counts.entry(ch).or_insert(0);
                if *e == 0 {
                    order.push(ch);
                }
                *e += 1;
            }
        }

        if order.is_empty() {
            continue;
        }

        // Pick majority; tie-break by first-seen order, then lexicographic.
        let mut best_char = order[0];
        let mut best = counts[&best_char];
        for &ch in order.iter().skip(1) {
            let c = counts[&ch];
            if c > best || (c == best && ch < best_char) {
                best = c;
                best_char = ch;
            }
        }
        result.push(best_char);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consensus_recovers_original_from_noisy_copies() {
        let original = "ACGTACGTACGTACGT";
        // Several noisy copies, each with a few flips at different spots.
        let copies = vec![
            "ACGTACGTACGTACGT".to_string(), // clean
            "TCGTACGTACGTACGT".to_string(), // pos 0 flipped
            "ACGTACGTACGAACGT".to_string(), // pos 11 flipped
            "ACGTACCTACGTACGT".to_string(), // pos 6 flipped
            "ACGTACGTACGTACGT".to_string(), // clean
        ];
        let refs: Vec<&str> = copies.iter().map(|s| s.as_str()).collect();
        let got = consensus(&refs);
        assert_eq!(got, original, "consensus failed to recover original");
    }

    #[test]
    fn consensus_empty() {
        let empty: Vec<&str> = Vec::new();
        assert_eq!(consensus(&empty), "");
    }

    #[test]
    fn consensus_handles_varying_lengths() {
        // Modal length is 4 (three reads of length 4).
        let reads = vec!["ACGT", "ACGT", "ACGT", "ACG", "ACGTACGT"];
        let got = consensus(&reads);
        assert_eq!(got, "ACGT");
    }

    #[test]
    fn cluster_groups_and_separates() {
        let reads = vec![
            "AAAAAAAAAA".to_string(),
            "AAAAAAAAAT".to_string(), // 1 diff from first
            "AAAAAAAATT".to_string(), // 2 diff from first
            "GGGGGGGGGG".to_string(), // clearly different
            "GGGGGGGGGC".to_string(), // near the G read
        ];
        let clusters = cluster(&reads, 0.3);

        // Find which cluster each index landed in.
        let find = |idx: usize| {
            clusters
                .iter()
                .position(|c| c.contains(&idx))
                .expect("index must be clustered")
        };

        // 0,1,2 together; 3,4 together; the two groups distinct.
        assert_eq!(find(0), find(1));
        assert_eq!(find(0), find(2));
        assert_eq!(find(3), find(4));
        assert_ne!(find(0), find(3));
    }

    #[test]
    fn cluster_empty() {
        let reads: Vec<String> = Vec::new();
        assert!(cluster(&reads, 0.5).is_empty());
    }
}
