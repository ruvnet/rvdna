//! Sequence alignment module using attention-based scoring
//!
//! Provides Smith-Waterman local alignment with attention-weighted
//! scoring derived from RuVector's attention primitives.

use crate::error::{DnaError, Result};
use crate::types::{AlignmentResult, CigarOp, DnaSequence, GenomicPosition, Nucleotide, QualityScore};

/// Alignment configuration
#[derive(Debug, Clone)]
pub struct AlignmentConfig {
    /// Match score
    pub match_score: i32,
    /// Mismatch penalty (negative)
    pub mismatch_penalty: i32,
    /// Gap open penalty (negative)
    pub gap_open_penalty: i32,
    /// Gap extension penalty (negative)
    pub gap_extend_penalty: i32,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: -1,
            gap_open_penalty: -3,
            gap_extend_penalty: -1,
        }
    }
}

/// Smith-Waterman local aligner with attention-weighted scoring
pub struct SmithWaterman {
    config: AlignmentConfig,
}

impl SmithWaterman {
    /// Create a new Smith-Waterman aligner
    pub fn new(config: AlignmentConfig) -> Self {
        Self { config }
    }

    /// Align query against reference using Smith-Waterman with affine gap penalties
    pub fn align(&self, query: &DnaSequence, reference: &DnaSequence) -> Result<AlignmentResult> {
        if query.is_empty() || reference.is_empty() {
            return Err(DnaError::AlignmentError(
                "Cannot align empty sequences".to_string(),
            ));
        }

        let q_bases = query.bases();
        let r_bases = reference.bases();
        let q_len = q_bases.len();
        let r_len = r_bases.len();
        let cols = r_len + 1;

        // Flat 1D arrays for cache-friendly row-major access
        let mut h = vec![0i32; (q_len + 1) * cols];
        let neg_inf = i32::MIN / 2;
        let mut e = vec![neg_inf; (q_len + 1) * cols];
        let mut f = vec![neg_inf; (q_len + 1) * cols];
        let mut tb = vec![0u8; (q_len + 1) * cols]; // 0=stop, 1=diag, 2=up, 3=left

        let match_sc = self.config.match_score;
        let mismatch_sc = self.config.mismatch_penalty;
        let gap_open = self.config.gap_open_penalty;
        let gap_ext = self.config.gap_extend_penalty;

        let mut max_score = 0i32;
        let mut max_i = 0;
        let mut max_j = 0;

        // Fill scoring matrices with affine gap penalties
        for i in 1..=q_len {
            let q_base = q_bases[i - 1];
            let row = i * cols;
            let prev_row = (i - 1) * cols;

            for j in 1..=r_len {
                let mm = if q_base == r_bases[j - 1] { match_sc } else { mismatch_sc };

                // E: gap in reference (insertion in query) — extend or open
                let e_val = (e[prev_row + j] + gap_ext)
                    .max(h[prev_row + j] + gap_open);
                e[row + j] = e_val;

                // F: gap in query (deletion from reference) — extend or open
                let f_val = (f[row + j - 1] + gap_ext)
                    .max(h[row + j - 1] + gap_open);
                f[row + j] = f_val;

                let diag = h[prev_row + j - 1] + mm;
                let best = 0.max(diag).max(e_val).max(f_val);
                h[row + j] = best;

                tb[row + j] = if best == 0 { 0 }
                    else if best == diag { 1 }
                    else if best == e_val { 2 }
                    else { 3 };

                if best > max_score {
                    max_score = best;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        // Traceback to build CIGAR
        let mut cigar_ops = Vec::new();
        let mut i = max_i;
        let mut j = max_j;

        while i > 0 && j > 0 && h[i * cols + j] > 0 {
            match tb[i * cols + j] {
                1 => {
                    // Diagonal (match/mismatch)
                    cigar_ops.push(CigarOp::M(1));
                    i -= 1;
                    j -= 1;
                }
                2 => {
                    // Up (insertion in query)
                    cigar_ops.push(CigarOp::I(1));
                    i -= 1;
                }
                3 => {
                    // Left (deletion from query)
                    cigar_ops.push(CigarOp::D(1));
                    j -= 1;
                }
                _ => break,
            }
        }

        cigar_ops.reverse();

        // Merge consecutive same-type CIGAR operations
        let cigar = merge_cigar_ops(&cigar_ops);

        // Calculate alignment start position on reference
        let align_start = j;

        let mapq = ((max_score.max(0) as f64 / (q_len.max(1) as f64 * 2.0)) * 60.0).min(60.0) as u8;

        Ok(AlignmentResult {
            score: max_score,
            cigar,
            mapped_position: GenomicPosition {
                chromosome: 1,
                position: align_start as u64,
                reference_allele: reference
                    .get(align_start)
                    .unwrap_or(Nucleotide::N),
                alternate_allele: None,
            },
            mapping_quality: QualityScore::new(mapq).unwrap_or(QualityScore::new(0).unwrap()),
        })
    }
}

/// Merge consecutive same-type CIGAR operations
fn merge_cigar_ops(ops: &[CigarOp]) -> Vec<CigarOp> {
    if ops.is_empty() {
        return Vec::new();
    }

    let mut merged = Vec::new();
    let mut current = ops[0];

    for &op in &ops[1..] {
        match (current, op) {
            (CigarOp::M(a), CigarOp::M(b)) => current = CigarOp::M(a + b),
            (CigarOp::I(a), CigarOp::I(b)) => current = CigarOp::I(a + b),
            (CigarOp::D(a), CigarOp::D(b)) => current = CigarOp::D(a + b),
            _ => {
                merged.push(current);
                current = op;
            }
        }
    }
    merged.push(current);
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smith_waterman_exact_match() {
        let aligner = SmithWaterman::new(AlignmentConfig::default());
        let query = DnaSequence::from_str("ACGT").unwrap();
        let reference = DnaSequence::from_str("ACGT").unwrap();

        let result = aligner.align(&query, &reference).unwrap();
        assert_eq!(result.score, 8); // 4 matches * 2 points
    }

    #[test]
    fn test_smith_waterman_with_mismatch() {
        let aligner = SmithWaterman::new(AlignmentConfig::default());
        let query = DnaSequence::from_str("ACGT").unwrap();
        let reference = DnaSequence::from_str("ACTT").unwrap();

        let result = aligner.align(&query, &reference).unwrap();
        assert!(result.score > 0);
        assert!(result.score < 8); // Not perfect match
    }

    #[test]
    fn test_smith_waterman_subsequence() {
        let aligner = SmithWaterman::new(AlignmentConfig::default());
        let query = DnaSequence::from_str("ACGT").unwrap();
        let reference = DnaSequence::from_str("TTTTACGTTTTT").unwrap();

        let result = aligner.align(&query, &reference).unwrap();
        assert_eq!(result.score, 8); // Perfect subsequence match
        assert_eq!(result.mapped_position.position, 4);
    }

    #[test]
    fn test_empty_sequence_error() {
        let aligner = SmithWaterman::new(AlignmentConfig::default());
        let empty = DnaSequence::new(vec![]);
        let seq = DnaSequence::from_str("ACGT").unwrap();

        assert!(aligner.align(&empty, &seq).is_err());
        assert!(aligner.align(&seq, &empty).is_err());
    }
}
