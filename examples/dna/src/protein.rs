//! Protein translation and amino acid analysis module
//!
//! Provides DNA to protein translation using the standard genetic code,
//! and amino acid property calculations.

use serde::{Deserialize, Serialize};

/// Amino acid representation with full names
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AminoAcid {
    /// Alanine
    Ala,
    /// Arginine
    Arg,
    /// Asparagine
    Asn,
    /// Aspartic acid
    Asp,
    /// Cysteine
    Cys,
    /// Glutamic acid
    Glu,
    /// Glutamine
    Gln,
    /// Glycine
    Gly,
    /// Histidine
    His,
    /// Isoleucine
    Ile,
    /// Leucine
    Leu,
    /// Lysine
    Lys,
    /// Methionine (start codon)
    Met,
    /// Phenylalanine
    Phe,
    /// Proline
    Pro,
    /// Serine
    Ser,
    /// Threonine
    Thr,
    /// Tryptophan
    Trp,
    /// Tyrosine
    Tyr,
    /// Valine
    Val,
    /// Stop codon
    Stop,
}

impl AminoAcid {
    /// Get single-letter code
    pub fn to_char(&self) -> char {
        match self {
            AminoAcid::Ala => 'A',
            AminoAcid::Arg => 'R',
            AminoAcid::Asn => 'N',
            AminoAcid::Asp => 'D',
            AminoAcid::Cys => 'C',
            AminoAcid::Glu => 'E',
            AminoAcid::Gln => 'Q',
            AminoAcid::Gly => 'G',
            AminoAcid::His => 'H',
            AminoAcid::Ile => 'I',
            AminoAcid::Leu => 'L',
            AminoAcid::Lys => 'K',
            AminoAcid::Met => 'M',
            AminoAcid::Phe => 'F',
            AminoAcid::Pro => 'P',
            AminoAcid::Ser => 'S',
            AminoAcid::Thr => 'T',
            AminoAcid::Trp => 'W',
            AminoAcid::Tyr => 'Y',
            AminoAcid::Val => 'V',
            AminoAcid::Stop => '*',
        }
    }

    /// Get Kyte-Doolittle hydrophobicity value
    pub fn hydrophobicity(&self) -> f32 {
        match self {
            AminoAcid::Ile => 4.5,
            AminoAcid::Val => 4.2,
            AminoAcid::Leu => 3.8,
            AminoAcid::Phe => 2.8,
            AminoAcid::Cys => 2.5,
            AminoAcid::Met => 1.9,
            AminoAcid::Ala => 1.8,
            AminoAcid::Gly => -0.4,
            AminoAcid::Thr => -0.7,
            AminoAcid::Ser => -0.8,
            AminoAcid::Trp => -0.9,
            AminoAcid::Tyr => -1.3,
            AminoAcid::Pro => -1.6,
            AminoAcid::His => -3.2,
            AminoAcid::Glu => -3.5,
            AminoAcid::Gln => -3.5,
            AminoAcid::Asp => -3.5,
            AminoAcid::Asn => -3.5,
            AminoAcid::Lys => -3.9,
            AminoAcid::Arg => -4.5,
            AminoAcid::Stop => 0.0,
        }
    }
}

/// Translate a DNA sequence to a vector of amino acids using the standard genetic code.
///
/// Translation proceeds in triplets (codons) from the start of the sequence.
/// Stop codons (TAA, TAG, TGA) terminate translation.
/// Incomplete codons at the end are ignored.
pub fn translate_dna(dna: &[u8]) -> Vec<AminoAcid> {
    let mut proteins = Vec::new();

    for chunk in dna.chunks(3) {
        if chunk.len() < 3 {
            break;
        }

        let codon = [
            chunk[0].to_ascii_uppercase(),
            chunk[1].to_ascii_uppercase(),
            chunk[2].to_ascii_uppercase(),
        ];

        let aa = match &codon {
            b"ATG" => AminoAcid::Met,
            b"TGG" => AminoAcid::Trp,
            b"TTT" | b"TTC" => AminoAcid::Phe,
            b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => AminoAcid::Leu,
            b"ATT" | b"ATC" | b"ATA" => AminoAcid::Ile,
            b"GTT" | b"GTC" | b"GTA" | b"GTG" => AminoAcid::Val,
            b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => AminoAcid::Ser,
            b"CCT" | b"CCC" | b"CCA" | b"CCG" => AminoAcid::Pro,
            b"ACT" | b"ACC" | b"ACA" | b"ACG" => AminoAcid::Thr,
            b"GCT" | b"GCC" | b"GCA" | b"GCG" => AminoAcid::Ala,
            b"TAT" | b"TAC" => AminoAcid::Tyr,
            b"CAT" | b"CAC" => AminoAcid::His,
            b"CAA" | b"CAG" => AminoAcid::Gln,
            b"AAT" | b"AAC" => AminoAcid::Asn,
            b"AAA" | b"AAG" => AminoAcid::Lys,
            b"GAT" | b"GAC" => AminoAcid::Asp,
            b"GAA" | b"GAG" => AminoAcid::Glu,
            b"TGT" | b"TGC" => AminoAcid::Cys,
            b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => AminoAcid::Arg,
            b"GGT" | b"GGC" | b"GGA" | b"GGG" => AminoAcid::Gly,
            b"TAA" | b"TAG" | b"TGA" => break, // Stop codons
            _ => continue, // Unknown codon, skip
        };

        proteins.push(aa);
    }

    proteins
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_basic() {
        let dna = b"ATGGCAGGT";
        let result = translate_dna(dna);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], AminoAcid::Met);
        assert_eq!(result[1], AminoAcid::Ala);
        assert_eq!(result[2], AminoAcid::Gly);
    }

    #[test]
    fn test_translate_stop_codon() {
        let dna = b"ATGGCATAA"; // Met-Ala-Stop
        let result = translate_dna(dna);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_hydrophobicity() {
        assert_eq!(AminoAcid::Ile.hydrophobicity(), 4.5);
        assert_eq!(AminoAcid::Arg.hydrophobicity(), -4.5);
    }
}
