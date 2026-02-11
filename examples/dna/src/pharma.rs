//! Pharmacogenomics module
//!
//! Provides CYP enzyme star allele calling and metabolizer phenotype
//! prediction for pharmacogenomic analysis.

use serde::{Deserialize, Serialize};

/// CYP2D6 star allele classification
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum StarAllele {
    /// *1 - Normal function (wild-type)
    Star1,
    /// *2 - Normal function
    Star2,
    /// *3 - No function (frameshift)
    Star3,
    /// *4 - No function (splicing defect)
    Star4,
    /// *5 - No function (gene deletion)
    Star5,
    /// *6 - No function (frameshift)
    Star6,
    /// *10 - Decreased function
    Star10,
    /// *17 - Decreased function
    Star17,
    /// *41 - Decreased function
    Star41,
    /// Unknown allele
    Unknown,
}

impl StarAllele {
    /// Get the activity score for this allele
    pub fn activity_score(&self) -> f64 {
        match self {
            StarAllele::Star1 | StarAllele::Star2 => 1.0,
            StarAllele::Star10 | StarAllele::Star17 | StarAllele::Star41 => 0.5,
            StarAllele::Star3
            | StarAllele::Star4
            | StarAllele::Star5
            | StarAllele::Star6 => 0.0,
            StarAllele::Unknown => 0.5,
        }
    }
}

/// Drug metabolizer phenotype
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum MetabolizerPhenotype {
    /// Ultra-rapid metabolizer (activity score > 2.0)
    UltraRapid,
    /// Normal metabolizer (1.0 <= activity score <= 2.0)
    Normal,
    /// Intermediate metabolizer (0.5 <= activity score < 1.0)
    Intermediate,
    /// Poor metabolizer (activity score < 0.5)
    Poor,
}

/// Pharmacogenomic variant for a specific gene
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PharmaVariant {
    /// Gene name (e.g., "CYP2D6")
    pub gene: String,
    /// Genomic position
    pub position: u64,
    /// Reference allele
    pub ref_allele: u8,
    /// Alternate allele
    pub alt_allele: u8,
    /// Clinical significance
    pub significance: String,
}

/// Call CYP2D6 star allele from observed variants
///
/// Uses a simplified lookup table based on key defining variants.
pub fn call_star_allele(variants: &[(u64, u8, u8)]) -> StarAllele {
    for &(pos, ref_allele, alt_allele) in variants {
        match (pos, ref_allele, alt_allele) {
            // *4: G>A at intron 3/exon 4 boundary (rs3892097)
            (42130692, b'G', b'A') => return StarAllele::Star4,
            // *5: whole gene deletion
            (42126611, b'T', b'-') => return StarAllele::Star5,
            // *3: frameshift (A deletion at rs35742686)
            (42127941, b'A', b'-') => return StarAllele::Star3,
            // *6: T deletion at rs5030655
            (42127803, b'T', b'-') => return StarAllele::Star6,
            // *10: C>T at rs1065852
            (42126938, b'C', b'T') => return StarAllele::Star10,
            _ => {}
        }
    }

    StarAllele::Star1 // Wild-type
}

/// Predict metabolizer phenotype from diplotype (two alleles)
pub fn predict_phenotype(allele1: &StarAllele, allele2: &StarAllele) -> MetabolizerPhenotype {
    let total_activity = allele1.activity_score() + allele2.activity_score();

    if total_activity > 2.0 {
        MetabolizerPhenotype::UltraRapid
    } else if total_activity >= 1.0 {
        MetabolizerPhenotype::Normal
    } else if total_activity >= 0.5 {
        MetabolizerPhenotype::Intermediate
    } else {
        MetabolizerPhenotype::Poor
    }
}

/// Drug recommendation based on metabolizer phenotype
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DrugRecommendation {
    /// Drug name
    pub drug: String,
    /// Gene involved
    pub gene: String,
    /// Recommendation text
    pub recommendation: String,
    /// Dosing adjustment factor (1.0 = standard dose)
    pub dose_factor: f64,
}

/// Get drug recommendations for a given phenotype
pub fn get_recommendations(
    gene: &str,
    phenotype: &MetabolizerPhenotype,
) -> Vec<DrugRecommendation> {
    match (gene, phenotype) {
        ("CYP2D6", MetabolizerPhenotype::Poor) => vec![
            DrugRecommendation {
                drug: "Codeine".to_string(),
                gene: gene.to_string(),
                recommendation: "Avoid codeine; use alternative analgesic".to_string(),
                dose_factor: 0.0,
            },
            DrugRecommendation {
                drug: "Tamoxifen".to_string(),
                gene: gene.to_string(),
                recommendation: "Consider alternative endocrine therapy".to_string(),
                dose_factor: 0.0,
            },
        ],
        ("CYP2D6", MetabolizerPhenotype::UltraRapid) => vec![DrugRecommendation {
            drug: "Codeine".to_string(),
            gene: gene.to_string(),
            recommendation: "Avoid codeine; risk of toxicity from rapid conversion".to_string(),
            dose_factor: 0.0,
        }],
        ("CYP2D6", MetabolizerPhenotype::Intermediate) => vec![DrugRecommendation {
            drug: "Codeine".to_string(),
            gene: gene.to_string(),
            recommendation: "Use lower dose or alternative".to_string(),
            dose_factor: 0.5,
        }],
        _ => vec![DrugRecommendation {
            drug: "Standard".to_string(),
            gene: gene.to_string(),
            recommendation: "Use standard dosing".to_string(),
            dose_factor: 1.0,
        }],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_star_allele_calling() {
        // Wild-type
        assert_eq!(call_star_allele(&[]), StarAllele::Star1);

        // *4 variant
        let star4 = call_star_allele(&[(42130692, b'G', b'A')]);
        assert_eq!(star4, StarAllele::Star4);
        assert_eq!(star4.activity_score(), 0.0);

        // *10 variant (decreased function)
        let star10 = call_star_allele(&[(42126938, b'C', b'T')]);
        assert_eq!(star10, StarAllele::Star10);
        assert_eq!(star10.activity_score(), 0.5);
    }

    #[test]
    fn test_phenotype_prediction() {
        assert_eq!(
            predict_phenotype(&StarAllele::Star1, &StarAllele::Star1),
            MetabolizerPhenotype::Normal
        );
        assert_eq!(
            predict_phenotype(&StarAllele::Star1, &StarAllele::Star4),
            MetabolizerPhenotype::Normal
        );
        assert_eq!(
            predict_phenotype(&StarAllele::Star4, &StarAllele::Star10),
            MetabolizerPhenotype::Intermediate
        );
        assert_eq!(
            predict_phenotype(&StarAllele::Star4, &StarAllele::Star4),
            MetabolizerPhenotype::Poor
        );
    }

    #[test]
    fn test_drug_recommendations() {
        let recs = get_recommendations("CYP2D6", &MetabolizerPhenotype::Poor);
        assert!(recs.len() >= 1);
        assert_eq!(recs[0].dose_factor, 0.0);

        let recs_normal = get_recommendations("CYP2D6", &MetabolizerPhenotype::Normal);
        assert_eq!(recs_normal[0].dose_factor, 1.0);
    }
}
