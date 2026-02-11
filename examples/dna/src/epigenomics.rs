//! Epigenomics analysis module
//!
//! Provides methylation profiling and epigenetic age prediction
//! using the Horvath clock model.

use serde::{Deserialize, Serialize};

/// A CpG site with methylation data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CpGSite {
    /// Chromosome number
    pub chromosome: u8,
    /// Genomic position
    pub position: u64,
    /// Methylation level (beta value, 0.0 to 1.0)
    pub methylation_level: f32,
}

/// Methylation profile containing CpG site measurements
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MethylationProfile {
    /// CpG sites with measured methylation levels
    pub sites: Vec<CpGSite>,
}

impl MethylationProfile {
    /// Create a methylation profile from position and beta value arrays
    pub fn from_beta_values(positions: Vec<(u8, u64)>, betas: Vec<f32>) -> Self {
        let sites = positions
            .into_iter()
            .zip(betas.into_iter())
            .map(|((chr, pos), beta)| CpGSite {
                chromosome: chr,
                position: pos,
                methylation_level: beta.clamp(0.0, 1.0),
            })
            .collect();

        Self { sites }
    }

    /// Calculate mean methylation across all sites
    pub fn mean_methylation(&self) -> f32 {
        if self.sites.is_empty() {
            return 0.0;
        }
        let sum: f32 = self.sites.iter().map(|s| s.methylation_level).sum();
        sum / self.sites.len() as f32
    }
}

/// Horvath epigenetic clock for biological age prediction
///
/// Uses a simplified linear model based on CpG site methylation levels
/// to predict biological age.
pub struct HorvathClock {
    /// Intercept term
    intercept: f64,
    /// Coefficient per CpG site bin
    coefficients: Vec<f64>,
    /// Number of bins to partition sites into
    num_bins: usize,
}

impl HorvathClock {
    /// Create the default Horvath clock model
    ///
    /// Uses a simplified model with binned methylation values.
    /// Real implementation would use 353 specific CpG sites.
    pub fn default_clock() -> Self {
        Self {
            intercept: 30.0,
            coefficients: vec![
                -15.0, // Low methylation bin (young)
                10.0,  // High methylation bin (age-associated)
                0.5,   // Neutral bin
            ],
            num_bins: 3,
        }
    }

    /// Predict biological age from a methylation profile
    pub fn predict_age(&self, profile: &MethylationProfile) -> f64 {
        if profile.sites.is_empty() {
            return self.intercept;
        }

        // Partition sites into bins and compute mean methylation per bin
        let bin_size = profile.sites.len() / self.num_bins.max(1);
        let mut age = self.intercept;

        for (bin_idx, coefficient) in self.coefficients.iter().enumerate() {
            let start = bin_idx * bin_size;
            let end = ((bin_idx + 1) * bin_size).min(profile.sites.len());

            if start >= profile.sites.len() {
                break;
            }

            let bin_sites = &profile.sites[start..end];
            if !bin_sites.is_empty() {
                let mean_meth: f64 = bin_sites
                    .iter()
                    .map(|s| s.methylation_level as f64)
                    .sum::<f64>()
                    / bin_sites.len() as f64;

                age += coefficient * mean_meth;
            }
        }

        age.max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_methylation_profile() {
        let positions = vec![(1, 1000), (1, 2000)];
        let betas = vec![0.3, 0.7];
        let profile = MethylationProfile::from_beta_values(positions, betas);

        assert_eq!(profile.sites.len(), 2);
        assert!((profile.mean_methylation() - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_horvath_clock() {
        let clock = HorvathClock::default_clock();
        let positions = vec![(1, 1000), (1, 2000), (1, 3000)];
        let betas = vec![0.5, 0.5, 0.5];
        let profile = MethylationProfile::from_beta_values(positions, betas);
        let age = clock.predict_age(&profile);
        assert!(age > 0.0);
    }
}
