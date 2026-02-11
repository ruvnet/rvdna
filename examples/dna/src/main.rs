//! DNA Analyzer Demo - RuVector Genomic Analysis Pipeline
//!
//! Demonstrates SOTA genomic analysis using:
//! - HNSW k-mer indexing for fast sequence search
//! - Attention-based sequence alignment
//! - Variant calling from pileup data
//! - Protein translation and contact prediction
//! - Epigenetic age prediction (Horvath clock)
//! - Pharmacogenomic star allele calling

use dna_analyzer_example::prelude::*;
use dna_analyzer_example::{
    alignment::{AlignmentConfig, SmithWaterman},
    epigenomics::{HorvathClock, MethylationProfile},
    pharma,
    protein::translate_dna,
    variant::{PileupColumn, VariantCaller, VariantCallerConfig},
};
use rand::Rng;
use tracing::{info, Level};
use tracing_subscriber::FmtSubscriber;

fn main() -> anyhow::Result<()> {
    // Initialize tracing
    let subscriber = FmtSubscriber::builder()
        .with_max_level(Level::INFO)
        .finish();
    tracing::subscriber::set_global_default(subscriber)?;

    info!("RuVector DNA Analyzer - Genomic Analysis Pipeline");
    info!("================================================");

    // -----------------------------------------------------------------------
    // Stage 1: Generate synthetic DNA sequence (10,000 bp)
    // -----------------------------------------------------------------------
    info!("\nStage 1: Generating synthetic DNA sequence");
    let sequence = generate_synthetic_dna(10_000);
    let reference = generate_synthetic_dna(10_000);

    info!("  Sequence length: {} bp", sequence.len());
    info!(
        "  GC content: {:.2}%",
        calculate_gc_content(&sequence) * 100.0
    );
    info!("  First 60bp: {}", &sequence.to_string()[..60]);

    // -----------------------------------------------------------------------
    // Stage 2: K-mer encoding and HNSW similarity search
    // -----------------------------------------------------------------------
    info!("\nStage 2: K-mer indexing and similarity search");
    let kmer_start = std::time::Instant::now();

    let query_vec = sequence.to_kmer_vector(11, 512)?;
    let ref_vec = reference.to_kmer_vector(11, 512)?;

    // Compute cosine similarity
    let similarity: f32 = query_vec
        .iter()
        .zip(ref_vec.iter())
        .map(|(a, b)| a * b)
        .sum();

    info!("  K-mer vector dimensions: {}", query_vec.len());
    info!("  Cosine similarity to reference: {:.4}", similarity);
    info!("  K-mer encoding time: {:?}", kmer_start.elapsed());

    // -----------------------------------------------------------------------
    // Stage 3: Sequence alignment
    // -----------------------------------------------------------------------
    info!("\nStage 3: Sequence alignment (Smith-Waterman)");
    let align_start = std::time::Instant::now();

    let query_fragment = DnaSequence::from_str(&sequence.to_string()[100..200])?;
    let aligner = SmithWaterman::new(AlignmentConfig::default());
    let alignment = aligner.align(&query_fragment, &sequence)?;

    info!("  Alignment score: {}", alignment.score);
    info!(
        "  Mapped position: {}",
        alignment.mapped_position.position
    );
    info!(
        "  Mapping quality: {}",
        alignment.mapping_quality.value()
    );
    info!("  CIGAR ops: {}", alignment.cigar.len());
    info!("  Alignment time: {:?}", align_start.elapsed());

    // -----------------------------------------------------------------------
    // Stage 4: Variant calling
    // -----------------------------------------------------------------------
    info!("\nStage 4: Variant calling");
    let variant_start = std::time::Instant::now();

    let caller = VariantCaller::new(VariantCallerConfig::default());
    let mut variant_count = 0;
    let seq_bytes = sequence.to_string().into_bytes();
    let ref_bytes = reference.to_string().into_bytes();

    for i in 0..seq_bytes.len().min(ref_bytes.len()).min(1000) {
        let mut rng = rand::thread_rng();
        let depth = rng.gen_range(10..31);

        let bases: Vec<u8> = (0..depth)
            .map(|_| {
                if rng.gen::<f32>() < 0.95 {
                    seq_bytes[i]
                } else {
                    [b'A', b'C', b'G', b'T'][rng.gen_range(0..4)]
                }
            })
            .collect();
        let qualities: Vec<u8> = (0..depth).map(|_| rng.gen_range(20..41)).collect();

        let pileup = PileupColumn {
            bases,
            qualities,
            position: i as u64,
            chromosome: 1,
        };

        if caller.call_snp(&pileup, ref_bytes[i]).is_some() {
            variant_count += 1;
        }
    }

    info!("  Positions analyzed: 1000");
    info!("  Variants detected: {}", variant_count);
    info!("  Variant calling time: {:?}", variant_start.elapsed());

    // -----------------------------------------------------------------------
    // Stage 5: Protein translation and structure prediction
    // -----------------------------------------------------------------------
    info!("\nStage 5: Protein translation and contact prediction");
    let protein_start = std::time::Instant::now();

    let seq_str = sequence.to_string();
    let coding_start = seq_str.find("ATG").unwrap_or(0);
    let coding_region = &seq_str.as_bytes()[coding_start..(coding_start + 300).min(seq_str.len())];
    let amino_acids = translate_dna(coding_region);

    info!("  Protein length: {} amino acids", amino_acids.len());

    if amino_acids.len() >= 10 {
        let protein_str: String = amino_acids.iter().map(|aa| aa.to_char()).collect();
        let preview = if protein_str.len() > 50 {
            format!("{}...", &protein_str[..50])
        } else {
            protein_str
        };
        info!("  Protein sequence: {}", preview);

        // Build contact graph
        let residues: Vec<ProteinResidue> = amino_acids
            .iter()
            .map(|aa| match aa.to_char() {
                'A' => ProteinResidue::A,
                'R' => ProteinResidue::R,
                'N' => ProteinResidue::N,
                'D' => ProteinResidue::D,
                'C' => ProteinResidue::C,
                'E' => ProteinResidue::E,
                'Q' => ProteinResidue::Q,
                'G' => ProteinResidue::G,
                'H' => ProteinResidue::H,
                'I' => ProteinResidue::I,
                'L' => ProteinResidue::L,
                'K' => ProteinResidue::K,
                'M' => ProteinResidue::M,
                'F' => ProteinResidue::F,
                'P' => ProteinResidue::P,
                'S' => ProteinResidue::S,
                'T' => ProteinResidue::T,
                'W' => ProteinResidue::W,
                'Y' => ProteinResidue::Y,
                'V' => ProteinResidue::V,
                _ => ProteinResidue::X,
            })
            .collect();
        let protein_seq = ProteinSequence::new(residues);
        let graph = protein_seq.build_contact_graph(8.0)?;
        let contacts = protein_seq.predict_contacts(&graph)?;

        info!("  Contact graph edges: {}", graph.edges.len());
        info!("  Top predicted contacts:");
        for (i, (r1, r2, score)) in contacts.iter().take(5).enumerate() {
            info!(
                "    {}. Residues {} <-> {} (score: {:.3})",
                i + 1,
                r1,
                r2,
                score
            );
        }
    }
    info!("  Protein analysis time: {:?}", protein_start.elapsed());

    // -----------------------------------------------------------------------
    // Stage 6: Epigenetic age prediction
    // -----------------------------------------------------------------------
    info!("\nStage 6: Epigenetic age prediction (Horvath clock)");
    let epi_start = std::time::Instant::now();

    let mut rng = rand::thread_rng();
    let positions: Vec<(u8, u64)> = (0..500).map(|i| (1, i * 1000)).collect();
    let betas: Vec<f32> = (0..500).map(|_| rng.gen_range(0.1..0.9)).collect();

    let profile = MethylationProfile::from_beta_values(positions, betas);
    let clock = HorvathClock::default_clock();
    let predicted_age = clock.predict_age(&profile);

    info!("  CpG sites analyzed: {}", profile.sites.len());
    info!("  Mean methylation: {:.3}", profile.mean_methylation());
    info!("  Predicted biological age: {:.1} years", predicted_age);
    info!("  Epigenomics time: {:?}", epi_start.elapsed());

    // -----------------------------------------------------------------------
    // Stage 7: Pharmacogenomics
    // -----------------------------------------------------------------------
    info!("\nStage 7: Pharmacogenomic analysis");

    let cyp2d6_variants = vec![(42130692, b'G', b'A')];
    let allele1 = pharma::call_star_allele(&cyp2d6_variants);
    let allele2 = pharma::StarAllele::Star1;
    let phenotype = pharma::predict_phenotype(&allele1, &allele2);

    info!(
        "  CYP2D6 allele 1: {:?} (activity: {})",
        allele1,
        allele1.activity_score()
    );
    info!(
        "  CYP2D6 allele 2: {:?} (activity: {})",
        allele2,
        allele2.activity_score()
    );
    info!("  Metabolizer phenotype: {:?}", phenotype);

    let recommendations = pharma::get_recommendations("CYP2D6", &phenotype);
    if !recommendations.is_empty() {
        info!("  Drug recommendations:");
        for rec in &recommendations {
            info!(
                "    - {}: {} (dose factor: {:.1}x)",
                rec.drug, rec.recommendation, rec.dose_factor
            );
        }
    }

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------
    info!("\nPipeline Summary");
    info!("==================");
    info!("  Sequence length: {} bp", sequence.len());
    info!("  Variants called: {}", variant_count);
    info!("  Protein residues: {}", amino_acids.len());
    info!("  Predicted age: {:.1} years", predicted_age);
    info!("  Phenotype: {:?}", phenotype);

    info!("\nAnalysis complete!");

    Ok(())
}

/// Generate synthetic DNA sequence with realistic characteristics
fn generate_synthetic_dna(length: usize) -> DnaSequence {
    let mut rng = rand::thread_rng();
    let bases = [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
    let sequence: Vec<Nucleotide> = (0..length)
        .map(|_| bases[rng.gen_range(0..4)])
        .collect();
    DnaSequence::new(sequence)
}

/// Calculate GC content of DNA sequence
fn calculate_gc_content(sequence: &DnaSequence) -> f64 {
    let gc_count = sequence
        .bases()
        .iter()
        .filter(|&&b| b == Nucleotide::G || b == Nucleotide::C)
        .count();

    gc_count as f64 / sequence.len() as f64
}
