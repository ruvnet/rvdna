//! # DNA Analyzer - State-of-the-Art Genomic Analysis with RuVector
//!
//! A comprehensive genomic analysis pipeline demonstrating RuVector's
//! vector computing capabilities applied to bioinformatics:
//!
//! - **K-mer HNSW Indexing**: Fast sequence similarity search via vector embeddings
//! - **Attention Alignment**: Smith-Waterman alignment with attention scoring
//! - **Variant Calling**: SNP/indel detection from pileup data
//! - **Protein Translation**: DNA-to-protein with contact graph prediction
//! - **Epigenomics**: Methylation profiling and Horvath clock age prediction
//! - **Pharmacogenomics**: CYP enzyme star allele calling and drug recommendations
//! - **Pipeline**: DAG-based orchestration of analysis stages
//! - **RVDNA Format**: AI-native binary file format with pre-computed tensors

#![warn(missing_docs)]
#![allow(clippy::all)]

pub mod error;
pub mod types;
pub mod kmer;
pub mod alignment;
pub mod variant;
pub mod protein;
pub mod epigenomics;
pub mod pharma;
pub mod pipeline;
pub mod rvdna;
pub mod real_data;

pub use error::{DnaError, Result};
pub use types::{
    AlignmentResult, AnalysisConfig, CigarOp, ContactGraph, DnaSequence, GenomicPosition,
    KmerIndex, Nucleotide, ProteinResidue, ProteinSequence, QualityScore, Variant,
};
pub use variant::{
    FilterStatus, Genotype, PileupColumn, VariantCall, VariantCaller, VariantCallerConfig,
};
pub use protein::{AminoAcid, translate_dna};
pub use epigenomics::{CpGSite, HorvathClock, MethylationProfile};
pub use alignment::{AlignmentConfig, SmithWaterman};
pub use pharma::{
    call_star_allele, get_recommendations, predict_phenotype, DrugRecommendation,
    MetabolizerPhenotype, PharmaVariant, StarAllele,
};
pub use rvdna::{
    Codec, RvdnaHeader, RvdnaReader, RvdnaWriter, RvdnaStats,
    SparseAttention, VariantTensor, KmerVectorBlock,
    encode_2bit, decode_2bit, fasta_to_rvdna,
};

pub use ruvector_core::{
    types::{DbOptions, DistanceMetric, HnswConfig, SearchQuery, SearchResult, VectorEntry},
    VectorDB,
};

/// Prelude module for common imports
pub mod prelude {
    pub use crate::alignment::*;
    pub use crate::epigenomics::*;
    pub use crate::error::{DnaError, Result};
    pub use crate::kmer::*;
    pub use crate::pharma::*;
    pub use crate::protein::*;
    pub use crate::types::*;
    pub use crate::variant::*;
}
