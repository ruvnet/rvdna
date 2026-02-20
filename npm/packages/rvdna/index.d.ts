/**
 * @ruvector/rvdna — AI-native genomic analysis and the .rvdna file format.
 *
 * Provides variant calling, protein translation, k-mer vector search,
 * and the compact .rvdna binary format via Rust NAPI-RS bindings.
 */

/**
 * Encode a DNA string to 2-bit packed bytes (4 bases per byte).
 * A=00, C=01, G=10, T=11. Ambiguous bases (N) map to A.
 */
export function encode2bit(sequence: string): Buffer;

/**
 * Decode 2-bit packed bytes back to a DNA string.
 * @param buffer - The 2-bit packed buffer
 * @param length - Number of bases to decode
 */
export function decode2bit(buffer: Buffer, length: number): string;

/**
 * Translate a DNA string to a protein amino acid string.
 * Uses the standard genetic code. Stops at the first stop codon.
 */
export function translateDna(sequence: string): string;

/**
 * Compute cosine similarity between two numeric arrays.
 * Returns a value between -1 and 1.
 */
export function cosineSimilarity(a: number[], b: number[]): number;

export interface RvdnaOptions {
  /** K-mer size (default: 11) */
  k?: number;
  /** Vector dimensions (default: 512) */
  dims?: number;
  /** Block size in bases (default: 500) */
  blockSize?: number;
}

/**
 * Convert a FASTA sequence string to .rvdna binary format.
 * Requires native bindings.
 */
export function fastaToRvdna(sequence: string, options?: RvdnaOptions): Buffer;

export interface RvdnaFile {
  /** Format version */
  version: number;
  /** Sequence length in bases */
  sequenceLength: number;
  /** Decoded DNA sequence */
  sequence: string;
  /** Pre-computed k-mer vector blocks */
  kmerVectors: Array<{
    k: number;
    dimensions: number;
    startPos: number;
    regionLen: number;
    vector: Float32Array;
  }>;
  /** Variant positions and genotype likelihoods */
  variants: Array<{
    position: number;
    refAllele: string;
    altAllele: string;
    likelihoods: [number, number, number];
    quality: number;
  }> | null;
  /** Metadata key-value pairs */
  metadata: Record<string, unknown> | null;
  /** File statistics */
  stats: {
    totalSize: number;
    bitsPerBase: number;
    compressionRatio: number;
  };
}

/**
 * Read a .rvdna file from a Buffer. Returns parsed sections.
 * Requires native bindings.
 */
export function readRvdna(buffer: Buffer): RvdnaFile;

/**
 * Check if native bindings are available for the current platform.
 */
export function isNativeAvailable(): boolean;

/**
 * Direct access to the native NAPI-RS module (null if not available).
 */
export const native: Record<string, Function> | null;

// -------------------------------------------------------------------
// 23andMe Genotyping Pipeline (v0.2.0)
// -------------------------------------------------------------------

/**
 * Normalize a genotype string: uppercase, trim, sort allele pair.
 * "ag" → "AG", "TC" → "CT", "DI" → "DI"
 */
export function normalizeGenotype(gt: string): string;

export interface Snp {
  rsid: string;
  chromosome: string;
  position: number;
  genotype: string;
}

export interface GenotypeData {
  snps: Record<string, Snp>;
  totalMarkers: number;
  noCalls: number;
  chrCounts: Record<string, number>;
  build: 'GRCh37' | 'GRCh38' | 'Unknown';
}

/**
 * Parse a 23andMe raw data file (v4/v5 tab-separated format).
 * Normalizes all genotype strings on load.
 */
export function parse23andMe(text: string): {
  snps: Map<string, Snp>;
  totalMarkers: number;
  noCalls: number;
  chrCounts: Map<string, number>;
  build: string;
};

export interface CypDiplotype {
  gene: string;
  allele1: string;
  allele2: string;
  activity: number;
  phenotype: 'UltraRapid' | 'Normal' | 'Intermediate' | 'Poor';
  confidence: 'Unsupported' | 'Weak' | 'Moderate' | 'Strong';
  rsidsGenotyped: number;
  rsidsMatched: number;
  rsidsTotal: number;
  notes: string[];
  details: string[];
}

/** Call CYP2D6 diplotype from a genotype map */
export function callCyp2d6(gts: Map<string, string>): CypDiplotype;

/** Call CYP2C19 diplotype from a genotype map */
export function callCyp2c19(gts: Map<string, string>): CypDiplotype;

export interface ApoeResult {
  genotype: string;
  rs429358: string;
  rs7412: string;
}

/** Determine APOE genotype from rs429358 + rs7412 */
export function determineApoe(gts: Map<string, string>): ApoeResult;

export interface AnalysisResult {
  data: GenotypeData;
  cyp2d6: CypDiplotype;
  cyp2c19: CypDiplotype;
  apoe: ApoeResult;
  homozygous: number;
  heterozygous: number;
  indels: number;
  hetRatio: number;
}

/**
 * Run the full 23andMe analysis pipeline.
 * @param text - Raw 23andMe file contents
 */
export function analyze23andMe(text: string): AnalysisResult;
