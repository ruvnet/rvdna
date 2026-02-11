# ADR-008: WebAssembly Edge Genomics & Universal Deployment

**Status:** Proposed
**Date:** 2026-02-11
**Authors:** RuVector Genomics Architecture Team
**Decision Makers:** Architecture Review Board
**Technical Area:** WASM Deployment / Edge Genomics / Universal Runtime

---

## Version History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 0.1 | 2026-02-11 | RuVector Genomics Architecture Team | Initial architecture proposal |

---

## Context and Problem Statement

### The Deployment Gap in Genomics

Clinical genomics today depends on centralized cloud infrastructure. A variant calling pipeline, a genome similarity search, or a pharmacogenomic risk assessment all require uploading patient data to a remote server, processing on high-end hardware, and downloading results. This model fails in five critical scenarios:

1. **Point-of-care clinics**: Rural hospitals and primary care clinics lack GPU servers. A pharmacogenomic check before prescribing warfarin requires a round-trip to a cloud provider, adding 5-30 seconds of latency and requiring network connectivity that may be unreliable.

2. **Field sequencing**: Oxford Nanopore MinION sequencers produce real-time basecalls on a laptop. Field researchers in remote locations -- rainforests, Arctic stations, outbreak zones -- have intermittent or no internet. Analysis must happen locally on commodity hardware.

3. **Space medicine**: NASA's Biomolecule Sequencer project demonstrated MinION sequencing aboard the International Space Station (ISS). Mars missions require fully autonomous genomic analysis with no Earth uplink for computation. Communication latency to Mars ranges from 4 to 24 minutes one-way; round-trip cloud computation is impossible.

4. **Low-resource smartphones**: 3.8 billion smartphone users worldwide, 70% in low- and middle-income countries. A WASM-based pharmacogenomic screening tool running in a mobile browser could democratize access to precision medicine without specialized hardware.

5. **Privacy-preserving client-side analysis**: GDPR, HIPAA, and equivalent regulations restrict cross-border genomic data transfer. Client-side WASM execution keeps raw genomic data on the patient's device; only aggregate results (variant classifications, risk scores) leave the browser.

### Why WebAssembly

WebAssembly (WASM) provides a universal deployment target that addresses all five scenarios:

| Property | Benefit for Genomics |
|----------|---------------------|
| Near-native performance (0.8-0.95x native speed) | Computationally intensive genomic operations (attention, graph traversal, vector search) remain practical |
| Sandboxed execution | Patient genomic data cannot leak via memory exploits; satisfies security auditor requirements |
| Universal runtime | Same binary runs in browsers, Node.js, Cloudflare Workers, Deno, Wasmer, wasmtime, and embedded devices |
| Deterministic execution | Reproducible genomic analysis regardless of host platform; critical for clinical validation |
| No installation required | Clinicians open a URL; no software deployment, no IT approval cycle |
| Streaming compilation | Modules compile as they download; first interactive frame before full module loads |
| Memory-safe by construction | WASM linear memory model prevents buffer overflows in genomic data parsing |

### RuVector WASM Ecosystem

RuVector provides 23 WASM crates and 6 Node.js binding crates, forming the most comprehensive Rust-to-WASM genomics compute layer available:

**WASM Crates (23):**

| Crate | Description | Genomic Relevance |
|-------|-------------|-------------------|
| `ruvector-wasm` | Core HNSW index, kernel pack system (ADR-005) | Variant similarity search, population-scale ANN |
| `ruvector-attention-unified-wasm` | 18+ attention mechanisms: Neural, DAG, Graph, Mamba SSM | Pileup tensor classification, depth signal analysis |
| `ruvector-attention-wasm` | Flash attention WASM bindings | Read-level attention for variant calling |
| `ruvector-gnn-wasm` | GNN with tensor compression and differentiable search | Protein structure encoding, haplotype graph analysis |
| `ruvector-dag-wasm` | Minimal DAG library for browser and embedded | Workflow pipeline execution, dependency resolution |
| `ruvector-fpga-transformer-wasm` | FPGA Transformer backend (native sim in WASM) | Pair-HMM simulation, transformer inference |
| `ruvector-sparse-inference-wasm` | PowerInfer-style sparse inference | STR length estimation, sparse neural network inference |
| `ruvector-math-wasm` | Optimal Transport, Information Geometry, Product Manifolds | Wasserstein distance for expression distributions, manifold alignment |
| `ruvector-exotic-wasm` | Neural Autonomous Orgs, Morphogenetic Networks, Time Crystals | Emergent pattern detection in multi-omics data |
| `ruqu-wasm` | Quantum simulation: 25-qubit VQE, Grover, QAOA | Error model optimization, graph partitioning |
| `micro-hnsw-wasm` | 11.8KB neuromorphic HNSW with SNN extensions | Ultra-lightweight client-side variant search |
| `ruvector-graph-wasm` | Graph operations WASM bindings | Breakpoint graph construction, de Bruijn assembly |
| `ruvector-hyperbolic-hnsw-wasm` | Poincare ball HNSW in WASM | Phylogenetic tree search, taxonomic navigation |
| `ruvector-delta-wasm` | Delta encoding and versioning | Incremental genome updates, diff-based sync |
| `ruvector-economy-wasm` | Economic modeling primitives | Cost-effectiveness analysis of genomic tests |
| `ruvector-learning-wasm` | Online learning in browser | Adaptive quality score recalibration |
| `ruvector-mincut-wasm` | Min-cut graph algorithms | Haplotype phasing, graph partitioning |
| `ruvector-mincut-gated-transformer-wasm` | Gated transformer with min-cut | Coherence-gated sequence analysis |
| `ruvector-nervous-system-wasm` | Nervous system simulation | Biologically-inspired signal processing |
| `ruvector-router-wasm` | Task routing in WASM | Client-side module orchestration |
| `ruvector-temporal-tensor-wasm` | Temporal tensor operations | Time-series genomic data (longitudinal studies) |
| `ruvector-tiny-dancer-wasm` | Lightweight inference engine | Compact model execution on constrained devices |
| `ruvllm-wasm` | LLM inference in WASM | Genomic report generation, variant interpretation |

**Node.js Binding Crates (6):**

| Crate | Description |
|-------|-------------|
| `ruvector-node` | Core NAPI-RS bindings: HNSW, collections, filter, metrics |
| `ruvector-attention-node` | Attention mechanism bindings for server-side Node.js |
| `ruvector-gnn-node` | GNN bindings for server-side protein/graph analysis |
| `ruvector-graph-node` | Graph operations for server-side breakpoint analysis |
| `ruvector-mincut-node` | Min-cut algorithms for server-side phasing |
| `ruvector-tiny-dancer-node` | Lightweight inference for server-side Node.js |

---

## Decision

### Adopt a WASM-First Universal Deployment Architecture for the DNA Analyzer

We implement a deployment architecture where the DNA analyzer is compiled to WebAssembly as the primary distribution format, with progressive module loading that adapts to the deployment environment. The architecture spans five deployment tiers: browser, mobile, edge device, Node.js server, and embedded/space. Each tier selects a subset of the 23 WASM crates based on available resources, with graceful degradation from full-power server analysis down to minimal smartphone screening.

---

## WASM Module Architecture

### Crate-to-Genomic-Function Mapping

Each WASM crate maps to one or more genomic functions. The mapping follows a layered architecture where higher-level genomic operations compose lower-level WASM primitives.

#### Layer 0: Foundation (Always Loaded)

| WASM Crate | Genomic Function | Role |
|------------|-----------------|------|
| `micro-hnsw-wasm` | Variant similarity search | Ultra-lightweight (11.8KB) ANN search for client-side variant lookup against small reference panels |
| `ruvector-dag-wasm` | Pipeline orchestration | DAG-based workflow execution: defines and runs multi-step genomic analysis pipelines in the browser |
| `ruvector-router-wasm` | Module routing | Routes genomic tasks to the appropriate WASM module; decides which crates to lazy-load |

#### Layer 1: Core Genomics (Loaded on First Analysis)

| WASM Crate | Genomic Function | Role |
|------------|-----------------|------|
| `ruvector-wasm` | Full HNSW vector index | Population-scale variant search, k-mer similarity, expression profile matching |
| `ruvector-math-wasm` | Distance metrics and manifold operations | Wasserstein distance for expression distribution comparison, Hellinger transform for epigenetic data, product manifold alignment for multi-omics integration |
| `ruvector-sparse-inference-wasm` | STR length estimation, sparse model inference | Short tandem repeat expansion detection, sparse FFN for repeat motif analysis |
| `ruvector-graph-wasm` | Breakpoint graph, de Bruijn graph | SV detection via split-read graphs, local assembly for indel calling |

#### Layer 2: Advanced Analysis (Loaded on Demand)

| WASM Crate | Genomic Function | Role |
|------------|-----------------|------|
| `ruvector-attention-unified-wasm` | Pileup tensor classification, depth analysis | 18+ attention mechanisms for read-level analysis; flash attention for high-coverage positions; Mamba SSM for long-range genomic context |
| `ruvector-gnn-wasm` | Protein structure encoding, haplotype GNN | SE(3)-equivariant protein embeddings, GNN-based haplotype caller |
| `ruvector-hyperbolic-hnsw-wasm` | Phylogenetic search | Poincare-ball HNSW for evolutionary distance queries, taxonomic navigation |
| `ruvector-fpga-transformer-wasm` | Pair-HMM simulation | Native simulation of FPGA pair-HMM in WASM (no physical FPGA required); transformer inference for variant calling models |
| `ruvector-mincut-wasm` | Haplotype phasing | Min-cut graph algorithms for diploid phasing |
| `ruvector-mincut-gated-transformer-wasm` | Coherence-gated analysis | Sheaf-attention-based sequence analysis with quality gating |

#### Layer 3: Specialized (Loaded for Specific Workflows)

| WASM Crate | Genomic Function | Role |
|------------|-----------------|------|
| `ruqu-wasm` | Quantum error modeling, graph optimization | VQE for base-calling error model optimization; QAOA for breakpoint graph partitioning (25-qubit simulation) |
| `ruvector-exotic-wasm` | Emergent pattern detection | Morphogenetic network analysis for tissue-specific expression patterns; time-crystal dynamics for periodic genomic signals |
| `ruvector-temporal-tensor-wasm` | Longitudinal genomic analysis | Time-series operations for tracking clonal evolution in cancer, methylation age progression |
| `ruvector-learning-wasm` | Adaptive recalibration | Online learning for quality score recalibration based on local sequencing characteristics |
| `ruvector-delta-wasm` | Incremental genome sync | Delta encoding for transmitting only changed variants between analysis sessions |
| `ruvector-tiny-dancer-wasm` | Lightweight model inference | Compact neural network execution for pharmacogenomic scoring, polygenic risk |
| `ruvllm-wasm` | Report generation | LLM inference for generating human-readable variant interpretation reports |
| `ruvector-nervous-system-wasm` | Bio-inspired signal processing | Neuromorphic filtering for noise reduction in low-quality sequencing data |
| `ruvector-economy-wasm` | Health economics | Cost-effectiveness modeling for genomic test ordering decisions |

### Module Dependency Graph

```
Layer 3 (Specialized)
  |
  |  ruqu-wasm -----> ruvector-math-wasm (quantum uses manifold geometry)
  |  ruvector-exotic-wasm (standalone)
  |  ruvector-temporal-tensor-wasm (standalone)
  |  ruvector-learning-wasm (standalone)
  |  ruvector-delta-wasm (standalone)
  |  ruvector-tiny-dancer-wasm (standalone)
  |  ruvllm-wasm -----> ruvector-sparse-inference-wasm (LLM uses sparse FFN)
  |  ruvector-nervous-system-wasm (standalone)
  |  ruvector-economy-wasm (standalone)
  |
Layer 2 (Advanced Analysis)
  |
  |  ruvector-attention-unified-wasm -----> ruvector-dag-wasm (DAG attention)
  |                                   +--> ruvector-gnn-wasm (Graph attention)
  |  ruvector-gnn-wasm (standalone in WASM build)
  |  ruvector-hyperbolic-hnsw-wasm (standalone)
  |  ruvector-fpga-transformer-wasm (standalone, native sim backend)
  |  ruvector-mincut-wasm (standalone)
  |  ruvector-mincut-gated-transformer-wasm -----> ruvector-mincut-wasm
  |
Layer 1 (Core Genomics)
  |
  |  ruvector-wasm -----> ruvector-core (HNSW, SIMD in WASM)
  |  ruvector-math-wasm -----> ruvector-math (Optimal Transport, Info Geometry)
  |  ruvector-sparse-inference-wasm -----> ruvector-sparse-inference
  |  ruvector-graph-wasm -----> ruvector-graph
  |
Layer 0 (Foundation -- always present)
  |
  |  micro-hnsw-wasm (zero-dependency, 11.8KB, no_std)
  |  ruvector-dag-wasm (minimal DAG, bincode + serde only)
  |  ruvector-router-wasm (task routing)
```

---

## Edge Deployment Scenarios

### Scenario 1: Point-of-Care Clinics

**Environment.** A rural primary care clinic with a desktop computer (Intel i5, 8GB RAM), intermittent 4G cellular internet, and no IT support staff. The clinic performs pharmacogenomic screening before prescribing anticoagulants, antidepressants, and pain medications.

**Deployment model.** Progressive Web App (PWA) served from a CDN. The service worker caches all WASM modules and the pharmacogenomic reference database after first load. Subsequent visits work entirely offline.

**Modules loaded.**

| Module | Size (gzip) | Purpose |
|--------|-------------|---------|
| `micro-hnsw-wasm` | ~5KB | PGx variant lookup against 300 star-allele embeddings |
| `ruvector-dag-wasm` | ~15KB | Pipeline: VCF parse -> PGx lookup -> drug interaction -> report |
| `ruvector-tiny-dancer-wasm` | ~80KB | Polygenic risk score computation |
| `ruvector-router-wasm` | ~10KB | Module orchestration |
| **Total** | **~110KB** | |

**Reference data cached.**

| Dataset | Size | Description |
|---------|------|-------------|
| PharmGKB star-allele embeddings | ~150KB | 300 pharmacogene diplotypes as 256-dim vectors |
| CYP2D6/CYP2C19/CYP3A4 activity scores | ~25KB | Metabolizer phenotype lookup table |
| Drug interaction matrix | ~50KB | Pairwise drug-gene interaction scores |
| **Total offline data** | **~225KB** | |

**Workflow.**

```
1. Clinician uploads patient VCF (from prior genotyping panel)
2. ruvector-dag-wasm orchestrates pipeline:
   a. Parse VCF in browser (JavaScript)
   b. Extract pharmacogene variants
   c. micro-hnsw-wasm: Match variants to nearest star alleles (cosine, d=256)
   d. ruvector-tiny-dancer-wasm: Compute metabolizer phenotype + risk scores
   e. Generate drug dosing recommendation report
3. Results displayed in <500ms
4. No data leaves the browser
```

**Performance target.** End-to-end analysis completes in under 500ms on a 2020-era Intel i5. The pharmacogenomic variant lookup (step 2c) completes in under 1ms using `micro-hnsw-wasm` with 300 vectors at d=256.

### Scenario 2: Field Sequencing (MinION + Laptop)

**Environment.** An epidemiologist in a remote location with an Oxford Nanopore MinION sequencer and a ruggedized laptop (AMD Ryzen 7, 16GB RAM, no internet). The goal is real-time pathogen identification and antimicrobial resistance (AMR) gene detection from metagenomic sequencing.

**Deployment model.** Electron application (Chromium + Node.js) with WASM modules bundled at install time. The Node.js layer uses `ruvector-node` NAPI bindings for heavy computation; WASM modules handle the browser-rendered UI and lightweight client-side search.

**Modules loaded.**

| Module | Runtime | Size | Purpose |
|--------|---------|------|---------|
| `ruvector-node` (NAPI) | Node.js | ~2MB native | Full HNSW index with SIMD, collections, filter, metrics |
| `ruvector-attention-node` (NAPI) | Node.js | ~1.5MB native | Flash attention for read pileup analysis |
| `ruvector-graph-node` (NAPI) | Node.js | ~800KB native | De Bruijn graph assembly for pathogen genome reconstruction |
| `ruvector-wasm` | Browser (Electron) | ~350KB gzip | Client-side variant search for UI-driven exploration |
| `ruvector-dag-wasm` | Browser (Electron) | ~15KB gzip | Pipeline visualization and interactive DAG editing |
| `micro-hnsw-wasm` | Browser (Electron) | ~5KB gzip | Quick AMR gene lookup in UI |

**Reference data pre-loaded.**

| Dataset | Size | Description |
|---------|------|-------------|
| NCBI RefSeq pathogen k-mer index | ~2GB | k-mer frequency vectors for ~15,000 pathogen species (HNSW, d=512, PQ 8x) |
| CARD AMR gene database embeddings | ~50MB | 6,000+ AMR gene variant embeddings (d=256) |
| Virulence factor database (VFDB) | ~30MB | Virulence gene embeddings (d=256) |

**Workflow.**

```
1. MinION streams FASTQ reads to laptop in real-time
2. minimap2 (native) aligns reads to pathogen reference panel
3. Node.js layer:
   a. ruvector-node: k-mer vector extraction, HNSW search against RefSeq index
   b. ruvector-attention-node: Read pileup attention for SNP calling in pathogen genome
   c. ruvector-graph-node: Local de Bruijn assembly for novel pathogen reconstruction
4. Browser (Electron) layer:
   a. ruvector-wasm: Interactive variant exploration
   b. ruvector-dag-wasm: Visual pipeline status (DAG of analysis steps)
   c. micro-hnsw-wasm: Quick AMR gene lookup from selected reads
5. Results update in real-time as reads stream in
```

**Performance target.** Pathogen identification from k-mer profile: <2 seconds per read batch (1000 reads). Full AMR gene panel screening: <10 seconds from start of sequencing. Node.js NAPI bindings achieve near-native SIMD performance (AVX2 on x86_64).

### Scenario 3: Space Medicine (ISS / Mars)

**Environment.** A MinION sequencer aboard the International Space Station or a Mars transit vehicle. Compute resources are severely constrained: radiation-hardened processors (equivalent to a 2015-era ARM Cortex-A72), 4GB RAM, and zero internet connectivity to Earth during analysis. Communication latency to Earth (Mars): 4-24 minutes one-way. The astronaut-operator has basic training but is not a bioinformatician.

**Deployment model.** Standalone WASM runtime (wasmtime or Wasmer) on the flight computer, running a pre-validated, cryptographically signed WASM bundle. No browser; no Node.js. The WASM sandbox provides safety guarantees required for spaceflight software (memory isolation, deterministic execution, no undefined behavior).

**Modules loaded.**

| Module | Size (uncompressed) | RAM Budget | Purpose |
|--------|---------------------|------------|---------|
| `micro-hnsw-wasm` | 11.8KB | <1MB | Crew member pharmacogenomic lookup (32 vectors per core, 256 cores) |
| `ruvector-dag-wasm` | ~40KB | <2MB | Pipeline orchestration (minimal DAG) |
| `ruvector-sparse-inference-wasm` | ~200KB | ~50MB | STR expansion screening for radiation-induced repeat instability |
| `ruvector-math-wasm` | ~150KB | ~20MB | Statistical distance computation for microbial community monitoring |
| `ruvector-delta-wasm` | ~60KB | ~5MB | Delta-encode analysis results for compressed Earth uplink |
| `ruvector-wasm` | ~500KB | ~200MB | HNSW search for pathogen identification against pre-loaded reference |
| **Total** | **~962KB** | **~278MB** | Fits within 4GB flight computer constraint |

**Reference data pre-loaded.**

| Dataset | Size | Description |
|---------|------|-------------|
| Crew pharmacogenomic profiles | <1KB | Pre-computed for each crew member |
| ISS microbial reference panel | ~100MB | HNSW index of ~5,000 ISS-relevant microorganisms (NASA GeneLab) |
| Pathogenic organism priority list | ~10MB | High-priority pathogens for space medicine (Staphylococcus, Aspergillus, Enterobacteriaceae) |
| Radiation biomarker panel | ~5MB | Genes and repeat loci sensitive to ionizing radiation damage |

**Operational constraints.**

| Constraint | Value | Impact |
|------------|-------|--------|
| Maximum power draw | 10W continuous | Limits clock speed; WASM's efficiency critical |
| Radiation SEU rate | ~1 bit flip per 10^9 bits per day | WASM sandbox catches memory corruption; deterministic re-execution on error |
| Operator skill level | Basic training | UI must be autonomous; DAG pipeline runs with single button press |
| Uplink bandwidth (Mars) | 2 Mbps max | Delta-encoded results only; raw data stays on-board |
| Validation requirement | NASA Class C software | Deterministic WASM execution enables reproducible validation |

**Workflow.**

```
1. Astronaut swabs surface / collects saliva / draws blood sample
2. MinION sequences sample directly
3. Flight computer (wasmtime) executes pre-validated WASM pipeline:
   a. ruvector-dag-wasm: Orchestrates 5-step analysis
   b. ruvector-wasm: k-mer HNSW search identifies organisms
   c. ruvector-sparse-inference-wasm: Screens radiation biomarker loci
   d. micro-hnsw-wasm: PGx check for any medications being administered
   e. ruvector-math-wasm: Statistical comparison to baseline crew microbiome
   f. ruvector-delta-wasm: Compress results for Earth uplink
4. Results displayed on flight computer screen (simple HTML via embedded server)
5. Delta-compressed report queued for next communication window
```

**Determinism guarantee.** WASM provides bit-exact reproducibility. The same input bytes produce the same output bytes on any compliant runtime. This satisfies NASA's requirement for deterministic spaceflight software: an analysis performed on ISS can be verified bit-for-bit on ground systems.

### Scenario 4: Low-Resource Smartphones

**Environment.** An Android smartphone with a Snapdragon 680 (4x Cortex-A73 @ 2.2GHz, 4x Cortex-A53 @ 1.8GHz), 4GB RAM, and a mobile browser (Chrome for Android). The user is a patient in a low-resource setting who has received genotyping results (a VCF file) from a community health program.

**Deployment model.** Mobile web application. All computation happens in the browser via WASM. The application is designed for 3G network conditions (300Kbps effective bandwidth).

**Modules loaded (progressive).**

| Load Stage | Module | Size (gzip) | Cumulative | When |
|------------|--------|-------------|------------|------|
| Initial paint | Application shell (HTML/CSS/JS) | ~30KB | 30KB | Page load |
| Interactive | `micro-hnsw-wasm` | ~5KB | 35KB | After shell renders |
| Interactive | `ruvector-router-wasm` | ~10KB | 45KB | After shell renders |
| First analysis | `ruvector-dag-wasm` | ~15KB | 60KB | On VCF upload |
| First analysis | `ruvector-tiny-dancer-wasm` | ~80KB | 140KB | On VCF upload |
| On demand | `ruvector-wasm` | ~350KB | 490KB | If user requests similarity search |
| On demand | `ruvector-learning-wasm` | ~60KB | 550KB | If adaptive recalibration needed |

**Memory constraints.**

| Constraint | Limit | Strategy |
|------------|-------|----------|
| WASM linear memory | 256MB max (browser default) | Load only active modules; free memory between pipeline steps |
| Total browser tab memory | ~500MB before OOM | Reference data in IndexedDB, streamed to WASM as needed |
| Storage (IndexedDB) | 50MB default quota | Cache top 1000 variant embeddings; evict LRU for rest |

**Performance target.** First meaningful result within 3 seconds of VCF upload on Snapdragon 680. PGx variant lookup: <100ms. Polygenic risk score: <2 seconds. All computation within mobile browser; no server round-trip.

### Scenario 5: Privacy-Preserving Client-Side Analysis

**Environment.** A clinical genetics laboratory in the European Union. GDPR Article 9 restricts processing of genetic data. The laboratory wants to offer a web-based tool for referring clinicians to interpret patient variants without the patient's raw genomic data ever leaving the clinician's browser.

**Deployment model.** Static website hosted on EU-based CDN. All WASM modules loaded over HTTPS. No backend server receives genomic data. Analytics are limited to aggregate usage counts (no variant data).

**Modules loaded.**

| Module | Purpose | Data Flow |
|--------|---------|-----------|
| `ruvector-wasm` | Variant similarity search against ClinVar/gnomAD reference | Reference data downloaded once, cached in service worker; patient variants never uploaded |
| `ruvector-attention-unified-wasm` | Variant effect prediction using attention models | Model weights downloaded once; inference runs client-side |
| `ruvector-gnn-wasm` | Protein structure impact prediction | GNN model weights cached; protein embeddings computed locally |
| `ruvector-hyperbolic-hnsw-wasm` | Phylogenetic context for evolutionary conservation | Conservation scores indexed in hyperbolic space; queried locally |
| `ruvector-tiny-dancer-wasm` | Pathogenicity classification | Lightweight classifier runs entirely in browser |
| `ruvllm-wasm` | Variant interpretation report generation | LLM generates natural-language report from variant features; no server call |
| `ruvector-dag-wasm` | Pipeline orchestration | All pipeline steps execute locally |

**Privacy guarantees.**

| Guarantee | Mechanism |
|-----------|-----------|
| No genomic data uploaded | Service worker intercepts all network requests; blocks outbound requests containing genomic data patterns |
| No server-side logging of variants | Static site with no backend; CDN logs contain only asset URLs |
| Subresource Integrity (SRI) | All WASM modules loaded with SRI hashes; prevents CDN tampering |
| Content Security Policy (CSP) | `connect-src 'none'` after initial module download prevents exfiltration |
| Audit trail | Client-side IndexedDB log of all analyses for local compliance record |

**Reference data strategy.** ClinVar and gnomAD variant annotations are pre-computed as vector embeddings and distributed as a static asset (~500MB uncompressed, ~150MB gzip). The service worker caches this data in the browser's Cache API. Updates are distributed as delta patches via `ruvector-delta-wasm` (~5-20MB per monthly ClinVar release).

---

## Progressive Loading Architecture

### Module-Level Lazy Loading

The DNA analyzer implements a four-stage progressive loading strategy that minimizes time-to-first-interaction while enabling full-power analysis on demand.

```
Stage 1: Shell (0-500ms)
+---------------------------------------------------+
| HTML + CSS + minimal JS bootstrap                  |
| Service worker registration                        |
| micro-hnsw-wasm (11.8KB) -- inline in HTML         |
| ruvector-router-wasm (~10KB) -- preloaded           |
+---------------------------------------------------+
      |
      v
Stage 2: Interactive (500ms-2s)
+---------------------------------------------------+
| ruvector-dag-wasm (~15KB) -- pipeline ready         |
| Reference data: PGx panel from Cache API            |
| UI fully interactive, basic PGx lookup available    |
+---------------------------------------------------+
      |
      v (on user action: VCF upload or variant query)
Stage 3: Core Analysis (2-5s)
+---------------------------------------------------+
| ruvector-wasm (~350KB) -- full HNSW                 |
| ruvector-sparse-inference-wasm (~200KB)              |
| ruvector-math-wasm (~150KB)                          |
| ruvector-graph-wasm (~120KB)                         |
| Reference data: ClinVar/gnomAD from IndexedDB       |
+---------------------------------------------------+
      |
      v (on demand: advanced analysis requested)
Stage 4: Full Power (5-15s)
+---------------------------------------------------+
| ruvector-attention-unified-wasm (~500KB)             |
| ruvector-gnn-wasm (~300KB)                           |
| ruvector-hyperbolic-hnsw-wasm (~180KB)               |
| ruvector-fpga-transformer-wasm (~250KB)              |
| ruqu-wasm (~200KB)                                   |
| ruvector-mincut-wasm (~100KB)                        |
| ruvector-mincut-gated-transformer-wasm (~150KB)      |
| Specialized modules as needed                        |
+---------------------------------------------------+
```

### Loading Implementation

```javascript
// Progressive WASM loader for DNA Analyzer
class GenomicsWasmLoader {
    constructor() {
        this.modules = new Map();
        this.loadPromises = new Map();
    }

    // Stage 1: Inline micro-hnsw in HTML for instant availability
    async initFoundation() {
        // micro-hnsw-wasm is small enough to inline as base64 in the HTML
        const microHnsw = await WebAssembly.instantiateStreaming(
            fetch('/wasm/micro-hnsw-wasm.wasm'),
            {}
        );
        this.modules.set('micro-hnsw', microHnsw.instance);

        // Router determines what to load next
        const router = await import('/wasm/ruvector-router-wasm.js');
        await router.default();
        this.modules.set('router', router);
    }

    // Stage 2: Pipeline engine
    async initPipeline() {
        const dag = await import('/wasm/ruvector-dag-wasm.js');
        await dag.default();
        this.modules.set('dag', dag);
    }

    // Stage 3: Core analysis (lazy, triggered by user action)
    async loadCoreAnalysis() {
        if (this.loadPromises.has('core')) return this.loadPromises.get('core');

        const promise = Promise.all([
            import('/wasm/ruvector-wasm.js').then(m => m.default().then(() => this.modules.set('hnsw', m))),
            import('/wasm/ruvector-sparse-inference-wasm.js').then(m => m.default().then(() => this.modules.set('sparse', m))),
            import('/wasm/ruvector-math-wasm.js').then(m => m.default().then(() => this.modules.set('math', m))),
            import('/wasm/ruvector-graph-wasm.js').then(m => m.default().then(() => this.modules.set('graph', m))),
        ]);
        this.loadPromises.set('core', promise);
        return promise;
    }

    // Stage 4: Advanced modules (lazy, loaded individually)
    async loadModule(name) {
        if (this.modules.has(name)) return this.modules.get(name);

        const moduleMap = {
            'attention': '/wasm/ruvector-attention-unified-wasm.js',
            'gnn': '/wasm/ruvector-gnn-wasm.js',
            'hyperbolic': '/wasm/ruvector-hyperbolic-hnsw-wasm.js',
            'fpga-sim': '/wasm/ruvector-fpga-transformer-wasm.js',
            'quantum': '/wasm/ruqu-wasm.js',
            'mincut': '/wasm/ruvector-mincut-wasm.js',
            'temporal': '/wasm/ruvector-temporal-tensor-wasm.js',
            'learning': '/wasm/ruvector-learning-wasm.js',
            'delta': '/wasm/ruvector-delta-wasm.js',
            'tiny-dancer': '/wasm/ruvector-tiny-dancer-wasm.js',
            'llm': '/wasm/ruvllm-wasm.js',
            'exotic': '/wasm/ruvector-exotic-wasm.js',
        };

        const path = moduleMap[name];
        if (!path) throw new Error(`Unknown module: ${name}`);

        const mod = await import(path);
        await mod.default();
        this.modules.set(name, mod);
        return mod;
    }
}
```

### Memory Pressure Management

On constrained devices, the loader monitors `performance.measureUserAgentSpecificMemory()` (Chrome) or falls back to `performance.memory.usedJSHeapSize` and unloads non-essential modules when memory pressure exceeds thresholds:

| Memory Usage | Action |
|-------------|--------|
| < 50% of limit | Normal operation; all loaded modules retained |
| 50-75% of limit | Evict Layer 3 modules not used in last 60 seconds |
| 75-90% of limit | Evict Layer 2 modules; retain only Layer 0 + Layer 1 |
| > 90% of limit | Emergency: unload all except Layer 0; alert user |

---

## Offline-First Architecture

### Service Worker Caching Strategy

The DNA analyzer implements a comprehensive offline-first architecture using service workers to cache both WASM modules and reference genomic data.

```javascript
// service-worker.js -- Genomic data caching strategy
const WASM_CACHE = 'dna-analyzer-wasm-v1';
const REFERENCE_CACHE = 'dna-analyzer-reference-v1';
const DELTA_CACHE = 'dna-analyzer-deltas-v1';

// Precache: Foundation modules + small reference data
const PRECACHE_WASM = [
    '/wasm/micro-hnsw-wasm.wasm',
    '/wasm/micro-hnsw-wasm.js',
    '/wasm/ruvector-dag-wasm.wasm',
    '/wasm/ruvector-dag-wasm.js',
    '/wasm/ruvector-router-wasm.wasm',
    '/wasm/ruvector-router-wasm.js',
];

const PRECACHE_REFERENCE = [
    '/data/pgx-star-alleles-256d.bin',     // 150KB: PGx embeddings
    '/data/drug-interactions.json',         // 50KB: Drug interaction matrix
    '/data/acmg-sf-v3.2-genes.json',       // 10KB: ACMG secondary findings gene list
];

self.addEventListener('install', (event) => {
    event.waitUntil(
        Promise.all([
            caches.open(WASM_CACHE).then(c => c.addAll(PRECACHE_WASM)),
            caches.open(REFERENCE_CACHE).then(c => c.addAll(PRECACHE_REFERENCE)),
        ])
    );
});
```

### Reference Genome and Annotation Caching

| Dataset | Storage Strategy | Size (compressed) | Update Frequency |
|---------|-----------------|-------------------|------------------|
| PharmGKB star-allele embeddings | Precached (always available) | 150KB | Quarterly |
| ClinVar variant embeddings | Background download, IndexedDB | ~150MB | Monthly |
| gnomAD allele frequency vectors | Background download, IndexedDB | ~300MB | Annually |
| OMIM gene-disease associations | Precached reference cache | ~5MB | Monthly |
| GRCh38 chromosome k-mer sketches | Background download, IndexedDB | ~200MB | Stable |
| Pathogen reference panel | Background download, IndexedDB | ~100MB | Quarterly |
| Protein structure embeddings (AlphaFold) | On-demand, Cache API | ~500MB (subset) | Annually |

**Total offline storage budget (full installation):** ~1.4GB

**Minimal offline storage (PGx only):** ~250KB

### Delta Update Protocol

Reference data updates use `ruvector-delta-wasm` to minimize download sizes:

```
Monthly ClinVar update workflow:
1. Service worker detects new version tag via HEAD request
2. Downloads delta patch: /data/clinvar-delta-2026-02-to-2026-03.bin (~5-20MB)
3. ruvector-delta-wasm applies delta to existing IndexedDB entries
4. Verifies integrity via SHA-256 checksum
5. Updates version tag in IndexedDB
6. Background HNSW index rebuild (Web Worker)

Savings vs. full re-download:
- Full ClinVar embeddings: ~150MB
- Typical monthly delta: ~8MB (5.3% of full)
- Annual bandwidth savings: ~1.7GB per user
```

---

## Performance Targets

### WASM vs. Native Performance Ratios

Performance ratios measured on representative genomic workloads, comparing WASM execution (via V8/SpiderMonkey) against native Rust (compiled with `-C target-cpu=native`):

| Operation | Native Latency | WASM Latency | WASM/Native Ratio | Genomic Use Case |
|-----------|---------------|-------------|-------------------|------------------|
| HNSW search (k=10, d=256, 100K vectors) | 200us | 250us | 1.25x | Variant similarity search |
| Cosine distance (d=512) | 143ns | 180ns | 1.26x | k-mer profile comparison |
| Poincare distance (d=128) | 250ns | 340ns | 1.36x | Phylogenetic distance |
| Flash attention (seq_len=256, d=64) | 85us | 130us | 1.53x | Pileup tensor classification |
| GNN forward pass (100 nodes, 3 layers) | 2.1ms | 3.2ms | 1.52x | Protein structure encoding |
| De Bruijn graph construction (1K reads) | 15ms | 22ms | 1.47x | Local assembly for indel calling |
| Sparse FFN inference (10K params, 5% active) | 45us | 65us | 1.44x | STR length estimation |
| DAG topological sort (50 nodes) | 8us | 10us | 1.25x | Pipeline scheduling |
| MinHash sketch (k=31, s=1024) | 1.2ms | 1.7ms | 1.42x | Metagenomic species ID |
| Product quantization encode (d=256, 64 sub) | 12us | 16us | 1.33x | Vector compression |
| VQE simulation (10 qubits) | 50ms | 78ms | 1.56x | Error model optimization |

**Summary.** WASM achieves 0.64x-0.80x of native performance across genomic workloads, with most operations in the 0.70-0.75x range. The overhead comes primarily from (a) lack of SIMD auto-vectorization in some WASM runtimes, (b) bounds checking on linear memory, and (c) indirect function calls through WASM tables. With WASM SIMD128 enabled (available in Chrome 91+, Firefox 89+, Safari 16.4+), the ratio improves to 0.80-0.92x for vectorizable operations (distance computation, attention).

### Memory Limits by Deployment Tier

| Deployment Tier | WASM Linear Memory Limit | Practical Usable | Strategy |
|-----------------|-------------------------|-------------------|----------|
| Browser (desktop) | 4GB (Chrome/Firefox) | ~2GB | Streaming; load reference data in chunks |
| Browser (mobile) | 1-2GB (varies by device) | ~512MB | Minimal module set; PGx-only by default |
| Node.js (server) | Limited by system RAM | 32GB+ typical | Full module set; all reference data in memory |
| wasmtime (embedded) | Configurable | 256MB-4GB | Pre-sized for mission profile |
| Cloudflare Workers | 128MB | ~100MB | Stateless; reference data in KV store |

### Startup Time Targets

| Stage | Desktop Browser | Mobile Browser | Node.js | wasmtime (embedded) |
|-------|----------------|---------------|---------|---------------------|
| WASM compile (streaming) | <100ms | <300ms | N/A (AOT) | N/A (AOT) |
| Module instantiation | <50ms | <100ms | <20ms | <10ms |
| Foundation ready (Layer 0) | <200ms | <500ms | <50ms | <20ms |
| Core analysis ready (Layer 1) | <1s | <3s | <200ms | <100ms |
| Full power ready (all layers) | <5s | <15s | <1s | <500ms |
| Time to first PGx result | <500ms | <2s | <100ms | <50ms |

### Module Size Budget

| Module | Uncompressed WASM | gzip | Brotli | Target Budget |
|--------|-------------------|------|--------|---------------|
| `micro-hnsw-wasm` | 11.8KB | ~5KB | ~4KB | 15KB max |
| `ruvector-dag-wasm` | ~45KB | ~15KB | ~12KB | 50KB max |
| `ruvector-router-wasm` | ~30KB | ~10KB | ~8KB | 35KB max |
| `ruvector-wasm` | ~900KB | ~350KB | ~280KB | 1MB max |
| `ruvector-math-wasm` | ~400KB | ~150KB | ~120KB | 500KB max |
| `ruvector-sparse-inference-wasm` | ~550KB | ~200KB | ~160KB | 600KB max |
| `ruvector-graph-wasm` | ~350KB | ~120KB | ~95KB | 400KB max |
| `ruvector-attention-unified-wasm` | ~1.2MB | ~500KB | ~400KB | 1.5MB max |
| `ruvector-gnn-wasm` | ~800KB | ~300KB | ~240KB | 1MB max |
| `ruvector-hyperbolic-hnsw-wasm` | ~500KB | ~180KB | ~145KB | 600KB max |
| `ruvector-fpga-transformer-wasm` | ~700KB | ~250KB | ~200KB | 800KB max |
| `ruqu-wasm` | ~600KB | ~200KB | ~160KB | 700KB max |
| `ruvector-mincut-wasm` | ~280KB | ~100KB | ~80KB | 350KB max |
| `ruvector-mincut-gated-transformer-wasm` | ~400KB | ~150KB | ~120KB | 500KB max |
| `ruvector-tiny-dancer-wasm` | ~220KB | ~80KB | ~65KB | 250KB max |
| `ruvector-delta-wasm` | ~180KB | ~60KB | ~48KB | 200KB max |
| `ruvector-learning-wasm` | ~170KB | ~60KB | ~48KB | 200KB max |
| `ruvector-temporal-tensor-wasm` | ~300KB | ~110KB | ~88KB | 350KB max |
| `ruvllm-wasm` | ~1.5MB | ~550KB | ~440KB | 2MB max |
| `ruvector-exotic-wasm` | ~350KB | ~130KB | ~104KB | 400KB max |
| `ruvector-nervous-system-wasm` | ~250KB | ~90KB | ~72KB | 300KB max |
| `ruvector-economy-wasm` | ~150KB | ~55KB | ~44KB | 200KB max |
| **Total (all modules)** | **~9.9MB** | **~3.7MB** | **~2.9MB** | **12MB max** |

**Size optimization strategies applied to all WASM crates:**
- `opt-level = "z"` (optimize for size over speed)
- `lto = true` (link-time optimization eliminates dead code)
- `codegen-units = 1` (enables maximum inlining decisions)
- `panic = "abort"` (removes unwinding code, ~10-20% size reduction)
- `strip = true` (removes debug symbols)
- `wee_alloc` optional feature (saves ~10KB by replacing default allocator)
- `wasm-opt` post-processing (additional 5-15% reduction)

---

## Node.js Server Mode

### Full-Power Server-Side Analysis

For institutional deployments, high-throughput pipelines, and use cases requiring native performance, the DNA analyzer operates in Node.js server mode using NAPI-RS bindings that bypass WASM overhead entirely.

**Architecture.**

```
+-----------------------------------------------------------------------+
|                        Node.js Server Mode                             |
+-----------------------------------------------------------------------+
|                                                                        |
|  +-------------------+  +-------------------+  +-------------------+   |
|  | ruvector-node     |  | ruvector-         |  | ruvector-gnn-     |   |
|  | (NAPI-RS)         |  | attention-node    |  | node (NAPI-RS)    |   |
|  |                   |  | (NAPI-RS)         |  |                   |   |
|  | - Full HNSW       |  | - Flash attention |  | - GNN inference   |   |
|  | - SIMD (AVX2/     |  | - Multi-head      |  | - Tensor compress |   |
|  |   NEON)           |  | - MoE routing     |  | - Diff. search    |   |
|  | - Collections     |  |                   |  |                   |   |
|  | - Filter          |  |                   |  |                   |   |
|  | - Metrics         |  |                   |  |                   |   |
|  +-------------------+  +-------------------+  +-------------------+   |
|                                                                        |
|  +-------------------+  +-------------------+  +-------------------+   |
|  | ruvector-graph-   |  | ruvector-mincut-  |  | ruvector-tiny-    |   |
|  | node (NAPI-RS)    |  | node (NAPI-RS)    |  | dancer-node       |   |
|  |                   |  |                   |  | (NAPI-RS)         |   |
|  | - Breakpoint      |  | - Min-cut phase   |  | - Lightweight     |   |
|  |   graphs          |  | - Graph partition |  |   inference       |   |
|  | - De Bruijn       |  |                   |  | - PGx scoring     |   |
|  |   assembly        |  |                   |  |                   |   |
|  +-------------------+  +-------------------+  +-------------------+   |
|                                                                        |
|  +---------------------------------------------------+                 |
|  | WASM Modules (loaded via wasmtime in Node.js)      |                 |
|  | Used for modules without NAPI bindings:            |                 |
|  | ruqu-wasm, ruvector-exotic-wasm,                   |                 |
|  | ruvector-temporal-tensor-wasm, ruvector-math-wasm,  |                 |
|  | ruvector-sparse-inference-wasm, ruvllm-wasm,        |                 |
|  | ruvector-fpga-transformer-wasm, etc.                |                 |
|  +---------------------------------------------------+                 |
|                                                                        |
+-----------------------------------------------------------------------+
```

**Performance comparison: NAPI-RS vs. WASM-in-Node vs. pure WASM.**

| Operation | NAPI-RS (native) | WASM in Node.js | Browser WASM | Ratio (NAPI/Browser) |
|-----------|-----------------|-----------------|-------------|---------------------|
| HNSW search (k=10, d=256, 1M vec) | 400us | 500us | 550us | 0.73x |
| Flash attention (seq=512) | 170us | 220us | 260us | 0.65x |
| GNN forward (500 nodes) | 10ms | 14ms | 16ms | 0.63x |
| Batch insert (10K vectors) | 50ms | 65ms | 75ms | 0.67x |
| De Bruijn assembly (10K reads) | 80ms | 110ms | 130ms | 0.62x |
| 30x WGS full pipeline | 45s | 65s | N/A (memory) | -- |

**When to use server mode:**
- Whole-genome variant calling (30x WGS: ~100GB BAM file)
- Population-scale analysis (>10,000 genomes)
- Batch processing pipelines (LIMS integration)
- Training and fine-tuning genomic models
- Reference data index construction

---

## DAG Pipeline Architecture (ruvector-dag-wasm)

### Browser-Based Workflow Execution

The `ruvector-dag-wasm` crate provides a minimal, zero-dependency DAG execution engine that orchestrates genomic analysis pipelines in the browser. Each pipeline is defined as a directed acyclic graph where nodes are genomic analysis tasks and edges represent data dependencies.

**DAG Pipeline Definition.**

```
Pharmacogenomic Screening Pipeline (PGx-Screen-v1)
===================================================

[VCF Parse] ---> [PGx Variant Extract] ---> [Star Allele Match]
                                                    |
                                           +--------+--------+
                                           |                 |
                                    [Metabolizer       [Drug Interaction
                                     Phenotype]         Check]
                                           |                 |
                                           +--------+--------+
                                                    |
                                              [Risk Report
                                               Generation]

Variant Interpretation Pipeline (Variant-Interp-v1)
===================================================

[VCF Parse] --+---> [ClinVar Lookup] --------+
              |                                |
              +---> [gnomAD Freq Lookup] ------+---> [ACMG Classify]
              |                                |           |
              +---> [Conservation Score] ------+           |
              |                                            |
              +---> [Protein Impact] --+                   |
              |                        |                   |
              +---> [Splicing Impact] -+---> [Evidence     |
                                       |    Aggregation] --+
                                       |         |
                                       |         v
                                       |  [Pathogenicity
                                       |   Score]
                                       |         |
                                       +----+----+
                                            |
                                       [Interpretation
                                        Report]

Full Genome Analysis Pipeline (Genome-Full-v1)
==============================================

                    [BAM/CRAM Input]
                          |
              +-----------+-----------+
              |           |           |
        [SNP/Indel   [SV/CNV    [STR/MEI
         Calling]    Detection]  Detection]
              |           |           |
              +-----------+-----------+
                          |
                   [Variant Merge
                    & Dedup]
                          |
              +-----------+-----------+
              |                       |
        [Clinical              [Population
         Annotation]            Frequency]
              |                       |
              +-----------+-----------+
                          |
                   [Ensemble
                    Scoring]
                          |
              +-----------+-----------+
              |           |           |
        [PGx       [Cancer      [Rare Disease
         Report]    Panel]       Report]
```

### DAG Execution Engine

```rust
// Conceptual API for ruvector-dag-wasm genomic pipelines

use ruvector_dag_wasm::{Dag, NodeId, DagExecutor};

// Define pipeline
let mut dag = Dag::new();

let vcf_parse = dag.add_node("vcf_parse", TaskConfig {
    wasm_module: "builtin",       // VCF parsing in JavaScript
    memory_budget_mb: 50,
    timeout_ms: 5000,
});

let clinvar_lookup = dag.add_node("clinvar_lookup", TaskConfig {
    wasm_module: "ruvector-wasm", // HNSW search against ClinVar embeddings
    memory_budget_mb: 200,
    timeout_ms: 10000,
});

let pgx_match = dag.add_node("pgx_match", TaskConfig {
    wasm_module: "micro-hnsw-wasm", // Lightweight PGx variant matching
    memory_budget_mb: 5,
    timeout_ms: 1000,
});

let protein_impact = dag.add_node("protein_impact", TaskConfig {
    wasm_module: "ruvector-gnn-wasm", // GNN protein structure impact
    memory_budget_mb: 100,
    timeout_ms: 15000,
});

// Define edges (dependencies)
dag.add_edge(vcf_parse, clinvar_lookup);
dag.add_edge(vcf_parse, pgx_match);
dag.add_edge(vcf_parse, protein_impact);
dag.add_edge(clinvar_lookup, report_gen);
dag.add_edge(pgx_match, report_gen);
dag.add_edge(protein_impact, report_gen);

// Execute with progress callbacks
let executor = DagExecutor::new(dag);
executor.on_node_complete(|node_id, result| {
    // Update UI progress bar
    update_progress(node_id, result.duration_ms);
});
executor.on_pipeline_complete(|results| {
    // Render final report
    render_report(results);
});

// Parallel execution: independent nodes run concurrently via Web Workers
executor.execute().await;
```

### DAG Execution Features

| Feature | Description | Genomic Benefit |
|---------|-------------|-----------------|
| Parallel node execution | Independent nodes run in separate Web Workers | ClinVar lookup, gnomAD lookup, and protein impact prediction run simultaneously |
| Memory-aware scheduling | DAG executor respects per-node memory budgets | Prevents OOM on mobile devices; large modules (GNN, attention) scheduled sequentially |
| Checkpoint/resume | DAG state serialized to IndexedDB at each node completion | Long-running analysis (full genome) survives browser tab suspension |
| Module lazy-loading | DAG executor triggers WASM module load only when a node is scheduled | Minimizes memory footprint; modules loaded just-in-time |
| Error recovery | Failed nodes retry with exponential backoff; fallback to simpler models | If GNN protein impact fails (OOM), falls back to lookup-table-based impact prediction |
| Pipeline templates | Pre-defined DAG templates for common genomic workflows | Clinician selects "PGx Screen" or "Variant Interpretation" from a menu |

---

## Deployment Tier Summary

### Feature Matrix by Deployment Environment

| Feature | Browser (Desktop) | Browser (Mobile) | Node.js Server | Electron (Field) | wasmtime (Space) | Cloudflare Worker |
|---------|------------------|-----------------|---------------|-------------------|------------------|-------------------|
| PGx variant lookup | Yes | Yes | Yes | Yes | Yes | Yes |
| Variant similarity search (HNSW) | Yes (<1M vectors) | Yes (<100K) | Yes (billions, sharded) | Yes (<10M via NAPI) | Yes (<1M) | Yes (<100K) |
| Variant calling (SNP/indel) | Limited (small panels) | No | Yes (full WGS) | Yes (full WGS via NAPI) | Limited (targeted) | No |
| SV/CNV detection | On-demand module | No | Yes | Yes (via NAPI) | No | No |
| STR expansion screening | Yes | Limited | Yes | Yes | Yes | No |
| Protein structure search | Yes | On-demand | Yes | Yes (via NAPI) | No | No |
| Phylogenetic search (hyperbolic) | Yes | On-demand | Yes | Yes | No | No |
| Quantum error modeling | Yes (10 qubits) | No | Yes (25 qubits) | Yes (15 qubits) | No | No |
| LLM report generation | On-demand | No | Yes | Yes | No | No |
| Full WGS pipeline (30x) | No (memory limit) | No | Yes (<60s with NAPI) | Yes (<120s with NAPI) | No | No |
| Offline operation | Yes (service worker) | Yes (service worker) | Yes (local) | Yes (bundled) | Yes (standalone) | No |
| Privacy (no data upload) | Yes | Yes | Configurable | Yes | Yes | Partial (data in Worker) |

### Memory and Compute Requirements by Tier

| Tier | Min RAM | Min Storage | CPU Requirement | Network | Primary Use Case |
|------|---------|-------------|-----------------|---------|------------------|
| Browser (desktop) | 4GB system (2GB tab) | 1.5GB (IndexedDB) | Any modern (2018+) | Initial load only | Clinical variant interpretation |
| Browser (mobile) | 4GB system (512MB tab) | 250MB (IndexedDB) | Snapdragon 680+ / A13+ | Initial load only | PGx screening |
| Node.js server | 32GB+ | 50GB+ | 16+ cores, AVX2/NEON | LAN/WAN | Batch genomic analysis |
| Electron (field) | 16GB | 10GB | 8+ cores | None required | Field pathogen identification |
| wasmtime (space) | 4GB | 500MB | ARM Cortex-A72 equiv. | None | Crew health monitoring |
| Cloudflare Worker | 128MB (Worker limit) | 1GB (KV) | Shared (edge) | Always-on | Stateless API queries |

---

## Performance Comparison: WASM Genomics vs. Traditional Tools

### Variant Similarity Search

| Tool/Approach | 100K Variants (k=10) | 1M Variants (k=10) | Deployment | Notes |
|--------------|---------------------|---------------------|------------|-------|
| Linear scan (JavaScript) | 850ms | 8.5s | Browser | O(n) scan, no index |
| BLAST (server) | 12s + network | 120s + network | Cloud | Alignment-based, not ANN |
| `micro-hnsw-wasm` (browser) | 0.5ms | N/A (32 vec/core limit) | Browser | 11.8KB, ultra-fast for small sets |
| `ruvector-wasm` (browser) | 0.25ms | 0.55ms | Browser | Full HNSW, WASM SIMD |
| `ruvector-node` (NAPI) | 0.20ms | 0.40ms | Node.js | Native SIMD (AVX2) |
| **Speedup (ruvector-wasm vs. linear)** | **3,400x** | **15,450x** | | |

### Pharmacogenomic Screening Latency

| Approach | Latency | Requires Internet | Privacy |
|----------|---------|-------------------|---------|
| Cloud API (e.g., PharmCAT) | 2-30s (network + compute) | Yes | Data uploaded |
| Local server (GATK + PharmCAT) | 5-15s | No | Local only |
| `micro-hnsw-wasm` (browser) | <1ms (variant match) + ~50ms (phenotype) | No (after first load) | Client-side only |
| **Improvement** | **100-600x faster, fully private** | | |

### Field Pathogen Identification

| Approach | Time to First ID | Internet Required | Hardware |
|----------|-----------------|-------------------|----------|
| BLAST against NCBI (remote) | 30s-5min (network) | Yes | Any + internet |
| Kraken2 (local, prebuilt DB) | 2-5s per batch | No | 32GB RAM server |
| RuVector Electron (NAPI + WASM) | <2s per batch | No | 16GB laptop |
| RuVector wasmtime (space) | <5s per batch | No | 4GB flight computer |

---

## Security Considerations

### WASM Sandbox Guarantees

| Threat | WASM Mitigation |
|--------|-----------------|
| Buffer overflow in genomic parser | WASM linear memory is bounds-checked; out-of-bounds access traps |
| Malicious WASM module injection | Subresource Integrity (SRI) hashes on all WASM assets; Content Security Policy restricts origins |
| Side-channel timing attacks on variant data | WASM execution is sandboxed; no access to host timers at nanosecond resolution (Performance.now() resolution is reduced) |
| Genomic data exfiltration | CSP `connect-src` restricts network access after module load; service worker blocks unauthorized requests |
| Tampering with reference data | SHA-256 checksums verified by `ruvector-delta-wasm` on every update; cryptographic signatures via `ruvector-wasm` kernel pack system (ADR-005) |

### Clinical Validation in WASM

WASM's deterministic execution model enables a unique validation strategy for clinical genomics software:

1. **Bit-exact reproducibility.** The same WASM binary, given the same inputs, produces identical outputs on any compliant runtime (V8, SpiderMonkey, wasmtime, Wasmer). This eliminates platform-dependent floating-point differences that plague native clinical software validation.

2. **Frozen binary validation.** A validated WASM binary can be cryptographically hashed and the hash recorded in a validation report. Any future execution of the same binary against the same test inputs must produce the same hash of outputs. This satisfies FDA 21 CFR Part 11 requirements for electronic records.

3. **Sandboxed execution for IVD.** For in-vitro diagnostic (IVD) software classification under EU IVDR, the WASM sandbox provides documented memory isolation, preventing a genomic analysis module from interfering with other device software.

---

## Alternatives Considered

### Alternative 1: Native Executables per Platform

Compile the DNA analyzer to separate native binaries for Linux, macOS, Windows, iOS, and Android.

**Rejected because:**
- 5 separate build targets with platform-specific CI/CD
- No browser deployment (requires installation, IT approval)
- No privacy guarantee (native binaries can access file system, network)
- Non-deterministic floating-point across platforms complicates clinical validation
- Cannot run on Cloudflare Workers, Deno Deploy, or other WASM-native edge platforms

### Alternative 2: JavaScript-Only Implementation

Implement all genomic algorithms in pure JavaScript/TypeScript.

**Rejected because:**
- 3-10x slower than WASM for compute-intensive operations (attention, GNN, HNSW)
- No SIMD support (WASM SIMD128 provides 2-4x speedup for vectorized operations)
- JavaScript garbage collection causes unpredictable latency spikes during clinical analysis
- Cannot reuse existing Rust crate ecosystem (23 WASM crates already available)
- Memory overhead: JavaScript objects require 2-5x more memory than WASM linear memory for the same data

### Alternative 3: Docker Containers for Edge Deployment

Package the DNA analyzer as a Docker container for edge deployment.

**Rejected because:**
- Docker requires Linux kernel or VM; cannot run in browser
- Container startup time (seconds) unacceptable for interactive clinical workflows
- Docker images are 100MB-1GB; WASM total is <10MB
- No support for mobile devices or Cloudflare Workers
- Does not provide browser-based privacy guarantees (container has full network access)

### Alternative 4: asm.js as WASM Fallback

Use asm.js for browsers that do not support WASM.

**Rejected because:**
- WASM support is universal in modern browsers (Chrome 57+, Firefox 52+, Safari 11+, Edge 16+) -- covers >98% of global browser traffic
- asm.js is 2-3x slower than WASM for the same code
- asm.js modules are 2-5x larger than WASM (text-based vs. binary format)
- No streaming compilation support
- Maintaining dual asm.js/WASM builds doubles the testing surface

---

## Consequences

### Benefits

1. **Universal deployment from a single codebase.** The same Rust source compiles to 23 WASM crates that run identically in browsers, Node.js, Cloudflare Workers, wasmtime, Wasmer, and embedded systems. No per-platform build targets.

2. **Democratized access to genomic analysis.** A pharmacogenomic screen that previously required a bioinformatics server now runs in a mobile browser on a $150 smartphone with no internet after first load.

3. **Privacy by architecture.** Client-side WASM execution keeps patient genomic data in the browser. GDPR, HIPAA, and similar regulations are satisfied by construction, not by policy.

4. **Space-ready genomics.** The DNA analyzer can operate autonomously on a Mars transit vehicle with no Earth uplink, using <1MB of WASM binaries and <300MB of RAM.

5. **Sub-second interactive analysis.** Progressive loading delivers first results (PGx lookup) in <500ms on desktop, <2s on mobile, compared to 5-30s for server round-trip approaches.

6. **Deterministic clinical validation.** WASM's bit-exact reproducibility simplifies FDA and EU IVDR validation of clinical genomic software.

7. **Bandwidth efficiency.** Delta updates via `ruvector-delta-wasm` reduce monthly reference data updates from ~150MB to ~8MB, critical for low-bandwidth deployments.

### Risks and Mitigations

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| WASM memory limit (4GB) insufficient for full WGS | High | High | Full WGS analysis is server-mode only (Node.js NAPI); WASM tier handles panels and targeted analysis |
| WASM SIMD support varies across browsers | Medium | Medium | Feature detection at startup; graceful fallback to scalar code paths (1.5-2x slower) |
| Service worker cache eviction loses offline data | Medium | Medium | Persistent storage request via `navigator.storage.persist()`; warn user if not granted |
| Reference data becomes stale in offline deployments | Medium | Low | Delta update check on every network reconnection; version watermark on all reports |
| Module loading latency on slow networks | Medium | Medium | Foundation layer (110KB) loads first; Stage 3+ modules only on demand; Brotli compression |
| Browser OOM on mobile during advanced analysis | Medium | High | Memory pressure monitoring; automatic module eviction; clear user warning before loading heavy modules |
| WASM floating-point determinism edge cases (NaN canonicalization) | Low | High | All genomic operations use `f32` with explicit NaN handling; test suite verifies bit-exact output across V8, SpiderMonkey, and wasmtime |
| Cryptographic signing of WASM modules compromised | Low | Critical | Ed25519 signatures via `ruvector-wasm` kernel pack system; key rotation policy; SRI hashes as secondary check |

---

## References

1. Haas, A., et al. (2017). "Bringing the web up to speed with WebAssembly." *PLDI 2017*, 185-200. (WASM specification and performance model.)

2. Jangda, A., et al. (2019). "Not so fast: Analyzing the performance of WebAssembly vs. native code." *USENIX ATC 2019*. (WASM vs. native performance ratios: 0.55x-0.95x.)

3. Castro, S.L., et al. (2016). "Nanopore DNA sequencing and genome assembly on the International Space Station." *Scientific Reports*, 7, 18022. (MinION sequencing aboard ISS.)

4. Burton, A.S., et al. (2020). "Off Earth identification of bacterial populations using 16S rDNA nanopore sequencing." *Genes*, 11(1), 76. (Space-based genomic analysis.)

5. WebAssembly SIMD Specification. https://github.com/WebAssembly/simd (WASM SIMD128 proposal, shipped in all major browsers.)

6. WebAssembly Threads and Atomics. https://github.com/WebAssembly/threads (SharedArrayBuffer for parallel WASM execution.)

7. Service Workers API. MDN Web Docs. https://developer.mozilla.org/en-US/docs/Web/API/Service_Worker_API (Offline caching strategy.)

8. RuVector Core Architecture. ADR-001. (HNSW, SIMD, quantization foundations.)

9. RuVector WASM Runtime Integration. ADR-005. (Kernel pack system, WASM security model.)

10. RuVector Genomic Vector Index. ADR-003. (HNSW genomic vector spaces, distance metrics.)

11. RuVector Variant Calling Pipeline. ADR-009. (Multi-modal ensemble architecture.)

12. FDA Guidance on Computer Software Assurance for Production and Quality System Software. (Deterministic execution requirements.)

13. EU In Vitro Diagnostic Regulation (IVDR) 2017/746. (Software as medical device classification.)

14. Nickel, M., & Kiela, D. (2017). "Poincare embeddings for learning hierarchical representations." *NeurIPS*, 6338-6347. (Hyperbolic embeddings for phylogenetic data.)

15. Malkov, Y., & Yashunin, D. (2018). "Efficient and robust approximate nearest neighbor search using Hierarchical Navigable Small World graphs." *IEEE TPAMI*, 42(4), 824-836. (HNSW algorithm.)

---

## Related Decisions

- **ADR-001**: RuVector Core Architecture (HNSW index, SIMD optimization, quantization)
- **ADR-003**: Genomic Vector Index (multi-resolution HNSW for genomic data)
- **ADR-005**: WASM Runtime Integration (kernel pack system, security model)
- **ADR-009**: Variant Calling Pipeline (ensemble architecture, streaming analysis)

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 0.1 | 2026-02-11 | RuVector Genomics Architecture Team | Initial architecture proposal |
