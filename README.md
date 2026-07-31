# rvDNA — AI-Native Genomic Analysis in Pure Rust + WebAssembly

[![Rust](https://img.shields.io/badge/Rust-2021-orange?logo=rust)](https://www.rust-lang.org)
[![npm](https://img.shields.io/badge/npm-%40ruvector%2Frvdna-cb3837?logo=npm)](https://www.npmjs.com/package/@ruvector/rvdna)
[![WebAssembly](https://img.shields.io/badge/WebAssembly-ready-654ff0?logo=webassembly&logoColor=white)](#webassembly--npm)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](#license)

**`rvDNA`** is an **AI-native genomics engine** in pure Rust (with WebAssembly + Node bindings):
biomarker risk scoring, 23andMe genotyping, pharmacogenomics, variant calling, protein prediction,
epigenomics, and **HNSW vector search** over genomic profiles — plus the **`.rvdna`** AI-native
cognitive-container format.

> 20-SNP biomarker risk scoring · streaming anomaly detection · 64-dim profile vectors · 23andMe
> genotyping · CYP2D6/CYP2C19 pharmacogenomics · variant calling · protein prediction · HNSW
> similarity search — no Python, runs on the edge via WASM.

## Why rvDNA

- **Pure Rust, edge-ready** — memory-safe, fast, compiles to native and `wasm32`; runs in the
  browser, on Node, or on-device.
- **Precision-medicine toolkit** — biomarker risk scores, pharmacogenomic star-alleles
  (**CYP2D6 / CYP2C19**), and 23andMe-style genotyping out of the box.
- **Whole pipeline** — alignment (Smith–Waterman), k-mer analysis, variant calling, protein
  prediction, and epigenomic/temporal modeling.
- **Vector-native** — 64-dimensional genomic profile vectors with **HNSW** approximate
  nearest-neighbor search for similarity and cohort matching.
- **`.rvdna` format** — an AI-native cognitive container for genomic data and derived models.

## Capabilities

| Module | Capability |
|--------|-----------|
| `biomarker`, `biomarker_stream` | 20-SNP **biomarker risk scoring** + streaming anomaly detection |
| `genotyping` | **23andMe** genotype parsing & calling |
| `pharma` | **pharmacogenomics** — CYP2D6 / CYP2C19 star-alleles, drug-response |
| `variant` | **variant calling** pipeline |
| `protein` | protein structure / function prediction |
| `alignment` | Smith–Waterman sequence **alignment** |
| `kmer`, `kmer_pagerank` | k-mer analysis & PageRank over genomic graphs |
| `epigenomics` | temporal **epigenomic** modeling |
| `health` | health & biomarker reporting |
| `archaic` | **archaic introgression & ghost-lineage recovery** — structured coalescent, HNSW-accelerated genealogy search, Darwin-mode tuning, flywheel refinement |

## Install

```bash
# Rust
cargo add rvdna

# Node / browser (NAPI + WASM wrapper)
npm install @ruvector/rvdna
```

## Quick start (Rust)

```rust
use rvdna::prelude::*;

// Score 20-SNP biomarker risk from a genotype profile.
let profile = GenomicProfile::from_23andme("genome.txt")?;
let risk = profile.biomarker_risk();
println!("risk vector (64-dim): {:?}", risk.vector());
```

## WebAssembly & npm

```js
import rvdna from "@ruvector/rvdna";
// biomarker scoring, genotyping and .rvdna I/O — in Node or the browser
```

## Build

```bash
cargo build --release                 # native lib + rvdna-cli
cargo test                            # test suite
cd npm/packages/rvdna && npm run build:napi   # native/WASM bindings
```

## Use cases

Consumer genomics & 23andMe analysis · precision medicine & pharmacogenomics · biomarker
risk scoring · variant calling pipelines · protein prediction · genomic similarity / cohort
search with HNSW · edge & in-browser genomics.

## Ghost lineages — a worked study

`examples/dna/src/archaic.rs` rebuilds, from first principles, the logic behind
[TRACE](https://news.berkeley.edu/2026/07/30/new-technique-pinpoints-human-dna-inherited-from-ghost-ancestors/)
(*Science*, 30 July 2026), which recovered two extinct hominin lineages from
present-day genomes alone. It simulates a hominin cohort under a structured
coalescent, finds introgressed segments by their coalescent depth using a
RuVector HNSW index for candidate retrieval, tunes itself with metaharness
**Darwin mode**, and refines with a **flywheel** that turns each confirmed ghost
segment into evidence for the next.

```bash
cargo build --release
./target/release/depth-hist artifacts/data   # the brute-force reference measurement
./target/release/trace-rv   artifacts/data   # simulate, evolve, spin, cluster
node artifacts/build.mjs                     # render the illustrated report
```

The rendered study, with animated diagrams, lives in [`artifacts/`](artifacts).
See [ADR-016](examples/dna/adr/ADR-016-archaic-introgression-ghost-lineages.md).

> The cohort is simulated against known truth, not real 1000 Genomes data. It
> validates the method; it makes no claim about any living person's ancestry.

## Layout

- [`examples/dna`](examples/dna) — the `rvdna` Rust crate (library + `rvdna-cli`)
- [`npm/packages/rvdna`](npm/packages/rvdna) — `@ruvector/rvdna` NAPI/WASM wrapper
- [`artifacts`](artifacts) — the ghost-lineage study: animated diagrams, generated data, illustrated report

## License

MIT © Ruvector Team. Part of the [ruvector](https://github.com/ruvnet/ruvector) ecosystem
(extracted per ADR-257). Built on the published `ruvector-*` crates (HNSW, attention, GNN, solver).

> Research/educational software — **not** a medical device and not for clinical diagnosis.
