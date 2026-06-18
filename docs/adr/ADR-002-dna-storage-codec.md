# ADR-002: DNA Storage Codec Simulator: architecture

- Status: accepted
- Date: 2026-06-18

## Context

DNA data storage offers extreme density and millennia-scale archival, but real
systems face synthesis/sequencing errors (substitutions and indels), strand
dropout, and biological sequence constraints (homopolymer runs, GC imbalance).
We want an end-to-end, wet-lab-free SIMULATOR inside the existing rvdna crate so
the full encode/channel/decode loop can be exercised, demoed, and tested without
any physical synthesis or sequencing.

## Decision

Add a new module `examples/dna/src/storage/` implementing a layered pipeline.

Encode path:

- file bytes -> CRC32 + manifest
- split into K source blocks
- Fountain (LT) OUTER rateless erasure code produces droplets
- per-strand header `[index|seed|degree]`
- Reed-Solomon (GF(256)) INNER code adds parity
- constraint-aware DNA mapping (homopolymer-free, GC-balanced)
- `DnaArchive { manifest, Vec<Strand> }`

Decode path reverses it:

- read pool -> cluster + majority consensus
- RS decode per strand
- Fountain peel
- reassemble -> CRC32 verify

Submodules: `constraints`, `gf256` (Reed-Solomon), `fountain` (LT), `channel`
(error injection), `consensus`, plus the `DnaStorageCodec` orchestrator and the
`EncodeParams` / `DnaArchive` / `DecodeReport` types in `mod.rs`. The module
reuses the crate's `DnaError`. The codec is exposed via a CLI subcommand and the
browser visualizer.

## Consequences

- The feature is additive, pure-Rust, and WASM-friendly; nothing in the existing
  crate depends on it.
- The two-layer code mirrors real DNA-storage systems (e.g., DNA Fountain):
  fountain handles erasures and whole-strand loss, RS handles in-strand
  substitutions.
- Indels are the hardest case. They are mitigated, not perfectly solved, via
  consensus over coverage plus strand over-provisioning, rather than by an
  alignment-based indel-correcting code.
