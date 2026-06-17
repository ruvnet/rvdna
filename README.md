# rvdna

Quantum/AI genomics engine and the `.rvdna` cognitive container format, in Rust + WebAssembly.

> Extracted from [ruvnet/ruvector](https://github.com/ruvnet/ruvector) per ADR-257.
> Depends on the published `ruvector-*` crates (crates.io). Builds standalone:
> `cargo build`. The npm wrapper lives in `npm/packages/rvdna` (`@ruvector/rvdna`).

## Layout
- `examples/dna` — the `rvdna` Rust crate (lib + `rvdna-cli`)
- `npm/packages/rvdna` — `@ruvector/rvdna` NAPI wrapper

## License
MIT © Ruvector Team
