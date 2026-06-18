# ADR-008: Browser visualizer & WASM strategy

- Status: accepted
- Date: 2026-06-18

## Context

The codec is most convincing when seen running: encode a real file to DNA,
mutate it, and watch the error correction recover it. We want that demo to work
with zero setup, while keeping the Rust crate as the single source of truth for
the codec.

## Decision

Ship a dependency-free single-page web demo (`web/`) that runs by opening
`index.html`:

- drop a file/PNG, encode to DNA, mutate via sliders, watch error correction
  recover it live

The demo reimplements the codec's CONCEPTS in pure JS for an instant, build-free
experience. The authoritative codec remains the Rust crate, which is
`wasm32`-ready, so a future build can swap the JS core for the compiled WASM with
no UI change.

This supports the "I stored my GitHub repo in DNA" hook.

## Consequences

- The demo opens with no build step, server, or dependencies — maximal
  approachability.
- There are two implementations to keep conceptually aligned (JS demo core and
  Rust codec); the JS core is illustrative, and the Rust crate is authoritative
  for correctness.
- The `wasm32` path is deliberately left open so the JS core can later be
  replaced by the compiled WASM behind the same UI.
