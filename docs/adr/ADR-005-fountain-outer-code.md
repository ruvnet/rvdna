# ADR-005: Luby-Transform rateless outer erasure code

- Status: accepted
- Date: 2026-06-18

## Context

Whole strands are lost in DNA storage: synthesis dropout, sequencing that never
covers a given molecule, or strands the inner RS code cannot repair. The outer
layer must recover the original data from an arbitrary subset of surviving
strands, without a fixed code rate chosen up front.

## Decision

Use an LT (Luby-Transform) fountain code as the OUTER layer so the system is
rateless:

- emit as many strands as desired (overhead factor, default ~1.8x)
- recover from arbitrary strand DROPOUT once roughly `K(1+ε)` good droplets
  arrive
- Robust Soliton degree distribution
- deterministic seed-driven neighbour selection (SplitMix64) shared by the
  encoder and the peeling decoder, so the decoder can reconstruct each droplet's
  neighbour set from its seed alone

This mirrors Erlich & Zielinski's "DNA Fountain."

## Consequences

- The system is rateless: more robustness is bought simply by emitting more
  strands, with no re-encode of the source blocks.
- Trade-off: a small overhead above K droplets, plus a tiny seed field in each
  strand header (`[index|seed|degree]`).
- Recovery is probabilistic in the number of received droplets; supplying
  meaningfully fewer than `K(1+ε)` good droplets may leave source blocks
  unrecovered.
