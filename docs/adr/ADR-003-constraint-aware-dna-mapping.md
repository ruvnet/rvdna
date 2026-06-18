# ADR-003: Homopolymer-free, GC-balanced byte<->DNA mapping

- Status: accepted
- Date: 2026-06-18

## Context

Mapping bytes to a 4-letter DNA alphabet is the layer that meets biology. Real
synthesis and sequencing degrade badly on long homopolymer runs (e.g.
`AAAAAA`) and on sequences whose GC content drifts far from 50%. The mapping
must therefore be constraint-aware, not just information-dense.

## Decision

Reject naive 2-bit packing (A=00, C=01, G=10, T=11): it produces long
homopolymer runs and uneven GC that hurt real synthesis/sequencing.

Instead use a reversible base-3 ("trit") rotating transform (Goldman-style):

- each byte -> 6 ternary digits
- each trit picks the next base as `cur = (prev + 1 + trit) % 4` over the
  alphabet A=0, C=1, G=2, T=3

Because `cur` is always at least `prev + 1 (mod 4)`, the emitted base ALWAYS
differs from the previous one. This guarantees a maximum homopolymer run of 1,
and GC stays near 50% statistically.

Cost: ~6 bases/byte (about 1.33 bits/base) versus the 2 bits/base theoretical
max. This is a deliberate density-for-robustness trade.

## Consequences

- Homopolymer runs are eliminated by construction (max run = 1), and GC balance
  is maintained without an explicit balancing pass.
- The running-`prev` chain introduces mild 2-byte error-locality: a corrupted
  base can perturb the decoding of the adjacent byte. The inner RS code absorbs
  this localized damage.
- Density is ~33% below the 2 bits/base ceiling — an accepted cost for
  synthesis/sequencing robustness.
