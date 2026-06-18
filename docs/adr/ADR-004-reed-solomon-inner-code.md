# ADR-004: Reed-Solomon (GF(256)) inner error-correction

- Status: accepted
- Date: 2026-06-18

## Context

After consensus, the dominant residual error within a strand is byte-level
substitution. The inner code must correct these substitutions on strand-sized
codewords before the outer fountain layer attempts to peel.

## Decision

Use systematic Reed-Solomon over GF(256):

- primitive polynomial 0x11D, generator α=2
- `nsym` parity bytes per strand correct up to `nsym/2` substitution errors
- decoder pipeline: Berlekamp-Massey (error locator) + Chien search (root
  finding) + Forney (error magnitudes)
- operates on codewords of at most 255 bytes (strand-sized)

Why GF(256): symbols are byte-aligned, the algorithm is mature and well
understood, and it is strong against the substitution-dominated residual that
remains after consensus.

## Consequences

- Each strand carries `nsym` bytes of parity overhead; correction capacity is
  `nsym/2` substituted bytes per strand, tunable via `EncodeParams`.
- Codewords are capped at 255 bytes, which sets the strand payload size and
  bounds per-strand work.
- RS corrects substitutions but not indels; surviving indels are handled
  upstream by consensus and absorbed as localized substitution-like damage where
  possible. Strands RS cannot fix are dropped and left to the fountain layer.
