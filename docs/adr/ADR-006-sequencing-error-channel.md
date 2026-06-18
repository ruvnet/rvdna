# ADR-006: Synthetic sequencing-error channel model

- Status: accepted
- Date: 2026-06-18

## Context

To validate the codes without a wet lab, we need a channel that injects the
error modes real DNA storage exhibits: per-base substitutions, insertions,
deletions, whole-strand dropout, and variable sequencing coverage. The channel
must be reproducible so demos and tests are deterministic.

## Decision

Provide a configurable `ErrorModel` with:

- per-base substitution / insertion / deletion probabilities
- whole-strand dropout probability
- sequencing coverage (reads per surviving strand)

The channel is deterministic given a seed, for reproducible demos and tests.
This lets the demo "mutate" DNA live and show the codes recovering it.

## Consequences

- Demos and tests are reproducible: the same seed yields the same corrupted read
  pool.
- Indels are modelled explicitly via insertion/deletion probabilities, so the
  hardest failure mode is exercised rather than assumed away.
- Coverage > 1 produces multiple noisy reads per surviving strand, which feeds
  the downstream consensus stage; coverage = 1 disables consensus's averaging
  benefit and stresses the inner code.
