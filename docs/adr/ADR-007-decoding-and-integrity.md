# ADR-007: Decoding pipeline & integrity verification

- Status: accepted
- Date: 2026-06-18

## Context

The decoder receives a noisy read pool — multiple, possibly corrupted reads per
surviving strand, with some strands missing entirely. It must turn this into the
original file bytes and report, honestly, whether recovery succeeded.

## Decision

Decode is a fixed pipeline:

- cluster noisy reads by similarity
- per-cluster majority-vote consensus (cancels random substitutions when
  coverage > 1)
- demap DNA -> bytes
- RS-decode each strand, recovering the true `[index|seed|degree]` header and
  payload, dropping strands RS cannot fix
- Fountain peeling to recover the K source blocks
- reassemble -> CRC32 check

A `DecodeReport` surfaces `strands_recovered` / `blocks_recovered` / `crc_ok`
for the visualizer's live stats.

## Consequences

- Consensus before RS lowers the per-strand error rate the inner code must
  absorb, raising the effective survivable noise.
- The report makes failure observable: a passing CRC confirms byte-exact
  recovery, while the recovered-counts expose how much margin remained.
- Clustering by similarity is the indel-tolerant step; misclustered reads
  degrade consensus rather than hard-failing, leaving residual damage to RS and
  the fountain layer.
