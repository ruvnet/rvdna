# ADR-016: Archaic introgression and ghost-lineage recovery (`TRACE-rv`)

- Status: accepted
- Date: 2026-07-31

## Context

On 30 July 2026 Zhang, Biddanda, Johnson, O'Dushlaine and Moorjani published
*Recovering signatures of archaic hominin introgression using ancestral
recombination graphs* in **Science**. Their method, **TRACE** (TRacking Archaic
Contributions via ARG Estimation), reconstructs ancestral recombination graphs
across present-day genomes and finds regions whose ancestry reaches unusually
far back in time. It recovered two previously unknown "ghost" hominin lineages:

- one that diverged from the modern human lineage ~800 ka, interbred inside
  Africa more than 50 ka, and contributes ~0.5–1% of every living person's
  genome;
- one ~1.8 Ma "super-archaic" lineage that interbred with **Denisovans**, and so
  reaches modern genomes only second-hand, concentrated in Oceania.

Neither lineage has ever been sequenced. The whole result rests on the shape of
genealogies among the living.

rvDNA already has the primitives this kind of analysis needs — k-mer profile
vectors, an HNSW vector index via `ruvector-core`, sequence alignment, and the
`.rvdna` container — but had no module that used them for population genetics.
There was also no worked example of the metaharness **Darwin mode** and
**flywheel** patterns operating on a real rvDNA workload.

## Decision

Add `examples/dna/src/archaic.rs`: a self-contained, first-principles
implementation of the *logic* behind ARG-based introgression detection, built on
rvDNA and validated against simulated ground truth.

### Structure

| Piece | Responsibility |
|---|---|
| `build_population_tree` | Hominin demography calibrated to the published divergence times (1.8 Ma / 800 ka / 650 ka / 400 ka / 200 ka / 60 ka). |
| `simulate_cohort` | Structured (multispecies) coalescent per window, Poisson mutations dropped on the resulting genealogy, introgression planted by reassigning a lineage to the donor deme. |
| `TraceEngine` | The detector. Embeds each window as an rvDNA k-mer vector, indexes it in a RuVector HNSW graph, retrieves candidates, computes exact Jukes–Cantor-corrected divergence against those only, and converts it to coalescent depth. |
| `TraceEngine::cluster_ghosts` | Two-means on log depth, with a separation test, to split unattributed deep segments into lineages. |
| `darwin` | Hyperparameter hill-climb: mutate one gene, evaluate on a train split, accept only past an explicit margin, report on a held-out split. |
| `flywheel` | Confirmed ghost segments become references for the next round; runs until two consecutive dry rounds. |

Two binaries drive it: `trace-rv` (the full study) and `depth-hist` (the
brute-force reference measurement the indexed detector is judged against).

### Why the division of labour between HNSW and exact divergence

Coalescent depth is the quantity that decides every call, and it must be exact —
an approximate distance in k-mer space is not a genealogy. But comparing every
segment against every other is quadratic and mostly wasted.

Because independent windows have independent genealogies, segments from the same
window cluster tightly in k-mer space and segments from different windows do not.
That makes approximate nearest-neighbour search a genealogical shortcut rather
than a heuristic: HNSW decides *what to look at*, exact divergence decides *the
answer*. The engine counts both the comparisons it ran and the all-pairs
equivalent, so the saving is a measured number rather than a claim.

### Why the flywheel needs a promotion gate

The first implementation let any GHOST call become a reference for the next
round, and rescued a segment whenever it fell within `ghost_link_ka` of one.
That collapsed: round 1 scored F1 0.84, round 2 dropped to 0.22 and round 3 to
0.19, as false positives were promoted to references and then validated more
false positives.

Two changes fix it, and both follow from the population genetics rather than from
tuning:

1. **Only primary calls may be promoted.** A segment that cleared the depth test
   on its own evidence can become a reference; a segment that was rescued by
   resembling a reference never can. Without this the loop compounds its own
   errors.
2. **The rescue rule is relative, not absolute.** Two copies of the same
   introgressed haplotype coalesce with each other *inside the donor population*
   — recently, in coalescent terms — while either coalesces with a living person
   only back past the donor's split. So the test is "closer to a known ghost than
   to anyone alive", not "close to a known ghost".

`ArchaicCall::primary` records which test a call passed, so the distinction is
visible in the output rather than implicit.

### Why divergence times are quoted from the lower tail

A segment's coalescent depth is always *older* than the split that produced it,
because the donor lineage had its own history to get through first. The mean
depth therefore overestimates a divergence time, sometimes badly. `GhostCluster`
reports `p10_split_ka` alongside the mean and median, and the report quotes the
p10 as the divergence-time estimate.

## Consequences

- The cohort is **simulated**, not real 1000 Genomes / HGDP data. The module
  reproduces the reasoning of the published method and validates it against
  known truth; it is not a reimplementation of TRACE and makes no claim about any
  living person's ancestry. `FinalReport::caveats` ships this in the output
  itself so it travels with the numbers.
- The per-base mutation rate is rescaled (~109×) so a 2 kb simulated window
  carries the segregating-site information of ~217 kb of real human sequence.
  All reported times use the true human clock, because the same rate is used in
  both directions. `trace-rv` verifies this before doing anything else and
  reports the correlation between true and estimated TMRCA.
- `archaic.rs` adds no new dependencies — `rand`, `rand_distr`, `serde` and
  `ruvector-core` were already in the manifest.
- Darwin's search dominates runtime (~8 min for 20 proposals over 240 windows),
  because most proposals change the embedding and force an HNSW rebuild. The
  index cache is bounded at three entries to keep memory flat; raising it trades
  memory for speed.
- Everything is seeded. The same checkout produces the same numbers, which is
  what makes the generated report and the visual story auditable.

## Addendum — real genomes (`real-dna`)

The original decision record covers a simulated cohort only. A later stage adds
`src/bin/real_dna.rs`, which runs on real data and makes two decisions worth
recording.

### Why the panel is loaded as bitsets and compared exhaustively

The real-data stage does not use HNSW. The whole point of the index is to avoid
comparisons that cannot matter, and the fidelity sweep already establishes it
costs nothing to use it — so on real data the more informative choice is to do
the work exactly and report the cost. 800 haplotypes over 19 windows is
6,072,400 exact pairwise comparisons, and by popcount over bit-packed haplotypes
that runs in 0.4 s. There is nothing to approximate away at this scale, and an
exact number is easier to defend.

Subsampling 400 of the 2,504 phase 3 samples does not bias a pairwise
comparison: any site at which two included haplotypes differ is by definition
polymorphic and therefore present in the VCF.

### Why the introgression claim rests on D and not on divergence

The natural comparison is divergence from each modern haplotype to a sequenced
archaic, contrasted between Africans and non-Africans. It is not sound. GRCh37 is
a European-weighted reference and the archaic genotypes were called against it,
so African haplotypes accumulate apparent mismatches to *any* archaic genome for
reasons unrelated to admixture. Running it gives a directionally correct,
uninterpretable number.

`real-dna` computes it anyway — as `archaic_affinity`, flagged in both the doc
comment and the report caveats — and rests the claim on Patterson's D instead.
ABBA and BABA both condition on the archaic carrying the derived allele, so bias
in the archaic call set is shared between the two counts and cancels in the
ratio.

The standard error is a leave-one-window-out block jackknife rather than a
binomial, because sites inside a 50 kb window are linked and are not independent
draws. With 19 blocks the resulting Z is 1.1 — the effect is correctly signed and
not significant, and the report says so in those words. The extrapolated
sequence length needed for Z = 3 (~7 Mb) is emitted as `mb_for_z3` so the limit
is a number in the output rather than a hedge in the prose.

### Consequences

- Archaic genotypes are the producers' raw snpAD calls with no extra GQ/DP
  filter. Damage and reference bias inflate absolute divergence; neither acts
  differently on African and non-African haplotypes, which is why the contrast
  and not the absolute value carries the result.
- Per-base divergence to an archaic needs a callable-base denominator, but the
  archaic table stores only variable sites. Total callable bases come from the
  producers' own per-individual count, apportioned across windows in proportion
  to local callable-site density. This affects only `archaic_affinity`; D is
  a ratio of site-pattern counts and does not use it.
- `data/real/` is checked in (~35 MB) so the stage is reproducible without
  network access. Re-fetching needs outbound HTTPS to EBI and MPI-EVA.
