//! # Archaic introgression & ghost-lineage recovery (`TRACE-rv`)
//!
//! A first-principles reimplementation, on top of the rvDNA engine and the
//! RuVector HNSW index, of the logic behind **TRACE** (*TRacking Archaic
//! Contributions via ARG Estimation*, Zhang, Biddanda, Moorjani et al.,
//! *Science*, 30 July 2026), which recovered two previously unknown "ghost"
//! hominin lineages from present-day human genomes alone.
//!
//! ## The idea in one sentence
//!
//! If you can reconstruct, window by window, the genealogy relating a set of
//! present-day haplotypes, then a haplotype that entered the population by
//! interbreeding with an extinct lineage betrays itself by *how deep in time it
//! has to go before it meets anybody else* — no ancient DNA required.
//!
//! ## What this module actually does
//!
//! 1. [`sim`] — simulates haplotypes under a **multispecies (structured)
//!    coalescent** on a hominin population tree calibrated to the published
//!    divergence times, with introgression pulses planted at known
//!    (haplotype, window) coordinates. This gives us ground truth.
//! 2. [`TraceEngine`] — the detector. Each window of each haplotype is embedded
//!    as an rvDNA k-mer profile vector and indexed in a RuVector **HNSW**
//!    graph. HNSW retrieves the handful of genealogically plausible relatives
//!    for a query segment in sublinear time; exact Jukes–Cantor-corrected
//!    divergence is then computed only against those, yielding a local
//!    coalescent depth. Segments whose depth to the *entire* modern panel is an
//!    outlier are archaic; those that additionally fail to match any *observed*
//!    archaic reference (Neanderthal, Denisovan) are **ghosts**.
//! 3. [`darwin`] — a Darwin-mode evolutionary search over the detector's own
//!    hyperparameters. Mutate, evaluate in a sandbox on a train split, keep
//!    only what measurably improves, report on a held-out test split.
//! 4. [`flywheel`] — each round's confirmed ghost segments are condensed into
//!    consensus references and written back into the HNSW index, so the next
//!    round can attribute borderline segments by similarity to a ghost that is
//!    no longer invisible. Runs until it goes dry.
//!
//! ## Honest framing
//!
//! The genomes analysed here are **simulated**, not real 1000 Genomes / HGDP
//! data — the module ships a demography calibrated to the published parameter
//! estimates so the pipeline can be validated against known truth. Divergence
//! times are reported on the real human clock; the per-base mutation rate is
//! rescaled (see [`DemographyConfig::mu`]) so a short simulated window carries
//! the same expected number of segregating sites as a realistically sized real
//! one. Nothing here is a clinical or population-genetic claim about any living
//! person.

use std::collections::{BTreeMap, HashMap, HashSet};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Exp, Poisson};
use serde::{Deserialize, Serialize};

use crate::error::{DnaError, Result};
use crate::types::DnaSequence;
use ruvector_core::types::{DbOptions, DistanceMetric, HnswConfig, SearchQuery, VectorEntry};
use ruvector_core::VectorDB;

// ---------------------------------------------------------------------------
// Population tree
// ---------------------------------------------------------------------------

/// Years per human generation, used for every generations <-> years conversion.
pub const GENERATION_YEARS: f64 = 29.0;

/// A deme in the hominin population tree.
#[derive(Debug, Clone)]
pub struct PopNode {
    /// Short label, e.g. `"YRI"` or `"NEA"`.
    pub name: &'static str,
    /// Index of the deme this one merges into, going backwards in time.
    pub parent: Option<usize>,
    /// Generations before present at which this deme merges into `parent`.
    pub split_gens: f64,
    /// Effective (diploid) population size for this branch.
    pub ne: f64,
    /// Broad grouping used for reporting (`"AFR"`, `"OCE"`, `"ARCHAIC"`, ...).
    pub group: &'static str,
}

/// Where a lineage sits in the tree, and how we treat it in the panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Role {
    /// A present-day human haplotype: both queried and usable as panel.
    Modern,
    /// A sequenced archaic genome available as a labelled reference.
    ArchaicReference,
}

/// One simulated haplotype.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Haplotype {
    pub id: String,
    pub population: String,
    pub group: String,
    pub role: Role,
}

/// The archaic sources that can donate segments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ArchaicSource {
    Neanderthal,
    Denisovan,
    /// The ~800 ka lineage that introgressed into modern humans inside Africa.
    GhostDeepAfrican,
    /// The ~1.8 Ma lineage that reached us indirectly, via Denisovans.
    SuperArchaic,
}

impl ArchaicSource {
    pub fn label(&self) -> &'static str {
        match self {
            ArchaicSource::Neanderthal => "Neanderthal",
            ArchaicSource::Denisovan => "Denisovan",
            ArchaicSource::GhostDeepAfrican => "Ghost-A (deep African)",
            ArchaicSource::SuperArchaic => "Ghost-B (super-archaic)",
        }
    }

    /// Is this lineage one we have actually sequenced?
    pub fn is_observed(&self) -> bool {
        matches!(self, ArchaicSource::Neanderthal | ArchaicSource::Denisovan)
    }
}

/// Index of every deme in the tree built by [`build_population_tree`].
#[derive(Debug, Clone)]
pub struct PopulationTree {
    pub nodes: Vec<PopNode>,
    pub by_name: HashMap<&'static str, usize>,
}

impl PopulationTree {
    pub fn idx(&self, name: &str) -> usize {
        *self
            .by_name
            .get(name)
            .unwrap_or_else(|| panic!("unknown deme {name}"))
    }
}

/// Builds the hominin population tree used throughout this module.
///
/// Divergence times follow the estimates reported for the *Science* 2026 ghost
/// lineage study and the archaic genomes literature it builds on:
///
/// | Split | Years BP | Generations |
/// |---|---|---|
/// | Super-archaic lineage | 1,800,000 | 62,069 |
/// | Ghost-A (deep African) | 800,000 | 27,586 |
/// | Neanderthal+Denisovan vs modern | 650,000 | 22,414 |
/// | Neanderthal vs Denisovan | 400,000 | 13,793 |
/// | Deep African structure | 200,000 | 6,897 |
/// | Out-of-Africa | 60,000 | 2,069 |
pub fn build_population_tree() -> PopulationTree {
    let g = |years: f64| years / GENERATION_YEARS;

    // (name, parent, split years BP, Ne, group)
    let spec: Vec<(&'static str, Option<&'static str>, f64, f64, &'static str)> = vec![
        ("ROOT", None, f64::INFINITY, 14_000.0, "ROOT"),
        ("SUPER", Some("ROOT"), 1_800_000.0, 2_000.0, "GHOST"),
        ("ANC_A", Some("ROOT"), 1_800_000.0, 12_000.0, "ANC"),
        ("GHOSTA", Some("ANC_A"), 800_000.0, 3_000.0, "GHOST"),
        ("ANC_B", Some("ANC_A"), 800_000.0, 12_000.0, "ANC"),
        ("ARCH", Some("ANC_B"), 650_000.0, 3_000.0, "ARCHAIC"),
        ("MODERN", Some("ANC_B"), 650_000.0, 12_000.0, "MODERN"),
        ("NEA", Some("ARCH"), 400_000.0, 2_500.0, "ARCHAIC"),
        ("DEN", Some("ARCH"), 400_000.0, 2_500.0, "ARCHAIC"),
        // Modern human structure.
        ("AFR_ANC", Some("MODERN"), 200_000.0, 12_000.0, "AFR"),
        ("SAN", Some("MODERN"), 200_000.0, 7_000.0, "AFR"),
        ("YRI", Some("AFR_ANC"), 90_000.0, 10_000.0, "AFR"),
        ("MSL", Some("AFR_ANC"), 90_000.0, 10_000.0, "AFR"),
        ("LWK", Some("AFR_ANC"), 75_000.0, 10_000.0, "AFR"),
        ("OOA", Some("AFR_ANC"), 60_000.0, 2_500.0, "OOA"),
        ("PNG", Some("OOA"), 48_000.0, 2_000.0, "OCE"),
        ("EURASIA", Some("OOA"), 48_000.0, 4_000.0, "EURASIA"),
        ("SAS", Some("EURASIA"), 42_000.0, 4_000.0, "SAS"),
        ("WEUR", Some("EURASIA"), 42_000.0, 4_000.0, "EURASIA"),
        ("CEU", Some("WEUR"), 35_000.0, 4_000.0, "EUR"),
        ("EAS_ANC", Some("WEUR"), 35_000.0, 4_000.0, "EAS"),
        ("CHB", Some("EAS_ANC"), 26_000.0, 4_000.0, "EAS"),
        ("AMR", Some("EAS_ANC"), 20_000.0, 2_000.0, "AMR"),
    ];

    let mut by_name = HashMap::new();
    for (i, s) in spec.iter().enumerate() {
        by_name.insert(s.0, i);
    }
    let nodes = spec
        .iter()
        .map(|(name, parent, years, ne, group)| PopNode {
            name,
            parent: parent.map(|p| by_name[p]),
            split_gens: if years.is_infinite() {
                f64::INFINITY
            } else {
                g(*years)
            },
            ne: *ne,
            group,
        })
        .collect();

    PopulationTree { nodes, by_name }
}

// ---------------------------------------------------------------------------
// Simulation configuration
// ---------------------------------------------------------------------------

/// Knobs for the simulated cohort.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DemographyConfig {
    /// Length of each independent window, in base pairs.
    pub window_bp: usize,
    /// How many independent windows (loci) to simulate.
    pub n_windows: usize,
    /// Per-base, per-generation mutation rate **as simulated**.
    ///
    /// This is deliberately rescaled upward from the human autosomal rate
    /// (~1.25e-8) so that a short simulated window carries the same expected
    /// number of segregating sites as a much longer real one. With the default
    /// `window_bp = 2000` and `mu = 1.36e-6`, one simulated window is
    /// information-equivalent to roughly 218 kb of real human sequence — a
    /// realistic scale for introgression-tract detection. All *times* reported
    /// by the pipeline are on the true human clock, because the same `mu` is
    /// used in both directions.
    pub mu: f64,
    /// Master RNG seed; the whole pipeline is deterministic given this.
    pub seed: u64,
    /// Introgression pulses to plant, as (source, recipient group, fraction).
    pub pulses: Vec<Pulse>,
}

/// One planted admixture pulse.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pulse {
    pub source: ArchaicSource,
    /// Which reporting group receives it (`"ALL"` means every modern group).
    pub recipient_group: String,
    /// Expected fraction of that group's windows carrying a donated segment.
    pub fraction: f64,
    /// Mean length, in consecutive windows, of each donated tract.
    ///
    /// `1.0` is the original behaviour: every planted segment is one isolated
    /// window, which models a panel of unlinked loci. Anything larger turns the
    /// windows into a **linked chromosome**, where a single interbreeding event
    /// leaves a contiguous run that recombination has not yet broken up. Tract
    /// length is the observable that carries *when* the interbreeding happened,
    /// as opposed to when the lineages split, and it is the signal a
    /// linkage-based detector would exploit.
    #[serde(default = "one")]
    pub tract_mean_windows: f64,
}

fn one() -> f64 {
    1.0
}

/// One planted introgression tract on a linked chromosome.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tract {
    pub haplotype: usize,
    pub start_window: usize,
    pub length: usize,
    pub source: ArchaicSource,
}

impl Default for DemographyConfig {
    fn default() -> Self {
        Self {
            window_bp: 2000,
            n_windows: 420,
            mu: 1.36e-6,
            seed: 0xA5CE_57A1_2026,
            pulses: vec![
                // Ghost-A entered the modern human lineage inside Africa, before
                // the last out-of-Africa migration, so every group carries it.
                Pulse {
                    source: ArchaicSource::GhostDeepAfrican,
                    recipient_group: "ALL".to_string(),
                    fraction: 0.0075,
                    tract_mean_windows: 1.0,
                },
                // Neanderthal into everyone outside Africa.
                Pulse {
                    source: ArchaicSource::Neanderthal,
                    recipient_group: "NONAFR".to_string(),
                    fraction: 0.019,
                    tract_mean_windows: 1.0,
                },
                // Denisovan, concentrated in Oceania.
                Pulse {
                    source: ArchaicSource::Denisovan,
                    recipient_group: "OCE".to_string(),
                    fraction: 0.042,
                    tract_mean_windows: 1.0,
                },
                Pulse {
                    source: ArchaicSource::Denisovan,
                    recipient_group: "EAS".to_string(),
                    fraction: 0.003,
                    tract_mean_windows: 1.0,
                },
                // The super-archaic lineage reached modern humans only through
                // Denisovans, so it is rare and Oceania-biased.
                Pulse {
                    source: ArchaicSource::SuperArchaic,
                    recipient_group: "OCE".to_string(),
                    fraction: 0.0055,
                    tract_mean_windows: 1.0,
                },
                Pulse {
                    source: ArchaicSource::SuperArchaic,
                    recipient_group: "EAS".to_string(),
                    fraction: 0.0006,
                    tract_mean_windows: 1.0,
                },
            ],
        }
    }
}

// ---------------------------------------------------------------------------
// The simulated cohort
// ---------------------------------------------------------------------------

/// A simulated panel of haplotypes with per-window sequences and known truth.
pub struct Cohort {
    pub config: DemographyConfig,
    pub tree: PopulationTree,
    pub haplotypes: Vec<Haplotype>,
    /// `sequences[window][haplotype]` — ASCII `ACGT`, `window_bp` long.
    pub sequences: Vec<Vec<Vec<u8>>>,
    /// Ground truth: `(window, haplotype) -> donating lineage`.
    pub truth: HashMap<(usize, usize), ArchaicSource>,
    /// Ground truth in tract form: one entry per planted contiguous run.
    ///
    /// With the default `tract_mean_windows = 1.0` every tract has
    /// `length == 1`, so this is simply `truth` re-expressed as a list.
    pub tracts: Vec<Tract>,
    /// True TMRCA in generations between every pair, per window — kept only for
    /// the calibration report, not used by the detector.
    pub calibration: Vec<CalibrationPoint>,
}

/// A (true TMRCA, observed divergence) pair used to sanity-check the clock.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalibrationPoint {
    pub true_tmrca_gens: f64,
    pub observed_divergence: f64,
    pub estimated_tmrca_gens: f64,
}

/// Jukes–Cantor correction: raw mismatch fraction -> expected substitutions/site.
pub fn jukes_cantor(p: f64) -> f64 {
    let capped = p.min(0.7449);
    -0.75 * (1.0 - (4.0 / 3.0) * capped).ln()
}

/// Convert a raw mismatch fraction into a TMRCA estimate, in generations.
pub fn divergence_to_tmrca_gens(p: f64, mu: f64) -> f64 {
    jukes_cantor(p) / (2.0 * mu)
}

/// Convert generations to thousands of years.
pub fn gens_to_ka(gens: f64) -> f64 {
    gens * GENERATION_YEARS / 1000.0
}

// ---------------------------------------------------------------------------
// Coalescent machinery
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct TreeNode {
    time: f64,
    children: Option<(usize, usize)>,
}

/// Simulate one window's genealogy under the structured coalescent, then drop
/// mutations on it and read off the leaf sequences.
fn simulate_window(
    tree: &PopulationTree,
    start_demes: &[usize],
    window_bp: usize,
    mu: f64,
    rng: &mut StdRng,
) -> (Vec<Vec<u8>>, Vec<TreeNode>, Vec<usize>) {
    let n = start_demes.len();

    // --- backward-in-time structured coalescent -------------------------
    let mut nodes: Vec<TreeNode> = (0..n)
        .map(|_| TreeNode {
            time: 0.0,
            children: None,
        })
        .collect();

    // deme index -> live lineage node ids
    let mut demes: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (leaf, &d) in start_demes.iter().enumerate() {
        demes.entry(d).or_default().push(leaf);
    }

    // Merge schedule: child deme folds into its parent at `split_gens`.
    let mut merges: Vec<(f64, usize, usize)> = tree
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(i, node)| node.parent.map(|p| (node.split_gens, i, p)))
        .collect();
    merges.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut next_merge = 0usize;

    let mut t = 0.0f64;
    let mut live: usize = n;

    while live > 1 {
        let merge_time = merges.get(next_merge).map(|m| m.0).unwrap_or(f64::INFINITY);

        // Total coalescence rate over all demes at the current time.
        let mut rates: Vec<(usize, f64)> = Vec::new();
        let mut total_rate = 0.0;
        for (&d, lin) in demes.iter() {
            let k = lin.len() as f64;
            if k >= 2.0 {
                let r = k * (k - 1.0) / 2.0 / (2.0 * tree.nodes[d].ne);
                rates.push((d, r));
                total_rate += r;
            }
        }

        let dt = if total_rate > 0.0 {
            Exp::new(total_rate)
                .map(|e| e.sample(rng))
                .unwrap_or(f64::INFINITY)
        } else {
            f64::INFINITY
        };

        if t + dt < merge_time {
            // A coalescence happens first.
            t += dt;
            let mut pick = rng.gen_range(0.0..total_rate);
            let mut chosen = rates[0].0;
            for (d, r) in &rates {
                if pick < *r {
                    chosen = *d;
                    break;
                }
                pick -= *r;
            }
            let lin = demes.get_mut(&chosen).unwrap();
            let i = rng.gen_range(0..lin.len());
            let a = lin.swap_remove(i);
            let j = rng.gen_range(0..lin.len());
            let b = lin.swap_remove(j);
            let parent = nodes.len();
            nodes.push(TreeNode {
                time: t,
                children: Some((a, b)),
            });
            lin.push(parent);
            live -= 1;
        } else if merge_time.is_finite() {
            // A deme merge happens first.
            t = merge_time;
            let (_, child, parent) = merges[next_merge];
            next_merge += 1;
            if let Some(moving) = demes.remove(&child) {
                demes.entry(parent).or_default().extend(moving);
            }
        } else {
            // No merges left and no coalescence possible: force everything into
            // the root deme. (Only reachable on a malformed tree.)
            let root = tree.idx("ROOT");
            let all: Vec<usize> = demes.values().flatten().copied().collect();
            demes.clear();
            demes.insert(root, all);
        }
    }

    let root_id = nodes.len() - 1;

    // --- drop mutations on the branches ---------------------------------
    let alphabet = [b'A', b'C', b'G', b'T'];
    let mut root_seq = vec![0u8; window_bp];
    for b in root_seq.iter_mut() {
        *b = alphabet[rng.gen_range(0..4)];
    }

    let mut seqs: HashMap<usize, Vec<u8>> = HashMap::new();
    seqs.insert(root_id, root_seq);

    // Preorder from the root so a node's parent sequence always exists.
    let mut stack = vec![root_id];
    while let Some(node_id) = stack.pop() {
        let Some((a, b)) = nodes[node_id].children else {
            continue;
        };
        let parent_time = nodes[node_id].time;
        let parent_seq = seqs.get(&node_id).cloned().unwrap();
        for child in [a, b] {
            let branch = (parent_time - nodes[child].time).max(0.0);
            let lambda = mu * window_bp as f64 * branch;
            let n_mut = if lambda > 0.0 {
                Poisson::new(lambda)
                    .map(|p| p.sample(rng) as usize)
                    .unwrap_or(0)
            } else {
                0
            };
            let mut child_seq = parent_seq.clone();
            for _ in 0..n_mut {
                let pos = rng.gen_range(0..window_bp);
                let cur = child_seq[pos];
                loop {
                    let nb = alphabet[rng.gen_range(0..4)];
                    if nb != cur {
                        child_seq[pos] = nb;
                        break;
                    }
                }
            }
            seqs.insert(child, child_seq);
            stack.push(child);
        }
        // Interior sequences are only needed to derive children.
        if node_id != root_id {
            seqs.remove(&node_id);
        }
    }

    let leaves = (0..n).map(|i| seqs.remove(&i).unwrap()).collect();
    let leaf_ids: Vec<usize> = (0..n).collect();
    (leaves, nodes, leaf_ids)
}

/// True TMRCA between two leaves, read off the simulated tree.
fn true_tmrca(nodes: &[TreeNode], a: usize, b: usize) -> f64 {
    // Walk up from each leaf collecting ancestors with their times.
    let parent_of = {
        let mut p = vec![usize::MAX; nodes.len()];
        for (i, n) in nodes.iter().enumerate() {
            if let Some((x, y)) = n.children {
                p[x] = i;
                p[y] = i;
            }
        }
        p
    };
    let mut anc_a = HashSet::new();
    let mut cur = a;
    loop {
        anc_a.insert(cur);
        if parent_of[cur] == usize::MAX {
            break;
        }
        cur = parent_of[cur];
    }
    let mut cur = b;
    loop {
        if anc_a.contains(&cur) {
            return nodes[cur].time;
        }
        if parent_of[cur] == usize::MAX {
            return nodes[cur].time;
        }
        cur = parent_of[cur];
    }
}

// ---------------------------------------------------------------------------
// Cohort construction
// ---------------------------------------------------------------------------

/// Which deme a haplotype's lineage starts in for a given window. Introgressed
/// windows start in the *donor* deme, which is exactly what admixture means
/// backwards in time.
fn source_deme(tree: &PopulationTree, src: ArchaicSource) -> usize {
    match src {
        ArchaicSource::Neanderthal => tree.idx("NEA"),
        ArchaicSource::Denisovan => tree.idx("DEN"),
        ArchaicSource::GhostDeepAfrican => tree.idx("GHOSTA"),
        ArchaicSource::SuperArchaic => tree.idx("SUPER"),
    }
}

/// Build the simulated cohort: haplotypes, per-window sequences, planted truth.
pub fn simulate_cohort(config: DemographyConfig) -> Cohort {
    let tree = build_population_tree();
    let mut rng = StdRng::seed_from_u64(config.seed);

    // Panel composition. Counts are deliberately uneven, as real reference
    // panels are, so the detector cannot lean on a balanced design.
    let panel: Vec<(&'static str, usize, Role)> = vec![
        ("YRI", 8, Role::Modern),
        ("MSL", 6, Role::Modern),
        ("LWK", 6, Role::Modern),
        ("SAN", 4, Role::Modern),
        ("CEU", 8, Role::Modern),
        ("SAS", 6, Role::Modern),
        ("CHB", 8, Role::Modern),
        ("PNG", 8, Role::Modern),
        ("AMR", 4, Role::Modern),
        ("NEA", 3, Role::ArchaicReference),
        ("DEN", 1, Role::ArchaicReference),
    ];

    let mut haplotypes = Vec::new();
    let mut home_deme = Vec::new();
    for (pop, count, role) in &panel {
        let idx = tree.idx(pop);
        for i in 0..*count {
            haplotypes.push(Haplotype {
                id: format!("{pop}_{:02}", i + 1),
                population: pop.to_string(),
                group: tree.nodes[idx].group.to_string(),
                role: *role,
            });
            home_deme.push(idx);
        }
    }
    let n_hap = haplotypes.len();

    // --- plant introgression --------------------------------------------
    //
    // At any given locus an archaic haplotype is rare, so each planted segment
    // goes into exactly one carrier. That keeps truth unambiguous and matches
    // the low per-locus frequency of real archaic tracts.
    let mut truth: HashMap<(usize, usize), ArchaicSource> = HashMap::new();
    let mut tracts: Vec<Tract> = Vec::new();
    for w in 0..config.n_windows {
        for pulse in &config.pulses {
            let eligible: Vec<usize> = (0..n_hap)
                .filter(|&h| {
                    if haplotypes[h].role != Role::Modern {
                        return false;
                    }
                    match pulse.recipient_group.as_str() {
                        "ALL" => true,
                        "NONAFR" => haplotypes[h].group != "AFR",
                        g => haplotypes[h].group == g,
                    }
                })
                .collect();
            if eligible.is_empty() {
                continue;
            }
            // Expected number of tracts *starting* at this window. Longer tracts
            // mean fewer starts for the same total burden, so the genome-wide
            // archaic fraction stays put as tract length changes.
            let mean_len = pulse.tract_mean_windows.max(1.0);
            let expected = pulse.fraction * eligible.len() as f64 / mean_len;
            let n_starts = {
                let floor = expected.floor();
                let frac = expected - floor;
                floor as usize + usize::from(rng.gen::<f64>() < frac)
            };
            for _ in 0..n_starts {
                let h = eligible[rng.gen_range(0..eligible.len())];
                // Recombination breaks a tract at a constant rate per window, so
                // its surviving length is geometric with the requested mean.
                let len = if mean_len <= 1.0 {
                    1
                } else {
                    let p = 1.0 / mean_len;
                    let u: f64 = rng.gen::<f64>().max(1e-12);
                    ((u.ln() / (1.0 - p).ln()).floor() as usize + 1).clamp(1, config.n_windows)
                };
                let end = (w + len).min(config.n_windows);
                if end <= w {
                    continue;
                }
                let mut placed = 0usize;
                for ww in w..end {
                    if truth.contains_key(&(ww, h)) {
                        break; // do not overwrite an earlier tract
                    }
                    truth.insert((ww, h), pulse.source);
                    placed += 1;
                }
                if placed > 0 {
                    tracts.push(Tract {
                        haplotype: h,
                        start_window: w,
                        length: placed,
                        source: pulse.source,
                    });
                }
            }
        }
    }

    // --- simulate each window -------------------------------------------
    let mut sequences = Vec::with_capacity(config.n_windows);
    let mut calibration = Vec::new();
    for w in 0..config.n_windows {
        let mut demes = home_deme.clone();
        for h in 0..n_hap {
            if let Some(src) = truth.get(&(w, h)) {
                demes[h] = source_deme(&tree, *src);
            }
        }
        let (leaves, nodes, _) =
            simulate_window(&tree, &demes, config.window_bp, config.mu, &mut rng);

        // Sample a few pairs per window for the clock calibration report.
        if w % 12 == 0 {
            for _ in 0..6 {
                let a = rng.gen_range(0..n_hap);
                let b = rng.gen_range(0..n_hap);
                if a == b {
                    continue;
                }
                let p = raw_divergence(&leaves[a], &leaves[b]);
                calibration.push(CalibrationPoint {
                    true_tmrca_gens: true_tmrca(&nodes, a, b),
                    observed_divergence: p,
                    estimated_tmrca_gens: divergence_to_tmrca_gens(p, config.mu),
                });
            }
        }
        sequences.push(leaves);
    }

    Cohort {
        config,
        tree,
        haplotypes,
        sequences,
        truth,
        tracts,
        calibration,
    }
}

/// Fraction of positions at which two equal-length sequences differ.
pub fn raw_divergence(a: &[u8], b: &[u8]) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 {
        return 0.0;
    }
    let mut d = 0usize;
    for i in 0..n {
        if a[i] != b[i] {
            d += 1;
        }
    }
    d as f64 / n as f64
}

// ---------------------------------------------------------------------------
// Detector parameters — the Darwin genome
// ---------------------------------------------------------------------------

/// The detector's tunable hyperparameters. Darwin mode mutates exactly these.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DetectorParams {
    /// k-mer size used to embed a window (rvDNA caps this at 15).
    pub k: usize,
    /// Dimensionality of the profile vector fed to HNSW.
    pub dims: usize,
    /// HNSW graph connectivity.
    pub hnsw_m: usize,
    /// HNSW search breadth.
    pub ef_search: usize,
    /// Neighbours retrieved per query before exact divergence is computed.
    pub top_k: usize,
    /// Minimum coalescent depth (ka) to the modern panel to call a segment archaic.
    pub tau_archaic_ka: f64,
    /// A segment is attributed to an observed archaic reference if its depth to
    /// that reference is below this (ka).
    pub match_margin_ka: f64,
    /// Two ghost segments join the same lineage if their mutual depth is below
    /// this (ka).
    pub ghost_link_ka: f64,
}

impl Default for DetectorParams {
    fn default() -> Self {
        Self {
            k: 8,
            dims: 256,
            hnsw_m: 16,
            ef_search: 96,
            top_k: 20,
            tau_archaic_ka: 620.0,
            match_margin_ka: 520.0,
            ghost_link_ka: 900.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Detection results
// ---------------------------------------------------------------------------

/// One archaic call made by the detector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArchaicCall {
    pub window: usize,
    pub haplotype: usize,
    pub haplotype_id: String,
    pub population: String,
    pub group: String,
    /// Coalescent depth to the nearest modern panel haplotype, in ka.
    pub depth_to_modern_ka: f64,
    /// Depth to the closest Neanderthal reference, in ka.
    pub depth_to_neanderthal_ka: f64,
    /// Depth to the closest Denisovan reference, in ka.
    pub depth_to_denisovan_ka: f64,
    /// What the detector concluded.
    pub attribution: String,
    /// Ghost lineage id, once clustering has run.
    pub ghost_cluster: Option<usize>,
    /// Which flywheel round produced this call.
    pub round: usize,
    /// True if the segment cleared the coalescent-depth test on its own
    /// evidence. False means it was rescued by resembling an already-confirmed
    /// ghost reference. Only primary calls are ever promoted to references —
    /// otherwise the flywheel would bootstrap on its own mistakes.
    pub primary: bool,
    /// Depth to the nearest confirmed ghost reference in this window, if any.
    pub depth_to_ghost_ref_ka: Option<f64>,
}

/// Precision / recall / F1 against planted truth.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Scores {
    pub true_positives: usize,
    pub false_positives: usize,
    pub false_negatives: usize,
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
}

impl Scores {
    pub fn compute(tp: usize, fp: usize, fneg: usize) -> Self {
        let precision = if tp + fp == 0 {
            0.0
        } else {
            tp as f64 / (tp + fp) as f64
        };
        let recall = if tp + fneg == 0 {
            0.0
        } else {
            tp as f64 / (tp + fneg) as f64
        };
        let f1 = if precision + recall == 0.0 {
            0.0
        } else {
            2.0 * precision * recall / (precision + recall)
        };
        Self {
            true_positives: tp,
            false_positives: fp,
            false_negatives: fneg,
            precision,
            recall,
            f1,
        }
    }
}

/// Everything one detector pass produced.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectionRun {
    pub params: DetectorParams,
    pub calls: Vec<ArchaicCall>,
    pub scores: Scores,
    /// Exact pairwise divergence computations actually performed.
    pub exact_comparisons: usize,
    /// What an all-pairs scan would have cost.
    pub bruteforce_comparisons: usize,
    pub elapsed_ms: u128,
}

// ---------------------------------------------------------------------------
// The TRACE-rv engine
// ---------------------------------------------------------------------------

/// A per-(k, dims, m) HNSW index over the cohort's window vectors.
struct IndexBundle {
    db: VectorDB,
    /// `vectors[window][haplotype]`, kept so the flywheel can add consensus refs.
    vectors: Vec<Vec<Vec<f32>>>,
    /// Maps a RuVector entry id back to (window, haplotype).
    lookup: HashMap<String, (usize, usize)>,
}

/// The detector. Owns the cohort, the HNSW indices, and the flywheel state.
pub struct TraceEngine<'a> {
    cohort: &'a Cohort,
    /// Cached HNSW indices keyed by the parameters that change the vectors.
    indices: HashMap<(usize, usize, usize), IndexBundle>,
    /// Ghost consensus references written back by the flywheel.
    ghost_refs: Vec<GhostReference>,
    storage_root: String,
    next_index_id: usize,
    /// Running tally of exact divergence computations.
    pub exact_comparisons: usize,
}

/// A ghost lineage made observable by the flywheel: a consensus of the segments
/// already confidently attributed to it.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GhostReference {
    pub cluster: usize,
    pub window: usize,
    pub haplotype: usize,
    pub round: usize,
}

impl<'a> TraceEngine<'a> {
    pub fn new(cohort: &'a Cohort, storage_root: impl Into<String>) -> Self {
        Self {
            cohort,
            indices: HashMap::new(),
            ghost_refs: Vec::new(),
            storage_root: storage_root.into(),
            next_index_id: 0,
            exact_comparisons: 0,
        }
    }

    /// Build (or reuse) the HNSW index for a given embedding configuration.
    fn index_for(&mut self, p: &DetectorParams) -> Result<&IndexBundle> {
        let key = (p.k, p.dims, p.hnsw_m);
        if self.indices.contains_key(&key) {
            return Ok(self.indices.get(&key).unwrap());
        }

        // Bound the cache: Darwin proposes many embedding configurations and
        // each index holds every window vector, so keeping them all would grow
        // without limit.
        const MAX_CACHED_INDICES: usize = 3;
        if self.indices.len() >= MAX_CACHED_INDICES {
            self.indices.clear();
            std::fs::remove_dir_all(&self.storage_root).ok();
        }

        let id = self.next_index_id;
        self.next_index_id += 1;
        let path = format!(
            "{}/hnsw_k{}_d{}_m{}_{}",
            self.storage_root, p.k, p.dims, p.hnsw_m, id
        );
        std::fs::create_dir_all(&self.storage_root).ok();

        let db = VectorDB::new(DbOptions {
            dimensions: p.dims,
            distance_metric: DistanceMetric::Cosine,
            storage_path: path,
            hnsw_config: Some(HnswConfig {
                m: p.hnsw_m,
                ef_construction: (p.hnsw_m * 12).max(120),
                ef_search: p.ef_search,
                max_elements: 1_000_000,
            }),
            quantization: None,
        })
        .map_err(DnaError::VectorDbError)?;

        let n_win = self.cohort.sequences.len();
        let n_hap = self.cohort.haplotypes.len();
        let mut vectors = Vec::with_capacity(n_win);
        let mut lookup = HashMap::new();
        let mut entries = Vec::with_capacity(n_win * n_hap);

        for w in 0..n_win {
            let mut row = Vec::with_capacity(n_hap);
            for h in 0..n_hap {
                let text = std::str::from_utf8(&self.cohort.sequences[w][h])
                    .map_err(|e| DnaError::InvalidSequence(e.to_string()))?;
                let seq = DnaSequence::from_str(text)?;
                let v = seq.to_kmer_vector(p.k, p.dims)?;
                let entry_id = format!("w{w}:h{h}");
                lookup.insert(entry_id.clone(), (w, h));
                let mut meta = HashMap::new();
                meta.insert("window".to_string(), serde_json::json!(w));
                meta.insert("haplotype".to_string(), serde_json::json!(h));
                meta.insert(
                    "role".to_string(),
                    serde_json::json!(match self.cohort.haplotypes[h].role {
                        Role::Modern => "modern",
                        Role::ArchaicReference => "archaic_ref",
                    }),
                );
                entries.push(VectorEntry {
                    id: Some(entry_id),
                    vector: v.clone(),
                    metadata: Some(meta),
                });
                row.push(v);
            }
            vectors.push(row);
        }

        db.insert_batch(entries).map_err(DnaError::VectorDbError)?;

        self.indices.insert(
            key,
            IndexBundle {
                db,
                vectors,
                lookup,
            },
        );
        Ok(self.indices.get(&key).unwrap())
    }

    /// Coalescent depth, in ka, from segment `(w, h)` to its nearest relative
    /// anywhere in the *modern* panel.
    ///
    /// This is the single measurement the whole method rests on, factored out so
    /// that [`TraceEngine::detect`] and [`TraceEngine::depth_grid`] cannot drift
    /// apart. HNSW proposes the handful of genealogically plausible relatives;
    /// exact Jukes–Cantor-corrected divergence is computed only against those. If
    /// HNSW returns nothing usable we fall back to the full panel rather than
    /// make a call on missing evidence.
    ///
    /// Returns `(depth_ka, exact_comparisons_performed)`. The caller owns the
    /// tally so that `detect`'s reported `exact_comparisons` stays exact.
    fn modern_depth(
        &self,
        p: &DetectorParams,
        key: (usize, usize, usize),
        w: usize,
        h: usize,
        modern: &[usize],
    ) -> Result<(f64, usize)> {
        let mu = self.cohort.config.mu;
        let mut exact = 0usize;

        // --- HNSW retrieval: who is even plausibly related? ---------
        let bundle = self.indices.get(&key).unwrap();
        let query = bundle.vectors[w][h].clone();
        let results = bundle
            .db
            .search(SearchQuery {
                vector: query,
                k: p.top_k,
                filter: None,
                ef_search: Some(p.ef_search),
            })
            .map_err(DnaError::VectorDbError)?;

        let mut candidates: Vec<usize> = Vec::new();
        for r in results {
            if let Some(&(rw, rh)) = bundle.lookup.get(&r.id) {
                // Only same-window haplotypes share a genealogy.
                if rw == w && rh != h && self.cohort.haplotypes[rh].role == Role::Modern {
                    candidates.push(rh);
                }
            }
        }
        candidates.sort_unstable();
        candidates.dedup();

        // --- exact divergence to the retrieved candidates -----------
        let seq_h = &self.cohort.sequences[w][h];
        let mut best_modern = f64::INFINITY;
        for &c in &candidates {
            let d = raw_divergence(seq_h, &self.cohort.sequences[w][c]);
            exact += 1;
            let t = gens_to_ka(divergence_to_tmrca_gens(d, mu));
            if t < best_modern {
                best_modern = t;
            }
        }
        if candidates.is_empty() {
            // HNSW returned nothing usable; fall back so we never make a
            // call on missing evidence.
            for &c in modern {
                if c == h {
                    continue;
                }
                let d = raw_divergence(seq_h, &self.cohort.sequences[w][c]);
                exact += 1;
                let t = gens_to_ka(divergence_to_tmrca_gens(d, mu));
                if t < best_modern {
                    best_modern = t;
                }
            }
        }

        Ok((best_modern, exact))
    }

    /// Coalescent depth to the nearest modern relative for **every** segment, as
    /// a dense grid.
    ///
    /// `grid[i][j]` is the depth in ka of window `windows[i]` in the `j`-th
    /// *modern* haplotype, where modern haplotypes are taken in ascending
    /// cohort index order (the same order [`TraceEngine::modern_haplotypes`]
    /// returns). Archaic reference genomes are neither rows nor columns: they are
    /// the thing a call is later attributed *against*, not part of the panel a
    /// segment is measured against.
    ///
    /// Thresholding this grid at `p.tau_archaic_ka` reproduces exactly the calls
    /// [`TraceEngine::detect`] makes on a fresh engine — the grid is the raw
    /// evidence, `detect` is one decision rule over it, and
    /// [`crate::linkage`] is another that also looks sideways at neighbouring
    /// windows.
    pub fn depth_grid(&mut self, p: &DetectorParams, windows: &[usize]) -> Result<Vec<Vec<f64>>> {
        self.index_for(p)?;
        let key = (p.k, p.dims, p.hnsw_m);
        let modern = self.modern_haplotypes();

        let mut grid = Vec::with_capacity(windows.len());
        let mut exact = 0usize;
        for &w in windows {
            let mut row = Vec::with_capacity(modern.len());
            for &h in &modern {
                let (depth, used) = self.modern_depth(p, key, w, h, &modern)?;
                exact += used;
                row.push(depth);
            }
            grid.push(row);
        }
        self.exact_comparisons += exact;
        Ok(grid)
    }

    /// Indices of the modern haplotypes, ascending. These are the columns of
    /// [`TraceEngine::depth_grid`].
    pub fn modern_haplotypes(&self) -> Vec<usize> {
        (0..self.cohort.haplotypes.len())
            .filter(|&h| self.cohort.haplotypes[h].role == Role::Modern)
            .collect()
    }

    /// The cohort this engine is reading.
    pub fn cohort(&self) -> &Cohort {
        self.cohort
    }

    /// Run one detection pass over the given windows.
    pub fn detect(
        &mut self,
        p: &DetectorParams,
        windows: &[usize],
        round: usize,
    ) -> Result<DetectionRun> {
        let start = std::time::Instant::now();

        // Build the index first so we can drop the borrow before mutating self.
        self.index_for(p)?;
        let key = (p.k, p.dims, p.hnsw_m);

        let mu = self.cohort.config.mu;
        let n_hap = self.cohort.haplotypes.len();
        let modern = self.modern_haplotypes();
        let nea_refs: Vec<usize> = (0..n_hap)
            .filter(|&h| self.cohort.haplotypes[h].population == "NEA")
            .collect();
        let den_refs: Vec<usize> = (0..n_hap)
            .filter(|&h| self.cohort.haplotypes[h].population == "DEN")
            .collect();

        let ghost_refs = self.ghost_refs.clone();
        let mut calls: Vec<ArchaicCall> = Vec::new();
        let mut exact = 0usize;
        let mut bruteforce = 0usize;

        for &w in windows {
            for &h in &modern {
                bruteforce += modern.len() - 1 + nea_refs.len() + den_refs.len();

                // --- coalescent depth to the nearest living relative ---------
                let (best_modern, used) = self.modern_depth(p, key, w, h, &modern)?;
                exact += used;
                let seq_h = &self.cohort.sequences[w][h];

                // The flywheel's payoff: once a ghost lineage has confirmed
                // segments, a borderline segment can be rescued by matching one
                // of them directly, without having to clear the depth threshold
                // on its own. Only same-window references are comparable.
                //
                // The rescue rule has to be strict or the loop eats itself. Two
                // copies of the same introgressed haplotype coalesce with each
                // other inside the *donor* population — recently, in coalescent
                // terms — while either one coalesces with a living person only
                // back past the donor's split. So the test is not "close to a
                // ghost" but "closer to a known ghost than to anyone alive".
                let mut best_ghost = f64::INFINITY;
                for g in ghost_refs.iter().filter(|g| g.window == w) {
                    if g.haplotype == h {
                        continue;
                    }
                    let d = raw_divergence(seq_h, &self.cohort.sequences[w][g.haplotype]);
                    exact += 1;
                    let t = gens_to_ka(divergence_to_tmrca_gens(d, mu));
                    if t < best_ghost {
                        best_ghost = t;
                    }
                }
                let rescued = best_ghost < best_modern && best_ghost < p.ghost_link_ka;

                let primary = best_modern >= p.tau_archaic_ka;
                if !primary && !rescued {
                    continue;
                }

                // --- attribute against the observed archaic genomes ---------
                let mut depth_nea = f64::INFINITY;
                for &r in &nea_refs {
                    let d = raw_divergence(seq_h, &self.cohort.sequences[w][r]);
                    exact += 1;
                    depth_nea = depth_nea.min(gens_to_ka(divergence_to_tmrca_gens(d, mu)));
                }
                let mut depth_den = f64::INFINITY;
                for &r in &den_refs {
                    let d = raw_divergence(seq_h, &self.cohort.sequences[w][r]);
                    exact += 1;
                    depth_den = depth_den.min(gens_to_ka(divergence_to_tmrca_gens(d, mu)));
                }

                let attribution = if depth_nea <= p.match_margin_ka && depth_nea <= depth_den {
                    "Neanderthal"
                } else if depth_den <= p.match_margin_ka {
                    "Denisovan"
                } else {
                    "GHOST"
                };

                calls.push(ArchaicCall {
                    window: w,
                    haplotype: h,
                    haplotype_id: self.cohort.haplotypes[h].id.clone(),
                    population: self.cohort.haplotypes[h].population.clone(),
                    group: self.cohort.haplotypes[h].group.clone(),
                    depth_to_modern_ka: best_modern,
                    depth_to_neanderthal_ka: depth_nea,
                    depth_to_denisovan_ka: depth_den,
                    attribution: attribution.to_string(),
                    ghost_cluster: None,
                    round,
                    primary,
                    depth_to_ghost_ref_ka: if best_ghost.is_finite() {
                        Some(best_ghost)
                    } else {
                        None
                    },
                });
            }
        }

        self.exact_comparisons += exact;

        let scores = self.score(&calls, windows);
        Ok(DetectionRun {
            params: *p,
            calls,
            scores,
            exact_comparisons: exact,
            bruteforce_comparisons: bruteforce,
            elapsed_ms: start.elapsed().as_millis(),
        })
    }

    /// Score calls against planted truth, restricted to the given windows.
    pub fn score(&self, calls: &[ArchaicCall], windows: &[usize]) -> Scores {
        let win: HashSet<usize> = windows.iter().copied().collect();
        let called: HashSet<(usize, usize)> =
            calls.iter().map(|c| (c.window, c.haplotype)).collect();
        let truth: HashSet<(usize, usize)> = self
            .cohort
            .truth
            .keys()
            .copied()
            .filter(|(w, h)| win.contains(w) && self.cohort.haplotypes[*h].role == Role::Modern)
            .collect();

        let tp = called.intersection(&truth).count();
        let fp = called.len() - tp;
        let fneg = truth.len() - tp;
        Scores::compute(tp, fp, fneg)
    }

    /// Group GHOST calls into lineages by mutual coalescent depth.
    ///
    /// Two ghost segments in the *same* window can be compared directly. Across
    /// windows we compare their *depth profiles* — a segment from a lineage that
    /// split 1.8 Ma sits at a systematically different depth from the modern
    /// panel than one from a lineage that split 800 ka — which is what lets a
    /// single-carrier-per-window design still resolve two distinct sources.
    pub fn cluster_ghosts(&self, calls: &mut [ArchaicCall]) -> Vec<GhostCluster> {
        let ghost_idx: Vec<usize> = (0..calls.len())
            .filter(|&i| calls[i].attribution == "GHOST")
            .collect();
        if ghost_idx.is_empty() {
            return Vec::new();
        }

        // Cluster on log depth-to-modern-panel. Depth is a direct read of the
        // donor lineage's split time, and it is log-normal-ish rather than
        // normal, because coalescent waiting times are exponential.
        let pts: Vec<(usize, f64)> = ghost_idx
            .iter()
            .map(|&i| (i, calls[i].depth_to_modern_ka.max(1.0).log10()))
            .collect();

        // Two-means, initialised at the extremes so the split is deterministic.
        let lo = pts.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let hi = pts.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
        let mut c = [lo, hi];
        let mut labels = vec![0usize; pts.len()];
        for _ in 0..40 {
            let mut moved = false;
            for (j, p) in pts.iter().enumerate() {
                let l = usize::from((p.1 - c[1]).abs() < (p.1 - c[0]).abs());
                if labels[j] != l {
                    labels[j] = l;
                    moved = true;
                }
            }
            for k in 0..2 {
                let members: Vec<f64> = pts
                    .iter()
                    .zip(&labels)
                    .filter(|(_, l)| **l == k)
                    .map(|(p, _)| p.1)
                    .collect();
                if !members.is_empty() {
                    c[k] = members.iter().sum::<f64>() / members.len() as f64;
                }
            }
            if !moved {
                break;
            }
        }

        // Accept two lineages only if the split is a real separation rather than
        // an arbitrary cut through one distribution: the centroids must be
        // further apart than the pooled within-cluster spread, and each side
        // must hold enough segments to mean anything.
        let mut ss = 0.0;
        for (j, p) in pts.iter().enumerate() {
            ss += (p.1 - c[labels[j]]).powi(2);
        }
        let pooled_sd = (ss / pts.len().max(1) as f64).sqrt();
        let separation = (c[1] - c[0]).abs();
        let n0 = labels.iter().filter(|l| **l == 0).count();
        let n1 = labels.len() - n0;
        let two_clusters = separation > 1.15 * pooled_sd && n0 >= 3 && n1 >= 3;

        let mut assignment: HashMap<usize, usize> = HashMap::new();
        for (j, p) in pts.iter().enumerate() {
            assignment.insert(p.0, if two_clusters { labels[j] } else { 0 });
        }

        for (&i, &c) in assignment.iter() {
            calls[i].ghost_cluster = Some(c);
        }

        let n_clusters = if two_clusters { 2 } else { 1 };
        (0..n_clusters)
            .map(|c| {
                let members: Vec<usize> = ghost_idx
                    .iter()
                    .copied()
                    .filter(|i| assignment.get(i) == Some(&c))
                    .collect();
                let mut depths: Vec<f64> = members
                    .iter()
                    .map(|&i| calls[i].depth_to_modern_ka)
                    .collect();
                depths.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let mean = depths.iter().sum::<f64>() / depths.len().max(1) as f64;
                let var = depths.iter().map(|d| (d - mean).powi(2)).sum::<f64>()
                    / depths.len().max(1) as f64;
                let q = |f: f64| -> f64 {
                    if depths.is_empty() {
                        return 0.0;
                    }
                    depths[(((depths.len() - 1) as f64) * f).round() as usize]
                };
                let mut by_group: BTreeMap<String, usize> = BTreeMap::new();
                for &i in &members {
                    *by_group.entry(calls[i].group.clone()).or_insert(0) += 1;
                }
                GhostCluster {
                    id: c,
                    n_segments: members.len(),
                    mean_split_ka: mean,
                    sd_split_ka: var.sqrt(),
                    // A segment's coalescent depth is always older than the
                    // split that produced it, so the lower tail — not the mean —
                    // is the usable divergence-time estimator.
                    p10_split_ka: q(0.10),
                    median_split_ka: q(0.50),
                    min_split_ka: depths.first().copied().unwrap_or(0.0),
                    max_split_ka: depths.last().copied().unwrap_or(0.0),
                    carriers_by_group: by_group,
                }
            })
            .collect()
    }

    /// Register confirmed ghost segments as references for the next round.
    ///
    /// Only *primary* calls are eligible. A segment that was itself rescued by
    /// resembling a reference can never become a reference, which is what stops
    /// the loop from compounding its own errors.
    pub fn add_ghost_references(&mut self, calls: &[ArchaicCall], round: usize) -> usize {
        let mut added = 0;
        for c in calls {
            if c.attribution != "GHOST" || !c.primary {
                continue;
            }
            let Some(cluster) = c.ghost_cluster else {
                continue;
            };
            let already = self
                .ghost_refs
                .iter()
                .any(|g| g.window == c.window && g.haplotype == c.haplotype);
            if already {
                continue;
            }
            self.ghost_refs.push(GhostReference {
                cluster,
                window: c.window,
                haplotype: c.haplotype,
                round,
            });
            added += 1;
        }
        added
    }

    pub fn ghost_reference_count(&self) -> usize {
        self.ghost_refs.len()
    }
}

/// A recovered ghost lineage.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GhostCluster {
    pub id: usize,
    pub n_segments: usize,
    /// Mean inferred coalescent depth of this lineage against modern humans (ka).
    /// Biased *older* than the true split — see [`GhostCluster::p10_split_ka`].
    pub mean_split_ka: f64,
    pub sd_split_ka: f64,
    /// 10th percentile of the depth distribution: the practical divergence-time
    /// estimator, since coalescence cannot happen before the split but can
    /// happen arbitrarily long after it.
    pub p10_split_ka: f64,
    pub median_split_ka: f64,
    pub min_split_ka: f64,
    pub max_split_ka: f64,
    /// How many segments each reporting group contributed.
    pub carriers_by_group: BTreeMap<String, usize>,
}

// ---------------------------------------------------------------------------
// Darwin mode
// ---------------------------------------------------------------------------

pub mod darwin {
    use super::*;

    /// One accepted or rejected mutation in the evolutionary record.
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct Mutation {
        pub generation: usize,
        pub gene: String,
        pub from: String,
        pub to: String,
        pub train_f1: f64,
        pub accepted: bool,
    }

    /// The result of a Darwin run.
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct Evolution {
        pub baseline: DetectorParams,
        pub baseline_train_f1: f64,
        pub baseline_test: Scores,
        pub evolved: DetectorParams,
        pub evolved_train_f1: f64,
        pub evolved_test: Scores,
        pub lineage: Vec<Mutation>,
        pub generations: usize,
        pub proposals: usize,
        pub accepted: usize,
    }

    fn mutate(p: &DetectorParams, rng: &mut StdRng) -> (DetectorParams, String, String, String) {
        let mut q = *p;
        let gene = rng.gen_range(0..8);
        let (name, from, to) = match gene {
            0 => {
                let from = q.k;
                let opts = [6usize, 7, 8, 9, 10, 11, 12];
                q.k = opts[rng.gen_range(0..opts.len())];
                ("k", from.to_string(), q.k.to_string())
            }
            1 => {
                let from = q.dims;
                let opts = [64usize, 128, 256, 384, 512];
                q.dims = opts[rng.gen_range(0..opts.len())];
                ("dims", from.to_string(), q.dims.to_string())
            }
            2 => {
                let from = q.hnsw_m;
                let opts = [8usize, 12, 16, 24, 32];
                q.hnsw_m = opts[rng.gen_range(0..opts.len())];
                ("hnsw_m", from.to_string(), q.hnsw_m.to_string())
            }
            3 => {
                let from = q.ef_search;
                let opts = [32usize, 64, 96, 128, 192, 256];
                q.ef_search = opts[rng.gen_range(0..opts.len())];
                ("ef_search", from.to_string(), q.ef_search.to_string())
            }
            4 => {
                let from = q.top_k;
                let opts = [8usize, 12, 16, 20, 28, 36, 48];
                q.top_k = opts[rng.gen_range(0..opts.len())];
                ("top_k", from.to_string(), q.top_k.to_string())
            }
            5 => {
                let from = q.tau_archaic_ka;
                q.tau_archaic_ka =
                    (q.tau_archaic_ka + rng.gen_range(-90.0..90.0)).clamp(300.0, 1400.0);
                (
                    "tau_archaic_ka",
                    format!("{from:.0}"),
                    format!("{:.0}", q.tau_archaic_ka),
                )
            }
            6 => {
                let from = q.match_margin_ka;
                q.match_margin_ka =
                    (q.match_margin_ka + rng.gen_range(-80.0..80.0)).clamp(250.0, 1200.0);
                (
                    "match_margin_ka",
                    format!("{from:.0}"),
                    format!("{:.0}", q.match_margin_ka),
                )
            }
            _ => {
                let from = q.ghost_link_ka;
                q.ghost_link_ka =
                    (q.ghost_link_ka + rng.gen_range(-120.0..120.0)).clamp(400.0, 2000.0);
                (
                    "ghost_link_ka",
                    format!("{from:.0}"),
                    format!("{:.0}", q.ghost_link_ka),
                )
            }
        };
        (q, name.to_string(), from, to)
    }

    /// Hill-climb the detector's hyperparameters against a train split, keeping
    /// only mutations that measurably improve F1, then report on a held-out
    /// test split so any overfitting is visible rather than hidden.
    pub fn evolve(
        engine: &mut TraceEngine<'_>,
        baseline: DetectorParams,
        train: &[usize],
        test: &[usize],
        generations: usize,
        proposals_per_gen: usize,
        seed: u64,
    ) -> Result<Evolution> {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut lineage = Vec::new();

        let base_train = engine.detect(&baseline, train, 0)?.scores.f1;
        let base_test = engine.detect(&baseline, test, 0)?.scores;

        let mut best = baseline;
        let mut best_f1 = base_train;
        let mut proposals = 0usize;
        let mut accepted = 0usize;

        // A mutation must clear this margin to count as a real improvement
        // rather than evaluation noise.
        const EPSILON: f64 = 0.004;

        for gen in 1..=generations {
            for _ in 0..proposals_per_gen {
                let (cand, gene, from, to) = mutate(&best, &mut rng);
                if cand == best {
                    continue;
                }
                proposals += 1;
                let f1 = engine.detect(&cand, train, 0)?.scores.f1;
                let take = f1 > best_f1 + EPSILON;
                lineage.push(Mutation {
                    generation: gen,
                    gene,
                    from,
                    to,
                    train_f1: f1,
                    accepted: take,
                });
                if take {
                    best = cand;
                    best_f1 = f1;
                    accepted += 1;
                }
            }
        }

        let evolved_test = engine.detect(&best, test, 0)?.scores;

        Ok(Evolution {
            baseline,
            baseline_train_f1: base_train,
            baseline_test: base_test,
            evolved: best,
            evolved_train_f1: best_f1,
            evolved_test,
            lineage,
            generations,
            proposals,
            accepted,
        })
    }
}

// ---------------------------------------------------------------------------
// Flywheel
// ---------------------------------------------------------------------------

pub mod flywheel {
    use super::*;

    /// One turn of the flywheel.
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct Round {
        pub round: usize,
        pub calls: usize,
        pub ghost_calls: usize,
        pub new_ghost_references: usize,
        pub cumulative_ghost_references: usize,
        pub scores: Scores,
        pub clusters: Vec<GhostCluster>,
        pub elapsed_ms: u128,
    }

    /// The whole flywheel run.
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct FlywheelRun {
        pub rounds: Vec<Round>,
        pub final_calls: Vec<ArchaicCall>,
        pub final_clusters: Vec<GhostCluster>,
        pub converged_after: usize,
    }

    /// Detect, condense confirmed ghosts into references, re-index, repeat.
    /// Stops after `dry_rounds` consecutive rounds that add nothing new.
    pub fn spin(
        engine: &mut TraceEngine<'_>,
        params: &DetectorParams,
        windows: &[usize],
        max_rounds: usize,
        dry_rounds: usize,
    ) -> Result<FlywheelRun> {
        let mut rounds = Vec::new();
        let mut dry = 0usize;
        let mut final_calls = Vec::new();
        let mut final_clusters = Vec::new();
        let mut converged_after = max_rounds;

        for r in 1..=max_rounds {
            let start = std::time::Instant::now();
            let run = engine.detect(params, windows, r)?;
            let mut calls = run.calls;
            let clusters = engine.cluster_ghosts(&mut calls);
            let ghost_calls = calls.iter().filter(|c| c.attribution == "GHOST").count();
            let added = engine.add_ghost_references(&calls, r);

            rounds.push(Round {
                round: r,
                calls: calls.len(),
                ghost_calls,
                new_ghost_references: added,
                cumulative_ghost_references: engine.ghost_reference_count(),
                scores: run.scores,
                clusters: clusters.clone(),
                elapsed_ms: start.elapsed().as_millis(),
            });

            final_calls = calls;
            final_clusters = clusters;

            if added == 0 {
                dry += 1;
                if dry >= dry_rounds {
                    converged_after = r;
                    break;
                }
            } else {
                dry = 0;
            }
        }

        Ok(FlywheelRun {
            rounds,
            final_calls,
            final_clusters,
            converged_after,
        })
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_config() -> DemographyConfig {
        DemographyConfig {
            window_bp: 800,
            n_windows: 24,
            ..Default::default()
        }
    }

    #[test]
    fn population_tree_is_well_formed() {
        let t = build_population_tree();
        // Every non-root node has a parent that splits deeper in time.
        for n in &t.nodes {
            if let Some(p) = n.parent {
                assert!(
                    t.nodes[p].split_gens >= n.split_gens,
                    "{} splits after its parent {}",
                    n.name,
                    t.nodes[p].name
                );
            }
        }
        assert!(t.nodes[t.idx("SUPER")].split_gens > t.nodes[t.idx("GHOSTA")].split_gens);
    }

    #[test]
    fn jukes_cantor_is_monotone_and_conservative() {
        assert!(jukes_cantor(0.0) == 0.0);
        assert!(jukes_cantor(0.1) > 0.1);
        assert!(jukes_cantor(0.2) > jukes_cantor(0.1));
    }

    #[test]
    fn simulated_cohort_has_planted_truth() {
        let c = simulate_cohort(tiny_config());
        assert_eq!(c.sequences.len(), 24);
        assert_eq!(c.sequences[0].len(), c.haplotypes.len());
        assert!(!c.truth.is_empty(), "no introgression was planted");
        // Every sequence is valid DNA of the right length.
        for w in &c.sequences {
            for s in w {
                assert_eq!(s.len(), 800);
                assert!(s.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')));
            }
        }
    }

    #[test]
    fn archaic_segments_are_deeper_than_ordinary_variation() {
        let c = simulate_cohort(tiny_config());
        let modern: Vec<usize> = (0..c.haplotypes.len())
            .filter(|&h| c.haplotypes[h].role == Role::Modern)
            .collect();

        let mut archaic_depths = Vec::new();
        let mut ordinary_depths = Vec::new();
        for w in 0..c.sequences.len() {
            for &h in &modern {
                let mut best = f64::INFINITY;
                for &o in &modern {
                    if o == h {
                        continue;
                    }
                    let d = raw_divergence(&c.sequences[w][h], &c.sequences[w][o]);
                    best = best.min(gens_to_ka(divergence_to_tmrca_gens(d, c.config.mu)));
                }
                if c.truth.contains_key(&(w, h)) {
                    archaic_depths.push(best);
                } else {
                    ordinary_depths.push(best);
                }
            }
        }
        assert!(!archaic_depths.is_empty());
        let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
        assert!(
            mean(&archaic_depths) > mean(&ordinary_depths) * 1.5,
            "introgressed segments should coalesce far deeper: {:.0} vs {:.0} ka",
            mean(&archaic_depths),
            mean(&ordinary_depths)
        );
    }

    /// Regression guard for the flywheel feedback collapse.
    ///
    /// The first implementation promoted *any* GHOST call to a reference and
    /// rescued anything within `ghost_link_ka` of one. False positives became
    /// references, validated more false positives, and precision fell from 0.84
    /// to 0.11 over three rounds. Two rules stop that: only primary calls are
    /// promoted, and a rescue requires the segment to be closer to a known ghost
    /// than to anyone living. This test fails if either is removed.
    #[test]
    fn flywheel_does_not_bootstrap_on_its_own_errors() {
        let cohort = simulate_cohort(DemographyConfig {
            window_bp: 900,
            n_windows: 30,
            ..Default::default()
        });
        let storage = std::env::temp_dir().join(format!(
            "trace_rv_test_{}_{}",
            std::process::id(),
            "flywheel"
        ));
        std::fs::remove_dir_all(&storage).ok();
        let mut engine = TraceEngine::new(&cohort, storage.to_string_lossy().to_string());

        let params = DetectorParams::default();
        let windows: Vec<usize> = (0..cohort.sequences.len()).collect();
        let run = flywheel::spin(&mut engine, &params, &windows, 4, 2).expect("flywheel");

        std::fs::remove_dir_all(&storage).ok();

        assert!(!run.rounds.is_empty());
        let first = run.rounds[0].scores.precision;
        for r in &run.rounds {
            assert!(
                r.scores.precision >= first - 0.15,
                "precision collapsed on round {}: {:.3} vs {:.3} in round 1 — the \
                 flywheel is promoting its own false positives",
                r.round,
                r.scores.precision,
                first
            );
        }

        // A rescued call must never have been promoted.
        let promoted_non_primary = run
            .final_calls
            .iter()
            .any(|c| !c.primary && c.attribution == "GHOST" && c.depth_to_ghost_ref_ka.is_none());
        assert!(
            !promoted_non_primary,
            "a non-primary ghost call exists with no ghost reference to explain it"
        );
    }

    #[test]
    fn scores_are_sane() {
        let s = Scores::compute(8, 2, 2);
        assert!((s.precision - 0.8).abs() < 1e-9);
        assert!((s.recall - 0.8).abs() < 1e-9);
        assert!((s.f1 - 0.8).abs() < 1e-9);
    }
}
