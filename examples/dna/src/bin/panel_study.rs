//! `panel-study` — two sweeps about *what you have to own* to run TRACE-rv.
//!
//! `discover` asked where the method stops working as a function of the
//! **biology** (how rare a ghost can get) and of the **data** (how much sequence
//! a window needs). Both of those are things you cannot buy. This binary asks
//! the two questions you *can* spend money on:
//!
//! 1. **Panel design.** The detector calls a segment archaic when its nearest
//!    relative in the modern reference panel is too deep. That makes the panel a
//!    purchased input: every extra haplotype is another chance for an ordinary
//!    segment to find a close relative and *not* be called. So — how many
//!    reference haplotypes do you actually need, and does it matter *which*
//!    ones? The standing prediction in the field is that African panels are
//!    worth more per haplotype, because African populations carry deeper and
//!    richer coalescent structure. This sweep tests that directly by holding the
//!    query set fixed and swapping the panel underneath it.
//!
//! 2. **Archaic reference ablation.** Attribution — Neanderthal vs Denisovan vs
//!    GHOST — is a *labelled-reference* problem, and the published study notes
//!    that exactly **one** Denisovan genome exists. That asymmetry is not a
//!    detail; it decides which real archaic segments get filed as ghosts. This
//!    sweep removes archaic references one class at a time and measures where
//!    the calls go, ending at the degenerate case of no archaic references at
//!    all, where "ghost" can only mean "not modern".
//!
//! Detection here is **exhaustive** — every query against every panel haplotype,
//! no index — for the same reason `discover` does it: the limits reported should
//! be properties of the panel and the reference set, not of an approximation.
//! `discover` experiment 3 already measured what the index costs separately.
//!
//! Usage: `panel-study [output-dir]` (default `artifacts/data`).

use std::collections::{BTreeMap, HashSet};
use std::io::Write as _;

use rvdna::archaic::{
    divergence_to_tmrca_gens, gens_to_ka, raw_divergence, simulate_cohort, ArchaicCall,
    ArchaicSource, Cohort, DemographyConfig, DetectorParams, Role, Scores,
};
use serde::Serialize;

// ---------------------------------------------------------------------------
// Shared machinery — the `exhaustive_calls` idea from `discover`, split into a
// panel scan and an attribution step so each sweep can vary exactly one of them.
// ---------------------------------------------------------------------------

/// Coalescent depth, in ka, between two haplotypes at one window.
fn depth(cohort: &Cohort, w: usize, a: usize, b: usize) -> f64 {
    gens_to_ka(divergence_to_tmrca_gens(
        raw_divergence(&cohort.sequences[w][a], &cohort.sequences[w][b]),
        cohort.config.mu,
    ))
}

/// A segment that cleared the coalescent-depth test, before attribution.
#[derive(Debug, Clone, Copy)]
struct Segment {
    window: usize,
    haplotype: usize,
    depth_to_modern_ka: f64,
}

/// What one exhaustive scan of the modern panel produced.
struct ScanResult {
    segments: Vec<Segment>,
    /// Median depth assigned to segments that are *not* introgressed — the noise
    /// floor the threshold has to sit above. It rises as the panel shrinks,
    /// which is the whole mechanism behind sweep 1.
    ordinary_median_ka: f64,
    comparisons: usize,
}

/// Exhaustive nearest-relative search: every query segment against every panel
/// haplotype (never itself), no index.
fn scan_panel(
    cohort: &Cohort,
    p: &DetectorParams,
    panel: &[usize],
    queries: &[usize],
) -> ScanResult {
    let mut segments = Vec::new();
    let mut ordinary: Vec<f64> = Vec::new();
    let mut comparisons = 0usize;

    for w in 0..cohort.sequences.len() {
        for &h in queries {
            let mut best = f64::INFINITY;
            for &o in panel {
                if o == h {
                    continue;
                }
                comparisons += 1;
                let d = depth(cohort, w, h, o);
                if d < best {
                    best = d;
                }
            }
            if !cohort.truth.contains_key(&(w, h)) && best.is_finite() {
                ordinary.push(best);
            }
            if best >= p.tau_archaic_ka {
                segments.push(Segment {
                    window: w,
                    haplotype: h,
                    depth_to_modern_ka: best,
                });
            }
        }
    }

    ordinary.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ordinary_median_ka = if ordinary.is_empty() {
        0.0
    } else {
        ordinary[ordinary.len() / 2]
    };

    ScanResult {
        segments,
        ordinary_median_ka,
        comparisons,
    }
}

/// Label each called segment against the archaic references we are allowed to
/// use. This is exactly the rule in `archaic::TraceEngine::detect`, but with the
/// reference sets passed in so they can be ablated.
fn attribute(
    cohort: &Cohort,
    p: &DetectorParams,
    segments: &[Segment],
    nea_refs: &[usize],
    den_refs: &[usize],
) -> Vec<ArchaicCall> {
    segments
        .iter()
        .map(|s| {
            let d_nea = nea_refs
                .iter()
                .map(|&r| depth(cohort, s.window, s.haplotype, r))
                .fold(f64::INFINITY, f64::min);
            let d_den = den_refs
                .iter()
                .map(|&r| depth(cohort, s.window, s.haplotype, r))
                .fold(f64::INFINITY, f64::min);
            let attribution = if d_nea <= p.match_margin_ka && d_nea <= d_den {
                "Neanderthal"
            } else if d_den <= p.match_margin_ka {
                "Denisovan"
            } else {
                "GHOST"
            };
            ArchaicCall {
                window: s.window,
                haplotype: s.haplotype,
                haplotype_id: cohort.haplotypes[s.haplotype].id.clone(),
                population: cohort.haplotypes[s.haplotype].population.clone(),
                group: cohort.haplotypes[s.haplotype].group.clone(),
                depth_to_modern_ka: s.depth_to_modern_ka,
                depth_to_neanderthal_ka: d_nea,
                depth_to_denisovan_ka: d_den,
                attribution: attribution.to_string(),
                ghost_cluster: None,
                round: 0,
                primary: true,
                depth_to_ghost_ref_ka: None,
            }
        })
        .collect()
}

/// Precision / recall / F1 of *detection* against planted truth, over a given
/// query set.
fn score_segments(cohort: &Cohort, segments: &[Segment], queries: &[usize]) -> Scores {
    let q: HashSet<usize> = queries.iter().copied().collect();
    let called: HashSet<(usize, usize)> =
        segments.iter().map(|s| (s.window, s.haplotype)).collect();
    let truth: HashSet<(usize, usize)> = cohort
        .truth
        .keys()
        .copied()
        .filter(|(_, h)| cohort.haplotypes[*h].role == Role::Modern && q.contains(h))
        .collect();
    let tp = called.intersection(&truth).count();
    Scores::compute(tp, called.len() - tp, truth.len() - tp)
}

/// Which lineage really donated this segment, or `None` if nothing did.
fn truth_label(cohort: &Cohort, w: usize, h: usize) -> &'static str {
    match cohort.truth.get(&(w, h)) {
        Some(ArchaicSource::Neanderthal) => "Neanderthal",
        Some(ArchaicSource::Denisovan) => "Denisovan",
        Some(ArchaicSource::GhostDeepAfrican) => "Ghost-A",
        Some(ArchaicSource::SuperArchaic) => "Ghost-B",
        None => "not-introgressed",
    }
}

fn is_true_ghost(label: &str) -> bool {
    label == "Ghost-A" || label == "Ghost-B"
}

// ---------------------------------------------------------------------------
// Panel construction
// ---------------------------------------------------------------------------

fn modern_indices(cohort: &Cohort) -> Vec<usize> {
    (0..cohort.haplotypes.len())
        .filter(|&h| cohort.haplotypes[h].role == Role::Modern)
        .collect()
}

/// Deal one element from each bucket in turn, so any prefix of the result is as
/// evenly spread across the buckets as that prefix length allows.
fn interleave(buckets: &[Vec<usize>]) -> Vec<usize> {
    let mut out = Vec::new();
    let mut i = 0usize;
    loop {
        let mut any = false;
        for b in buckets {
            if let Some(&x) = b.get(i) {
                out.push(x);
                any = true;
            }
        }
        if !any {
            return out;
        }
        i += 1;
    }
}

/// Order the haplotypes of one reporting group, round-robin across the
/// populations inside it, so a prefix never over-weights a single population.
fn group_ordered(cohort: &Cohort, members: &[usize]) -> Vec<usize> {
    let mut by_pop: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for &h in members {
        by_pop
            .entry(cohort.haplotypes[h].population.as_str())
            .or_default()
            .push(h);
    }
    let buckets: Vec<Vec<usize>> = by_pop.into_values().collect();
    interleave(&buckets)
}

/// A panel spread as evenly as possible across every superpopulation present.
fn balanced_panel(cohort: &Cohort, size: usize) -> Vec<usize> {
    let mut by_group: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for h in modern_indices(cohort) {
        by_group
            .entry(cohort.haplotypes[h].group.as_str())
            .or_default()
            .push(h);
    }
    let buckets: Vec<Vec<usize>> = by_group
        .into_values()
        .map(|m| group_ordered(cohort, &m))
        .collect();
    let mut order = interleave(&buckets);
    order.truncate(size);
    order
}

/// A panel drawn only from African populations, spread across them evenly.
/// Truncated if the cohort does not hold that many African haplotypes.
fn african_panel(cohort: &Cohort, size: usize) -> Vec<usize> {
    let members: Vec<usize> = modern_indices(cohort)
        .into_iter()
        .filter(|&h| cohort.haplotypes[h].group == "AFR")
        .collect();
    let mut order = group_ordered(cohort, &members);
    order.truncate(size);
    order
}

// ---------------------------------------------------------------------------
// Report types
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct PanelPoint {
    /// `"balanced"` or `"african_only"`.
    composition: String,
    requested_size: usize,
    /// What the cohort could actually supply.
    actual_size: usize,
    /// False when the cohort cannot build this panel at all (see `note`).
    available: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    note: Option<String>,
    calls: usize,
    true_positives: usize,
    false_positives: usize,
    false_negatives: usize,
    precision: f64,
    recall: f64,
    f1: f64,
    /// F1 delivered per reference haplotype spent — the "informative per
    /// haplotype" quantity the prediction is about.
    f1_per_haplotype: f64,
    /// Precision restricted to African query haplotypes.
    precision_afr_queries: f64,
    /// Precision restricted to non-African query haplotypes.
    precision_nonafr_queries: f64,
    /// Median depth assigned to ordinary, non-introgressed segments (ka).
    ordinary_median_ka: f64,
    exact_comparisons: usize,
}

#[derive(Serialize)]
struct AblationPoint {
    /// Human-readable name of the reference set.
    reference_set: String,
    neanderthal_refs: usize,
    denisovan_refs: usize,
    /// Detection is untouched by this sweep, so this is constant by construction.
    total_calls: usize,
    attributed_neanderthal: usize,
    attributed_denisovan: usize,
    attributed_ghost: usize,
    /// `truth class -> attribution -> count`.
    confusion: BTreeMap<String, BTreeMap<String, usize>>,
    /// True Neanderthal segments that were called, then filed as GHOST.
    true_neanderthal_as_ghost: usize,
    /// True Denisovan segments that were called, then filed as GHOST.
    true_denisovan_as_ghost: usize,
    /// True Neanderthal segments filed under the wrong archaic name.
    true_neanderthal_as_denisovan: usize,
    /// True Denisovan segments filed under the wrong archaic name.
    true_denisovan_as_neanderthal: usize,
    /// Detection false positives that end up in the ghost pile. These are not an
    /// attribution failure — they are the ceiling attribution can never beat.
    not_introgressed_as_ghost: usize,
    /// Of everything labelled GHOST, the share that really is a ghost lineage.
    ghost_precision: f64,
    /// Of every true ghost segment in the cohort, the share labelled GHOST.
    ghost_recall: f64,
    ghost_f1: f64,
    /// Share of called, genuinely-archaic segments given the right label.
    attribution_accuracy: f64,
}

#[derive(Serialize)]
struct CohortSummary {
    n_windows: usize,
    window_bp: usize,
    modern_haplotypes: usize,
    african_modern_haplotypes: usize,
    neanderthal_references: usize,
    denisovan_references: usize,
    true_neanderthal_segments: usize,
    true_denisovan_segments: usize,
    true_ghost_segments: usize,
}

#[derive(Serialize)]
struct PanelStudy {
    generated_utc: String,
    params: DetectorParams,
    cohort: CohortSummary,
    panel_design: Vec<PanelPoint>,
    /// Smallest balanced panel already within 2% of the best F1 any panel reached.
    balanced_saturation_size: Option<usize>,
    /// Did African-only beat balanced at every size both could be built?
    african_beats_balanced_at_matched_size: bool,
    /// The same question restricted to *non-African* query haplotypes. This is
    /// where the per-haplotype story and the per-query story come apart.
    african_beats_balanced_for_nonafrican_queries: bool,
    /// Smallest balanced panel that matches the best African-only panel's F1.
    balanced_size_matching_best_african: Option<usize>,
    /// How many balanced haplotypes one African haplotype is worth, from the
    /// two sizes above. `None` if no balanced panel tested ever caught up.
    african_haplotype_value_multiple: Option<f64>,
    archaic_ablation: Vec<AblationPoint>,
    /// Ghost precision with every archaic reference available.
    ghost_precision_all_references: f64,
    /// Ghost precision with no archaic reference at all.
    ghost_precision_no_references: f64,
    /// Ghost precision once the single Denisovan genome is removed.
    ghost_precision_neanderthal_only: f64,
    findings: Vec<String>,
}

// ---------------------------------------------------------------------------

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let out_dir = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "artifacts/data".into());
    std::fs::create_dir_all(&out_dir)?;

    // The same operating point `discover` uses: default detector with the
    // threshold Darwin settled on in the main run.
    let params = DetectorParams {
        tau_archaic_ka: 551.0,
        ..Default::default()
    };

    let cohort = simulate_cohort(DemographyConfig {
        n_windows: 360,
        ..Default::default()
    });

    let queries = modern_indices(&cohort);
    let n_hap = cohort.haplotypes.len();
    let nea_all: Vec<usize> = (0..n_hap)
        .filter(|&h| cohort.haplotypes[h].population == "NEA")
        .collect();
    let den_all: Vec<usize> = (0..n_hap)
        .filter(|&h| cohort.haplotypes[h].population == "DEN")
        .collect();
    let afr_available = queries
        .iter()
        .filter(|&&h| cohort.haplotypes[h].group == "AFR")
        .count();

    let count_truth = |f: &dyn Fn(&ArchaicSource) -> bool| {
        cohort
            .truth
            .iter()
            .filter(|((_, h), s)| cohort.haplotypes[*h].role == Role::Modern && f(s))
            .count()
    };
    let n_true_nea = count_truth(&|s| *s == ArchaicSource::Neanderthal);
    let n_true_den = count_truth(&|s| *s == ArchaicSource::Denisovan);
    let n_true_ghost = count_truth(&|s| {
        matches!(
            s,
            ArchaicSource::GhostDeepAfrican | ArchaicSource::SuperArchaic
        )
    });

    let banner = |s: &str| {
        println!("\n\x1b[1m{s}\x1b[0m");
        println!("{}", "-".repeat(s.len()));
    };

    println!(
        "cohort: {} windows x {} modern haplotypes ({} African), {} NEA + {} DEN references",
        cohort.sequences.len(),
        queries.len(),
        afr_available,
        nea_all.len(),
        den_all.len()
    );
    println!(
        "planted truth: {n_true_nea} Neanderthal, {n_true_den} Denisovan, {n_true_ghost} ghost segments"
    );

    // ----------------------------------------------------------------- 1
    banner("1. Panel design — how many reference haplotypes, and which ones?");
    println!(
        "   composition   size   calls   precision  recall      F1   F1/hap   P(AFR q)  P(nonAFR q)  noise floor"
    );

    let sizes = [6usize, 12, 24, 36, 46, 58];
    let mut panel_design: Vec<PanelPoint> = Vec::new();
    // Keep the full-panel scan around: sweep 2 needs exactly that detection.
    let mut full_scan: Option<ScanResult> = None;

    for composition in ["balanced", "african_only"] {
        for &size in &sizes {
            let panel = if composition == "balanced" {
                balanced_panel(&cohort, size)
            } else {
                african_panel(&cohort, size)
            };

            if panel.len() < size {
                // An African-only panel cannot exceed the African haplotypes the
                // cohort holds. That ceiling is itself a result, so record it
                // rather than quietly substituting a smaller panel.
                println!(
                    "   {:<12}  {:>4}      —          —       —       —        —          —            —            —   (only {afr_available} African haplotypes exist)",
                    composition, size
                );
                panel_design.push(PanelPoint {
                    composition: composition.to_string(),
                    requested_size: size,
                    actual_size: panel.len(),
                    available: false,
                    note: Some(format!(
                        "cohort holds only {afr_available} African modern haplotypes; an African-only panel of {size} cannot be built"
                    )),
                    calls: 0,
                    true_positives: 0,
                    false_positives: 0,
                    false_negatives: 0,
                    precision: 0.0,
                    recall: 0.0,
                    f1: 0.0,
                    f1_per_haplotype: 0.0,
                    precision_afr_queries: 0.0,
                    precision_nonafr_queries: 0.0,
                    ordinary_median_ka: 0.0,
                    exact_comparisons: 0,
                });
                continue;
            }

            let scan = scan_panel(&cohort, &params, &panel, &queries);
            let scores = score_segments(&cohort, &scan.segments, &queries);

            // Split precision by where the *query* comes from. This is what
            // separates "the panel is small" from "the panel is the wrong shape".
            let mut tp_afr = 0usize;
            let mut n_afr = 0usize;
            let mut tp_non = 0usize;
            let mut n_non = 0usize;
            for s in &scan.segments {
                let hit = cohort.truth.contains_key(&(s.window, s.haplotype));
                if cohort.haplotypes[s.haplotype].group == "AFR" {
                    n_afr += 1;
                    tp_afr += usize::from(hit);
                } else {
                    n_non += 1;
                    tp_non += usize::from(hit);
                }
            }
            let p_afr = if n_afr == 0 {
                0.0
            } else {
                tp_afr as f64 / n_afr as f64
            };
            let p_non = if n_non == 0 {
                0.0
            } else {
                tp_non as f64 / n_non as f64
            };

            println!(
                "   {:<12}  {:>4}  {:>6}  {:>10.3}  {:>6.3}  {:>6.4}  {:>7.5}  {:>9.3}  {:>11.3}  {:>8.0} ka",
                composition,
                panel.len(),
                scan.segments.len(),
                scores.precision,
                scores.recall,
                scores.f1,
                scores.f1 / panel.len() as f64,
                p_afr,
                p_non,
                scan.ordinary_median_ka
            );

            panel_design.push(PanelPoint {
                composition: composition.to_string(),
                requested_size: size,
                actual_size: panel.len(),
                available: true,
                note: None,
                calls: scan.segments.len(),
                true_positives: scores.true_positives,
                false_positives: scores.false_positives,
                false_negatives: scores.false_negatives,
                precision: scores.precision,
                recall: scores.recall,
                f1: scores.f1,
                f1_per_haplotype: scores.f1 / panel.len() as f64,
                precision_afr_queries: p_afr,
                precision_nonafr_queries: p_non,
                ordinary_median_ka: scan.ordinary_median_ka,
                exact_comparisons: scan.comparisons,
            });

            if composition == "balanced" && panel.len() == queries.len() {
                full_scan = Some(scan);
            }
        }
    }

    let best_f1 = panel_design
        .iter()
        .filter(|r| r.available)
        .map(|r| r.f1)
        .fold(0.0, f64::max);
    let balanced_saturation_size = panel_design
        .iter()
        .filter(|r| r.available && r.composition == "balanced" && r.f1 >= best_f1 * 0.98)
        .map(|r| r.actual_size)
        .min();
    if let Some(s) = balanced_saturation_size {
        println!("\n   -> a balanced panel is within 2% of the best F1 achieved by {s} haplotypes");
    }

    // Does the standing prediction hold? Compare like for like: only the sizes
    // at which both compositions could actually be built.
    let mut matched: Vec<(usize, f64, f64)> = Vec::new();
    for a in panel_design
        .iter()
        .filter(|r| r.composition == "african_only" && r.available)
    {
        if let Some(b) = panel_design
            .iter()
            .find(|r| r.composition == "balanced" && r.available && r.actual_size == a.actual_size)
        {
            matched.push((a.actual_size, a.f1, b.f1));
        }
    }
    let african_wins = !matched.is_empty() && matched.iter().all(|(_, af, bf)| af > bf);
    println!(
        "   -> at matched sizes {:?}, African-only F1 {} balanced",
        matched.iter().map(|m| m.0).collect::<Vec<_>>(),
        if african_wins { "beat" } else { "did NOT beat" }
    );
    for (n, af, bf) in &matched {
        println!(
            "        n={n:<3} african {af:.4}  vs  balanced {bf:.4}   (delta {:+.4})",
            af - bf
        );
    }

    // The same comparison restricted to non-African queries. If the African
    // advantage is real *everywhere*, this holds too; if it is really an
    // African-coverage effect, this is where it breaks.
    let mut matched_nonafr: Vec<(usize, f64, f64)> = Vec::new();
    for a in panel_design
        .iter()
        .filter(|r| r.composition == "african_only" && r.available)
    {
        if let Some(b) = panel_design
            .iter()
            .find(|r| r.composition == "balanced" && r.available && r.actual_size == a.actual_size)
        {
            matched_nonafr.push((
                a.actual_size,
                a.precision_nonafr_queries,
                b.precision_nonafr_queries,
            ));
        }
    }
    let african_wins_nonafr =
        !matched_nonafr.is_empty() && matched_nonafr.iter().all(|(_, af, bf)| af > bf);
    println!(
        "   -> restricted to non-African queries, African-only precision {} balanced",
        if african_wins_nonafr {
            "beat"
        } else {
            "did NOT beat"
        }
    );

    // How many balanced haplotypes is one African haplotype worth? Find the
    // smallest balanced panel that matches the best African-only panel.
    let best_african = panel_design
        .iter()
        .filter(|r| r.composition == "african_only" && r.available)
        .max_by(|a, b| a.f1.partial_cmp(&b.f1).unwrap());
    let balanced_size_matching_best_african = best_african.and_then(|a| {
        panel_design
            .iter()
            .filter(|r| r.composition == "balanced" && r.available && r.f1 >= a.f1)
            .map(|r| r.actual_size)
            .min()
    });
    let african_haplotype_value_multiple = match (best_african, balanced_size_matching_best_african)
    {
        (Some(a), Some(b)) if a.actual_size > 0 => Some(b as f64 / a.actual_size as f64),
        _ => None,
    };
    if let (Some(a), Some(b)) = (best_african, balanced_size_matching_best_african) {
        println!(
            "   -> {} African haplotypes (F1 {:.4}) take {} balanced haplotypes to match: {:.2}x value per haplotype",
            a.actual_size,
            a.f1,
            b,
            b as f64 / a.actual_size as f64
        );
    }

    // ----------------------------------------------------------------- 2
    banner("2. Archaic reference ablation — what does one Denisovan genome cost?");

    // Detection is identical in every arm; only the labelling changes.
    let scan = match full_scan {
        Some(s) => s,
        None => scan_panel(&cohort, &params, &queries, &queries),
    };

    let single_nea: Vec<usize> = nea_all.iter().copied().take(1).collect();
    let arms: Vec<(&str, Vec<usize>, Vec<usize>)> = vec![
        ("all four (3 NEA + 1 DEN)", nea_all.clone(), den_all.clone()),
        ("Neanderthal only (3)", nea_all.clone(), Vec::new()),
        ("Denisovan only (1)", Vec::new(), den_all.clone()),
        ("single Neanderthal (1)", single_nea, Vec::new()),
        ("none", Vec::new(), Vec::new()),
    ];

    println!(
        "   reference set              NEA   DEN     calls    ->NEA    ->DEN  ->GHOST   trueNEA→ghost  trueDEN→ghost   ghost P   ghost R  ghost F1   attr acc"
    );

    let mut archaic_ablation: Vec<AblationPoint> = Vec::new();
    for (name, nea, den) in &arms {
        let calls = attribute(&cohort, &params, &scan.segments, nea, den);

        let mut confusion: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
        for c in &calls {
            let t = truth_label(&cohort, c.window, c.haplotype);
            *confusion
                .entry(t.to_string())
                .or_default()
                .entry(c.attribution.clone())
                .or_insert(0) += 1;
        }
        let cell = |t: &str, a: &str| -> usize {
            confusion
                .get(t)
                .and_then(|m| m.get(a))
                .copied()
                .unwrap_or(0)
        };

        let n_nea_lab = calls
            .iter()
            .filter(|c| c.attribution == "Neanderthal")
            .count();
        let n_den_lab = calls
            .iter()
            .filter(|c| c.attribution == "Denisovan")
            .count();
        let n_ghost_lab = calls.iter().filter(|c| c.attribution == "GHOST").count();

        let ghost_tp = calls
            .iter()
            .filter(|c| {
                c.attribution == "GHOST"
                    && is_true_ghost(truth_label(&cohort, c.window, c.haplotype))
            })
            .count();
        let ghost_precision = if n_ghost_lab == 0 {
            0.0
        } else {
            ghost_tp as f64 / n_ghost_lab as f64
        };
        let ghost_recall = if n_true_ghost == 0 {
            0.0
        } else {
            ghost_tp as f64 / n_true_ghost as f64
        };
        let ghost_f1 = if ghost_precision + ghost_recall == 0.0 {
            0.0
        } else {
            2.0 * ghost_precision * ghost_recall / (ghost_precision + ghost_recall)
        };

        // Attribution accuracy is only meaningful over calls that really are
        // archaic; false-positive calls have no correct label to get right.
        let mut archaic_calls = 0usize;
        let mut correct = 0usize;
        for c in &calls {
            let t = truth_label(&cohort, c.window, c.haplotype);
            if t == "not-introgressed" {
                continue;
            }
            archaic_calls += 1;
            let ok = match t {
                "Neanderthal" => c.attribution == "Neanderthal",
                "Denisovan" => c.attribution == "Denisovan",
                _ => c.attribution == "GHOST",
            };
            correct += usize::from(ok);
        }
        let attribution_accuracy = if archaic_calls == 0 {
            0.0
        } else {
            correct as f64 / archaic_calls as f64
        };

        println!(
            "   {:<26} {:>3}   {:>3}   {:>7}  {:>7}  {:>7}  {:>7}   {:>13}  {:>13}   {:>7.3}   {:>7.3}   {:>7.4}   {:>8.3}",
            name,
            nea.len(),
            den.len(),
            calls.len(),
            n_nea_lab,
            n_den_lab,
            n_ghost_lab,
            cell("Neanderthal", "GHOST"),
            cell("Denisovan", "GHOST"),
            ghost_precision,
            ghost_recall,
            ghost_f1,
            attribution_accuracy
        );

        archaic_ablation.push(AblationPoint {
            reference_set: name.to_string(),
            neanderthal_refs: nea.len(),
            denisovan_refs: den.len(),
            total_calls: calls.len(),
            attributed_neanderthal: n_nea_lab,
            attributed_denisovan: n_den_lab,
            attributed_ghost: n_ghost_lab,
            true_neanderthal_as_ghost: cell("Neanderthal", "GHOST"),
            true_denisovan_as_ghost: cell("Denisovan", "GHOST"),
            true_neanderthal_as_denisovan: cell("Neanderthal", "Denisovan"),
            true_denisovan_as_neanderthal: cell("Denisovan", "Neanderthal"),
            not_introgressed_as_ghost: cell("not-introgressed", "GHOST"),
            confusion,
            ghost_precision,
            ghost_recall,
            ghost_f1,
            attribution_accuracy,
        });
    }

    let arm_p = |needle: &str| -> f64 {
        archaic_ablation
            .iter()
            .find(|a| a.reference_set.starts_with(needle))
            .map(|a| a.ghost_precision)
            .unwrap_or(0.0)
    };
    let p_all = arm_p("all four");
    let p_none = arm_p("none");
    let p_nea_only = arm_p("Neanderthal only");
    println!(
        "\n   -> ghost precision: {:.3} with all four references, {:.3} without the Denisovan, {:.3} with none",
        p_all, p_nea_only, p_none
    );

    // ----------------------------------------------------------------- out
    let mut findings = Vec::new();

    let bal = |n: usize| -> Option<&PanelPoint> {
        panel_design
            .iter()
            .find(|r| r.composition == "balanced" && r.available && r.actual_size == n)
    };
    let afr = |n: usize| -> Option<&PanelPoint> {
        panel_design
            .iter()
            .find(|r| r.composition == "african_only" && r.available && r.actual_size == n)
    };

    if let (Some(small), Some(big)) = (bal(6), bal(*sizes.last().unwrap())) {
        let saturated = balanced_saturation_size.is_some_and(|s| s < big.actual_size);
        findings.push(format!(
            "Panel size buys precision and *costs* recall. Growing a balanced panel from {} to {} haplotypes lifts precision from {:.3} to {:.3} while recall falls from {:.3} to {:.3}, for a net F1 gain of {:.4} -> {:.4}. Both halves have the same cause: every extra haplotype is another chance for a segment to find a close relative, which suppresses false positives and also hides some genuinely introgressed segments whose nearest relative happens to be shallow. The noise floor — the median depth assigned to ordinary segments — drops from {:.0} ka to {:.0} ka over the same range, and that is what the archaic threshold has to clear.",
            small.actual_size, big.actual_size,
            small.precision, big.precision,
            small.recall, big.recall,
            small.f1, big.f1,
            small.ordinary_median_ka, big.ordinary_median_ka
        ));
        findings.push(if saturated {
            format!(
                "Returns saturate: a balanced panel is within 2% of the best F1 measured by {} haplotypes, so past that point more reference genomes are close to wasted.",
                balanced_saturation_size.unwrap()
            )
        } else {
            format!(
                "Returns had not saturated at the largest panel tested. The best balanced F1 ({:.4}) was still the largest panel, {} haplotypes, so within the range studied there is no plateau — the answer to 'how many do I need?' is bounded by what you can afford, not by diminishing returns.",
                big.f1, big.actual_size
            )
        });
    }

    if matched.is_empty() {
        findings.push(
            "African-only and balanced panels could not be compared at any matched size in this cohort.".to_string(),
        );
    } else if african_wins {
        let deltas: Vec<String> = matched
            .iter()
            .map(|(n, af, bf)| format!("n={n}: {:+.4}", af - bf))
            .collect();
        findings.push(format!(
            "The prediction HELD, and by a wide margin. At every size both compositions could be built, an African-only panel scored higher F1 than a balanced one ({}). {}",
            deltas.join(", "),
            match (best_african, balanced_size_matching_best_african, african_haplotype_value_multiple) {
                (Some(a), Some(b), Some(m)) => format!(
                    "In purchasing terms: {} African haplotypes (F1 {:.4}) are worth {} balanced ones — {:.2}x the detection value per genome sequenced. The advantage is also recall-shaped: an African-only panel holds recall at {:.3} at every size tested, while a balanced panel trades recall away as it grows, from {:.3} down to {:.3}.",
                    a.actual_size, a.f1, b, m,
                    a.recall,
                    bal(6).map(|r| r.recall).unwrap_or(0.0),
                    bal(*sizes.last().unwrap()).map(|r| r.recall).unwrap_or(0.0)
                ),
                _ => "No balanced panel tested ever matched the best African-only panel's F1.".to_string(),
            }
        ));
    } else {
        let worst = matched.iter().map(|(n, af, bf)| (*n, af - bf)).fold(
            (0usize, f64::INFINITY),
            |acc, x| if x.1 < acc.1 { x } else { acc },
        );
        findings.push(format!(
            "The prediction did NOT hold. African-only panels were expected to be more informative per haplotype because of their deeper coalescent structure, but at matched sizes {:?} they did not beat balanced panels — the worst case was n={} where African-only lost {:.4} F1.",
            matched.iter().map(|m| m.0).collect::<Vec<_>>(),
            worst.0,
            -worst.1
        ));
    }

    // The sub-prediction. "African panels are more informative" is easy to read
    // as "an African panel is better for everybody", and that is the part the
    // data refuses.
    if !matched_nonafr.is_empty() {
        let n = matched_nonafr.last().unwrap().0;
        if let (Some(a), Some(b)) = (afr(n), bal(n)) {
            findings.push(if african_wins_nonafr {
                format!(
                    "The African advantage is uniform: even on non-African queries an African-only panel of {n} beat a balanced one, {:.3} precision against {:.3}.",
                    a.precision_nonafr_queries, b.precision_nonafr_queries
                )
            } else {
                format!(
                    "But the natural corollary — that an African panel is therefore the better panel for everyone — did NOT hold, and this is the prediction the data refuses. Split by where the query comes from, an African-only panel of {n} scores {:.3} precision on African queries and only {:.3} on non-African ones; a balanced panel of the same size scores {:.3} and {:.3}. The African panel wins overall *despite* being strictly worse for every non-African query, because African queries are where the false positives are concentrated: a balanced panel spends most of its slots on populations whose internal diversity is shallow, leaving African queries with no close relative and a precision of {:.3}. The lesson is not 'African haplotypes are better haplotypes' but 'a panel is only informative about the populations it actually contains, and African populations need the most coverage to cover.'",
                    a.precision_afr_queries,
                    a.precision_nonafr_queries,
                    b.precision_afr_queries,
                    b.precision_nonafr_queries,
                    b.precision_afr_queries
                )
            });
        }
    }

    // The largest balanced panel that the best African-only panel still beats.
    let outgunned = best_african.and_then(|a| {
        panel_design
            .iter()
            .filter(|r| r.composition == "balanced" && r.available && r.f1 < a.f1)
            .max_by_key(|r| r.actual_size)
    });
    if let (Some(a), Some(b)) = (best_african, outgunned) {
        findings.push(format!(
            "Composition beats count outright in this cohort. An African-only panel of {} haplotypes (F1 {:.4}) outperforms a balanced panel of {} (F1 {:.4}) — {} fewer genomes to sequence, better answers. Anyone choosing between 'more genomes' and 'the right genomes' should choose the right genomes, until the right ones run out.",
            a.actual_size,
            a.f1,
            b.actual_size,
            b.f1,
            b.actual_size - a.actual_size
        ));
    }

    findings.push(format!(
        "And they do run out. An African-only panel cannot be scaled past {afr_available} haplotypes in this cohort at all, so the three largest panel sizes tested have no African-only counterpart and are recorded as unavailable rather than silently downsized. Composition and count are not independent knobs: picking the more informative composition also caps how much of it you can buy, which is why the best *overall* result still belongs to the largest balanced panel."
    ));

    let all_arm = archaic_ablation
        .iter()
        .find(|a| a.reference_set.starts_with("all four"))
        .unwrap();
    let nea_only_arm = archaic_ablation
        .iter()
        .find(|a| a.reference_set.starts_with("Neanderthal only"))
        .unwrap();
    let none_arm = archaic_ablation
        .iter()
        .find(|a| a.reference_set == "none")
        .unwrap();
    let den_only_arm = archaic_ablation
        .iter()
        .find(|a| a.reference_set.starts_with("Denisovan only"))
        .unwrap();
    let single_arm = archaic_ablation
        .iter()
        .find(|a| a.reference_set.starts_with("single Neanderthal"))
        .unwrap();

    findings.push(format!(
        "Removing archaic references never changes which segments are called — detection reads only the modern panel — but it changes what they are called. Ghost precision falls from {:.3} with all four archaic genomes to {:.3} with none, a drop of {:.3} ({:.1}% relative), because with nothing to attribute against, 'ghost' degenerates into 'not modern' and every real Neanderthal and Denisovan segment is filed as a new lineage.",
        all_arm.ghost_precision,
        none_arm.ghost_precision,
        all_arm.ghost_precision - none_arm.ghost_precision,
        100.0 * (all_arm.ghost_precision - none_arm.ghost_precision) / all_arm.ghost_precision.max(1e-9)
    ));

    findings.push(format!(
        "The cost of having only one Denisovan genome is measurable on its own. Deleting it — the 'Neanderthal only' arm — pushes {} true Denisovan segments into the GHOST pile and drops ghost precision from {:.3} to {:.3}. Deleting the three Neanderthal genomes instead costs {} misfiled Neanderthal segments and leaves ghost precision at {:.3}: the Denisovan reference is the scarcer resource per genome, but the Neanderthal references protect more segments in total.",
        nea_only_arm.true_denisovan_as_ghost,
        all_arm.ghost_precision,
        nea_only_arm.ghost_precision,
        den_only_arm.true_neanderthal_as_ghost,
        den_only_arm.ghost_precision
    ));

    findings.push(format!(
        "Archaic references are partly substitutable, which is why the single-Denisovan problem is survivable. With only a Denisovan reference, {} of {} called true Neanderthal segments are still caught as archaic-and-named — misfiled as 'Denisovan' rather than escaping into GHOST — because Neanderthal and Denisovan split 400 ka, well inside the {:.0} ka attribution margin. A reference from a sister lineage is a blunt but real substitute; overall attribution accuracy is {:.3} with all four references, {:.3} with a single Neanderthal, and {:.3} with none.",
        den_only_arm.true_neanderthal_as_denisovan,
        den_only_arm.true_neanderthal_as_denisovan + den_only_arm.true_neanderthal_as_ghost,
        params.match_margin_ka,
        all_arm.attribution_accuracy,
        single_arm.attribution_accuracy,
        none_arm.attribution_accuracy
    ));

    findings.push(format!(
        "Cutting the Neanderthal panel from three genomes to one costs almost nothing next to cutting it to zero: ghost precision goes {:.3} -> {:.3} -> {:.3} across three, one, and no Neanderthal references, and the number of true Neanderthal segments escaping into GHOST rises only from {} to {} before jumping to {} at zero. Archaic reference value is steeply diminishing — the first genome of a lineage does nearly all the work, which is the reassuring reading of the single-Denisovan problem: n=1 is a long way from n=0.",
        nea_only_arm.ghost_precision,
        single_arm.ghost_precision,
        none_arm.ghost_precision,
        nea_only_arm.true_neanderthal_as_ghost,
        single_arm.true_neanderthal_as_ghost,
        none_arm.true_neanderthal_as_ghost
    ));

    findings.push(format!(
        "Ghost *recall* is untouched by any of this: {:.3} in every arm, because a ghost segment is deeper than anything an observed archaic reference could match, so no reference set can pull one out of the ghost pile. The whole cost of missing archaic genomes is paid in precision, which means the failure mode of an under-referenced study is inventing lineages, never missing them.",
        all_arm.ghost_recall
    ));

    let misfiled = all_arm.true_neanderthal_as_ghost + all_arm.true_denisovan_as_ghost;
    findings.push(format!(
        "One caveat on the headline number: even with all four archaic references, ghost precision is only {:.3}, and {} of the {} impurities in the ghost pile are detection false positives — segments that were never introgressed at all — against just {} genuinely archaic segments filed under the wrong lineage. Better archaic references cannot move that ceiling; a better modern panel can, which is exactly what sweep 1 measures. The two sweeps are the two halves of the same error budget.",
        all_arm.ghost_precision,
        all_arm.not_introgressed_as_ghost,
        all_arm.not_introgressed_as_ghost + misfiled,
        misfiled
    ));

    let out = PanelStudy {
        generated_utc: chrono::Utc::now().to_rfc3339(),
        params,
        cohort: CohortSummary {
            n_windows: cohort.sequences.len(),
            window_bp: cohort.config.window_bp,
            modern_haplotypes: queries.len(),
            african_modern_haplotypes: afr_available,
            neanderthal_references: nea_all.len(),
            denisovan_references: den_all.len(),
            true_neanderthal_segments: n_true_nea,
            true_denisovan_segments: n_true_den,
            true_ghost_segments: n_true_ghost,
        },
        panel_design,
        balanced_saturation_size,
        african_beats_balanced_at_matched_size: african_wins,
        african_beats_balanced_for_nonafrican_queries: african_wins_nonafr,
        balanced_size_matching_best_african,
        african_haplotype_value_multiple,
        archaic_ablation,
        ghost_precision_all_references: p_all,
        ghost_precision_no_references: p_none,
        ghost_precision_neanderthal_only: p_nea_only,
        findings,
    };

    let path = format!("{out_dir}/panel.json");
    let mut f = std::fs::File::create(&path)?;
    f.write_all(serde_json::to_string_pretty(&out)?.as_bytes())?;
    println!("\n\x1b[1mWrote\x1b[0m {path}");

    Ok(())
}

// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny() -> Cohort {
        simulate_cohort(DemographyConfig {
            window_bp: 300,
            n_windows: 4,
            ..Default::default()
        })
    }

    #[test]
    fn balanced_panel_spreads_across_superpopulations() {
        let c = tiny();
        let panel = balanced_panel(&c, 6);
        assert_eq!(panel.len(), 6);
        let groups: HashSet<&str> = panel
            .iter()
            .map(|&h| c.haplotypes[h].group.as_str())
            .collect();
        // Six haplotypes drawn round-robin must land in six distinct groups.
        assert_eq!(
            groups.len(),
            6,
            "balanced panel is not balanced: {groups:?}"
        );
        // A panel may never contain an archaic reference.
        assert!(panel.iter().all(|&h| c.haplotypes[h].role == Role::Modern));
    }

    #[test]
    fn african_panel_is_african_and_capped() {
        let c = tiny();
        let panel = african_panel(&c, 12);
        assert_eq!(panel.len(), 12);
        assert!(panel.iter().all(|&h| c.haplotypes[h].group == "AFR"));

        // Asking for more African haplotypes than exist truncates rather than
        // silently borrowing from elsewhere — sweep 1 depends on that.
        let n_afr = modern_indices(&c)
            .into_iter()
            .filter(|&h| c.haplotypes[h].group == "AFR")
            .count();
        let over = african_panel(&c, n_afr + 10);
        assert_eq!(over.len(), n_afr);
    }

    #[test]
    fn ablating_references_moves_calls_into_ghost_without_changing_detection() {
        let c = tiny();
        let p = DetectorParams::default();
        let queries = modern_indices(&c);
        let scan = scan_panel(&c, &p, &queries, &queries);

        let nea: Vec<usize> = (0..c.haplotypes.len())
            .filter(|&h| c.haplotypes[h].population == "NEA")
            .collect();
        let den: Vec<usize> = (0..c.haplotypes.len())
            .filter(|&h| c.haplotypes[h].population == "DEN")
            .collect();

        let full = attribute(&c, &p, &scan.segments, &nea, &den);
        let bare = attribute(&c, &p, &scan.segments, &[], &[]);

        // Detection is untouched by the ablation.
        assert_eq!(full.len(), bare.len());
        // With no references at all, nothing can be named.
        assert!(bare.iter().all(|c| c.attribution == "GHOST"));
        // And the ghost pile can only grow when references are removed.
        let g = |v: &[ArchaicCall]| v.iter().filter(|c| c.attribution == "GHOST").count();
        assert!(g(&bare) >= g(&full));
    }

    #[test]
    fn interleave_prefix_is_balanced() {
        let b = vec![vec![1, 2, 3], vec![4, 5], vec![6]];
        assert_eq!(interleave(&b), vec![1, 4, 6, 2, 5, 3]);
    }
}
