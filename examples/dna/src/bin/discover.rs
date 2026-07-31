//! `discover` — three questions the main study left open.
//!
//! The `trace-rv` run answers "does this work on one cohort?". These sweeps ask
//! where it stops working, which is the more useful thing to know.
//!
//! 1. **Resolution limit.** A ghost lineage is only resolvable as a *lineage*
//!    once enough of its segments survive to form a mode. How rare can it get?
//! 2. **Information limit.** Coalescent depth is estimated from segregating
//!    sites, so a window that is too short gives a noisy estimate. How much
//!    sequence does a window actually need?
//! 3. **Index fidelity.** The main pipeline lets an approximate nearest-
//!    neighbour index decide which comparisons to run. Does that approximation
//!    cost any calls against an exhaustive scan?
//!
//! Experiments 1 and 2 deliberately use an **exhaustive** detector, so the
//! limits they report are properties of the biology and the method, not of the
//! index. Experiment 3 is the one that isolates the index.
//!
//! Usage: `discover [output-dir]` (default `artifacts/data`).

use std::collections::HashSet;
use std::io::Write as _;

use rvdna::archaic::{
    divergence_to_tmrca_gens, gens_to_ka, raw_divergence, simulate_cohort, ArchaicCall,
    ArchaicSource, Cohort, DemographyConfig, DetectorParams, Role, Scores, TraceEngine,
};
use serde::Serialize;

/// Run the detector with no index at all: every modern segment against every
/// other modern segment, then against each archaic reference.
fn exhaustive_calls(cohort: &Cohort, p: &DetectorParams) -> Vec<ArchaicCall> {
    let mu = cohort.config.mu;
    let n_hap = cohort.haplotypes.len();
    let modern: Vec<usize> = (0..n_hap)
        .filter(|&h| cohort.haplotypes[h].role == Role::Modern)
        .collect();
    let nea: Vec<usize> = (0..n_hap)
        .filter(|&h| cohort.haplotypes[h].population == "NEA")
        .collect();
    let den: Vec<usize> = (0..n_hap)
        .filter(|&h| cohort.haplotypes[h].population == "DEN")
        .collect();

    let depth = |a: &[u8], b: &[u8]| gens_to_ka(divergence_to_tmrca_gens(raw_divergence(a, b), mu));

    let mut calls = Vec::new();
    for w in 0..cohort.sequences.len() {
        for &h in &modern {
            let seq = &cohort.sequences[w][h];
            let mut best = f64::INFINITY;
            for &o in &modern {
                if o != h {
                    best = best.min(depth(seq, &cohort.sequences[w][o]));
                }
            }
            if best < p.tau_archaic_ka {
                continue;
            }
            let d_nea = nea
                .iter()
                .map(|&r| depth(seq, &cohort.sequences[w][r]))
                .fold(f64::INFINITY, f64::min);
            let d_den = den
                .iter()
                .map(|&r| depth(seq, &cohort.sequences[w][r]))
                .fold(f64::INFINITY, f64::min);
            let attribution = if d_nea <= p.match_margin_ka && d_nea <= d_den {
                "Neanderthal"
            } else if d_den <= p.match_margin_ka {
                "Denisovan"
            } else {
                "GHOST"
            };
            calls.push(ArchaicCall {
                window: w,
                haplotype: h,
                haplotype_id: cohort.haplotypes[h].id.clone(),
                population: cohort.haplotypes[h].population.clone(),
                group: cohort.haplotypes[h].group.clone(),
                depth_to_modern_ka: best,
                depth_to_neanderthal_ka: d_nea,
                depth_to_denisovan_ka: d_den,
                attribution: attribution.to_string(),
                ghost_cluster: None,
                round: 0,
                primary: true,
                depth_to_ghost_ref_ka: None,
            });
        }
    }
    calls
}

fn score_against_truth(cohort: &Cohort, calls: &[ArchaicCall]) -> Scores {
    let called: HashSet<(usize, usize)> = calls.iter().map(|c| (c.window, c.haplotype)).collect();
    let truth: HashSet<(usize, usize)> = cohort
        .truth
        .keys()
        .copied()
        .filter(|(_, h)| cohort.haplotypes[*h].role == Role::Modern)
        .collect();
    let tp = called.intersection(&truth).count();
    Scores::compute(tp, called.len() - tp, truth.len() - tp)
}

/// Scale every super-archaic pulse by `factor`, leaving the rest alone.
fn config_with_ghostb_scale(factor: f64, seed: u64) -> DemographyConfig {
    let mut c = DemographyConfig {
        n_windows: 360,
        seed,
        ..Default::default()
    };
    for p in c.pulses.iter_mut() {
        if p.source == ArchaicSource::SuperArchaic {
            p.fraction *= factor;
        }
    }
    c
}

// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct ResolutionPoint {
    scale: f64,
    /// Segments of the super-archaic lineage actually planted.
    ghost_b_segments: usize,
    /// Its share of all modern segments, as a percentage.
    ghost_b_percent: f64,
    /// Did clustering split the ghost calls into two modes at all?
    two_modes_found: bool,
    /// Share of all Ghost-B segments that landed in the deepest mode.
    captured_by_deep_mode: f64,
    /// Median coalescent depth of the deepest mode, in ka.
    deep_mode_median_ka: f64,
    /// Oceanian share of the deep mode against the shallow one.
    oce_enrichment: f64,
    /// Was the deep mode identifiable as super-archaic at all?
    resolved: bool,
}

#[derive(Serialize)]
struct InformationPoint {
    window_bp: usize,
    equivalent_real_kb: usize,
    /// Expected segregating sites between two haplotypes at a 650 ka split.
    expected_sites_at_650ka: f64,
    precision: f64,
    recall: f64,
    f1: f64,
    /// Median depth the detector assigns to ordinary segments — the noise floor.
    ordinary_median_ka: f64,
}

#[derive(Serialize)]
struct FidelityReport {
    exhaustive_calls: usize,
    indexed_calls: usize,
    agreed: usize,
    missed_by_index: usize,
    added_by_index: usize,
    /// Share of the exhaustive detector's calls the indexed one reproduced.
    call_recall: f64,
    exhaustive_f1: f64,
    indexed_f1: f64,
    exact_comparisons_indexed: usize,
    exact_comparisons_exhaustive: usize,
    work_ratio: f64,
}

#[derive(Serialize)]
struct Discoveries {
    generated_utc: String,
    params: DetectorParams,
    resolution: Vec<ResolutionPoint>,
    resolution_threshold_segments: Option<usize>,
    information: Vec<InformationPoint>,
    information_knee_bp: Option<usize>,
    fidelity: FidelityReport,
    findings: Vec<String>,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let out_dir = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "artifacts/data".into());
    std::fs::create_dir_all(&out_dir)?;

    // The parameters Darwin settled on in the main run.
    let params = DetectorParams {
        tau_archaic_ka: 551.0,
        ..Default::default()
    };
    let real_mu = 1.25e-8f64;

    let banner = |s: &str| {
        println!("\n\x1b[1m{s}\x1b[0m");
        println!("{}", "-".repeat(s.len()));
    };

    // ----------------------------------------------------------------- 1
    banner("1. Resolution limit — how rare can a ghost lineage get?");
    println!("   scale  segments   %seg  modes  captured  deep median  OCE enrich  resolved");
    let mut resolution = Vec::new();
    for scale in [0.15f64, 0.3, 0.55, 1.0, 1.8, 3.0, 5.0] {
        let cohort = simulate_cohort(config_with_ghostb_scale(scale, 0xA5CE_57A1_2026));
        let n_modern = cohort
            .haplotypes
            .iter()
            .filter(|h| h.role == Role::Modern)
            .count();
        let gb: Vec<(usize, usize)> = cohort
            .truth
            .iter()
            .filter(|(_, s)| **s == ArchaicSource::SuperArchaic)
            .map(|(k, _)| *k)
            .collect();

        let mut calls = exhaustive_calls(&cohort, &params);
        let engine = TraceEngine::new(&cohort, std::env::temp_dir().to_string_lossy().to_string());
        let clusters = engine.cluster_ghosts(&mut calls);

        // The deepest mode is the candidate super-archaic one.
        let deep = clusters
            .iter()
            .enumerate()
            .max_by(|a, b| {
                a.1.median_split_ka
                    .partial_cmp(&b.1.median_split_ka)
                    .unwrap()
            })
            .map(|(i, _)| i);
        let (captured, deep_median, oce_enrich) = match deep {
            Some(di) => {
                let gb_set: HashSet<(usize, usize)> = gb.iter().copied().collect();
                let in_deep = calls
                    .iter()
                    .filter(|c| {
                        c.ghost_cluster == Some(clusters[di].id)
                            && gb_set.contains(&(c.window, c.haplotype))
                    })
                    .count();
                let share = |id: usize| {
                    let m: Vec<&ArchaicCall> = calls
                        .iter()
                        .filter(|c| c.ghost_cluster == Some(id))
                        .collect();
                    if m.is_empty() {
                        return 0.0;
                    }
                    m.iter().filter(|c| c.group == "OCE").count() as f64 / m.len() as f64
                };
                let deep_share = share(clusters[di].id);
                let other: f64 = clusters
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != di)
                    .map(|(_, c)| share(c.id))
                    .fold(0.0, f64::max);
                (
                    in_deep as f64 / gb.len().max(1) as f64,
                    clusters[di].median_split_ka,
                    if other > 0.0 { deep_share / other } else { 0.0 },
                )
            }
            None => (0.0, 0.0, 0.0),
        };

        // "Resolved" means: two modes exist, the deep one holds most of the
        // super-archaic segments, and it is Oceania-enriched — the same three
        // conditions the main study's conclusion rests on.
        let two_modes = clusters.len() >= 2;
        let resolved = two_modes && captured >= 0.5 && oce_enrich >= 1.5;

        println!(
            "   {:>5.2}  {:>8}  {:>5.3}  {:>5}  {:>7.0}%  {:>10.0}  {:>10.1}x  {}",
            scale,
            gb.len(),
            100.0 * gb.len() as f64 / (cohort.sequences.len() * n_modern) as f64,
            clusters.len(),
            captured * 100.0,
            deep_median,
            oce_enrich,
            if resolved { "yes" } else { "no" }
        );

        resolution.push(ResolutionPoint {
            scale,
            ghost_b_segments: gb.len(),
            ghost_b_percent: 100.0 * gb.len() as f64 / (cohort.sequences.len() * n_modern) as f64,
            two_modes_found: two_modes,
            captured_by_deep_mode: captured,
            deep_mode_median_ka: deep_median,
            oce_enrichment: oce_enrich,
            resolved,
        });
    }
    let resolution_threshold_segments = resolution
        .iter()
        .filter(|r| r.resolved)
        .map(|r| r.ghost_b_segments)
        .min();
    if let Some(t) = resolution_threshold_segments {
        println!("\n   -> a ghost lineage becomes resolvable at about {t} recovered segments");
    } else {
        println!("\n   -> no tested frequency resolved the lineage");
    }

    // ----------------------------------------------------------------- 2
    banner("2. Information limit — how much sequence does a window need?");
    println!("   window   real-eq   sites@650ka   precision  recall      F1   noise floor");
    let mut information = Vec::new();
    for window_bp in [250usize, 500, 1000, 2000, 4000] {
        let config = DemographyConfig {
            n_windows: 360,
            window_bp,
            ..Default::default()
        };
        let mu = config.mu;
        let cohort = simulate_cohort(config);
        let calls = exhaustive_calls(&cohort, &params);
        let scores = score_against_truth(&cohort, &calls);

        // The noise floor: what depth does the detector assign to ordinary
        // segments? If that climbs to meet the threshold, everything is archaic.
        let modern: Vec<usize> = (0..cohort.haplotypes.len())
            .filter(|&h| cohort.haplotypes[h].role == Role::Modern)
            .collect();
        let mut ordinary: Vec<f64> = Vec::new();
        for w in (0..cohort.sequences.len()).step_by(6) {
            for &h in &modern {
                if cohort.truth.contains_key(&(w, h)) {
                    continue;
                }
                let mut best = f64::INFINITY;
                for &o in &modern {
                    if o != h {
                        let d = raw_divergence(&cohort.sequences[w][h], &cohort.sequences[w][o]);
                        best = best.min(gens_to_ka(divergence_to_tmrca_gens(d, mu)));
                    }
                }
                ordinary.push(best);
            }
        }
        ordinary.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let ordinary_median = ordinary[ordinary.len() / 2];

        // 650 ka split -> 2 x 22,414 generations of independent mutation.
        let sites = 2.0 * (650_000.0 / 29.0) * mu * window_bp as f64;
        let eq_kb = (window_bp as f64 * mu / real_mu / 1000.0) as usize;

        println!(
            "   {:>6}  {:>6} kb  {:>11.0}  {:>10.3}  {:>6.3}  {:>6.4}  {:>8.0} ka",
            window_bp, eq_kb, sites, scores.precision, scores.recall, scores.f1, ordinary_median
        );

        information.push(InformationPoint {
            window_bp,
            equivalent_real_kb: eq_kb,
            expected_sites_at_650ka: sites,
            precision: scores.precision,
            recall: scores.recall,
            f1: scores.f1,
            ordinary_median_ka: ordinary_median,
        });
    }
    // The knee: the smallest window still within 5% of the best F1 achieved.
    let best_f1 = information.iter().map(|i| i.f1).fold(0.0, f64::max);
    let information_knee_bp = information
        .iter()
        .filter(|i| i.f1 >= best_f1 * 0.95)
        .map(|i| i.window_bp)
        .min();
    if let Some(k) = information_knee_bp {
        let eq = information
            .iter()
            .find(|i| i.window_bp == k)
            .map(|i| i.equivalent_real_kb)
            .unwrap_or(0);
        println!("\n   -> {k} bp (~{eq} kb of real sequence) is enough; more buys almost nothing");
    }

    // ----------------------------------------------------------------- 3
    banner("3. Index fidelity — what does the approximation cost?");
    let cohort = simulate_cohort(DemographyConfig {
        n_windows: 360,
        ..Default::default()
    });
    let storage = std::env::temp_dir()
        .join("discover_hnsw")
        .to_string_lossy()
        .to_string();
    std::fs::remove_dir_all(&storage).ok();

    let exhaustive = exhaustive_calls(&cohort, &params);
    let ex_scores = score_against_truth(&cohort, &exhaustive);

    let mut engine = TraceEngine::new(&cohort, storage.clone());
    let all: Vec<usize> = (0..cohort.sequences.len()).collect();
    let indexed_run = engine.detect(&params, &all, 0)?;
    std::fs::remove_dir_all(&storage).ok();

    let ex_set: HashSet<(usize, usize)> =
        exhaustive.iter().map(|c| (c.window, c.haplotype)).collect();
    let ix_set: HashSet<(usize, usize)> = indexed_run
        .calls
        .iter()
        .map(|c| (c.window, c.haplotype))
        .collect();
    let agreed = ex_set.intersection(&ix_set).count();

    let n_modern = cohort
        .haplotypes
        .iter()
        .filter(|h| h.role == Role::Modern)
        .count();
    let exhaustive_comparisons = cohort.sequences.len() * n_modern * (n_modern - 1 + 4);

    let fidelity = FidelityReport {
        exhaustive_calls: ex_set.len(),
        indexed_calls: ix_set.len(),
        agreed,
        missed_by_index: ex_set.len() - agreed,
        added_by_index: ix_set.len() - agreed,
        call_recall: agreed as f64 / ex_set.len().max(1) as f64,
        exhaustive_f1: ex_scores.f1,
        indexed_f1: indexed_run.scores.f1,
        exact_comparisons_indexed: indexed_run.exact_comparisons,
        exact_comparisons_exhaustive: exhaustive_comparisons,
        work_ratio: exhaustive_comparisons as f64 / indexed_run.exact_comparisons.max(1) as f64,
    };
    println!(
        "   exhaustive {} calls (F1 {:.4}) vs indexed {} calls (F1 {:.4})",
        fidelity.exhaustive_calls,
        fidelity.exhaustive_f1,
        fidelity.indexed_calls,
        fidelity.indexed_f1
    );
    println!(
        "   agreement {:.2}% | missed by index {} | added by index {}",
        fidelity.call_recall * 100.0,
        fidelity.missed_by_index,
        fidelity.added_by_index
    );
    println!(
        "   work: {} comparisons vs {} exhaustive ({:.1}x less)",
        fidelity.exact_comparisons_indexed,
        fidelity.exact_comparisons_exhaustive,
        fidelity.work_ratio
    );

    // ----------------------------------------------------------------- out
    let mut findings = Vec::new();
    if let Some(t) = resolution_threshold_segments {
        findings.push(format!(
            "A ghost lineage stops being resolvable as a distinct lineage below roughly {t} recovered segments. Above that it separates into its own mode, keeps most of its segments, and shows the carrier-geography enrichment that identifies its route; below it, the deep tail of the more common ghost swallows it entirely."
        ));
    }
    if let Some(k) = information_knee_bp {
        let pt = information.iter().find(|i| i.window_bp == k).unwrap();
        findings.push(format!(
            "Detection saturates at about {k} bp per window ({} kb of real-sequence information, ~{:.0} expected segregating sites at a 650 ka split). Longer windows buy almost nothing; shorter ones fail fast, because the coalescent-depth estimate turns to noise before the threshold does.",
            pt.equivalent_real_kb, pt.expected_sites_at_650ka
        ));
    }
    findings.push(format!(
        "Approximate nearest-neighbour retrieval reproduces {:.1}% of the calls an exhaustive scan makes, for {:.1}x less comparison work, and the two detectors land within {:.4} F1 of each other. The index is choosing what to compare, not what to conclude — and this is the measurement that shows it.",
        fidelity.call_recall * 100.0,
        fidelity.work_ratio,
        (fidelity.exhaustive_f1 - fidelity.indexed_f1).abs()
    ));

    let out = Discoveries {
        generated_utc: chrono::Utc::now().to_rfc3339(),
        params,
        resolution,
        resolution_threshold_segments,
        information,
        information_knee_bp,
        fidelity,
        findings,
    };

    let path = format!("{out_dir}/discoveries.json");
    let mut f = std::fs::File::create(&path)?;
    f.write_all(serde_json::to_string_pretty(&out)?.as_bytes())?;
    println!("\n\x1b[1mWrote\x1b[0m {path}");

    Ok(())
}
