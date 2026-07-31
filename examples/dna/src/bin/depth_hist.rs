//! `depth-hist` — the distribution the whole method rests on.
//!
//! For every modern segment in the simulated cohort, computes the coalescent
//! depth to its nearest relative anywhere in the modern panel, and buckets those
//! depths by what the segment actually is. This is the picture that shows *why*
//! a single threshold can separate archaic ancestry from ordinary variation, and
//! where the two distributions unavoidably overlap.
//!
//! Deliberately brute-force: this is the reference measurement the HNSW-indexed
//! detector is judged against, so it must not use the same shortcut.
//!
//! Usage: `depth-hist [output-dir]` (default `artifacts/data`).

use std::collections::BTreeMap;
use std::io::Write as _;

use rvdna::archaic::{
    divergence_to_tmrca_gens, gens_to_ka, raw_divergence, simulate_cohort, ArchaicSource,
    DemographyConfig, Role,
};
use serde::Serialize;

#[derive(Serialize)]
struct Series {
    label: String,
    /// Counts per bin, bins defined by `bin_edges_ka`.
    counts: Vec<usize>,
    n: usize,
    mean_ka: f64,
    median_ka: f64,
    p05_ka: f64,
    p95_ka: f64,
}

#[derive(Serialize)]
struct Overlap {
    /// Threshold in ka at which the two distributions are best separated.
    best_threshold_ka: f64,
    best_f1: f64,
    /// Fraction of ordinary segments that sit above the best threshold.
    ordinary_above: f64,
    /// Fraction of archaic segments that sit below it.
    archaic_below: f64,
}

/// One haplotype, in panel order, so the visualisation can lay the cohort out
/// the same way the engine indexes it.
#[derive(Serialize)]
struct RosterEntry {
    index: usize,
    id: String,
    population: String,
    group: String,
    role: String,
    /// How many of this haplotype's windows are archaic in truth.
    archaic_windows: usize,
}

#[derive(Serialize)]
struct Output {
    bin_edges_ka: Vec<f64>,
    series: Vec<Series>,
    overlap: Overlap,
    n_segments: usize,
    n_windows: usize,
    roster: Vec<RosterEntry>,
    /// `truth[window][haplotype]` as a source label, empty string for none.
    /// Encoded as one string per window with one character per haplotype:
    /// `.` none, `N` Neanderthal, `D` Denisovan, `A` Ghost-A, `B` Ghost-B.
    truth_map: Vec<String>,
    /// `depth_map[window][i]` — coalescent depth in ka, rounded, for the i-th
    /// *modern* haplotype in roster order. This is the raw measurement behind
    /// every figure, kept so the narrative can paint it directly rather than
    /// re-deriving it.
    depth_map: Vec<Vec<u32>>,
    /// Indices into `roster` that `depth_map` columns correspond to.
    depth_map_haplotypes: Vec<usize>,
    note: String,
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let i = ((sorted.len() - 1) as f64 * q).round() as usize;
    sorted[i]
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let out_dir = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "artifacts/data".into());
    std::fs::create_dir_all(&out_dir)?;

    // Same seed and config as `trace-rv`, so this is the same cohort.
    let config = DemographyConfig {
        n_windows: 360,
        ..Default::default()
    };
    let mu = config.mu;
    println!("simulating cohort (identical seed to trace-rv)...");
    let cohort = simulate_cohort(config);

    let modern: Vec<usize> = (0..cohort.haplotypes.len())
        .filter(|&h| cohort.haplotypes[h].role == Role::Modern)
        .collect();

    println!(
        "measuring nearest-relative depth for {} segments (brute force)...",
        cohort.sequences.len() * modern.len()
    );

    let mut by_label: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut depth_map: Vec<Vec<u32>> = Vec::with_capacity(cohort.sequences.len());
    let t0 = std::time::Instant::now();
    for w in 0..cohort.sequences.len() {
        let mut row: Vec<u32> = Vec::with_capacity(modern.len());
        for &h in &modern {
            let mut best = f64::INFINITY;
            for &o in &modern {
                if o == h {
                    continue;
                }
                let d = raw_divergence(&cohort.sequences[w][h], &cohort.sequences[w][o]);
                let t = gens_to_ka(divergence_to_tmrca_gens(d, mu));
                if t < best {
                    best = t;
                }
            }
            row.push(best.max(0.0).round() as u32);
            let label = cohort
                .truth
                .get(&(w, h))
                .map(|s| s.label().to_string())
                .unwrap_or_else(|| "No archaic ancestry".to_string());
            by_label.entry(label).or_default().push(best);
        }
        depth_map.push(row);
    }
    println!("   done in {:.1}s", t0.elapsed().as_secs_f64());

    // Bins from 0 to 2600 ka in 50 ka steps.
    let bin_w = 50.0;
    let n_bins = 52usize;
    let bin_edges_ka: Vec<f64> = (0..=n_bins).map(|i| i as f64 * bin_w).collect();

    let mut series = Vec::new();
    let mut ordinary: Vec<f64> = Vec::new();
    let mut archaic: Vec<f64> = Vec::new();

    for (label, mut vals) in by_label.into_iter() {
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut counts = vec![0usize; n_bins];
        for v in &vals {
            let b = ((v / bin_w) as usize).min(n_bins - 1);
            counts[b] += 1;
        }
        let n = vals.len();
        let mean = vals.iter().sum::<f64>() / n.max(1) as f64;
        if label == "No archaic ancestry" {
            ordinary.extend(&vals);
        } else {
            archaic.extend(&vals);
        }
        series.push(Series {
            label,
            counts,
            n,
            mean_ka: mean,
            median_ka: quantile(&vals, 0.5),
            p05_ka: quantile(&vals, 0.05),
            p95_ka: quantile(&vals, 0.95),
        });
    }

    // Where does a single depth threshold separate them best?
    let mut best_threshold = 0.0;
    let mut best_f1 = 0.0;
    let mut best_ord_above = 0.0;
    let mut best_arc_below = 0.0;
    let mut tau = 200.0;
    while tau <= 1600.0 {
        let tp = archaic.iter().filter(|d| **d >= tau).count();
        let fp = ordinary.iter().filter(|d| **d >= tau).count();
        let fneg = archaic.len() - tp;
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
        if f1 > best_f1 {
            best_f1 = f1;
            best_threshold = tau;
            best_ord_above = fp as f64 / ordinary.len().max(1) as f64;
            best_arc_below = fneg as f64 / archaic.len().max(1) as f64;
        }
        tau += 10.0;
    }

    println!(
        "   best single threshold {:.0} ka -> F1 {:.4}; {:.2}% of ordinary segments sit above it, {:.1}% of archaic below",
        best_threshold,
        best_f1,
        best_ord_above * 100.0,
        best_arc_below * 100.0
    );
    for s in &series {
        println!(
            "   {:<26} n={:<5} median {:>6.0} ka  [p05 {:>5.0} .. p95 {:>6.0}]",
            s.label, s.n, s.median_ka, s.p05_ka, s.p95_ka
        );
    }

    let roster: Vec<RosterEntry> = cohort
        .haplotypes
        .iter()
        .enumerate()
        .map(|(i, h)| RosterEntry {
            index: i,
            id: h.id.clone(),
            population: h.population.clone(),
            group: h.group.clone(),
            role: match h.role {
                Role::Modern => "modern".into(),
                _ => "archaic_reference".into(),
            },
            archaic_windows: (0..cohort.sequences.len())
                .filter(|&w| cohort.truth.contains_key(&(w, i)))
                .count(),
        })
        .collect();

    let truth_map: Vec<String> = (0..cohort.sequences.len())
        .map(|w| {
            (0..cohort.haplotypes.len())
                .map(|h| match cohort.truth.get(&(w, h)) {
                    Some(ArchaicSource::Neanderthal) => 'N',
                    Some(ArchaicSource::Denisovan) => 'D',
                    Some(ArchaicSource::GhostDeepAfrican) => 'A',
                    Some(ArchaicSource::SuperArchaic) => 'B',
                    None => '.',
                })
                .collect()
        })
        .collect();

    let out = Output {
        bin_edges_ka,
        n_segments: series.iter().map(|s| s.n).sum(),
        n_windows: cohort.sequences.len(),
        roster,
        truth_map,
        depth_map,
        depth_map_haplotypes: modern.clone(),
        series,
        overlap: Overlap {
            best_threshold_ka: best_threshold,
            best_f1,
            ordinary_above: best_ord_above,
            archaic_below: best_arc_below,
        },
        note: "Coalescent depth to the nearest relative anywhere in the modern panel, computed by exhaustive comparison rather than by the indexed detector, so it can serve as the reference the detector is judged against.".into(),
    };

    let path = format!("{out_dir}/depth-distribution.json");
    let mut f = std::fs::File::create(&path)?;
    f.write_all(serde_json::to_string_pretty(&out)?.as_bytes())?;
    println!("\nWrote {path}");
    Ok(())
}
