//! `trace-rv` — ghost-lineage recovery over the rvDNA engine.
//!
//! Runs the whole study end to end and writes machine-readable results:
//!
//! 1. simulate a hominin cohort under a structured coalescent with planted
//!    introgression from four archaic sources, two of which are "ghosts" that
//!    are never observed directly;
//! 2. verify the molecular clock the detector will rely on;
//! 3. evolve the detector's hyperparameters with Darwin mode (train split),
//!    and report on a held-out test split;
//! 4. spin the flywheel over the full genome until it goes dry;
//! 5. cluster the unattributed deep segments into ghost lineages and compare
//!    the recovered divergence times against the planted truth.
//!
//! Usage: `trace-rv [output-dir]` (default `artifacts/data`).

use std::collections::BTreeMap;
use std::io::Write as _;

use rvdna::archaic::{
    darwin, flywheel, gens_to_ka, simulate_cohort, DemographyConfig, DetectorParams, Role,
    TraceEngine,
};
use serde::Serialize;

#[derive(Serialize)]
struct ClockReport {
    n_points: usize,
    /// Pearson correlation between true and estimated TMRCA.
    correlation: f64,
    /// Mean relative error of the estimate.
    mean_relative_error: f64,
    median_true_ka: f64,
    median_estimated_ka: f64,
}

#[derive(Serialize)]
struct CohortReport {
    n_haplotypes: usize,
    n_modern: usize,
    n_archaic_references: usize,
    n_windows: usize,
    window_bp: usize,
    /// One simulated window carries as many segregating sites as this much real
    /// human sequence.
    equivalent_real_bp: usize,
    populations: BTreeMap<String, usize>,
    planted_segments: BTreeMap<String, usize>,
    planted_total: usize,
    archaic_fraction_of_segments: f64,
    /// Per reporting group, the fraction of that group's segments that are
    /// archaic in truth.
    truth_fraction_by_group: BTreeMap<String, f64>,
}

#[derive(Serialize)]
struct RetrievalReport {
    exact_comparisons: usize,
    bruteforce_comparisons: usize,
    speedup: f64,
}

#[derive(Serialize)]
struct LineageReport {
    cluster: usize,
    label: String,
    n_segments: usize,
    /// Lower-tail estimator of the divergence time — the one to quote.
    inferred_split_ka: f64,
    /// Mean coalescent depth, which is biased older than the split.
    mean_depth_ka: f64,
    median_depth_ka: f64,
    sd_ka: f64,
    truth_split_ka: f64,
    error_ka: f64,
    carriers_by_group: BTreeMap<String, usize>,
    /// What the planted truth says these segments actually were.
    truth_composition: BTreeMap<String, usize>,
    purity: f64,
}

#[derive(Serialize)]
struct AttributionRow {
    attribution: String,
    truth: String,
    count: usize,
}

#[derive(Serialize)]
struct FinalReport {
    title: String,
    generated_utc: String,
    source_study: SourceStudy,
    cohort: CohortReport,
    clock: ClockReport,
    darwin: darwin::Evolution,
    flywheel: flywheel::FlywheelRun,
    retrieval: RetrievalReport,
    lineages: Vec<LineageReport>,
    attribution_matrix: Vec<AttributionRow>,
    rvdna_container: RvdnaContainerReport,
    caveats: Vec<String>,
}

#[derive(Serialize)]
struct SourceStudy {
    title: String,
    method: String,
    journal: String,
    published: String,
    authors: Vec<String>,
    institutions: Vec<String>,
}

#[derive(Serialize)]
struct RvdnaContainerReport {
    haplotype: String,
    window: usize,
    sequence_bp: usize,
    container_bytes: usize,
    bits_per_base: f64,
    kmer_blocks: usize,
    roundtrip_exact: bool,
}

fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for i in 0..xs.len() {
        let dx = xs[i] - mx;
        let dy = ys[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 {
        0.0
    } else {
        sxy / (sxx.sqrt() * syy.sqrt())
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let out_dir = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "artifacts/data".into());
    std::fs::create_dir_all(&out_dir)?;

    let banner = |s: &str| {
        println!("\n\x1b[1m{s}\x1b[0m");
        println!("{}", "-".repeat(s.len()));
    };

    // ---------------------------------------------------------------- 1
    banner("1. Simulating hominin cohort (structured coalescent)");
    let config = DemographyConfig {
        n_windows: 360,
        ..Default::default()
    };
    let real_mu = 1.25e-8f64;
    let equivalent_real_bp = (config.window_bp as f64 * config.mu / real_mu) as usize;
    let t0 = std::time::Instant::now();
    let cohort = simulate_cohort(config);
    println!(
        "   {} haplotypes x {} windows x {} bp in {:.1}s",
        cohort.haplotypes.len(),
        cohort.sequences.len(),
        cohort.config.window_bp,
        t0.elapsed().as_secs_f64()
    );
    println!(
        "   mutation rate rescaled {:.0}x -> each window ~= {} kb of real human sequence",
        cohort.config.mu / real_mu,
        equivalent_real_bp / 1000
    );

    let n_modern = cohort
        .haplotypes
        .iter()
        .filter(|h| h.role == Role::Modern)
        .count();
    let mut populations: BTreeMap<String, usize> = BTreeMap::new();
    for h in &cohort.haplotypes {
        *populations.entry(h.population.clone()).or_insert(0) += 1;
    }

    let mut planted: BTreeMap<String, usize> = BTreeMap::new();
    for src in cohort.truth.values() {
        *planted.entry(src.label().to_string()).or_insert(0) += 1;
    }
    let planted_total = cohort.truth.len();
    let total_segments = cohort.sequences.len() * n_modern;

    let mut group_totals: BTreeMap<String, usize> = BTreeMap::new();
    let mut group_archaic: BTreeMap<String, usize> = BTreeMap::new();
    for (i, h) in cohort.haplotypes.iter().enumerate() {
        if h.role != Role::Modern {
            continue;
        }
        *group_totals.entry(h.group.clone()).or_insert(0) += cohort.sequences.len();
        let n = (0..cohort.sequences.len())
            .filter(|&w| cohort.truth.contains_key(&(w, i)))
            .count();
        *group_archaic.entry(h.group.clone()).or_insert(0) += n;
    }
    let truth_fraction_by_group: BTreeMap<String, f64> = group_totals
        .iter()
        .map(|(g, t)| {
            (
                g.clone(),
                *group_archaic.get(g).unwrap_or(&0) as f64 / *t as f64,
            )
        })
        .collect();

    for (k, v) in &planted {
        println!("   planted {v:>4} segments from {k}");
    }
    println!(
        "   total archaic burden: {:.2}% of all modern segments",
        100.0 * planted_total as f64 / total_segments as f64
    );

    // ---------------------------------------------------------------- 2
    banner("2. Verifying the molecular clock");
    let truths: Vec<f64> = cohort
        .calibration
        .iter()
        .map(|c| gens_to_ka(c.true_tmrca_gens))
        .collect();
    let ests: Vec<f64> = cohort
        .calibration
        .iter()
        .map(|c| gens_to_ka(c.estimated_tmrca_gens))
        .collect();
    let corr = pearson(&truths, &ests);
    let mre = truths
        .iter()
        .zip(&ests)
        .filter(|(t, _)| **t > 1.0)
        .map(|(t, e)| ((e - t) / t).abs())
        .sum::<f64>()
        / truths.iter().filter(|t| **t > 1.0).count().max(1) as f64;
    let clock = ClockReport {
        n_points: truths.len(),
        correlation: corr,
        mean_relative_error: mre,
        median_true_ka: median(truths.clone()),
        median_estimated_ka: median(ests.clone()),
    };
    println!(
        "   r = {:.4} over {} pairs; mean relative error {:.1}%",
        clock.correlation,
        clock.n_points,
        clock.mean_relative_error * 100.0
    );
    println!(
        "   median true TMRCA {:.0} ka vs estimated {:.0} ka",
        clock.median_true_ka, clock.median_estimated_ka
    );

    // ---------------------------------------------------------------- 3
    let storage = std::env::temp_dir()
        .join("trace_rv_hnsw")
        .to_string_lossy()
        .to_string();
    std::fs::remove_dir_all(&storage).ok();
    let mut engine = TraceEngine::new(&cohort, storage.clone());

    let n_win = cohort.sequences.len();
    let train: Vec<usize> = (0..n_win).filter(|w| w % 3 != 0).collect();
    let test: Vec<usize> = (0..n_win).filter(|w| w % 3 == 0).collect();

    banner("3. Darwin mode — evolving the detector");
    println!(
        "   train {} windows / test {} windows (held out)",
        train.len(),
        test.len()
    );
    let t0 = std::time::Instant::now();
    let evolution = darwin::evolve(
        &mut engine,
        DetectorParams::default(),
        &train,
        &test,
        6,
        4,
        0xDA_2E_11_60,
    )?;
    println!(
        "   {} proposals, {} accepted, {:.1}s",
        evolution.proposals,
        evolution.accepted,
        t0.elapsed().as_secs_f64()
    );
    println!(
        "   baseline  train F1 {:.4} | test F1 {:.4} (P {:.3} / R {:.3})",
        evolution.baseline_train_f1,
        evolution.baseline_test.f1,
        evolution.baseline_test.precision,
        evolution.baseline_test.recall
    );
    println!(
        "   evolved   train F1 {:.4} | test F1 {:.4} (P {:.3} / R {:.3})",
        evolution.evolved_train_f1,
        evolution.evolved_test.f1,
        evolution.evolved_test.precision,
        evolution.evolved_test.recall
    );
    for m in evolution.lineage.iter().filter(|m| m.accepted) {
        println!(
            "   gen {} kept {}: {} -> {} (F1 {:.4})",
            m.generation, m.gene, m.from, m.to, m.train_f1
        );
    }

    // ---------------------------------------------------------------- 4
    banner("4. Flywheel — detect, condense, re-index, repeat");
    let all: Vec<usize> = (0..n_win).collect();
    let wheel = flywheel::spin(&mut engine, &evolution.evolved, &all, 6, 2)?;
    for r in &wheel.rounds {
        println!(
            "   round {} | calls {:>4} ghost {:>4} | new refs {:>3} | F1 {:.4} (P {:.3} R {:.3}) | {} ms",
            r.round,
            r.calls,
            r.ghost_calls,
            r.new_ghost_references,
            r.scores.f1,
            r.scores.precision,
            r.scores.recall,
            r.elapsed_ms
        );
    }
    println!("   converged after round {}", wheel.converged_after);

    // ---------------------------------------------------------------- 5
    banner("5. Recovered ghost lineages");
    let mut calls = wheel.final_calls.clone();
    let clusters = engine.cluster_ghosts(&mut calls);

    // Which planted source dominates each recovered cluster?
    let truth_label = |w: usize, h: usize| -> String {
        cohort
            .truth
            .get(&(w, h))
            .map(|s| s.label().to_string())
            .unwrap_or_else(|| "none (false positive)".to_string())
    };

    let mut lineages = Vec::new();
    for c in &clusters {
        let mut composition: BTreeMap<String, usize> = BTreeMap::new();
        for call in calls
            .iter()
            .filter(|x| x.attribution == "GHOST" && x.ghost_cluster == Some(c.id))
        {
            *composition
                .entry(truth_label(call.window, call.haplotype))
                .or_insert(0) += 1;
        }
        let (dominant, dom_n) = composition
            .iter()
            .max_by_key(|(_, n)| **n)
            .map(|(k, n)| (k.clone(), *n))
            .unwrap_or_else(|| ("unknown".to_string(), 0));
        let purity = dom_n as f64 / c.n_segments.max(1) as f64;
        let truth_split_ka = if dominant.contains("super-archaic") {
            1800.0
        } else if dominant.contains("deep African") {
            800.0
        } else if dominant.contains("Neanderthal") || dominant.contains("Denisovan") {
            650.0
        } else {
            f64::NAN
        };
        println!(
            "   cluster {} | {:>4} segments | split {:.0} ka (p10) · depth {:.0} median / {:.0} mean | truth {:.0} ka | {} ({:.0}% pure)",
            c.id,
            c.n_segments,
            c.p10_split_ka,
            c.median_split_ka,
            c.mean_split_ka,
            truth_split_ka,
            dominant,
            purity * 100.0
        );
        print!("      carriers:");
        for (g, n) in &c.carriers_by_group {
            print!(" {g}={n}");
        }
        println!();
        lineages.push(LineageReport {
            cluster: c.id,
            label: dominant.clone(),
            n_segments: c.n_segments,
            inferred_split_ka: c.p10_split_ka,
            mean_depth_ka: c.mean_split_ka,
            median_depth_ka: c.median_split_ka,
            sd_ka: c.sd_split_ka,
            truth_split_ka,
            error_ka: (c.p10_split_ka - truth_split_ka).abs(),
            carriers_by_group: c.carriers_by_group.clone(),
            truth_composition: composition,
            purity,
        });
    }

    // Attribution confusion matrix.
    let mut matrix: BTreeMap<(String, String), usize> = BTreeMap::new();
    for call in &calls {
        let key = (
            call.attribution.clone(),
            truth_label(call.window, call.haplotype),
        );
        *matrix.entry(key).or_insert(0) += 1;
    }
    banner("6. Attribution vs planted truth");
    for ((a, t), n) in &matrix {
        println!("   called {a:<12} truth {t:<26} {n:>5}");
    }
    let attribution_matrix: Vec<AttributionRow> = matrix
        .into_iter()
        .map(|((attribution, truth), count)| AttributionRow {
            attribution,
            truth,
            count,
        })
        .collect();

    // ---------------------------------------------------------------- 7
    banner("7. HNSW retrieval efficiency");
    let final_run = engine.detect(&evolution.evolved, &all, 99)?;
    let speedup =
        final_run.bruteforce_comparisons as f64 / final_run.exact_comparisons.max(1) as f64;
    println!(
        "   {} exact divergence computations vs {} for an all-pairs scan ({:.1}x fewer)",
        final_run.exact_comparisons, final_run.bruteforce_comparisons, speedup
    );
    let retrieval = RetrievalReport {
        exact_comparisons: final_run.exact_comparisons,
        bruteforce_comparisons: final_run.bruteforce_comparisons,
        speedup,
    };

    // ---------------------------------------------------------------- 8
    banner("8. Packing a ghost carrier into the .rvdna container");
    let ghost_call = calls
        .iter()
        .find(|c| c.attribution == "GHOST")
        .cloned()
        .expect("no ghost segments were recovered");
    let seq_bytes = &cohort.sequences[ghost_call.window][ghost_call.haplotype];
    let seq_text = std::str::from_utf8(seq_bytes)?;
    let container = rvdna::rvdna::fasta_to_rvdna(seq_text, 11, 512, 500)?;
    let reader = rvdna::rvdna::RvdnaReader::from_bytes(container.clone())?;
    let restored = reader.read_sequence()?;
    let blocks = reader.read_kmer_vectors()?;
    let stats = reader.stats();
    let roundtrip_exact = restored.to_string() == seq_text;
    println!(
        "   {} window {} -> {} bytes ({:.2} bits/base), {} k-mer blocks, roundtrip {}",
        ghost_call.haplotype_id,
        ghost_call.window,
        container.len(),
        stats.bits_per_base,
        blocks.len(),
        if roundtrip_exact { "exact" } else { "LOSSY" }
    );

    // ---------------------------------------------------------------- write
    let report = FinalReport {
        title: "TRACE-rv: ghost hominin lineage recovery on the rvDNA engine".into(),
        generated_utc: chrono::Utc::now().to_rfc3339(),
        source_study: SourceStudy {
            title: "Recovering signatures of archaic hominin introgression using ancestral recombination graphs".into(),
            method: "TRACE (TRacking Archaic Contributions via ARG Estimation)".into(),
            journal: "Science".into(),
            published: "2026-07-30".into(),
            authors: vec![
                "Yulin Zhang".into(),
                "Arjun Biddanda".into(),
                "Sarah Johnson".into(),
                "Colm O'Dushlaine".into(),
                "Priya Moorjani".into(),
            ],
            institutions: vec![
                "UC Berkeley".into(),
                "Johns Hopkins University".into(),
                "54Gene, Inc.".into(),
            ],
        },
        cohort: CohortReport {
            n_haplotypes: cohort.haplotypes.len(),
            n_modern,
            n_archaic_references: cohort.haplotypes.len() - n_modern,
            n_windows: cohort.sequences.len(),
            window_bp: cohort.config.window_bp,
            equivalent_real_bp,
            populations,
            planted_segments: planted,
            planted_total,
            archaic_fraction_of_segments: planted_total as f64 / total_segments as f64,
            truth_fraction_by_group,
        },
        clock,
        darwin: evolution,
        flywheel: wheel,
        retrieval,
        lineages,
        attribution_matrix,
        rvdna_container: RvdnaContainerReport {
            haplotype: ghost_call.haplotype_id.clone(),
            window: ghost_call.window,
            sequence_bp: seq_bytes.len(),
            container_bytes: container.len(),
            bits_per_base: stats.bits_per_base,
            kmer_blocks: blocks.len(),
            roundtrip_exact,
        },
        caveats: vec![
            "Genomes are simulated under a coalescent model calibrated to published divergence times, not real 1000 Genomes / HGDP data.".into(),
            "The per-base mutation rate is rescaled so a short simulated window carries the segregating-site information of a much longer real one; all reported times use the true human clock.".into(),
            "Each planted archaic segment has exactly one carrier per window, matching the low per-locus frequency of real archaic tracts and keeping ground truth unambiguous.".into(),
            "This reproduces the logic of the published ARG-based method; it is not a reimplementation of TRACE itself and makes no claim about any living person's ancestry.".into(),
        ],
    };

    let path = format!("{out_dir}/trace-rv-report.json");
    let json = serde_json::to_string_pretty(&report)?;
    let mut f = std::fs::File::create(&path)?;
    f.write_all(json.as_bytes())?;

    // A compact per-call table for the visualisation layer.
    let calls_path = format!("{out_dir}/archaic-calls.json");
    let mut f = std::fs::File::create(&calls_path)?;
    f.write_all(serde_json::to_string(&calls)?.as_bytes())?;

    println!("\n\x1b[1mWrote\x1b[0m {path}");
    println!("\x1b[1mWrote\x1b[0m {calls_path}");
    std::fs::remove_dir_all(&storage).ok();
    Ok(())
}
