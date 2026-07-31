//! `linkage-study` — does looking sideways along the chromosome help?
//!
//! The per-window detector in [`rvdna::archaic`] judges each segment alone. This
//! binary asks whether the *linkage* structure of introgression — the fact that
//! donated ancestry arrives in contiguous tracts, not as a scatter of
//! independent loci — is worth exploiting, and by how much.
//!
//! The experiment:
//!
//! 1. Simulate a cohort with `tract_mean_windows` in the 4–6 range on every
//!    pulse, so the planted truth actually has tract structure to find. (The
//!    default cohort has `tract_mean_windows = 1.0`, i.e. unlinked loci, where a
//!    linkage detector has nothing to work with by construction.)
//! 2. Compute the coalescent depth grid once, with the HNSW-indexed retrieval
//!    path, and run **both** decision rules over the *same* grid. Any difference
//!    in the numbers is then attributable to the decision rule alone, not to
//!    different evidence or a different random draw.
//! 3. Score both against planted truth; reconstruct called tracts; compare the
//!    recovered tract-length distribution against the planted one; and invert
//!    mean tract length into a time since admixture.
//! 4. Hold out half the windows, tune the linkage rule on the other half, and
//!    report the held-out numbers separately — so the headline claim is not a
//!    number that was fitted on the data it is quoted against.
//!
//! Usage: `linkage-study [output-dir]` (default `artifacts/data`).
//! Writes `<output-dir>/linkage.json`.

use std::collections::BTreeMap;
use std::io::Write as _;

use rvdna::archaic::{
    simulate_cohort, ArchaicSource, DemographyConfig, DetectorParams, Role, Scores, TraceEngine,
};
use rvdna::linkage::{
    admixture_time_from_tracts, call_per_window, compare_tract_lengths, planted_tract_lengths,
    planted_tracts_in, reconstruct_tracts, run_linkage, score_calls, AdmixtureTimeEstimate,
    LinkageParams, MorganScale, TractComparison, TractLengthStats,
};
use serde::Serialize;

/// Mean planted tract length, in windows, applied to every pulse.
const TRACT_MEAN_WINDOWS: f64 = 5.0;

#[derive(Serialize)]
struct CohortReport {
    n_haplotypes: usize,
    n_modern: usize,
    n_windows: usize,
    window_bp: usize,
    tract_mean_windows_requested: f64,
    planted_tracts: usize,
    planted_archaic_windows: usize,
    archaic_fraction_of_segments: f64,
    planted_by_source: BTreeMap<String, usize>,
    planted_tract_lengths: TractLengthStats,
}

#[derive(Serialize)]
struct DetectorReport {
    name: String,
    rule: String,
    scores: Scores,
    n_called_windows: usize,
    n_called_tracts: usize,
    called_tract_lengths: TractLengthStats,
    tract_comparison: TractComparison,
    admixture_time: Option<AdmixtureTimeEstimate>,
}

#[derive(Serialize)]
struct Delta {
    precision: f64,
    recall: f64,
    f1: f64,
    /// True if the linkage rule's F1 is higher than the per-window rule's.
    linkage_wins: bool,
    /// Windows the linkage rule added on top of the per-window call set.
    rescued_windows: usize,
    rescued_true: usize,
    rescued_precision: f64,
}

#[derive(Serialize)]
struct SweepRow {
    tau_rescue_ka: f64,
    tau_context_ka: f64,
    half_width: usize,
    decay: f64,
    train_f1: f64,
    test_precision: f64,
    test_recall: f64,
    test_f1: f64,
}

#[derive(Serialize)]
struct HeldOut {
    n_train_windows: usize,
    n_test_windows: usize,
    tuned_params: LinkageParams,
    per_window_test: Scores,
    linkage_test: Scores,
    linkage_wins_on_held_out: bool,
    sweep: Vec<SweepRow>,
    note: String,
}

#[derive(Serialize)]
struct ControlReport {
    note: String,
    per_window: Scores,
    linkage: Scores,
    linkage_wins: bool,
}

#[derive(Serialize)]
struct Reproducibility {
    cohort_deterministic: bool,
    hnsw_deterministic_across_processes: bool,
    note: String,
}

#[derive(Serialize)]
struct Report {
    title: String,
    generated_utc: String,
    question: String,
    cohort: CohortReport,
    detector_params: DetectorParams,
    linkage_params: LinkageParams,
    morgan_scale: MorganScale,
    coordinate_convention: String,
    per_window: DetectorReport,
    linkage: DetectorReport,
    delta: Delta,
    admixture_time_truth: Option<AdmixtureTimeEstimate>,
    held_out: HeldOut,
    unlinked_control: ControlReport,
    findings: Vec<String>,
    caveats: Vec<String>,
    reproducibility: Reproducibility,
}

fn tracted_config(n_windows: usize, mean: f64) -> DemographyConfig {
    let mut c = DemographyConfig {
        n_windows,
        ..Default::default()
    };
    for p in c.pulses.iter_mut() {
        p.tract_mean_windows = mean;
    }
    c
}

fn source_label(s: &ArchaicSource) -> String {
    s.label().to_string()
}

/// Build a full [`DetectorReport`] from a boolean call grid.
fn report_for(
    name: &str,
    rule: &str,
    called: &[Vec<bool>],
    windows: &[usize],
    modern: &[usize],
    cohort: &rvdna::archaic::Cohort,
    planted_lengths: &[usize],
    scale: &MorganScale,
) -> DetectorReport {
    let scores = score_calls(called, windows, modern, &cohort.truth);
    let tracts = reconstruct_tracts(called, windows, modern);
    let lengths: Vec<usize> = tracts.iter().map(|t| t.length).collect();
    let n_called: usize = called
        .iter()
        .map(|r| r.iter().filter(|b| **b).count())
        .sum();

    DetectorReport {
        name: name.to_string(),
        rule: rule.to_string(),
        scores,
        n_called_windows: n_called,
        n_called_tracts: tracts.len(),
        called_tract_lengths: TractLengthStats::from_lengths(&lengths),
        tract_comparison: compare_tract_lengths(planted_lengths, &lengths),
        admixture_time: admixture_time_from_tracts(&lengths, scale, name),
    }
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
    banner("1. Simulating a cohort with linked introgression tracts");
    let config = tracted_config(360, TRACT_MEAN_WINDOWS);
    let scale = MorganScale::from_config(&config);
    let t0 = std::time::Instant::now();
    let cohort = simulate_cohort(config);
    let n_windows = cohort.sequences.len();
    let modern_count = cohort
        .haplotypes
        .iter()
        .filter(|h| h.role == Role::Modern)
        .count();
    println!(
        "   {} haplotypes x {} windows in {:.1}s; {} tracts planted covering {} windows",
        cohort.haplotypes.len(),
        n_windows,
        t0.elapsed().as_secs_f64(),
        cohort.tracts.len(),
        cohort.truth.len()
    );
    println!("   {}", scale.convention());

    let mut planted_by_source: BTreeMap<String, usize> = BTreeMap::new();
    for t in &cohort.tracts {
        *planted_by_source
            .entry(source_label(&t.source))
            .or_insert(0) += t.length;
    }

    let all_windows: Vec<usize> = (0..n_windows).collect();
    let planted_lengths = planted_tract_lengths(&cohort, &all_windows);
    let planted_stats = TractLengthStats::from_lengths(&planted_lengths);
    println!(
        "   planted tract length: mean {:.2} windows, median {:.0}, max {}",
        planted_stats.mean_windows, planted_stats.median_windows, planted_stats.max_windows
    );

    // ---------------------------------------------------------------- 2
    banner("2. Depth grid (one HNSW-indexed pass, shared by both detectors)");
    let storage = std::env::temp_dir().join(format!("trace_rv_linkage_{}", std::process::id()));
    std::fs::remove_dir_all(&storage).ok();
    let mut engine = TraceEngine::new(&cohort, storage.to_string_lossy().to_string());
    let dp = DetectorParams::default();
    let t1 = std::time::Instant::now();
    let grid = engine.depth_grid(&dp, &all_windows)?;
    let modern = engine.modern_haplotypes();
    println!(
        "   {} x {} depths in {:.1}s ({} exact divergence computations)",
        grid.len(),
        modern.len(),
        t1.elapsed().as_secs_f64(),
        engine.exact_comparisons
    );

    // Consistency guard: thresholding the grid must reproduce `detect()`.
    let det = engine.detect(&dp, &all_windows, 0)?;
    let base_called = call_per_window(&grid, dp.tau_archaic_ka);
    let base_scores = score_calls(&base_called, &all_windows, &modern, &cohort.truth);
    assert_eq!(
        det.scores.true_positives, base_scores.true_positives,
        "depth_grid thresholding disagrees with detect(); the two detectors would \
         not be reading the same evidence"
    );
    println!("   grid agrees with detect(): {} calls", det.calls.len());
    std::fs::remove_dir_all(&storage).ok();

    // ---------------------------------------------------------------- 3
    banner("3. Per-window vs linkage-aware, on the same grid");
    let lp = LinkageParams::from_detector(&dp);
    let linkage = run_linkage(&grid, &all_windows, &modern, &cohort.truth, &lp);
    let linked_called = rvdna::linkage::call_linkage(&grid, &all_windows, &lp);

    let per_window_report = report_for(
        "per-window",
        "call iff this window's own coalescent depth >= tau_archaic_ka",
        &base_called,
        &all_windows,
        &modern,
        &cohort,
        &planted_lengths,
        &scale,
    );
    let linkage_report = report_for(
        "linkage-aware",
        "call iff own depth >= tau_archaic_ka, OR own depth >= tau_rescue_ka and the \
         decay-weighted geometric mean depth of the neighbouring windows (self excluded) \
         >= tau_context_ka",
        &linked_called,
        &all_windows,
        &modern,
        &cohort,
        &planted_lengths,
        &scale,
    );

    let pw = per_window_report.scores;
    let lk = linkage_report.scores;
    println!(
        "   per-window    P {:.4}  R {:.4}  F1 {:.4}   (tp {} fp {} fn {})",
        pw.precision, pw.recall, pw.f1, pw.true_positives, pw.false_positives, pw.false_negatives
    );
    println!(
        "   linkage-aware P {:.4}  R {:.4}  F1 {:.4}   (tp {} fp {} fn {})",
        lk.precision, lk.recall, lk.f1, lk.true_positives, lk.false_positives, lk.false_negatives
    );
    println!(
        "   rescued {} borderline windows, {} of them true ({:.1}% precision on the rescues)",
        linkage.rescued_windows,
        linkage.rescued_true,
        100.0 * linkage.rescued_true as f64 / linkage.rescued_windows.max(1) as f64
    );

    let delta = Delta {
        precision: lk.precision - pw.precision,
        recall: lk.recall - pw.recall,
        f1: lk.f1 - pw.f1,
        linkage_wins: lk.f1 > pw.f1,
        rescued_windows: linkage.rescued_windows,
        rescued_true: linkage.rescued_true,
        rescued_precision: linkage.rescued_true as f64 / linkage.rescued_windows.max(1) as f64,
    };

    // ---------------------------------------------------------------- 4
    banner("4. Tract lengths and the admixture clock");
    let truth_time = admixture_time_from_tracts(&planted_lengths, &scale, "planted truth");
    for (label, est) in [
        ("planted truth", &truth_time),
        ("per-window", &per_window_report.admixture_time),
        ("linkage-aware", &linkage_report.admixture_time),
    ] {
        if let Some(e) = est {
            println!(
                "   {label:<14} mean tract {:.2} windows = {:.4e} M -> t = {:.0} generations \
                 ({:.1} ka)  [95% {:.0}-{:.0} gen, n = {}]",
                e.mean_tract_windows,
                e.mean_tract_morgans,
                e.generations,
                e.ka,
                e.generations_ci95.0,
                e.generations_ci95.1,
                e.n_tracts
            );
        }
    }
    println!(
        "   fragmentation: per-window {:.2}x, linkage {:.2}x tracts per planted tract",
        per_window_report.tract_comparison.fragmentation,
        linkage_report.tract_comparison.fragmentation
    );

    // ---------------------------------------------------------------- 5
    banner("5. Held-out check (tune on train windows, report on test windows)");
    let train: Vec<usize> = (0..n_windows / 2).collect();
    let test: Vec<usize> = (n_windows / 2..n_windows).collect();
    let sub = |ws: &[usize]| -> Vec<Vec<f64>> { ws.iter().map(|&w| grid[w].clone()).collect() };
    let train_grid = sub(&train);
    let test_grid = sub(&test);

    let mut sweep: Vec<SweepRow> = Vec::new();
    let mut best: Option<(f64, LinkageParams)> = None;
    for rescue_frac in [0.60f64, 0.66, 0.72, 0.80, 0.88] {
        for context_frac in [0.40f64, 0.50, 0.55, 0.65, 0.75] {
            for half_width in [1usize, 2, 3, 5] {
                for decay in [0.4f64, 0.6, 0.85] {
                    let mut p = LinkageParams::from_detector(&dp);
                    p.tau_rescue_ka = rescue_frac * dp.tau_archaic_ka;
                    p.tau_context_ka = context_frac * dp.tau_archaic_ka;
                    p.half_width = half_width;
                    p.decay = decay;

                    let tr = run_linkage(&train_grid, &train, &modern, &cohort.truth, &p);
                    let te = run_linkage(&test_grid, &test, &modern, &cohort.truth, &p);
                    sweep.push(SweepRow {
                        tau_rescue_ka: p.tau_rescue_ka,
                        tau_context_ka: p.tau_context_ka,
                        half_width,
                        decay,
                        train_f1: tr.scores.f1,
                        test_precision: te.scores.precision,
                        test_recall: te.scores.recall,
                        test_f1: te.scores.f1,
                    });
                    if best.as_ref().map_or(true, |(f, _)| tr.scores.f1 > *f) {
                        best = Some((tr.scores.f1, p));
                    }
                }
            }
        }
    }
    let tuned = best.map(|(_, p)| p).unwrap_or(lp);
    let pw_test = score_calls(
        &call_per_window(&test_grid, dp.tau_archaic_ka),
        &test,
        &modern,
        &cohort.truth,
    );
    let lk_test = run_linkage(&test_grid, &test, &modern, &cohort.truth, &tuned).scores;
    println!(
        "   tuned on {} train windows: tau_rescue {:.0} ka, tau_context {:.0} ka, \
         half_width {}, decay {:.2}",
        train.len(),
        tuned.tau_rescue_ka,
        tuned.tau_context_ka,
        tuned.half_width,
        tuned.decay
    );
    println!(
        "   held-out per-window    P {:.4}  R {:.4}  F1 {:.4}",
        pw_test.precision, pw_test.recall, pw_test.f1
    );
    println!(
        "   held-out linkage-aware P {:.4}  R {:.4}  F1 {:.4}",
        lk_test.precision, lk_test.recall, lk_test.f1
    );

    // ---------------------------------------------------------------- 6
    banner("6. Control: the same rule on an UNLINKED cohort (tract_mean = 1.0)");
    // If the linkage rule helps on a tracted cohort it must *not* help on a
    // cohort with no tract structure, or the gain is an artefact of the
    // threshold relaxation rather than of linkage.
    let ctrl_cohort = simulate_cohort(DemographyConfig {
        n_windows: 180,
        ..Default::default()
    });
    let ctrl_storage =
        std::env::temp_dir().join(format!("trace_rv_linkage_ctrl_{}", std::process::id()));
    std::fs::remove_dir_all(&ctrl_storage).ok();
    let mut ctrl_engine =
        TraceEngine::new(&ctrl_cohort, ctrl_storage.to_string_lossy().to_string());
    let ctrl_windows: Vec<usize> = (0..ctrl_cohort.sequences.len()).collect();
    let ctrl_grid = ctrl_engine.depth_grid(&dp, &ctrl_windows)?;
    let ctrl_modern = ctrl_engine.modern_haplotypes();
    std::fs::remove_dir_all(&ctrl_storage).ok();

    let ctrl_pw = score_calls(
        &call_per_window(&ctrl_grid, dp.tau_archaic_ka),
        &ctrl_windows,
        &ctrl_modern,
        &ctrl_cohort.truth,
    );
    let ctrl_lk = run_linkage(
        &ctrl_grid,
        &ctrl_windows,
        &ctrl_modern,
        &ctrl_cohort.truth,
        &lp,
    )
    .scores;
    println!(
        "   unlinked per-window    P {:.4}  R {:.4}  F1 {:.4}",
        ctrl_pw.precision, ctrl_pw.recall, ctrl_pw.f1
    );
    println!(
        "   unlinked linkage-aware P {:.4}  R {:.4}  F1 {:.4}",
        ctrl_lk.precision, ctrl_lk.recall, ctrl_lk.f1
    );

    // ---------------------------------------------------------------- 7
    banner("7. Findings");
    let mut findings: Vec<String> = Vec::new();

    findings.push(format!(
        "On a cohort where introgression arrives in tracts of mean {:.1} windows, the \
         linkage-aware rule scored precision {:.3}, recall {:.3}, F1 {:.3}, against the \
         per-window detector's precision {:.3}, recall {:.3}, F1 {:.3} — a change in F1 of \
         {:+.4}.",
        planted_stats.mean_windows,
        lk.precision,
        lk.recall,
        lk.f1,
        pw.precision,
        pw.recall,
        pw.f1,
        delta.f1
    ));

    findings.push(if delta.f1 > 0.0 {
        format!(
            "Looking sideways along the chromosome helped: F1 rose by {:+.4}. The rule is \
             monotone by construction, so recall could only rise ({:+.4}); the question was \
             whether precision would survive the extra calls, and it moved {:+.4}.",
            delta.f1, delta.recall, delta.precision
        )
    } else if delta.f1 == 0.0 {
        format!(
            "Looking sideways along the chromosome changed nothing: F1 was identical to four \
             decimal places. The rescue band added {} windows.",
            delta.rescued_windows
        )
    } else {
        format!(
            "NEGATIVE RESULT: looking sideways along the chromosome did not pay. F1 fell by \
             {:.4}. Recall rose {:+.4}, as it must for a rule that only adds calls, but the \
             {} windows it added were only {:.1}% true, and the precision cost ({:+.4}) \
             outweighed the recall gain.",
            -delta.f1,
            delta.recall,
            delta.rescued_windows,
            100.0 * delta.rescued_precision,
            delta.precision
        )
    });

    findings.push(format!(
        "Of the {} borderline windows the linkage rule rescued, {} were genuinely \
         introgressed ({:.1}%). For comparison, the per-window detector's overall precision \
         is {:.1}%, so the rescued windows are {} than the calls it was already making.",
        delta.rescued_windows,
        delta.rescued_true,
        100.0 * delta.rescued_precision,
        100.0 * pw.precision,
        if delta.rescued_precision > pw.precision {
            "cleaner"
        } else {
            "dirtier"
        }
    ));

    findings.push(format!(
        "Neither detector recovers tracts at their planted length. The planted tracts average \
         {:.2} windows; the per-window detector recovers {:.2} and the linkage-aware rule \
         {:.2}. Both fragment real tracts rather than fusing them ({:.2}x and {:.2}x as many \
         tracts as were planted), because a single missed window in the middle of a run \
         splits it in two.",
        planted_stats.mean_windows,
        per_window_report.called_tract_lengths.mean_windows,
        linkage_report.called_tract_lengths.mean_windows,
        per_window_report.tract_comparison.fragmentation,
        linkage_report.tract_comparison.fragmentation
    ));

    if let (Some(t), Some(l)) = (&truth_time, &linkage_report.admixture_time) {
        findings.push(format!(
            "Inverting mean tract length through t = 1/E[L] Morgans gives {:.0} generations \
             ({:.1} ka) from the planted tracts and {:.0} generations ({:.1} ka) from the \
             linkage-aware calls. The recovered estimate is biased {} because fragmentation \
             shortens the measured tracts, and a shorter mean tract always reads as an older \
             pulse.",
            t.generations,
            t.ka,
            l.generations,
            l.ka,
            if l.generations > t.generations {
                "older"
            } else {
                "younger"
            }
        ));
    }

    findings.push(format!(
        "The admixture-time numbers are only as good as the windows-to-Morgans convention, \
         which is asserted rather than measured: {}",
        scale.convention()
    ));

    findings.push(format!(
        "On a held-out half of the genome, with the linkage parameters tuned only on the \
         other half, the linkage rule scored F1 {:.4} against the per-window rule's {:.4} \
         ({:+.4}). {}",
        lk_test.f1,
        pw_test.f1,
        lk_test.f1 - pw_test.f1,
        if lk_test.f1 > pw_test.f1 {
            "The advantage survives tuning on data it is not quoted against."
        } else {
            "The advantage does not survive being tuned on one half and quoted on the other."
        }
    ));

    findings.push(format!(
        "Control: on an UNLINKED cohort (tract_mean_windows = 1.0, the published default) the \
         same rule moved F1 from {:.4} to {:.4} ({:+.4}). {}",
        ctrl_pw.f1,
        ctrl_lk.f1,
        ctrl_lk.f1 - ctrl_pw.f1,
        if ctrl_lk.f1 >= ctrl_pw.f1 {
            "It does not hurt there, but with no tract structure to exploit there is nothing \
             for it to gain either — any gain seen on the tracted cohort that also appears \
             here is threshold relaxation, not linkage."
        } else {
            "It actively hurts there, which is the expected sign: with no tract structure, \
             relaxing the threshold on the strength of neighbouring windows can only add noise."
        }
    ));

    for f in &findings {
        println!("   - {f}");
    }

    // ---------------------------------------------------------------- 8
    let report = Report {
        title: "Linkage-aware archaic introgression detection".to_string(),
        generated_utc: chrono::Utc::now().to_rfc3339(),
        question: "Introgressed ancestry arrives in contiguous tracts. Does a detector that \
                   smooths coalescent depth across neighbouring windows before thresholding \
                   beat one that judges every window in isolation?"
            .to_string(),
        cohort: CohortReport {
            n_haplotypes: cohort.haplotypes.len(),
            n_modern: modern_count,
            n_windows,
            window_bp: cohort.config.window_bp,
            tract_mean_windows_requested: TRACT_MEAN_WINDOWS,
            planted_tracts: planted_tracts_in(&cohort, &all_windows).len(),
            planted_archaic_windows: cohort.truth.len(),
            archaic_fraction_of_segments: cohort.truth.len() as f64
                / (n_windows * modern_count) as f64,
            planted_by_source,
            planted_tract_lengths: planted_stats,
        },
        detector_params: dp,
        linkage_params: lp,
        morgan_scale: scale,
        coordinate_convention: scale.convention(),
        per_window: per_window_report,
        linkage: linkage_report,
        delta,
        admixture_time_truth: truth_time,
        held_out: HeldOut {
            n_train_windows: train.len(),
            n_test_windows: test.len(),
            tuned_params: tuned,
            per_window_test: pw_test,
            linkage_test: lk_test,
            linkage_wins_on_held_out: lk_test.f1 > pw_test.f1,
            sweep,
            note: "Linkage parameters were chosen by maximising F1 on the train windows only. \
                   The per-window detector was not tuned at all, so this comparison is if \
                   anything generous to the linkage rule."
                .to_string(),
        },
        unlinked_control: ControlReport {
            note: "The same linkage rule applied to a cohort with tract_mean_windows = 1.0, \
                   where introgressed windows are isolated by construction. A rule that helps \
                   here is not exploiting linkage."
                .to_string(),
            per_window: ctrl_pw,
            linkage: ctrl_lk,
            linkage_wins: ctrl_lk.f1 > ctrl_pw.f1,
        },
        findings,
        caveats: vec![
            "The genomes are simulated, not real. Nothing here is a claim about any living \
             person or population."
                .to_string(),
            "Windows are simulated as independent loci; linkage is imposed afterwards by \
             planting tracts as geometric runs of consecutive windows. There is no \
             recombination map in the simulation, so the windows-to-Morgans conversion is an \
             asserted convention, not a measurement."
                .to_string(),
            "The single-pulse tract-length clock (L ~ Exp(t), t = 1/E[L]) assumes one pulse \
             and no subsequent gene flow. This cohort carries six pulses from four sources, \
             so the inferred time is a blend across all of them, not any one event."
                .to_string(),
            "With 360 windows the whole simulated chromosome is under a Morgan, and a mean \
             tract of 5 windows is ~1 Mb — an order of magnitude longer than real Neanderthal \
             tracts. The implied admixture time is correspondingly far more recent than the \
             true Neanderthal or Denisovan pulses. That is a property of the simulator's \
             length budget, not a finding about hominin history."
                .to_string(),
            "Both detectors read the same depth grid, so the comparison isolates the decision \
             rule. It does not isolate the choice of thresholds, which interact."
                .to_string(),
            "The cohort is fully deterministic in the seed, but the HNSW index build is not \
             reproducible across processes, so the depth grid — and with it both detectors' \
             absolute scores — move by roughly +-0.003 F1 from run to run. Both detectors read \
             the same grid within a run, so the *difference* between them is far more stable \
             than either absolute number. Quote the delta, not the third decimal place."
                .to_string(),
        ],
        reproducibility: Reproducibility {
            cohort_deterministic: true,
            hnsw_deterministic_across_processes: false,
            note: "Rerunning this binary reproduces the cohort exactly (same seed, same planted \
                   tracts) but rebuilds the HNSW graph, whose construction is not deterministic \
                   across processes. Observed run-to-run spread on the headline F1 is a few \
                   thousandths for both detectors, with the linkage-vs-per-window delta holding \
                   its sign and rough magnitude."
                .to_string(),
        },
    };

    let path = format!("{out_dir}/linkage.json");
    let mut f = std::fs::File::create(&path)?;
    f.write_all(serde_json::to_string_pretty(&report)?.as_bytes())?;
    println!("\nWrote {path}");
    Ok(())
}
