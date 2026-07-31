//! # Linkage-aware introgression detection
//!
//! The per-window detector in [`crate::archaic`] judges every segment in
//! isolation: measure the coalescent depth from this window of this haplotype to
//! its nearest living relative, and call it archaic if that depth clears a
//! threshold. That throws away the single strongest structural fact about
//! introgression.
//!
//! ## The fact it throws away
//!
//! Ancestry does not arrive one locus at a time. A single interbreeding event
//! hands over whole chromosomes; recombination then chews them down over the
//! following generations. What survives into a present-day genome is a set of
//! **contiguous runs** — tracts — not a scatter of independent deep windows. So
//! the depth grid is not `n_windows x n_haplotypes` independent measurements. It
//! is a set of measurements with strong spatial autocorrelation *along the
//! chromosome*, and a detector that ignores that is leaving evidence unused.
//!
//! Concretely: a window whose depth sits just under the threshold is a coin
//! flip on its own evidence. The same window with four deep neighbours on either
//! side is almost certainly the middle of a real tract, because ordinary
//! variation has no mechanism that would make five adjacent windows deep at
//! once. That asymmetry is what this module exploits.
//!
//! ## What it does, and what it deliberately does not do
//!
//! [`call_linkage`] is **monotone**: it starts from the per-window call set and
//! can only *add* to it. A borderline window is rescued when the smoothed depth
//! of its neighbourhood — computed excluding the window itself, so the context
//! is independent evidence — is itself deep. Because it only adds, recall can
//! never fall below the per-window detector's; the open question, which
//! `linkage-study` answers empirically rather than by assertion, is whether the
//! windows it adds are mostly true (precision holds, F1 rises) or mostly noise.
//!
//! The optional [`LinkageParams::prune_isolated`] mode breaks that guarantee on
//! purpose: it drops called runs shorter than `min_run`, trading recall for
//! precision. It is off by default and reported separately.
//!
//! ## Tract length as a clock
//!
//! Tract length carries information the depth grid does not. Depth tells you
//! *when the donor lineage split* from ours — hundreds of thousands of years.
//! Tract length tells you *when the interbreeding happened* — a far more recent
//! event, and one no amount of per-window depth measurement can reach.
//!
//! After a single pulse `t` generations ago, the ancestry tract length `L` is
//! exponentially distributed with mean `1/t` Morgans (each generation of
//! recombination breaks a tract at rate 1 per Morgan per generation, so the
//! surviving length after `t` generations is `Exp(t)`). Inverting,
//! `t = 1 / E[L]`. See [`admixture_time_from_tracts`] — and read
//! [`MorganScale`] before believing any number it returns, because converting
//! *windows* to *Morgans* depends on the simulator's rescaled length scale.

use std::collections::{HashMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::archaic::{
    ArchaicSource, Cohort, DemographyConfig, DetectorParams, Role, Scores, Tract, GENERATION_YEARS,
};

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Knobs for the linkage-aware decision rule.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LinkageParams {
    /// Primary threshold, in ka. Identical in meaning to
    /// [`DetectorParams::tau_archaic_ka`]: any window at or above this is called
    /// on its own evidence, exactly as the per-window detector would.
    pub tau_archaic_ka: f64,
    /// Relaxed threshold, in ka, for windows that are *borderline*. A window
    /// between this and `tau_archaic_ka` is a candidate for rescue, but is only
    /// called if its neighbourhood backs it up. Must be `<= tau_archaic_ka`;
    /// setting it equal disables rescue entirely.
    pub tau_rescue_ka: f64,
    /// How deep, in ka, the smoothed neighbourhood must be for a rescue.
    pub tau_context_ka: f64,
    /// Neighbourhood radius, in windows, on each side.
    pub half_width: usize,
    /// Geometric weight decay per window of separation. A neighbour `d` windows
    /// away contributes with weight `decay^d`. `1.0` is a flat box filter; small
    /// values shrink the effective neighbourhood towards the immediate
    /// neighbours.
    pub decay: f64,
    /// If true, additionally *drop* called runs shorter than `min_run`. This
    /// breaks the monotonicity guarantee and can lower recall; off by default.
    pub prune_isolated: bool,
    /// Minimum run length kept when `prune_isolated` is set.
    pub min_run: usize,
}

impl Default for LinkageParams {
    fn default() -> Self {
        Self::from_detector(&DetectorParams::default())
    }
}

impl LinkageParams {
    /// Derive linkage parameters from a per-window detector, so the two share a
    /// primary threshold and the comparison is like-for-like.
    ///
    /// The rescue band is set to 72% of the primary threshold and the required
    /// neighbourhood depth to 55% of it. Both are well above the bulk of the
    /// ordinary-variation depth distribution (which sits near 100–250 ka for
    /// this demography) while sitting below the archaic mode, which is the band
    /// where a window's own evidence genuinely is ambiguous.
    pub fn from_detector(p: &DetectorParams) -> Self {
        Self {
            tau_archaic_ka: p.tau_archaic_ka,
            tau_rescue_ka: 0.72 * p.tau_archaic_ka,
            tau_context_ka: 0.55 * p.tau_archaic_ka,
            half_width: 3,
            decay: 0.6,
            prune_isolated: false,
            min_run: 2,
        }
    }
}

// ---------------------------------------------------------------------------
// The smoother
// ---------------------------------------------------------------------------

/// Weighted-geometric-mean smoother over neighbouring windows, **excluding the
/// window itself**.
///
/// `grid[i][j]` is the coalescent depth in ka of `windows[i]` in modern
/// haplotype `j`, as returned by [`crate::archaic::TraceEngine::depth_grid`].
/// The result has the same shape, and `out[i][j]` is the smoothed depth of the
/// *context around* `(i, j)`.
///
/// Two deliberate choices:
///
/// * **Smoothing happens in log space.** Coalescent waiting times are
///   exponential, so depth is log-normal-ish rather than normal; a single very
///   deep neighbour would otherwise drag an arithmetic mean up on its own.
/// * **The centre window is excluded.** The point of the context score is to be
///   evidence the window does not already have. Including the centre would let a
///   deep isolated window vouch for itself, and the rescue rule would collapse
///   into thresholding at `tau_rescue_ka`.
///
/// Adjacency is computed from the *window indices in `windows`*, not from
/// positions in the slice, so a non-contiguous window set (a train/test split,
/// say) does not silently glue unrelated windows together. Windows further than
/// `half_width` apart are not neighbours. A window with no neighbours at all
/// gets a context of `0.0`, which fails every threshold — no neighbours means no
/// support, not free support.
pub fn smooth_depth_grid(grid: &[Vec<f64>], windows: &[usize], p: &LinkageParams) -> Vec<Vec<f64>> {
    let n_win = grid.len();
    let n_hap = grid.first().map(|r| r.len()).unwrap_or(0);
    let mut out = vec![vec![0.0f64; n_hap]; n_win];
    if n_win == 0 || n_hap == 0 {
        return out;
    }

    // Precompute, for each row, which other rows are neighbours and with what
    // weight. Window sets are small and this is done once per grid.
    let mut neighbours: Vec<Vec<(usize, f64)>> = Vec::with_capacity(n_win);
    for i in 0..n_win {
        let mut ns = Vec::new();
        // Rows are in the order the caller asked for; scan outwards but bail as
        // soon as the window distance exceeds the half width in both directions.
        for i2 in 0..n_win {
            if i2 == i {
                continue;
            }
            let d = windows[i].abs_diff(windows[i2]);
            if d == 0 || d > p.half_width {
                continue;
            }
            ns.push((i2, p.decay.powi(d as i32)));
        }
        neighbours.push(ns);
    }

    for i in 0..n_win {
        for j in 0..n_hap {
            let mut wsum = 0.0;
            let mut acc = 0.0;
            for &(i2, w) in &neighbours[i] {
                // `max(1.0)` keeps the log defined; a depth under 1 ka is
                // indistinguishable from zero for our purposes anyway.
                acc += w * grid[i2][j].max(1.0).ln();
                wsum += w;
            }
            out[i][j] = if wsum > 0.0 { (acc / wsum).exp() } else { 0.0 };
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Decision rules
// ---------------------------------------------------------------------------

/// The per-window rule, for comparison: call a segment archaic iff its own depth
/// clears `tau`. This is exactly what [`crate::archaic::TraceEngine::detect`]
/// does on an engine with no ghost references, and the equivalence is asserted
/// by `linkage-study` and by unit test.
pub fn call_per_window(grid: &[Vec<f64>], tau_ka: f64) -> Vec<Vec<bool>> {
    grid.iter()
        .map(|row| row.iter().map(|&d| d >= tau_ka).collect())
        .collect()
}

/// The linkage-aware rule.
///
/// A window is called if **either**
///
/// 1. its own depth clears `tau_archaic_ka` — the per-window rule, unchanged; or
/// 2. its own depth clears the relaxed `tau_rescue_ka` **and** the smoothed
///    depth of its neighbourhood clears `tau_context_ka`.
///
/// Rule 1 alone is the baseline, so with `prune_isolated == false` the call set
/// is a superset of the baseline's and recall cannot decrease. With
/// `prune_isolated == true` a third step removes runs shorter than `min_run`,
/// which can and does remove true calls.
pub fn call_linkage(grid: &[Vec<f64>], windows: &[usize], p: &LinkageParams) -> Vec<Vec<bool>> {
    let context = smooth_depth_grid(grid, windows, p);
    let n_win = grid.len();
    let n_hap = grid.first().map(|r| r.len()).unwrap_or(0);

    let mut called = vec![vec![false; n_hap]; n_win];
    for i in 0..n_win {
        for j in 0..n_hap {
            let own = grid[i][j];
            let primary = own >= p.tau_archaic_ka;
            let rescued = own >= p.tau_rescue_ka && context[i][j] >= p.tau_context_ka;
            called[i][j] = primary || rescued;
        }
    }

    if p.prune_isolated && p.min_run > 1 {
        prune_short_runs(&mut called, windows, p.min_run);
    }
    called
}

/// Drop maximal runs of called windows shorter than `min_run`.
fn prune_short_runs(called: &mut [Vec<bool>], windows: &[usize], min_run: usize) {
    let n_win = called.len();
    let n_hap = called.first().map(|r| r.len()).unwrap_or(0);
    for j in 0..n_hap {
        let mut i = 0usize;
        while i < n_win {
            if !called[i][j] {
                i += 1;
                continue;
            }
            let mut end = i + 1;
            while end < n_win && called[end][j] && windows[end] == windows[end - 1] + 1 {
                end += 1;
            }
            if end - i < min_run {
                for k in i..end {
                    called[k][j] = false;
                }
            }
            i = end;
        }
    }
}

// ---------------------------------------------------------------------------
// Tract reconstruction
// ---------------------------------------------------------------------------

/// A maximal run of consecutive called windows in one haplotype.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CalledTract {
    /// Cohort haplotype index (not the column index in the grid).
    pub haplotype: usize,
    pub start_window: usize,
    pub length: usize,
}

impl CalledTract {
    /// One past the last window in the tract.
    pub fn end_window(&self) -> usize {
        self.start_window + self.length
    }
}

/// Reconstruct called tracts from a boolean call grid.
///
/// A tract is a **maximal** run of called windows that are consecutive both in
/// the call grid and in true window coordinates: a gap in `windows` breaks a run
/// even if both sides are called, because we have no evidence about what sits in
/// the gap. `modern` maps grid columns back to cohort haplotype indices.
pub fn reconstruct_tracts(
    called: &[Vec<bool>],
    windows: &[usize],
    modern: &[usize],
) -> Vec<CalledTract> {
    let n_win = called.len();
    let mut out = Vec::new();
    for (j, &h) in modern.iter().enumerate() {
        let mut i = 0usize;
        while i < n_win {
            if !called[i][j] {
                i += 1;
                continue;
            }
            let start = i;
            let mut end = i + 1;
            while end < n_win && called[end][j] && windows[end] == windows[end - 1] + 1 {
                end += 1;
            }
            out.push(CalledTract {
                haplotype: h,
                start_window: windows[start],
                length: end - start,
            });
            i = end;
        }
    }
    out.sort_by_key(|t| (t.haplotype, t.start_window));
    out
}

// ---------------------------------------------------------------------------
// Tract-length distributions
// ---------------------------------------------------------------------------

/// Summary of a tract-length distribution, in windows.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TractLengthStats {
    pub n: usize,
    pub mean_windows: f64,
    pub sd_windows: f64,
    pub median_windows: f64,
    pub p90_windows: f64,
    pub max_windows: usize,
    /// `histogram[i]` is the number of tracts of length `i + 1` windows.
    pub histogram: Vec<usize>,
    /// Total windows covered, i.e. the sum of all lengths.
    pub total_windows: usize,
}

impl TractLengthStats {
    pub fn from_lengths(lengths: &[usize]) -> Self {
        let n = lengths.len();
        if n == 0 {
            return Self {
                n: 0,
                mean_windows: 0.0,
                sd_windows: 0.0,
                median_windows: 0.0,
                p90_windows: 0.0,
                max_windows: 0,
                histogram: Vec::new(),
                total_windows: 0,
            };
        }
        let mut sorted: Vec<usize> = lengths.to_vec();
        sorted.sort_unstable();
        let total: usize = sorted.iter().sum();
        let mean = total as f64 / n as f64;
        let var = sorted
            .iter()
            .map(|&l| (l as f64 - mean).powi(2))
            .sum::<f64>()
            / n as f64;
        let q = |f: f64| -> f64 { sorted[(((n - 1) as f64) * f).round() as usize] as f64 };
        let max = *sorted.last().unwrap();
        let mut histogram = vec![0usize; max];
        for &l in &sorted {
            if l >= 1 {
                histogram[l - 1] += 1;
            }
        }
        Self {
            n,
            mean_windows: mean,
            sd_windows: var.sqrt(),
            median_windows: q(0.5),
            p90_windows: q(0.9),
            max_windows: max,
            histogram,
            total_windows: total,
        }
    }
}

/// Planted versus recovered tract-length distributions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TractComparison {
    pub planted: TractLengthStats,
    pub recovered: TractLengthStats,
    /// `recovered.mean / planted.mean`. Below 1 means the detector is
    /// fragmenting real tracts; above 1 means it is fusing them (or inventing
    /// long false ones).
    pub mean_length_ratio: f64,
    /// Recovered tracts per planted tract. Above 1 is fragmentation.
    pub fragmentation: f64,
    /// Two-sample Kolmogorov–Smirnov statistic between the two length
    /// distributions: the largest gap between their empirical CDFs, in [0, 1].
    /// 0 is identical distributions.
    pub ks_statistic: f64,
}

/// Compare a recovered tract-length distribution against the planted truth.
///
/// Both are taken as plain length lists in windows; the KS statistic compares
/// their shapes without assuming either is exponential.
pub fn compare_tract_lengths(planted: &[usize], recovered: &[usize]) -> TractComparison {
    let p = TractLengthStats::from_lengths(planted);
    let r = TractLengthStats::from_lengths(recovered);

    let ks = {
        let max_len = p.max_windows.max(r.max_windows);
        let mut d: f64 = 0.0;
        if p.n > 0 && r.n > 0 {
            let (mut cp, mut cr) = (0usize, 0usize);
            for l in 1..=max_len {
                cp += planted.iter().filter(|&&x| x == l).count();
                cr += recovered.iter().filter(|&&x| x == l).count();
                let gap = (cp as f64 / p.n as f64 - cr as f64 / r.n as f64).abs();
                if gap > d {
                    d = gap;
                }
            }
        }
        d
    };

    TractComparison {
        mean_length_ratio: if p.mean_windows > 0.0 {
            r.mean_windows / p.mean_windows
        } else {
            0.0
        },
        fragmentation: if p.n > 0 {
            r.n as f64 / p.n as f64
        } else {
            0.0
        },
        ks_statistic: ks,
        planted: p,
        recovered: r,
    }
}

// ---------------------------------------------------------------------------
// Windows -> Morgans, and the admixture clock
// ---------------------------------------------------------------------------

/// Human autosomal per-base per-generation mutation rate, the real one.
pub const HUMAN_MU_PER_BP: f64 = 1.25e-8;

/// Genome-average human recombination rate, in Morgans per base pair
/// (≈ 1 cM/Mb — the standard autosomal average; real rates vary by an order of
/// magnitude between hotspots and cold deserts).
pub const HUMAN_RECOMBINATION_M_PER_BP: f64 = 1.0e-8;

/// The coordinate convention that turns *windows* into *Morgans*.
///
/// **Read this before using any number derived from it.** The simulator does not
/// have a recombination map. Its windows are simulated as independent loci and
/// linkage is imposed afterwards, by planting tracts as geometric runs of
/// consecutive windows. So there is no physical distance in the simulation to
/// read off; a genetic length scale has to be *asserted*, and every
/// tract-length-based time estimate inherits whatever we assert.
///
/// The assertion made here, in three steps:
///
/// 1. **A simulated window is not its own length in base pairs.** The simulator
///    rescales the mutation rate upward (see [`DemographyConfig::mu`]) so a
///    short window carries the segregating-site count of a much longer real one.
///    The rescaling factor is `mu_simulated / mu_real`, and the *effective*
///    length of a window is `window_bp * mu_simulated / mu_real`. With the
///    defaults (`window_bp = 2000`, `mu = 1.36e-6`) that is ≈ 218 kb.
/// 2. **Effective base pairs convert to Morgans at the genome-average rate**,
///    `1e-8 M/bp`. So one default window ≈ `2.18e-3` Morgans.
/// 3. **Adjacent windows are assumed adjacent on one chromosome**, with no gap.
///    A tract of `L` windows is therefore `L * morgans_per_window` Morgans.
///
/// Step 1 is the load-bearing one and it is where the honesty lives: if you
/// changed `mu` without changing `window_bp`, every window would silently become
/// a different genetic length, and every admixture-time estimate would move even
/// though nothing about the genealogy changed. Step 2 is a genome-wide average
/// that no individual locus obeys. Step 3 is true by construction here and false
/// in any real genome, where windows are separated by unsequenced gaps.
///
/// The practical consequence for this study: with 420 windows the whole
/// simulated chromosome is only ~0.9 Morgans, and a mean tract of 4–6 windows is
/// ~1 Mb — an order of magnitude longer than real Neanderthal tracts (~50–80 kb).
/// So the implied admixture time comes out far *more recent* than the true
/// Neanderthal or Denisovan pulses. That is a statement about the simulator's
/// length budget, not a finding about hominin history.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct MorganScale {
    pub window_bp: usize,
    /// The simulator's rescaled mutation rate.
    pub mu_simulated: f64,
    /// The real human rate the simulation is rescaled against.
    pub mu_real: f64,
    /// `mu_simulated / mu_real`.
    pub rescale_factor: f64,
    /// `window_bp * rescale_factor` — how much real sequence a window stands in
    /// for, in base pairs.
    pub effective_window_bp: f64,
    pub recombination_m_per_bp: f64,
    /// The number this whole struct exists to produce.
    pub morgans_per_window: f64,
}

impl MorganScale {
    /// Derive the scale from a demography config.
    pub fn from_config(c: &DemographyConfig) -> Self {
        Self::new(
            c.window_bp,
            c.mu,
            HUMAN_MU_PER_BP,
            HUMAN_RECOMBINATION_M_PER_BP,
        )
    }

    pub fn new(window_bp: usize, mu_simulated: f64, mu_real: f64, recomb: f64) -> Self {
        let rescale = mu_simulated / mu_real;
        let eff = window_bp as f64 * rescale;
        Self {
            window_bp,
            mu_simulated,
            mu_real,
            rescale_factor: rescale,
            effective_window_bp: eff,
            recombination_m_per_bp: recomb,
            morgans_per_window: eff * recomb,
        }
    }

    /// Human-readable statement of the convention, for embedding in reports so a
    /// reader never sees the derived time without the assumption behind it.
    pub fn convention(&self) -> String {
        format!(
            "1 window = {} bp simulated at mu = {:.3e}, which is {:.0}x the real human rate \
             {:.2e}, so it stands in for {:.0} kb of real sequence; at {:.1} cM/Mb that is \
             {:.3e} Morgans per window. Adjacent windows are assumed contiguous on one \
             chromosome. Every tract-length-derived time below inherits these three \
             assumptions.",
            self.window_bp,
            self.mu_simulated,
            self.rescale_factor,
            self.mu_real,
            self.effective_window_bp / 1000.0,
            self.recombination_m_per_bp * 1.0e8,
            self.morgans_per_window
        )
    }
}

/// An admixture time inferred from mean tract length.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdmixtureTimeEstimate {
    /// What this estimate was computed from, e.g. `"planted"` / `"recovered"`.
    pub source: String,
    pub n_tracts: usize,
    pub mean_tract_windows: f64,
    pub mean_tract_morgans: f64,
    /// `t = 1 / E[L]`, in generations before present.
    pub generations: f64,
    pub years: f64,
    pub ka: f64,
    /// 95% interval on `generations`, from the sampling error of an exponential
    /// mean (`SE(L) = L / sqrt(n)`), propagated through `t = 1/L`. It covers
    /// sampling noise only — it says nothing about whether the coordinate
    /// convention in [`MorganScale`] is right.
    pub generations_ci95: (f64, f64),
    pub scale: MorganScale,
    pub convention: String,
}

/// Estimate time since admixture from a set of tract lengths, in windows.
///
/// Uses the standard single-pulse result: `t` generations after a pulse, an
/// ancestry tract has length `L ~ Exp(t)` in Morgans, i.e. `E[L] = 1/t`, so
/// `t = 1 / E[L]`. This assumes a *single* pulse, no subsequent gene flow, and
/// tract boundaries observed without error — none of which hold exactly here:
/// this cohort carries six pulses from four sources, and the recovered tracts
/// are measured with error.
///
/// Returns `None` if there are no tracts, or if the mean length is zero.
pub fn admixture_time_from_tracts(
    lengths: &[usize],
    scale: &MorganScale,
    source: impl Into<String>,
) -> Option<AdmixtureTimeEstimate> {
    if lengths.is_empty() {
        return None;
    }
    let n = lengths.len();
    let mean_windows = lengths.iter().sum::<usize>() as f64 / n as f64;
    let mean_morgans = mean_windows * scale.morgans_per_window;
    if mean_morgans <= 0.0 {
        return None;
    }
    let t = 1.0 / mean_morgans;

    // SE of an exponential mean is mean/sqrt(n); push the +-1.96 SE bounds on
    // the mean through the reciprocal, which flips their order.
    let rel = 1.96 / (n as f64).sqrt();
    let lo_mean = mean_morgans * (1.0 + rel);
    let hi_mean = mean_morgans * (1.0 - rel).max(1e-9);

    Some(AdmixtureTimeEstimate {
        source: source.into(),
        n_tracts: n,
        mean_tract_windows: mean_windows,
        mean_tract_morgans: mean_morgans,
        generations: t,
        years: t * GENERATION_YEARS,
        ka: t * GENERATION_YEARS / 1000.0,
        generations_ci95: (1.0 / lo_mean, 1.0 / hi_mean),
        scale: *scale,
        convention: scale.convention(),
    })
}

// ---------------------------------------------------------------------------
// Scoring
// ---------------------------------------------------------------------------

/// Score a boolean call grid against planted truth, restricted to `windows` and
/// to modern haplotypes — the same accounting
/// [`crate::archaic::TraceEngine::score`] uses, so the two detectors' numbers
/// are directly comparable.
pub fn score_calls(
    called: &[Vec<bool>],
    windows: &[usize],
    modern: &[usize],
    truth: &HashMap<(usize, usize), ArchaicSource>,
) -> Scores {
    let win: HashSet<usize> = windows.iter().copied().collect();
    let modern_set: HashSet<usize> = modern.iter().copied().collect();

    let mut called_set: HashSet<(usize, usize)> = HashSet::new();
    for (i, &w) in windows.iter().enumerate() {
        for (j, &h) in modern.iter().enumerate() {
            if called[i][j] {
                called_set.insert((w, h));
            }
        }
    }
    let truth_set: HashSet<(usize, usize)> = truth
        .keys()
        .copied()
        .filter(|(w, h)| win.contains(w) && modern_set.contains(h))
        .collect();

    let tp = called_set.intersection(&truth_set).count();
    Scores::compute(tp, called_set.len() - tp, truth_set.len() - tp)
}

/// Planted tract lengths, restricted to tracts that lie wholly inside `windows`
/// and belong to a modern haplotype.
///
/// Tracts that run off the edge of the analysed window set are dropped rather
/// than truncated: a truncated tract has a censored length, and averaging
/// censored lengths in with complete ones biases the mean down, which would bias
/// the admixture-time estimate *older* for no reason but bookkeeping.
pub fn planted_tract_lengths(cohort: &Cohort, windows: &[usize]) -> Vec<usize> {
    let win: HashSet<usize> = windows.iter().copied().collect();
    cohort
        .tracts
        .iter()
        .filter(|t| cohort.haplotypes[t.haplotype].role == Role::Modern)
        .filter(|t| (t.start_window..t.start_window + t.length).all(|w| win.contains(&w)))
        .map(|t| t.length)
        .collect()
}

/// Planted tracts as [`CalledTract`]s, for like-for-like comparison against
/// reconstructed ones.
pub fn planted_tracts_in(cohort: &Cohort, windows: &[usize]) -> Vec<CalledTract> {
    let win: HashSet<usize> = windows.iter().copied().collect();
    let mut out: Vec<CalledTract> = cohort
        .tracts
        .iter()
        .filter(|t: &&Tract| cohort.haplotypes[t.haplotype].role == Role::Modern)
        .filter(|t| (t.start_window..t.start_window + t.length).all(|w| win.contains(&w)))
        .map(|t| CalledTract {
            haplotype: t.haplotype,
            start_window: t.start_window,
            length: t.length,
        })
        .collect();
    out.sort_by_key(|t| (t.haplotype, t.start_window));
    out
}

// ---------------------------------------------------------------------------
// One-shot driver
// ---------------------------------------------------------------------------

/// Everything one linkage-aware pass produced.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinkageRun {
    pub params: LinkageParams,
    pub scores: Scores,
    pub tracts: Vec<CalledTract>,
    pub tract_lengths: TractLengthStats,
    /// Windows called by the linkage rule that the per-window rule missed.
    pub rescued_windows: usize,
    /// Of those, how many were true.
    pub rescued_true: usize,
}

/// Run the linkage rule over a depth grid and score it.
pub fn run_linkage(
    grid: &[Vec<f64>],
    windows: &[usize],
    modern: &[usize],
    truth: &HashMap<(usize, usize), ArchaicSource>,
    p: &LinkageParams,
) -> LinkageRun {
    let base = call_per_window(grid, p.tau_archaic_ka);
    let called = call_linkage(grid, windows, p);

    let mut rescued = 0usize;
    let mut rescued_true = 0usize;
    for (i, &w) in windows.iter().enumerate() {
        for (j, &h) in modern.iter().enumerate() {
            if called[i][j] && !base[i][j] {
                rescued += 1;
                if truth.contains_key(&(w, h)) {
                    rescued_true += 1;
                }
            }
        }
    }

    let tracts = reconstruct_tracts(&called, windows, modern);
    let lengths: Vec<usize> = tracts.iter().map(|t| t.length).collect();

    LinkageRun {
        params: *p,
        scores: score_calls(&called, windows, modern, truth),
        tract_lengths: TractLengthStats::from_lengths(&lengths),
        tracts,
        rescued_windows: rescued,
        rescued_true,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archaic::{simulate_cohort, Pulse};

    fn tracted_config(mean: f64) -> DemographyConfig {
        let mut c = DemographyConfig {
            window_bp: 800,
            n_windows: 60,
            ..Default::default()
        };
        for p in c.pulses.iter_mut() {
            p.tract_mean_windows = mean;
        }
        c
    }

    /// Planted tracts must be contiguous runs in `truth`, not a scatter.
    #[test]
    fn tract_planting_produces_contiguous_runs() {
        let cohort = simulate_cohort(tracted_config(5.0));
        assert!(!cohort.tracts.is_empty(), "no tracts were planted");

        for t in &cohort.tracts {
            assert!(t.length >= 1);
            for w in t.start_window..t.start_window + t.length {
                assert_eq!(
                    cohort.truth.get(&(w, t.haplotype)),
                    Some(&t.source),
                    "tract {:?} claims window {w} but truth disagrees",
                    t
                );
            }
        }

        // Every truth cell is covered by exactly one tract, and vice versa.
        let mut covered: HashSet<(usize, usize)> = HashSet::new();
        for t in &cohort.tracts {
            for w in t.start_window..t.start_window + t.length {
                assert!(
                    covered.insert((w, t.haplotype)),
                    "tracts overlap at ({w}, {})",
                    t.haplotype
                );
            }
        }
        let truth_cells: HashSet<(usize, usize)> = cohort.truth.keys().copied().collect();
        assert_eq!(
            covered, truth_cells,
            "tracts and truth cover different cells"
        );

        // The whole point: with mean 5 the runs are actually long.
        let mean_len = cohort.tracts.iter().map(|t| t.length).sum::<usize>() as f64
            / cohort.tracts.len() as f64;
        assert!(
            mean_len > 2.0,
            "mean planted tract length {mean_len:.2} is not tract-like"
        );
    }

    /// The published-numbers guard: `tract_mean_windows = 1.0` must reproduce the
    /// original single-window behaviour exactly.
    #[test]
    fn tract_mean_one_reproduces_single_window_behaviour() {
        let cohort = simulate_cohort(DemographyConfig {
            window_bp: 700,
            n_windows: 40,
            ..Default::default()
        });

        // The default config must still be all-ones, or every published number
        // silently moves.
        for p in &DemographyConfig::default().pulses {
            assert_eq!(
                p.tract_mean_windows, 1.0,
                "default pulse {:?} is no longer single-window",
                p.source
            );
        }

        assert!(!cohort.tracts.is_empty());
        for t in &cohort.tracts {
            assert_eq!(t.length, 1, "tract {:?} is longer than one window", t);
        }
        // `tracts` is then exactly `truth` re-expressed as a list.
        assert_eq!(cohort.tracts.len(), cohort.truth.len());
        for t in &cohort.tracts {
            assert_eq!(
                cohort.truth.get(&(t.start_window, t.haplotype)),
                Some(&t.source)
            );
        }
    }

    /// Serde back-compat: an old config with no `tract_mean_windows` field must
    /// still deserialise, and must default to the single-window behaviour.
    #[test]
    fn missing_tract_field_defaults_to_one() {
        let json = r#"{"source":"Neanderthal","recipient_group":"NONAFR","fraction":0.019}"#;
        let p: Pulse = serde_json::from_str(json).expect("legacy pulse should deserialise");
        assert_eq!(p.tract_mean_windows, 1.0);
    }

    /// The monotonicity guarantee, checked on a synthetic grid: the linkage rule
    /// calls a superset of the per-window rule.
    #[test]
    fn linkage_never_lowers_recall_on_a_synthetic_grid() {
        // 12 windows, 3 haplotypes. Haplotype 0 carries a deep run with one
        // borderline dip in the middle; haplotype 1 is ordinary; haplotype 2 has
        // one isolated borderline window with no support.
        let deep = 1200.0;
        let borderline = 500.0;
        let ordinary = 150.0;
        let mut grid = vec![vec![ordinary; 3]; 12];
        for w in 3..9 {
            grid[w][0] = deep;
        }
        grid[6][0] = borderline;
        grid[5][2] = borderline;

        let windows: Vec<usize> = (0..12).collect();
        let p = LinkageParams::default();
        let base = call_per_window(&grid, p.tau_archaic_ka);
        let linked = call_linkage(&grid, &windows, &p);

        for w in 0..12 {
            for h in 0..3 {
                assert!(
                    linked[w][h] || !base[w][h],
                    "linkage dropped a per-window call at ({w}, {h})"
                );
            }
        }
        // The dip inside the run gets rescued...
        assert!(
            !base[6][0] && linked[6][0],
            "the supported dip was not rescued"
        );
        // ...but the unsupported isolated window does not.
        assert!(
            !linked[5][2],
            "an isolated borderline window was rescued with no neighbourhood support"
        );
    }

    /// The same guarantee on a real tracted cohort, using real depths rather
    /// than a hand-built grid.
    #[test]
    fn linkage_smoother_never_lowers_recall_on_a_tracted_cohort() {
        use crate::archaic::TraceEngine;

        let cohort = simulate_cohort(tracted_config(5.0));
        let storage = std::env::temp_dir().join(format!(
            "trace_rv_linkage_test_{}_{}",
            std::process::id(),
            "recall"
        ));
        std::fs::remove_dir_all(&storage).ok();
        let mut engine = TraceEngine::new(&cohort, storage.to_string_lossy().to_string());

        let dp = DetectorParams::default();
        let windows: Vec<usize> = (0..cohort.sequences.len()).collect();
        let grid = engine.depth_grid(&dp, &windows).expect("depth grid");
        let modern = engine.modern_haplotypes();
        std::fs::remove_dir_all(&storage).ok();

        let lp = LinkageParams::from_detector(&dp);
        let base = call_per_window(&grid, dp.tau_archaic_ka);
        let base_scores = score_calls(&base, &windows, &modern, &cohort.truth);
        let run = run_linkage(&grid, &windows, &modern, &cohort.truth, &lp);

        assert!(
            run.scores.recall >= base_scores.recall,
            "linkage lowered recall: {:.4} vs per-window {:.4}",
            run.scores.recall,
            base_scores.recall
        );
        assert!(
            run.scores.true_positives >= base_scores.true_positives,
            "linkage lost true positives: {} vs {}",
            run.scores.true_positives,
            base_scores.true_positives
        );
    }

    /// Thresholding the depth grid must reproduce `detect()` exactly, or the two
    /// detectors are not being compared on the same evidence.
    #[test]
    fn depth_grid_agrees_with_detect() {
        use crate::archaic::TraceEngine;

        let cohort = simulate_cohort(DemographyConfig {
            window_bp: 700,
            n_windows: 24,
            ..Default::default()
        });
        let storage = std::env::temp_dir().join(format!(
            "trace_rv_linkage_test_{}_{}",
            std::process::id(),
            "grid"
        ));
        std::fs::remove_dir_all(&storage).ok();
        let mut engine = TraceEngine::new(&cohort, storage.to_string_lossy().to_string());

        let dp = DetectorParams::default();
        let windows: Vec<usize> = (0..cohort.sequences.len()).collect();
        let run = engine.detect(&dp, &windows, 0).expect("detect");
        let grid = engine.depth_grid(&dp, &windows).expect("grid");
        let modern = engine.modern_haplotypes();
        std::fs::remove_dir_all(&storage).ok();

        let from_detect: HashSet<(usize, usize)> =
            run.calls.iter().map(|c| (c.window, c.haplotype)).collect();
        let mut from_grid: HashSet<(usize, usize)> = HashSet::new();
        for (i, &w) in windows.iter().enumerate() {
            for (j, &h) in modern.iter().enumerate() {
                if grid[i][j] >= dp.tau_archaic_ka {
                    from_grid.insert((w, h));
                }
            }
        }
        assert_eq!(from_detect, from_grid, "depth_grid and detect disagree");

        // And the depths themselves must match, not just the calls.
        for c in &run.calls {
            let i = windows.iter().position(|&w| w == c.window).unwrap();
            let j = modern.iter().position(|&h| h == c.haplotype).unwrap();
            assert!(
                (grid[i][j] - c.depth_to_modern_ka).abs() < 1e-9,
                "depth mismatch at ({}, {}): {} vs {}",
                c.window,
                c.haplotype,
                grid[i][j],
                c.depth_to_modern_ka
            );
        }
    }

    #[test]
    fn reconstruct_tracts_finds_maximal_runs() {
        //          w0     w1     w2     w3     w4
        let called = vec![
            vec![true, false],
            vec![true, false],
            vec![false, true],
            vec![true, true],
            vec![true, false],
        ];
        let windows = vec![0usize, 1, 2, 3, 4];
        let modern = vec![10usize, 20];
        let t = reconstruct_tracts(&called, &windows, &modern);
        assert_eq!(
            t,
            vec![
                CalledTract {
                    haplotype: 10,
                    start_window: 0,
                    length: 2
                },
                CalledTract {
                    haplotype: 10,
                    start_window: 3,
                    length: 2
                },
                CalledTract {
                    haplotype: 20,
                    start_window: 2,
                    length: 2
                },
            ]
        );
    }

    /// A gap in the analysed window set must break a run: we have no evidence
    /// about the windows we did not look at.
    #[test]
    fn a_window_gap_breaks_a_run() {
        let called = vec![vec![true], vec![true], vec![true]];
        let windows = vec![0usize, 1, 7];
        let t = reconstruct_tracts(&called, &windows, &[0]);
        assert_eq!(t.len(), 2);
        assert_eq!(t[0].length, 2);
        assert_eq!(t[1].start_window, 7);
        assert_eq!(t[1].length, 1);
    }

    #[test]
    fn morgan_scale_matches_the_documented_convention() {
        let s = MorganScale::from_config(&DemographyConfig::default());
        // 1.36e-6 / 1.25e-8 = 108.8x; 2000 bp * 108.8 = 217.6 kb.
        assert!((s.rescale_factor - 108.8).abs() < 1e-6);
        assert!((s.effective_window_bp - 217_600.0).abs() < 1e-3);
        // 217600 bp * 1e-8 M/bp = 2.176e-3 M.
        assert!((s.morgans_per_window - 2.176e-3).abs() < 1e-9);
    }

    /// The clock inverts: feed it tracts of a known mean and it returns the
    /// generation count that would have produced them.
    #[test]
    fn admixture_clock_inverts_exponential_mean() {
        let scale = MorganScale::from_config(&DemographyConfig::default());
        // 100 tracts of exactly 5 windows -> mean 5 * 2.176e-3 = 1.088e-2 M.
        let lengths = vec![5usize; 100];
        let est = admixture_time_from_tracts(&lengths, &scale, "test").unwrap();
        assert!((est.mean_tract_morgans - 1.088e-2).abs() < 1e-9);
        assert!((est.generations - 1.0 / 1.088e-2).abs() < 1e-6);
        // Shorter tracts mean an older pulse.
        let older = admixture_time_from_tracts(&vec![2usize; 100], &scale, "test").unwrap();
        assert!(older.generations > est.generations);
        // The interval brackets the point estimate.
        assert!(est.generations_ci95.0 < est.generations);
        assert!(est.generations_ci95.1 > est.generations);
        assert!(admixture_time_from_tracts(&[], &scale, "test").is_none());
    }

    #[test]
    fn tract_comparison_detects_fragmentation() {
        let planted = vec![6usize; 20];
        let fragmented: Vec<usize> = vec![2usize; 60];
        let c = compare_tract_lengths(&planted, &fragmented);
        assert!((c.mean_length_ratio - 2.0 / 6.0).abs() < 1e-9);
        assert!((c.fragmentation - 3.0).abs() < 1e-9);
        assert!(
            c.ks_statistic > 0.9,
            "KS should be near 1 for disjoint supports"
        );

        let identical = compare_tract_lengths(&planted, &planted);
        assert!(identical.ks_statistic < 1e-12);
        assert!((identical.mean_length_ratio - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pruning_can_lower_recall_and_is_off_by_default() {
        assert!(!LinkageParams::default().prune_isolated);

        let mut grid = vec![vec![150.0; 1]; 6];
        grid[2][0] = 1200.0; // an isolated deep window
        let windows: Vec<usize> = (0..6).collect();
        let mut p = LinkageParams::default();
        assert!(call_linkage(&grid, &windows, &p)[2][0]);
        p.prune_isolated = true;
        p.min_run = 2;
        assert!(
            !call_linkage(&grid, &windows, &p)[2][0],
            "pruning should drop a singleton run"
        );
    }
}
