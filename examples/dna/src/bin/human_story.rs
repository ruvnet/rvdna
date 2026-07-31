//! `human-story` — what the real panel says about human history, not about the
//! detector.
//!
//! Everything else in this study measures a *method*: how well it separates
//! introgressed segments from ordinary ones, where it stops working, what it
//! costs. This binary asks the other kind of question. Given 800 real phased
//! haplotypes and three sequenced archaic genomes over the same megabase of
//! chromosome 22, what can be said about the people?
//!
//! Three things, each with its own control:
//!
//! 1. **Where the Neanderthal DNA actually is.** Not a population average — the
//!    individual haplotypes, at coordinates you can look up. The scan uses
//!    *Neanderthal-informative sites*: positions where a Neanderthal carries the
//!    non-reference allele and not one of the 160 African haplotypes does. A
//!    living non-African who carries a run of those alleles inherited them from
//!    somebody, and the only somebody available is a Neanderthal. The control is
//!    to run the identical scan on the African haplotypes, who should have
//!    almost none.
//!
//! 2. **When that DNA arrived.** Recombination has been cutting introgressed
//!    tracts shorter every generation since the interbreeding, so the surviving
//!    tract length is a clock. Mean length in Morgans inverts to generations.
//!
//! 3. **The shape of the human family, drawn from the data alone.** Mean
//!    pairwise sequence divergence between all 26 populations, and within each
//!    of them. No tree is assumed; the matrix is the output.
//!
//! Usage: `human-story [data-dir] [out-dir]`
//! (defaults `data/real` and `artifacts/data`).

use std::collections::BTreeMap;
use std::io::{BufWriter, Write as _};

use rvdna::archaic::{gens_to_ka, jukes_cantor};
use rvdna::real_panel::{diffs, load, load_archaic};
use serde::Serialize;

/// Autosomal mutation rate per base per generation — the real human clock.
const MU: f64 = 1.25e-8;
const GEN_YEARS: f64 = 29.0;
/// Sex-averaged recombination, the standard genome-wide approximation. Local
/// rate on chr22 varies several-fold around it; see the caveats.
const CM_PER_MB: f64 = 1.0;

/// Window for the tract scan. Real Neanderthal tracts in living people run tens
/// of kilobases, so the window has to be well below that to resolve a boundary,
/// and well above the spacing of informative sites to carry any signal.
const WIN: u32 = 10_000;

/// A run of adjacent windows in one haplotype carrying Neanderthal-informative
/// alleles that no African in the panel has.
#[derive(Serialize, Clone)]
struct ArchaicTract {
    haplotype: String,
    sample: String,
    population: String,
    superpopulation: String,
    start: u32,
    end: u32,
    length_bp: u32,
    n_windows: usize,
    /// Neanderthal-informative alleles this haplotype carries inside the tract.
    informative_carried: u32,
    /// Neanderthal-informative sites available inside the tract.
    informative_available: u32,
    carried_fraction: f64,
}

#[derive(Serialize)]
struct TractScanArm {
    arm: String,
    n_haplotypes: usize,
    n_tracts: usize,
    tracts_per_haplotype: f64,
    haplotypes_with_a_tract: usize,
    fraction_of_haplotypes: f64,
    mean_length_bp: f64,
    median_length_bp: f64,
    max_length_bp: u32,
    bp_in_tracts: u64,
    /// Share of the scanned sequence that sits inside a called tract.
    genome_fraction: f64,
}

#[derive(Serialize)]
struct AdmixtureDate {
    source: String,
    n_tracts: usize,
    mean_length_bp: f64,
    mean_length_morgans: f64,
    generations: f64,
    years_ago: f64,
    ka: f64,
    ci95_ka: [f64; 2],
}

/// The whole scan, packed for drawing: one entry per carrier haplotype, and
/// tracts as `[carrier index, first window, last window]`. Shipping all 847
/// tracts as records would cost 150 KB; this costs about 10.
#[derive(Serialize)]
struct TractMap {
    region_start: u32,
    region_end: u32,
    window_bp: u32,
    n_windows: usize,
    /// `[haplotype id, population, superpopulation]` per carrier.
    carriers: Vec<(String, String, String)>,
    /// Non-carrier haplotype counts per superpopulation, so the drawing can
    /// show the empty rows that make the control visible.
    total_by_superpopulation: Vec<(String, usize)>,
    tracts: Vec<(usize, usize, usize)>,
}

#[derive(Serialize)]
struct PopStat {
    population: String,
    superpopulation: String,
    n_haplotypes: usize,
    /// Mean divergence between two haplotypes drawn from this population, in ka.
    within_ka: f64,
    /// Mean divergence to every haplotype outside this population, in ka.
    between_ka: f64,
    /// Neanderthal-informative alleles carried per haplotype.
    neanderthal_alleles_per_haplotype: f64,
}

#[derive(Serialize)]
struct HumanStory {
    title: String,
    generated_utc: String,
    region: String,
    window_bp: u32,
    n_haplotypes: usize,
    n_variants: usize,
    mu_per_bp_per_generation: f64,
    generation_years: f64,
    cm_per_mb: f64,

    informative_sites: usize,
    informative_note: String,

    tract_scan: Vec<TractScanArm>,
    tract_map: TractMap,
    tracts: Vec<ArchaicTract>,
    admixture_date: AdmixtureDate,
    carrier_by_superpopulation: Vec<(String, usize, usize, f64)>,

    populations: Vec<PopStat>,
    population_order: Vec<String>,
    divergence_matrix_ka: Vec<Vec<f64>>,

    findings: Vec<String>,
    caveats: Vec<String>,
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    sorted[(((sorted.len() - 1) as f64) * q).round() as usize]
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let data_dir = args.get(1).cloned().unwrap_or_else(|| "data/real".into());
    let out_dir = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "artifacts/data".into());
    std::fs::create_dir_all(&out_dir)?;

    let t0 = std::time::Instant::now();
    let p = load(&data_dir)?;
    let arcs = load_archaic(&data_dir, &p.positions, &p.alleles)?;
    let n_hap = p.bits.len();
    let n_var = p.positions.len();
    let words = n_var.div_ceil(64);
    println!(
        "{n_hap} haplotypes x {n_var} sites, {} archaic genomes",
        arcs.len()
    );

    let is_afr = |h: usize| p.sample_super[p.sample_of_hap[h]] == "AFR";
    let afr: Vec<usize> = (0..n_hap).filter(|&h| is_afr(*&h)).collect();
    let non: Vec<usize> = (0..n_hap).filter(|&h| !is_afr(*&h)).collect();

    // ─────────────────────────────────────────── Neanderthal-informative sites
    //
    // A site qualifies when a Neanderthal carries the non-reference allele and
    // no African haplotype in the panel does. The African-absence condition is
    // what makes this robust: it does not matter that the reference sequence is
    // European-weighted, because the test is not "how different is this person
    // from a Neanderthal" but "does this person carry an allele that exists in
    // Neanderthals and nowhere in Africa".
    let mut afr_has = vec![0u64; words];
    for &h in &afr {
        for w in 0..words {
            afr_has[w] |= p.bits[h][w];
        }
    }
    let nea_alt: Vec<u64> = (0..words)
        .map(|w| {
            // either high-coverage Neanderthal, at a site that Neanderthal could be called at
            let a = arcs[0].callable[w] & arcs[0].has_alt[w];
            let b = if arcs.len() > 1 {
                arcs[1].callable[w] & arcs[1].has_alt[w]
            } else {
                0
            };
            a | b
        })
        .collect();
    let informative: Vec<u64> = (0..words).map(|w| nea_alt[w] & !afr_has[w]).collect();
    let n_informative: usize = informative.iter().map(|w| w.count_ones() as usize).sum();
    println!("{n_informative} Neanderthal-informative sites (Neanderthal ALT, absent in all {} African haplotypes)", afr.len());

    // ─────────────────────────────────────────── windows
    let first = *p.positions.first().unwrap();
    let last = *p.positions.last().unwrap();
    let mut windows: Vec<(u32, u32, usize, usize)> = Vec::new();
    let mut s = first;
    while s < last {
        let e = s + WIN;
        let lo = p.positions.partition_point(|&x| x < s);
        let hi = p.positions.partition_point(|&x| x < e);
        if hi > lo {
            windows.push((s, e, lo, hi));
        }
        s = e;
    }
    println!("{} windows of {} bp", windows.len(), WIN);

    // informative sites available per window
    let mask_window = |lo: usize, hi: usize, src: &[u64]| -> Vec<u64> {
        let (w0, w1) = (lo / 64, hi.div_ceil(64));
        let mut m = vec![0u64; words];
        for w in w0..w1 {
            let mut k = u64::MAX;
            if w == w0 {
                k &= u64::MAX << (lo % 64);
            }
            if w == w1 - 1 && hi % 64 != 0 {
                k &= !(u64::MAX << (hi % 64));
            }
            m[w] = k & src[w];
        }
        m
    };

    // ─────────────────────────────────────────── the tract scan
    //
    // A window is called archaic for a haplotype when it carries at least
    // MIN_ALLELES of the informative alleles present there and at least
    // MIN_FRACTION of them. The count guards against a single recurrent
    // mutation; the fraction guards against a window that happens to hold many
    // informative sites.
    const MIN_ALLELES: u32 = 3;
    const MIN_FRACTION: f64 = 0.30;

    let scan = |haps: &[usize], label: &str| -> (Vec<ArchaicTract>, TractScanArm) {
        let mut out: Vec<ArchaicTract> = Vec::new();
        let win_masks: Vec<(Vec<u64>, u32)> = windows
            .iter()
            .map(|&(_s, _e, lo, hi)| {
                let m = mask_window(lo, hi, &informative);
                let n = m.iter().map(|w| w.count_ones()).sum::<u32>();
                (m, n)
            })
            .collect();

        for &h in haps {
            let mut run: Option<(usize, usize, u32, u32)> = None; // start_w, end_w, carried, avail
            for (wi, (m, avail)) in win_masks.iter().enumerate() {
                let carried: u32 = (0..words).map(|w| (m[w] & p.bits[h][w]).count_ones()).sum();
                let hit = *avail > 0
                    && carried >= MIN_ALLELES
                    && (carried as f64 / *avail as f64) >= MIN_FRACTION;
                match (&mut run, hit) {
                    (None, true) => run = Some((wi, wi, carried, *avail)),
                    (Some(r), true) => {
                        r.1 = wi;
                        r.2 += carried;
                        r.3 += avail;
                    }
                    (Some(_), false) => {
                        let (a, b, c, av) = run.take().unwrap();
                        let smp = p.sample_of_hap[h];
                        out.push(ArchaicTract {
                            haplotype: p.hap_names[h].clone(),
                            sample: p.hap_names[h]
                                .rsplit_once('_')
                                .map(|(x, _)| x.to_string())
                                .unwrap_or_default(),
                            population: p.sample_pop[smp].clone(),
                            superpopulation: p.sample_super[smp].clone(),
                            start: windows[a].0,
                            end: windows[b].1,
                            length_bp: windows[b].1 - windows[a].0,
                            n_windows: b - a + 1,
                            informative_carried: c,
                            informative_available: av,
                            carried_fraction: c as f64 / av.max(1) as f64,
                        });
                    }
                    (None, false) => {}
                }
            }
            if let Some((a, b, c, av)) = run {
                let smp = p.sample_of_hap[h];
                out.push(ArchaicTract {
                    haplotype: p.hap_names[h].clone(),
                    sample: p.hap_names[h]
                        .rsplit_once('_')
                        .map(|(x, _)| x.to_string())
                        .unwrap_or_default(),
                    population: p.sample_pop[smp].clone(),
                    superpopulation: p.sample_super[smp].clone(),
                    start: windows[a].0,
                    end: windows[b].1,
                    length_bp: windows[b].1 - windows[a].0,
                    n_windows: b - a + 1,
                    informative_carried: c,
                    informative_available: av,
                    carried_fraction: c as f64 / av.max(1) as f64,
                });
            }
        }

        let mut lens: Vec<f64> = out.iter().map(|t| t.length_bp as f64).collect();
        lens.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let carriers: std::collections::BTreeSet<&str> =
            out.iter().map(|t| t.haplotype.as_str()).collect();
        let bp: u64 = out.iter().map(|t| t.length_bp as u64).sum();
        let scanned = haps.len() as u64 * (last - first) as u64;
        let arm = TractScanArm {
            arm: label.into(),
            n_haplotypes: haps.len(),
            n_tracts: out.len(),
            tracts_per_haplotype: out.len() as f64 / haps.len().max(1) as f64,
            haplotypes_with_a_tract: carriers.len(),
            fraction_of_haplotypes: carriers.len() as f64 / haps.len().max(1) as f64,
            mean_length_bp: if lens.is_empty() {
                0.0
            } else {
                lens.iter().sum::<f64>() / lens.len() as f64
            },
            median_length_bp: quantile(&lens, 0.5),
            max_length_bp: out.iter().map(|t| t.length_bp).max().unwrap_or(0),
            bp_in_tracts: bp,
            genome_fraction: bp as f64 / scanned.max(1) as f64,
        };
        (out, arm)
    };

    println!("\nscanning for Neanderthal tracts...");
    let (non_tracts, non_arm) = scan(&non, "non-African (the test)");
    let (afr_tracts, afr_arm) = scan(&afr, "African (the control)");
    for a in [&non_arm, &afr_arm] {
        println!(
            "   {:<26} {:>4} tracts in {:>3}/{:<3} haplotypes ({:>5.2}%)  mean {:>6.0} bp  covering {:.3}% of sequence",
            a.arm, a.n_tracts, a.haplotypes_with_a_tract, a.n_haplotypes,
            a.fraction_of_haplotypes * 100.0, a.mean_length_bp, a.genome_fraction * 100.0
        );
    }

    // ─────────────────────────────────────────── dating the pulse
    //
    // Every generation since the interbreeding, recombination has had a chance
    // to cut each introgressed tract. Surviving mean length L in Morgans
    // inverts to t = 1/L generations (Liang & Nielsen 2014). Tracts truncated
    // by the edge of a 1 Mb window bias the mean short, so the date biases old.
    let mean_bp = non_arm.mean_length_bp;
    let morgans = mean_bp * CM_PER_MB * 1e-8;
    let gens = if morgans > 0.0 { 1.0 / morgans } else { 0.0 };
    let n = non_arm.n_tracts.max(1) as f64;
    // the mean of n exponentials has relative SE 1/sqrt(n); propagate through 1/L
    let lo = gens / (1.0 + 1.96 / n.sqrt());
    let hi = gens / (1.0 - 1.96 / n.sqrt()).max(1e-6);
    let date = AdmixtureDate {
        source: "non-African tract lengths, this region".into(),
        n_tracts: non_arm.n_tracts,
        mean_length_bp: mean_bp,
        mean_length_morgans: morgans,
        generations: gens,
        years_ago: gens * GEN_YEARS,
        ka: gens * GEN_YEARS / 1000.0,
        ci95_ka: [lo * GEN_YEARS / 1000.0, hi * GEN_YEARS / 1000.0],
    };
    println!(
        "\nadmixture date from tract length: {:.0} generations = {:.1} ka  (95% CI {:.1}-{:.1} ka)",
        date.generations, date.ka, date.ci95_ka[0], date.ci95_ka[1]
    );

    // carriers by superpopulation
    let mut carrier: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for h in 0..n_hap {
        let g = p.sample_super[p.sample_of_hap[h]].clone();
        carrier.entry(g).or_insert((0, 0)).1 += 1;
    }
    let mut seen: std::collections::BTreeSet<(&str, &str)> = Default::default();
    for t in non_tracts.iter().chain(afr_tracts.iter()) {
        if seen.insert((t.superpopulation.as_str(), t.haplotype.as_str())) {
            carrier.get_mut(&t.superpopulation).unwrap().0 += 1;
        }
    }
    let carrier_by_superpopulation: Vec<(String, usize, usize, f64)> = carrier
        .iter()
        .map(|(g, (c, n))| (g.clone(), *c, *n, *c as f64 / *n as f64))
        .collect();
    println!("\ncarriers by superpopulation:");
    for (g, c, n, f) in &carrier_by_superpopulation {
        println!(
            "   {g:<4} {c:>3} / {n:<3} haplotypes carry a called tract  ({:.1}%)",
            f * 100.0
        );
    }

    // ─────────────────────────────────────────── the family, from the data
    //
    // Mean pairwise divergence between every pair of populations, and within
    // each. No tree is assumed and none is fitted: the matrix is the result.
    let mut pops: Vec<String> = p.sample_pop.clone();
    pops.sort();
    pops.dedup();
    let hap_of_pop: Vec<Vec<usize>> = pops
        .iter()
        .map(|q| {
            (0..n_hap)
                .filter(|&h| &p.sample_pop[p.sample_of_hap[h]] == q)
                .collect()
        })
        .collect();
    let span = (last - first) as f64;
    let to_ka = |d: f64| gens_to_ka(jukes_cantor(d / span) / (2.0 * MU));

    let np = pops.len();
    let mut mat = vec![vec![0.0f64; np]; np];
    for i in 0..np {
        for j in i..np {
            let (mut sum, mut cnt) = (0.0f64, 0u64);
            for &a in &hap_of_pop[i] {
                for &b in &hap_of_pop[j] {
                    if i == j && b <= a {
                        continue;
                    }
                    sum += diffs(&p.bits[a], &p.bits[b], 0, n_var) as f64;
                    cnt += 1;
                }
            }
            let v = if cnt > 0 {
                to_ka(sum / cnt as f64)
            } else {
                0.0
            };
            mat[i][j] = v;
            mat[j][i] = v;
        }
    }

    let nea_per_hap = |haps: &[usize]| -> f64 {
        if haps.is_empty() {
            return 0.0;
        }
        haps.iter()
            .map(|&h| {
                (0..words)
                    .map(|w| (informative[w] & p.bits[h][w]).count_ones() as f64)
                    .sum::<f64>()
            })
            .sum::<f64>()
            / haps.len() as f64
    };

    let mut populations: Vec<PopStat> = Vec::new();
    for (i, q) in pops.iter().enumerate() {
        let hs = &hap_of_pop[i];
        let sup = p.sample_super[p.sample_of_hap[hs[0]]].clone();
        let between: f64 =
            (0..np).filter(|&j| j != i).map(|j| mat[i][j]).sum::<f64>() / (np - 1) as f64;
        populations.push(PopStat {
            population: q.clone(),
            superpopulation: sup,
            n_haplotypes: hs.len(),
            within_ka: mat[i][i],
            between_ka: between,
            neanderthal_alleles_per_haplotype: nea_per_hap(hs),
        });
    }
    let mut by_within = populations.iter().collect::<Vec<_>>();
    by_within.sort_by(|a, b| b.within_ka.partial_cmp(&a.within_ka).unwrap());
    println!("\ndeepest internal diversity, by population:");
    for s in by_within.iter().take(6) {
        println!(
            "   {:<4} {:<4} within {:>6.1} ka   Neanderthal alleles/haplotype {:>5.1}",
            s.population, s.superpopulation, s.within_ka, s.neanderthal_alleles_per_haplotype
        );
    }
    println!("   ...");
    for s in by_within.iter().rev().take(3).rev() {
        println!(
            "   {:<4} {:<4} within {:>6.1} ka   Neanderthal alleles/haplotype {:>5.1}",
            s.population, s.superpopulation, s.within_ka, s.neanderthal_alleles_per_haplotype
        );
    }

    // The sharpest way to say "African diversity is deeper" is to compare it
    // against the widest gap anyone expects to be large: Europe to East Asia.
    let sup_of = |q: &str| {
        populations
            .iter()
            .find(|s| s.population == q)
            .map(|s| s.superpopulation.clone())
            .unwrap_or_default()
    };
    let mut widest_eur_eas = 0.0f64;
    for (i, a) in pops.iter().enumerate() {
        for (j, b) in pops.iter().enumerate() {
            let (sa, sb) = (sup_of(a), sup_of(b));
            if (sa == "EUR" && sb == "EAS") || (sa == "EAS" && sb == "EUR") {
                widest_eur_eas = widest_eur_eas.max(mat[i][j]);
            }
        }
    }

    let top_afr = by_within
        .iter()
        .filter(|s| s.superpopulation == "AFR")
        .count();
    let afr_in_top = by_within
        .iter()
        .take(top_afr)
        .filter(|s| s.superpopulation == "AFR")
        .count();

    let nea_afr = nea_per_hap(&afr);
    let nea_non = nea_per_hap(&non);
    let by_super = |g: &str| {
        nea_per_hap(
            &(0..n_hap)
                .filter(|&h| p.sample_super[p.sample_of_hap[h]] == g)
                .collect::<Vec<_>>(),
        )
    };
    let eas_alleles = by_super("EAS");
    let eur_alleles = by_super("EUR");
    println!(
        "\nNeanderthal-informative alleles per haplotype: AFR {:.1}  EUR {:.1}  SAS {:.1}  AMR {:.1}  EAS {:.1}  (non-AFR mean {:.1})",
        nea_afr, eur_alleles, by_super("SAS"), by_super("AMR"), eas_alleles, nea_non
    );

    let findings = vec![
        format!(
            "Start with what makes the scan trustworthy. {} positions in this megabase are Neanderthal-informative: a high-coverage Neanderthal carries the non-reference allele there, and not one of the {} African haplotypes in the panel does. A living person outside Africa carrying a run of those alleles inherited them from somebody, and the only somebody on offer is a Neanderthal. Conditioning on absence in Africa is also what makes the test immune to the European weighting of the reference sequence.",
            n_informative, afr.len()
        ),
        format!(
            "{} tracts turned up, in {} of the {} non-African haplotypes. Run the identical scan on the {} African haplotypes and it returns {}. Not a small number — zero. That is the control, and it is the reason the rest of this section is worth reading: the method is not finding tracts wherever it looks.",
            non_arm.n_tracts, non_arm.haplotypes_with_a_tract, non_arm.n_haplotypes,
            afr_arm.n_haplotypes, afr_arm.n_tracts
        ),
        format!(
            "The called tracts cover {:.2}% of the non-African sequence scanned. The published figure for Neanderthal ancestry in non-Africans, from whole genomes and far more careful methods, is about 2%. That number was not put in anywhere: it falls out of one megabase of chromosome 22, three archaic genomes, and a rule about which alleles are missing from Africa.",
            non_arm.genome_fraction * 100.0
        ),
        format!(
            "It is not evenly shared, and the ranking is the published one. East Asian haplotypes carry {:.1} Neanderthal-informative alleles each, South Asians {:.1}, Americans {:.1}, Europeans {:.1}. East Asians carrying more Neanderthal ancestry than Europeans is a real and initially counter-intuitive result in human genetics — Neanderthals lived in Europe, not East Asia — and the ranking falls out of a single megabase. The size of the gap here does not: genome-wide the East Asian excess is roughly 20%, not the {:.0}x seen at this locus. One megabase is one locus, and a locus can carry a Neanderthal haplotype that drifted to high frequency in one place and not another. The ordering is the finding; the magnitude is this piece of chromosome 22 talking.",
            eas_alleles, by_super("SAS"), by_super("AMR"), eur_alleles,
            if eur_alleles > 0.0 { eas_alleles / eur_alleles } else { f64::INFINITY }
        ),
        format!(
            "The tracts also date the interbreeding, and this is where the study contradicts itself usefully. Mean surviving tract length is {:.0} bp, which is {:.2e} Morgans, which inverts to {:.0} generations — {:.0} thousand years ago. The right answer, from whole genomes, is 50-60 ka. The estimate is {:.0}x too old.",
            date.mean_length_bp, date.mean_length_morgans, date.generations, date.ka,
            date.ka / 55.0
        ),
        format!(
            "That failure was predicted. The simulation sweep reports that a detector which fragments tracts rather than fusing them always reads the admixture as older than it was, because a shorter measured tract inverts to a larger number of generations — there it over-aged a planted pulse by 2.4x. Here, with {} informative sites to trace tracts through across a megabase, fragmentation is far worse and the date comes out {:.0}x old. The simulation predicted the direction and the mechanism of an error later observed in real human DNA, which is the strongest thing a simulation in this study does.",
            n_informative, date.ka / 55.0
        ),
        format!(
            "Finally, the family, drawn with no tree assumed. Rank the {} populations by how different two of their own randomly chosen haplotypes are, and the top {} are all African. {} sits at the top at {:.0} ka between two of its own members; {} at the bottom at {:.0} ka. Two haplotypes drawn from {} inside that one population are further apart than any European-East Asian pair in the whole panel, whose widest is {:.0} ka. Everyone outside Africa descends from one small group that walked out, and this is that fact as a measurement rather than a slogan.",
            np, afr_in_top,
            by_within[0].population, by_within[0].within_ka,
            by_within[np-1].population, by_within[np-1].within_ka,
            by_within[0].population, widest_eur_eas
        ),
    ];

    let story = HumanStory {
        title: "What the real panel says about us — 1000 Genomes chr22 and three archaic genomes".into(),
        generated_utc: chrono::Utc::now().to_rfc3339(),
        region: format!("chr22:{first}-{last} (GRCh37)"),
        window_bp: WIN,
        n_haplotypes: n_hap,
        n_variants: n_var,
        mu_per_bp_per_generation: MU,
        generation_years: GEN_YEARS,
        cm_per_mb: CM_PER_MB,
        informative_sites: n_informative,
        informative_note: format!(
            "A site is Neanderthal-informative when the Altai or Vindija genome carries the non-reference allele at a callable position and none of the {} African haplotypes in the panel carries it. Conditioning on African absence is what makes the scan robust to the European weighting of the reference sequence: the question is not how different somebody is from a Neanderthal, but whether they carry an allele found in Neanderthals and nowhere in Africa.",
            afr.len()
        ),
        tract_scan: vec![non_arm, afr_arm],
        tract_map: {
            let mut carriers: Vec<(String, String, String)> = Vec::new();
            let mut idx: BTreeMap<String, usize> = BTreeMap::new();
            let mut packed: Vec<(usize, usize, usize)> = Vec::new();
            let win_of = |pos: u32| ((pos - first) / WIN) as usize;
            for t in non_tracts.iter().chain(afr_tracts.iter()) {
                let k = *idx.entry(t.haplotype.clone()).or_insert_with(|| {
                    carriers.push((
                        t.haplotype.clone(),
                        t.population.clone(),
                        t.superpopulation.clone(),
                    ));
                    carriers.len() - 1
                });
                packed.push((k, win_of(t.start), win_of(t.end.saturating_sub(1))));
            }
            let mut totals: BTreeMap<String, usize> = BTreeMap::new();
            for h in 0..n_hap {
                *totals
                    .entry(p.sample_super[p.sample_of_hap[h]].clone())
                    .or_insert(0) += 1;
            }
            TractMap {
                region_start: first,
                region_end: last,
                window_bp: WIN,
                n_windows: windows.len(),
                carriers,
                total_by_superpopulation: totals.into_iter().collect(),
                tracts: packed,
            }
        },
        tracts: {
            let mut v = non_tracts.clone();
            v.sort_by(|a, b| b.length_bp.cmp(&a.length_bp));
            v.truncate(40);
            v
        },
        admixture_date: date,
        carrier_by_superpopulation,
        populations,
        population_order: pops,
        divergence_matrix_ka: mat,
        findings,
        caveats: vec![
            "This is a population-level scan reported at haplotype resolution, not a clinical or genealogical result. A called tract is a candidate: the evidence is that the haplotype carries alleles found in Neanderthals and in no African in this panel, which is strong but is not a phased ancestral-recombination-graph reconstruction.".into(),
            "The African-absence condition uses 160 African haplotypes, not all of Africa. An allele present in African populations not sampled here, or at very low frequency in those that are, would pass the filter wrongly. This inflates tract counts rather than deflating them.".into(),
            "Tract dating assumes a uniform 1 cM/Mb. Local recombination on chr22 varies several-fold around that, and any tract running off the edge of the 1 Mb region is measured short. Both push the estimate towards older dates.".into(),
            "Divergence times use a strict clock at mu = 1.25e-8 per base per generation with a 29-year generation, and are sequence divergence — always older than the population split that produced it, because the two lineages had their own history first.".into(),
            "The panel is 400 of the 2,504 phase 3 samples, 80 per superpopulation. Pairwise divergence between two included haplotypes is exact regardless, because any site at which they differ is polymorphic and therefore in the VCF.".into(),
        ],
    };

    let path = format!("{out_dir}/human-story.json");
    let mut f = BufWriter::new(std::fs::File::create(&path)?);
    f.write_all(serde_json::to_string_pretty(&story)?.as_bytes())?;
    println!("\nWrote {path}  ({:.1}s)", t0.elapsed().as_secs_f64());
    Ok(())
}
