//! Loading and comparing the real genome panel.
//!
//! Shared by every binary that works on real data rather than on a simulated
//! cohort: the 1000 Genomes phased haplotype table, the archaic genotype table
//! aligned onto the same variant index, and the two bit-parallel comparison
//! primitives everything else is built out of.
//!
//! Haplotypes are stored as one bit per variant per haplotype, so a pairwise
//! comparison over a window is a masked XOR and a popcount rather than a walk
//! over a genotype array. That is what makes an exhaustive all-pairs scan the
//! cheap option on this data instead of the expensive one.

use std::collections::BTreeMap;
use std::io::BufRead;

pub struct Panel {
    pub positions: Vec<u32>,
    /// `(ref, alt)` per variant — needed to check that the archaic call set
    /// agrees with the panel about what the alleles at a site are.
    pub alleles: Vec<(String, String)>,
    /// `bits[hap]` — one bit per variant, set when that haplotype carries ALT.
    pub bits: Vec<Vec<u64>>,
    pub hap_names: Vec<String>,
    pub sample_of_hap: Vec<usize>,
    pub sample_pop: Vec<String>,
    pub sample_super: Vec<String>,
}

pub fn load(dir: &str) -> anyhow::Result<Panel> {
    // --- samples ---------------------------------------------------------
    let mut sample_pop = Vec::new();
    let mut sample_super = Vec::new();
    let mut sample_name = Vec::new();
    let f = std::fs::File::open(format!("{dir}/chr22_region.samples.tsv"))?;
    for (i, line) in std::io::BufReader::new(f).lines().enumerate() {
        let line = line?;
        if i == 0 {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 3 {
            continue;
        }
        sample_name.push(c[0].to_string());
        sample_pop.push(c[1].to_string());
        sample_super.push(c[2].to_string());
    }

    // --- haplotypes ------------------------------------------------------
    let f = std::fs::File::open(format!("{dir}/chr22_region.haplotypes.tsv"))?;
    let mut rdr = std::io::BufReader::with_capacity(1 << 20, f);
    let mut header = String::new();
    rdr.read_line(&mut header)?;
    let hap_names: Vec<String> = header
        .trim_end()
        .split('\t')
        .skip(3)
        .map(|s| s.to_string())
        .collect();
    let n_hap = hap_names.len();

    // Map each haplotype column back to its sample row.
    let index_of: BTreeMap<&str, usize> = sample_name
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let sample_of_hap: Vec<usize> = hap_names
        .iter()
        .map(|h| {
            let base = h.rsplit_once('_').map(|(a, _)| a).unwrap_or(h);
            *index_of.get(base).unwrap_or(&0)
        })
        .collect();

    let mut positions: Vec<u32> = Vec::new();
    let mut alleles: Vec<(String, String)> = Vec::new();
    let mut rows: Vec<Vec<bool>> = Vec::new();
    let mut line = String::new();
    while {
        line.clear();
        rdr.read_line(&mut line)? > 0
    } {
        let t = line.trim_end();
        if t.is_empty() {
            continue;
        }
        let mut it = t.split('\t');
        let pos: u32 = match it.next().and_then(|s| s.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        let r = it.next().unwrap_or("").to_string();
        let a = it.next().unwrap_or("").to_string();
        let row: Vec<bool> = it.map(|g| g == "1").collect();
        if row.len() != n_hap {
            continue;
        }
        positions.push(pos);
        alleles.push((r, a));
        rows.push(row);
    }

    let n_var = positions.len();
    let words = n_var.div_ceil(64);
    let mut bits = vec![vec![0u64; words]; n_hap];
    for (v, row) in rows.iter().enumerate() {
        let (w, b) = (v / 64, v % 64);
        for (h, &alt) in row.iter().enumerate() {
            if alt {
                bits[h][w] |= 1u64 << b;
            }
        }
    }

    Ok(Panel {
        positions,
        alleles,
        bits,
        hap_names,
        sample_of_hap,
        sample_pop,
        sample_super,
    })
}

/// One sequenced archaic genome, aligned 1:1 onto the panel's variant index.
///
/// snpAD genotype calls are diploid, so a site is stored as *which alleles this
/// individual carries* rather than as a dosage: a heterozygote carries both and
/// therefore matches a modern haplotype whichever allele that haplotype has.
/// Sites the producers could not call are excluded from both numerator and
/// denominator — an all-sites VCF means absence is real uncallability, not a
/// missing record.
pub struct Archaic {
    pub name: String,
    pub has_ref: Vec<u64>,
    pub has_alt: Vec<u64>,
    pub callable: Vec<u64>,
    /// Callable bases in the whole region, from the producers' own count. The
    /// union table only stores *variable* sites, so the denominator for a
    /// per-base divergence has to come from here.
    pub callable_bp_region: f64,
}

/// Loads the archaic genotype table and aligns it to `positions`.
///
/// Rows whose REF/ALT pair disagrees with the panel are dropped rather than
/// reconciled: a mismatched pair means the two call sets disagree about what
/// the alleles at that site even are, and guessing would silently invent
/// divergence.
pub fn load_archaic(
    dir: &str,
    positions: &[u32],
    panel_alleles: &[(String, String)],
) -> anyhow::Result<Vec<Archaic>> {
    let path = format!("{dir}/archaic_chr22_region.tsv");
    let f = std::fs::File::open(&path)?;
    let mut rdr = std::io::BufReader::with_capacity(1 << 20, f);
    let mut header = String::new();
    rdr.read_line(&mut header)?;
    let names: Vec<String> = header
        .trim_end()
        .split('\t')
        .skip(3)
        .map(String::from)
        .collect();

    let words = positions.len().div_ceil(64);
    let mut arcs: Vec<Archaic> = names
        .iter()
        .map(|n| Archaic {
            name: n.clone(),
            has_ref: vec![0u64; words],
            has_alt: vec![0u64; words],
            callable: vec![0u64; words],
            callable_bp_region: 0.0,
        })
        .collect();

    let index_of: BTreeMap<u32, usize> =
        positions.iter().enumerate().map(|(i, &p)| (p, i)).collect();

    let mut line = String::new();
    while {
        line.clear();
        rdr.read_line(&mut line)? > 0
    } {
        let t = line.trim_end();
        if t.is_empty() {
            continue;
        }
        let mut it = t.split('\t');
        let pos: u32 = match it.next().and_then(|s| s.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        let Some(&v) = index_of.get(&pos) else {
            continue;
        };
        let (r, a) = (it.next().unwrap_or(""), it.next().unwrap_or(""));
        if panel_alleles[v].0 != r || panel_alleles[v].1 != a {
            continue;
        }
        let (w, b) = (v / 64, 1u64 << (v % 64));
        for (k, g) in it.enumerate() {
            if k >= arcs.len() {
                break;
            }
            match g {
                "0" => {
                    arcs[k].has_ref[w] |= b;
                    arcs[k].callable[w] |= b;
                }
                "1" => {
                    arcs[k].has_ref[w] |= b;
                    arcs[k].has_alt[w] |= b;
                    arcs[k].callable[w] |= b;
                }
                "2" => {
                    arcs[k].has_alt[w] |= b;
                    arcs[k].callable[w] |= b;
                }
                _ => {}
            }
        }
    }
    Ok(arcs)
}

/// Sites at which a modern haplotype carries an allele the archaic individual
/// does not, over a half-open variant range, counting only callable sites.
///
/// Returns `(mismatches, callable_sites)`.
#[inline]
pub fn arc_diffs(hap: &[u64], arc: &Archaic, lo: usize, hi: usize) -> (u32, u32) {
    let (w0, w1) = (lo / 64, hi.div_ceil(64));
    let (mut d, mut c) = (0u32, 0u32);
    for w in w0..w1 {
        let mut mask = u64::MAX;
        if w == w0 {
            mask &= u64::MAX << (lo % 64);
        }
        if w == w1 - 1 && hi % 64 != 0 {
            mask &= !(u64::MAX << (hi % 64));
        }
        let cal = arc.callable[w] & mask;
        let a = hap[w];
        // carries ALT but the archaic has no ALT · carries REF but the archaic has no REF
        d += (cal & ((a & !arc.has_alt[w]) | (!a & !arc.has_ref[w]))).count_ones();
        c += cal.count_ones();
    }
    (d, c)
}

/// Differing sites between two haplotypes over a half-open variant range.
#[inline]
pub fn diffs(a: &[u64], b: &[u64], lo: usize, hi: usize) -> u32 {
    let (w0, w1) = (lo / 64, hi.div_ceil(64));
    let mut n = 0u32;
    for w in w0..w1 {
        let mut x = a[w] ^ b[w];
        if w == w0 {
            let sh = lo % 64;
            x &= u64::MAX << sh;
        }
        if w == w1 - 1 && hi % 64 != 0 {
            x &= !(u64::MAX << (hi % 64));
        }
        n += x.count_ones();
    }
    n
}
