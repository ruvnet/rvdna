"""Fetch REAL archaic hominin genotypes for chr22:20,000,000-21,000,000 (GRCh37).

Source: MPI-EVA cdna.eva.mpg.de, snpAD genotype calls (Prufer et al. 2017)
        https://cdna.eva.mpg.de/neandertal/Vindija/VCF/{Altai,Vindija33.19,Denisova}/

Method: no tabix/bcftools available, so this parses the .tbi tabix index in pure
Python, converts the virtual offsets for the target region into byte offsets,
pulls only those bytes with HTTP Range requests, and inflates the BGZF blocks
with zlib. Stops as soon as POS passes the end of the region.
"""
import json, os, sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-rvdna/db34c107-a32d-5b63-a37a-dfc1856a1402/scratchpad')
from tbx import parse_tbi, chunks_for, http_range, BGZFStream

BASE = 'https://cdna.eva.mpg.de/neandertal/Vindija/VCF'
# (output label, server directory)
SAMPLES = [('AltaiNeandertal', 'Altai'),
           ('Vindija33.19', 'Vindija33.19'),
           ('Denisova3', 'Denisova')]
BEG, END = 20000000, 21000000          # 1-based inclusive
SCR = '/tmp/claude-0/-home-user-rvdna/db34c107-a32d-5b63-a37a-dfc1856a1402/scratchpad'
OUT = '/home/user/rvdna/data/real'
MODERN = os.path.join(OUT, 'chr22_region.haplotypes.tsv')
SLICE = 4 << 20


def url_for(d):
    return '%s/%s/chr22_mq25_mapab100.vcf.gz' % (BASE, d)


def plan(d):
    idx = parse_tbi(os.path.join(SCR, '%s.chr22.tbi' % d))
    ri = idx['names'].index('22')
    ch = chunks_for(idx, ri, BEG - 1, END)
    return min(c[0] >> 16 for c in ch), max(c[1] >> 16 for c in ch)


def load_modern_positions():
    """Positions (and REF/ALT) of the modern 1000G panel already fetched."""
    pos = {}
    if not os.path.exists(MODERN):
        return pos
    with open(MODERN) as fh:
        fh.readline()
        for line in fh:
            p, r, a = line.split('\t', 3)[:3]
            pos[int(p)] = (r, a)
    return pos


def scan(d, label, stats):
    """Stream the region for one individual; return {pos: (ref, alt, gt_str, gq)}."""
    url = url_for(d)
    lo, hi = plan(d)
    stream = BGZFStream()
    tail = b''
    got = {}
    downloaded = 0
    off = lo
    done = False
    while off <= hi and not done:
        end = min(off + SLICE - 1, hi)
        raw = http_range(url, off, end)
        if not raw:
            break
        downloaded += len(raw)
        off += len(raw)
        text = tail + stream.feed(raw)
        nl = text.rfind(b'\n')
        if nl < 0:
            tail = text
            continue
        tail = text[nl + 1:]
        for line in text[:nl].split(b'\n'):
            if not line or line[0:1] == b'#':
                continue
            f = line.split(b'\t')
            if len(f) < 10:
                continue
            p = int(f[1])
            if p < BEG:
                continue
            if p > END:
                done = True
                break
            fmt = f[8].split(b':')
            val = f[9].split(b':')
            rec = dict(zip(fmt, val))
            gt = rec.get(b'GT', b'./.').decode()
            gq = rec.get(b'GQ', b'.').decode()
            filt = f[6].decode()
            got[p] = (f[3].decode(), f[4].decode(), gt, gq, filt)
    stats[label] = {'url': url, 'byte_range': [lo, hi], 'bytes_downloaded': downloaded,
                    'sites_in_region': len(got)}
    sys.stderr.write('%s: %d sites, %.2f MB downloaded\n' % (label, len(got), downloaded / 1e6))
    return got


def main():
    t0 = time.time()
    stats = {}
    data = {}
    for label, d in SAMPLES:
        data[label] = scan(d, label, stats)

    modern = load_modern_positions()
    labels = [s[0] for s in SAMPLES]

    # union of: any archaic site with a called ALT, plus every modern panel position
    interesting = set(modern)
    for lab in labels:
        for p, (r, a, gt, gq, filt) in data[lab].items():
            if a != '.':
                interesting.add(p)

    n_third = 0
    n_skipped_indel = 0
    rows = []
    for p in sorted(interesting):
        # resolve REF from whichever source has it
        ref = None
        for lab in labels:
            if p in data[lab]:
                ref = data[lab][p][0]
                break
        if ref is None and p in modern:
            ref = modern[p][0]
        if ref is None:
            continue
        # resolve ALT: prefer modern panel allele, else archaic-observed allele
        alt = None
        if p in modern and modern[p][1] not in ('.', ''):
            alt = modern[p][1]
        if alt is None:
            for lab in labels:
                if p in data[lab] and data[lab][p][1] != '.':
                    alt = data[lab][p][1].split(',')[0]
                    break
        if alt is None:
            continue
        if len(ref) != 1 or len(alt) != 1 or ',' in alt:
            n_skipped_indel += 1
            continue

        cols = []
        for lab in labels:
            if p not in data[lab]:
                cols.append('.')
                continue
            r, a, gt, gq, filt = data[lab][p]
            if r != ref or gt in ('./.', '.', '.|.') or filt not in ('.', 'PASS'):
                cols.append('.')
                continue
            alleles = [r] + ([] if a == '.' else a.split(','))
            sep = '|' if '|' in gt else '/'
            try:
                idxs = [int(x) for x in gt.split(sep)]
            except ValueError:
                cols.append('.')
                continue
            if any(i >= len(alleles) for i in idxs):
                cols.append('.')
                continue
            obs = [alleles[i] for i in idxs]
            if any(b != ref and b != alt for b in obs):
                n_third += 1
                cols.append('.')
                continue
            cols.append(str(sum(1 for b in obs if b == alt)))
        if all(c == '.' for c in cols):
            continue
        rows.append((p, ref, alt, cols))

    tsv = os.path.join(OUT, 'archaic_chr22_region.tsv')
    with open(tsv, 'w') as fh:
        fh.write('pos\tref\talt\t' + '\t'.join(labels) + '\n')
        for p, ref, alt, cols in rows:
            fh.write('%d\t%s\t%s\t%s\n' % (p, ref, alt, '\t'.join(cols)))

    # summary counts per individual
    per = {}
    for i, lab in enumerate(labels):
        c = {'0': 0, '1': 0, '2': 0, '.': 0}
        for _, _, _, cols in rows:
            c[cols[i]] = c.get(cols[i], 0) + 1
        per[lab] = c

    meta = {
        'individuals': labels,
        'source_urls': [stats[l]['url'] for l in labels] +
                       ['%s/README' % BASE] +
                       ['%s.tbi' % stats[l]['url'] for l in labels],
        'region': {'chrom': '22', 'start': BEG, 'end': END,
                   'string': 'chr22:%d-%d' % (BEG, END), 'coordinates': '1-based inclusive'},
        'build': 'GRCh37/hg19',
        'coverage_note': (
            'snpAD genotype calls (Prufer et al. 2017, Science 358:655) from MPI-EVA. '
            'Altai Neandertal ~52x, Vindija 33.19 ~30x, Denisova 3 ~30x nuclear coverage. '
            'Calls restricted by the producers to mappable sites (Heng Li 35mer filter) with MQ>=25. '
            'These are all-sites VCFs, so a position absent from a file is genuinely uncallable in '
            'that individual and is reported as "." (no-call). Genotypes are the raw snpAD GT '
            'converted to alt-allele dosage; NO extra GQ/DP filter was applied here, so low-GQ '
            'calls are included as reported by the producers. Rows are the union of (a) every '
            'position in the modern 1000G panel file chr22_region.haplotypes.tsv and (b) every '
            'position where at least one archaic carries a called ALT allele; rows where all three '
            'individuals are no-call are dropped. Only biallelic SNVs are emitted. Where an archaic '
            'carries a third allele not matching the chosen REF/ALT pair, that individual is set to '
            '"." (%d such individual-site calls).' % n_third),
        'per_individual_dosage_counts': per,
        'n_sites': len(rows),
        'n_sites_skipped_non_snv': n_skipped_indel,
        'retrieval_method': ('pure-Python tabix (.tbi) index parse -> virtual offsets -> HTTP Range '
                             'requests -> BGZF block inflate with zlib; no tabix/bcftools available'),
        'bytes_downloaded': {l: stats[l]['bytes_downloaded'] for l in labels},
        'total_bytes_downloaded': sum(stats[l]['bytes_downloaded'] for l in labels),
        'sites_callable_in_region': {l: stats[l]['sites_in_region'] for l in labels},
        'retrieved_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    }
    with open(os.path.join(OUT, 'archaic_chr22_region.meta.json'), 'w') as fh:
        json.dump(meta, fh, indent=2)
    sys.stderr.write('rows=%d  elapsed=%.0fs\n' % (len(rows), time.time() - t0))


if __name__ == '__main__':
    main()
