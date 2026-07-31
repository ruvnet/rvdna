#!/usr/bin/env node
/**
 * Assembles artifacts/index.html from the template, the analysis output, and
 * the animated diagrams.
 *
 * Page weight is the binding constraint here — the artifact frame has to parse
 * everything synchronously — so this script does not simply inline the JSON the
 * Rust binaries emit:
 *
 *   - the 20,880 per-segment coalescent depths ship as a base64 Uint16Array
 *     rather than a JSON number array (253 KB -> ~56 KB);
 *   - archaic calls are reduced to the six fields the page actually reads;
 *   - the diagrams do NOT get a copy of the mono webfont each. Inlining it three
 *     times cost 126 KB to change the labels on three static illustrations.
 *
 * Usage: node artifacts/build.mjs
 */

import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const dataDir = join(here, 'data');

const need = (p) => {
  if (!existsSync(p)) {
    console.error(`missing: ${p}\nrun the Rust binaries first:\n  ./target/release/depth-hist artifacts/data\n  ./target/release/trace-rv  artifacts/data\n  ./target/release/discover  artifacts/data`);
    process.exit(1);
  }
  return p;
};
const readJson = (name) => JSON.parse(readFileSync(need(join(dataDir, name)), 'utf8'));

const report = readJson('trace-rv-report.json');
const calls = readJson('archaic-calls.json');
const depth = readJson('depth-distribution.json');
const disc = readJson('discoveries.json');

// The report embeds the full call list under `flywheel.final_calls`; the calls
// file is the same data. Carry it once.
delete report.flywheel.final_calls;

// ---------------------------------------------------------------- particles
//
// One particle per (window, modern haplotype). Everything the 3D narrative
// needs about a segment is packed into three parallel typed arrays.

const roster = depth.roster.filter((r) => r.role === 'modern');
const hapCol = new Map(depth.depth_map_haplotypes.map((h, i) => [h, i]));
const nW = depth.n_windows;
const nH = roster.length;
const N = nW * nH;

const depths = new Uint16Array(N);        // coalescent depth in ka, capped
const klass = new Uint8Array(N);          // what the detector concluded
const truthK = new Uint8Array(N);         // what the segment actually is

// 0 nothing · 1 Neanderthal · 2 Denisovan · 3 ghost mode 0 · 4 ghost mode 1
const CALL_CLASS = { Neanderthal: 1, Denisovan: 2 };
const TRUTH_CLASS = { '.': 0, N: 1, D: 2, A: 3, B: 4 };

const callAt = new Map();
for (const c of calls) callAt.set(c.window * 4096 + c.haplotype, c);

for (let w = 0; w < nW; w++) {
  const truthRow = depth.truth_map[w] || '';
  for (let r = 0; r < nH; r++) {
    const i = w * nH + r;
    const h = roster[r].index;
    depths[i] = Math.min(65535, depth.depth_map[w][hapCol.get(h)] | 0);
    truthK[i] = TRUTH_CLASS[truthRow[h]] ?? 0;
    const c = callAt.get(w * 4096 + h);
    klass[i] = !c ? 0 : (CALL_CLASS[c.attribution] ?? (c.ghost_cluster === 1 ? 4 : 3));
  }
}

const b64 = (buf) => Buffer.from(buf.buffer, buf.byteOffset, buf.byteLength).toString('base64');

const particles = {
  n: N, windows: nW, haplotypes: nH,
  depths: b64(depths),
  calls: b64(klass),
  truth: b64(truthK),
  roster: roster.map((r) => ({ id: r.id, pop: r.population, grp: r.group })),
};

// The charts still need per-call detail, but only these fields.
const callsSlim = calls.map((c) => ({
  w: c.window, h: c.haplotype, id: c.haplotype_id, g: c.group,
  a: c.attribution, d: Math.round(c.depth_to_modern_ka), c: c.ghost_cluster,
}));

// The depth histogram keeps its series and overlap; the bulky maps are now in
// `particles`.
const depthSlim = {
  bin_edges_ka: depth.bin_edges_ka,
  series: depth.series,
  overlap: depth.overlap,
  n_segments: depth.n_segments,
};

// ---------------------------------------------------------------- assets
const font = (name) => {
  const buf = readFileSync(need(join(here, 'fonts', name)));
  return `data:font/woff2;base64,${buf.toString('base64')}`;
};

const diagram = (name) => {
  const svg = readFileSync(need(join(here, 'diagrams', name)), 'utf8');
  return `data:image/svg+xml;base64,${Buffer.from(svg, 'utf8').toString('base64')}`;
};

// ---------------------------------------------------------------- assemble
let html = readFileSync(join(here, 'story.template.html'), 'utf8');

const subs = {
  __TRACE_DATA__: JSON.stringify(report),
  __CALLS_DATA__: JSON.stringify(callsSlim),
  __DEPTH_DATA__: JSON.stringify(depthSlim),
  __DISC_DATA__: JSON.stringify(disc),
  __PARTICLE_DATA__: JSON.stringify(particles),
  __FONT_OUTFIT__: font('Outfit.woff2'),
  __FONT_MONO__: font('JetBrainsMono.woff2'),
  __DIAGRAM_01__: diagram('01-deep-time-tree.svg'),
  __DIAGRAM_04__: diagram('04-darwin-flywheel.svg'),
};

for (const [k, v] of Object.entries(subs)) {
  if (!html.includes(k)) {
    console.error(`template has no placeholder ${k}`);
    process.exit(1);
  }
  // Function form so `$&` / `$1` inside the payload are not treated as
  // replacement patterns — the JSON is full of `$`-adjacent text.
  html = html.replaceAll(k, () => v);
}

const leftover = html.match(/__[A-Z0-9_]+__/g);
if (leftover) {
  console.error(`unresolved placeholders: ${[...new Set(leftover)].join(', ')}`);
  process.exit(1);
}

const out = join(here, 'index.html');
writeFileSync(out, html);

const kb = (n) => (n / 1024).toFixed(0) + ' KB';
const total = Buffer.byteLength(html);
console.log(`wrote ${out}  ${kb(total)}`);
console.log(`  particles ${kb(JSON.stringify(particles).length)} (${N} segments)`);
console.log(`  report ${kb(subs.__TRACE_DATA__.length)} · calls ${kb(subs.__CALLS_DATA__.length)} · depth ${kb(subs.__DEPTH_DATA__.length)} · discoveries ${kb(subs.__DISC_DATA__.length)}`);
console.log(`  fonts ${kb(subs.__FONT_OUTFIT__.length + subs.__FONT_MONO__.length)} · diagrams ${kb(subs.__DIAGRAM_01__.length + subs.__DIAGRAM_04__.length)}`);
if (total > 420 * 1024) console.warn(`  WARNING: ${kb(total)} is heavy for a single artifact page`);
