#!/usr/bin/env node
/**
 * Assembles artifacts/index.html from the template, the analysis output, and
 * the animated diagrams.
 *
 * The diagrams are inlined as base64 data URIs inside <img> elements rather
 * than as inline <svg>. That keeps each diagram's stylesheet in its own
 * document — they all use short class names like `.bg` and `.note` that would
 * otherwise collide with each other and with the page — while still animating,
 * because CSS animations run in img-referenced SVG.
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
    console.error(`missing: ${p}\nrun the Rust binaries first:\n  ./target/release/depth-hist artifacts/data\n  ./target/release/trace-rv  artifacts/data`);
    process.exit(1);
  }
  return p;
};

const reportObj = JSON.parse(readFileSync(need(join(dataDir, 'trace-rv-report.json')), 'utf8'));
const calls = readFileSync(need(join(dataDir, 'archaic-calls.json')), 'utf8');
const depth = readFileSync(need(join(dataDir, 'depth-distribution.json')), 'utf8');
const disc = readFileSync(need(join(dataDir, 'discoveries.json')), 'utf8');

// The report embeds the full call list under `flywheel.final_calls`, and
// archaic-calls.json is the same data. Inline it once.
delete reportObj.flywheel.final_calls;
const report = JSON.stringify(reportObj);

// Outfit and JetBrains Mono, both SIL OFL, embedded as data URIs. The Artifact
// CSP blocks font CDNs, and a silent fallback would lose the typography the
// design system is built on.
const font = (name) => {
  const buf = readFileSync(need(join(here, 'fonts', name)));
  return `data:font/woff2;base64,${buf.toString('base64')}`;
};

// An <img>-referenced SVG is its own document: its stylesheet cannot reach the
// page's @font-face, and every diagram uses short class names that would collide
// if inlined. So each diagram carries its own copy of the mono face.
const monoFace = `@font-face{font-family:'JetBrains Mono';font-style:normal;font-weight:100 800;src:url(${font('JetBrainsMono.woff2')}) format('woff2');}`;

const diagram = (name) => {
  let svg = readFileSync(need(join(here, 'diagrams', name)), 'utf8');
  if (!svg.includes('<style>')) throw new Error(`${name} has no <style> block to inject the face into`);
  svg = svg.replace('<style>', `<style>\n${monoFace}\n`);
  return `data:image/svg+xml;base64,${Buffer.from(svg, 'utf8').toString('base64')}`;
};

let html = readFileSync(join(here, 'story.template.html'), 'utf8');

const subs = {
  __TRACE_DATA__: report,
  __CALLS_DATA__: calls,
  __DEPTH_DATA__: depth,
  __DISC_DATA__: disc,
  __FONT_OUTFIT__: font('Outfit.woff2'),
  __FONT_MONO__: font('JetBrainsMono.woff2'),
  __DIAGRAM_01__: diagram('01-deep-time-tree.svg'),
  __DIAGRAM_03__: diagram('03-hnsw-retrieval.svg'),
  __DIAGRAM_04__: diagram('04-darwin-flywheel.svg'),
};

for (const [k, v] of Object.entries(subs)) {
  if (!html.includes(k)) {
    console.error(`template has no placeholder ${k}`);
    process.exit(1);
  }
  // Function form so `$&`, `$1` etc. in the payload are not treated as
  // replacement patterns — the JSON contains plenty of `$`-adjacent text.
  html = html.replaceAll(k, () => v);
}

const leftover = html.match(/__[A-Z0-9_]+__/g);
if (leftover) {
  console.error(`unresolved placeholders: ${[...new Set(leftover)].join(', ')}`);
  process.exit(1);
}

const out = join(here, 'index.html');
writeFileSync(out, html);

const kb = (n) => (n / 1024).toFixed(0);
console.log(`wrote ${out}  (${kb(Buffer.byteLength(html))} KB)`);
console.log(`  report ${kb(report.length)} KB · calls ${kb(calls.length)} KB · depth ${kb(depth.length)} KB`);
