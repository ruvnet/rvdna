// =============================================================================
// DNA Storage Codec Simulator — UI wiring, rendering & visualization
// =============================================================================
// Glue between the pure codec (codec.js) and the DOM. No frameworks, no network.

import {
  BASES, encodeBytes, decodeBytes,
  crc32, applyChannel,
  DnaStorageCodec, DEFAULT_PARAMS,
  totalBases, bitsPerBase, meanGc, meanMaxHomopolymer,
  selfTest,
} from './codec.js';

// A tiny 32x32 colorful PNG embedded as a data URI so "Use sample" works offline.
const SAMPLE_PNG_B64 =
  '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';
// NOTE: the long b64 above is a literal 32x32 RGB PNG (built offline) so the
// "Use sample PNG" button works with zero network access.

// -----------------------------------------------------------------------------
// App state
// -----------------------------------------------------------------------------
const state = {
  fileBytes: null,     // Uint8Array of the loaded file
  fileName: '',
  mimeType: '',        // best-effort
  isImage: false,
  archive: null,       // result of codec.encode()
  reads: null,         // mutated read pool
  report: null,        // result of codec.decode()
  codec: new DnaStorageCodec(DEFAULT_PARAMS),
};

const $ = (id) => document.getElementById(id);
const fmt = (n) => n.toLocaleString('en-US');
const pct = (x) => (x * 100).toFixed(1) + '%';

const BASE_COLOR = { A: '#27d6a4', C: '#4aa8ff', G: '#ffb454', T: '#ff5d8f' };

// -----------------------------------------------------------------------------
// File loading
// -----------------------------------------------------------------------------
function detectImage(name, bytes) {
  const lower = name.toLowerCase();
  if (lower.endsWith('.png')) return 'image/png';
  if (lower.endsWith('.jpg') || lower.endsWith('.jpeg')) return 'image/jpeg';
  if (lower.endsWith('.gif')) return 'image/gif';
  if (lower.endsWith('.webp')) return 'image/webp';
  // magic-byte sniff
  if (bytes.length > 8 && bytes[0] === 0x89 && bytes[1] === 0x50) return 'image/png';
  if (bytes.length > 3 && bytes[0] === 0xff && bytes[1] === 0xd8) return 'image/jpeg';
  return '';
}

function loadBytes(name, bytes) {
  state.fileBytes = bytes;
  state.fileName = name;
  state.mimeType = detectImage(name, bytes);
  state.isImage = state.mimeType.startsWith('image/');
  // reset downstream state
  state.archive = state.reads = state.report = null;
  $('encodeBtn').disabled = false;
  $('mutateBtn').disabled = true;
  $('recoverBtn').disabled = true;
  $('fileInfo').textContent =
    `Loaded "${name}" — ${fmt(bytes.length)} bytes${state.isImage ? ' (image)' : ''}`;
  $('stats').innerHTML = '';
  $('strandPreview').innerHTML = '';
  $('mutateInfo').innerHTML = '';
  $('recoverBanner').className = 'banner hidden';
  $('recoverInfo').innerHTML = '';
  renderImagePanels(); // show "before" if image
  drawStrandAnimation(bytes);
}

function readFile(file) {
  const reader = new FileReader();
  reader.onload = () => loadBytes(file.name, new Uint8Array(reader.result));
  reader.readAsArrayBuffer(file);
}

function loadSamplePng() {
  const bin = atob(SAMPLE_PNG_B64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  loadBytes('sample.png', bytes);
}

// -----------------------------------------------------------------------------
// Encode
// -----------------------------------------------------------------------------
function doEncode() {
  state.codec = new DnaStorageCodec(readParams());
  state.archive = state.codec.encode(state.fileName, state.fileBytes);
  state.reads = null;
  state.report = null;
  renderStats();
  renderStrandPreview();
  $('mutateBtn').disabled = false;
  $('recoverBtn').disabled = true;
  renderImagePanels();
}

function readParams() {
  return {
    ...DEFAULT_PARAMS,
    blockSize: +$('blockSize').value,
    rsParity: +$('rsParity').value,
    overhead: +$('overhead').value,
  };
}

function renderStats() {
  const a = state.archive;
  const bases = totalBases(a);
  const strandLen = a.strands.length ? a.strands[0].sequence.length : 0;
  const rows = [
    ['Original size', `${fmt(a.byteLen)} bytes`],
    ['CRC32', '0x' + a.crc32.toString(16).padStart(8, '0').toUpperCase()],
    ['Source blocks', fmt(a.numBlocks)],
    ['Strands (oligos)', fmt(a.strands.length)],
    ['Strand length', `${fmt(strandLen)} bases`],
    ['Total bases', fmt(bases)],
    ['Density', bitsPerBase(a).toFixed(3) + ' bits/base'],
    ['Mean GC', pct(meanGc(a))],
    ['Max homopolymer run', String(meanMaxHomopolymer(a))],
    ['RS parity', `${a.params.rsParity} bytes (corrects ${Math.floor(a.params.rsParity / 2)}/strand)`],
    ['Redundancy', `${Math.ceil(a.params.overhead)}× (dropout tolerance)`],
  ];
  $('stats').innerHTML = rows
    .map(([k, v]) => `<div class="stat"><span class="k">${k}</span><span class="v">${v}</span></div>`)
    .join('');
}

// Render a sequence as colored base spans. `diffSeq` optionally highlights diffs.
function renderSeq(seq, diffSeq) {
  let html = '';
  for (let i = 0; i < seq.length; i++) {
    const b = seq[i];
    const changed = diffSeq && (i >= diffSeq.length || diffSeq[i] !== b);
    const cls = changed ? 'base diff' : 'base';
    const color = BASE_COLOR[b] || '#888';
    html += `<span class="${cls}" style="--bc:${color}">${b}</span>`;
  }
  return html;
}

function renderStrandPreview() {
  const a = state.archive;
  const n = Math.min(6, a.strands.length);
  let html = '';
  for (let i = 0; i < n; i++) {
    const s = a.strands[i];
    html += `<div class="strand-row">
      <span class="strand-idx">#${s.index}</span>
      <span class="strand-seq">${renderSeq(s.sequence)}</span>
    </div>`;
  }
  $('strandPreview').innerHTML =
    `<div class="legend">${BASES.map(b => `<span class="lg" style="--bc:${BASE_COLOR[b]}">${b}</span>`).join('')}
      <span class="muted">showing ${n} of ${fmt(a.strands.length)} strands</span></div>` + html;
}

// -----------------------------------------------------------------------------
// Mutate (channel)
// -----------------------------------------------------------------------------
function doMutate() {
  const model = {
    pSub: +$('pSub').value / 100,
    pIns: +$('pIns').value / 100,
    pDel: +$('pDel').value / 100,
    pDrop: +$('pDrop').value / 100,
    coverage: +$('coverage').value,
  };
  const seed = (Date.now() & 0xffffff) ^ 0xc0ffee;
  state.reads = applyChannel(state.archive.strands, model, seed);
  state.report = null;

  // measure damage: compare each read to its true strand (by index)
  const byIdx = new Map(state.archive.strands.map(s => [s.index, s.sequence]));
  let basesHit = 0, strandsHit = 0, totalReadBases = 0;
  for (const r of state.reads) {
    const orig = byIdx.get(r.index) || '';
    totalReadBases += r.sequence.length;
    let hit = 0;
    const m = Math.min(orig.length, r.sequence.length);
    for (let i = 0; i < m; i++) if (orig[i] !== r.sequence[i]) hit++;
    hit += Math.abs(orig.length - r.sequence.length); // indel-induced shifts
    if (hit > 0) strandsHit++;
    basesHit += hit;
  }
  const dropped = state.archive.strands.length * Math.max(1, model.coverage) - state.reads.length;

  $('mutateInfo').innerHTML = `
    <div class="stat"><span class="k">Reads produced</span><span class="v">${fmt(state.reads.length)}</span></div>
    <div class="stat"><span class="k">Reads with errors</span><span class="v">${fmt(strandsHit)}</span></div>
    <div class="stat"><span class="k">Base mismatches</span><span class="v">${fmt(basesHit)} / ${fmt(totalReadBases)}</span></div>
    <div class="stat"><span class="k">Dropped reads</span><span class="v">${fmt(Math.max(0, dropped))}</span></div>`;

  renderMutatedPreview(byIdx);
  $('recoverBtn').disabled = false;
  renderImagePanels(); // updates the MUTATED-raw panel
}

function renderMutatedPreview(byIdx) {
  const reads = state.reads;
  const n = Math.min(6, reads.length);
  let html = `<div class="legend"><span class="muted">corrupted bases highlighted (vs original strand)</span></div>`;
  for (let i = 0; i < n; i++) {
    const r = reads[i];
    const orig = byIdx.get(r.index) || '';
    html += `<div class="strand-row">
      <span class="strand-idx">#${r.index}</span>
      <span class="strand-seq">${renderSeq(r.sequence, orig)}</span>
    </div>`;
  }
  $('mutatedPreview').innerHTML = html;
}

// -----------------------------------------------------------------------------
// Recover
// -----------------------------------------------------------------------------
function doRecover() {
  const rep = state.codec.decode(state.archive, state.reads);
  state.report = rep;

  const banner = $('recoverBanner');
  banner.className = 'banner ' + (rep.crcOk ? 'ok' : 'fail');
  banner.textContent = rep.crcOk
    ? '✓ RECOVERED — byte-perfect, CRC32 matches the original'
    : '✗ RECOVERY INCOMPLETE — too much damage for the current redundancy';

  const crcRecovered = crc32(rep.bytes);
  $('recoverInfo').innerHTML = `
    <div class="stat"><span class="k">Strands recovered</span><span class="v">${fmt(rep.strandsRecovered)}</span></div>
    <div class="stat"><span class="k">Blocks recovered</span><span class="v">${fmt(rep.blocksRecovered)} / ${fmt(rep.numBlocks)}</span></div>
    <div class="stat"><span class="k">RS errors corrected</span><span class="v">${fmt(rep.errorsCorrected)} symbols</span></div>
    <div class="stat"><span class="k">Reads consumed</span><span class="v">${fmt(rep.readsIn)}</span></div>
    <div class="stat"><span class="k">Length match</span><span class="v">${rep.bytes.length === state.fileBytes.length ? 'yes' : 'no'}</span></div>
    <div class="stat"><span class="k">CRC32 match</span><span class="v">${
      '0x' + crcRecovered.toString(16).padStart(8, '0').toUpperCase()} ${rep.crcOk ? '✓' : '✗ (0x' + state.archive.crc32.toString(16).padStart(8,'0').toUpperCase() + ')'}</span></div>`;

  renderImagePanels();
}

// -----------------------------------------------------------------------------
// Image panels: BEFORE / MUTATED-raw (no ECC) / AFTER (error-corrected)
// -----------------------------------------------------------------------------
function renderImagePanels() {
  const wrap = $('imagePanels');
  if (!state.isImage || !state.fileBytes) { wrap.classList.add('hidden'); return; }
  wrap.classList.remove('hidden');

  // BEFORE — the original file bytes
  setImage('imgBefore', state.fileBytes, state.mimeType, 'Original');

  // MUTATED-raw — decode the corrupted DNA directly with NO error correction.
  if (state.reads) {
    const raw = decodeNoEcc();
    setImage('imgMutated', raw, state.mimeType, 'No error correction');
  } else {
    clearImage('imgMutated', 'mutate first');
  }

  // AFTER — error-corrected recovery
  if (state.report) {
    setImage('imgAfter', state.report.bytes, state.mimeType,
      state.report.crcOk ? 'Recovered (clean)' : 'Recovered (partial)');
  } else {
    clearImage('imgAfter', 'recover first');
  }
}

// Decode corrupted reads with NO RS / consensus — pure base-3 inverse of the
// first surviving read per logical block. Shows the visible glitches.
function decodeNoEcc() {
  const a = state.archive;
  const { blockSize, rsParity } = a.params;
  const frameLen = 2 + blockSize;
  const codewordBases = (frameLen + rsParity) * 6;
  const out = new Uint8Array(a.numBlocks * blockSize);
  const seen = new Set();
  // group by physical index, take ONE read each, strip header+parity blindly
  const firstRead = new Map();
  for (const r of state.reads) if (!firstRead.has(r.index)) firstRead.set(r.index, r.sequence);

  for (const [idx, seqRaw] of firstRead) {
    let seq = seqRaw;
    if (seq.length > codewordBases) seq = seq.slice(0, codewordBases);
    else seq = seq + 'A'.repeat(Math.max(0, codewordBases - seq.length));
    let bytes;
    try { bytes = decodeBytes(seq, true); } catch { continue; }
    const blk = (bytes[0] << 8) | bytes[1];
    if (blk < 0 || blk >= a.numBlocks || seen.has(blk)) continue;
    seen.add(blk);
    const payload = bytes.slice(2, 2 + blockSize);
    out.set(payload, blk * blockSize);
  }
  return out.slice(0, a.byteLen);
}

function setImage(id, bytes, mime, caption) {
  const fig = $(id);
  const img = fig.querySelector('img');
  const cap = fig.querySelector('figcaption .sub');
  try {
    const blob = new Blob([bytes], { type: mime || 'image/png' });
    const url = URL.createObjectURL(blob);
    if (img.dataset.url) URL.revokeObjectURL(img.dataset.url);
    img.dataset.url = url;
    img.src = url;
    img.style.visibility = 'visible';
  } catch {
    img.style.visibility = 'hidden';
  }
  if (cap) cap.textContent = caption || '';
}
function clearImage(id, caption) {
  const fig = $(id);
  const img = fig.querySelector('img');
  img.removeAttribute('src');
  img.style.visibility = 'hidden';
  const cap = fig.querySelector('figcaption .sub');
  if (cap) cap.textContent = caption || '';
}

// -----------------------------------------------------------------------------
// Lightweight DNA double-helix animation on a canvas (decorative)
// -----------------------------------------------------------------------------
let animState = { seed: 1, raf: 0 };
function drawStrandAnimation(bytes) {
  const canvas = $('helix');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  // derive a base sequence from the file (or a default) to "flow" through
  let seq = bytes && bytes.length ? encodeBytes(bytes.slice(0, 64)) : encodeBytes(new TextEncoder().encode('DNA STORAGE'));
  if (seq.length < 80) seq = (seq + seq + seq).slice(0, 80);

  cancelAnimationFrame(animState.raf);
  let t = 0;
  function frame() {
    ctx.clearRect(0, 0, W, H);
    const n = 40;
    const amp = H * 0.32, midY = H / 2;
    for (let i = 0; i < n; i++) {
      const x = (i / (n - 1)) * (W - 20) + 10;
      const phase = (i * 0.5) - t;
      const y1 = midY + Math.sin(phase) * amp;
      const y2 = midY - Math.sin(phase) * amp;
      // rungs
      ctx.strokeStyle = 'rgba(120,140,180,0.25)';
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x, y1); ctx.lineTo(x, y2); ctx.stroke();
      // backbone nodes colored by the flowing base sequence
      const b = seq[(i + Math.floor(t * 4)) % seq.length];
      ctx.fillStyle = BASE_COLOR[b] || '#88a';
      const r = 3 + Math.cos(phase) * 1.5;
      ctx.beginPath(); ctx.arc(x, y1, Math.abs(r), 0, 7); ctx.fill();
      ctx.fillStyle = BASE_COLOR[BASES[(BASES.indexOf(b) + 2) % 4]];
      ctx.beginPath(); ctx.arc(x, y2, Math.abs(r), 0, 7); ctx.fill();
    }
    t += 0.04;
    animState.raf = requestAnimationFrame(frame);
  }
  frame();
}

// -----------------------------------------------------------------------------
// Slider value labels
// -----------------------------------------------------------------------------
function bindSlider(id, fmtFn) {
  const el = $(id), out = $(id + 'Val');
  const upd = () => { out.textContent = fmtFn(el.value); };
  el.addEventListener('input', upd); upd();
}

// -----------------------------------------------------------------------------
// Wire up
// -----------------------------------------------------------------------------
function init() {
  // self-test in the console so the demo is auditable; surfaces any regression.
  try {
    const log = selfTest();
    console.log('%cDNA codec self-test: ALL PASSED', 'color:#27d6a4;font-weight:bold');
    console.log(log.join('\n'));
  } catch (e) {
    console.error('DNA codec self-test FAILED:', e);
  }

  $('fileInput').addEventListener('change', (e) => {
    if (e.target.files[0]) readFile(e.target.files[0]);
  });
  $('sampleBtn').addEventListener('click', loadSamplePng);
  $('encodeBtn').addEventListener('click', doEncode);
  $('mutateBtn').addEventListener('click', doMutate);
  $('recoverBtn').addEventListener('click', doRecover);

  // drag & drop
  const dz = $('dropzone');
  ['dragover', 'dragenter'].forEach(ev => dz.addEventListener(ev, (e) => {
    e.preventDefault(); dz.classList.add('drag');
  }));
  ['dragleave', 'drop'].forEach(ev => dz.addEventListener(ev, (e) => {
    e.preventDefault(); dz.classList.remove('drag');
  }));
  dz.addEventListener('drop', (e) => {
    if (e.dataTransfer.files[0]) readFile(e.dataTransfer.files[0]);
  });
  dz.addEventListener('click', () => $('fileInput').click());

  bindSlider('pSub', v => v + '%');
  bindSlider('pIns', v => v + '%');
  bindSlider('pDel', v => v + '%');
  bindSlider('pDrop', v => v + '%');
  bindSlider('coverage', v => v + '×');
  bindSlider('blockSize', v => v + ' B');
  bindSlider('rsParity', v => v + ' B');
  bindSlider('overhead', v => v + '×');

  drawStrandAnimation(null);
}

window.addEventListener('DOMContentLoaded', init);
