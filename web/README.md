# 🧬 DNA Storage Codec Simulator — browser demo

The live demo behind the hook **"I stored my GitHub repo in DNA."**

Turn any file into a pool of synthetic DNA strands ("oligos"), push them through a
noisy sequencing channel (substitutions, insertions, deletions, whole-strand
dropout), then reconstruct the original bytes — demonstrating end-to-end error
correction with **no wet lab, no build step, and no network**.

## How to open it

Just open `index.html` in any modern browser:

```
# from the repo root
xdg-open web/index.html      # Linux
open    web/index.html       # macOS
# or double-click web/index.html in a file manager
```

That's it. No `npm`, no bundler, no server, no CDN. Everything is vanilla
HTML/CSS/JS using ES modules. The sample image is embedded as a base64 PNG so
the **"Use sample PNG"** button works fully offline.

> If your browser blocks `file://` ES-module imports (some do for security),
> serve the folder with any static server, e.g. `python3 -m http.server` inside
> `web/` and open <http://localhost:8000>.

## What you can do

1. **Load a file** — drag-and-drop, file picker, or the built-in sample PNG.
   Any file type works; images get a before/after visualization.
2. **Encode to DNA** — see stats (size, # strands, strand length, total bases,
   bits/base density, mean GC%, max homopolymer run) and a colored ACGT preview
   of the first strands.
3. **Mutate** — apply the sequencing channel with sliders for substitution,
   insertion, deletion, strand-dropout %, and coverage (reads/strand). Corrupted
   bases are highlighted vs. the original.
4. **Recover** — run consensus + Reed–Solomon decode + reassembly. A banner
   reports success/failure with # errors corrected, # strands/blocks recovered,
   and a byte-level CRC32/length match.

For images you get three panels side by side:

| BEFORE | MUTATED — raw | AFTER — corrected |
|--------|---------------|-------------------|
| the original | the corrupted DNA decoded **with no error correction** (visible glitches) | the same data after consensus + RS recovery (clean) |

The contrast between the glitched middle image and the clean recovery is the
whole point.

Open the browser console to see the codec **self-test** run on load (CRC vector,
byte round-trip, homopolymer-free check, GF(256) inverse, RS error correction,
end-to-end recovery). You can also run it from a terminal: `node web/codec.js`.

## Files

| File | Role |
|------|------|
| `index.html` | the page / layout |
| `styles.css` | dark "lab" aesthetic, responsive |
| `codec.js`   | the pure-JS codec (the core) + self-test |
| `app.js`     | UI wiring, rendering, visualization, embedded sample PNG |
| `README.md`  | this file |

## Mapping to the Rust crate (`examples/dna/src/storage/`)

`codec.js` is a faithful, dependency-free JavaScript mirror of the Rust storage
codec. It implements the **same LT fountain + Reed–Solomon + base-3 mapping
pipeline** as the Rust crate — including a bit-for-bit port of the fountain
code's SplitMix64 PRNG and seed→neighbour selection, so JS and Rust agree
exactly on which source blocks each droplet XORs. Each section maps 1:1 to a
Rust submodule:

| Rust (`examples/dna/src/storage/`) | JavaScript (`web/codec.js`) | Notes |
|---|---|---|
| `constraints.rs` — `encode_bytes` / `decode_bytes` / `gc_content` / `max_homopolymer_run` | `encodeBytes` / `decodeBytes` / `gcContent` / `maxHomopolymerRun` | Identical homopolymer-free base-3 transform: `A=0,C=1,G=2,T=3`; `prev` init 0; each byte → 6 big-endian ternary digits `d[i]=⌊b/3^(5-i)⌋%3`; per trit `cur=(prev+1+t)%4`; inverse `t=(cur+4-prev-1)%4`. Guarantees max homopolymer run = 1. |
| `gf256.rs` — `ReedSolomon { nsym }` | `ReedSolomon` class | GF(256) with primitive polynomial `0x11D`, generator `2`. Systematic encoder appends `nsym` parity bytes; decoder is the classic syndrome → Berlekamp–Massey → Chien search → Forney pipeline and corrects up to `nsym/2` substitutions. Default `nsym = 8` (corrects 4). |
| `channel.rs` — `ErrorModel`, `apply` | `DEFAULT_MODEL`, `applyChannel`, `Rng` | Per-base substitution / insertion / deletion, whole-strand dropout, and `coverage` reads per surviving strand. Deterministic via a seeded `mulberry32` PRNG; returns a shuffled read pool. |
| `consensus.rs` — `cluster`, `consensus` | `cluster`, `consensus` + the salvage pass in `DnaStorageCodec.decode` | Reads that fail per-strand RS are greedily clustered by normalized Hamming similarity, then a position-wise majority vote (up to the modal read length) cancels random substitutions when coverage > 1, and RS is retried. |
| `fountain.rs` — Luby-Transform (LT) outer code | `SplitMix64`, `robustSolitonPmf`, `sampleDegree`, `neighbours`, `LtEncoder`, `LtDecoder` | Faithful bit-for-bit port: SplitMix64 PRNG (BigInt 64-bit state, constants `0x9E3779B97F4A7C15` / `0xBF58476D1CE4E5B9` / `0x94D049BB133111EB`), robust-soliton degree distribution (`c=0.03`, `delta=0.05`), deterministic seed→neighbour selection (degree draw + partial Fisher–Yates, sorted distinct indices), and belief-propagation peeling decoder. Encode and decode agree on every droplet's neighbour set because both recompute it from the droplet `seed`. |
| `mod.rs` — `crc32`, `EncodeParams`, `Strand`, `DnaArchive`, `DecodeReport`, `DnaStorageCodec` | `crc32`, `DEFAULT_PARAMS`, archive objects, `DnaStorageCodec`, `totalBases` / `bitsPerBase` / `meanGc` | IEEE CRC32 (reflected, poly `0xEDB88320`) — matches the Rust known-answer test `crc32("123456789") == 0xCBF43926`. Same archive fields and derived stats. |

### Per-strand wire format

Mirrors the Rust pipeline exactly — `header [index|seed] → LT droplet → RS parity → DNA mapping`:

```
[ index(4) | seed(4) | droplet payload(block_size) | RS parity(rs_parity) ]  →  encodeBytes(...)
   u32 BE strand id  u32 BE fountain seed   LT droplet bytes    GF(256) parity      ACGT string
```

The 8-byte header (strand index + fountain seed) lives **inside** the RS
codeword, so the exact `seed` is error-corrected before it is trusted. The
decoder recomputes each droplet's neighbour set from that seed, dedupes droplets
by seed, and peels them with the LT decoder — identical to the Rust contract.

### Honest differences from the Rust crate

The demo now matches the Rust outer code; the remaining simplifications are:

- **Indels are the hard frontier.** Substitutions and whole-strand dropout are
  fully recovered at the default settings. Insertions/deletions *frame-shift* a
  strand and desynchronize the differential base-3 decode, so the demo ships
  with the indel sliders at 0 by default and lets you turn them up to *watch*
  recovery degrade — an honest depiction of the real open problem in DNA storage.
- **Determinism.** The channel is seeded for reproducible-ish runs in the
  browser; the Rust `apply` takes an explicit `seed`.

Everything else — the LT fountain outer code, the base-3 constraint mapping,
GF(256) Reed–Solomon, majority-vote consensus, CRC32 integrity check, and the
archive/stat definitions — matches the Rust contracts directly.
