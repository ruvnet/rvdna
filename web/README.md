# 🧬 DNA Storage Codec — Browser Visualizer

A dependency-free, single-page demo of the [rvDNA DNA Storage
Codec](../examples/dna/src/storage). Drop in a file (or use the built-in sample
PNG), **encode it to synthetic DNA**, **mutate** it with a noisy sequencing
channel, and watch **error correction recover the original — live**.

> Viral hook: *"I stored my GitHub repo in DNA."*

## Run it

No build step, no npm, no network. Just open the page:

```bash
# from the repo root
open web/index.html            # macOS
xdg-open web/index.html        # Linux
# …or serve it:
python3 -m http.server -d web 8000   # then visit http://localhost:8000
```

(The files are ES modules, so loading over `file://` works in most browsers; if
yours blocks module loading from disk, use the `http.server` option above.)

## What you can do

1. **Load** any file by drag-and-drop, or click **Use sample PNG**.
2. **Encode to DNA** — see strand count, strand length, total bases, density
   (bits/base), mean **GC%**, and max homopolymer run, plus a colored ACGT
   preview of the strands.
3. **Mutate** — sliders for substitution / insertion / deletion / strand-dropout
   rates and sequencing **coverage**. Corrupted bases are highlighted.
4. **Recover** — runs consensus + Reed–Solomon decode + fountain peeling and
   reports strands recovered, errors corrected, and a CRC32 / length match. For
   images you get a **before / mutated-raw / error-corrected** comparison so the
   glitch-vs-clean contrast is obvious.

## How it maps to the Rust codec

`codec.js` faithfully re-implements the same algorithms as the authoritative
Rust crate so the demo is honest, not a mock:

| `web/codec.js` | Rust `examples/dna/src/storage/` |
|----------------|----------------------------------|
| `encodeBytes` / `decodeBytes` (base-3, homopolymer-free) | `constraints.rs` |
| `ReedSolomon` (GF(256), prim `0x11D`, α=2) | `gf256.rs` |
| `applyChannel` (sub/ins/del/dropout + coverage) | `channel.rs` |
| `consensus` (majority vote) | `consensus.rs` |
| `DnaStorageCodec`, `crc32` | `mod.rs` |

A fountain/LT layer is modeled in the Rust codec (`fountain.rs`); the browser
demo focuses on the per-strand Reed–Solomon + consensus recovery that is most
visually instructive. The Rust crate is `wasm32`-ready, so a future build can
swap this JS core for the compiled WASM with no UI change.

`codec.js` exposes a `selfTest()` that checks the CRC32 vector, the
homopolymer-free guarantee, full byte round-trips, GF(256) inverses, RS error
correction, and end-to-end recovery under noise.

## Files

- `index.html` — page layout and controls
- `styles.css` — dark "lab" styling
- `codec.js` — the pure-JS codec (constraints, Reed–Solomon, channel, consensus)
- `app.js` — UI wiring, rendering, and the DNA visualization
