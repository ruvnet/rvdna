const { platform, arch } = process;

// Platform-specific native binary packages
const platformMap = {
  'linux': {
    'x64': '@ruvector/rvdna-linux-x64-gnu',
    'arm64': '@ruvector/rvdna-linux-arm64-gnu'
  },
  'darwin': {
    'x64': '@ruvector/rvdna-darwin-x64',
    'arm64': '@ruvector/rvdna-darwin-arm64'
  },
  'win32': {
    'x64': '@ruvector/rvdna-win32-x64-msvc'
  }
};

function loadNativeModule() {
  const platformPackage = platformMap[platform]?.[arch];

  if (!platformPackage) {
    throw new Error(
      `Unsupported platform: ${platform}-${arch}\n` +
      `@ruvector/rvdna native bindings are available for:\n` +
      `- Linux (x64, ARM64)\n` +
      `- macOS (x64, ARM64)\n` +
      `- Windows (x64)\n\n` +
      `For other platforms, use the WASM build: npm install @ruvector/rvdna-wasm`
    );
  }

  try {
    return require(platformPackage);
  } catch (error) {
    if (error.code === 'MODULE_NOT_FOUND') {
      throw new Error(
        `Native module not found for ${platform}-${arch}\n` +
        `Please install: npm install ${platformPackage}\n` +
        `Or reinstall @ruvector/rvdna to get optional dependencies`
      );
    }
    throw error;
  }
}

// Try native first, fall back to pure JS shim with basic functionality
let nativeModule;
try {
  nativeModule = loadNativeModule();
} catch (e) {
  // Native bindings not available — provide JS shim for basic operations
  nativeModule = null;
}

// -------------------------------------------------------------------
// Public API — wraps native bindings or provides JS fallbacks
// -------------------------------------------------------------------

/**
 * Encode a DNA string to 2-bit packed bytes (4 bases per byte).
 * A=00, C=01, G=10, T=11. Returns a Buffer.
 */
function encode2bit(sequence) {
  if (nativeModule?.encode2bit) return nativeModule.encode2bit(sequence);

  // JS fallback
  const map = { A: 0, C: 1, G: 2, T: 3, N: 0 };
  const len = sequence.length;
  const buf = Buffer.alloc(Math.ceil(len / 4));
  for (let i = 0; i < len; i++) {
    const byteIdx = i >> 2;
    const bitOff = 6 - (i & 3) * 2;
    buf[byteIdx] |= (map[sequence[i]] || 0) << bitOff;
  }
  return buf;
}

/**
 * Decode 2-bit packed bytes back to a DNA string.
 */
function decode2bit(buffer, length) {
  if (nativeModule?.decode2bit) return nativeModule.decode2bit(buffer, length);

  const bases = ['A', 'C', 'G', 'T'];
  let result = '';
  for (let i = 0; i < length; i++) {
    const byteIdx = i >> 2;
    const bitOff = 6 - (i & 3) * 2;
    result += bases[(buffer[byteIdx] >> bitOff) & 3];
  }
  return result;
}

/**
 * Translate a DNA string to a protein amino acid string.
 */
function translateDna(sequence) {
  if (nativeModule?.translateDna) return nativeModule.translateDna(sequence);

  // JS fallback — standard genetic code
  const codons = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
  };
  let protein = '';
  for (let i = 0; i + 2 < sequence.length; i += 3) {
    const codon = sequence.slice(i, i + 3).toUpperCase();
    const aa = codons[codon] || 'X';
    if (aa === '*') break;
    protein += aa;
  }
  return protein;
}

/**
 * Compute cosine similarity between two numeric arrays.
 */
function cosineSimilarity(a, b) {
  if (nativeModule?.cosineSimilarity) return nativeModule.cosineSimilarity(a, b);

  let dot = 0, magA = 0, magB = 0;
  for (let i = 0; i < a.length; i++) {
    dot += a[i] * b[i];
    magA += a[i] * a[i];
    magB += b[i] * b[i];
  }
  magA = Math.sqrt(magA);
  magB = Math.sqrt(magB);
  return (magA && magB) ? dot / (magA * magB) : 0;
}

/**
 * Convert a FASTA sequence string to .rvdna binary format.
 * Returns a Buffer with the complete .rvdna file contents.
 */
function fastaToRvdna(sequence, options = {}) {
  if (nativeModule?.fastaToRvdna) {
    return nativeModule.fastaToRvdna(sequence, options.k || 11, options.dims || 512, options.blockSize || 500);
  }
  throw new Error('fastaToRvdna requires native bindings. Install the platform-specific package.');
}

/**
 * Read a .rvdna file from a Buffer. Returns parsed sections.
 */
function readRvdna(buffer) {
  if (nativeModule?.readRvdna) return nativeModule.readRvdna(buffer);
  throw new Error('readRvdna requires native bindings. Install the platform-specific package.');
}

/**
 * Check if native bindings are available.
 */
function isNativeAvailable() {
  return nativeModule !== null;
}

module.exports = {
  encode2bit,
  decode2bit,
  translateDna,
  cosineSimilarity,
  fastaToRvdna,
  readRvdna,
  isNativeAvailable,

  // Re-export native module for advanced use
  native: nativeModule,
};
