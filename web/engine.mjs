// Discriminase engine, browser edition. A faithful port of the Python core
// (pack.py / genome.py / index.py). It runs entirely client-side.
//
// Key trick: a guide of <= 26 letters packs into an integer < 2^52, which a JS
// Number (float64) represents EXACTLY. So packing is `v = v*4 + code` (not bit
// shifts, which JS caps at 32 bits), sorting/binary-search use plain Numbers, and
// only the tiny per-candidate mismatch count uses BigInt. Guides are stored
// PAM-proximal (seed) first, so sorting groups guides by seed -> binary search.

const CODE = new Int8Array(256);           // ASCII -> 2-bit code, default 0 (A)
CODE[67] = 1; CODE[71] = 2; CODE[84] = 3;  // C G T
CODE[99] = 1; CODE[103] = 2; CODE[116] = 3; // c g t
const BASE = "ACGT";

const IUPAC = {
  A: [0], C: [1], G: [2], T: [3],
  R: [0, 2], Y: [1, 3], S: [1, 2], W: [0, 3], K: [2, 3], M: [0, 1],
  B: [1, 2, 3], D: [0, 2, 3], H: [0, 1, 3], V: [0, 1, 2], N: [0, 1, 2, 3],
};

// --- packing ---------------------------------------------------------------
export function packGuide(seq) {
  let v = 0;
  for (let i = 0; i < seq.length; i++) v = v * 4 + CODE[seq.charCodeAt(i) & 0xff];
  return v;
}
export function unpackGuide(value, L) {
  const out = new Array(L);
  for (let i = L - 1; i >= 0; i--) { out[i] = BASE[value % 4]; value = Math.floor(value / 4); }
  return out.join("");
}
export function seedOf(value, L, seedLen) {
  return Math.floor(value / 4 ** (L - seedLen));
}

// --- Hamming over packed guides (BigInt; used only on small candidate blocks)
const FOLD = 0x5555555555555555n;
function popcountBig(x) { let c = 0; while (x) { x &= x - 1n; c++; } return c; }
export function hamming(a, b) {
  const x = BigInt(a) ^ BigInt(b);
  return popcountBig((x | (x >> 1n)) & FOLD);
}

export function seedNeighbors(seed, seedLen, maxMm) {
  if (maxMm <= 0) return [seed];
  const results = new Set([seed]);
  let frontier = [seed];
  for (let step = 0; step < maxMm; step++) {
    const next = [];
    for (const s of frontier) {
      for (let pos = 0; pos < seedLen; pos++) {
        const place = 4 ** pos;
        const cur = Math.floor(s / place) % 4;
        for (let c = 0; c < 4; c++) {
          if (c === cur) continue;
          const nb = s + (c - cur) * place;
          if (!results.has(nb)) { results.add(nb); next.push(nb); }
        }
      }
    }
    frontier = next;
  }
  return [...results];
}

// --- genome encoding + PAM extraction --------------------------------------
export function encode(seq) {
  const a = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) a[i] = CODE[seq.charCodeAt(i) & 0xff];
  return a;
}
export function reverseComplement(codes) {
  const n = codes.length, out = new Uint8Array(n);
  for (let i = 0; i < n; i++) out[i] = 3 - codes[n - 1 - i];
  return out;
}
function pamSets(pam) {
  return [...pam.toUpperCase()].map((b) => {
    if (!(b in IUPAC)) throw new Error(`invalid PAM letter ${b}`);
    const a = IUPAC[b];
    return a.length === 4 ? null : a;       // null == any base
  });
}
function pamPositions(codes, sets) {
  const n = codes.length, p = sets.length, out = [];
  for (let i = 0; i + p <= n; i++) {
    let ok = true;
    for (let j = 0; j < p && ok; j++) {
      const allowed = sets[j];
      if (allowed === null) continue;
      if (!allowed.includes(codes[i + j])) ok = false;
    }
    if (ok) out.push(i);
  }
  return out;
}
function packWindow(codes, start, L, reverse) {
  let v = 0;
  if (!reverse) for (let j = 0; j < L; j++) v = v * 4 + codes[start + j];
  else for (let j = L - 1; j >= 0; j--) v = v * 4 + codes[start + j];
  return v;
}
// A/C/G/T validity: anything else (N, IUPAC ambiguity, gaps) is not a definite base.
const IS_ACGT = new Uint8Array(256);
for (const c of [65, 67, 71, 84, 97, 99, 103, 116]) IS_ACGT[c] = 1; // ACGT acgt
function invArray(seq) {
  const n = seq.length, a = new Uint8Array(n);
  for (let i = 0; i < n; i++) a[i] = IS_ACGT[seq.charCodeAt(i) & 0xff] ? 0 : 1;
  return a;
}
function prefixSum(inv) {
  const cum = new Int32Array(inv.length + 1);
  for (let i = 0; i < inv.length; i++) cum[i + 1] = cum[i] + inv[i];
  return cum;
}

function strandGuides(codes, cumInv, L, sets, gap, side) {
  const n = codes.length, pamlen = sets.length, starts = [], packed = [];
  for (const p of pamPositions(codes, sets)) {
    let start, reverse;
    if (side === "5prime") { start = p + pamlen + gap; reverse = false; }
    else if (side === "3prime") { start = p - gap - L; reverse = true; }
    else throw new Error(`pam_side must be 5prime/3prime`);
    // skip windows that run off the end or contain any ambiguous/invalid base
    if (start >= 0 && start + L <= n && cumInv[start + L] - cumInv[start] === 0) {
      starts.push(start); packed.push(packWindow(codes, start, L, reverse));
    }
  }
  return { starts, packed };
}

export function extractGuides(seq, L, pam = "TTT", gap = 1, side = "5prime") {
  const codes = encode(seq), sets = pamSets(pam), inv = invArray(seq);
  const f = strandGuides(codes, prefixSum(inv), L, sets, gap, side);
  const r = strandGuides(reverseComplement(codes), prefixSum(inv.slice().reverse()), L, sets, gap, side);
  return f.packed.concat(r.packed);
}

// Target extraction keeps provenance (forward position + strand) for the output.
export function extractTargetGuides(seq, L, pam = "TTT", gap = 1, side = "5prime") {
  const codes = encode(seq), n = codes.length, sets = pamSets(pam), inv = invArray(seq);
  const f = strandGuides(codes, prefixSum(inv), L, sets, gap, side);
  const r = strandGuides(reverseComplement(codes), prefixSum(inv.slice().reverse()), L, sets, gap, side);
  const packed = f.packed.concat(r.packed);
  const starts = f.starts.concat(r.starts.map((s) => n - s - L));   // rc -> forward coord
  const strands = f.starts.map(() => 0).concat(r.starts.map(() => 1));
  return { packed, starts, strands };
}

export function gcFraction(value, L) {
  let gc = 0, v = value;
  for (let k = 0; k < L; k++) { const f = v % 4; if (f === 1 || f === 2) gc++; v = Math.floor(v / 4); }
  return gc / L;
}

// --- the index -------------------------------------------------------------
function lowerBound(arr, x) {            // first index with arr[i] >= x
  let lo = 0, hi = arr.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (arr[mid] < x) lo = mid + 1; else hi = mid; }
  return lo;
}

export class GuideIndex {
  constructor(guides, guideLength, seedLen, orgIds = null, organisms = []) {
    this.guides = guides;               // sorted Float64Array
    this.orgIds = orgIds;
    this.guideLength = guideLength;
    this.seedLen = seedLen;
    this.organisms = organisms;
    this.shiftDiv = 4 ** (guideLength - seedLen);
  }

  static fromPacked(values, guideLength, seedLen, orgIds = null, organisms = []) {
    // sort by guide, dedup, keep first org per unique guide
    const idx = Array.from(values.keys()).sort((a, b) => values[a] - values[b]);
    const g = [], o = orgIds ? [] : null;
    let prev = NaN;
    for (const i of idx) {
      if (values[i] !== prev) { g.push(values[i]); if (o) o.push(orgIds[i]); prev = values[i]; }
    }
    return new GuideIndex(Float64Array.from(g), guideLength, seedLen,
      o ? Uint16Array.from(o) : null, organisms);
  }

  get length() { return this.guides.length; }

  // nearest colliding commensal under the seed-anchored model, or null
  query(g, totalMm, seedMm) {
    const seeds = seedNeighbors(seedOf(g, this.guideLength, this.seedLen), this.seedLen, seedMm);
    let best = null;
    for (const s of seeds) {
      const lo = lowerBound(this.guides, s * this.shiftDiv);
      const hi = lowerBound(this.guides, (s + 1) * this.shiftDiv);
      for (let i = lo; i < hi; i++) {
        const d = hamming(g, this.guides[i]);
        if (d <= totalMm && (best === null || d < best.dist)) {
          best = { dist: d, org: this.orgIds ? this.orgIds[i] : -1 };
          if (d === 0) return best;
        }
      }
    }
    return best;
  }

  isSpared(g, totalMm, seedMm) { return this.query(g, totalMm, seedMm) === null; }
}
