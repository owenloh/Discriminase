// Build a commensal index in the browser, one genome at a time.
// Storage-agnostic: the caller supplies `fetchSeq`; IndexedDB caching lives in
// store.mjs and is layered on top by the UI.

import { extractGuides, extractTargetGuides, gcFraction, unpackGuide, GuideIndex } from "./engine.mjs";

// Fetch a source's sequence, retrying transient failures with backoff. NCBI throttles
// hard (especially without an email), so a single attempt drops genomes silently and a
// missing commensal means its cut-sites go UNPROTECTED. Retrying keeps the index complete.
async function fetchSeqWithRetry(fetchSeq, src, { attempts = 4, onProgress, i, n } = {}) {
  let err;
  for (let a = 0; a < attempts; a++) {
    try { return await fetchSeq(src); }
    catch (e) {
      err = e;
      if (a < attempts - 1) {
        onProgress?.({ i, n, name: src.name, status: "retrying", attempt: a + 1, error: String(e) });
        await new Promise((r) => setTimeout(r, 500 * 2 ** a));   // 0.5s, 1s, 2s
      }
    }
  }
  throw err;
}

// params: {guideLength, seedLen, pam, side}
// sources: [{name, seq?} | {name, accession?} | {name, taxid?}]
// fetchSeq(src) -> {seq}   (only called when src has no inline seq)
export async function buildIndex(sources, params, { fetchSeq, onProgress } = {}) {
  const guides = [];
  const orgIds = [];
  const organisms = [];
  const failed = [];

  for (let i = 0; i < sources.length; i++) {
    const src = sources[i];
    onProgress?.({ i, n: sources.length, name: src.name, status: "fetching" });
    let seq;
    try {
      seq = src.seq ?? (await fetchSeqWithRetry(fetchSeq, src, { onProgress, i, n: sources.length })).seq;
      if (!seq || seq.length < params.guideLength) throw new Error("empty/short sequence");
    } catch (e) {
      failed.push({ name: src.name, accession: src.accession, taxid: src.taxid, error: String(e) });
      onProgress?.({ i, n: sources.length, name: src.name, status: "failed", error: String(e) });
      continue;
    }
    const g = extractGuides(seq, params.guideLength, params.pam, params.side);
    const orgIndex = organisms.length;
    organisms.push({
      name: src.name, accession: src.accession || null, taxid: src.taxid || null,
      length_bp: seq.length, n_guides: g.length,
    });
    for (const v of g) { guides.push(v); orgIds.push(orgIndex); }
    onProgress?.({ i, n: sources.length, name: src.name, status: "done", guides: g.length });
  }

  if (!guides.length) throw new Error("no genomes ingested; cannot build index");
  const index = GuideIndex.fromPacked(
    Float64Array.from(guides), params.guideLength, params.seedLen,
    Uint16Array.from(orgIds), organisms);
  return { index, failed };
}

// Collect the target guides shared by EVERY target genome (exact-match intersection),
// the way the commensal panel is built one genome at a time. With one target this is
// just that genome's guides; with several it's the guides that cut all of them.
//
// sources: same shape as buildIndex's. params: {guideLength, pam, side}.
// region:  optional {start, end} applied ONLY when there is a single target (its
//          coordinates are meaningless across multiple genomes).
// Returns {packed: Float64Array, starts, strands, total, organisms, failed}. Position
// provenance (start/strand) is carried from the FIRST target the guide appears in.
export async function collectCommonTargetGuides(sources, params, { fetchSeq, onProgress, region } = {}) {
  let common = null;                          // Map<value, {start, strand}> running intersection
  const organisms = [];
  const failed = [];
  for (let i = 0; i < sources.length; i++) {
    const src = sources[i];
    onProgress?.({ i, n: sources.length, name: src.name, status: "fetching" });
    let seq;
    try {
      seq = src.seq ?? (await fetchSeqWithRetry(fetchSeq, src, { onProgress, i, n: sources.length })).seq;
      if (!seq || seq.length < params.guideLength) throw new Error("empty/short sequence");
    } catch (e) {
      failed.push({ name: src.name, accession: src.accession, taxid: src.taxid, error: String(e) });
      onProgress?.({ i, n: sources.length, name: src.name, status: "failed", error: String(e) });
      continue;
    }
    let sub = seq, off = 0;
    if (region && sources.length === 1) {     // region only makes sense for one genome
      sub = seq.slice(region.start, region.end); off = region.start;
    }
    const { packed, starts, strands } = extractTargetGuides(sub, params.guideLength, params.pam, params.side);
    organisms.push({ name: src.name, accession: src.accession || null, taxid: src.taxid || null,
                     length_bp: seq.length, n_guides: packed.length });
    if (common === null) {                    // first good target seeds the set (+ provenance)
      common = new Map();
      for (let k = 0; k < packed.length; k++)
        if (!common.has(packed[k])) common.set(packed[k], { start: starts[k] + off, strand: strands[k] });
    } else {                                  // later targets prune to the shared subset
      const here = new Set(packed);
      for (const v of common.keys()) if (!here.has(v)) common.delete(v);
    }
    onProgress?.({ i, n: sources.length, name: src.name, status: "done", guides: packed.length, common: common.size });
    if (common.size === 0) break;             // nothing shared -> no point fetching the rest
  }
  if (common === null) throw new Error("no target genomes ingested; add at least one target");
  const packed = [], starts = [], strands = [];
  for (const [v, prov] of common) { packed.push(v); starts.push(prov.start); strands.push(prov.strand); }
  return { packed: Float64Array.from(packed), starts, strands, total: packed.length, organisms, failed };
}

// Screen a precomputed target-guide set against a commensal index.
//   target    = {packed, starts, strands} from collectCommonTargetGuides
//   total     = how many target guides there are (the universe)
//   scanned   = how many we actually screened (< total if the cap hit or cancelled)
//   rows      = the commensal-sparing keepers (partial if cancelled)
//   cancelled = true if `shouldStop()` asked us to bail mid-run
// Async + time-sliced: every ~`tickMs` of work it reports progress and yields to the
// event loop so the browser can repaint the bar/ETA (and a Stop click can land).
export async function screenTargetGuides(target, index, params,
  { onProgress, shouldStop, maxGuides = 1000, tickMs = 80 } = {}) {
  const { packed, starts, strands } = target;
  const total = packed.length;
  const flip = params.side === "3prime";
  const d = params.totalMm, s = params.seedMm;
  const cap = maxGuides || total;
  const rows = [];
  const now = () => (typeof performance !== "undefined" ? performance.now() : Date.now());
  let lastTick = now();
  let scanned = total;                       // a full pass, unless the cap/cancel breaks us out early
  let cancelled = false;
  for (let i = 0; i < total; i++) {
    if (onProgress && now() - lastTick >= tickMs) {
      onProgress(i, total, rows.length);
      await new Promise((r) => setTimeout(r, 0));   // let the browser paint + handle a Stop click
      if (shouldStop && shouldStop()) { scanned = i; cancelled = true; break; }
      lastTick = now();
    }
    const g = packed[i];
    const frac = gcFraction(g, params.guideLength);
    if (frac < params.minGc || frac > params.maxGc) continue;
    if (index.query(g, d, s) !== null) continue;
    let seq = unpackGuide(g, params.guideLength);
    if (flip) seq = [...seq].reverse().join("");
    rows.push({ position: starts[i], strand: strands[i] ? "-" : "+", guide_sequence: seq, gc: Math.round(frac * 1000) / 1000 });
    if (rows.length >= cap) { scanned = i + 1; break; }
  }
  onProgress?.(scanned, total, rows.length);
  return { rows, scanned, total, cancelled };
}
