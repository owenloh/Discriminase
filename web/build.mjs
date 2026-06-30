// Build a commensal index in the browser, one genome at a time.
// Storage-agnostic: the caller supplies `fetchSeq`; IndexedDB caching lives in
// store.mjs and is layered on top by the UI.

import { extractGuides, extractTargetGuides, gcFraction, unpackGuide, GuideIndex } from "./engine.mjs";

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
      seq = src.seq ?? (await fetchSeq(src)).seq;
      if (!seq || seq.length < params.guideLength) throw new Error("empty/short sequence");
    } catch (e) {
      failed.push({ ...src, error: String(e) });
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

// Screen a target sequence against an index -> array of result rows.
// Async + time-sliced: every ~`tickMs` of work it reports progress and yields to
// the event loop so the browser can repaint the progress bar / ETA. Without the
// yield the whole loop runs in one synchronous task and the bar only "appears"
// once everything is already done.
export async function findSparingGuides(targetSeq, index, params,
  { onProgress, maxGuides = 1000, positionOffset = 0, tickMs = 80 } = {}) {
  const { packed, starts, strands } = extractTargetGuides(
    targetSeq, params.guideLength, params.pam, params.side);
  const flip = params.side === "3prime";
  const d = params.totalMm, s = params.seedMm;
  const cap = maxGuides || packed.length;
  const rows = [];
  const now = () => (typeof performance !== "undefined" ? performance.now() : Date.now());
  let lastTick = now();
  for (let i = 0; i < packed.length; i++) {
    if (onProgress && now() - lastTick >= tickMs) {
      onProgress(i, packed.length, rows.length);
      await new Promise((r) => setTimeout(r, 0));   // let the browser paint
      lastTick = now();
    }
    const g = packed[i];
    const frac = gcFraction(g, params.guideLength);
    if (frac < params.minGc || frac > params.maxGc) continue;
    if (index.query(g, d, s) !== null) continue;
    let seq = unpackGuide(g, params.guideLength);
    if (flip) seq = [...seq].reverse().join("");
    rows.push({ position: starts[i] + positionOffset, strand: strands[i] ? "-" : "+", guide_sequence: seq, gc: Math.round(frac * 1000) / 1000 });
    if (rows.length >= cap) break;
  }
  onProgress?.(packed.length, packed.length, rows.length);
  return rows;
}
