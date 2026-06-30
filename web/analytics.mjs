// Optional, privacy-bounded usage analytics via self-hosted Umami.
// Dormant until web/config.js is filled in (after Umami is deployed) -- every call is
// a no-op when unconfigured, so the app works identically with or without it.
//
// What we send: organism names/accessions, parameters, result counts, and SHORT pasted
// DNA. What we never send: uploaded FASTA contents or any large sequence (only its
// length + filename). Keep it that way.

const cfg = (typeof window !== "undefined" && window.DISC_ANALYTICS) || {};
const MAX_SEQ = 120;   // pasted DNA at/under this many bases may be logged verbatim

export function initAnalytics() {
  if (!cfg.websiteId || !cfg.src) return;            // not configured -> off
  const s = document.createElement("script");
  s.async = true;
  s.src = cfg.src;
  s.setAttribute("data-website-id", cfg.websiteId);
  if (cfg.host) s.setAttribute("data-host-url", cfg.host);
  s.setAttribute("data-auto-track", "true");
  document.head.appendChild(s);
}

export function track(event, data = {}) {
  try { if (typeof window !== "undefined" && window.umami) window.umami.track(event, data); }
  catch (e) { /* never let analytics break the app */ }
}

// Summarize a chosen target for logging without ever including a big sequence.
export function targetSummary(info, seq) {
  const base = { type: info.type, length: seq ? seq.length : undefined };
  if (info.type === "accession") base.accession = info.value;
  else if (info.type === "name") { base.accession = info.value; base.name = info.name; }
  else if (info.type === "upload") base.file = info.name;       // filename only, not contents
  else if (info.type === "paste") { if (seq && seq.length <= MAX_SEQ) base.seq = seq; }
  return base;
}
