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

// Flat, primitive-only target fields (Umami stores primitives cleanly; nested
// objects/arrays do not). Never includes a big sequence -- uploads log filename+length.
export function targetSummary(info, seq) {
  const f = { target_type: info.type, target_len: seq ? seq.length : undefined };
  if (info.type === "accession") f.target = info.value;
  else if (info.type === "name") { f.target = info.value; f.target_name = info.name; }
  else if (info.type === "upload") f.target_file = info.name;     // filename only, not contents
  else if (info.type === "paste" && seq && seq.length <= MAX_SEQ) f.target_seq = seq;
  return f;
}
