// UI controller. Holds no algorithm logic itself — it wires the DOM to the
// verified engine (engine.mjs / build.mjs / ncbi.mjs) and persistence (store.mjs).

import { GuideIndex } from "./engine.mjs";
import { buildIndex, findSparingGuides } from "./build.mjs";
import * as ncbi from "./ncbi.mjs";
import * as store from "./store.mjs";
import { initAnalytics, track, targetSummary } from "./analytics.mjs";

const $ = (id) => document.getElementById(id);
const PRESETS = {
  cas12a: { pam: "TTTV", side: "5prime", guideLength: 23, seedLen: 10 },
  spcas9: { pam: "NGG", side: "3prime", guideLength: 20, seedLen: 10 },
};
const PARAM_IDS = ["pam", "side", "guideLength", "seedLen", "seedMm",
                   "totalMm", "minGc", "maxGc", "maxGuides", "email"];

const state = { panel: [], index: null, indexMeta: null, target: null, results: null };

// --- parameters ------------------------------------------------------------
function params() {
  return {
    pam: $("pam").value.trim().toUpperCase(), side: $("side").value,
    guideLength: +$("guideLength").value, seedLen: +$("seedLen").value,
    seedMm: +$("seedMm").value, totalMm: +$("totalMm").value,
    minGc: +$("minGc").value, maxGc: +$("maxGc").value, maxGuides: +$("maxGuides").value,
    email: $("email").value.trim(),
  };
}
function persistParams() {
  const v = {}; PARAM_IDS.forEach((id) => (v[id] = $(id).value));
  store.saveSetting("params", v);
}
function restoreParams() {
  const v = store.loadSetting("params", null);
  if (v) PARAM_IDS.forEach((id) => { if (v[id] != null) $(id).value = v[id]; });
}

// --- commensal panel -------------------------------------------------------
function srcLabel(s) { return s.accession ? `acc ${s.accession}` : s.taxid ? `txid ${s.taxid}` : s.uploaded ? "uploaded" : "—"; }
function renderPanel() {
  const tb = $("panel-list");
  tb.innerHTML = "";
  if (!state.panel.length) { tb.innerHTML = `<tr><td colspan="4" class="muted">No commensals yet. Add some above.</td></tr>`; return; }
  state.panel.forEach((s, i) => {
    const tr = document.createElement("tr");
    const status = s.status ? `<span class="pill ${s.status === "failed" ? "bad" : "ok"}">${s.status}${s.guides ? " · " + s.guides : ""}</span>` : "";
    tr.innerHTML = `<td>${s.name}</td><td class="muted">${srcLabel(s)}</td><td>${status}</td><td><button class="ghost" data-i="${i}">✕</button></td>`;
    tr.querySelector("button").onclick = () => { state.panel.splice(i, 1); renderPanel(); savePanelState(); };
    tb.appendChild(tr);
  });
}
function addSource(s) {
  const key = s.accession || s.taxid || ("up:" + s.name);
  if (state.panel.some((x) => (x.accession || x.taxid || ("up:" + x.name)) === key)) return;
  state.panel.push(s); renderPanel(); savePanelState();
}
function savePanelState() {
  // persist metadata only (not uploaded sequences — too big for localStorage)
  store.saveSetting("panel", state.panel.filter((s) => !s.uploaded).map((s) => ({ name: s.name, accession: s.accession, taxid: s.taxid })));
}

async function searchInto(query, container, onPick, single = false) {
  container.style.display = "block";
  container.innerHTML = `<p class="muted" style="padding:8px">Searching…</p>`;
  try {
    const cands = await ncbi.searchCandidates(query, { email: params().email, retmax: 20 });
    if (!cands.length) { container.innerHTML = `<p class="muted" style="padding:8px">No matches.</p>`; return; }
    const t = document.createElement("table");
    t.innerHTML = `<thead><tr><th>accession</th><th>length</th><th>description</th><th></th></tr></thead><tbody></tbody>`;
    const tbody = t.querySelector("tbody");
    cands.forEach((c) => {
      const tr = document.createElement("tr");
      tr.innerHTML = `<td>${c.accession}</td><td>${c.length_bp ? c.length_bp.toLocaleString() + " bp" : "?"}</td><td>${c.title.slice(0, 64)}</td><td><button class="ghost">pick</button></td>`;
      const btn = tr.querySelector("button");
      btn.onclick = () => {
        if (single) tbody.querySelectorAll("tr.picked").forEach((r) => { r.classList.remove("picked"); r.querySelector("button").textContent = "pick"; });
        tr.classList.add("picked");
        btn.textContent = single ? "✓ picked" : "✓ added";
        onPick(c);
      };
      tbody.appendChild(tr);
    });
    container.innerHTML = ""; container.appendChild(t);
  } catch (e) { container.innerHTML = `<p class="bad" style="padding:8px">${e}</p>`; }
}

async function readFasta(file) {
  const text = await file.text();
  return text.split("\n").filter((l) => !l.startsWith(">")).join("").replace(/\s/g, "");
}

// --- build / load index ----------------------------------------------------
function panelKey(p) {
  const geom = { L: p.guideLength, pam: p.pam, side: p.side, seed: p.seedLen };
  const srcs = state.panel.map((s) => s.accession || s.taxid || ("up:" + s.name + ":" + (s.seq ? s.seq.length : 0))).sort();
  return JSON.stringify({ geom, srcs });
}
async function buildOrLoad() {
  if (!state.panel.length) return setStatus("build-status", "Add commensals first.", true);
  const p = params();
  const key = panelKey(p);
  setStatus("build-status", "Checking cache…");
  const cached = await store.loadIndex(key).catch(() => null);
  if (cached) {
    state.index = new GuideIndex(cached.guides, cached.guideLength, cached.seedLen, cached.orgIds, cached.organisms);
    state.indexMeta = cached;
    setStatus("build-status", `Loaded cached index: ${state.index.length.toLocaleString()} guides from ${cached.organisms.length} organisms.`);
    return updateRunnable();
  }
  $("build-btn").disabled = true;
  try {
    const { index, failed } = await buildIndex(state.panel, p, {
      fetchSeq: (s) => s.accession ? ncbi.fetchSequence(s.accession, { email: p.email }) : ncbi.fetchSequenceByTaxid(s.taxid, { email: p.email }),
      onProgress: (ev) => {
        const src = state.panel.find((s) => s.name === ev.name);
        if (src) { src.status = ev.status; src.guides = ev.guides; }
        renderPanel();
        setStatus("build-status", `${ev.status} ${ev.name} (${ev.i + 1}/${ev.n})…`);
      },
    });
    state.index = index; state.indexMeta = { organisms: index.organisms };
    await store.saveIndex(key, { guides: index.guides, orgIds: index.orgIds, guideLength: index.guideLength, seedLen: index.seedLen, organisms: index.organisms }).catch(() => {});
    const note = failed.length ? ` (${failed.length} failed)` : "";
    setStatus("build-status", `Index ready: ${index.length.toLocaleString()} guides from ${index.organisms.length} organisms${note}.`);
    track("build", { preset: $("preset").value, pam: p.pam, side: p.side, guideLength: p.guideLength,
      n_organisms: index.organisms.length, guides: index.length, failed: failed.length,
      organisms: index.organisms.map((o) => o.name).join(", ").slice(0, 480) });  // primitive string for Umami
  } catch (e) {
    setStatus("build-status", String(e), true);
  } finally { $("build-btn").disabled = false; updateRunnable(); }
}

// --- panel sources: examples, prebuilt indexes, saved panels, uploaded CSV --
// Examples + Saved + CSV give you a LIST of organisms that you still Build.
// Prebuilt is already Built (an index) and loads instantly.
function _fillOptgroup(id, items, mkValue, mkLabel) {
  const og = $(id); og.innerHTML = "";
  items.forEach((it) => { const o = document.createElement("option"); o.value = mkValue(it); o.textContent = mkLabel(it); og.appendChild(o); });
  og.hidden = items.length === 0;          // hide empty groups so the menu stays clean
}
async function loadExampleList() {
  try { _fillOptgroup("og-example", await (await fetch("./examples/index.json")).json(),
    (p) => "example:" + p.file, (p) => p.name); } catch (e) { $("og-example").hidden = true; }
}
async function loadPrebuiltList() {
  try { _fillOptgroup("og-prebuilt", await (await fetch("./panels/index.json")).json(),
    (p) => "prebuilt:" + p.prefix, (p) => p.name + " (instant)"); } catch (e) { $("og-prebuilt").hidden = true; }
}
function refreshSavedPanels() {
  _fillOptgroup("og-saved", store.loadSetting("panelNames", []) || [], (n) => "saved:" + n, (n) => n);
}
function saveNamedPanel() {
  if (!state.panel.length) return setStatus("build-status", "Add commensals before saving.", true);
  const name = prompt("Save this panel in your browser as:");
  if (!name) return;
  const names = store.loadSetting("panelNames", []);
  if (!names.includes(name)) names.push(name);
  store.saveSetting("panelNames", names);
  store.saveSetting("panel:" + name, state.panel.filter((s) => !s.uploaded).map((s) => ({ name: s.name, accession: s.accession, taxid: s.taxid })));
  refreshSavedPanels();
  setStatus("build-status", `Saved “${name}” in this browser.`);
}

function loadPanelSources(sources, label) {
  if (!sources.length) return setStatus("build-status", "That panel had no usable organisms.", true);
  state.panel = sources; renderPanel(); savePanelState();
  setStatus("build-status", `Loaded ${sources.length} organisms${label ? " from " + label : ""}. Now press “Build / load index”.`);
}
async function loadExamplePanel(file) {
  setStatus("build-status", "Loading example panel…");
  try { loadPanelSources(parsePanelCsv(await (await fetch("./examples/" + file)).text()), file); }
  catch (e) { setStatus("build-status", String(e), true); }
}
function loadSavedPanel(name) {
  loadPanelSources((store.loadSetting("panel:" + name, []) || []).map((s) => ({ ...s })), name);
}

// CSV parser for uploaded / example panels (columns: taxid / accession / organism_strain)
function _csvLine(line) {
  const out = []; let cur = "", q = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
    else if (c === '"') q = true; else if (c === ",") { out.push(cur); cur = ""; } else cur += c;
  }
  out.push(cur); return out;
}
function parsePanelCsv(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim());
  if (!lines.length) return [];
  const header = _csvLine(lines[0]).map((h) => h.trim().toLowerCase());
  const find = (names) => { for (const n of names) { const i = header.indexOf(n); if (i >= 0) return i; } return -1; };
  let iTax = find(["taxid", "tax_id", "txid"]), iAcc = find(["accession", "acc"]),
      iName = find(["organism_strain", "organism", "name", "strain"]), start = 1;
  if (iTax < 0 && iAcc < 0 && iName < 0) { iTax = 0; iName = 1; start = 0; }   // headerless: taxid,name
  const out = [];
  for (let r = start; r < lines.length; r++) {
    const f = _csvLine(lines[r]);
    const taxid = iTax >= 0 ? (f[iTax] || "").trim() : "";
    const accession = iAcc >= 0 ? (f[iAcc] || "").trim() : "";
    const name = (iName >= 0 ? (f[iName] || "").trim() : "") || accession || (taxid ? "txid" + taxid : "");
    if (!taxid && !accession) continue;          // need something fetchable from NCBI
    out.push({ name, taxid: taxid || undefined, accession: accession || undefined });
  }
  return out;
}

async function loadPrebuilt(prefix) {
  setStatus("build-status", "Loading prebuilt index…");
  try {
    const man = await (await fetch("./" + prefix + ".web.json")).json();
    const buf = await (await fetch("./" + prefix + ".guides.f64")).arrayBuffer();
    const guides = new Float64Array(buf);
    state.index = new GuideIndex(guides, man.guide_length, man.seed_len, null, man.organisms || []);
    state.indexMeta = man;
    setStatus("build-status", `Loaded ${guides.length.toLocaleString()} guides from ${(man.organisms || []).length} organisms (prebuilt).`);
    updateRunnable();
  } catch (e) { setStatus("build-status", String(e), true); }
}

// --- target ----------------------------------------------------------------
function setTarget(seq, label, info = { type: "unknown" }) {
  state.target = { seq, label, info };
  const box = $("t-status");
  box.classList.add("has");
  box.innerHTML = `✓ Target: <b>${label}</b><br><span class="muted">${seq.length.toLocaleString()} bp · ${info.type}</span>`;
  track("target", targetSummary(info, seq));
  updateRunnable();
}

// --- run -------------------------------------------------------------------
function updateRunnable() { $("run-btn").disabled = !(state.index && state.target); }
const _now = () => (typeof performance !== "undefined" ? performance.now() : Date.now());
function fmtDuration(s) {                      // 7s · 1m 03s · 1h 04m
  s = Math.max(0, Math.round(s));
  if (s < 60) return `${s}s`;
  const m = Math.floor(s / 60), r = s % 60;
  if (m < 60) return `${m}m ${String(r).padStart(2, "0")}s`;
  return `${Math.floor(m / 60)}h ${String(m % 60).padStart(2, "0")}m`;
}
async function run() {
  const btn = $("run-btn");
  if (btn.disabled) return;
  const p = params();
  // optional region: restrict to bases [start, end) of the target
  let seq = state.target.seq, offset = 0;
  const rs = parseInt($("region-start").value, 10), re = parseInt($("region-end").value, 10);
  const start = Number.isFinite(rs) ? Math.max(0, rs) : 0;
  const end = Number.isFinite(re) ? Math.min(seq.length, re) : seq.length;
  if (start > 0 || end < seq.length) {
    if (end <= start) return setStatus("run-status", "Region end must be greater than start.", true);
    seq = seq.slice(start, end); offset = start;
  }
  btn.disabled = true;
  const bar = $("run-progress");
  bar.style.display = "block"; bar.value = 0; bar.max = 1;
  const t0 = _now();
  let rows;
  try {
    rows = await findSparingGuides(seq, state.index, p, {
      maxGuides: p.maxGuides, positionOffset: offset,
      onProgress: (done, total, kept) => {
        bar.max = total || 1; bar.value = done;
        const pct = total ? Math.round((done / total) * 100) : 0;
        // ETA extrapolated from the pace so far (early ticks set it from the first guides)
        let eta = "";
        if (done > 0 && done < total) {
          const remaining = ((_now() - t0) / 1000) * (total - done) / done;
          eta = ` · ~${fmtDuration(remaining)} left`;
        }
        setStatus("run-status",
          `Screening guide ${done.toLocaleString()} / ${total.toLocaleString()} (${pct}%) — ${kept.toLocaleString()} kept${eta}`);
      },
    });
  } catch (e) {
    bar.style.display = "none"; btn.disabled = false;
    return setStatus("run-status", String(e), true);
  }
  bar.style.display = "none"; btn.disabled = false;
  state.results = rows;
  setStatus("run-status",
    `Done in ${fmtDuration((_now() - t0) / 1000)} — ${rows.length.toLocaleString()} commensal-sparing guides.`);
  renderResults(rows);
  track("run", { preset: $("preset").value, pam: p.pam, side: p.side, seedMm: p.seedMm,
    totalMm: p.totalMm, minGc: p.minGc, maxGc: p.maxGc, maxGuides: p.maxGuides,
    region_start: offset || undefined, region_len: seq.length,
    ...targetSummary(state.target.info, state.target.seq), kept: rows.length });
}
function renderResults(rows) {
  $("results-wrap").style.display = "block";
  $("results-count").textContent = `${rows.length} commensal-sparing guides`;
  const tb = $("results"); tb.innerHTML = "";
  rows.slice(0, 2000).forEach((r, i) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${i + 1}</td><td>${r.position}</td><td>${r.strand}</td><td><code>${r.guide_sequence}</code></td><td>${r.gc}</td>`;
    tb.appendChild(tr);
  });
}
function downloadCSV() {
  const rows = state.results || [];
  const csv = ["position,strand,guide_sequence,gc", ...rows.map((r) => `${r.position},${r.strand},${r.guide_sequence},${r.gc}`)].join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  a.download = `${(state.target?.label || "target").replace(/[^A-Za-z0-9]+/g, "_")}_guides.csv`;
  a.click(); URL.revokeObjectURL(a.href);
}

function csvCell(v) { v = String(v ?? ""); return /[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v; }
function downloadPanel() {
  if (!state.panel.length) return setStatus("build-status", "Panel is empty.", true);
  const rows = state.panel.map((s) => `${csvCell(s.name)},${s.taxid || ""},${s.accession || ""}`);
  const csv = ["organism_strain,taxid,accession", ...rows].join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  a.download = "discriminase_panel.csv";
  a.click(); URL.revokeObjectURL(a.href);
}

function setStatus(id, msg, bad) { const e = $(id); e.textContent = msg; e.classList.toggle("bad", !!bad); }

// --- init ------------------------------------------------------------------
// The four fields a preset pins. Editing any of them re-derives the dropdown.
const PRESET_GEOM = ["pam", "side", "guideLength", "seedLen"];
function applyPreset() {
  const p = PRESETS[$("preset").value];
  if (!p) return;                          // "custom" pins nothing
  $("pam").value = p.pam; $("side").value = p.side;
  $("guideLength").value = p.guideLength; $("seedLen").value = p.seedLen; persistParams();
}
// Reflect the live params back onto the dropdown: a known nuclease if they match
// one exactly, otherwise "custom". Keeps the label honest after manual edits.
function syncPresetFromParams() {
  const cur = { pam: $("pam").value.trim().toUpperCase(), side: $("side").value,
                guideLength: +$("guideLength").value, seedLen: +$("seedLen").value };
  const hit = Object.keys(PRESETS).find((k) => PRESET_GEOM.every((f) => PRESETS[k][f] === cur[f]));
  $("preset").value = hit || "custom";
}

function init() {
  restoreParams();
  // restore last panel (metadata only)
  (store.loadSetting("panel", []) || []).forEach((s) => state.panel.push(s));
  renderPanel(); refreshSavedPanels(); loadPrebuiltList(); loadExampleList(); initAnalytics();

  // onboarding guide: a right-side drawer. Opens on first visit, reopen via header.
  const guide = $("guide"), overlay = $("guide-overlay");
  const openGuide = (v) => { guide.classList.toggle("open", v); overlay.classList.toggle("open", v); };
  if (!store.loadSetting("guideSeen", false)) openGuide(true);
  const dismiss = () => { openGuide(false); store.saveSetting("guideSeen", true); };
  $("guide-toggle").onclick = () => openGuide(!guide.classList.contains("open"));
  $("guide-close").onclick = dismiss;
  overlay.onclick = dismiss;

  $("load-panel").onchange = (e) => {
    const v = e.target.value; e.target.value = "";
    if (!v) return;
    const i = v.indexOf(":"), type = v.slice(0, i), id = v.slice(i + 1);
    if (type === "example") loadExamplePanel(id);
    else if (type === "prebuilt") loadPrebuilt(id);
    else if (type === "saved") loadSavedPanel(id);
  };
  $("preset").onchange = applyPreset;
  PARAM_IDS.forEach((id) => $(id).addEventListener("change", persistParams));
  // flip the preset label to "custom" (or a matching nuclease) on manual edits
  PRESET_GEOM.forEach((id) => $(id).addEventListener("input", syncPresetFromParams));
  syncPresetFromParams();                  // and reflect restored params on load

  $("cp-search").onclick = () => $("cp-name").value.trim() && searchInto($("cp-name").value.trim(), $("cp-results"),
    (c) => addSource({ name: c.organism || c.title.slice(0, 40) || c.accession, accession: c.accession, taxid: c.taxid }));
  $("cp-add-acc").onclick = () => { const a = $("cp-acc").value.trim(); if (a) { addSource({ name: a, accession: a }); $("cp-acc").value = ""; } };
  $("cp-add-tax").onclick = () => { const t = $("cp-tax").value.trim(); if (t) { addSource({ name: "txid" + t, taxid: t }); $("cp-tax").value = ""; } };
  $("cp-file").onchange = async (e) => { for (const f of e.target.files) addSource({ name: f.name, seq: await readFasta(f), uploaded: true }); e.target.value = ""; };
  $("build-btn").onclick = buildOrLoad;
  $("cp-save").onclick = saveNamedPanel;
  $("panel-csv").onchange = async (e) => {
    const f = e.target.files[0]; if (!f) return;
    loadPanelSources(parsePanelCsv(await f.text()), f.name); e.target.value = "";
  };

  $("t-search").onclick = () => { const q = $("t-name").value.trim(); if (!q) return; track("search_target", { query: q }); searchInto(q, $("t-results"),
    async (c) => { setStatus("t-status", `Fetching ${c.accession}…`); try { const { seq } = await ncbi.fetchSequence(c.accession, { email: params().email }); setTarget(seq, c.accession, { type: "name", value: c.accession, name: q }); } catch (err) { setStatus("t-status", String(err), true); } }, true); };
  $("t-use-acc").onclick = async () => { const a = $("t-acc").value.trim(); if (!a) return; setStatus("t-status", `Fetching ${a}…`); try { const { seq } = await ncbi.fetchSequence(a, { email: params().email }); setTarget(seq, a, { type: "accession", value: a }); } catch (e) { setStatus("t-status", String(e), true); } };
  $("t-file").onchange = async (e) => { const f = e.target.files[0]; if (f) setTarget(await readFasta(f), f.name, { type: "upload", name: f.name }); };
  $("t-use-paste").onclick = () => { const s = $("t-paste").value.replace(/\s/g, ""); if (s) setTarget(s, "pasted_sequence", { type: "paste" }); };

  $("run-btn").onclick = run;
  $("dl-btn").onclick = downloadCSV;
  $("cp-download").onclick = downloadPanel;
  updateRunnable();
}

init();
