// UI controller. Holds no algorithm logic itself — it wires the DOM to the
// verified engine (engine.mjs / build.mjs / ncbi.mjs) and persistence (store.mjs).

import { GuideIndex } from "./engine.mjs";
import { buildIndex, collectCommonTargetGuides, screenTargetGuides } from "./build.mjs";
import * as ncbi from "./ncbi.mjs";
import * as store from "./store.mjs";
import { initAnalytics, track } from "./analytics.mjs";

const $ = (id) => document.getElementById(id);
const PRESETS = {
  cas12a: { pam: "TTTV", side: "5prime", guideLength: 23, seedLen: 10 },
  spcas9: { pam: "NGG", side: "3prime", guideLength: 20, seedLen: 10 },
};
const PARAM_IDS = ["pam", "side", "guideLength", "seedLen", "seedMm",
                   "totalMm", "minGc", "maxGc", "maxGuides", "email"];

const state = { panel: [], index: null, indexMeta: null,
                targets: [], targetGuides: null, results: null, cancel: false, pasteCount: 0 };

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

// --- source tables (shared by the commensal panel and the target panel) ----
function srcLabel(s) { return s.accession ? `acc ${s.accession}` : s.taxid ? `txid ${s.taxid}` : s.uploaded ? "uploaded" : "—"; }
function renderSourceTable(tbodyId, list, emptyMsg, onChange) {
  const tb = $(tbodyId);
  tb.innerHTML = "";
  if (!list.length) { tb.innerHTML = `<tr><td colspan="4" class="muted">${emptyMsg}</td></tr>`; return; }
  list.forEach((s, i) => {
    const tr = document.createElement("tr");
    const status = s.status ? `<span class="pill ${s.status === "failed" ? "bad" : "ok"}">${s.status}${s.guides != null ? " · " + s.guides : ""}</span>` : "";
    tr.innerHTML = `<td>${s.name}</td><td class="muted">${srcLabel(s)}</td><td>${status}</td><td><button class="ghost">✕</button></td>`;
    tr.querySelector("button").onclick = () => { list.splice(i, 1); onChange(); };
    tb.appendChild(tr);
  });
}
function sourceKey(s) { return s.accession || s.taxid || ("up:" + s.name); }

// --- commensal panel -------------------------------------------------------
function renderPanel() {
  renderSourceTable("panel-list", state.panel, "No commensals yet. Add some on the right.",
    () => { renderPanel(); savePanelState(); });
}
function addSource(s) {
  if (state.panel.some((x) => sourceKey(x) === sourceKey(s))) return;
  state.panel.push(s); renderPanel(); savePanelState();
}
function savePanelState() {
  // persist metadata only (not uploaded sequences — too big for localStorage)
  store.saveSetting("panel", state.panel.filter((s) => !s.uploaded).map((s) => ({ name: s.name, accession: s.accession, taxid: s.taxid })));
}

// --- target panel (mirrors the commensal panel) ----------------------------
function renderTargets() {
  renderSourceTable("target-list", state.targets, "No targets yet. Add some on the right.",
    () => { state.targetGuides = null; renderTargets(); saveTargetState(); updateRunnable(); });
}
function addTarget(s) {
  if (state.targets.some((x) => sourceKey(x) === sourceKey(s))) return;
  state.targets.push(s); state.targetGuides = null;   // panel changed -> recompute needed
  renderTargets(); saveTargetState(); updateRunnable();
}
function saveTargetState() {
  store.saveSetting("targets", state.targets.filter((s) => !s.uploaded).map((s) => ({ name: s.name, accession: s.accession, taxid: s.taxid })));
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

function fetchSeqFor(s, email) {
  return s.accession ? ncbi.fetchSequence(s.accession, { email }) : ncbi.fetchSequenceByTaxid(s.taxid, { email });
}

// --- build / load commensal index ------------------------------------------
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
      fetchSeq: (s) => fetchSeqFor(s, p.email),
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
  setStatus("build-status", `Loaded ${sources.length} organisms${label ? " from " + label : ""}. Now press “Build commensal index”.`);
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

// --- find common target guides (Step 2 action) -----------------------------
// The region inputs, parsed; only honoured when there is a single target.
function regionFromInputs() {
  const rs = parseInt($("region-start").value, 10), re = parseInt($("region-end").value, 10);
  if (!Number.isFinite(rs) && !Number.isFinite(re)) return null;
  return { start: Number.isFinite(rs) ? Math.max(0, rs) : 0,
           end: Number.isFinite(re) ? re : Infinity };
}
async function findCommonTargetGuides() {
  if (!state.targets.length) return setStatus("target-status", "Add at least one target first.", true);
  const p = params();
  state.targets.forEach((s) => { s.status = undefined; s.guides = undefined; });   // clear stale row pills
  renderTargets();
  let region = null;
  if (state.targets.length === 1) {
    region = regionFromInputs();
    if (region && region.end <= region.start) return setStatus("target-status", "Region end must be greater than start.", true);
  }
  $("target-btn").disabled = true;
  try {
    const res = await collectCommonTargetGuides(state.targets, p, {
      fetchSeq: (s) => fetchSeqFor(s, p.email),
      onProgress: (ev) => {
        const src = state.targets.find((s) => s.name === ev.name);
        if (src) { src.status = ev.status; src.guides = ev.guides; }
        renderTargets();
        setStatus("target-status", `${ev.status} ${ev.name} (${ev.i + 1}/${ev.n})…`);
      },
      region,
    });
    state.targetGuides = res;
    const failNote = res.failed.length ? ` · ${res.failed.length} failed` : "";
    const scope = res.organisms.length > 1 ? `shared by all ${res.organisms.length} targets` : `in ${res.organisms[0]?.name ?? "target"}`;
    const capNote = p.maxGuides && res.total > p.maxGuides ? ` — Run screens the first ${p.maxGuides.toLocaleString()} (Max guides)` : "";
    setStatus("target-status", `${res.total.toLocaleString()} common target guides ${scope}${failNote}${capNote}.`, res.total === 0);
    track("find_targets", { preset: $("preset").value, pam: p.pam, side: p.side,
      n_targets: res.organisms.length, common: res.total, failed: res.failed.length });
  } catch (e) {
    setStatus("target-status", String(e), true);
  } finally { $("target-btn").disabled = false; updateRunnable(); }
}

// --- run -------------------------------------------------------------------
function updateRunnable() {
  $("run-btn").disabled = !(state.index && state.targetGuides && state.targetGuides.total > 0);
}
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
  const tg = state.targetGuides;
  btn.disabled = true;
  const cancelBtn = $("cancel-btn");
  state.cancel = false;
  cancelBtn.style.display = ""; cancelBtn.disabled = false;
  const bar = $("run-progress");
  bar.style.display = "block"; bar.value = 0; bar.max = 1;
  const t0 = _now();
  let result;
  try {
    result = await screenTargetGuides(tg, state.index, p, {
      maxGuides: p.maxGuides,
      shouldStop: () => state.cancel,
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
          `Screening target guide ${done.toLocaleString()} / ${total.toLocaleString()} (${pct}%) — ${kept.toLocaleString()} kept${eta}`);
      },
    });
  } catch (e) {
    bar.style.display = "none"; btn.disabled = false; cancelBtn.style.display = "none";
    return setStatus("run-status", String(e), true);
  }
  const { rows, scanned, total, cancelled } = result;
  bar.style.display = "none"; btn.disabled = false; cancelBtn.style.display = "none";
  state.results = rows;
  const capped = !cancelled && scanned < total;
  const dur = fmtDuration((_now() - t0) / 1000);
  setStatus("run-status", cancelled
    ? `Stopped after ${dur} — showing ${rows.length.toLocaleString()} guides from the first ${scanned.toLocaleString()} / ${total.toLocaleString()} target guides screened.`
    : `Done in ${dur} — screened ${scanned.toLocaleString()} / ${total.toLocaleString()} target guides.`);
  renderResults(rows, { scanned, total, capped, cancelled, maxGuides: p.maxGuides });
  track("run", { preset: $("preset").value, pam: p.pam, side: p.side, seedMm: p.seedMm,
    totalMm: p.totalMm, minGc: p.minGc, maxGc: p.maxGc, maxGuides: p.maxGuides,
    n_targets: tg.organisms.length, target_guides: total, scanned, capped, cancelled, kept: rows.length });
}
function renderResults(rows, meta = {}) {
  $("results-wrap").style.display = "block";
  const { scanned, total, capped, cancelled, maxGuides } = meta;
  let summary = `${rows.length.toLocaleString()} commensal-sparing guides`;
  if (total != null) summary += ` · screened ${scanned.toLocaleString()} / ${total.toLocaleString()} target guides`;
  $("results-count").textContent = summary;
  // explain when the scan didn't cover the whole target guide universe
  const note = $("results-note");
  if (cancelled) {
    note.textContent = `Stopped early — the ${rows.length.toLocaleString()} guides above are what was found in the first ${scanned.toLocaleString()} of ${total.toLocaleString()} target guides. Press Run to screen the whole target from the start.`;
    note.style.display = "block";
  } else if (capped) {
    note.textContent = `Stopped at the Max guides cap (${maxGuides.toLocaleString()}) — ${(total - scanned).toLocaleString()} target guides weren't screened. Raise “Max guides”, or narrow to a region (Step 2) to search a different part of the genome.`;
    note.style.display = "block";
  } else {
    note.style.display = "none";
  }
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
  const orgs = state.targetGuides?.organisms || [];
  const base = orgs.length > 1 ? "common_targets" : (orgs[0]?.name || "target");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  a.download = `${base.replace(/[^A-Za-z0-9]+/g, "_")}_guides.csv`;
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
  invalidateGuides();
}
// Reflect the live params back onto the dropdown: a known nuclease if they match
// one exactly, otherwise "custom". Keeps the label honest after manual edits.
function syncPresetFromParams() {
  const cur = { pam: $("pam").value.trim().toUpperCase(), side: $("side").value,
                guideLength: +$("guideLength").value, seedLen: +$("seedLen").value };
  const hit = Object.keys(PRESETS).find((k) => PRESET_GEOM.every((f) => PRESETS[k][f] === cur[f]));
  $("preset").value = hit || "custom";
}
// Changing geometry or region makes any computed common-target-guide set stale.
function invalidateGuides() {
  if (state.targetGuides) { state.targetGuides = null; setStatus("target-status", "Geometry changed — press “Find common target guides” again."); }
  updateRunnable();
}

function init() {
  restoreParams();
  // restore last panel + targets (metadata only)
  (store.loadSetting("panel", []) || []).forEach((s) => state.panel.push(s));
  (store.loadSetting("targets", []) || []).forEach((s) => state.targets.push(s));
  renderPanel(); renderTargets();
  refreshSavedPanels(); loadPrebuiltList(); loadExampleList(); initAnalytics();

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
  // geometry / region edits invalidate the common-target-guide set
  [...PRESET_GEOM, "region-start", "region-end"].forEach((id) =>
    $(id).addEventListener("change", invalidateGuides));
  syncPresetFromParams();                  // and reflect restored params on load

  // commensal panel
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

  // target panel (mirrors the commensal panel; adds the paste-DNA option)
  $("t-search").onclick = () => { const q = $("t-name").value.trim(); if (!q) return; track("search_target", { query: q }); searchInto(q, $("t-results"),
    (c) => addTarget({ name: c.organism || c.title.slice(0, 40) || c.accession, accession: c.accession, taxid: c.taxid })); };
  $("t-add-acc").onclick = () => { const a = $("t-acc").value.trim(); if (a) { addTarget({ name: a, accession: a }); $("t-acc").value = ""; } };
  $("t-add-tax").onclick = () => { const t = $("t-tax").value.trim(); if (t) { addTarget({ name: "txid" + t, taxid: t }); $("t-tax").value = ""; } };
  $("t-file").onchange = async (e) => { for (const f of e.target.files) addTarget({ name: f.name, seq: await readFasta(f), uploaded: true }); e.target.value = ""; };
  $("t-add-paste").onclick = () => { const s = $("t-paste").value.replace(/\s/g, ""); if (s) { addTarget({ name: `pasted #${++state.pasteCount} (${s.length.toLocaleString()} bp)`, seq: s, uploaded: true }); $("t-paste").value = ""; } };
  $("target-btn").onclick = findCommonTargetGuides;

  // run
  $("run-btn").onclick = run;
  $("cancel-btn").onclick = () => { state.cancel = true; $("cancel-btn").disabled = true; setStatus("run-status", "Stopping — showing results so far…"); };
  $("dl-btn").onclick = downloadCSV;
  $("cp-download").onclick = downloadPanel;
  updateRunnable();
}

init();
