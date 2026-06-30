// UI controller. Holds no algorithm logic itself — it wires the DOM to the
// verified engine (engine.mjs / build.mjs / ncbi.mjs) and persistence (store.mjs).

import { GuideIndex } from "./engine.mjs";
import { buildIndex, findSparingGuides } from "./build.mjs";
import * as ncbi from "./ncbi.mjs";
import * as store from "./store.mjs";
import { initAnalytics, track, targetSummary } from "./analytics.mjs";

const $ = (id) => document.getElementById(id);
const PRESETS = {
  cas12a: { pam: "TTTV", side: "5prime", gap: 1, guideLength: 23, seedLen: 10 },
  spcas9: { pam: "NGG", side: "3prime", gap: 0, guideLength: 20, seedLen: 10 },
};
const PARAM_IDS = ["pam", "side", "gap", "guideLength", "seedLen", "seedMm",
                   "totalMm", "minGc", "maxGc", "maxGuides", "email"];

const state = { panel: [], index: null, indexMeta: null, target: null, results: null };

// --- parameters ------------------------------------------------------------
function params() {
  return {
    pam: $("pam").value.trim().toUpperCase(), side: $("side").value,
    gap: +$("gap").value, guideLength: +$("guideLength").value, seedLen: +$("seedLen").value,
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

async function searchInto(query, container, onPick) {
  container.style.display = "block";
  container.innerHTML = `<p class="muted" style="padding:8px">Searching…</p>`;
  try {
    const cands = await ncbi.searchCandidates(query, { email: params().email, retmax: 20 });
    if (!cands.length) { container.innerHTML = `<p class="muted" style="padding:8px">No matches.</p>`; return; }
    const t = document.createElement("table");
    t.innerHTML = `<thead><tr><th>accession</th><th>length</th><th>description</th><th></th></tr></thead><tbody></tbody>`;
    cands.forEach((c) => {
      const tr = document.createElement("tr");
      tr.innerHTML = `<td>${c.accession}</td><td>${c.length_bp ? c.length_bp.toLocaleString() + " bp" : "?"}</td><td>${c.title.slice(0, 64)}</td><td><button class="ghost">pick</button></td>`;
      tr.querySelector("button").onclick = () => onPick(c);
      t.querySelector("tbody").appendChild(tr);
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
  const geom = { L: p.guideLength, pam: p.pam, side: p.side, gap: p.gap, seed: p.seedLen };
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
      organisms: index.organisms.map((o) => o.name), n: index.organisms.length, guides: index.length, failed: failed.length });
  } catch (e) {
    setStatus("build-status", String(e), true);
  } finally { $("build-btn").disabled = false; updateRunnable(); }
}

// --- saved panels ----------------------------------------------------------
function refreshSavedPanels() {
  const names = store.loadSetting("panelNames", []);
  const sel = $("saved-panels");
  sel.innerHTML = `<option value="">— saved panels —</option>` + names.map((n) => `<option>${n}</option>`).join("");
}
function saveNamedPanel() {
  const name = prompt("Save this panel as:");
  if (!name) return;
  const names = store.loadSetting("panelNames", []);
  if (!names.includes(name)) names.push(name);
  store.saveSetting("panelNames", names);
  store.saveSetting("panel:" + name, state.panel.filter((s) => !s.uploaded).map((s) => ({ name: s.name, accession: s.accession, taxid: s.taxid })));
  refreshSavedPanels();
}

// --- prebuilt panels (static files shipped with the site) ------------------
async function loadPrebuiltList() {
  try {
    const list = await (await fetch("./panels/index.json")).json();
    const sel = $("prebuilt");
    list.forEach((p) => { const o = document.createElement("option"); o.value = p.prefix; o.textContent = p.name; sel.appendChild(o); });
  } catch (e) { /* no prebuilt panels shipped — fine */ }
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
  $("t-status").innerHTML = `Target: <b>${label}</b> (${seq.length.toLocaleString()} bp)`;
  track("target", targetSummary(info, seq));
  updateRunnable();
}

// --- run -------------------------------------------------------------------
function updateRunnable() { $("run-btn").disabled = !(state.index && state.target); }
function run() {
  const p = params();
  $("run-progress").style.display = "block";
  const rows = findSparingGuides(state.target.seq, state.index, p, {
    maxGuides: p.maxGuides,
    onProgress: (done, total, kept) => {
      $("run-progress").max = total; $("run-progress").value = done;
      setStatus("run-status", `screening ${done.toLocaleString()}/${total.toLocaleString()} — ${kept} kept`);
    },
  });
  state.results = rows;
  $("run-progress").style.display = "none";
  renderResults(rows);
  track("run", { preset: $("preset").value, pam: p.pam, side: p.side, seedMm: p.seedMm,
    totalMm: p.totalMm, minGc: p.minGc, maxGc: p.maxGc, maxGuides: p.maxGuides,
    target: targetSummary(state.target.info, state.target.seq), kept: rows.length });
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

function setStatus(id, msg, bad) { const e = $(id); e.textContent = msg; e.className = "muted " + (bad ? "bad" : ""); }

// --- init ------------------------------------------------------------------
function applyPreset() {
  const p = PRESETS[$("preset").value];
  if (!p) return;
  $("pam").value = p.pam; $("side").value = p.side; $("gap").value = p.gap;
  $("guideLength").value = p.guideLength; $("seedLen").value = p.seedLen; persistParams();
}

function init() {
  restoreParams();
  // restore last panel (metadata only)
  (store.loadSetting("panel", []) || []).forEach((s) => state.panel.push(s));
  renderPanel(); refreshSavedPanels(); loadPrebuiltList(); initAnalytics();

  // onboarding guide: show until dismissed, reopen via the header button
  const guide = $("guide");
  if (store.loadSetting("guideSeen", false)) guide.style.display = "none";
  $("guide-toggle").onclick = () => { guide.style.display = guide.style.display === "none" ? "block" : "none"; };
  $("guide-close").onclick = () => { guide.style.display = "none"; store.saveSetting("guideSeen", true); };

  $("prebuilt").onchange = (e) => { if (e.target.value) loadPrebuilt(e.target.value); };
  $("preset").onchange = applyPreset;
  PARAM_IDS.forEach((id) => $(id).addEventListener("change", persistParams));

  $("cp-search").onclick = () => $("cp-name").value.trim() && searchInto($("cp-name").value.trim(), $("cp-results"),
    (c) => addSource({ name: c.organism || c.title.slice(0, 40) || c.accession, accession: c.accession, taxid: c.taxid }));
  $("cp-add-acc").onclick = () => { const a = $("cp-acc").value.trim(); if (a) { addSource({ name: a, accession: a }); $("cp-acc").value = ""; } };
  $("cp-add-tax").onclick = () => { const t = $("cp-tax").value.trim(); if (t) { addSource({ name: "txid" + t, taxid: t }); $("cp-tax").value = ""; } };
  $("cp-file").onchange = async (e) => { for (const f of e.target.files) addSource({ name: f.name, seq: await readFasta(f), uploaded: true }); e.target.value = ""; };
  $("build-btn").onclick = buildOrLoad;
  $("cp-save").onclick = saveNamedPanel;
  $("saved-panels").onchange = (e) => {
    if (!e.target.value) return;
    state.panel = (store.loadSetting("panel:" + e.target.value, []) || []).map((s) => ({ ...s }));
    renderPanel(); savePanelState();
  };

  $("t-search").onclick = () => { const q = $("t-name").value.trim(); if (!q) return; track("search_target", { query: q }); searchInto(q, $("t-results"),
    async (c) => { setStatus("t-status", `Fetching ${c.accession}…`); try { const { seq } = await ncbi.fetchSequence(c.accession, { email: params().email }); setTarget(seq, c.accession, { type: "name", value: c.accession, name: q }); } catch (err) { setStatus("t-status", String(err), true); } }); };
  $("t-use-acc").onclick = async () => { const a = $("t-acc").value.trim(); if (!a) return; setStatus("t-status", `Fetching ${a}…`); try { const { seq } = await ncbi.fetchSequence(a, { email: params().email }); setTarget(seq, a, { type: "accession", value: a }); } catch (e) { setStatus("t-status", String(e), true); } };
  $("t-file").onchange = async (e) => { const f = e.target.files[0]; if (f) setTarget(await readFasta(f), f.name, { type: "upload", name: f.name }); };
  $("t-use-paste").onclick = () => { const s = $("t-paste").value.replace(/\s/g, ""); if (s) setTarget(s, "pasted_sequence", { type: "paste" }); };

  $("run-btn").onclick = run;
  $("dl-btn").onclick = downloadCSV;
  updateRunnable();
}

init();
