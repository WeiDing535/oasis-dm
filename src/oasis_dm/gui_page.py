"""The single HTML page served by oasis_dm.gui (no external assets).

``__TOKEN__`` is replaced with the per-session token at serve time; every
API call sends it in the X-Oasis-Token header so other local processes or
cross-origin pages cannot drive the engine.
"""
PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>oasis-dm</title>
<style>
  :root {
    --ink: #1c2733; --mut: #64748b; --line: #d8dee6; --card: #ffffff;
    --acc: #0369a1; --ok: #059669; --bad: #b91c1c; --bg: #eef1f5;
  }
  * { box-sizing: border-box; }
  body { margin: 0; background: var(--bg); color: var(--ink);
         font: 14px/1.45 -apple-system, "Segoe UI", Roboto, sans-serif; }
  .wrap { max-width: 960px; margin: 0 auto; padding: 18px 16px 60px; }
  h1 { font-size: 19px; margin: 6px 0 2px; }
  h1 small { color: var(--mut); font-weight: 400; margin-left: 8px; }
  .sub { color: var(--mut); margin: 0 0 14px; }
  .card { background: var(--card); border: 1px solid var(--line);
          border-radius: 8px; padding: 12px 14px; margin-bottom: 12px; }
  .card > h2 { font-size: 13px; text-transform: uppercase; letter-spacing: .04em;
               color: var(--mut); margin: 0 0 10px; }
  .row { display: flex; flex-wrap: wrap; gap: 8px 14px; align-items: center; }
  .row + .row { margin-top: 8px; }
  .f { display: flex; align-items: center; gap: 6px; flex: 1 1 220px; }
  .f label { color: var(--mut); white-space: nowrap; }
  input[type=text], input[type=number], select {
    font: inherit; padding: 4px 7px; border: 1px solid var(--line);
    border-radius: 5px; background: #fff; min-width: 0; flex: 1; }
  input[type=number] { width: 76px; flex: 0; }
  .cols input[type=text], .cols select { flex: 0 1 130px; }
  button { font: inherit; padding: 5px 13px; border: 1px solid var(--line);
           border-radius: 5px; background: #f6f8fa; cursor: pointer; }
  button:hover { background: #eef2f6; }
  button.primary { background: var(--acc); border-color: var(--acc);
                   color: #fff; font-weight: 600; padding: 7px 26px; }
  button.primary:hover { background: #075e8e; }
  button.danger { background:#fee2e2; border-color:#fca5a5; color: var(--bad); }
  button:disabled { opacity: .45; cursor: default; }
  .info { color: var(--mut); font-size: 12.5px; }
  .cmd { font: 12px/1.5 ui-monospace, Menlo, Consolas, monospace; color: #334;
         background: #f4f6f9; border: 1px solid var(--line); border-radius: 6px;
         padding: 8px 10px; white-space: pre-wrap; word-break: break-all; }
  #status { margin-left: auto; font-weight: 600; }
  #status.ok { color: var(--ok); } #status.bad { color: var(--bad); }
  #log { font: 12px/1.55 ui-monospace, Menlo, Consolas, monospace;
         background: #10151c; color: #cfe3f2; border-radius: 6px;
         padding: 10px 12px; height: 300px; overflow-y: auto;
         white-space: pre-wrap; word-break: break-all; }
  #chart { width: 100%; height: 120px; background: #fff; border: 1px solid
           var(--line); border-radius: 6px; margin-bottom: 10px; }
  .dlg { position: fixed; inset: 0; background: rgba(15,23,32,.45);
         display: none; align-items: center; justify-content: center; }
  .dlg .box { background: #fff; border-radius: 10px; width: min(560px, 92vw);
              padding: 14px; box-shadow: 0 18px 50px rgba(0,0,0,.25); }
  .dlg .list { border: 1px solid var(--line); border-radius: 6px;
               height: 320px; overflow-y: auto; margin: 10px 0; }
  .dlg .list div { padding: 5px 10px; cursor: pointer; display: flex;
                   gap: 8px; }
  .dlg .list div:hover { background: #e8f2fa; }
  .dlg .list .dir { color: var(--acc); }
  .dlg .list .mtz { font-weight: 600; }
  .dlg .list .dim { color: var(--mut); margin-left: auto; }
  label.radio { display: flex; gap: 5px; align-items: center; }
  .problems { color: var(--bad); margin: 6px 0 0; }
</style>
</head>
<body>
<div class="wrap">
  <h1>oasis-dm <small>density modification</small></h1>
  <p class="sub">Solvent flattening, histogram matching and phase extension
     with omit-based figure-of-merit estimation. The command line equivalent
     of every run is shown before you launch it.</p>

  <div class="card">
    <h2>Files</h2>
    <div class="row">
      <div class="f"><label>input MTZ</label>
        <input type="text" id="mtzin" spellcheck="false">
        <button id="b_in">Browse…</button></div>
    </div>
    <div class="row">
      <div class="f"><label>output MTZ</label>
        <input type="text" id="mtzout" spellcheck="false">
        <button id="b_out">Browse…</button></div>
    </div>
    <div class="row"><span class="info" id="mtzinfo">no input selected yet</span></div>
  </div>

  <div class="card cols">
    <h2>Columns (auto-detected from the input)</h2>
    <div class="row">
      <div class="f"><label>Fo</label><select id="fo"></select></div>
      <div class="f"><label>σ(Fo)</label><select id="sigf"></select></div>
      <div class="f"><label>PHIB</label><select id="phib"></select></div>
      <div class="f"><label>FOM</label><select id="fom"></select></div>
      <div class="f"><label>free set</label><select id="free"></select></div>
      <div class="f"><label>truth</label><select id="truth"></select></div>
    </div>
  </div>

  <div class="card">
    <h2>Parameters</h2>
    <div class="row">
      <div class="f"><label>solvent content</label>
        <input type="number" id="solc" value="0.5" step="0.01" min="0.05" max="0.95"></div>
      <div class="f"><label>max cycles</label>
        <input type="number" id="cycles" value="8" step="1" min="1" max="100"></div>
      <div class="f"><label>free-set early stop</label>
        <input type="checkbox" id="auto_stop" checked></div>
    </div>
    <div class="row">
      <div class="f"><label>mask</label>
        <select id="mask"><option>variance</option><option>wang</option></select></div>
      <div class="f"><label>histogram target</label>
        <select id="hist_mode"><option>cell</option><option>protein</option></select></div>
      <div class="f"><label>free metric</label>
        <select id="free_metric"><option>cc</option><option>rfree</option></select></div>
    </div>
  </div>

  <div class="card">
    <h2>Options</h2>
    <div class="row">
      <div class="f"><label>phase extension</label>
        <input type="checkbox" id="extend" checked></div>
      <div class="f"><label>anisotropic correction</label>
        <input type="checkbox" id="aniso"></div>
      <div class="f"><label>per-shell D in output</label>
        <input type="checkbox" id="per_shell_d"></div>
    </div>
  </div>

  <div class="card">
    <h2>Non-crystallographic symmetry — expert option: averaging only pays
       with tight operators and several copies</h2>
    <div class="row">
      <label class="radio"><input type="radio" name="ncs" value="none" checked>none</label>
      <label class="radio"><input type="radio" name="ncs" value="file">operator file</label>
      <label class="radio"><input type="radio" name="ncs" value="auto">auto-detect (free-set gated)</label>
    </div>
    <div class="row">
      <div class="f"><label>operator file</label>
        <input type="text" id="ncs_file" spellcheck="false" disabled>
        <button id="b_ncs">Browse…</button></div>
    </div>
  </div>

  <div class="card">
    <h2>Advanced</h2>
    <div class="row">
      <div class="f"><label>SigmaQ source</label>
        <select id="sigmaq_source">
          <option>blend</option><option>kfold</option><option>kcal</option>
          <option>omit</option><option>work</option></select></div>
    </div>
    <div class="row">
      <div class="f"><label>hist.lib</label>
        <input type="text" id="hist_lib" spellcheck="false">
        <button id="b_hist">Browse…</button></div>
      <div class="f"><label>no histogram library</label>
        <input type="checkbox" id="no_hist"></div>
    </div>
  </div>

  <div class="card">
    <h2>Command</h2>
    <div class="cmd" id="cmd">oasis-dm --mtzin … --mtzout …</div>
    <div class="problems" id="problems"></div>
  </div>

  <div class="row" style="margin-bottom: 12px;">
    <button class="primary" id="run">Run</button>
    <button class="danger" id="stop" disabled>Stop</button>
    <span id="status"></span>
  </div>

  <svg id="chart" preserveAspectRatio="none"></svg>
  <div id="log">The engine log will appear here once you press Run.</div>
</div>

<div class="dlg" id="dlg">
  <div class="box">
    <div class="row">
      <b id="dlg_title">Choose a file</b>
      <button id="dlg_up" style="margin-left:auto">⬑ up</button>
      <button id="dlg_home">home</button>
    </div>
    <div class="list" id="dlg_list"></div>
    <div class="row">
      <input type="text" id="dlg_path" spellcheck="false" style="flex:1">
      <button id="dlg_pick">Use this path</button>
      <button id="dlg_cancel">Cancel</button>
    </div>
  </div>
</div>

<script>
"use strict";
const TOKEN = "__TOKEN__";
const $ = id => document.getElementById(id);
const api = {};

async function call(method, route, body) {
  const opt = {method, headers: {"X-Oasis-Token": TOKEN}};
  if (body !== undefined) {
    opt.headers["Content-Type"] = "application/json";
    opt.body = JSON.stringify(body);
  }
  const r = await fetch(route, opt);
  return r.json();
}
api.open   = p    => call("POST", "/api/open", {path: p});
api.list   = p    => call("POST", "/api/list", {path: p});
api.preview= st   => call("POST", "/api/preview", {state: st});
api.run    = st   => call("POST", "/api/run", {state: st});
api.stop   = ()   => call("POST", "/api/stop", {});
api.poll   = seq  => call("GET", "/api/poll?seq=" + seq);

// ---------- state ------------------------------------------------------------
const IDS = ["mtzin","mtzout","fo","sigf","phib","fom","free","truth","solc",
             "cycles","auto_stop","mask","hist_mode","free_metric","extend",
             "aniso","per_shell_d","ncs_file","no_hist","hist_lib",
             "sigmaq_source"];
function state() {
  const st = {};
  for (const id of IDS) st[id] = $(id).value.trim();
  st.solc = parseFloat(st.solc); if (isNaN(st.solc)) st.solc = 0.5;
  st.cycles = parseInt(st.cycles, 10); if (isNaN(st.cycles)) st.cycles = 8;
  for (const k of ["auto_stop","extend","aniso","per_shell_d","no_hist"])
    st[k] = $(k).checked;
  st.ncs = (document.querySelector("input[name=ncs]:checked")||{}).value
           || "none";
  if (st.ncs !== "file") st.ncs_file = "";
  if (st.no_hist) st.hist_lib = "";
  st.tag = st.mtzin ? (st.mtzin.split(/[\\/]/).pop()||"").replace(/\.mtz$/i,"")
                    : "dm";
  return st;
}

// ---------- command preview --------------------------------------------------
let previewTimer = null;
function schedulePreview() {
  clearTimeout(previewTimer);
  previewTimer = setTimeout(async () => {
    const r = await api.preview(state());
    if (r.ok) $("cmd").textContent = "$ " + r.cmd;
  }, 250);
}
for (const id of IDS) $(id).addEventListener("input", schedulePreview);
for (const el of document.querySelectorAll("input[type=checkbox],input[name=ncs],select"))
  el.addEventListener("change", schedulePreview);
document.querySelectorAll("input[name=ncs]").forEach(el =>
  el.addEventListener("change", () => { $("ncs_file").disabled =
    state().ncs !== "file"; }));
$("no_hist").addEventListener("change", () => {
  if ($("no_hist").checked) $("hist_lib").value = ""; });

// ---------- input file -------------------------------------------------------
$("mtzin").addEventListener("change", async () => {
  const p = $("mtzin").value.trim();
  if (!p) return;
  const r = await api.open(p);
  if (!r.ok) { $("mtzinfo").textContent = "cannot read: " + r.error; return; }
  const wants = {fo:["FOBS","FP","F"], sigf:["SIGFOBS","SIGFP","SIGF"],
                 phib:["PHIB","PHI"], fom:["FOM","WT"],
                 free:["FreeR","FREE"], truth:["PHIC"]};
  for (const [key, cand] of Object.entries(wants)) {
    const sel = $(key);
    sel.innerHTML = "";
    for (const lab of r.labels.concat(["(none)"])) {
      const o = document.createElement("option"); o.textContent = lab;
      sel.appendChild(o);
    }
    const hit = cand.find(c => r.labels.includes(c));
    sel.value = hit || (r.labels[0] || "(none)");
  }
  $("mtzinfo").textContent =
    (p.split(/[\\/]/).pop()) + " · " + r.spacegroup + " · " +
    r.nref + " reflections · d_min " + r.dmin.toFixed(2) + " Å";
  if (!$("mtzout").value.trim() && p) {
    const m = p.match(/^(.*?)(\.[^.\\/]*)?$/);
    $("mtzout").value = m[1] + "_dm" + (m[2] || ".mtz");
  }
  schedulePreview();
});

// ---------- file browser dialog ---------------------------------------------
const dlg = {target: null, filter: null, cwd: ""};
function openDlg(targetId, title, filter) {
  dlg.target = targetId; dlg.filter = filter;
  $("dlg_title").textContent = title;
  $("dlg").style.display = "flex";
  $("dlg_path").value = $(targetId).value.trim() ||
    (dlg.cwd || (location.hostname ? "." : "."));
  loadDir(dirname($("dlg_path").value) || ".");
}
function dirname(p) {
  const i = Math.max(p.lastIndexOf("/"), p.lastIndexOf("\\"));
  return i > 0 ? p.slice(0, i) : "";
}
async function loadDir(dir) {
  const r = await api.list(dir || ".");
  if (!r.ok) { $("dlg_path").value = dir; return; }
  dlg.cwd = r.cwd;
  $("dlg_list").innerHTML = "";
  for (const e of r.entries) {
    if (!e.dir && dlg.filter && !dlg.filter(e)) continue;
    const div = document.createElement("div");
    div.innerHTML = (e.dir ? "📁 " : "") +
      '<span class="' + (e.dir ? "dir" : (e.mtz ? "mtz" : "")) + '">' +
      e.name + '</span>' + (e.dir ? '<span class="dim">dir</span>' : "");
    div.onclick = () => {
      const path = dlg.cwd.replace(/[\\/]$/, "") + "/" + e.name;
      if (e.dir) { loadDir(path); }
      else { $("dlg_path").value = path; pick(); }
    };
    $("dlg_list").appendChild(div);
  }
}
function pick() {
  const v = $("dlg_path").value.trim();
  if (!v) return;
  $(dlg.target).value = v;
  $("dlg").style.display = "none";
  if (dlg.target === "mtzin") $("mtzin").dispatchEvent(
    new Event("change"));
  schedulePreview();
}
$("b_in").onclick   = () => openDlg("mtzin",   "Choose input MTZ",
                                    e => e.mtz);
$("b_out").onclick  = () => openDlg("mtzout",  "Choose output MTZ",
                                    e => e.mtz);
$("b_ncs").onclick  = () => { document.querySelector(
  'input[name=ncs][value=file]').checked = true;
  $("ncs_file").disabled = false; openDlg("ncs_file", "Choose NCS operator file",
                                    e => /\.txt$/i.test(e.name) || e.dir); };
$("b_hist").onclick = () => openDlg("hist_lib", "Choose hist.lib",
                                    e => /hist/i.test(e.name) || e.dir);
$("dlg_up").onclick    = () => loadDir(dirname(dlg.cwd) || dlg.cwd);
$("dlg_home").onclick  = () => loadDir("~");
$("dlg_pick").onclick  = pick;
$("dlg_cancel").onclick = () => $("dlg").style.display = "none";
$("dlg_path").addEventListener("keydown",
  e => { if (e.key === "Enter") pick(); });

// ---------- run / poll -------------------------------------------------------
let polling = null, nlines = 0, cycles = [];
$("run").onclick = async () => {
  const r = await api.run(state());
  $("problems").textContent = "";
  if (!r.ok) { $("problems").textContent = "· " + r.problems.join(" · ");
               return; }
  nlines = 0; cycles = [];
  $("log").textContent = "$ " + r.cmd + "\n";
  drawChart();
  $("run").disabled = true; $("stop").disabled = false;
  setStatus("running…", "");
  pollLoop();
};
$("stop").onclick = async () => { await api.stop(); };

async function pollLoop() {
  const r = await api.poll(nlines);
  if (r.idle) return;
  if (r.lines.length) {
    $("log").textContent += r.lines.join("");
    $("log").scrollTop = $("log").scrollHeight;
    nlines = r.total;
  }
  if (r.cycles.length !== cycles.length) { cycles = r.cycles; drawChart(); }
  if (r.running) {
    setStatus("running… cycle " + (cycles.length ? cycles[cycles.length-1][0] : 0)
              + " · " + Math.round(r.elapsed) + " s", "");
    setTimeout(pollLoop, 600);
  } else {
    $("run").disabled = false; $("stop").disabled = true;
    if (r.rc === 0) setStatus("done in " + Math.round(r.elapsed) + " s", "ok");
    else setStatus("FAILED (exit " + r.rc + ") — see the log", "bad");
  }
}
function setStatus(text, cls) {
  $("status").textContent = text;
  $("status").className = cls;
}

// ---------- free-cc chart (inline SVG) --------------------------------------
function drawChart() {
  const svg = $("chart"), W = svg.clientWidth || 900, H = 120;
  svg.setAttribute("viewBox", "0 0 " + W + " " + H);
  svg.innerHTML = "";
  const ns = "http://www.w3.org/2000/svg";
  const el = (t, attrs, text) => {
    const e = document.createElementNS(ns, t);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text !== undefined) e.textContent = text;
    svg.appendChild(e); return e;
  };
  if (!cycles.length) {
    el("text", {x: W/2, y: H/2, "text-anchor": "middle", fill: "#999",
                "font-size": 12}, "free-set cc per cycle");
    return;
  }
  const pad = 40, padR = 12, padT = 10, padB = 20;
  const ys = cycles.map(c => c[2]).concat([0]);
  let lo = Math.min(...ys), hi = Math.max(...ys);
  if (hi - lo < 1e-6) hi = lo + 1e-6;
  const xy = (i, y) => [
    pad + (W - pad - padR) * (i / Math.max(cycles.length - 1, 1)),
    H - padB - (H - padT - padB) * ((y - lo) / (hi - lo))];
  el("line", {x1: pad-6, y1: H-padB, x2: W-padR, y2: H-padB,
              stroke: "#ccc"});
  el("line", {x1: pad-6, y1: H-padB, x2: pad-6, y2: padT, stroke: "#ccc"});
  el("text", {x: 2, y: padT+8, fill: "#888", "font-size": 10},
     hi.toFixed(2));
  el("text", {x: 2, y: H-padB, fill: "#888", "font-size": 10},
     lo.toFixed(2));
  const pts = cycles.map((c, i) => xy(i, c[2]));
  el("polyline", {points: pts.map(p => p.join(",")).join(" "),
                  fill: "none", stroke: "#08c", "stroke-width": 2});
  cycles.forEach((c, i) => {
    const [x, y] = xy(i, c[2]);
    el("circle", {cx: x, cy: y, r: 3, fill: "#08c"});
    el("text", {x: x, y: H-6, "text-anchor": "middle", fill: "#888",
                "font-size": 9}, c[0]);
  });
}
window.addEventListener("resize", drawChart);
drawChart();
schedulePreview();
</script>
</body>
</html>
"""
