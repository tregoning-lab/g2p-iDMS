"""
Build the interactive mutational-scan site from the tidy tables in data/.

    python build_site.py
    python build_site.py --doi 10.1234/abcd.2026.5678 --citation "Gerardi et al. (2026)"
    python build_site.py --offline --output standalone.html

Every scan listed in datasets.csv becomes one selectable heatmap. The page is a
single self-contained HTML file: pick an antibody and antigen, switch between
fold change and K_d, search for a mutation, hide everything below an effect
threshold, read the sortable table, and download the data.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent

TITLE = "In silico deep mutational scan"
CAPTION = (
    "Every amino acid substitution was modelled at each residue of interest in the "
    "antigen, and binding affinity (K<sub>d</sub>) to the antibody predicted with PRODIGY. "
    "Fold change is K<sub>d</sub> mutant / K<sub>d</sub> wild type, so 10× is binding ten "
    "times weaker than wild type and 0.1× ten times tighter. Both views are logarithmic and "
    "centred on wild type, which sits at the pale midpoint: blue is tighter binding than "
    "wild type, red is weaker. Grey cells are the wild-type residue."
)

# One diverging ramp, blue to red through a neutral midpoint, used for both
# views. Each scan's range is symmetric about wild type, so the pale midpoint
# always means "no change" and the same colour means the same change from wild
# type in every scan. Dark mode has its own steps rather than an automatic
# flip, so that the midpoint sits on the dark surface colour.
BLUE_RED_LIGHT = ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f4a9a8", "#e34948", "#a32d2c"]
BLUE_RED_DARK = ["#cde2fb", "#6da7ec", "#256abf", "#15304f", "#383835",
                 "#4a2a28", "#a32d2c", "#e66767", "#f4a9a8"]

THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#000000", "muted": "#6b6a66", "line": "#9b9a96",
              "hover": "#ffffff", "wt": "#c8c7c2", "dim": "#edecea", "mark": "#000000",
              "ramp": BLUE_RED_LIGHT},
    "dark": {"surface": "#1a1a19", "ink": "#f2f1ee", "muted": "#a3a29d", "line": "#6f6e69",
             "hover": "#262624", "wt": "#6a6963", "dim": "#2a2a27", "mark": "#f2f1ee",
             "ramp": BLUE_RED_DARK},
}


def ramp(colors: list[str]) -> list[list]:
    n = len(colors) - 1
    return [[i / n, c] for i, c in enumerate(colors)]


def build_dataset(spec: dict) -> dict:
    df = pd.read_csv(HERE / "data" / f"{spec['key']}.csv")
    positions = list(dict.fromkeys(df["position"]))
    subs = list(dict.fromkeys(df["substitution"]))

    def grid(column):
        frame = (df.pivot_table(index="substitution", columns="position", values=column, sort=False)
                   .reindex(index=subs, columns=positions))
        return [[None if pd.isna(v) else float(v) for v in row] for row in frame.values]

    fold, kd_M = grid("fold_change"), grid("kd_M")
    kd_nM = [[None if v is None else v * 1e9 for v in row] for row in kd_M]
    wt_mask = [[1 if v is None else None for v in row] for row in fold]

    flat_fold = [v for row in fold for v in row if v is not None]
    limit = math.ceil(max(abs(math.log10(min(flat_fold))), math.log10(max(flat_fold))) * 10) / 10

    return {
        "key": spec["key"],
        "antibody": spec["antibody"],
        "antigen": spec["antigen"],
        "wtKd": float(spec["wt_kd_M"]) * 1e9,
        "x": positions,
        "y": subs,
        "fold": fold,
        "logFold": [[None if v is None else math.log10(v) for v in row] for row in fold],
        "kd": kd_nM,
        "logKd": [[None if v is None else math.log10(v) for v in row] for row in kd_nM],
        "wt": wt_mask,
        "foldLimit": limit,
        "n": len(flat_fold),
        "csv": (HERE / "data" / f"{spec['key']}.csv").read_text(encoding="utf-8"),
    }


PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%%PAGE_TITLE%%</title>
%%PLOTLY%%
<style>
  :root { --surface:#fcfcfb; --page:#f5f5f2; --ink:#1f1f1e; --muted:#6b6a66;
          --line:#e6e5e1; --accent:#256abf; --chip:#ffffff; }
  @media (prefers-color-scheme: dark) {
    :root { --surface:#1a1a19; --page:#121211; --ink:#f2f1ee; --muted:#a3a29d;
            --line:#33332f; --accent:#6da7ec; --chip:#232321; }
  }
  * { box-sizing: border-box; }
  body { margin:0; background:var(--page); color:var(--ink); line-height:1.55;
         font-family: Arial, Helvetica, sans-serif; }
  main { max-width: 1180px; margin: 0 auto; padding: 36px 20px 64px; }
  h1 { font-size: 1.5rem; margin: 0 0 10px; letter-spacing: -0.01em; }
  .caption { color: var(--muted); font-size: .9rem; margin: 0 0 22px; }

  .controls { display:flex; flex-wrap:wrap; gap:16px; align-items:flex-end; margin-bottom:14px; }
  .group { display:flex; flex-direction:column; gap:6px; }
  .group-label { font-size:.7rem; text-transform:uppercase; letter-spacing:.07em; color:var(--muted); }
  select, input[type="search"] {
    appearance:none; font:inherit; font-size:.88rem; color:var(--ink); background:var(--chip);
    border:1px solid var(--line); border-radius:8px; padding:7px 12px; min-height:36px; }
  select { padding-right:30px;
    background-image:linear-gradient(45deg,transparent 50%,currentColor 50%),linear-gradient(135deg,currentColor 50%,transparent 50%);
    background-position:calc(100% - 15px) 16px, calc(100% - 10px) 16px;
    background-size:5px 5px, 5px 5px; background-repeat:no-repeat; }
  input[type="search"] { min-width:190px; }
  select:focus-visible, input:focus-visible, button:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
  .seg { display:flex; border:1px solid var(--line); border-radius:8px; overflow:hidden; background:var(--chip); min-height:36px; }
  .seg button { appearance:none; border:0; background:transparent; color:var(--muted);
                font:inherit; font-size:.88rem; padding:0 15px; cursor:pointer; }
  .seg button + button { border-left:1px solid var(--line); }
  .seg button[aria-pressed="true"] { background:var(--accent); color:#fff; font-weight:600; }
  .spacer { flex:1 1 auto; }
  .btn { appearance:none; border:1px solid var(--line); background:var(--chip); color:var(--ink);
         font:inherit; font-size:.85rem; padding:8px 14px; border-radius:8px; cursor:pointer; min-height:36px; }
  .btn:hover { border-color:var(--accent); color:var(--accent); }

  .meta { color:var(--muted); font-size:.83rem; margin:0 2px 12px; }
  .meta b { color:var(--ink); font-weight:600; }
  .found { color:var(--accent); font-size:.85rem; margin:0 2px 12px; min-height:1.2em; }
  .card { background:var(--surface); border:1px solid var(--line); border-radius:12px; padding:14px 14px 6px; }
  #chart { width:100%; height:%%HEIGHT%%px; }
  .hint { color:var(--muted); font-size:.82rem; margin:12px 4px 0; }

  h2 { font-size:1.05rem; margin:34px 0 10px; }
  .table-wrap { background:var(--surface); border:1px solid var(--line); border-radius:12px;
                overflow:auto; max-height:460px; }
  table { border-collapse:collapse; width:100%; font-size:.84rem; font-variant-numeric:tabular-nums; }
  th, td { padding:7px 14px; border-bottom:1px solid var(--line); text-align:right; white-space:nowrap; }
  th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) { text-align:left; }
  thead th { position:sticky; top:0; background:var(--surface); color:var(--muted);
             font-weight:600; cursor:pointer; user-select:none; z-index:1; }
  thead th:hover { color:var(--accent); }
  thead th[aria-sort] { color:var(--ink); }
  tbody tr:hover { background:color-mix(in srgb, var(--accent) 7%, transparent); }
  .up::after { content:" ▲"; font-size:.7em; }
  .down::after { content:" ▼"; font-size:.7em; }
  footer { color:var(--muted); font-size:.82rem; margin-top:30px; border-top:1px solid var(--line); padding-top:14px; }
  footer a { color:var(--accent); }
</style>
</head>
<body>
<main>
  <h1>%%TITLE%%</h1>
  <p class="caption">%%CAPTION%%</p>

  <div class="controls">
    <div class="group">
      <label class="group-label" for="antibody">Antibody</label>
      <select id="antibody"></select>
    </div>
    <div class="group">
      <label class="group-label" for="antigen">Antigen</label>
      <select id="antigen"></select>
    </div>
    <div class="group">
      <span class="group-label" id="lbl-metric">Shown</span>
      <div class="seg" role="group" aria-labelledby="lbl-metric" id="metric-switch">
        <button type="button" data-metric="fold" aria-pressed="true">Fold change</button>
        <button type="button" data-metric="kd" aria-pressed="false">K<sub>d</sub></button>
      </div>
    </div>
    <div class="group">
      <label class="group-label" for="search">Find a mutation</label>
      <input type="search" id="search" placeholder="e.g. N487G or 487" autocomplete="off">
    </div>
    <div class="group">
      <label class="group-label" for="threshold">Effect at least</label>
      <select id="threshold">
        <option value="1">show everything</option>
        <option value="2">2-fold</option>
        <option value="5">5-fold</option>
        <option value="10">10-fold</option>
        <option value="25">25-fold</option>
        <option value="100">100-fold</option>
      </select>
    </div>
    <span class="spacer"></span>
    <button type="button" class="btn" id="download">Download CSV</button>
  </div>

  <p class="meta" id="meta"></p>
  <p class="found" id="found"></p>

  <div class="card"><div id="chart"></div></div>

  <h2 id="table-heading">All mutations</h2>
  <div class="table-wrap">
    <table>
      <thead><tr>
        <th data-sort="mutation">Mutation</th>
        <th data-sort="position">Position</th>
        <th data-sort="substitution">To</th>
        <th data-sort="kd">K<sub>d</sub> (nM)</th>
        <th data-sort="fold">Fold change</th>
      </tr></thead>
      <tbody id="tbody"></tbody>
    </table>
  </div>

  <footer>
    <!-- DOI: rebuild with  python build_site.py --doi 10.xxxx/xxxxx  or edit the line below. -->
    %%DOI%%
    <div>Predicted with PRODIGY · page built with Python and Plotly.</div>
    <div>Made with Claude Code.</div>
  </footer>
</main>

<script>
const DATA = %%DATA%%;
const THEMES = %%THEMES%%;
const T = (window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches)
          ? THEMES.dark : THEMES.light;

const state = { key: DATA[0].key, metric: "fold", query: "", threshold: 1,
                sort: { col: "fold", dir: -1 } };
const current = () => DATA.find(d => d.key === state.key);

const sig = (v, n) => v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v >= 1 ? v.toFixed(2) : v.toPrecision(n);
const fmtFold = v => v === null ? "—" : sig(v, 2) + "×";
const fmtKd = v => v === null ? "—" : sig(v, 2);
// How far from wild type, as a plain multiple: 0.1x and 10x are both 10-fold.
const effect = v => v === null ? 0 : (v >= 1 ? v : 1 / v);

/* ---------- the rows behind the table and the search ---------- */
function rowsOf(d) {
  const out = [];
  d.y.forEach((aa, i) => d.x.forEach((pos, j) => {
    if (d.fold[i][j] === null) return;
    out.push({ mutation: pos + aa, position: pos, substitution: aa,
               fold: d.fold[i][j], kd: d.kd[i][j], i, j });
  }));
  return out;
}

function matches(row, q) {
  if (!q) return true;
  return row.mutation.toLowerCase().startsWith(q) ||
         row.position.toLowerCase().startsWith(q) ||
         row.position.slice(1).startsWith(q);
}

function visible() {
  const q = state.query.trim().toLowerCase();
  return rowsOf(current()).filter(r => matches(r, q) && effect(r.fold) >= state.threshold);
}

/* ---------- the heatmap ---------- */
function logTicks(lo, hi, suffix) {
  const steps = (hi - lo) < 2.5 ? [1, 3] : [1];
  const vals = [], text = [];
  for (let e = Math.floor(lo); e <= Math.ceil(hi); e++) {
    for (const m of steps) {
      const n = m * Math.pow(10, e), v = Math.log10(n);
      if (v < lo || v > hi) continue;
      vals.push(v);
      text.push((n >= 1 ? String(Math.round(n)) : n.toFixed(n >= 0.1 ? 1 : 2)) + suffix);
    }
  }
  return { vals, text };
}

function traces(shown) {
  const d = current(), isFold = state.metric === "fold";
  const keep = new Set(shown.map(r => r.i + "," + r.j));
  const live = state.threshold > 1 || state.query.trim() !== "";
  const z = (isFold ? d.logFold : d.logKd)
    .map((row, i) => row.map((v, j) => (!live || keep.has(i + "," + j)) ? v : null));
  const dimmed = d.fold.map((row, i) => row.map((v, j) =>
    (live && v !== null && !keep.has(i + "," + j)) ? 1 : null));
  const custom = d.fold.map((row, i) => row.map((v, j) => [fmtFold(v), fmtKd(d.kd[i][j])]));
  // Both views are centred on wild type: 1x for fold change, the wild-type
  // K_d itself for K_d, so the pale midpoint always means "no change".
  const mid = isFold ? 0 : Math.log10(d.wtKd);
  const lo = mid - d.foldLimit, hi = mid + d.foldLimit;
  const ticks = logTicks(lo, hi, isFold ? "×" : "");
  const hover = "<b>%{x} → %{y}</b><br>Fold change: %{customdata[0]}<br>K<sub>d</sub>: %{customdata[1]} nM<extra></extra>";

  const flat = (color, zz, tip) => ({
    type: "heatmap", x: d.x, y: d.y, z: zz, xgap: 2, ygap: 2,
    colorscale: [[0, color], [1, color]], showscale: false, hoverongaps: false,
    customdata: custom, hovertemplate: tip
  });

  return [
    { type: "heatmap", x: d.x, y: d.y, z, xgap: 2, ygap: 2, hoverongaps: false,
      customdata: custom, hovertemplate: hover,
      colorscale: T.ramp,
      zmin: lo, zmax: hi,
      colorbar: {
        title: { text: isFold ? "Fold change in K<sub>d</sub>" : "K<sub>d</sub> (nM)",
                 side: "right", font: { color: T.ink, size: 15 } },
        thickness: 14, outlinewidth: 0, len: 0.92,
        tickfont: { color: T.ink, size: 13 },
        ticks: "outside", ticklen: 5, tickcolor: T.ink,
        tickmode: "array", tickvals: ticks.vals, ticktext: ticks.text }
    },
    flat(T.dim, dimmed, hover),
    flat(T.wt, d.wt, "<b>%{x} → %{y}</b><br>Wild-type residue<extra></extra>")
  ];
}

function outlines(shown) {
  // Ring the cells matching the search, so they are findable in a wall of colour.
  if (!state.query.trim()) return [];
  return shown.slice(0, 60).map(r => ({
    type: "rect", xref: "x", yref: "y",
    x0: r.j - 0.5, x1: r.j + 0.5, y0: r.i - 0.5, y1: r.i + 0.5,
    line: { color: T.mark, width: 2 }, fillcolor: "rgba(0,0,0,0)", layer: "above"
  }));
}

function layout(shown) {
  const axis = { showgrid: false, zeroline: false, type: "category",
                 ticks: "outside", ticklen: 6, tickwidth: 1, tickcolor: T.ink,
                 showline: true, linecolor: T.ink, linewidth: 1,
                 tickfont: { color: T.ink, size: 13 },
                 title: { font: { color: T.ink, size: 16 } } };
  return {
    paper_bgcolor: T.surface, plot_bgcolor: T.surface,
    font: { family: "Arial, Helvetica, sans-serif", color: T.ink },
    margin: { l: 90, r: 16, t: 46, b: 104 },
    xaxis: { ...axis, title: { ...axis.title, text: "Antigen Residue Position", standoff: 34 },
             side: "bottom", tickangle: -45 },
    yaxis: { ...axis, title: { ...axis.title, text: "Amino acid substitution", standoff: 18 },
             autorange: "reversed" },
    hoverlabel: { bgcolor: T.hover, bordercolor: T.line, font: { color: T.ink, size: 13 } },
    shapes: outlines(shown), dragmode: "zoom"
  };
}

/* ---------- the table ---------- */
function renderTable(shown) {
  const { col, dir } = state.sort;
  const sorted = [...shown].sort((a, b) =>
    (typeof a[col] === "string" ? a[col].localeCompare(b[col]) : a[col] - b[col]) * dir);
  document.getElementById("tbody").innerHTML = sorted.map(r =>
    `<tr><td>${r.mutation}</td><td>${r.position}</td><td>${r.substitution}</td>` +
    `<td>${fmtKd(r.kd)}</td><td>${fmtFold(r.fold)}</td></tr>`).join("");
  document.querySelectorAll("thead th").forEach(th => {
    const on = th.dataset.sort === col;
    th.className = on ? (dir > 0 ? "up" : "down") : "";
    if (on) th.setAttribute("aria-sort", dir > 0 ? "ascending" : "descending");
    else th.removeAttribute("aria-sort");
  });
}

/* ---------- the page ---------- */
function fillAntigens() {
  const ab = document.getElementById("antibody").value;
  const sel = document.getElementById("antigen");
  const options = DATA.filter(d => d.antibody === ab);
  sel.innerHTML = options.map(d => `<option value="${d.key}">${d.antigen}</option>`).join("");
  if (!options.some(d => d.key === state.key)) state.key = options[0].key;
  sel.value = state.key;
}

function draw() {
  const d = current(), shown = visible(), all = rowsOf(d).length;

  document.getElementById("meta").innerHTML =
    `<b>${d.antibody}</b> vs <b>${d.antigen}</b> · wild-type K<sub>d</sub> ${sig(d.wtKd, 2)} nM · ` +
    `${d.x.length} positions · ${all} mutations modelled`;

  const q = state.query.trim();
  const note = document.getElementById("found");
  if (q && shown.length === 1) {
    const r = shown[0];
    note.innerHTML = `<b>${r.mutation}</b>: ${fmtFold(r.fold)} · K<sub>d</sub> ${fmtKd(r.kd)} nM`;
  } else if (q && shown.length) {
    note.textContent = `${shown.length} mutations match “${q}”.`;
  } else if (q) {
    note.textContent = `Nothing matches “${q}” in this scan.`;
  } else if (state.threshold > 1) {
    note.textContent = `${shown.length} of ${all} mutations change binding at least ${state.threshold}-fold.`;
  } else {
    note.textContent = "";
  }

  document.getElementById("table-heading").textContent =
    shown.length === all ? "All mutations" : `Mutations shown (${shown.length} of ${all})`;

  Plotly.react("chart", traces(shown), layout(shown), {
    responsive: true, displaylogo: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d", "sendDataToCloud"],
    toImageButtonOptions: { filename: `${state.key}_${state.metric}`, scale: 3 }
  });
  renderTable(shown);

  const hash = new URLSearchParams({ scan: state.key, show: state.metric });
  history.replaceState(null, "", "#" + hash.toString());
}

/* ---------- wiring ---------- */
const antibodySel = document.getElementById("antibody");
antibodySel.innerHTML = [...new Set(DATA.map(d => d.antibody))]
  .map(a => `<option value="${a}">${a}</option>`).join("");

(function restore() {
  const p = new URLSearchParams(location.hash.slice(1));
  const wanted = DATA.find(d => d.key === p.get("scan"));
  if (wanted) state.key = wanted.key;
  if (p.get("show") === "kd") state.metric = "kd";
  antibodySel.value = current().antibody;
  fillAntigens();
  document.querySelectorAll("#metric-switch button").forEach(b =>
    b.setAttribute("aria-pressed", String(b.dataset.metric === state.metric)));
})();

antibodySel.addEventListener("change", () => { fillAntigens(); draw(); });
document.getElementById("antigen").addEventListener("change", e => { state.key = e.target.value; draw(); });
document.getElementById("metric-switch").addEventListener("click", e => {
  const b = e.target.closest("button"); if (!b) return;
  state.metric = b.dataset.metric;
  e.currentTarget.querySelectorAll("button").forEach(x =>
    x.setAttribute("aria-pressed", String(x.dataset.metric === state.metric)));
  draw();
});
document.getElementById("search").addEventListener("input", e => { state.query = e.target.value; draw(); });
document.getElementById("threshold").addEventListener("change", e => { state.threshold = +e.target.value; draw(); });
document.querySelectorAll("thead th").forEach(th => th.addEventListener("click", () => {
  const col = th.dataset.sort;
  state.sort = { col, dir: state.sort.col === col ? -state.sort.dir : (col === "fold" || col === "kd" ? -1 : 1) };
  renderTable(visible());
}));
document.getElementById("download").addEventListener("click", () => {
  const d = current();
  const url = URL.createObjectURL(new Blob([d.csv], { type: "text/csv" }));
  const a = document.createElement("a");
  a.href = url; a.download = `${d.key}.csv`; a.click();
  URL.revokeObjectURL(url);
});

draw();
</script>
</body>
</html>
"""


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", "-o", type=Path, default=HERE / "docs" / "index.html")
    p.add_argument("--title", default=TITLE)
    p.add_argument("--doi", default="", help="DOI of the paper, e.g. 10.1234/abcd.2026.5678")
    p.add_argument("--citation", default="", help="citation line shown next to the DOI")
    p.add_argument("--offline", action="store_true",
                   help="embed Plotly in the file so it works with no internet connection")
    args = p.parse_args()

    manifest = list(csv.DictReader((HERE / "datasets.csv").open()))
    data = [build_dataset(spec) for spec in manifest]

    if args.offline:
        import plotly
        js = (Path(plotly.__file__).parent / "package_data" / "plotly.min.js").read_text(encoding="utf-8")
        plotly_tag = f"<script>{js}</script>"
    else:
        plotly_tag = '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>'

    bits = []
    if args.citation:
        bits.append(html.escape(args.citation))
    if args.doi:
        doi = args.doi.replace("https://doi.org/", "").strip()
        bits.append(f'<a href="https://doi.org/{html.escape(doi)}">doi:{html.escape(doi)}</a>')

    page = PAGE
    for token, value in {
        "%%PAGE_TITLE%%": html.escape(args.title),
        "%%TITLE%%": html.escape(args.title),
        "%%CAPTION%%": CAPTION,
        "%%PLOTLY%%": plotly_tag,
        "%%HEIGHT%%": str(max(520, 120 + len(data[0]["y"]) * 26)),
        "%%DOI%%": f'<div>{" · ".join(bits)}</div>' if bits else "",
        "%%DATA%%": json.dumps(data),
        "%%THEMES%%": json.dumps({k: {**v, "ramp": ramp(v["ramp"])}
                                  for k, v in THEMES.items()}),
    }.items():
        page = page.replace(token, value)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(page, encoding="utf-8")
    print(f"Wrote {args.output}  ({args.output.stat().st_size / 1024:,.0f} KB, "
          f"{len(data)} scans: {', '.join(d['key'] for d in data)})")


if __name__ == "__main__":
    main()
