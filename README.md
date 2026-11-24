# g2p-heatmaps

Interactive heatmaps of computational mutational scans of antibody–antigen binding.

Every amino acid substitution was modelled at each residue of interest in the antigen, and binding affinity (K<sub>d</sub>) to the antibody predicted with [PRODIGY](https://rascar.science.uu.nl/prodigy/). Fold change is K<sub>d</sub> mutant / K<sub>d</sub> wild type, so 10× is binding ten times weaker than wild type and 0.1× ten times tighter.

![Heatmap preview](assets/preview.png)

## The scans

| Antibody | Antigen | Positions | Mutations | Wild-type K<sub>d</sub> |
|---|---|---|---|---|
| Nirsevimab | RSV-A F | 27 | 513 | 2.0 nM |
| Nirsevimab | RSV-B F | 27 | 513 | 1.5 nM |
| RB1 (Clesrovimab) | RSV F | 23 | 437 | 3.3 nM |
| Casirivimab | SARS-CoV-2 Delta RBD | 23 | 437 | 6.5 nM |

## The page

Pick an antibody and antigen from the dropdowns, then:

- switch between **fold change** and **K<sub>d</sub>**
- **find a mutation** — type `N487G` for one cell, or `487` for a whole position; matches are ringed on the heatmap and the table narrows to them
- **effect at least** — hide everything that changes binding by less than 2-, 5-, 10-, 25- or 100-fold, in either direction, to see only the residues that matter
- read the **sortable table** of every mutation, by fold change, K<sub>d</sub>, position or substitution
- hover any cell for its exact values, drag to zoom, save a PNG, download the CSV

The address bar tracks which scan you are looking at, so a link points at a particular view.

Open `docs/index.html` in a browser (it loads Plotly from a CDN, so it needs an internet connection). For a single file that works offline and can be emailed as one attachment:

```bash
python build_site.py --offline --output mutational-scans.html
```

## Rebuilding

Requires Python 3.9+.

```bash
pip install -r requirements.txt
python prepare_data.py     # data/source/*.csv  ->  data/<scan>.csv
python build_site.py       # data/<scan>.csv    ->  docs/index.html
```

### Adding a DOI

Once the paper is out, rebuild with the DOI and it appears in the page footer:

```bash
python build_site.py --doi 10.1234/abcd.2026.5678 --citation "Gerardi et al. (2026)"
```

### Adding another scan

1. Put the raw PRODIGY output in `data/source/`.
2. Add one row to `datasets.csv`.
3. Run `python prepare_data.py && python build_site.py`.

The new antibody and antigen appear in the dropdowns automatically; nothing else needs editing.

`datasets.csv` columns:

| Column | Meaning |
|---|---|
| `key` | file name stem for `data/<key>.csv` |
| `antibody`, `antigen` | names shown in the dropdowns |
| `source` | path to the raw file |
| `format` | `long` (one mutation per row, with a fold change), `long_log` (one per row, fold change as a log<sub>10</sub>) or `grid` (a matrix) |
| `wt_kd_M` | wild-type K<sub>d</sub> in molar — the reference for every fold change |
| `relabel` | optional `OLD=NEW` position fixes |

## Data

`data/source/` holds the raw output exactly as produced. `data/<scan>.csv` holds the tidy table the site is built from: position, substitution, wild-type flag, K<sub>d</sub> (M), fold change, log<sub>10</sub> fold change. The wild-type residue at each position carries no values.

`prepare_data.py` checks every scan before writing it, and stops rather than produce a figure that is quietly wrong:

- the gap in each column must fall on that position's own wild-type residue
- any fold change quoted in the source must match K<sub>d</sub> / K<sub>d</sub>(wild type)

All 1,881 fold changes across the four scans pass.

The clesrovimab source is the `Figure6_E` sheet of the supplied workbook, exported as a CSV. The workbook's other sheets (contact distances, GISAID accessions with isolation dates, locations and patient ages, and the publication tables) are not part of the scan and are not included here.

### One correction applied

In the raw RSV-B file, position 197 is labelled **S197** and the substitution list skips S. Residue 197 is N in both RSV subtypes, so the list should skip N. The 19 predicted values are in alphabetical order of substitution and are unaffected, but treating S as wild type shifted the values for N, P, Q and R by one row. The `relabel` entry in `datasets.csv` moves the gap back to N, which reproduces the original source spreadsheet exactly. Nothing else in any file was changed.

The RSV wild-type residue differs between subtypes at five of the 27 positions: 67 (N in A, T in B), 200 (D / N), 201 (K / N), 207 (V / L) and 209 (K / Q).

## Project structure

```
g2p-heatmaps/
├── datasets.csv         # the manifest: one row per scan
├── prepare_data.py      # raw output -> tidy CSVs, with checks
├── build_site.py        # tidy CSVs -> docs/index.html
├── requirements.txt
├── data/
│   ├── source/          # raw output, unmodified
│   └── *.csv            # one tidy table per scan
├── docs/index.html      # the site
└── assets/preview.png
```

## Colour

Both views use the same blue–white–red diverging ramp on a log axis: **blue is tighter binding than wild type, red is weaker**, and the pale midpoint is wild type itself — 1× for fold change, the wild-type K<sub>d</sub> for K<sub>d</sub>.

Each scan's range is symmetric about that midpoint, so ten times weaker and ten times tighter sit the same distance either side, and the same colour means the same change from wild type in every scan. The two views therefore show the same pattern in different units, which is expected: fold change is K<sub>d</sub> divided by a constant.

Grey cells mark the wild-type residue, and mutations hidden by the effect filter stay visible as faint cells rather than disappearing. Light and dark mode use their own steps rather than one being a flip of the other, so the midpoint always sits on the surface colour.
