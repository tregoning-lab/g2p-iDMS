"""
Turn the raw PRODIGY output listed in datasets.csv into one tidy table per scan.

    python prepare_data.py

Every scan is one row of `datasets.csv`:

    key        file name stem used for data/<key>.csv
    antibody   antibody name, shown in the site's filter
    antigen    antigen name, shown in the site's filter
    source     path to the raw file
    format     long     - one mutation per row: Mutation, Predicted KD (M), Kd FC
               long_log - one mutation per row, with the fold change given as
                          a log10: Mutation, PRODIGY Kd M, Log10(Mut Kd/ WT Kd)
               grid     - a matrix: positions across the top, substitutions
                          down the side, optionally followed by a block of
                          log10 fold changes, which is checked but not used
    wt_kd_M    wild-type K_d in molar, the reference for every fold change
    relabel    optional OLD=NEW position fixes, separated by spaces

Each scan is checked before it is written: the gap in every column must fall on
that position's own wild-type residue, and any fold change in the source must
match K_d / K_d(wild type). The script stops rather than write a table that
fails either check.

Known correction
----------------
In the RSV-B source, position 197 is labelled S197 and the substitution list
skips S. Residue 197 is N in both RSV subtypes, so the list should skip N. The
19 predicted values are in alphabetical order of substitution and are
unaffected, but treating S as wild type shifted the values for N, P, Q and R by
one row. `relabel = S197=N197` moves the gap back to N, which reproduces the
original source spreadsheet exactly.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
AA = "ACDEFGHIKLMNPQRSTVWY"
TOLERANCE = 0.02          # relative, on fold changes quoted to 3-4 figures


# --- readers ---------------------------------------------------------------
def read_long(path: Path) -> tuple[dict[str, dict[str, float]], dict[str, float]]:
    """One mutation per row: Mutation, Predicted KD (M), Kd FC."""
    kd: dict[str, dict[str, float]] = {}
    quoted: dict[str, float] = {}
    for row in csv.DictReader(path.open()):
        mutation = row["Mutation"].strip()
        position, substitution = mutation[:-1], mutation[-1]
        kd.setdefault(position, {})[substitution] = float(row["Predicted KD (M)"])
        if row.get("Kd FC", "").strip():
            quoted[mutation] = float(row["Kd FC"])
    return kd, quoted


def read_long_log(path: Path) -> tuple[dict[str, dict[str, float]], dict[str, float]]:
    """One mutation per row, fold change given as log10. A 'WT' row is skipped."""
    kd: dict[str, dict[str, float]] = {}
    quoted: dict[str, float] = {}
    for row in csv.DictReader(path.open()):
        mutation = row["Mutation"].strip()
        if mutation.upper() == "WT":
            continue
        position, substitution = mutation[:-1], mutation[-1]
        kd.setdefault(position, {})[substitution] = float(row["PRODIGY Kd M"])
        log = row.get("Log10(Mut Kd/ WT Kd)", "").strip()
        if log:
            quoted[mutation] = 10 ** float(log)
    return kd, quoted


def read_grid(path: Path) -> tuple[dict[str, dict[str, float]], dict[str, float]]:
    """A matrix of K_d, optionally followed by a block of log10 fold changes."""
    rows = list(csv.reader(path.open()))
    positions = [p.strip() for p in rows[0][1:] if p.strip()]

    kd: dict[str, dict[str, float]] = {p: {} for p in positions}
    log_block: list[list[str]] = []
    seen = []
    for row in rows[1:]:
        label = row[0].strip()
        if len(label) == 1 and label in AA and label not in seen:   # a K_d row
            seen.append(label)
            for p, v in zip(positions, row[1:]):
                if v.strip():
                    kd[p][label] = float(v)
        elif not label and any(c.strip() for c in row[1:]):   # a log fold-change row
            log_block.append(row[1:])

    if len(seen) != 20:
        raise SystemExit(f"{path.name}: found {len(seen)} substitution rows, expected 20")

    # The trailing block is a redundant copy; use it as a check, not as data.
    quoted: dict[str, float] = {}
    for aa, row in zip(seen, log_block):
        for p, v in zip(positions, row):
            v = v.strip()
            if v and v != "#NUM!":
                quoted[f"{p}{aa}"] = 10 ** float(v)
    return kd, quoted


READERS = {"long": read_long, "long_log": read_long_log, "grid": read_grid}


# --- checks ----------------------------------------------------------------
def relabel(subs: dict[str, float], new_wt: str) -> dict[str, float]:
    """Move the wild-type gap to `new_wt`, keeping the values in their order."""
    values = [subs[a] for a in AA if a in subs]
    return dict(zip([a for a in AA if a != new_wt], values))


def check(key: str, kd: dict[str, dict[str, float]], quoted: dict[str, float],
          wt_kd: float) -> None:
    for position, subs in kd.items():
        gap = set(AA) - set(subs)
        if gap != {position[0]}:
            raise SystemExit(
                f"{key} {position}: wild type should be {position[0]}, "
                f"but the gap is at {sorted(gap) or 'nothing'}")
    bad = [(m, f, kd[m[:-1]][m[-1]] / wt_kd) for m, f in quoted.items()
           if abs(kd[m[:-1]][m[-1]] / wt_kd - f) > TOLERANCE * max(f, 1e-9)]
    if bad:
        lines = "\n".join(f"    {m}: source says {f:.4g}x, K_d implies {c:.4g}x"
                          for m, f, c in bad[:10])
        raise SystemExit(f"{key}: {len(bad)} fold changes disagree with K_d\n{lines}")


# --- main ------------------------------------------------------------------
def main() -> None:
    manifest = list(csv.DictReader((HERE / "datasets.csv").open()))
    for spec in manifest:
        key, wt_kd = spec["key"], float(spec["wt_kd_M"])
        kd, quoted = READERS[spec["format"]](HERE / spec["source"])

        for pair in spec.get("relabel", "").split():
            old, new = pair.split("=")
            kd[new] = relabel(kd.pop(old), new[0])
            quoted = {m: f for m, f in quoted.items() if not m.startswith(old)}

        check(key, kd, quoted, wt_kd)
        positions = sorted(kd, key=lambda p: int(p[1:]))

        rows = []
        for aa in AA:
            for position in positions:
                value = kd[position].get(aa)
                if value is None:
                    rows.append([position, aa, "yes", "", "", ""])
                else:
                    fc = value / wt_kd
                    rows.append([position, aa, "no", f"{value:.3g}",
                                 f"{fc:.4g}", f"{math.log10(fc):.4g}"])

        out = HERE / "data" / f"{key}.csv"
        with out.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["position", "substitution", "wild_type", "kd_M",
                             "fold_change", "log10_fold_change"])
            writer.writerows(rows)
        print(f"{out.relative_to(HERE)}: {len(positions)} positions, "
              f"{sum(r[2] == 'no' for r in rows)} mutations, "
              f"wild-type Kd {wt_kd * 1e9:g} nM"
              + (f", {len(quoted)} fold changes verified" if quoted else ""))


if __name__ == "__main__":
    main()
