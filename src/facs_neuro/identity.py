"""Necessary, not sufficient, Figure 2 consistency checks for bead assignments.

Miltenyi MACSPlex EV Kit Neuro manual 140-006-716.02, printed p.7, shows
FITC increasing with the first bead-ID digit and PE with the second. The
schematic supplies ordinal levels, not Fortessa fluorescence coordinates.

This audit uses observed FITC/PE medians of uniquely assigned, pre-gated events.
It can expose an axis reversal or inconsistent labels; passing cannot establish
true bead identity, separation, suitable bead counts, or assay validity. It never
relabels events, fills a missing population, or uses APC to choose an identity.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import EXPECTED_MARKERS
from .gating import GateError, gate_events


# Figure 2, bottom to top, left to right. These are the full 39 IDs, including
# two isotypes. The diagram does not imply equal intensity spacing or alignment.
FIGURE2_ROWS = (
    ("22", "32", "42", "52"),
    ("23", "33", "43", "53", "63"),
    ("24", "34", "44", "54", "64", "74"),
    ("35", "45", "55", "65", "75", "85"),
    ("46", "56", "66", "76", "86", "96"),
    ("57", "67", "77", "87", "97"),
    ("68", "78", "88", "98"),
    ("79", "89", "99"),
)
FIGURE2_LEVELS = {
    bead: (int(bead[0]), int(bead[1]))
    for row in FIGURE2_ROWS for bead in row
}


def _reference_table(panel: pd.DataFrame) -> pd.DataFrame:
    if not {"bead_id", "marker", "role"}.issubset(panel.columns):
        raise GateError("Identity audit requires bead_id, marker and role columns")
    rows = panel.loc[:, ["bead_id", "marker", "role"]].copy()
    rows["bead_id"] = rows["bead_id"].astype(str)
    if len(rows) != 39 or rows["bead_id"].duplicated().any() or set(rows["bead_id"]) != set(FIGURE2_LEVELS):
        raise GateError("Identity audit requires all 39 canonical Figure 2 bead IDs")
    for _, row in rows.iterrows():
        bead = row["bead_id"]
        expected_role = "isotype_control" if bead in {"97", "99"} else "capture_marker"
        if row["marker"] != EXPECTED_MARKERS[bead] or row["role"] != expected_role:
            raise GateError(f"Canonical Table 1 identity mismatch for bead {bead}")
    rows["figure2_fitc_level"] = rows["bead_id"].map(lambda b: FIGURE2_LEVELS[b][0])
    rows["figure2_pe_level"] = rows["bead_id"].map(lambda b: FIGURE2_LEVELS[b][1])
    return rows


def audit_identity(data, setup: dict, panel: pd.DataFrame) -> pd.DataFrame:
    """Audit only adjacent within-row FITC and within-column PE order.

    Medians use linear analysis values after the declared compensation policy,
    not transformed gate coordinates or polygon centroids. Missing or ambiguous
    assignments remain absent. Counts and separation require their own QC.
    ``layout_status=consistent`` means this ordinal check passed, not that the
    instrument's biological bead assignments have been independently validated.
    """
    result = _reference_table(panel).reset_index(drop=True)
    values, _, membership = gate_events(data, setup, panel)
    assigned_counts = np.sum(np.column_stack(list(membership.values())), axis=1)
    fitc_index = data.channels.index(setup["channels"]["fitc"])
    pe_index = data.channels.index(setup["channels"]["pe"])
    counts, fitc, pe = [], [], []
    flags = {bead: [] for bead in result["bead_id"]}
    for bead in result["bead_id"]:
        selected = membership[bead] & (assigned_counts == 1)
        count = int(selected.sum())
        counts.append(count)
        fitc.append(float(np.median(values[selected, fitc_index])) if count else np.nan)
        pe.append(float(np.median(values[selected, pe_index])) if count else np.nan)
        if not count:
            flags[bead].append("missing_population")
    result["n_events"] = counts
    result["observed_fitc_median"] = fitc
    result["observed_pe_median"] = pe
    lookup = result.set_index("bead_id")

    def check_adjacent(beads, column, axis):
        for lower, upper in zip(beads, beads[1:]):
            a, b = lookup.loc[lower, column], lookup.loc[upper, column]
            if not np.isfinite(a) or not np.isfinite(b):
                continue
            if not a < b:
                flags[lower].append(f"figure2_{axis}_order_vs_{upper}")
                flags[upper].append(f"figure2_{axis}_order_vs_{lower}")

    for row in FIGURE2_ROWS:
        check_adjacent(row, "observed_fitc_median", "fitc")
    for level in range(2, 10):
        column = sorted((b for b in FIGURE2_LEVELS if FIGURE2_LEVELS[b][0] == level), key=lambda b: FIGURE2_LEVELS[b][1])
        check_adjacent(column, "observed_pe_median", "pe")
    result["layout_status"] = [
        "missing" if "missing_population" in flags[b]
        else "inconsistent" if flags[b] else "consistent"
        for b in result["bead_id"]
    ]
    result["layout_flags"] = [";".join(flags[b]) for b in result["bead_id"]]
    return result


def write_identity_reference(output_dir: str | Path, panel: pd.DataFrame) -> None:
    """Write the canonical Figure 2/Table 1 reference, never numerical gates.

    Both filenames must be unused. The plot shows explicit ordinal code levels;
    its coordinates must not be transferred to instrument data or gate files.
    """
    reference = _reference_table(panel).sort_values("bead_id").reset_index(drop=True)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "bead-identity-reference.csv"
    png_path = out / "bead-identity-reference.png"
    if csv_path.exists() or png_path.exists():
        raise FileExistsError("Bead identity reference output already exists")
    reference.to_csv(csv_path, index=False)

    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    fig = Figure(figsize=(12, 9))
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    fig.subplots_adjust(left=.08, right=.99, top=.91, bottom=.13)
    for row in reference.to_dict("records"):
        x, y = row["figure2_fitc_level"], row["figure2_pe_level"]
        color = "#fff0ce" if row["role"] == "isotype_control" else "#e7eff9"
        ax.add_patch(Rectangle((x - .46, y - .37), .92, .74, facecolor=color, edgecolor="#607186", linewidth=.8))
        label = row["marker"].replace(" control", "\ncontrol")
        if label == "Ganglioside GD2":
            label = "Ganglioside\nGD2"
        ax.text(x, y, f"{row['bead_id']}\n{label}", ha="center", va="center", fontsize=8)
    ax.set(xlim=(1.35, 9.65), ylim=(1.25, 9.75), xticks=range(2, 10), yticks=range(2, 10),
           xlabel="MACSPlex-FITC code level (first digit; ordinal)",
           ylabel="MACSPlex-PE code level (second digit; ordinal)",
           title="Schematic ordinal levels — NOT Fortessa fluorescence coordinates\nMACSPlex EV Kit Neuro: Figure 2 layout with Table 1 markers")
    ax.grid(alpha=.15)
    ax.set_axisbelow(True)
    fig.text(.5, .018, "37 capture markers + 2 isotype controls (orange). Reference: 140-006-716.02, printed pp.7 and 10.\nObserved setup-bead clusters and assignments still require independent review.", ha="center", fontsize=8)
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.clear()
