#!/usr/bin/env python3
"""Derive and name the 39 MACSPlex EV Neuro bead populations from a setup-bead FCS file.

The script writes an UNREVIEWED setup template for ``config/analysis.json`` together
with its evidence: a per-population table and an annotated FITC-versus-PE overlay.
Identification uses FITC/PE geometry only. Density peaks of gated singlet beads are
matched to the manufacturer's Figure 2 staircase (``identity.FIGURE2_LEVELS``) and
named from Table 1 (``config.EXPECTED_MARKERS``). APC is never read for identity, a
missing population is never filled, and the template is never marked reviewed: a
person must inspect the overlay, set ``min_events_per_population``, ``gate_version``
and ``review_note``, and only then set ``reviewed`` to true.

Example (Fortessa export, asinh display with cofactor 150)::

    .venv/bin/python scripts/identify_setup_beads.py data/raw/20261029_Beads_5000.fcs \
        --setup-id fortessa_20261029 --output results/identify-20261029
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from facs_neuro.config import EXPECTED_MARKERS  # noqa: E402
from facs_neuro.fcs import read_fcs  # noqa: E402
from facs_neuro.gating import polygon_mask, summarize, transform  # noqa: E402
from facs_neuro.identity import FIGURE2_LEVELS, FIGURE2_ROWS, audit_identity  # noqa: E402
from facs_neuro.pipeline import new_output  # noqa: E402


def robust_interval(values: np.ndarray, k: float = 6.0) -> list[float]:
    median = float(np.median(values))
    sigma = 1.4826 * float(np.median(np.abs(values - median))) or 1.0
    return [median - k * sigma, median + k * sigma]


def rectangle(x_interval, y_interval) -> list[list[float]]:
    (x0, x1), (y0, y1) = x_interval, y_interval
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]


def density_peaks(points: np.ndarray, sigma: float, bins: int = 400) -> np.ndarray:
    """Local maxima of a smoothed 2-D histogram, in data coordinates.

    ``sigma`` is the Gaussian smoothing in histogram bins. Peaks below a tenth of
    the highest density are debris, not bead populations, and are discarded.
    """
    low, high = points.min(axis=0), points.max(axis=0)
    pad = 0.05 * (high - low)
    hist, xe, ye = np.histogram2d(points[:, 0], points[:, 1], bins=bins,
                                  range=[[low[0] - pad[0], high[0] + pad[0]], [low[1] - pad[1], high[1] + pad[1]]])
    reach = int(np.ceil(3 * sigma))
    kernel = np.exp(-0.5 * (np.arange(-reach, reach + 1) / sigma) ** 2)
    kernel /= kernel.sum()
    smooth = np.apply_along_axis(lambda row: np.convolve(row, kernel, mode="same"), 0, hist)
    smooth = np.apply_along_axis(lambda row: np.convolve(row, kernel, mode="same"), 1, smooth)
    radius = 5
    padded = np.pad(smooth, radius, mode="constant", constant_values=-1)
    windows = np.lib.stride_tricks.sliding_window_view(padded, (2 * radius + 1, 2 * radius + 1))
    peaks = (smooth == windows.max(axis=(2, 3))) & (smooth > 0.1 * smooth.max())
    ix, iy = np.nonzero(peaks)
    return np.column_stack([(xe[ix] + xe[ix + 1]) / 2, (ye[iy] + ye[iy + 1]) / 2])


def typical_gap(centers: np.ndarray) -> float:
    gaps = np.linalg.norm(centers[:, None] - centers[None], axis=2)
    np.fill_diagonal(gaps, np.inf)
    return float(np.median(gaps.min(axis=1)))


def cluster(points: np.ndarray, centers: np.ndarray, keep_fraction: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """Assign points to the nearest centre within a fraction of the typical inter-centre gap."""
    cutoff = keep_fraction * typical_gap(centers)
    labels = np.full(len(points), -1)
    for _ in range(3):
        distance = np.linalg.norm(points[:, None] - centers[None], axis=2)
        labels = distance.argmin(axis=1)
        labels[distance.min(axis=1) > cutoff] = -1
        centers = np.array([np.median(points[labels == k], axis=0) if (labels == k).sum() >= 5 else centers[k]
                            for k in range(len(centers))])
    sizes = np.bincount(labels[labels >= 0], minlength=len(centers))
    keep = np.flatnonzero(sizes >= 10)
    remap = {old: new for new, old in enumerate(keep)}
    return centers[keep], np.array([remap.get(label, -1) for label in labels])


def merge_split_peaks(points: np.ndarray, centers: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Merge peaks closer than 1.5 median within-population SDs: one elongated population, two peaks."""
    spread = [np.sqrt(np.linalg.eigvalsh(np.cov(points[labels == k].T)).max()) for k in range(len(centers))]
    threshold = 1.5 * float(np.median(spread))
    sizes = np.bincount(labels[labels >= 0], minlength=len(centers))
    kept: list[int] = []
    for k in np.argsort(-sizes):
        if all(np.linalg.norm(centers[k] - centers[j]) > threshold for j in kept):
            kept.append(int(k))
    return centers[sorted(kept)]


def hungarian(cost: np.ndarray) -> dict[int, int]:
    """Minimum-cost one-to-one assignment of rows to columns (rows <= columns)."""
    n, m = cost.shape
    u, v, p, way = [0.0] * (n + 1), [0.0] * (m + 1), [0] * (m + 1), [0] * (m + 1)
    for i in range(1, n + 1):
        p[0], j0 = i, 0
        minv, used = [np.inf] * (m + 1), [False] * (m + 1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], np.inf, 0
            for j in range(1, m + 1):
                if not used[j]:
                    current = cost[i0 - 1, j - 1] - u[i0] - v[j]
                    if current < minv[j]:
                        minv[j], way[j] = current, j0
                    if minv[j] < delta:
                        delta, j1 = minv[j], j
            for j in range(m + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]
            p[j0], j0 = p[j1], j1
    return {p[j] - 1: j - 1 for j in range(1, m + 1) if p[j]}


def decode_lattice(centers: np.ndarray) -> tuple[dict[str, int], dict[str, float], dict[str, float]]:
    """Match centres to Figure 2 code levels through an iteratively refitted affine map."""
    ids = sorted(FIGURE2_LEVELS)
    uv = np.array([FIGURE2_LEVELS[b] for b in ids], dtype=float)
    dimmest = int(np.argmin(centers.sum(axis=1)))
    distance = np.linalg.norm(centers - centers[dimmest], axis=1)
    distance[dimmest] = np.inf
    first, second = np.argsort(distance)[:2]
    angle = [np.arctan2(*(centers[j] - centers[dimmest])[::-1]) for j in (first, second)]
    right, up = (first, second) if angle[0] < angle[1] else (second, first)

    def fit(source, target):
        return np.linalg.lstsq(np.column_stack([source, np.ones(len(source))]), target, rcond=None)[0]

    coefficients = fit(np.array([[2, 2], [3, 2], [2, 3]], float), centers[[dimmest, right, up]])
    assignment: dict[str, int] = {}
    for _ in range(20):
        predicted = np.column_stack([uv, np.ones(len(uv))]) @ coefficients
        cost = np.linalg.norm(predicted[:, None] - centers[None], axis=2)
        if len(centers) >= len(ids):
            assignment = {ids[row]: column for row, column in hungarian(cost).items()}
        else:
            assignment = {ids[row]: column for column, row in hungarian(cost.T).items()}
        rows = [ids.index(b) for b in assignment]
        refit = fit(uv[rows], centers[[assignment[ids[r]] for r in rows]])
        if np.allclose(refit, coefficients):
            break
        coefficients = refit
    residual = {b: float(cost[ids.index(b), k]) for b, k in assignment.items()}
    margin = {}
    for b, k in assignment.items():
        others = np.delete(cost[ids.index(b)], k)
        margin[b] = float(others.min() - cost[ids.index(b), k])
    return assignment, residual, margin


def clipped_ellipse(points: np.ndarray, centre_index: int, centers: np.ndarray, scale: float = 3.5) -> np.ndarray:
    """A 3.5-SD robust ellipse, cut back to 45% of the way towards every other centre."""
    mu = np.median(points, axis=0)
    covariance = np.cov(points.T)
    mahalanobis = np.sqrt(np.einsum("ij,jk,ik->i", points - mu, np.linalg.inv(covariance), points - mu))
    core = points[mahalanobis <= 2.5]
    mu, covariance = np.median(core, axis=0), np.cov(core.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    angles = np.linspace(0, 2 * np.pi, 17)[:-1]
    polygon = (eigenvectors @ (np.sqrt(eigenvalues)[:, None] * np.vstack([np.cos(angles), np.sin(angles)]) * scale)).T + mu
    for other, centre in enumerate(centers):
        if other == centre_index:
            continue
        normal = centre - mu
        cut = mu + 0.45 * normal
        kept = []
        for i in range(len(polygon)):
            a, b = polygon[i], polygon[(i + 1) % len(polygon)]
            ta, tb = float(np.dot(a - cut, normal)), float(np.dot(b - cut, normal))
            if ta <= 0:
                kept.append(a)
            if (ta <= 0) != (tb <= 0):
                kept.append(a + (b - a) * (ta / (ta - tb)))
        polygon = np.array(kept)
    return polygon


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("fcs")
    parser.add_argument("--output", required=True, help="New directory for template, table and overlay")
    parser.add_argument("--setup-id", default="setup")
    parser.add_argument("--panel", default=str(Path(__file__).resolve().parents[1] / "config" / "panel.csv"))
    parser.add_argument("--fitc", default="FITC-A")
    parser.add_argument("--pe", default="PE-A")
    parser.add_argument("--apc", default="APC-A")
    parser.add_argument("--transform", choices=("asinh", "linear"), default="asinh")
    parser.add_argument("--cofactor", type=float, default=150.0)
    parser.add_argument("--smoothing", type=float, default=2.5,
                        help="Density smoothing in histogram bins (1/400 of the range); raise if populations split")
    parser.add_argument("--singlet-axis", default=None, help="FSC-W (default when present) or FSC-H")
    parser.add_argument("--beads-gate", help="JSON FSC-A/SSC-A polygon; default: robust rectangle")
    parser.add_argument("--singlets-gate", help="JSON FSC-A/singlet-axis polygon; default: robust band")
    args = parser.parse_args(argv)

    data = read_fcs(args.fcs)
    index = {name: i for i, name in enumerate(data.channels)}
    for name in (args.fitc, args.pe, args.apc, "FSC-A", "SSC-A"):
        if name not in index:
            raise SystemExit(f"Channel {name} is absent; channels are {data.channels}")
    singlet_axis = args.singlet_axis or ("FSC-W" if "FSC-W" in index else "FSC-H")
    if singlet_axis not in index:
        raise SystemExit(f"Singlet axis {singlet_axis} is absent")
    events = data.events
    fsc, ssc, pulse = events[:, index["FSC-A"]], events[:, index["SSC-A"]], events[:, index[singlet_axis]]

    fsc_interval = robust_interval(fsc)
    beads_gate = json.loads(args.beads_gate) if args.beads_gate else rectangle(fsc_interval, robust_interval(ssc))
    beads = polygon_mask(np.column_stack([fsc, ssc]), beads_gate)
    if args.singlets_gate:
        singlets_gate = json.loads(args.singlets_gate)
    elif singlet_axis.upper().endswith("-W"):
        singlets_gate = rectangle(fsc_interval, robust_interval(pulse[beads]))
    else:
        low, high = robust_interval(pulse[beads] / fsc[beads], k=5.0)
        x0, x1 = fsc_interval
        singlets_gate = [[x0, low * x0], [x1, low * x1], [x1, high * x1], [x0, high * x0]]
    singlets = beads & polygon_mask(np.column_stack([fsc, pulse]), singlets_gate)

    spec = {"kind": "asinh", "cofactor": args.cofactor} if args.transform == "asinh" else {"kind": "linear"}
    x = transform(events[singlets, index[args.fitc]], spec)
    y = transform(events[singlets, index[args.pe]], spec)
    points = np.column_stack([x, y])
    peaks, labels = cluster(points, density_peaks(points, args.smoothing))
    peaks, labels = cluster(points, merge_split_peaks(points, peaks, labels))
    problems = []
    if len(peaks) != 39:
        problems.append(f"{len(peaks)} populations resolved after peak merging; 39 are expected (try --smoothing)")
    out = new_output(args.output)
    if len(peaks) < 3:
        raise SystemExit("Fewer than three populations; check the pre-gates and transform")

    assignment, residual, margin = decode_lattice(peaks)
    # Keep only the matched peaks; events of any extra (split) peak rejoin the nearest matched one.
    order = sorted(assignment, key=lambda b: FIGURE2_LEVELS[b])
    centers = peaks[[assignment[b] for b in order]]
    assignment = {b: i for i, b in enumerate(order)}
    centers, labels = cluster(points, centers)
    neighbours = []
    for bead, k in assignment.items():
        distance = np.linalg.norm(centers - centers[k], axis=1)
        distance[k] = np.inf
        nearest = int(distance.argmin())
        other = next((b for b, j in assignment.items() if j == nearest), None)
        if other is None or max(abs(int(bead[0]) - int(other[0])), abs(int(bead[1]) - int(other[1]))) > 1:
            neighbours.append(f"{bead} nearest cluster is not a Figure 2 neighbour")
    problems += neighbours

    linear = {name: events[singlets, index[name]] for name in (args.fitc, args.pe, args.apc)}
    rows = []
    polygons = {}
    for bead in sorted(assignment, key=lambda b: FIGURE2_LEVELS[b]):
        member = labels == assignment[bead]
        polygons[bead] = clipped_ellipse(points[member], assignment[bead], centers)
        rows.append({"bead_id": bead, "marker": EXPECTED_MARKERS[bead],
                     "figure2_fitc_level": FIGURE2_LEVELS[bead][0], "figure2_pe_level": FIGURE2_LEVELS[bead][1],
                     "n_events": int(member.sum()),
                     "fitc_median": float(np.median(linear[args.fitc][member])),
                     "pe_median": float(np.median(linear[args.pe][member])),
                     "apc_median": float(np.median(linear[args.apc][member])),
                     "lattice_residual": round(residual[bead], 3), "assignment_margin": round(margin[bead], 3)})
    table = pd.DataFrame(rows)
    lookup = table.set_index("bead_id")
    for row in FIGURE2_ROWS:
        for a, b in zip(row, row[1:], strict=False):
            if a in lookup.index and b in lookup.index and not lookup.loc[a, "fitc_median"] < lookup.loc[b, "fitc_median"]:
                problems.append(f"FITC order {a} < {b} violated")
    for level in range(2, 10):
        column = sorted((b for b in FIGURE2_LEVELS if FIGURE2_LEVELS[b][0] == level), key=lambda b: FIGURE2_LEVELS[b][1])
        for a, b in zip(column, column[1:], strict=False):
            if a in lookup.index and b in lookup.index and not lookup.loc[a, "pe_median"] < lookup.loc[b, "pe_median"]:
                problems.append(f"PE order {a} < {b} violated")

    setup = {
        "reviewed": False, "gate_version": "", "review_note": "",
        "channels": {"fitc": args.fitc, "pe": args.pe, "apc": args.apc},
        "compensation": {"input_state": "uncompensated", "mode": "none", "reason": ""},
        "bead_space": {"x_transform": spec, "y_transform": spec},
        "pre_gates": [
            {"name": "beads", "x": "FSC-A", "y": "SSC-A", "vertices": np.round(beads_gate, 2).tolist()},
            {"name": "singlets", "x": "FSC-A", "y": singlet_axis, "vertices": np.round(singlets_gate, 2).tolist()}],
        "min_events_per_population": None,
        "populations": [{"bead_id": b, "vertices": np.round(polygons[b], 4).tolist()} for b in sorted(polygons)],
        "derived_from": {"file": Path(args.fcs).name, "sha256": data.sha256, "script": Path(__file__).name,
                         "instrument": data.metadata.get("cyt", ""), "acquired": data.metadata.get("date", "")},
    }
    # Dry-run the pipeline's own gating and audit with the template to prove it reproduces the result.
    pipeline_check = {}
    if len(assignment) == 39:
        panel = pd.read_csv(args.panel, dtype=str, keep_default_na=False)
        trial = {**setup, "reviewed": True, "gate_version": "dry-run", "review_note": "dry-run", "min_events_per_population": 1,
                 "compensation": {**setup["compensation"], "reason": "dry-run; original detector values"}}
        measured, qc = summarize(data, trial, panel)
        audit = audit_identity(data, trial, panel)
        pipeline_check = {"ambiguous_events": qc["ambiguous_events"], "unassigned_pregated_events": qc["unassigned_events"],
                          "populations_valid": int(measured.measurement_valid.sum()),
                          "audit_consistent": int(audit.layout_status.eq("consistent").sum())}
        if qc["ambiguous_events"] or not measured.measurement_valid.all() or not audit.layout_status.eq("consistent").all():
            problems.append("pipeline dry-run flagged populations; see summary.json")

    table.to_csv(out / "populations.csv", index=False)
    (out / "setup-template.json").write_text(json.dumps({args.setup_id: setup}, indent=1) + "\n")
    summary = {"file": str(Path(args.fcs).resolve()), "sha256": data.sha256, "total_events": int(len(events)),
               "beads": int(beads.sum()), "singlets": int(singlets.sum()), "density_peaks": int(len(peaks)),
               "populations_resolved": int(len(centers)),
               "unassigned_singlets": int((labels == -1).sum()), "pipeline_dry_run": pipeline_check,
               "problems": problems, "status": "identified_unreviewed" if not problems else "review_required"}
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    axes[0].scatter(fsc, ssc, s=2, alpha=.3, c="gray")
    axes[0].scatter(fsc[beads], ssc[beads], s=2, alpha=.4)
    axes[0].plot(*np.vstack([beads_gate, beads_gate[0]]).T, c="#b45309")
    axes[0].set(xlabel="FSC-A", ylabel="SSC-A", title=f"beads: {int(beads.sum())} of {len(events)}")
    axes[1].scatter(fsc[beads], pulse[beads], s=2, alpha=.3, c="gray")
    axes[1].scatter(fsc[singlets], pulse[singlets], s=2, alpha=.4)
    axes[1].plot(*np.vstack([singlets_gate, singlets_gate[0]]).T, c="#b45309")
    axes[1].set(xlabel="FSC-A", ylabel=singlet_axis, title=f"singlets: {int(singlets.sum())}")
    axes[2].scatter(x, y, s=3, alpha=.25, c="lightgray")
    for bead, k in assignment.items():
        member = labels == k
        axes[2].scatter(x[member], y[member], s=3, alpha=.6)
        axes[2].plot(*np.vstack([polygons[bead], polygons[bead][0]]).T, lw=.7, c="k")
        axes[2].annotate(f"{bead}\n{EXPECTED_MARKERS[bead]}", centers[k] + 0.03 * np.ptp(points, axis=0), fontsize=6.5)
    label = f"{args.transform}" + (f" cofactor {args.cofactor:g}" if args.transform == "asinh" else "")
    axes[2].set(xlabel=f"{args.fitc} ({label})", ylabel=f"{args.pe} ({label})",
                title=f"{len(assignment)} populations decoded against Figure 2, named from Table 1 — UNREVIEWED")
    fig.tight_layout()
    fig.savefig(out / "identified.png", dpi=110)
    print(json.dumps(summary, indent=1))
    print(table.to_string(index=False))
    return 0 if not problems else 2


if __name__ == "__main__":
    raise SystemExit(main())
