"""Reviewed, endpoint-independent bead gating on documented linear FCS values.

Polygon coordinates are expressed after the configured display/gating transform;
reported APC medians always use untransformed analysis values. Shared polygon
boundaries and overlapping gates are deliberately ambiguous, never first-match.
This module provides no automatic biological gate or count threshold.

This first version supports Fortessa-style pre-gate axes named FSC-*, SSC-*,
and Time (case insensitive). Bead gates require FSC versus SSC; singlet gates
require area versus height/width for the same scatter detector. Arbitrarily
renamed detector axes require a separately reviewed implementation, not an alias
that could accidentally select on the APC endpoint.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


class GateError(ValueError):
    """An acquisition/gating policy cannot be applied safely."""


def transform(values, spec=None):
    """Apply a linear or asinh gating transform without changing source events."""
    spec = {} if spec is None else spec
    if not isinstance(spec, dict):
        raise GateError("Transform specifications must be objects")
    kind = spec.get("kind", "linear")
    values = np.asarray(values, dtype=float)
    if kind == "linear":
        return values.copy()
    if kind == "asinh":
        cofactor = spec.get("cofactor")
        if isinstance(cofactor, bool) or not isinstance(cofactor, (float, int)) or not np.isfinite(cofactor) or cofactor <= 0:
            raise GateError("asinh transforms require a finite positive cofactor")
        return np.arcsinh(values / cofactor)
    raise GateError(f"Unsupported transform kind: {kind!r}")


def _orientation(a, b, c):
    delta1, delta2 = b - a, c - a
    return delta1[0] * delta2[1] - delta1[1] * delta2[0]


def _on_segment(a, b, c, tolerance):
    return abs(_orientation(a, b, c)) <= tolerance and np.all(c >= np.minimum(a, b) - tolerance) and np.all(c <= np.maximum(a, b) + tolerance)


def _segments_intersect(a, b, c, d, tolerance):
    o1, o2 = _orientation(a, b, c), _orientation(a, b, d)
    o3, o4 = _orientation(c, d, a), _orientation(c, d, b)
    if ((o1 > tolerance and o2 < -tolerance) or (o1 < -tolerance and o2 > tolerance)) and ((o3 > tolerance and o4 < -tolerance) or (o3 < -tolerance and o4 > tolerance)):
        return True
    return (_on_segment(a, b, c, tolerance) or _on_segment(a, b, d, tolerance) or _on_segment(c, d, a, tolerance) or _on_segment(c, d, b, tolerance))


def validate_polygon(vertices):
    """Return a finite simple polygon, allowing a repeated closing vertex."""
    try:
        polygon = np.asarray(vertices, dtype=float)
    except (TypeError, ValueError) as exc:
        raise GateError("Polygon coordinates must be numeric") from exc
    if polygon.ndim != 2 or polygon.shape[1] != 2 or not np.isfinite(polygon).all():
        raise GateError("Polygon coordinates must be finite [x, y] pairs")
    if len(polygon) > 1 and np.array_equal(polygon[0], polygon[-1]):
        polygon = polygon[:-1]
    if len(polygon) < 3 or len(np.unique(polygon, axis=0)) != len(polygon):
        raise GateError("A polygon needs at least three distinct vertices")
    scale = max(1.0, float(np.ptp(polygon, axis=0).max()))
    tolerance = np.finfo(float).eps * scale * scale * 32
    area = float(np.sum(polygon[:, 0] * np.roll(polygon[:, 1], -1) - polygon[:, 1] * np.roll(polygon[:, 0], -1)))
    if abs(area) <= tolerance:
        raise GateError("Polygon has zero or numerically degenerate area")
    n = len(polygon)
    for i in range(n):
        for j in range(i + 1, n):
            if j == i + 1 or (i == 0 and j == n - 1):
                continue
            if _segments_intersect(polygon[i], polygon[(i + 1) % n], polygon[j], polygon[(j + 1) % n], tolerance):
                raise GateError("Self-intersecting polygons are not supported")
    return polygon


def polygon_mask(points, vertices):
    """Return point membership, including edges; caller resolves overlaps."""
    points = np.asarray(points, dtype=float)
    polygon = validate_polygon(vertices)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise GateError("Points must be finite [x, y] pairs")
    x, y = points[:, 0], points[:, 1]
    inside = np.zeros(len(points), dtype=bool)
    boundary = np.zeros(len(points), dtype=bool)
    scale = max(1.0, float(np.ptp(polygon, axis=0).max()))
    tolerance = np.finfo(float).eps * scale * scale * 32
    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        delta = b - a
        cross = delta[0] * (y - a[1]) - delta[1] * (x - a[0])
        boundary |= (np.abs(cross) <= tolerance) & (x >= min(a[0], b[0]) - tolerance) & (x <= max(a[0], b[0]) + tolerance) & (y >= min(a[1], b[1]) - tolerance) & (y <= max(a[1], b[1]) + tolerance)
        crossing = (a[1] > y) != (b[1] > y)
        if delta[1] != 0:
            x_intercept = a[0] + (y - a[1]) * delta[0] / delta[1]
            inside ^= crossing & (x < x_intercept)
    return inside | boundary


def _channel_indices(data):
    if len(set(data.channels)) != len(data.channels):
        raise GateError("FCS channel names must be unique")
    return {name: i for i, name in enumerate(data.channels)}


def _required_channels(data, setup):
    if not isinstance(setup, dict):
        raise GateError("Acquisition setup must be an object")
    channels = setup.get("channels", {})
    if not isinstance(channels, dict):
        raise GateError("Configured channels must be an object")
    names = [channels.get(role) for role in ("fitc", "pe", "apc")]
    if any(not isinstance(name, str) or not name for name in names) or len(set(names)) != 3:
        raise GateError("Specify three distinct FITC, PE and APC channel names")
    indices = _channel_indices(data)
    missing = [name for name in names if name not in indices]
    if missing:
        raise GateError(f"FCS is missing configured channels: {missing}")
    return channels, indices


def _embedded_matrix(data):
    raw = data.metadata.get("spillover") or data.metadata.get("spill")
    if not raw:
        raise GateError("Embedded compensation requested but no SPILLOVER/SPILL metadata exists")
    if not isinstance(raw, str):
        raise GateError("Embedded compensation metadata must be a string")
    parts = [s.strip() for s in raw.split(",")]
    try:
        n = int(parts[0])
        if n < 1 or len(parts) != 1 + n + n * n:
            raise ValueError
        return parts[1:1 + n], np.asarray(parts[1 + n:], dtype=float).reshape(n, n)
    except (ValueError, IndexError) as exc:
        raise GateError("Malformed embedded compensation matrix") from exc


def prepare_events(data, setup):
    """Validate explicit compensation policy and return linear analysis events.

    Matrix convention: measured = true @ spillover. Channel order is declared
    by the matrix, and may differ from FCS parameter order. Already compensated
    exports must use mode=none. Embedded metadata alone never implies processing.
    """
    _, indices = _required_channels(data, setup)
    values = np.asarray(data.events, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(data.channels) or not np.isfinite(values).all():
        raise GateError("Event array must be finite and match FCS channel metadata")
    values = values.copy()
    policy = setup.get("compensation", {})
    if not isinstance(policy, dict):
        raise GateError("Compensation policy must be an object")
    mode, state = policy.get("mode"), policy.get("input_state")
    if not isinstance(state, str) or state not in {"uncompensated", "compensated"}:
        raise GateError("Declare compensation.input_state as uncompensated or compensated from acquisition records")
    if not isinstance(mode, str) or mode not in {"none", "embedded", "explicit"}:
        raise GateError("Declare compensation.mode as none, embedded or explicit")
    if state == "compensated" and mode != "none":
        raise GateError("Already compensated inputs cannot be compensated a second time")
    if mode == "none":
        if not isinstance(policy.get("reason"), str) or not policy["reason"].strip():
            raise GateError("Compensation mode none requires a documented reason")
        return values
    if mode == "embedded":
        names, matrix = _embedded_matrix(data)
    else:
        names = policy.get("channels", [])
        try:
            matrix = np.asarray(policy.get("matrix", []), dtype=float)
        except (TypeError, ValueError) as exc:
            raise GateError("Compensation matrix must be numeric") from exc
    if not isinstance(names, (list, tuple)) or not names or not all(isinstance(c, str) for c in names) or len(set(names)) != len(names):
        raise GateError("Compensation matrix channels must be distinct names")
    if any(c not in indices for c in names):
        raise GateError("Compensation matrix references channels absent from FCS")
    if any(name.upper().startswith(("FSC-", "SSC-")) or name.upper() == "TIME" for name in names):
        raise GateError("Fluorescence compensation cannot modify scatter or Time channels")
    n = len(names)
    if matrix.shape != (n, n) or not np.isfinite(matrix).all() or not np.allclose(np.diag(matrix), 1.0, atol=1e-8, rtol=0):
        raise GateError("Compensation requires a finite square spillover matrix with unit diagonal")
    if np.linalg.cond(matrix) > 1e10:
        raise GateError("Compensation matrix is singular or numerically unstable")
    cols = [indices[name] for name in names]
    try:
        values[:, cols] = np.linalg.solve(matrix.T, values[:, cols].T).T
    except np.linalg.LinAlgError as exc:
        raise GateError("Compensation matrix cannot be solved") from exc
    if not np.isfinite(values).all():
        raise GateError("Compensation produced nonfinite values")
    return values


def gate_events(data, setup, panel, *, stages=None):
    """Return compensated events, sequential pre-gate mask and memberships.

    Memberships include pre-gate acceptance, but do not resolve overlap. The
    public helper supports audit plots using exactly the analysis definitions.
    If supplied, ``stages`` receives cumulative event counts after each gate;
    the return value remains the same three-item tuple.
    """
    channels, indices = _required_channels(data, setup)
    if setup.get("reviewed") is not True or any(not isinstance(setup.get(field), str) or not setup[field].strip() for field in ("gate_version", "review_note")):
        raise GateError("A reviewed gate template with gate_version and review_note is required")
    if stages is not None and not isinstance(stages, list):
        raise GateError("stages must be a list when supplied")
    minimum = setup.get("min_events_per_population")
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
        raise GateError("Set a positive integer min_events_per_population from the validated acquisition protocol")
    values = prepare_events(data, setup)
    pregates = setup.get("pre_gates", [])
    if not isinstance(pregates, list):
        raise GateError("pre_gates must be a list")
    if any(not isinstance(g, dict) for g in pregates):
        raise GateError("Every pre-gate must be an object")
    gate_names = [g.get("name") for g in pregates]
    if any(not isinstance(name, str) or not name.strip() for name in gate_names):
        raise GateError("Every pre-gate needs a nonempty name")
    if len(gate_names) != len(set(gate_names)) or not {"beads", "singlets"}.issubset(gate_names):
        raise GateError("Distinct pre_gates named beads and singlets are required")
    mask = np.ones(len(values), dtype=bool)
    for gate in pregates:
        names = [gate.get("x"), gate.get("y")]
        if any(not isinstance(name, str) or not name for name in names):
            raise GateError("Pre-gate axes must be nonempty channel names")
        if any(name not in indices for name in names):
            raise GateError(f"Pre-gate {gate.get('name')} references an absent channel")
        if names[0] == names[1]:
            raise GateError("Pre-gates must use two distinct channels")
        if any(name == channels["apc"] or "apc" in name.casefold() for name in names):
            raise GateError("APC detection signal cannot be used to select beads or singlets")
        conventional = [name.upper() for name in names]
        if any(not (name.startswith("FSC-") or name.startswith("SSC-") or name == "TIME") for name in conventional):
            raise GateError("Pre-gate axes must use conventional FSC-*, SSC-* or Time names")
        if gate["name"] == "beads" and {name.split("-", 1)[0] for name in conventional} != {"FSC", "SSC"}:
            raise GateError("The beads gate requires FSC versus SSC")
        if gate["name"] == "singlets":
            valid_pulses = any(
                set(conventional) in ({f"{detector}-A", f"{detector}-H"}, {f"{detector}-A", f"{detector}-W"})
                for detector in ("FSC", "SSC")
            )
            if not valid_pulses:
                raise GateError("The singlets gate requires same-detector scatter area versus height or width")
        points = np.column_stack([transform(values[:, indices[names[0]]], gate.get("x_transform")), transform(values[:, indices[names[1]]], gate.get("y_transform"))])
        mask &= polygon_mask(points, gate.get("vertices", []))
        if stages is not None:
            stages.append({"gate": gate["name"], "events": int(mask.sum())})
    populations = setup.get("populations", [])
    if not isinstance(populations, list):
        raise GateError("populations must be a list")
    if any(not isinstance(population, dict) for population in populations):
        raise GateError("Every population gate must be an object")
    population_ids = [str(pop.get("bead_id", "")) for pop in populations]
    panel_ids = list(panel["bead_id"].astype(str))
    if len(set(panel_ids)) != len(panel_ids):
        raise GateError("Panel bead IDs must be unique")
    if len(set(population_ids)) != len(population_ids) or set(population_ids) != set(panel_ids):
        raise GateError("Population gates must map to every panel bead ID exactly once")
    space = setup.get("bead_space", {})
    if not isinstance(space, dict):
        raise GateError("bead_space must be an object")
    points = np.column_stack([transform(values[:, indices[channels["fitc"]]], space.get("x_transform")), transform(values[:, indices[channels["pe"]]], space.get("y_transform"))])
    membership = {str(pop["bead_id"]): polygon_mask(points, pop.get("vertices", [])) & mask for pop in populations}
    return values, mask, membership


def summarize(data, setup, panel):
    """Summarize every panel bead ID, preserving failed measurements for audit."""
    stages = []
    values, pregated, membership = gate_events(data, setup, panel, stages=stages)
    apc_index = data.channels.index(setup["channels"]["apc"])
    counts = np.sum(np.column_stack(list(membership.values())), axis=1)
    ambiguous = pregated & (counts > 1)
    threshold = setup.get("apc_saturation_threshold")
    if threshold is not None:
        if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not np.isfinite(threshold) or threshold <= 0:
            raise GateError("apc_saturation_threshold must be finite and positive")
    else:
        upper = data.ranges.get(setup["channels"]["apc"])
        if upper is None or not np.isfinite(upper) or upper <= 1:
            raise GateError("APC saturation requires a valid parameter range or explicit apc_saturation_threshold")
        threshold = upper - 1
    is_original = setup["compensation"]["input_state"] == "uncompensated"
    review_note = setup.get("saturation_review_note")
    manually_reviewed = isinstance(review_note, str) and bool(review_note.strip())
    saturation_unverified = not is_original and not manually_reviewed
    saturated = np.zeros(len(values), dtype=bool)
    if is_original:
        # Saturation is detected before any computational compensation. Negative
        # values are not interpreted as a lower clipping limit. For already
        # compensated exports these original values are unavailable: only an
        # explicit review of the acquisition QC can certify saturation status.
        saturated = np.asarray(data.events)[:, apc_index] >= threshold
    rows = []
    for _, pop in panel.iterrows():
        bead_id = str(pop["bead_id"])
        selected = membership[bead_id] & ~ambiguous
        n = int(selected.sum())
        flags = []
        if n == 0:
            flags.append("missing_population")
        elif n < setup["min_events_per_population"]:
            flags.append("low_bead_count")
        if np.any(selected & saturated):
            flags.append("apc_saturation")
        if np.any(membership[bead_id] & ambiguous):
            flags.append("ambiguous_population_events")
        if saturation_unverified:
            flags.append("acquisition_saturation_unverified")
        rows.append({"bead_id": bead_id, "marker": pop["marker"], "role": pop["role"], "n_events": n, "raw_median_apc": float(np.median(values[selected, apc_index])) if n else np.nan, "measurement_valid": not flags, "qc_flags": ";".join(flags)})
    qc = {"total_events": len(values), "pregated_events": int(pregated.sum()), "assigned_events": int((pregated & (counts == 1)).sum()), "ambiguous_events": int(ambiguous.sum()), "unassigned_events": int((pregated & (counts == 0)).sum()), "apc_saturated_events": int((pregated & saturated).sum()) if is_original else None, "apc_saturation_threshold": float(threshold) if is_original else None, "apc_saturation_source": "original_fcs" if is_original else "manual_acquisition_review" if manually_reviewed else "unverified", "acquisition_saturation_unverified": saturation_unverified, "gate_counts": stages}
    return pd.DataFrame(rows), qc
