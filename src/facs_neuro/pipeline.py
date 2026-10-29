"""Run orchestration, provenance and gate diagnostics."""

from __future__ import annotations

from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil

import numpy as np
import pandas as pd

from . import __version__
from .config import ConfigurationError, load_config, load_inputs, readiness_errors, file_hash, select_samples
from .fcs import read_fcs
from .gating import gate_events, prepare_events, summarize, transform
from .identity import audit_identity, write_identity_reference
from .annotations import ANNOTATION_VERSION, annotate_measurements, load_annotations, load_sources
from .quantify import quantify


def save_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def new_output(path: str | Path) -> Path:
    path = Path(path).resolve()
    if path.exists():
        raise ConfigurationError(f"Output already exists; choose a new directory: {path}")
    path.mkdir(parents=True)
    return path


def inspect(config_path: str | Path) -> dict:
    config = load_config(config_path)
    all_samples, controls, panel = load_inputs(config)
    samples, excluded = select_samples(config, all_samples)
    errors = readiness_errors(config, samples, controls)
    inventory = []
    for row in samples + controls:
        entry = {"id": row.get("acquisition_id", row.get("control_id")), "path": row["_file"]}
        if row["_file"] and Path(row["_file"]).is_file():
            try:
                data = read_fcs(row["_file"])
                entry.update(events=len(data.events), channels=list(data.channels), sha256=data.sha256,
                             metadata=data.metadata, ranges=data.ranges)
            except (ValueError, OSError) as exc:
                entry["error"] = str(exc)
                errors.append(f"{entry['id']}: {exc}")
        else:
            entry["error"] = "FCS path not supplied or file missing"
        inventory.append(entry)
    # Gate review is separate from file/metadata readiness.
    gate_ready = all(s.get("reviewed") is True for s in config["setups"].values())
    exclusions = [{"sample_id": r["sample_id"], "acquisition_id": r["acquisition_id"],
                   "exclusion_reason": r["exclusion_reason"]} for r in excluded]
    inventory.extend({"id": r["acquisition_id"], "path": r["_file"], "excluded": True,
                      "exclusion_reason": r["exclusion_reason"]} for r in excluded)
    return {"dataset_kind": config["dataset_kind"], "acquisitions": len(samples),
            "input_acquisitions": len(all_samples), "excluded_acquisitions": len(excluded),
            "sample_exclusions": exclusions,
            "donors": len({r.get('subject_id') for r in samples if r.get('subject_id')}),
            "capture_markers": int((panel.role == "capture_marker").sum()),
            "isotype_controls": int((panel.role == "isotype_control").sum()), "inputs_ready": not errors,
            "gates_marked_reviewed": gate_ready, "errors": errors, "inventory": inventory}


def diagnostic_plot(data, setup: dict, panel: pd.DataFrame, path: Path, *, reviewed: bool) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    values = prepare_events(data, setup)
    ch = setup["channels"]
    space = setup.get("bead_space", {})
    x = transform(values[:, data.channels.index(ch["fitc"])], space.get("x_transform"))
    y = transform(values[:, data.channels.index(ch["pe"])], space.get("y_transform"))
    subset = np.linspace(0, len(values) - 1, min(len(values), 20000), dtype=int) if len(values) else np.empty(0, dtype=int)
    pregates = setup.get("pre_gates", []) if reviewed else []
    panels = 1 + len(pregates) if pregates else 2
    fig, axes = plt.subplots(1, panels, figsize=(6 * panels, 5))
    axes[0].scatter(x[subset], y[subset], s=1, alpha=.3, c="gray")
    if reviewed:
        _, _, membership = gate_events(data, setup, panel)
        unique = np.sum(np.column_stack(list(membership.values())), axis=1) == 1
        for bead_id, mask in membership.items():
            indices = np.flatnonzero(mask & unique)
            if len(indices):
                select = indices[np.linspace(0, len(indices) - 1, min(len(indices), 400), dtype=int)]
                axes[0].scatter(x[select], y[select], s=2)
        for pop in setup["populations"]:
            polygon = np.array(pop["vertices"])
            closed = np.vstack([polygon, polygon[0]])
            axes[0].plot(closed[:, 0], closed[:, 1], linewidth=.6)
            axes[0].text(*polygon.mean(axis=0), pop["bead_id"], fontsize=6)
    axes[0].set(xlabel=f"{ch['fitc']} ({space.get('x_transform', {}).get('kind', 'linear')})",
                ylabel=f"{ch['pe']} ({space.get('y_transform', {}).get('kind', 'linear')})",
                title="Bead identities: all events in gray")
    cumulative = np.ones(len(values), dtype=bool)
    if pregates:
        from .gating import polygon_mask
        for axis, gate in zip(axes[1:], pregates, strict=True):
            gx = transform(values[:, data.channels.index(gate["x"])], gate.get("x_transform"))
            gy = transform(values[:, data.channels.index(gate["y"])], gate.get("y_transform"))
            cumulative &= polygon_mask(np.column_stack([gx, gy]), gate["vertices"])
            axis.scatter(gx[subset], gy[subset], s=1, alpha=.2, c="gray")
            accepted = np.asarray(subset)[cumulative[subset]]
            axis.scatter(gx[accepted], gy[accepted], s=1, alpha=.5, c="#2563eb")
            polygon = np.array(gate["vertices"])
            closed = np.vstack([polygon, polygon[0]])
            axis.plot(closed[:, 0], closed[:, 1], c="#b45309", linewidth=1)
            axis.set_xlabel(gate["x"], parse_math=False)
            axis.set_ylabel(gate["y"], parse_math=False)
            axis.set_title(f"{gate['name']}: {int(cumulative.sum())} accepted", parse_math=False)
    elif "FSC-A" in data.channels and "SSC-A" in data.channels:
        axes[1].scatter(values[subset, data.channels.index("FSC-A")],
                        values[subset, data.channels.index("SSC-A")], s=1, alpha=.3)
        axes[1].set(xlabel="FSC-A (linear)", ylabel="SSC-A (linear)", title="Scatter overview")
    else:
        axes[1].text(.1, .5, "Conventional FSC-A / SSC-A unavailable")
        axes[1].set_axis_off()
    fig.suptitle("Reviewed gate template" if reviewed else "UNREVIEWED overview — not a gate calibration")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def calibrate(config_path: str | Path, output: str | Path) -> Path:
    """Produce reference overviews and an unreviewed template, never auto-label clusters."""
    config = load_config(config_path)
    _, controls, panel = load_inputs(config)
    references = [r for r in controls if r["control_type"] == "setup_beads"]
    if not references:
        raise ConfigurationError("Add setup_beads controls before preparing gate overviews")
    prepared = []
    for row in references:
        if row.get("setup_id") not in config["setups"] or not row["_file"]:
            raise ConfigurationError(f"Setup control {row['control_id']} needs file and valid setup_id")
        if not Path(row["_file"]).is_file():
            raise ConfigurationError(f"Setup control {row['control_id']}: FCS file missing")
        data = read_fcs(row["_file"])
        setup = config["setups"][row["setup_id"]]
        values = prepare_events(data, setup)
        space = setup.get("bead_space", {})
        if not isinstance(space, dict):
            raise ConfigurationError("bead_space must be an object")
        for role, axis in (("fitc", "x"), ("pe", "y")):
            transform(values[:, data.channels.index(setup["channels"][role])], space.get(f"{axis}_transform"))
        prepared.append((row, data))
    out = new_output(output)
    try:
        write_identity_reference(out, panel)
        templates = {}
        for index, (row, data) in enumerate(prepared):
            setup = config["setups"][row["setup_id"]]
            diagnostic_plot(data, setup, panel, out / f"reference-{index + 1:03d}.png", reviewed=False)
            template = json.loads(json.dumps(setup))
            template.update(reviewed=False, gate_version="", review_note="",
                            populations=[{"bead_id": b, "vertices": []} for b in panel.bead_id])
            templates[row["setup_id"]] = template
        save_json(out / "gate-templates.json", templates)
        (out / "README.txt").write_text(
            "These overviews and empty polygons are NOT calibrated gates. Review scatter, singlets, "
            "the manufacturer's bead identity map, compensation and the chosen transform/cofactor. "
            "Fill polygons and pilot-derived QC thresholds; copy the reviewed setup into your config. "
            "Use bead-identity-reference.png/csv to check the Figure 2 layout and Table 1 names. "
            "These reference levels are ordinal, not Fortessa fluorescence coordinates. "
            "The sample FCS, assay-blank FCS and setup-bead FCS are the required input types; "
            "a separate supplier gate/workspace file is not required. Derived gates must be "
            "checked against the setup-bead data before being marked reviewed.\n")
    except Exception as exc:
        save_json(out / "failure.json", {"status": "incomplete", "error": str(exc)})
        raise
    return out


def _instrument_signature(data, setup: dict) -> tuple:
    signature = []
    names = [setup["channels"][c] for c in ("fitc", "pe", "apc")]
    gates = setup.get("pre_gates", [])
    if not isinstance(gates, list):
        raise ConfigurationError("pre_gates must be a list")
    for gate in gates:
        if not isinstance(gate, dict):
            raise ConfigurationError("Every pre-gate must be an object")
        for axis in ("x", "y"):
            if not isinstance(gate.get(axis), str) or not gate[axis].strip():
                raise ConfigurationError("Pre-gate axes must be nonempty channel names")
            names.append(gate[axis])
    for name in dict.fromkeys(names):
        if name not in data.channels:
            raise ConfigurationError(f"Required detector {name} absent in FCS")
        i = data.channels.index(name) + 1
        vals = []
        for suffix in ("r", "v", "e", "g", "f", "l"):
            raw = data.metadata.get(f"p{i}{suffix}", "")
            try:
                value = tuple(float(v) for v in raw.split(",")) if raw else ()
            except ValueError:
                value = (raw,)
            vals.append(value)
        signature.append((name, *vals))
    # Different compensation matrices under one setup are not interchangeable.
    signature.append(("spillover", data.metadata.get("spillover", data.metadata.get("spill", ""))))
    signature.append(("instrument", data.metadata.get("cyt", ""), data.metadata.get("cytsn", "")))
    return tuple(signature)


def analyze(config_path: str | Path, output: str | Path) -> Path:
    config = load_config(config_path)
    all_samples, controls, panel = load_inputs(config)
    samples, excluded = select_samples(config, all_samples)
    errors = readiness_errors(config, samples, controls)
    if errors:
        raise ConfigurationError("Inputs are incomplete:\n" + "\n".join(errors[:30]) +
                                 (f"\n... {len(errors) - 30} additional issues; run inspect" if len(errors) > 30 else ""))
    normal = config.get("normalization", {})
    background = config.get("background", {})
    background_mode = background.get("mode", "assay_blank")
    background_control_type = "setup_beads" if background_mode == "setup_baseline" else "assay_blank"
    background_ids = {r["blank_id"] for r in samples}
    enabled = normal.get("enabled", False)
    if not isinstance(enabled, bool):
        raise ConfigurationError("normalization.enabled must be true or false")
    floor = normal.get("minimum_factor")
    if enabled and (isinstance(floor, bool) or not isinstance(floor, (float, int)) or not np.isfinite(floor) or floor < 0):
        raise ConfigurationError("Normalization requires an explicit nonnegative minimum_factor validated on pilot data")
    markers = normal.get("markers", ["CD9", "CD63", "CD81"])
    if not isinstance(markers, list) or not markers or not all(isinstance(m, str) for m in markers) or len(set(markers)) != len(markers) or not set(markers) <= {"CD9", "CD63", "CD81"}:
        raise ConfigurationError("Declare a unique, fixed set of CD9/CD63/CD81 normalizers")
    if enabled:
        for row in samples:
            blocked = {v.strip() for v in row.get("affinity_target", "").replace("+", ",").split(",")}
            if blocked & set(markers):
                raise ConfigurationError(f"Sample {row['sample_id']}: affinity enrichment blocks a selected normalizer; review the fixed set")
    out = new_output(output)
    (out / "diagnostics").mkdir()
    provenance = {
        "dataset_kind": config["dataset_kind"],
        "validation_status": "synthetic_only" if config["dataset_kind"] == "synthetic" else "unvalidated_assay",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_sha256": config["_sha256"], "python": platform.python_version(),
        "software": {"facs-neuro": __version__, **{
            x: importlib.metadata.version(x) for x in ("numpy", "pandas", "matplotlib")}},
        "manifest_hashes": {k: file_hash(Path(config[k])) for k in ("samples", "controls", "panel")},
        "configuration_source": config["_path"],
        "source_hashes": {p.name: file_hash(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
        "inputs": [], "gate_versions": {k: v.get("gate_version") for k, v in config["setups"].items()},
        "statistics": "Descriptive per-acquisition summaries; no group tests or automatic technical-replicate aggregation",
        "background": {
            "mode": background_mode,
            "rationale": background.get("rationale", ""),
            "limitation": "setup_baseline_only_reagent_background_unmeasured" if background_mode == "setup_baseline" else "",
        },
        "sample_exclusions": [{"sample_id": r["sample_id"], "acquisition_id": r["acquisition_id"],
                               "exclusion_reason": r["exclusion_reason"]} for r in excluded],
    }
    tables, blank_tables, control_tables, qc, identity_audits = [], [], [], [], []
    signatures, hashes = {}, {}
    try:
        for k in ("samples", "controls", "panel"):
            shutil.copyfile(config[k], out / f"{k}_input.csv")
        shutil.copyfile(config["_path"], out / "analysis_input.json")
        excluded_records = [{k: v for k, v in r.items() if not k.startswith("_")} for r in excluded]
        if excluded_records:
            pd.DataFrame(excluded_records).to_csv(out / "excluded_samples.csv", index=False)
        else:
            pd.DataFrame(columns=["sample_id", "acquisition_id", "exclusion_reason"]).to_csv(out / "excluded_samples.csv", index=False)
        write_identity_reference(out, panel)
        load_annotations().to_csv(out / "marker_annotations.csv", index=False)
        save_json(out / "marker_annotation_sources.json", load_sources())
        provenance["annotation_version"] = ANNOTATION_VERSION
        provenance["annotation_resource_hashes"] = {
            name: file_hash(out / name)
            for name in ("marker_annotations.csv", "marker_annotation_sources.json")
        }
        # Reference controls are validated first. Failed recognition blocks analysis.
        ordered = [(False, r) for r in sorted(controls, key=lambda c: c["control_type"] != "setup_beads")]
        ordered += [(True, r) for r in samples]
        for index, (is_sample, row) in enumerate(ordered):
            rid = row["acquisition_id"] if is_sample else row["control_id"]
            data = read_fcs(row["_file"])
            if data.sha256 in hashes:
                raise ConfigurationError(f"Identical FCS content assigned to {hashes[data.sha256]} and {rid}")
            hashes[data.sha256] = rid
            setup = config["setups"][row["setup_id"]]
            signature = _instrument_signature(data, setup)
            if row["setup_id"] in signatures and signatures[row["setup_id"]] != signature:
                raise ConfigurationError(f"FCS acquisition metadata differs within setup {row['setup_id']} ({rid}); split/review setups")
            signatures[row["setup_id"]] = signature
            table, quality = summarize(data, setup, panel)
            for key, value in row.items():
                if not key.startswith("_"):
                    table[key] = value
            record = {"record_id": rid, "file_type": "sample" if is_sample else row["control_type"], **quality}
            record["diagnostic_file"] = f"diagnostics/{index + 1:03d}.png"
            qc.append(record)
            provenance["inputs"].append({"id": rid, "type": record["file_type"],
                                         "path": row["_file"], "sha256": data.sha256})
            diagnostic_plot(data, setup, panel, out / "diagnostics" / f"{index + 1:03d}.png", reviewed=True)
            if not is_sample:
                control_tables.append(table)
                if row["control_type"] == "setup_beads":
                    audit = audit_identity(data, setup, panel)
                    audit["record_id"] = rid
                    audit["setup_id"] = row["setup_id"]
                    identity_audits.append(audit)
                    pd.concat(identity_audits, ignore_index=True).to_csv(out / "bead_identity_audit.csv", index=False)
                    if not audit.layout_status.eq("consistent").all():
                        raise ConfigurationError(
                            f"Setup reference {rid} disagrees with Figure 2 FITC/PE ordering or has missing populations; "
                            "check axis mapping, gate IDs, compensation and bead_identity_audit.csv. "
                            "Do not relabel ambiguous clusters by guesswork.")
                if row["control_type"] == "setup_beads" and not table.measurement_valid.all():
                    raise ConfigurationError(f"Setup reference {rid} failed population QC; review gates and acquisition")
                if row["control_type"] == background_control_type and rid in background_ids:
                    blank_tables.append(table)
            else:
                tables.append(table)
        measured = quantify(pd.concat(tables, ignore_index=True), pd.concat(blank_tables, ignore_index=True),
                            normalizers=tuple(markers), min_normalization_factor=float(floor or 0), normalize=enabled,
                            background_mode=background_mode)
        measured = annotate_measurements(measured)
        measured.to_csv(out / "measurements_long.csv", index=False)
        pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")} for r in samples]).to_csv(out / "sample_metadata.csv", index=False)
        annotate_measurements(pd.concat(control_tables, ignore_index=True)).to_csv(out / "control_measurements.csv", index=False)
        measured.loc[measured.role == "isotype_control"].to_csv(out / "isotype_controls.csv", index=False)
        for col, filename in (("raw_median_apc", "raw_median_apc.csv"), ("corrected_apc", "background_corrected.csv"),
                              ("normalized_apc", "normalized.csv")):
            measured[measured.role == "capture_marker"].pivot(index="acquisition_id", columns="marker", values=col).to_csv(out / filename)
        views = out / "marker_views"
        views.mkdir()
        for group, rows in measured.loc[measured.role == "capture_marker"].groupby("biological_group", sort=False):
            for col, scale in (("corrected_apc", "background_corrected"), ("normalized_apc", "normalized")):
                rows.pivot(index="acquisition_id", columns="marker", values=col).to_csv(views / f"{group}_{scale}.csv")
        measured[["acquisition_id", "sample_id", "subject_id", "bead_id", "n_events", "measurement_valid",
                  "correction_status", "detection_state", "normalization_status", "qc_flags"]].to_csv(out / "qc.csv", index=False)
        failed = ~measured.measurement_valid | measured.corrected_apc.isna()
        if enabled:
            failed |= measured.normalized_apc.isna()
        provenance["status"] = "completed_with_qc_failures" if failed.any() else "completed"
        provenance["failed_measurements"] = int(failed.sum())
        save_json(out / "qc.json", qc)
        save_json(out / "provenance.json", provenance)
        from .report import render_report
        render_report(out)
    except Exception as exc:
        save_json(out / "failure.json", {"status": "incomplete", "error": str(exc)})
        save_json(out / "qc.json", qc)
        provenance["status"] = "incomplete"
        save_json(out / "provenance.json", provenance)
        raise
    return out
