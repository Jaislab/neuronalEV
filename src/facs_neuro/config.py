"""Configuration and manifest validation; all paths have explicit anchors."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pandas as pd


class ConfigurationError(ValueError):
    """Actionable invalid or incomplete input configuration."""


EXPECTED_IDS = {
    "22", "23", "24", "32", "33", "34", "35", "42", "43", "44", "45", "46",
    "52", "53", "54", "55", "56", "57", "63", "64", "65", "66", "67", "68",
    "74", "75", "76", "77", "78", "79", "85", "86", "87", "88", "89", "96",
    "97", "98", "99",
}
EXPECTED_MARKERS = dict(zip(
    "22 23 24 32 33 34 35 42 43 44 45 46 52 53 54 55 56 57 63 64 65 66 67 68 74 75 76 77 78 79 85 86 87 88 89 96 97 98 99".split(),
    "ErbB-2,CD49a,CD140a,CD171,CD56,CD13,CD31,CD222,CD11b,CD24,CD45,CD133/1,O4,GLAST,CD90,CX3CR1,EGFR,Podoplanin,A2B5,CD29,CSPG4,CD9,CD63,CD81,CD47,PSA-NCAM,CD44,CD38,CD49f,CD106,CD49e,CD54,CD36,CD107a,CD68,VGLUT2,REA control,Ganglioside GD2,mIgG1 control".split(","),
    strict=True,
))
RESERVED_COLUMNS = {
    "bead_id", "marker", "role", "n_events", "raw_median_apc", "measurement_valid", "qc_flags",
    "blank_median_apc", "corrected_apc", "detection_state", "correction_status", "normalization_factor",
    "normalized_apc", "normalization_status",
    "biological_group", "association_tags", "gene_or_epitope", "interpretation_note",
    "annotation_source_ids", "annotation_version",
}
MATCH_FIELDS = (
    "batch_id", "setup_id", "instrument_settings_id", "kit_lot",
    "detection_reagent", "incubation_protocol",
)


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_config(path: str | Path) -> dict:
    path = Path(path).resolve()
    try:
        config = json.loads(path.read_text(), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"Cannot read configuration {path}: {exc}") from exc
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise ConfigurationError("Configuration must be an object with schema_version=1")
    if config.get("dataset_kind") not in {"experimental", "synthetic"}:
        raise ConfigurationError("dataset_kind must be experimental or synthetic")
    for key in ("samples", "controls", "panel"):
        if not isinstance(config.get(key), str) or not config[key].strip():
            raise ConfigurationError(f"Missing {key} manifest path")
        config[key] = str((path.parent / config[key]).resolve())
    if not isinstance(config.get("setups"), dict) or not config["setups"]:
        raise ConfigurationError("Define at least one acquisition setup")
    if not isinstance(config.get("normalization", {}), dict):
        raise ConfigurationError("normalization must be an object")
    for key, setup in config["setups"].items():
        if not isinstance(key, str) or not key.strip() or not isinstance(setup, dict):
            raise ConfigurationError(f"Setup {key} must be an object")
        channels = setup.get("channels", {})
        if not isinstance(channels, dict) or any(not isinstance(channels.get(c), str) or not channels[c].strip() for c in ("fitc", "pe", "apc")):
            raise ConfigurationError(f"Setup {key} needs FITC, PE and APC channel mappings")
        if len({channels[c] for c in ("fitc", "pe", "apc")}) != 3:
            raise ConfigurationError(f"Setup {key} requires three distinct fluorescence channels")
    config["_path"] = str(path)
    config["_sha256"] = file_hash(path)
    return config


def read_manifest(path: str | Path, required: tuple[str, ...]) -> list[dict]:
    path = Path(path)
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            headers = reader.fieldnames or []
            if any(h.startswith("_") for h in headers):
                raise ConfigurationError(f"{path}: underscore-prefixed headers are reserved")
            if len(set(headers)) != len(headers) or any(c not in headers for c in required):
                raise ConfigurationError(f"{path}: unique headers required, including {required}")
            rows = []
            for line, row in enumerate(reader, 2):
                if None in row or any(v is None for v in row.values()):
                    raise ConfigurationError(f"{path}:{line}: malformed CSV row")
                row = {k: v.strip() for k, v in row.items()}
                if not any(row.values()):
                    continue
                row["_manifest"] = str(path.resolve())
                row["_line"] = line
                row["_file"] = str((path.parent / row["fcs_path"]).resolve()) if row["fcs_path"] else ""
                rows.append(row)
            return rows
    except OSError as exc:
        raise ConfigurationError(f"Cannot read manifest {path}: {exc}") from exc


def load_inputs(config: dict) -> tuple[list[dict], list[dict], pd.DataFrame]:
    samples = read_manifest(config["samples"], ("sample_id", "fcs_path"))
    controls = read_manifest(config["controls"], ("control_id", "control_type", "fcs_path"))
    for label, rows, extra in (("samples", samples, {"control_id", "control_type"}),
                               ("controls", controls, {"sample_id", "acquisition_id", "subject_id"})):
        for row in rows:
            reserved = set(row) & (RESERVED_COLUMNS | extra)
            private = {key for key in row if key.startswith("_") and key not in {"_file", "_manifest", "_line"}}
            if reserved or private:
                raise ConfigurationError(f"{label} manifest contains reserved columns: {sorted(reserved | private)}")
    if not samples:
        raise ConfigurationError("Sample manifest is empty")
    seen = set()
    for sample in samples:
        if not sample["sample_id"]:
            raise ConfigurationError("Every sample needs a nonempty sample_id")
        sample["acquisition_id"] = sample.get("acquisition_id") or sample["sample_id"]
        if sample["acquisition_id"] in seen:
            raise ConfigurationError("Duplicate acquisition_id; give technical acquisitions unique IDs")
        seen.add(sample["acquisition_id"])
    cids = [r["control_id"] for r in controls]
    if any(not x for x in cids) or len(cids) != len(set(cids)):
        raise ConfigurationError("Every control needs a unique nonempty control_id")
    try:
        panel = pd.read_csv(config["panel"], dtype=str, keep_default_na=False)
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"Cannot read panel: {exc}") from exc
    if not {"bead_id", "marker", "role"} <= set(panel.columns):
        raise ConfigurationError("Panel requires bead_id, marker and role")
    if len(panel) != 39 or set(panel.bead_id) != EXPECTED_IDS or panel.marker.duplicated().any():
        raise ConfigurationError("Panel must identify all 39 unique Neuro bead IDs and marker names")
    for row in panel.to_dict("records"):
        expected = "isotype_control" if row["bead_id"] in {"97", "99"} else "capture_marker"
        if row["role"] != expected or not row["marker"]:
            raise ConfigurationError(f"Incorrect role or missing marker for bead {row['bead_id']}")
    actual = panel.set_index("bead_id").marker
    if any(actual[i] != name for i, name in EXPECTED_MARKERS.items()):
        raise ConfigurationError("Panel bead identity mismatch; use canonical marker names and keep aliases in a separate column")
    return samples, controls, panel


def readiness_errors(config: dict, samples: list[dict], controls: list[dict]) -> list[str]:
    errors = []
    by_id = {c["control_id"]: c for c in controls}
    files = {}
    for kind, rows in (("sample", samples), ("control", controls)):
        for row in rows:
            rid = row["acquisition_id"] if kind == "sample" else row["control_id"]
            prefix = f"{kind} {rid}"
            fpath = row["_file"]
            if not fpath or not Path(fpath).is_file():
                errors.append(f"{prefix}: FCS file missing ({row['fcs_path'] or 'path not supplied'})")
            elif fpath in files:
                errors.append(f"{prefix}: same FCS path already mapped to {files[fpath]}")
            else:
                files[fpath] = prefix
            if row.get("setup_id") not in config["setups"]:
                errors.append(f"{prefix}: missing or unknown setup_id")
            keys = MATCH_FIELDS if kind == "sample" or row.get("control_type") == "assay_blank" else (
                "batch_id", "setup_id", "instrument_settings_id", "kit_lot",
            )
            for field in keys:
                if not row.get(field) or row[field].lower() in {"unknown", "pending", "not_applicable"}:
                    errors.append(f"{prefix}: {field} must be documented")
            if kind == "sample":
                if not row.get("subject_id") or row["subject_id"].lower() in {"unknown", "pending", "not_applicable"}:
                    errors.append(f"{prefix}: subject_id required to preserve donor identity")
                for field, expected_type in (("blank_id", "assay_blank"), ("setup_control_id", "setup_beads")):
                    ref = row.get(field)
                    control = by_id.get(ref)
                    if control is None or control.get("control_type") != expected_type:
                        errors.append(f"{prefix}: {field} must reference {expected_type}")
                        continue
                    match = MATCH_FIELDS if expected_type == "assay_blank" else MATCH_FIELDS[:4]
                    for key in match:
                        if row.get(key) != control.get(key):
                            errors.append(f"{prefix}: {key} differs from referenced {expected_type} {ref}")
    return errors
