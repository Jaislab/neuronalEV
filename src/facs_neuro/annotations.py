"""Curated display annotations for the 39 MACSPlex EV Neuro bead identities.

The editorial groups and overlapping tags are interpretation aids, not validated
cell assignments. They never change quantitative values, create a biological
score, or estimate cell-origin fractions. The manufacturer's table establishes
assay identities only. See the packaged source registry for evidence scopes.
"""

from __future__ import annotations

import csv
from importlib.resources import files
import io
import json
import re

import pandas as pd

from .config import EXPECTED_MARKERS


ANNOTATION_VERSION = "1.0.0"
ANNOTATION_COLUMNS = (
    "biological_group", "association_tags", "gene_or_epitope",
    "interpretation_note", "annotation_source_ids", "annotation_version",
)
GROUP_LABELS = {
    "neuron_associated": "Neuron-associated",
    "glial_progenitor_associated": "Glia / progenitor-associated",
    "shared_neural_progenitor": "Shared neural / progenitor",
    "immune_myeloid": "Immune / myeloid",
    "vascular_stromal_adhesion": "Vascular / stromal / adhesion",
    "broad_receptor_signaling": "Broad receptors / signaling",
    "lysosomal_trafficking": "Lysosomal / trafficking",
    "tetraspanin_reference": "Tetraspanin references",
    "isotype_control": "Isotype controls",
}
_GROUP_IDS = {
    "neuron_associated": {"32", "33", "75", "96"},
    "glial_progenitor_associated": {"24", "52", "53", "63", "65"},
    "shared_neural_progenitor": {"46", "98"},
    "immune_myeloid": {"34", "43", "45", "55", "77", "87", "89"},
    "vascular_stromal_adhesion": {"23", "35", "54", "57", "64", "76", "78", "79", "85", "86"},
    "broad_receptor_signaling": {"22", "42", "44", "56", "74"},
    "lysosomal_trafficking": {"88"},
    "tetraspanin_reference": {"66", "67", "68"},
    "isotype_control": {"97", "99"},
}
_GROUP_FOR_ID = {bead: group for group, beads in _GROUP_IDS.items() for bead in beads}
_IDENTITY_COLUMNS = ("bead_id", "marker", "role")
_TAG_PATTERN = re.compile(r"[a-z][a-z0-9_]*\Z")


def _read_resource(name: str) -> str:
    return files("facs_neuro").joinpath("data", name).read_text(encoding="utf-8")


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Annotation source registry has duplicate key {key!r}")
        result[key] = value
    return result


def load_sources() -> dict:
    """Return the packaged registry, including explicit evidence-scope caveats."""
    try:
        registry = json.loads(
            _read_resource("annotation_sources.json"), object_pairs_hook=_unique_json_object
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read annotation source registry: {exc}") from exc
    if not isinstance(registry, dict) or registry.get("annotation_version") != ANNOTATION_VERSION:
        raise ValueError("Annotation source registry version mismatch")
    if not isinstance(registry.get("curation_notice"), str) or not registry["curation_notice"].strip():
        raise ValueError("Annotation source registry requires a curation notice")
    sources = registry.get("sources")
    if not isinstance(sources, dict) or not sources:
        raise ValueError("Annotation source registry requires sources")
    for source_id, source in sources.items():
        if not isinstance(source_id, str) or not source_id or ";" in source_id:
            raise ValueError("Annotation source IDs must be nonempty and cannot contain semicolons")
        if not isinstance(source, dict) or any(
            not isinstance(source.get(key), str) or not source[key].strip()
            for key in ("title", "evidence_scope", "verification_status")
        ):
            raise ValueError(f"Annotation source {source_id!r} lacks title, evidence scope or verification status")
        if not isinstance(source.get("url", ""), str):
            raise ValueError(f"Annotation source {source_id!r} URL must be text")
    return registry


def _parts(value: str, label: str) -> list[str]:
    parts = value.split(";")
    if any(not item or item != item.strip() for item in parts) or len(set(parts)) != len(parts):
        raise ValueError(f"{label} must contain unique nonempty semicolon-separated values")
    return parts


def load_annotations() -> pd.DataFrame:
    """Load and validate all 39 canonical bead identities and display annotations.

    The registry and CSV are package resources. A missing, duplicate, mislabeled,
    or version-inconsistent annotation fails validation instead of leaving a
    silently unannotated biological measurement.
    """
    try:
        reader = csv.DictReader(io.StringIO(_read_resource("marker_annotations.csv")))
        headers = reader.fieldnames or []
        expected_columns = {*_IDENTITY_COLUMNS, *ANNOTATION_COLUMNS}
        if len(headers) != len(set(headers)) or set(headers) != expected_columns:
            raise ValueError("Annotation CSV must contain the exact identity and annotation columns")
        rows = list(reader)
    except (OSError, csv.Error) as exc:
        raise ValueError(f"Cannot read marker annotations: {exc}") from exc
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("Annotation CSV contains malformed rows")
    table = pd.DataFrame(rows, columns=headers)
    if len(table) != 39 or table["bead_id"].duplicated().any() or set(table["bead_id"]) != set(EXPECTED_MARKERS):
        raise ValueError("Annotations must cover all 39 canonical bead IDs exactly once")
    sources = load_sources()["sources"]
    for row in rows:
        bead = row["bead_id"]
        expected_role = "isotype_control" if bead in {"97", "99"} else "capture_marker"
        if row["marker"] != EXPECTED_MARKERS[bead] or row["role"] != expected_role:
            raise ValueError(f"Annotation canonical identity mismatch for bead {bead}")
        if row["biological_group"] != _GROUP_FOR_ID[bead]:
            raise ValueError(f"Annotation display group mismatch for bead {bead}")
        if row["annotation_version"] != ANNOTATION_VERSION:
            raise ValueError(f"Annotation version mismatch for bead {bead}")
        for column in ("gene_or_epitope", "interpretation_note"):
            if not row[column].strip():
                raise ValueError(f"Annotation {column} missing for bead {bead}")
        tags = _parts(row["association_tags"], "Association tags")
        if any(not _TAG_PATTERN.fullmatch(tag) for tag in tags):
            raise ValueError(f"Invalid association tag for bead {bead}")
        if any(tag.endswith("_specific") or tag.startswith("cell_origin") for tag in tags):
            raise ValueError(f"Unsupported specificity or cell-origin claim for bead {bead}")
        refs = _parts(row["annotation_source_ids"], "Annotation source IDs")
        if "kit_table1" not in refs or "curated_interpretation_v1" not in refs:
            raise ValueError(f"Annotation identity/curation provenance missing for bead {bead}")
        if any(ref not in sources for ref in refs):
            raise ValueError(f"Unknown annotation source for bead {bead}")
        if expected_role == "isotype_control" and tags != ["assay_control", "isotype_background"]:
            raise ValueError("Isotype controls cannot receive biological association tags")
    return table


def annotate_measurements(frame: pd.DataFrame) -> pd.DataFrame:
    """Append display annotations without changing existing measurement values.

    Canonical identity is checked on every row, including controls. Original
    columns, values, dtypes, row order and index (including duplicate index
    labels) are preserved. Matching preexisting annotations are accepted for
    idempotent reporting; conflicting preexisting values are rejected.
    """
    if not isinstance(frame, pd.DataFrame):
        raise ValueError("Measurements must be a pandas DataFrame")
    if frame.columns.duplicated().any() or any(column not in frame for column in _IDENTITY_COLUMNS):
        raise ValueError("Measurements require unique columns including bead_id, marker and role")
    annotations = load_annotations().set_index("bead_id", verify_integrity=True)
    keys = frame["bead_id"].map(str).tolist()
    for position, (bead, marker, role) in enumerate(zip(keys, frame["marker"], frame["role"], strict=True)):
        if bead not in annotations.index:
            raise ValueError(f"Unknown canonical bead ID at measurement row {position}: {bead!r}")
        expected = annotations.loc[bead]
        if not isinstance(marker, str) or marker != expected["marker"] or not isinstance(role, str) or role != expected["role"]:
            raise ValueError(f"Measurement canonical identity mismatch at row {position} (bead {bead})")
    result = frame.copy(deep=True)
    for column in ANNOTATION_COLUMNS:
        values = [annotations.at[bead, column] for bead in keys]
        if column in result:
            existing = result[column].tolist()
            if any(not isinstance(actual, str) or actual != expected for actual, expected in zip(existing, values, strict=True)):
                raise ValueError(f"Conflicting preexisting annotation column: {column}")
        else:
            # Assign a list, not an indexed Series: duplicate row labels are safe.
            result[column] = values
    return result
