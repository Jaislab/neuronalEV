"""Auditable background correction and optional tetraspanin normalization.

Inputs contain one *acquisition* per acquisition/control identifier and one row per
bead population. Controls are matched by explicit ID; no controls, events, or
biological replicates are implicitly pooled. This module does not establish a
detection limit: a positive assay-blank-corrected value is labelled
``unclassified_positive``. Explicit setup-baseline subtraction is a limited
alternative: it does not measure detection-reagent background and never makes
detection/nondetection claims.
"""

from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np
import pandas as pd


_MEASUREMENT_FIELDS = {
    "batch_id", "setup_id", "bead_id", "marker", "role", "n_events",
    "raw_median_apc", "measurement_valid",
}
_SAMPLE_FIELDS = _MEASUREMENT_FIELDS | {"sample_id", "subject_id", "blank_id"}
_BLANK_FIELDS = _MEASUREMENT_FIELDS | {"control_id"}
SETUP_BASELINE_LIMITATION = "setup_baseline_only_reagent_background_unmeasured"


def _identifier(value: object) -> str:
    """Use exact textual identifiers, with absent values represented as empty."""
    if value is None or pd.isna(value):
        return ""
    return str(value)


def _require_columns(frame: pd.DataFrame, fields: set[str], label: str) -> None:
    missing = sorted(fields.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
    if frame.columns.duplicated().any():
        raise ValueError(f"{label} has duplicate column names")


def _validate_frame(frame: pd.DataFrame, id_column: str, label: str) -> None:
    for column in (id_column, "bead_id", "batch_id", "setup_id", "marker", "role"):
        if frame[column].map(_identifier).str.strip().eq("").any():
            raise ValueError(f"{label}.{column} must contain nonempty identifiers")
    keys = pd.DataFrame({
        id_column: frame[id_column].map(_identifier),
        "bead_id": frame["bead_id"].map(_identifier),
    })
    if keys.duplicated().any():
        raise ValueError(f"{label} has duplicate {id_column}/bead_id rows")
    if not frame["measurement_valid"].map(lambda x: isinstance(x, (bool, np.bool_))).all():
        raise ValueError(f"{label}.measurement_valid must contain booleans, not strings")
    consistent_fields = ["batch_id", "setup_id"]
    if id_column in {"sample_id", "acquisition_id"}:
        consistent_fields += ["subject_id", "blank_id"]
        if id_column == "acquisition_id":
            consistent_fields += ["sample_id"]
    for identifier, group in frame.groupby(frame[id_column].map(_identifier), sort=False):
        for column in consistent_fields:
            if group[column].map(_identifier).nunique() > 1:
                raise ValueError(f"{label} {identifier!r} has inconsistent {column}")


def _valid_measurement(row: pd.Series) -> bool:
    try:
        count = float(row["n_events"])
        raw = float(row["raw_median_apc"])
    except (TypeError, ValueError):
        return False
    return bool(row["measurement_valid"]) and math.isfinite(raw) and (
        math.isfinite(count) and count > 0 and count.is_integer()
    )


def _join_flags(existing: object, added: Iterable[str]) -> str:
    flags = [] if existing is None or pd.isna(existing) else str(existing).split(";")
    return ";".join(dict.fromkeys(flag for flag in [*flags, *added] if flag))


def quantify(
    samples: pd.DataFrame,
    blanks: pd.DataFrame,
    *,
    normalizers: tuple[str, ...] = ("CD9", "CD63", "CD81"),
    min_normalization_factor: float = 0.0,
    normalize: bool = False,
    background_mode: str = "assay_blank",
) -> pd.DataFrame:
    """Subtract exact, compatible per-bead reference medians.

    ``measurement_valid`` is an upstream gate/QC decision, not a truthy string.
    Counts must also be positive integers and medians finite. No count threshold
    is invented here. Invalid measurements and absent, mismatched, or invalid
    blanks yield NaN corrected values and explicit QC reasons. Existing columns,
    row order, and input indexes are retained; input frames are never modified.

    If present, ``acquisition_id`` distinguishes technical acquisitions of the
    same ``sample_id``. Otherwise ``sample_id`` is the acquisition identifier.

    When requested, normalization uses the arithmetic mean of **all** declared
    markers for that acquisition. Every normalizer must be a unique valid
    capture-marker measurement with positive corrected signal, and their factor
    must exceed the configured nonnegative floor. No per-sample marker fallback
    is used. Isotypes retain their role and may be corrected/normalized for
    diagnostics, but cannot be selected as normalizers.

    The default ``assay_blank`` mode retains the original full-reagent-blank
    semantics. Explicit ``setup_baseline`` mode requires setup-bead references,
    records its limitation on every row (including normalized measurements),
    and labels values only as below/at/above the setup baseline. It does not
    establish positivity, LOD, or removal of nonspecific detection-reagent signal.
    """
    _require_columns(samples, _SAMPLE_FIELDS, "samples")
    acquisition_column = "acquisition_id" if "acquisition_id" in samples else "sample_id"
    _validate_frame(samples, acquisition_column, "samples")
    for column in ("sample_id", "subject_id"):
        if samples[column].map(_identifier).str.strip().eq("").any():
            raise ValueError(f"samples.{column} must contain nonempty identifiers")
    for sample_id, group in samples.groupby(samples["sample_id"].map(_identifier), sort=False):
        if group["subject_id"].map(_identifier).nunique() > 1:
            raise ValueError(f"sample {sample_id!r} maps to multiple subject_id values")
    if not isinstance(background_mode, str) or background_mode not in {"assay_blank", "setup_baseline"}:
        raise ValueError("background_mode must be assay_blank or setup_baseline")
    if not isinstance(normalize, (bool, np.bool_)):
        raise ValueError("normalize must be a boolean")
    if isinstance(min_normalization_factor, (bool, np.bool_)):
        raise ValueError("min_normalization_factor must be finite and nonnegative, not boolean")
    try:
        floor = float(min_normalization_factor)
    except (TypeError, ValueError) as exc:
        raise ValueError("min_normalization_factor must be finite and nonnegative") from exc
    if not math.isfinite(floor) or floor < 0:
        raise ValueError("min_normalization_factor must be finite and nonnegative")
    if normalize and (
        isinstance(normalizers, str) or not normalizers
        or any(not isinstance(marker, str) or not marker.strip() for marker in normalizers)
        or len(set(normalizers)) != len(normalizers)
    ):
        raise ValueError("normalizers must be a nonempty sequence of unique marker names")

    # A truly empty frame is accepted for a workflow with no blank files yet.
    blank_lookup: dict[tuple[str, str], pd.Series] = {}
    if len(blanks):
        _require_columns(blanks, _BLANK_FIELDS, "blanks")
        if background_mode == "setup_baseline" and "control_type" not in blanks:
            raise ValueError("setup_baseline reference rows require control_type='setup_beads'")
        _validate_frame(blanks, "control_id", "blanks")
        blank_lookup = {
            (_identifier(row["control_id"]), _identifier(row["bead_id"])): row
            for _, row in blanks.iterrows()
        }

    result = samples.copy(deep=True)
    if "qc_flags" not in result:
        result["qc_flags"] = ""
    result["background_mode"] = background_mode
    result["background_reference_id"] = result["blank_id"].copy()
    result["background_limitation"] = (
        SETUP_BASELINE_LIMITATION if background_mode == "setup_baseline" else ""
    )
    if background_mode == "setup_baseline":
        result["qc_flags"] = result["qc_flags"].map(
            lambda value: _join_flags(value, ["setup_baseline_only"])
        )
    result["blank_median_apc"] = np.nan
    result["corrected_apc"] = np.nan
    result["correction_status"] = "unavailable"
    result["detection_state"] = "unavailable"
    result["normalization_factor"] = np.nan
    result["normalized_apc"] = np.nan
    result["normalization_status"] = "not_requested"

    # Positional updates preserve arbitrary, including non-unique, input indexes.
    columns = {name: result.columns.get_loc(name) for name in result.columns}

    def put(position: int, column: str, value: object) -> None:
        result.iat[position, columns[column]] = value

    for position, (_, sample) in enumerate(samples.iterrows()):
        reasons = []
        if not _valid_measurement(sample):
            reasons.append("invalid_sample_measurement")
        blank = blank_lookup.get((
            _identifier(sample["blank_id"]), _identifier(sample["bead_id"])
        ))
        if blank is None:
            reasons.append(
                "missing_setup_baseline" if background_mode == "setup_baseline" else "missing_assay_blank"
            )
        else:
            blank_reasons = []
            if "control_type" in blank.index:
                expected_type = "setup_beads" if background_mode == "setup_baseline" else "assay_blank"
                if _identifier(blank["control_type"]) != expected_type:
                    blank_reasons.append("background_reference_type_mismatch")
            for column in ("batch_id", "setup_id", "marker", "role"):
                if _identifier(sample[column]) != _identifier(blank[column]):
                    blank_reasons.append(f"blank_{column.removesuffix('_id')}_mismatch")
            if not _valid_measurement(blank):
                blank_reasons.append("invalid_blank_measurement")
            if not blank_reasons:
                put(position, "blank_median_apc", float(blank["raw_median_apc"]))
            reasons.extend(blank_reasons)
        if reasons:
            put(position, "correction_status", ";".join(reasons))
            put(position, "qc_flags", _join_flags(result.iloc[position]["qc_flags"], reasons))
            continue
        corrected = float(sample["raw_median_apc"]) - float(blank["raw_median_apc"])
        if not math.isfinite(corrected):
            reason = "nonfinite_corrected_value"
            put(position, "correction_status", reason)
            put(position, "qc_flags", _join_flags(result.iloc[position]["qc_flags"], [reason]))
            continue
        put(position, "corrected_apc", corrected)
        put(position, "correction_status", "valid")
        if background_mode == "setup_baseline":
            state = (
                "below_setup_baseline" if corrected < 0 else
                "at_setup_baseline" if corrected == 0 else "above_setup_baseline"
            )
        else:
            state = (
                "not_detected" if corrected < 0 else
                "at_background" if corrected == 0 else "unclassified_positive"
            )
        put(position, "detection_state", state)

    if not normalize:
        return result

    acquisition_ids = result[acquisition_column].map(_identifier).to_numpy()
    for acquisition_id in dict.fromkeys(acquisition_ids):
        positions = np.flatnonzero(acquisition_ids == acquisition_id)
        group = result.iloc[positions]
        selected = []
        status = "valid"
        for marker in normalizers:
            candidates = group.loc[group["marker"] == marker]
            if candidates.empty:
                status = "missing_normalizer"
                break
            if len(candidates) != 1:
                status = "ambiguous_normalizer"
                break
            row = candidates.iloc[0]
            if row["role"] != "capture_marker" or row["correction_status"] != "valid":
                status = "invalid_normalizer"
                break
            value = float(row["corrected_apc"])
            if not math.isfinite(value) or value <= 0:
                status = "nonpositive_normalizer"
                break
            selected.append(value)
        factor = math.fsum(value / len(selected) for value in selected) if status == "valid" else np.nan
        if status == "valid" and (not math.isfinite(factor) or factor <= floor):
            status = "factor_below_minimum"
            factor = np.nan
        for position in positions:
            put(position, "normalization_status", status)
            if status != "valid":
                put(position, "qc_flags", _join_flags(
                    result.iloc[position]["qc_flags"], [f"normalization_{status}"]
                ))
                continue
            put(position, "normalization_factor", factor)
            if result.iloc[position]["correction_status"] != "valid":
                put(position, "normalization_status", "invalid_measurement")
                continue
            normalized = float(result.iloc[position]["corrected_apc"]) / factor
            if math.isfinite(normalized):
                put(position, "normalized_apc", normalized)
            else:
                put(position, "normalization_status", "nonfinite_normalized_value")
                put(position, "qc_flags", _join_flags(
                    result.iloc[position]["qc_flags"], ["nonfinite_normalized_value"]
                ))
    return result
