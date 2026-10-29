"""Self-contained, descriptive reports for MACSPlex bead-assay runs.

No inferential statistics or assay-validation claims are generated here. Numeric
tables retain the analysis scale; asinh is used only to display signed signals.
"""

from __future__ import annotations

import base64
from html import escape
from io import BytesIO
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


_NUMBER_COLUMNS = (
    "raw_median_apc", "corrected_apc", "normalized_apc", "n_events"
)
_REQUIRED_COLUMNS = {
    "sample_id", "acquisition_id", "bead_id", "marker", "role",
    "raw_median_apc", "corrected_apc", "normalized_apc", "n_events",
}


def _natural_key(value: object) -> tuple:
    return tuple((1, int(part)) if part.isdigit() else (0, part.lower())
                 for part in re.split(r"(\d+)", str(value)))


def _figure_data_url(figure: Figure) -> str:
    buffer = BytesIO()
    FigureCanvasAgg(figure).print_png(buffer)
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def _table(frame: pd.DataFrame, *, max_rows: int | None = None) -> str:
    if frame.empty:
        return "<p>No records.</p>"
    display = frame if max_rows is None else frame.head(max_rows)
    table = display.to_html(index=False, escape=True, na_rep="unavailable", border=0,
                            float_format=lambda x: f"{x:.6g}")
    if max_rows is not None and len(frame) > max_rows:
        table += (f"<p>Showing {max_rows} of {len(frame)} rows; "
                  "download measurements_long.csv for every measurement.</p>")
    return '<div class="table-scroll">' + table + "</div>"


def _read_json(path: Path, default: object) -> object:
    if not path.exists():
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read report input {path.name}: {exc}") from exc


def _heatmap(frame: pd.DataFrame, value: str, title: str) -> str:
    capture = frame.loc[frame["role"] == "capture_marker"].copy()
    if capture.empty:
        return f"<h3>{escape(title)}</h3><p>No capture-marker measurements.</p>"
    if value == "raw_median_apc" and "measurement_valid" in capture:
        valid = capture["measurement_valid"].str.lower().isin({"true", "1"})
        capture.loc[~valid, value] = np.nan
    if capture.duplicated(["acquisition_id", "bead_id"]).any():
        raise ValueError("Duplicate acquisition/bead measurements cannot be plotted without aggregation.")
    # A single marker identity is required per bead ID; never relabel by row order.
    names = capture[["bead_id", "marker"]].drop_duplicates()
    if names["bead_id"].duplicated().any():
        raise ValueError("Inconsistent marker names for the same bead ID.")
    names = names.set_index("bead_id")["marker"].to_dict()
    matrix = capture.pivot(index="acquisition_id", columns="bead_id", values=value)
    matrix = matrix.reindex(index=sorted(matrix.index, key=_natural_key),
                            columns=sorted(matrix.columns, key=_natural_key))
    values = matrix.to_numpy(dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        detail = ("Normalization is unavailable for every capture-marker measurement. "
                  "See normalization_status and QC flags in the numerical tables."
                  if value == "normalized_apc" else
                  "No valid finite measurements are available on this scale.")
        return f"<h3>{escape(title)}</h3><p>{escape(detail)}</p>"
    displayed = np.ma.masked_invalid(np.arcsinh(values))
    figure = Figure(figsize=(max(8.0, 0.33 * len(matrix.columns) + 3),
                             max(3.0, 0.29 * len(matrix.index) + 2.5)),
                    layout="constrained")
    axis = figure.subplots()
    signed = bool((values[finite] < 0).any())
    from matplotlib import colormaps
    cmap = colormaps["coolwarm" if signed else "viridis"].copy()
    cmap.set_bad("#d1d5db")
    limits = {}
    if signed:
        limit = float(np.max(np.abs(displayed.compressed())))
        limits = {"vmin": -limit, "vmax": limit}
    image = axis.imshow(displayed, interpolation="nearest", aspect="auto", cmap=cmap, **limits)
    axis.set_xticks(range(len(matrix.columns)),
                    [f"{names[col]} [{col}]" for col in matrix.columns],
                    rotation=70, ha="right", fontsize=8, parse_math=False)
    axis.set_yticks(range(len(matrix.index)), matrix.index, fontsize=8, parse_math=False)
    axis.set_ylabel("Acquisition ID (not independent bead events)")
    axis.set_title(title)
    colorbar = figure.colorbar(image, ax=axis, shrink=0.75)
    colorbar.set_label("asinh(signal); display transform only")
    caption = ("Gray cells are unavailable or QC-invalid; they are not biological zeros. "
               "Colors use asinh(value) without imputation or per-marker scaling. "
               "Numerical CSV files retain signed values on the analysis scale.")
    return (f"<h3>{escape(title)}</h3><img class=\"heatmap\" alt=\"{escape(title, quote=True)}\" "
            f"src=\"{_figure_data_url(figure)}\"><p class=\"caption\">{caption}</p>")


def render_report(run_dir: Path | str) -> Path:
    """Render a run's CSV/JSON results to a standalone ``report.html``.

    Identifiers and metadata are escaped. Diagnostic PNGs are embedded rather
    than referenced externally. Existing report.html is refreshed in place;
    input measurement/provenance files are never modified.
    """
    root = Path(run_dir)
    csv_path = root / "measurements_long.csv"
    try:
        frame = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ValueError(f"Cannot read measurements_long.csv: {exc}") from exc
    missing = sorted(_REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"Report measurements missing columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError("No measurements are available to report.")
    for column in _NUMBER_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    provenance = _read_json(root / "provenance.json", {})
    if not isinstance(provenance, dict):
        raise ValueError("provenance.json must contain an object.")
    if provenance.get("status") == "incomplete":
        raise ValueError("Cannot report an incomplete run")
    qc = _read_json(root / "qc.json", [])
    if not isinstance(qc, list) or any(not isinstance(row, dict) for row in qc):
        raise ValueError("qc.json must contain a list of objects.")
    synthetic = provenance.get("dataset_kind") == "synthetic"
    banner = ("SYNTHETIC DATA — software demonstration only; no donor measurements."
              if synthetic else
              "EXPERIMENTAL DATA — assay performance and gates require reviewed real-file validation.")
    validation = provenance.get("validation_status", "not provided")
    subjects = (frame["subject_id"].replace("", pd.NA).dropna().nunique()
                if "subject_id" in frame else None)
    n_markers = frame.loc[frame["role"] == "capture_marker", "bead_id"].nunique()
    n_controls = frame.loc[frame["role"] == "isotype_control", "bead_id"].nunique()
    pieces = ["<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
              "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
              "<title>MACSPlex EV Neuro analysis report</title>",
              "<style>body{font:16px/1.5 system-ui,sans-serif;color:#172033;max-width:1400px;"
              "margin:2rem auto;padding:0 1.2rem}h1,h2,h3{line-height:1.2}"
              ".banner{padding:1rem;border:2px solid #b45309;background:#fffbeb;font-weight:700}"
              ".summary{padding:1rem;background:#eff6ff}table{border-collapse:collapse;font-size:.86rem}"
              "th,td{padding:.4rem .65rem;border:1px solid #d1d5db;text-align:left}"
              "th{background:#f3f4f6}.table-scroll{overflow:auto;max-height:34rem}"
              "img{max-width:100%;height:auto}.caption{color:#4b5563;font-size:.9rem}"
              "pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f3f4f6;padding:1rem}"
              "details{margin:1rem 0}a{color:#1d4ed8}</style></head><body>",
              "<h1>MACSPlex EV Neuro analysis report</h1>",
              f"<p class=\"banner\">{escape(banner)}</p>",
              f"<p>Validation status recorded by the run: <strong>{escape(str(validation))}</strong>.</p>",
              "<div class=\"summary\"><p>"
              f"{frame['sample_id'].nunique()} sample IDs; {frame['acquisition_id'].nunique()} acquisitions; "
              f"{subjects if subjects is not None else 'unknown'} donor/subject IDs; "
              f"{n_markers} capture-marker and {n_controls} isotype-control populations.</p>"
              "<p>Samples are planned as independent donors. Technical acquisitions are displayed separately; "
              "bead-event counts measure assay sampling and do not increase biological sample size.</p></div>",
              "<h2>Interpretation</h2><p>FITC/PE identify bead populations; APC is the detection endpoint. "
              "The CD9/CD63/CD81 cocktail produces a conditional capture/detection signal. "
              "These measurements do not establish absolute EV counts, EV purity, neuronal origin, "
              "or single-EV co-expression. ExoQuick enrichment can retain non-EV plasma material.</p>",
              "<p>Signed background-corrected values are retained. Negative corrected values are labeled "
              "not detected; zero is at background. A positive corrected value alone is not a validated "
              "detection call. No validated LOD/LOQ or inferential statistics are supplied by this report. "
              "Missing/invalid measurements are unavailable, not zero. Tetraspanin normalization is optional "
              "and changes the biological question.</p>",
              "<h2>Marker profiles</h2>"]
    pieces.extend(_heatmap(frame, key, title) for key, title in (
        ("raw_median_apc", "Raw APC medians — QC-valid measurements"),
        ("corrected_apc", "Background-corrected APC"),
        ("normalized_apc", "Optional tetraspanin-normalized signal"),
    ))
    pieces.append("<h2>Quality control</h2>")
    pieces.append(_table(pd.DataFrame(qc)) if qc else "<p>No acquisition-level QC records were supplied.</p>")
    flags = frame.get("qc_flags", pd.Series("", index=frame.index)).fillna("").astype(str)
    invalid = (frame["measurement_valid"].str.lower().isin({"false", "0"})
               if "measurement_valid" in frame else pd.Series(False, index=frame.index))
    issues = frame.loc[(flags.str.strip() != "") | invalid]
    issue_columns = [c for c in ("sample_id", "acquisition_id", "bead_id", "marker", "n_events",
                                "measurement_valid", "qc_flags", "correction_status",
                                "normalization_status") if c in frame]
    pieces.append(f"<p>{len(issues)} measurement rows have QC flags or invalid measurements.</p>")
    pieces.append(_table(issues[issue_columns]))
    pieces.append("<h2>Isotype controls</h2><p>Isotype controls diagnose nonspecific binding. "
                  "They do not replace the marker-specific assay blank and are not automatically "
                  "subtracted a second time.</p>")
    control_columns = [c for c in ("sample_id", "acquisition_id", "bead_id", "marker", "n_events",
                                  "raw_median_apc", "corrected_apc", "qc_flags") if c in frame]
    pieces.append(_table(frame.loc[frame["role"] == "isotype_control", control_columns]))
    pieces.append("<h2>Acquisition and gating diagnostics</h2>")
    diagnostics = sorted((root / "diagnostics").glob("*.png"), key=lambda p: _natural_key(p.name))
    if not diagnostics:
        pieces.append("<p>No diagnostic PNGs were supplied.</p>")
    for path in diagnostics:
        # Fixed data:image MIME prevents the filename or metadata from introducing markup.
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        reference = next((str(row.get("record_id", path.stem)) for row in qc
                          if row.get("diagnostic_file") == f"diagnostics/{path.name}"), path.stem)
        pieces.append(f"<details><summary>{escape(reference)} ({escape(path.stem)})</summary>"
                      f"<img alt=\"{escape(path.stem, quote=True)}\" "
                      f"src=\"data:image/png;base64,{encoded}\"></details>")
    pieces.append("<h2>Numerical outputs</h2><p>Signal values in CSV files remain on the analysis scale. "
                  "These download links require the run folder alongside the standalone report.</p><ul>")
    for name in ("measurements_long.csv", "raw_median_apc.csv", "background_corrected.csv",
                 "normalized.csv", "isotype_controls.csv", "qc.csv", "sample_metadata.csv"):
        if (root / name).is_file():
            pieces.append(f"<li><a href=\"{escape(name, quote=True)}\">{escape(name)}</a></li>")
    pieces.append("</ul><h2>Run provenance</h2><pre>")
    pieces.append(escape(json.dumps(provenance, indent=2, ensure_ascii=False)))
    pieces.append("</pre></body></html>")
    output = root / "report.html"
    output.write_text("\n".join(pieces), encoding="utf-8")
    return output
