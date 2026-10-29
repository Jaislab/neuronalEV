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
from urllib.parse import urlsplit

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
_BIOLOGICAL_FOCUS = (
    "neuron_associated", "glial_progenitor_associated", "shared_neural_progenitor",
)
_BIOLOGICAL_GROUP_ORDER = _BIOLOGICAL_FOCUS + (
    "immune_myeloid", "vascular_stromal_adhesion", "broad_receptor_signaling",
    "lysosomal_trafficking", "tetraspanin_reference",
)
_ANNOTATION_COLUMNS = (
    "biological_group", "association_tags", "gene_or_epitope", "interpretation_note",
    "annotation_source_ids", "annotation_version",
)


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


def _background_context(frame: pd.DataFrame, provenance: dict) -> tuple[str, str]:
    """Use explicit run semantics and reject conflicting measurement labels."""
    background = provenance.get("background", {})
    if not isinstance(background, dict):
        raise ValueError("Run background provenance must be an object")
    modes = set(frame["background_mode"]) if "background_mode" in frame else set()
    if len(modes) > 1 or modes - {"assay_blank", "setup_baseline"}:
        raise ValueError("Report measurements contain inconsistent background modes")
    mode = background.get("mode", next(iter(modes), "assay_blank"))
    if mode not in {"assay_blank", "setup_baseline"} or (modes and modes != {mode}):
        raise ValueError("Measurements and run provenance disagree on background mode")
    return mode, str(background.get("rationale", ""))


def _exclusion_section(root: Path, frame: pd.DataFrame, provenance: dict) -> str:
    path = root / "excluded_samples.csv"
    if path.is_file():
        try:
            excluded = pd.read_csv(path, dtype=str, keep_default_na=False)
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
            raise ValueError(f"Cannot read excluded_samples.csv: {exc}") from exc
    else:
        records = provenance.get("sample_exclusions", [])
        if not isinstance(records, list) or any(not isinstance(row, dict) for row in records):
            raise ValueError("sample_exclusions provenance must be a list of records")
        excluded = pd.DataFrame(records)
    if excluded.empty:
        return ""
    if not {"sample_id", "exclusion_reason"}.issubset(excluded.columns):
        raise ValueError("Excluded sample audit requires sample_id and exclusion_reason")
    if excluded["sample_id"].isin(frame["sample_id"]).any():
        raise ValueError("Excluded sample IDs remain in report measurements")
    if excluded["sample_id"].eq("").any() or excluded["exclusion_reason"].str.strip().eq("").any():
        raise ValueError("Excluded sample audit requires nonempty sample IDs and reasons")
    columns = [name for name in ("sample_id", "acquisition_id", "subject_id", "exclusion_reason")
               if name in excluded]
    return ("<h2>Explicit sample exclusions</h2>"
            f"<p>{excluded['sample_id'].nunique()} excluded sample IDs; {len(excluded)} excluded acquisition "
            "records. These records are separate from the active-sample counts above. Excluded samples "
            "have no rows in the numerical results and are not represented by fabricated zeros.</p>"
            + _table(excluded[columns]))


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


def _identity_section(root: Path) -> str:
    """Show reference evidence without treating a relative-order check as identity proof."""
    sections = [
        "<h2>Bead identity reference and setup-control audit</h2>",
        "<p>The manufacturer's Figure 2 assigns 39 bead IDs by their relative FITC/PE layout. "
        "The schematic uses ordinal levels, not fluorescence coordinates or calibrated gate boundaries. "
        "Compare it with the measured setup-control diagnostics. Manual reference validation of the "
        "derived bead identities and gate placement uses that setup-control data: a consistent relative-order audit does not "
        "prove true bead identity.</p>",
    ]
    reference = root / "bead-identity-reference.png"
    if reference.is_file():
        encoded = base64.b64encode(reference.read_bytes()).decode("ascii")
        sections.append(
            '<img alt="Figure 2 ordinal bead identity reference, not fluorescence gate coordinates" '
            f'src="data:image/png;base64,{encoded}">'
        )
    audit_path = root / "bead_identity_audit.csv"
    if audit_path.is_file():
        try:
            audit = pd.read_csv(audit_path, dtype=str, keep_default_na=False)
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
            raise ValueError(f"Cannot read bead_identity_audit.csv: {exc}") from exc
        needed = {"bead_id", "marker", "figure2_fitc_level", "figure2_pe_level", "n_events",
                  "observed_fitc_median", "observed_pe_median", "layout_status", "layout_flags"}
        missing = sorted(needed.difference(audit.columns))
        if missing:
            raise ValueError("Bead identity audit missing columns: " + ", ".join(missing))
        unexpected = set(audit["layout_status"]) - {"consistent", "inconsistent", "missing"}
        if unexpected:
            raise ValueError("Bead identity audit has unsupported layout_status values")
        flagged = int((audit["layout_status"] != "consistent").sum())
        sections.append(
            f"<p>{len(audit)} setup-control population records; {flagged} have inconsistent or missing "
            "relative-layout evidence. This check evaluates relative-order plausibility only.</p>"
        )
        sections.append(_table(audit))
    else:
        sections.append("<p>No setup-control relative-order audit was supplied for this run.</p>")
    return "\n".join(sections)


def _source_link(url: object) -> str:
    """Make a source URL clickable only when it is a conventional HTTPS URL."""
    value = str(url)
    try:
        parsed = urlsplit(value)
        safe = (parsed.scheme == "https" and bool(parsed.hostname)
                and parsed.username is None and parsed.password is None
                and not re.search(r"[\x00-\x20\x7f\\]", value))
    except ValueError:
        safe = False
    label = escape(value)
    if safe:
        return f'<a href="{escape(value, quote=True)}" rel="noopener noreferrer">{label}</a>'
    return label


def _annotation_sources(root: Path) -> str:
    path = root / "marker_annotation_sources.json"
    if not path.is_file():
        return "<p>No biological-annotation source registry was supplied for this run.</p>"
    registry = _read_json(path, {})
    if not isinstance(registry, dict) or not isinstance(registry.get("sources"), dict):
        raise ValueError("marker_annotation_sources.json must contain a sources object")
    content = ["<h3>Biological-annotation sources</h3>"]
    if registry.get("annotation_version"):
        content.append(f"<p>Annotation version: {escape(str(registry['annotation_version']))}.</p>")
    if registry.get("curation_notice"):
        content.append(f"<p>{escape(str(registry['curation_notice']))}</p>")
    content.append("<ul>")
    for identifier, source in registry["sources"].items():
        if not isinstance(source, dict):
            raise ValueError("Each biological-annotation source must be an object")
        content.append(f"<li><strong>{escape(str(identifier))}</strong>: "
                       f"{escape(str(source.get('title', 'Untitled source')))}")
        if source.get("url"):
            content.append(" — " + _source_link(source["url"]))
        for key, label in (("doi", "DOI"), ("evidence_scope", "Evidence scope"),
                           ("verification_status", "Source verification")):
            if source.get(key):
                content.append(f"<p>{label}: {escape(str(source[key]))}</p>")
        content.append("</li>")
    content.append("</ul>")
    return "\n".join(content)


def _biological_section(root: Path, frame: pd.DataFrame, *, background_mode: str = "assay_blank") -> str:
    """Organize individual marker columns using a saved descriptive annotation snapshot."""
    path = root / "marker_annotations.csv"
    if not path.is_file():
        return ("<h2>Biological marker views</h2><p>No biological-annotation snapshot was supplied "
                "for this run. Full marker profiles remain available above.</p>")
    try:
        annotations = pd.read_csv(path, dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ValueError(f"Cannot read marker_annotations.csv: {exc}") from exc
    required = {"bead_id", "marker", "role", *_ANNOTATION_COLUMNS}
    missing = sorted(required.difference(annotations.columns))
    if missing:
        raise ValueError("Biological annotations missing columns: " + ", ".join(missing))
    from .config import EXPECTED_MARKERS
    if (len(annotations) != 39 or annotations.bead_id.duplicated().any()
            or set(annotations.bead_id) != set(EXPECTED_MARKERS)):
        raise ValueError("Biological annotation snapshot must identify all 39 unique Neuro bead IDs")
    for row in annotations.to_dict("records"):
        control = row["bead_id"] in {"97", "99"}
        if row["marker"] != EXPECTED_MARKERS[row["bead_id"]]:
            raise ValueError("Biological annotation snapshot has a bead/marker identity mismatch")
        if row["role"] != ("isotype_control" if control else "capture_marker"):
            raise ValueError("Biological annotation snapshot has a bead/role mismatch")
        if not row["biological_group"] or (row["biological_group"] == "isotype_control") != control:
            raise ValueError("Isotype controls cannot be assigned to a biological marker view")
    snapshot = annotations.set_index("bead_id")
    enriched = frame.copy(deep=True)
    for row in enriched[["bead_id", "marker", "role"]].drop_duplicates().to_dict("records"):
        if row["bead_id"] not in snapshot.index or any(
            row[field] != snapshot.loc[row["bead_id"], field] for field in ("marker", "role")
        ):
            raise ValueError("Measurements do not match their biological annotation snapshot")
    # Read the saved snapshot, not a newer installed resource, when rerendering.
    # Existing per-measurement annotations must agree rather than be silently replaced.
    for column in _ANNOTATION_COLUMNS:
        expected = enriched["bead_id"].map(snapshot[column])
        if column in enriched and not enriched[column].fillna("").astype(str).equals(expected.astype(str)):
            raise ValueError(f"Measurements and biological annotation snapshot disagree on {column}")
        enriched[column] = expected
    try:
        from .annotations import GROUP_LABELS
    except ImportError:
        GROUP_LABELS = {}
    corrected_label = "Setup-baseline-subtracted APC" if background_mode == "setup_baseline" else "background-corrected APC"
    sections = [
        "<h2>Biological marker views</h2>",
        "<p>These descriptive associations organize individual markers for review. Association tags "
        "can overlap across cell types and do not identify an exclusive cell of origin. Each heatmap "
        f"column remains a single capture marker with signed {corrected_label}; no group sums, "
        "averages, EV fractions, or additional detection thresholds are calculated. The full 37-marker "
        "profiles remain available above. Isotype controls are excluded from these biological views.</p>",
        "<p>CD171/L1CAM can occur as soluble protein in human plasma, so its capture signal does not "
        "establish neuronal EVs. CD56 is also expressed by NK cells. CD140a/PDGFRA and CSPG4 can "
        "associate with mesenchymal or pericyte populations. CD11b, CX3CR1, and CD68 also occur in "
        "peripheral myeloid cells; their plasma signals do not count microglial EVs.</p>",
    ]
    capture = enriched.loc[enriched["role"] == "capture_marker"]
    present = set(annotations.loc[annotations.role == "capture_marker", "biological_group"])
    order = list(_BIOLOGICAL_GROUP_ORDER) + sorted(present - set(_BIOLOGICAL_GROUP_ORDER))
    for group in order:
        if group not in present:
            continue
        selected = capture.loc[capture["biological_group"] == group]
        label = str(GROUP_LABELS.get(group, group.replace("_", " ").capitalize()))
        ids = set(annotations.loc[(annotations.role == "capture_marker") &
                                   (annotations.biological_group == group), "bead_id"])
        heading = f"{label} — individual markers ({len(ids)})"
        if background_mode == "setup_baseline":
            heading = f"{label} — Setup-baseline-subtracted APC — individual markers ({len(ids)})"
        view = _heatmap(selected, "corrected_apc", heading)
        downloads = []
        if group in _BIOLOGICAL_GROUP_ORDER:
            for scale in ("background_corrected", "normalized"):
                filename = f"marker_views/{group}_{scale}.csv"
                if (root / filename).is_file():
                    scale_label = ("Setup-baseline-subtracted APC" if background_mode == "setup_baseline"
                                   and scale == "background_corrected" else scale)
                    downloads.append(f'<a href="{escape(filename, quote=True)}">{escape(scale_label)} CSV</a>')
        if downloads:
            view += "<p>Individual-marker tables: " + " · ".join(downloads) + "</p>"
        if group in _BIOLOGICAL_FOCUS:
            sections.append(view)
        else:
            sections.append(f"<details><summary>{escape(heading)}</summary>{view}</details>")
    sections.append("<h3>Marker annotation snapshot</h3><p>All 39 populations are listed, including "
                    "the two diagnostic isotype controls. These associations annotate measurements "
                    "without changing their numerical values.</p>")
    sections.append(_table(annotations))
    sections.append(_annotation_sources(root))
    return "\n".join(sections)


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
    background_mode, background_rationale = _background_context(frame, provenance)
    setup_baseline = background_mode == "setup_baseline"
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
    interpretation = (
        "<p>Setup-baseline-subtracted APC is the sample median minus the chosen kit setup-bead median. "
        "Kit setup beads alone do not measure detection-reagent or cocktail background. "
        "Detection-reagent background is unmeasured by this reference, so below/at/above setup baseline "
        "describes the numerical comparison only and does not establish detection or nondetection. "
        "Negative signed values are retained. Optional tetraspanin normalization does not remove this "
        "limitation; raw APC values remain available. No validated LOD/LOQ or inferential statistics "
        "are supplied by this report. Missing/invalid measurements remain unavailable, not zero.</p>"
        if setup_baseline else
        "<p>Signed background-corrected values are retained. Negative corrected values are labeled "
        "not detected; zero is at background. A positive corrected value alone is not a validated "
        "detection call. No validated LOD/LOQ or inferential statistics are supplied by this report. "
        "Missing/invalid measurements are unavailable, not zero. Tetraspanin normalization is optional "
        "and changes the biological question.</p>"
    )
    baseline_banner = (
        '<p class="banner">SETUP-BEAD BASELINE ONLY — detection-reagent background is unmeasured; '
        'the subtracted values do not establish detection or nondetection.</p>'
        f"<p>Reference-selection rationale: {escape(background_rationale or 'Not recorded')}.</p>"
        if setup_baseline else ""
    )
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
              baseline_banner,
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
              interpretation,
              "<h2>Marker profiles</h2>"]
    pieces.extend(_heatmap(frame, key, title) for key, title in (
        ("raw_median_apc", "Raw APC medians — QC-valid measurements"),
        ("corrected_apc", "Setup-baseline-subtracted APC" if setup_baseline else "Background-corrected APC"),
        ("normalized_apc", "Setup-baseline-subtracted APC — optional tetraspanin normalization"
         if setup_baseline else "Optional tetraspanin-normalized signal"),
    ))
    pieces.append(_biological_section(root, frame, background_mode=background_mode))
    pieces.append(_exclusion_section(root, frame, provenance))
    pieces.append(_identity_section(root))
    pieces.append("<h2>Quality control</h2>")
    pieces.append(
        "<p>An acquisition target of 5,000 total bead events per sample is distributed across all 39 "
        "populations: about 128 per population only if evenly distributed, before filtering losses. "
        "Actual accepted counts are recorded per population in n_events and through the gate stages. "
        "The total does not guarantee sufficient counts in every population. No fixed count cutoff "
        "is inferred from that target; use documented, pilot-defined QC thresholds.</p>"
    )
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
    pieces.append("<h2>Isotype controls</h2><p>Isotype controls diagnose nonspecific binding. " + (
        "Their displayed values use the same setup-bead-only reference; detection-reagent "
        "background remains unmeasured.</p>" if setup_baseline else
        "They do not replace the marker-specific assay blank and are not automatically "
        "subtracted a second time.</p>"))
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
                 "normalized.csv", "isotype_controls.csv", "qc.csv", "sample_metadata.csv",
                 "bead_identity_audit.csv", "bead-identity-reference.csv", "marker_annotations.csv",
                 "marker_annotation_sources.json", "excluded_samples.csv"):
        if (root / name).is_file():
            label = "Setup-baseline-subtracted APC CSV (background_corrected.csv)" if setup_baseline and name == "background_corrected.csv" else name
            pieces.append(f"<li><a href=\"{escape(name, quote=True)}\">{escape(label)}</a></li>")
    pieces.append("</ul><h2>Run provenance</h2><pre>")
    pieces.append(escape(json.dumps(provenance, indent=2, ensure_ascii=False)))
    pieces.append("</pre></body></html>")
    output = root / "report.html"
    output.write_text("\n".join(pieces), encoding="utf-8")
    return output
