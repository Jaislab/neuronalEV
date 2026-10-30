"""Report checks exercise plotting, missing values, and untrusted identifiers."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from facs_neuro.report import _biological_section, render_report


EXPECTED_GROUPS = {
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


class ReportTests(unittest.TestCase):
    def _write_annotations(self, root: Path) -> pd.DataFrame:
        # Synthetic annotation metadata exercises the report independently of
        # the packaged curation resource; memberships follow the stated design.
        from facs_neuro.config import EXPECTED_MARKERS
        rows = []
        for group, ids in EXPECTED_GROUPS.items():
            for bead in sorted(ids):
                rows.append({
                    "bead_id": bead, "marker": EXPECTED_MARKERS[bead],
                    "role": "isotype_control" if bead in {"97", "99"} else "capture_marker",
                    "biological_group": group, "association_tags": "neuron;immune_nk" if bead == "33" else group,
                    "gene_or_epitope": EXPECTED_MARKERS[bead],
                    "interpretation_note": "Synthetic annotation <script>bad()</script>",
                    "annotation_source_ids": "safe;unsafe", "annotation_version": "test-v1",
                })
        annotations = pd.DataFrame(rows)
        annotations.to_csv(root / "marker_annotations.csv", index=False)
        (root / "marker_annotation_sources.json").write_text(json.dumps({
            "annotation_version": "test-v1", "curation_notice": "TEST associations only <script>notice</script>",
            "sources": {
                "safe": {"title": "Source <script>title</script>",
                         "url": 'https://example.org/article?note="quoted"',
                         "evidence_scope": "Not exclusive <img src=x onerror=bad()>",
                         "verification_status": "test-only"},
                "unsafe": {"title": "Unsafe URL", "url": "javascript:alert(1)"},
            },
        }))
        return annotations

    def _write_run(self, root: Path, *, synthetic: bool = True) -> bytes:
        rows = []
        for index, sample in enumerate(("<script>alert(1)</script>", "2")):
            for bead, marker, role in (("66", "CD9", "capture_marker"),
                                       ("67", "CD63", "capture_marker"),
                                       ("97", "REA control", "isotype_control")):
                rows.append({
                    "sample_id": sample, "acquisition_id": f"a{index + 1}",
                    "subject_id": f"donor_{index + 1}", "bead_id": bead,
                    "marker": marker, "role": role, "n_events": 120,
                    "raw_median_apc": 20 + index, "corrected_apc": -0.5 if index == 0 else 2.5,
                    "normalized_apc": np.nan, "measurement_valid": True,
                    "qc_flags": "flag <img src=x onerror=alert(1)>" if bead == "97" else "",
                    "correction_status": "corrected", "detection_state": "not_detected",
                    "normalization_status": "disabled",
                })
        csv_file = root / "measurements_long.csv"
        pd.DataFrame(rows).to_csv(csv_file, index=False)
        (root / "qc.json").write_text(json.dumps([
            {"record_id": "a1<script>", "file_type": "sample", "total_events": 500,
             "pregated_events": 480, "assigned_events": 400, "ambiguous_events": 0,
             "unassigned_events": 80}
        ]))
        (root / "provenance.json").write_text(json.dumps({
            "dataset_kind": "synthetic" if synthetic else "experimental",
            "validation_status": "synthetic_only" if synthetic else "unvalidated_assay",
            "config_sha256": "example hash", "notes": "<script>provenance</script>",
        }))
        return csv_file.read_bytes()

    def test_synthetic_report_plots_and_preserves_signed_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = self._write_run(root)
            output = render_report(root)
            html = output.read_text()
            self.assertEqual(output, root / "report.html")
            self.assertIn("SYNTHETIC DATA", html)
            self.assertIn("software demonstration only", html)
            self.assertIn("Normalization is unavailable", html)
            self.assertIn("data:image/png;base64,", html)
            self.assertIn("No validated LOD/LOQ or inferential statistics", html)
            self.assertIn("not biological zeros", html)
            self.assertIn("-0.5", html)
            self.assertEqual((root / "measurements_long.csv").read_bytes(), original)
            self.assertNotIn("src=\"http", html)
            self.assertNotIn("<script>", html)
            self.assertIn("&lt;script&gt;alert(1)&lt;/script&gt;", html)
            self.assertIn("&lt;img src=x onerror=alert(1)&gt;", html)
            self.assertIn("&lt;script&gt;provenance&lt;/script&gt;", html)
            self.assertIn("No biological-annotation snapshot was supplied", html)

    def test_experimental_report_cannot_claim_assay_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root, synthetic=False)
            html = render_report(root).read_text()
            self.assertIn("EXPERIMENTAL DATA", html)
            self.assertIn("require reviewed real-file validation", html)
            self.assertNotIn("SYNTHETIC DATA", html)

    def test_duplicate_measurements_are_not_silently_averaged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            frame = pd.read_csv(root / "measurements_long.csv")
            frame = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
            frame.to_csv(root / "measurements_long.csv", index=False)
            with self.assertRaisesRegex(ValueError, "Duplicate acquisition/bead"):
                render_report(root)

    def test_all_missing_signal_generates_explanations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            frame = pd.read_csv(root / "measurements_long.csv")
            for name in ("raw_median_apc", "corrected_apc", "normalized_apc"):
                frame[name] = np.nan
            frame["measurement_valid"] = False
            frame.to_csv(root / "measurements_long.csv", index=False)
            html = render_report(root).read_text()
            self.assertIn("No valid finite measurements", html)
            self.assertIn("Normalization is unavailable", html)

    def test_missing_required_schema_is_clear(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pd.DataFrame({"sample_id": ["1"]}).to_csv(root / "measurements_long.csv", index=False)
            with self.assertRaisesRegex(ValueError, "missing columns"):
                render_report(root)

    def test_identifiers_are_not_parsed_as_math(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            frame = pd.read_csv(root / "measurements_long.csv")
            frame["acquisition_id"] = frame["acquisition_id"].replace({"a1": r"$\invalidcommand$"})
            frame.loc[frame["marker"] == "CD9", "marker"] = r"$\invalidmarker$"
            frame.to_csv(root / "measurements_long.csv", index=False)
            self.assertIn("data:image/png;base64,", render_report(root).read_text())

    def test_direct_api_rejects_incomplete_run_without_replacing_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            path = root / "provenance.json"
            provenance = json.loads(path.read_text())
            provenance["status"] = "incomplete"
            path.write_text(json.dumps(provenance))
            existing = root / "report.html"
            existing.write_text("existing report must be preserved")
            with self.assertRaisesRegex(ValueError, "Cannot report an incomplete run"):
                render_report(root)
            self.assertEqual(existing.read_text(), "existing report must be preserved")

    def test_optional_identity_audit_and_schematic_are_embedded_with_caveats(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            audit = pd.DataFrame([{
                "record_id": "setup<script>", "setup_id": "fortessa", "bead_id": "22",
                "marker": "ErbB-2", "role": "capture_marker", "figure2_fitc_level": 2,
                "figure2_pe_level": 2, "n_events": 128, "observed_fitc_median": 1500,
                "observed_pe_median": 2500, "layout_status": "inconsistent",
                "layout_flags": "flag <img src=x onerror=alert(1)>",
            }])
            audit.to_csv(root / "bead_identity_audit.csv", index=False)
            audit[["bead_id", "marker", "figure2_fitc_level", "figure2_pe_level"]].to_csv(
                root / "bead-identity-reference.csv", index=False)
            from matplotlib.figure import Figure
            figure = Figure(figsize=(2, 1))
            figure.subplots().text(.1, .5, "Ordinal schematic")
            figure.savefig(root / "bead-identity-reference.png")
            html = render_report(root).read_text()
            self.assertIn("Figure 2 ordinal bead identity reference", html)
            self.assertIn("1 setup-control population records; 1 have inconsistent or missing", html)
            self.assertIn("relative-order audit does not prove true bead identity", html)
            self.assertIn("Manual reference validation", html)
            self.assertIn("not fluorescence coordinates or calibrated gate boundaries", html)
            self.assertIn("5,000 total bead events", html)
            self.assertIn("128 per population only if evenly distributed", html)
            self.assertIn("No fixed count cutoff", html)
            self.assertIn('href="bead_identity_audit.csv"', html)
            self.assertIn('href="bead-identity-reference.csv"', html)
            self.assertIn("setup&lt;script&gt;", html)
            self.assertNotIn("<script>", html)

    def test_identity_audit_missing_schema_is_not_presented_as_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            pd.DataFrame({"bead_id": ["22"]}).to_csv(root / "bead_identity_audit.csv", index=False)
            with self.assertRaisesRegex(ValueError, "Bead identity audit missing columns"):
                render_report(root)

    def test_biological_views_select_exact_individual_markers_and_exclude_isotypes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            annotations = self._write_annotations(root)
            frame = annotations[["bead_id", "marker", "role"]].copy()
            frame["sample_id"] = "1"
            frame["acquisition_id"] = "1"
            frame["corrected_apc"] = 4.5
            frame.loc[frame.bead_id == "32", "corrected_apc"] = np.nan
            frame.loc[frame.bead_id == "24", "corrected_apc"] = -1.25
            original = frame.copy(deep=True)
            views = {}

            def inspect_columns(selected, value, title):
                groups = set(selected.biological_group)
                self.assertEqual(len(groups), 1)
                group = next(iter(groups))
                self.assertEqual(value, "corrected_apc")
                self.assertEqual(set(selected.role), {"capture_marker"})
                self.assertNotIn("97", set(selected.bead_id))
                self.assertNotIn("99", set(selected.bead_id))
                views[group] = selected.copy(deep=True)
                return f"<h3>{title}</h3>"

            with patch("facs_neuro.report._heatmap", side_effect=inspect_columns):
                html = _biological_section(root, frame)
            self.assertEqual(set(views), set(EXPECTED_GROUPS) - {"isotype_control"})
            for group, selected in views.items():
                self.assertEqual(set(selected.bead_id), EXPECTED_GROUPS[group])
            self.assertEqual(sum(len(view) for view in views.values()), 37)
            neuron = views["neuron_associated"].set_index("bead_id")
            glial = views["glial_progenitor_associated"].set_index("bead_id")
            self.assertTrue(pd.isna(neuron.loc["32", "corrected_apc"]))
            self.assertEqual(glial.loc["24", "corrected_apc"], -1.25)
            self.assertEqual(neuron.loc["33", "association_tags"], "neuron;immune_nk")
            pd.testing.assert_frame_equal(frame, original)
            self.assertIn("can overlap across cell types", html)
            self.assertIn("CD171/L1CAM can occur as soluble protein", html)
            self.assertIn("CD56 is also expressed by NK cells", html)
            self.assertIn("do not count microglial EVs", html)
            self.assertIn("no group sums", html)
            self.assertIn("All 39 populations", html)
            self.assertIn("&lt;script&gt;bad()&lt;/script&gt;", html)

    def test_biological_sources_are_safe_and_numerical_csv_is_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = self._write_run(root)
            self._write_annotations(root)
            html = render_report(root).read_text()
            self.assertEqual((root / "measurements_long.csv").read_bytes(), original)
            self.assertIn("Biological marker views", html)
            self.assertIn("Biological-annotation sources", html)
            self.assertIn("Source &lt;script&gt;title&lt;/script&gt;", html)
            self.assertIn('href="https://example.org/article?note=&quot;quoted&quot;"', html)
            self.assertNotIn('href="javascript:', html)
            self.assertNotIn("<script>", html)
            self.assertIn('href="marker_annotations.csv"', html)
            self.assertIn('href="marker_annotation_sources.json"', html)

    def test_isotype_cannot_be_mislabeled_as_neuron_associated(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            annotations = self._write_annotations(root)
            annotations.loc[annotations.bead_id == "97", "biological_group"] = "neuron_associated"
            annotations.to_csv(root / "marker_annotations.csv", index=False)
            with self.assertRaisesRegex(ValueError, "Isotype controls"):
                render_report(root)

    def test_measurement_annotation_must_match_saved_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            self._write_annotations(root)
            frame = pd.read_csv(root / "measurements_long.csv", dtype=str, keep_default_na=False)
            frame["biological_group"] = "neuron_associated"
            frame.to_csv(root / "measurements_long.csv", index=False)
            with self.assertRaisesRegex(ValueError, "snapshot disagree on biological_group"):
                render_report(root)

    def test_setup_only_baseline_keeps_limitation_in_titles_and_interpretation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            self._write_annotations(root)
            path = root / "measurements_long.csv"
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            frame["background_mode"] = "setup_baseline"
            frame["background_reference_id"] = "setup"
            frame["background_limitation"] = "setup_baseline_only_reagent_background_unmeasured"
            frame["qc_flags"] = "setup_baseline_only"
            frame["detection_state"] = np.where(frame.corrected_apc.astype(float) < 0,
                                                  "below_setup_baseline", "above_setup_baseline")
            frame.to_csv(path, index=False)
            original = path.read_bytes()
            provenance_path = root / "provenance.json"
            provenance = json.loads(provenance_path.read_text())
            provenance["background"] = {"mode": "setup_baseline", "rationale": "Setup only <script>reason</script>"}
            provenance_path.write_text(json.dumps(provenance))
            html = render_report(root).read_text()
            self.assertIn("SETUP-BEAD BASELINE ONLY", html)
            self.assertIn("Setup-baseline-subtracted APC", html)
            self.assertIn("detection-reagent background is unmeasured", html)
            self.assertIn("Optional tetraspanin normalization does not remove", html)
            self.assertIn("does not establish detection or nondetection", html)
            self.assertIn("&lt;script&gt;reason&lt;/script&gt;", html)
            self.assertNotIn("Negative corrected values are labeled not detected", html)
            self.assertNotIn("<h3>Background-corrected APC</h3>", html)
            self.assertNotIn("single capture marker with signed background-corrected APC", html)
            self.assertEqual(path.read_bytes(), original)

    def test_background_mode_conflicts_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            path = root / "measurements_long.csv"
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            frame["background_mode"] = "setup_baseline"
            frame.to_csv(path, index=False)
            provenance_path = root / "provenance.json"
            provenance = json.loads(provenance_path.read_text())
            provenance["background"] = {"mode": "assay_blank"}
            provenance_path.write_text(json.dumps(provenance))
            with self.assertRaisesRegex(ValueError, "disagree on background mode"):
                render_report(root)

    def test_explicit_exclusions_are_separate_escaped_audit_records(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            path = root / "measurements_long.csv"
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            frame = frame[frame.sample_id != "2"]
            frame.to_csv(path, index=False)
            original = path.read_bytes()
            pd.DataFrame([{"sample_id": "2", "acquisition_id": "a2", "subject_id": "donor_2",
                           "exclusion_reason": "Bead loss <script>reason</script>"}]).to_csv(
                root / "excluded_samples.csv", index=False)
            html = render_report(root).read_text()
            self.assertIn("Explicit sample exclusions", html)
            self.assertIn("1 excluded sample IDs; 1 excluded acquisition records", html)
            self.assertIn("Bead loss &lt;script&gt;reason&lt;/script&gt;", html)
            self.assertIn("not represented by fabricated zeros", html)
            self.assertIn('href="excluded_samples.csv"', html)
            self.assertIn("1 sample IDs; 1 acquisitions; 1 donor/subject IDs", html)
            self.assertEqual(path.read_bytes(), original)

    def test_excluded_samples_cannot_remain_in_measurement_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._write_run(root)
            pd.DataFrame([{"sample_id": "2", "exclusion_reason": "Bead loss"}]).to_csv(
                root / "excluded_samples.csv", index=False)
            with self.assertRaisesRegex(ValueError, "Excluded sample IDs remain"):
                render_report(root)


if __name__ == "__main__":
    unittest.main()
