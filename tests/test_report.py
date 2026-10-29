"""Report checks exercise plotting, missing values, and untrusted identifiers."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from facs_neuro.report import render_report


class ReportTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
