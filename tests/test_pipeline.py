"""End-to-end checks use independently generated, explicitly synthetic FCS data."""

from __future__ import annotations

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from io import StringIO
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from facs_neuro.cli import main
from facs_neuro.config import ConfigurationError, file_hash
from facs_neuro.demo import generate_demo, write_fcs
from facs_neuro.fcs import read_fcs
from facs_neuro.pipeline import analyze, inspect


class PipelineIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.fixture = generate_demo(cls.root / "fixture", donors=2, run=False)
        cls.results = analyze(cls.fixture / "analysis.json", cls.root / "results")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @contextmanager
    def altered(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / "fixture"
            shutil.copytree(self.fixture, fixture)
            yield fixture, root / "results"

    def test_two_donor_analysis_preserves_identity_signed_values_and_counts(self):
        measurements = pd.read_csv(self.results / "measurements_long.csv",
                                   dtype={"sample_id": str, "subject_id": str,
                                          "acquisition_id": str, "bead_id": str})
        self.assertEqual(len(measurements), 2 * 39)
        self.assertEqual(set(measurements.sample_id), {"1", "2"})
        self.assertEqual(set(measurements.subject_id), {"1", "2"})
        self.assertTrue(measurements.groupby("acquisition_id").size().eq(39).all())
        self.assertTrue(measurements.n_events.eq(80).all())
        self.assertTrue(measurements.measurement_valid.all())
        self.assertEqual(len(measurements[measurements.role == "capture_marker"]), 2 * 37)
        first = measurements[(measurements.sample_id == "1") & (measurements.bead_id == "22")].iloc[0]
        self.assertEqual(first.marker, "ErbB-2")
        self.assertLess(first.corrected_apc, 0)
        self.assertEqual(first.detection_state, "not_detected")
        self.assertLess(first.normalized_apc, 0)

        # Independently select the known synthetic bead-22 interval. This checks
        # identity and numerical agreement, not merely the shape of output files.
        source_medians = []
        for filename in ("1.fcs", "blank.fcs"):
            source = read_fcs(self.fixture / "fcs" / filename)
            fitc = source.events[:, source.channels.index("FITC-A")]
            pe = source.events[:, source.channels.index("PE-A")]
            selected = ((fitc > 1750) & (fitc < 2250) & (pe > 1750) & (pe < 2250))
            self.assertEqual(int(selected.sum()), 80)
            source_medians.append(float(np.median(source.events[selected, source.channels.index("APC-A")])))
        self.assertAlmostEqual(first.raw_median_apc, source_medians[0], places=8)
        self.assertAlmostEqual(first.blank_median_apc, source_medians[1], places=8)
        self.assertAlmostEqual(first.corrected_apc, source_medians[0] - source_medians[1], places=8)
        for filename in ("raw_median_apc.csv", "background_corrected.csv", "normalized.csv"):
            matrix = pd.read_csv(self.results / filename, index_col=0)
            self.assertEqual(matrix.shape, (2, 37))
            self.assertNotIn("REA control", matrix.columns)
            self.assertNotIn("mIgG1 control", matrix.columns)
        provenance = json.loads((self.results / "provenance.json").read_text())
        self.assertEqual(provenance["status"], "completed")
        self.assertEqual(provenance["dataset_kind"], "synthetic")
        self.assertEqual(provenance["validation_status"], "synthetic_only")
        self.assertEqual(provenance["config_sha256"], file_hash(self.fixture / "analysis.json"))
        self.assertEqual(len(provenance["inputs"]), 4)
        for entry in provenance["inputs"]:
            self.assertEqual(entry["sha256"], file_hash(Path(entry["path"])))
        self.assertEqual(len(list((self.results / "diagnostics").glob("*.png"))), 4)
        self.assertIn("SYNTHETIC DATA", (self.results / "report.html").read_text())
        self.assertFalse((self.results / "failure.json").exists())

    def test_missing_blank_rejected_before_output_is_created(self):
        with self.altered() as (fixture, output):
            controls = pd.read_csv(fixture / "controls.csv", dtype=str)
            controls = controls[controls.control_type != "assay_blank"]
            controls.to_csv(fixture / "controls.csv", index=False)
            with self.assertRaisesRegex(ConfigurationError, "blank_id.*assay_blank"):
                analyze(fixture / "analysis.json", output)
            self.assertFalse(output.exists())

    def test_blank_kit_lot_mismatch_rejected_before_output_is_created(self):
        with self.altered() as (fixture, output):
            controls = pd.read_csv(fixture / "controls.csv", dtype=str)
            controls.loc[controls.control_type == "assay_blank", "kit_lot"] = "OTHER-LOT"
            controls.to_csv(fixture / "controls.csv", index=False)
            with self.assertRaisesRegex(ConfigurationError, "kit_lot differs"):
                analyze(fixture / "analysis.json", output)
            self.assertFalse(output.exists())

    def test_manifest_metadata_cannot_override_computed_measurements(self):
        for filename, column, value in (("samples.csv", "raw_median_apc", "0"),
                                         ("samples.csv", "corrected_apc", "100000"),
                                         ("controls.csv", "measurement_valid", "True")):
            with self.subTest(manifest=filename, column=column):
                with self.altered() as (fixture, output):
                    frame = pd.read_csv(fixture / filename, dtype=str, keep_default_na=False)
                    frame[column] = value
                    frame.to_csv(fixture / filename, index=False)
                    with self.assertRaisesRegex(ConfigurationError, "reserved columns"):
                        analyze(fixture / "analysis.json", output)
                    self.assertFalse(output.exists())

    def test_nonnormalizer_marker_label_swap_cannot_change_bead_identity(self):
        with self.altered() as (fixture, output):
            path = fixture / "panel.csv"
            panel = pd.read_csv(path, dtype=str, keep_default_na=False)
            original = panel.set_index("bead_id").marker
            panel.loc[panel.bead_id == "22", "marker"] = original["23"]
            panel.loc[panel.bead_id == "23", "marker"] = original["22"]
            self.assertEqual(panel.marker.nunique(), 39)
            self.assertEqual(panel.loc[panel.bead_id == "66", "marker"].iloc[0], "CD9")
            panel.to_csv(path, index=False)
            with self.assertRaisesRegex(ConfigurationError, "bead identity mismatch"):
                analyze(fixture / "analysis.json", output)
            self.assertFalse(output.exists())

    def test_unreviewed_gates_cannot_generate_a_completed_report(self):
        with self.altered() as (fixture, output):
            path = fixture / "analysis.json"
            config = json.loads(path.read_text())
            config["setups"]["synthetic_setup"]["reviewed"] = False
            path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, "review"):
                analyze(path, output)
            if output.exists():
                self.assertEqual(json.loads((output / "provenance.json").read_text())["status"], "incomplete")
                self.assertTrue((output / "failure.json").is_file())
                self.assertFalse((output / "report.html").exists())

    def test_changed_detector_metadata_requires_separate_setup(self):
        with self.altered() as (fixture, output):
            path = fixture / "fcs" / "1.fcs"
            content = path.read_bytes()
            token = b"$P1R|262144|"
            self.assertIn(token, content)
            path.write_bytes(content.replace(token, b"$P1R|262143|", 1))
            with self.assertRaisesRegex(ConfigurationError, "metadata differs within setup"):
                analyze(fixture / "analysis.json", output)
            self.assertEqual(json.loads((output / "provenance.json").read_text())["status"], "incomplete")
            self.assertFalse((output / "report.html").exists())

    def test_technical_acquisition_keeps_same_donor_identity(self):
        with self.altered() as (fixture, _):
            source = read_fcs(fixture / "fcs" / "1.fcs")
            events = source.events.copy()
            events[:, source.channels.index("APC-A")] += 0.25
            write_fcs(fixture / "fcs" / "1-tech.fcs", events, label="SYNTHETIC-TECHNICAL-ACQUISITION")
            samples = pd.read_csv(fixture / "samples.csv", dtype=str, keep_default_na=False)
            technical = samples.iloc[[0]].copy()
            technical["acquisition_id"] = "1-tech"
            technical["fcs_path"] = "fcs/1-tech.fcs"
            pd.concat([samples, technical], ignore_index=True).to_csv(fixture / "samples.csv", index=False)
            inventory = inspect(fixture / "analysis.json")
            self.assertTrue(inventory["inputs_ready"])
            self.assertEqual(inventory["acquisitions"], 3)
            self.assertEqual(inventory["donors"], 2)

    def test_existing_result_directory_is_preserved(self):
        original = (self.results / "report.html").read_bytes()
        with self.assertRaisesRegex(ConfigurationError, "Output already exists"):
            analyze(self.fixture / "analysis.json", self.results)
        self.assertEqual((self.results / "report.html").read_bytes(), original)

    def test_real_template_inspection_reports_24_unready_donors(self):
        real_config = Path(__file__).resolve().parents[1] / "config" / "analysis.json"
        inventory = inspect(real_config)
        self.assertEqual(inventory["dataset_kind"], "experimental")
        self.assertEqual(inventory["acquisitions"], 24)
        self.assertEqual(inventory["donors"], 24)
        self.assertEqual(inventory["capture_markers"], 37)
        self.assertFalse(inventory["inputs_ready"])
        self.assertFalse(inventory["gates_marked_reviewed"])
        with redirect_stdout(StringIO()) as stdout, redirect_stderr(StringIO()):
            code = main(["inspect", "--config", str(real_config)])
        self.assertEqual(code, 2)
        self.assertFalse(json.loads(stdout.getvalue())["inputs_ready"])

    def test_cli_refuses_report_for_incomplete_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "provenance.json").write_text(json.dumps({"status": "incomplete"}))
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()) as stderr:
                code = main(["report", "--run", str(root)])
            self.assertEqual(code, 2)
            self.assertIn("Cannot report an incomplete run", stderr.getvalue())
            self.assertFalse((root / "report.html").exists())


if __name__ == "__main__":
    unittest.main()
