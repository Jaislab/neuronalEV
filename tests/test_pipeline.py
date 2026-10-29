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
from facs_neuro.pipeline import analyze, calibrate, inspect


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
        audit = pd.read_csv(self.results / "bead_identity_audit.csv", dtype={"bead_id": str})
        self.assertEqual(len(audit), 39)
        self.assertEqual(set(audit.record_id), {"setup"})
        self.assertEqual(set(audit.layout_status), {"consistent"})
        self.assertTrue((self.results / "bead-identity-reference.png").is_file())
        self.assertIn("SYNTHETIC DATA", (self.results / "report.html").read_text())
        self.assertFalse((self.results / "failure.json").exists())

    def test_annotation_snapshot_covers_all_markers_and_preserves_per_acquisition_rows(self):
        annotations = pd.read_csv(self.results / "marker_annotations.csv", dtype=str, keep_default_na=False)
        measurements = pd.read_csv(self.results / "measurements_long.csv", dtype=str, keep_default_na=False)
        metadata = ["biological_group", "association_tags", "gene_or_epitope", "interpretation_note",
                    "annotation_source_ids", "annotation_version"]
        self.assertEqual(len(annotations), 39)
        self.assertEqual(annotations.bead_id.nunique(), 39)
        self.assertEqual(len(measurements), 78)
        self.assertTrue(set(metadata).issubset(measurements.columns))
        expected_groups = {
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
        self.assertEqual(set(annotations.biological_group), set(expected_groups))
        for group, beads in expected_groups.items():
            self.assertEqual(set(annotations.loc[annotations.biological_group == group, "bead_id"]), beads)
        self.assertEqual(sum(annotations.role == "capture_marker"), 37)
        self.assertEqual(set(annotations.loc[annotations.role == "isotype_control", "biological_group"]),
                         {"isotype_control"})
        snapshot = annotations.set_index("bead_id")
        for column in metadata:
            self.assertTrue(measurements[column].equals(measurements.bead_id.map(snapshot[column])))
        self.assertEqual(measurements.acquisition_id.nunique(), 2)
        registry = json.loads((self.results / "marker_annotation_sources.json").read_text())
        self.assertTrue(registry["sources"])
        provenance = json.loads((self.results / "provenance.json").read_text())
        self.assertEqual(provenance["annotation_version"], registry["annotation_version"])
        self.assertEqual(set(annotations.annotation_version), {registry["annotation_version"]})
        for filename in ("marker_annotations.csv", "marker_annotation_sources.json"):
            self.assertEqual(provenance["annotation_resource_hashes"][filename], file_hash(self.results / filename))
        views = self.results / "marker_views"
        self.assertEqual(len(list(views.glob("*.csv"))), 16)
        total_columns = 0
        for group, beads in expected_groups.items():
            if group == "isotype_control":
                self.assertFalse(list(views.glob("isotype_control_*.csv")))
                continue
            markers = set(snapshot.loc[list(beads), "marker"])
            for scale, overall in (("background_corrected", "background_corrected.csv"),
                                   ("normalized", "normalized.csv")):
                view = pd.read_csv(views / f"{group}_{scale}.csv", index_col=0)
                full = pd.read_csv(self.results / overall, index_col=0)
                self.assertEqual(view.shape, (2, len(beads)))
                self.assertEqual(set(view.columns), markers)
                pd.testing.assert_frame_equal(view.sort_index(axis=1), full[view.columns].sort_index(axis=1))
                if scale == "background_corrected":
                    total_columns += len(view.columns)
        self.assertEqual(total_columns, 37)
        html = (self.results / "report.html").read_text()
        self.assertIn("Biological marker views", html)
        self.assertIn("Association tags can overlap", html)
        self.assertIn("full 37-marker profiles remain available", html)
        self.assertIn('href="marker_annotations.csv"', html)
        self.assertIn('href="marker_annotation_sources.json"', html)

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
                                         ("samples.csv", "biological_group", "neuron_associated"),
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

    def test_populations_with_swapped_gate_positions_fail_figure2_audit(self):
        with self.altered() as (fixture, output):
            path = fixture / "analysis.json"
            config = json.loads(path.read_text())
            populations = {row["bead_id"]: row for row in config["setups"]["synthetic_setup"]["populations"]}
            # Labels remain canonical, and every population still has 80 events.
            # Swapping the actual gates must be caught independently of counts.
            populations["22"]["vertices"], populations["32"]["vertices"] = (
                populations["32"]["vertices"], populations["22"]["vertices"])
            path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ConfigurationError, "Figure 2"):
                analyze(path, output)
            self.assertTrue((output / "failure.json").is_file())
            self.assertEqual(json.loads((output / "provenance.json").read_text())["status"], "incomplete")
            self.assertFalse((output / "report.html").exists())
            audit = pd.read_csv(output / "bead_identity_audit.csv", dtype={"bead_id": str})
            bad_ids = set(audit.loc[audit.layout_status != "consistent", "bead_id"])
            self.assertTrue({"22", "32"}.issubset(bad_ids))
            self.assertTrue(audit.n_events.eq(80).all())

    def test_calibration_outputs_unreviewed_gates_and_ordinal_reference(self):
        with self.altered() as (fixture, output):
            calibrate(fixture / "analysis.json", output)
            templates = json.loads((output / "gate-templates.json").read_text())
            setup = templates["synthetic_setup"]
            self.assertFalse(setup["reviewed"])
            self.assertEqual(setup["gate_version"], "")
            self.assertEqual(len(setup["populations"]), 39)
            self.assertTrue(all(population["vertices"] == [] for population in setup["populations"]))
            reference = pd.read_csv(output / "bead-identity-reference.csv", dtype={"bead_id": str})
            self.assertEqual(len(reference), 39)
            self.assertEqual(set(reference.bead_id), {p["bead_id"] for p in setup["populations"]})
            expected = reference.set_index("bead_id")
            self.assertEqual(expected.loc["22", "figure2_pe_level"], expected.loc["32", "figure2_pe_level"])
            self.assertLess(expected.loc["22", "figure2_fitc_level"], expected.loc["32", "figure2_fitc_level"])
            self.assertEqual(expected.loc["22", "figure2_fitc_level"], expected.loc["23", "figure2_fitc_level"])
            self.assertLess(expected.loc["22", "figure2_pe_level"], expected.loc["23", "figure2_pe_level"])
            self.assertEqual(expected.loc["22", "marker"], "ErbB-2")
            self.assertEqual(expected.loc["32", "marker"], "CD171")
            self.assertEqual(sum(reference.role == "isotype_control"), 2)
            self.assertTrue((output / "bead-identity-reference.png").is_file())

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
