"""Fail-closed configuration and manifest checks run before event processing."""

from copy import deepcopy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from facs_neuro.config import (
    ConfigurationError, CONTROL_TYPES, EXPECTED_MARKERS,
    load_config, load_inputs, readiness_errors,
)


class ConfigurationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.base = {
            "schema_version": 1,
            "dataset_kind": "synthetic",
            "samples": "samples.csv",
            "controls": "controls.csv",
            "panel": "panel.csv",
            "normalization": {
                "enabled": False, "minimum_factor": None,
                "markers": ["CD9", "CD63", "CD81"],
            },
            "setups": {"fortessa": {
                "channels": {"fitc": "FITC-A", "pe": "PE-A", "apc": "APC-A"},
            }},
        }

    def write_config(self, value=None):
        path = self.root / "analysis.json"
        path.write_text(json.dumps(self.base if value is None else value))
        return path

    def readiness_fixture(self):
        shared = {
            "batch_id": "batch1", "setup_id": "fortessa",
            "instrument_settings_id": "settings1", "kit_lot": "lot1",
            "detection_reagent": "CD9+CD63+CD81", "incubation_protocol": "short",
        }

        def file_fields(name):
            path = self.root / f"{name}.fcs"
            # Readiness checks file existence, not FCS event content.
            path.write_bytes(b"not parsed by metadata preflight")
            return {"fcs_path": path.name, "_file": str(path)}

        samples = [{
            **shared, **file_fields("sample1"), "sample_id": "1",
            "subject_id": "donor1", "acquisition_id": "1-a",
            "blank_id": "blank1", "setup_control_id": "setup1",
        }]
        controls = [
            {**shared, **file_fields("blank1"), "control_id": "blank1", "control_type": "assay_blank"},
            {**shared, **file_fields("setup1"), "control_id": "setup1", "control_type": "setup_beads"},
        ]
        return samples, controls

    def write_manifests(self, samples, controls):
        for filename, rows in (("samples.csv", samples), ("controls.csv", controls)):
            fields = list(dict.fromkeys(key for row in rows for key in row if not key.startswith("_")))
            with (self.root / filename).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows({key: row.get(key, "") for key in fields} for row in rows)
        with (self.root / "panel.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["bead_id", "marker", "role"])
            writer.writeheader()
            writer.writerows({
                "bead_id": bead, "marker": marker,
                "role": "isotype_control" if bead in {"97", "99"} else "capture_marker",
            } for bead, marker in EXPECTED_MARKERS.items())

    def test_disabled_normalization_accepts_none_without_guessing_a_floor(self):
        result = load_config(self.write_config())
        self.assertFalse(result["normalization"]["enabled"])
        self.assertIsNone(result["normalization"]["minimum_factor"])
        self.assertEqual(result["normalization"]["markers"], ["CD9", "CD63", "CD81"])
        omitted = deepcopy(self.base)
        del omitted["normalization"]
        self.assertNotIn("normalization", load_config(self.write_config(omitted)))

    def test_numeric_nonnegative_floors_are_accepted_in_both_modes(self):
        for enabled in (False, True):
            for floor in (0, 0.0, 2, 1.25):
                with self.subTest(enabled=enabled, floor=floor):
                    config = deepcopy(self.base)
                    config["normalization"].update(enabled=enabled, minimum_factor=floor)
                    result = load_config(self.write_config(config))
                    self.assertEqual(result["normalization"]["minimum_factor"], floor)

    def test_enabled_normalization_requires_explicit_nonnull_floor(self):
        for missing in (False, True):
            with self.subTest(missing=missing):
                config = deepcopy(self.base)
                config["normalization"]["enabled"] = True
                if missing:
                    del config["normalization"]["minimum_factor"]
                with self.assertRaisesRegex(ConfigurationError, "minimum_factor is required"):
                    load_config(self.write_config(config))

    def test_supplied_bad_floor_is_rejected_even_when_disabled(self):
        invalid = (-1, -0.1, True, False, "0", "", [], {}, 10 ** 400)
        for enabled in (False, True):
            for floor in invalid:
                with self.subTest(enabled=enabled, floor=str(floor)[:30]):
                    config = deepcopy(self.base)
                    config["normalization"].update(enabled=enabled, minimum_factor=floor)
                    with self.assertRaisesRegex(ConfigurationError, "minimum_factor.*finite nonnegative"):
                        load_config(self.write_config(config))

    def test_normalization_enabled_requires_a_real_boolean(self):
        for enabled in (None, 0, 1, "false", "true", [], {}):
            with self.subTest(enabled=enabled):
                config = deepcopy(self.base)
                config["normalization"]["enabled"] = enabled
                with self.assertRaisesRegex(ConfigurationError, "enabled must be true or false"):
                    load_config(self.write_config(config))

    def test_marker_configuration_is_validated_without_activating_normalization(self):
        for markers in (None, [], "CD9", ["CD9", "CD9"], ["CD56"], [["CD9"]], [1]):
            with self.subTest(markers=markers):
                config = deepcopy(self.base)
                config["normalization"]["markers"] = markers
                with self.assertRaisesRegex(ConfigurationError, "normalization.markers"):
                    load_config(self.write_config(config))
        config = deepcopy(self.base)
        config["normalization"]["markers"] = ["CD9", "CD81"]
        result = load_config(self.write_config(config))
        self.assertFalse(result["normalization"]["enabled"])
        self.assertEqual(result["normalization"]["markers"], ["CD9", "CD81"])

    def test_json_nonfinite_constants_have_explicit_configuration_errors(self):
        for constant in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(constant=constant):
                path = self.write_config()
                payload = path.read_text().replace('"minimum_factor": null', f'"minimum_factor": {constant}')
                path.write_text(payload)
                with self.assertRaisesRegex(ConfigurationError, f"Non-finite JSON constant.*{constant}"):
                    load_config(path)
        # JSON may encode a syntactically valid number that overflows a float.
        path = self.write_config()
        path.write_text(path.read_text().replace('"minimum_factor": null', '"minimum_factor": 1e999'))
        with self.assertRaisesRegex(ConfigurationError, "minimum_factor.*finite nonnegative"):
            load_config(path)

    def test_real_template_retains_unknown_acquisition_details_and_disabled_none(self):
        path = Path(__file__).resolve().parents[1] / "config" / "analysis.json"
        original_text = path.read_text()
        original = json.loads(original_text)
        result = load_config(path)
        self.assertFalse(result["normalization"]["enabled"])
        self.assertIsNone(result["normalization"]["minimum_factor"])
        self.assertEqual(result["setups"], original["setups"])
        self.assertEqual(path.read_text(), original_text)

    def test_known_control_types_are_all_accepted_when_unused(self):
        expected = {
            "setup_beads", "assay_blank", "reference_sample", "process_blank",
            "compensation_control", "unstained_control", "instrument_qc",
        }
        self.assertEqual(CONTROL_TYPES, expected)
        for control_type in expected:
            with self.subTest(control_type=control_type):
                samples, controls = self.readiness_fixture()
                extra = deepcopy(controls[0])
                path = self.root / "unused-control.fcs"
                path.write_bytes(b"unused control")
                extra.update(control_id="unused", control_type=control_type, _file=str(path), fcs_path=path.name)
                self.assertEqual(readiness_errors(self.base, samples, [*controls, extra]), [])

    def test_unknown_or_mistyped_unused_control_is_an_explicit_error(self):
        for control_type in ("assay blank", "setup_bead", "ASSAY_BLANK", "", "custom"):
            with self.subTest(control_type=control_type):
                samples, controls = self.readiness_fixture()
                extra = deepcopy(controls[0])
                path = self.root / "unused-control.fcs"
                path.write_bytes(b"unused control")
                extra.update(control_id="unused", control_type=control_type, _file=str(path), fcs_path=path.name)
                errors = readiness_errors(self.base, samples, [*controls, extra])
                self.assertEqual(len(errors), 1)
                self.assertIn("control unused: unknown control_type", errors[0])
                self.assertIn(repr(control_type), errors[0])

    def test_same_sample_mapped_to_different_donors_fails_metadata_preflight(self):
        samples, controls = self.readiness_fixture()
        technical = deepcopy(samples[0])
        path = self.root / "sample1-technical.fcs"
        path.write_bytes(b"not parsed by metadata preflight")
        technical.update(acquisition_id="1-b", subject_id="donor2", _file=str(path), fcs_path=path.name)
        original = deepcopy([*samples, technical])
        errors = readiness_errors(self.base, [*samples, technical], controls)
        self.assertEqual(len(errors), 1)
        self.assertIn("sample '1': sample_id maps to multiple subject_id values", errors[0])
        self.assertIn("donor1, donor2", errors[0])
        self.assertEqual([*samples, technical], original)

    def test_technical_acquisitions_of_same_donor_remain_valid(self):
        samples, controls = self.readiness_fixture()
        technical = deepcopy(samples[0])
        path = self.root / "sample1-technical.fcs"
        path.write_bytes(b"not parsed by metadata preflight")
        technical.update(acquisition_id="1-b", _file=str(path), fcs_path=path.name)
        self.assertEqual(readiness_errors(self.base, [*samples, technical], controls), [])

    def test_placeholder_donor_and_required_metadata_are_still_rejected(self):
        for value in ("unknown", "pending", "not_applicable", " UNKNOWN ", ""):
            for field in ("subject_id", "kit_lot"):
                with self.subTest(value=value, field=field):
                    samples, controls = self.readiness_fixture()
                    samples[0][field] = value
                    original = deepcopy(samples)
                    errors = readiness_errors(self.base, samples, controls)
                    message = "subject_id required" if field == "subject_id" else "kit_lot must be documented"
                    self.assertTrue(any(message in error for error in errors))
                    self.assertEqual(samples, original)

    def test_annotation_columns_remain_reserved_manifest_outputs(self):
        for label in ("sample", "control"):
            for column in ("biological_group", "association_tags", "annotation_source_ids", "annotation_version"):
                with self.subTest(label=label, column=column):
                    samples, controls = self.readiness_fixture()
                    (samples if label == "sample" else controls)[0][column] = "forged annotation"
                    self.write_manifests(samples, controls)
                    config = load_config(self.write_config())
                    with self.assertRaisesRegex(ConfigurationError, "reserved columns"):
                        load_inputs(config)


if __name__ == "__main__":
    unittest.main()
