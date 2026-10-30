"""Fail-closed configuration and manifest checks run before event processing."""

from copy import deepcopy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from facs_neuro.config import (
    ConfigurationError, CONTROL_TYPES, EXPECTED_MARKERS,
    load_config, load_inputs, readiness_errors, select_samples,
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

    def test_real_config_keeps_study_normalization_and_is_not_mutated(self):
        # The 2026-10-29 study enables tetraspanin normalization with a floor of 500.
        path = Path(__file__).resolve().parents[1] / "config" / "analysis.json"
        original_text = path.read_text()
        original = json.loads(original_text)
        result = load_config(path)
        self.assertIs(result["normalization"]["enabled"], True)
        self.assertEqual(result["normalization"]["minimum_factor"], 500)
        self.assertEqual(result["normalization"]["markers"], ["CD9", "CD63", "CD81"])
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

    def test_background_defaults_remain_full_assay_blank_without_inventing_a_choice(self):
        loaded = load_config(self.write_config())
        self.assertNotIn("background", loaded)
        explicit = deepcopy(self.base)
        explicit["background"] = {"mode": "assay_blank"}
        self.assertEqual(load_config(self.write_config(explicit))["background"], {"mode": "assay_blank"})
        samples, controls = self.readiness_fixture()
        samples[0]["blank_id"] = "setup1"
        errors = readiness_errors(self.base, samples, controls)
        self.assertTrue(any("blank_id must reference assay_blank" in error for error in errors))

    def test_setup_baseline_configuration_requires_explicit_mode_and_rationale(self):
        config = deepcopy(self.base)
        config["background"] = {
            "mode": "setup_baseline",
            "rationale": "User selected setup beads because the stained blank was improperly washed.",
        }
        result = load_config(self.write_config(config))
        self.assertEqual(result["background"], config["background"])
        for background in (
            None, [], "setup_baseline", {"mode": "baseline"}, {"mode": []},
            {"mode": "setup_beads"}, {"mode": ""}, {"mode": True},
            {"mode": "setup_baseline"}, {"mode": "setup_baseline", "rationale": ""},
            {"mode": "setup_baseline", "rationale": "   "},
            {"mode": "setup_baseline", "rationale": None},
            {"mode": "assay_blank", "rationale": 1},
        ):
            with self.subTest(background=background):
                invalid = deepcopy(self.base)
                invalid["background"] = background
                with self.assertRaisesRegex(ConfigurationError, "background"):
                    load_config(self.write_config(invalid))

    def test_one_setup_control_can_be_both_baseline_and_gate_reference_without_staining_metadata(self):
        config = deepcopy(self.base)
        config["background"] = {"mode": "setup_baseline", "rationale": "Selected unstained setup reference."}
        samples, controls = self.readiness_fixture()
        setup = next(control for control in controls if control["control_type"] == "setup_beads")
        setup["detection_reagent"] = ""
        setup["incubation_protocol"] = ""
        samples[0]["blank_id"] = setup["control_id"]
        self.assertEqual(samples[0]["blank_id"], samples[0]["setup_control_id"])
        self.assertEqual(readiness_errors(config, samples, [setup]), [])
        self.write_manifests(samples, [setup])
        loaded = load_config(self.write_config(config))
        actual_samples, actual_controls, _ = load_inputs(loaded)
        self.assertEqual(len(actual_controls), 1)
        self.assertEqual(readiness_errors(loaded, actual_samples, actual_controls), [])

    def test_baseline_still_requires_setup_type_and_matching_acquisition_metadata(self):
        config = deepcopy(self.base)
        config["background"] = {"mode": "setup_baseline", "rationale": "Selected unstained setup reference."}
        samples, controls = self.readiness_fixture()
        errors = readiness_errors(config, samples, controls)
        self.assertTrue(any("blank_id must reference setup_beads" in error for error in errors))
        for field in ("batch_id", "setup_id", "instrument_settings_id", "kit_lot"):
            with self.subTest(field=field):
                samples, controls = self.readiness_fixture()
                setup = next(control for control in controls if control["control_type"] == "setup_beads")
                samples[0]["blank_id"] = setup["control_id"]
                setup[field] = "different"
                errors = readiness_errors(config, samples, [setup])
                self.assertTrue(any(f"{field} differs from referenced setup_beads" in error for error in errors))

    def test_background_output_provenance_columns_are_reserved(self):
        for label in ("sample", "control"):
            for column in ("background_mode", "background_reference_id", "background_limitation"):
                with self.subTest(label=label, column=column):
                    samples, controls = self.readiness_fixture()
                    (samples if label == "sample" else controls)[0][column] = "forged background metadata"
                    self.write_manifests(samples, controls)
                    config = load_config(self.write_config())
                    with self.assertRaisesRegex(ConfigurationError, "reserved columns"):
                        load_inputs(config)

    def test_sample_exclusions_default_preserves_all_samples_without_inventing_config(self):
        config = load_config(self.write_config())
        self.assertNotIn("sample_exclusions", config)
        samples, _ = self.readiness_fixture()
        before = deepcopy(samples)
        active, excluded = select_samples(config, samples)
        self.assertEqual(active, samples)
        self.assertEqual(excluded, [])
        self.assertEqual(samples, before)
        self.assertIsNot(active[0], samples[0])

    def test_sample_exclusions_mapping_requires_nonempty_string_reasons(self):
        for exclusions in (None, [], "10", True, {"": "bead loss"}, {"  ": "bead loss"},
                           {"10": ""}, {"10": "  "}, {"10": None}, {"10": True},
                           {"10": False}, {"10": 1}, {"10": []}, {"10": {}}):
            with self.subTest(exclusions=exclusions):
                config = deepcopy(self.base)
                config["sample_exclusions"] = exclusions
                with self.assertRaisesRegex(ConfigurationError, "sample_exclusions"):
                    load_config(self.write_config(config))
        # JSON itself only admits string object keys. The helper also validates
        # direct Python callers so bool/int keys cannot masquerade as sample IDs.
        samples = [{"sample_id": "10"}, {"sample_id": "11"}]
        for key in (True, False, 10):
            with self.subTest(key=key), self.assertRaisesRegex(ConfigurationError, "keys.*strings"):
                select_samples({"sample_exclusions": {key: "bead loss"}}, samples)

    def test_sample_exclusions_remove_all_acquisitions_and_preserve_exact_reason_and_order(self):
        samples = [
            {"sample_id": "1", "acquisition_id": "1-a", "subject_id": "donor1"},
            {"sample_id": "10", "acquisition_id": "10-a", "subject_id": "donor10"},
            {"sample_id": "2", "acquisition_id": "2-a", "subject_id": "donor2"},
            {"sample_id": "10", "acquisition_id": "10-b", "subject_id": "donor10"},
        ]
        config = deepcopy(self.base)
        reason = "Confirmed bead loss in sample 10; exclude all acquisitions."
        config["sample_exclusions"] = {"10": reason}
        loaded = load_config(self.write_config(config))
        original_samples, original_config = deepcopy(samples), deepcopy(loaded)
        active, excluded = select_samples(loaded, samples)
        self.assertEqual([row["acquisition_id"] for row in active], ["1-a", "2-a"])
        self.assertEqual([row["acquisition_id"] for row in excluded], ["10-a", "10-b"])
        self.assertTrue(all(row["exclusion_reason"] == reason for row in excluded))
        self.assertFalse(any("exclusion_reason" in row for row in active))
        self.assertEqual(samples, original_samples)
        self.assertEqual(loaded, original_config)
        active[0]["subject_id"] = "changed copy"
        excluded[0]["exclusion_reason"] = "changed copy"
        self.assertEqual(samples, original_samples)
        self.assertEqual(loaded, original_config)

    def test_unknown_exclusion_ids_fail_instead_of_silently_skipping_a_typo(self):
        samples = [{"sample_id": "1"}, {"sample_id": "10"}]
        for unknown in ("11", "10.0", " 10"):
            with self.subTest(unknown=unknown), self.assertRaisesRegex(ConfigurationError, "unknown sample_id"):
                select_samples({"sample_exclusions": {unknown: "bead loss"}}, samples)

    def test_exclusions_produce_23_active_donors_from_a_24_sample_manifest(self):
        samples = [{"sample_id": str(i), "acquisition_id": str(i), "subject_id": str(i)} for i in range(1, 25)]
        active, excluded = select_samples({"sample_exclusions": {"10": "Confirmed bead loss"}}, samples)
        self.assertEqual(len(active), 23)
        self.assertEqual(len(excluded), 1)
        self.assertEqual(excluded[0]["sample_id"], "10")
        self.assertNotIn("10", {row["sample_id"] for row in active})

    def test_excluding_every_sample_is_an_error(self):
        samples = [{"sample_id": "10", "acquisition_id": "10-a"},
                   {"sample_id": "10", "acquisition_id": "10-b"}]
        with self.assertRaisesRegex(ConfigurationError, "No active samples remain"):
            select_samples({"sample_exclusions": {"10": "Confirmed bead loss"}}, samples)


if __name__ == "__main__":
    unittest.main()
