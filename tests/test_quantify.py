"""Tests for assay-level arithmetic and failure semantics, without clinical data."""

import unittest

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from facs_neuro.quantify import quantify


def make_inputs():
    markers = [("66", "CD9"), ("67", "CD63"), ("68", "CD81"), ("33", "CD56"), ("97", "REA control")]
    sample_rows, blank_rows = [], []
    for bead_id, marker in markers:
        shared = {
            "batch_id": "batch1", "setup_id": "fortessa1", "bead_id": bead_id,
            "marker": marker, "role": "isotype_control" if bead_id == "97" else "capture_marker",
            "n_events": 150, "measurement_valid": True, "qc_flags": "",
        }
        sample_rows.append({
            **shared, "sample_id": "1", "subject_id": "donor1", "blank_id": "blank1",
            "raw_median_apc": {"CD9": 20, "CD63": 30, "CD81": 40, "CD56": 5, "REA control": 10}[marker],
        })
        blank_rows.append({**shared, "control_id": "blank1", "raw_median_apc": 10})
    return pd.DataFrame(sample_rows), pd.DataFrame(blank_rows)


class QuantificationTests(unittest.TestCase):
    def test_signed_correction_and_detection_are_not_binary_positivity(self):
        samples, blanks = make_inputs()
        result = quantify(samples, blanks).set_index("marker")
        self.assertEqual(result.loc["CD56", "corrected_apc"], -5)
        self.assertEqual(result.loc["CD56", "detection_state"], "not_detected")
        self.assertEqual(result.loc["REA control", "detection_state"], "at_background")
        self.assertEqual(result.loc["CD9", "detection_state"], "unclassified_positive")
        self.assertTrue(result["normalized_apc"].isna().all())
        self.assertTrue(result["normalization_status"].eq("not_requested").all())

    def test_normalization_arithmetic_and_no_input_mutation(self):
        samples, blanks = make_inputs()
        original_samples, original_blanks = samples.copy(deep=True), blanks.copy(deep=True)
        result = quantify(samples, blanks, normalize=True).set_index("marker")
        self.assertTrue(result["normalization_factor"].eq(20).all())
        self.assertAlmostEqual(result.loc["CD9", "normalized_apc"], 0.5)
        self.assertAlmostEqual(result.loc["CD56", "normalized_apc"], -0.25)
        self.assertEqual(result.loc["REA control", "role"], "isotype_control")
        assert_frame_equal(samples, original_samples)
        assert_frame_equal(blanks, original_blanks)

    def test_missing_blank_is_missing_not_zero(self):
        samples, _ = make_inputs()
        result = quantify(samples, pd.DataFrame(), normalize=True)
        self.assertTrue(result["corrected_apc"].isna().all())
        self.assertTrue(result["detection_state"].eq("unavailable").all())
        self.assertTrue(result["qc_flags"].str.contains("missing_assay_blank").all())

    def test_blanks_cannot_cross_batch_or_setup(self):
        for field, status in [("batch_id", "blank_batch_mismatch"), ("setup_id", "blank_setup_mismatch")]:
            with self.subTest(field=field):
                samples, blanks = make_inputs()
                blanks[field] = "different"
                result = quantify(samples, blanks)
                self.assertTrue(result["corrected_apc"].isna().all())
                self.assertTrue(result["correction_status"].eq(status).all())

    def test_wrong_panel_identity_is_rejected(self):
        samples, blanks = make_inputs()
        blanks.loc[0, "marker"] = "different_marker"
        result = quantify(samples, blanks)
        self.assertEqual(result.loc[0, "correction_status"], "blank_marker_mismatch")
        self.assertTrue(np.isnan(result.loc[0, "corrected_apc"]))

    def test_flagged_low_counts_do_not_produce_numbers(self):
        samples, blanks = make_inputs()
        samples.loc[3, ["measurement_valid", "qc_flags"]] = [False, "low_bead_count"]
        result = quantify(samples, blanks, normalize=True)
        self.assertTrue(np.isnan(result.loc[3, "corrected_apc"]))
        self.assertTrue(np.isnan(result.loc[3, "normalized_apc"]))
        self.assertIn("low_bead_count", result.loc[3, "qc_flags"])
        self.assertEqual(result.loc[3, "normalization_status"], "invalid_measurement")
        self.assertEqual(result.loc[3, "normalization_factor"], 20)

    def test_zero_counts_and_invalid_blank_fail(self):
        samples, blanks = make_inputs()
        samples.loc[3, "n_events"] = 0
        blanks.loc[4, "measurement_valid"] = False
        result = quantify(samples, blanks)
        self.assertTrue(np.isnan(result.loc[3, "corrected_apc"]))
        self.assertEqual(result.loc[4, "correction_status"], "invalid_blank_measurement")

    def test_missing_or_nonpositive_normalizer_never_triggers_fallback(self):
        for invalid in ["missing", "nonpositive", "invalid"]:
            with self.subTest(invalid=invalid):
                samples, blanks = make_inputs()
                expected = {"missing": "missing_normalizer", "nonpositive": "nonpositive_normalizer", "invalid": "invalid_normalizer"}[invalid]
                if invalid == "missing":
                    samples = samples[samples["marker"] != "CD63"]
                elif invalid == "nonpositive":
                    samples.loc[1, "raw_median_apc"] = 9
                else:
                    samples.loc[1, "measurement_valid"] = False
                result = quantify(samples, blanks, normalize=True)
                self.assertTrue(result["normalization_status"].eq(expected).all())
                self.assertTrue(result["normalized_apc"].isna().all())

    def test_predeclared_subset_is_fixed_and_floor_is_strict(self):
        samples, blanks = make_inputs()
        samples.loc[1, "measurement_valid"] = False
        result = quantify(samples, blanks, normalize=True, normalizers=("CD9", "CD81"))
        self.assertEqual(result.loc[0, "normalization_factor"], 20)
        result = quantify(samples, blanks, normalize=True, normalizers=("CD9", "CD81"), min_normalization_factor=20)
        self.assertTrue(result["normalization_status"].eq("factor_below_minimum").all())

    def test_duplicate_rows_and_inconsistent_acquisitions_raise(self):
        samples, blanks = make_inputs()
        with self.assertRaisesRegex(ValueError, "duplicate sample_id/bead_id"):
            quantify(pd.concat([samples, samples.iloc[[0]]]), blanks)
        with self.assertRaisesRegex(ValueError, "duplicate control_id/bead_id"):
            quantify(samples, pd.concat([blanks, blanks.iloc[[0]]]))
        samples.loc[0, "blank_id"] = "another"
        with self.assertRaisesRegex(ValueError, "inconsistent blank_id"):
            quantify(samples, blanks)

    def test_explicit_acquisition_ids_keep_donors_independent(self):
        samples, blanks = make_inputs()
        second = samples.copy()
        second["sample_id"] = "2"
        second["subject_id"] = "donor2"
        second["raw_median_apc"] *= 2
        result = quantify(pd.concat([samples, second], ignore_index=True), blanks, normalize=True)
        self.assertEqual(result.loc[0, "normalization_factor"], 20)
        self.assertEqual(result.loc[5, "normalization_factor"], 50)

    def test_shuffled_rows_and_duplicate_dataframe_indexes_do_not_change_identity(self):
        samples, blanks = make_inputs()
        samples = samples.iloc[[3, 1, 4, 0, 2]].copy()
        samples.index = [0] * len(samples)
        result = quantify(samples, blanks.iloc[::-1], normalize=True)
        self.assertEqual(result["marker"].tolist(), samples["marker"].tolist())
        self.assertAlmostEqual(result.iloc[0]["normalized_apc"], -0.25)

    def test_technical_acquisitions_share_sample_but_not_normalization(self):
        samples, blanks = make_inputs()
        samples["acquisition_id"] = "sample1_run1"
        second = samples.copy()
        second["acquisition_id"] = "sample1_run2"
        second["raw_median_apc"] *= 2
        result = quantify(pd.concat([samples, second], ignore_index=True), blanks, normalize=True)
        self.assertTrue(result["sample_id"].eq("1").all())
        self.assertTrue(result["subject_id"].eq("donor1").all())
        self.assertEqual(result.loc[0, "normalization_factor"], 20)
        self.assertEqual(result.loc[5, "normalization_factor"], 50)
        with self.assertRaisesRegex(ValueError, "duplicate acquisition_id/bead_id"):
            quantify(pd.concat([samples, samples.iloc[[0]]]), blanks)

    def test_strings_cannot_masquerade_as_boolean_validity(self):
        samples, blanks = make_inputs()
        samples["measurement_valid"] = "False"
        with self.assertRaisesRegex(ValueError, "must contain booleans"):
            quantify(samples, blanks)

    def test_configuration_validation(self):
        samples, blanks = make_inputs()
        for normalizers in [(), ("CD9", "CD9"), "CD9"]:
            with self.subTest(normalizers=normalizers), self.assertRaises(ValueError):
                quantify(samples, blanks, normalize=True, normalizers=normalizers)
        for floor in [-1, float("nan"), float("inf"), True, False, np.bool_(True)]:
            with self.subTest(floor=floor), self.assertRaises(ValueError):
                quantify(samples, blanks, min_normalization_factor=floor)

    def test_a_biological_sample_cannot_change_donor_across_acquisitions(self):
        samples, blanks = make_inputs()
        samples["acquisition_id"] = "first"
        second = samples.copy()
        second["acquisition_id"] = "second"
        second["subject_id"] = "different_donor"
        with self.assertRaisesRegex(ValueError, "multiple subject_id"):
            quantify(pd.concat([samples, second]), blanks)

    def test_isotype_or_ambiguous_marker_cannot_be_a_normalizer(self):
        samples, blanks = make_inputs()
        result = quantify(samples, blanks, normalize=True, normalizers=("REA control",))
        self.assertTrue(result["normalization_status"].eq("invalid_normalizer").all())
        samples.loc[3, "marker"] = "CD9"
        blanks.loc[3, "marker"] = "CD9"
        result = quantify(samples, blanks, normalize=True)
        self.assertTrue(result["normalization_status"].eq("ambiguous_normalizer").all())


if __name__ == "__main__":
    unittest.main()
