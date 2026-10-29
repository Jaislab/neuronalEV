import copy
import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd

from facs_neuro.gating import GateError, gate_events, polygon_mask, prepare_events, summarize, transform


def rectangle(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]


def fixture():
    # Deliberately noncanonical channel order, values include signed APC.
    channels = ("APC-A", "PE-A", "FSC-H", "FITC-A", "SSC-A", "FSC-A")
    events = np.array([[-3, 10, 100, 10, 100, 100], [-1, 10, 100, 10, 100, 100], [9, 30, 100, 30, 100, 100], [11, 30, 100, 30, 100, 100], [500, 10, 900, 10, 100, 100]], dtype=float)
    data = SimpleNamespace(events=events, channels=channels, metadata={"datatype": "F"}, ranges={"APC-A": 1000.0})
    setup = {"reviewed": True, "gate_version": "synthetic-v1", "review_note": "Synthetic tests only", "channels": {"fitc": "FITC-A", "pe": "PE-A", "apc": "APC-A"}, "compensation": {"mode": "none", "input_state": "uncompensated", "reason": "Synthetic independent channels"}, "pre_gates": [{"name": "beads", "x": "FSC-A", "y": "SSC-A", "vertices": rectangle(50, 50, 150, 150)}, {"name": "singlets", "x": "FSC-A", "y": "FSC-H", "vertices": rectangle(50, 50, 150, 150)}], "bead_space": {}, "populations": [{"bead_id": "99", "vertices": rectangle(25, 25, 35, 35)}, {"bead_id": "22", "vertices": rectangle(5, 5, 15, 15)}], "min_events_per_population": 2}
    panel = pd.DataFrame([{"bead_id": "22", "marker": "ErbB-2", "role": "capture_marker"}, {"bead_id": "99", "marker": "mIgG1 control", "role": "isotype_control"}])
    return data, setup, panel


class GatingTests(unittest.TestCase):
    def test_shuffled_channels_population_ids_signed_medians_and_singlet_exclusion(self):
        data, setup, panel = fixture()
        result, qc = summarize(data, setup, panel)
        self.assertEqual(result.bead_id.tolist(), ["22", "99"])
        self.assertEqual(result.raw_median_apc.tolist(), [-2.0, 10.0])
        self.assertTrue(result.measurement_valid.all())
        self.assertEqual(qc["total_events"], 5)
        self.assertEqual(qc["pregated_events"], 4)
        self.assertEqual(qc["assigned_events"], 4)

    def test_missing_population_is_nan_and_low_count_is_invalid(self):
        data, setup, panel = fixture()
        data.events = data.events[:1]
        result, _ = summarize(data, setup, panel)
        self.assertEqual(result.loc[0, "raw_median_apc"], -3)
        self.assertFalse(result.loc[0, "measurement_valid"])
        self.assertEqual(result.loc[0, "qc_flags"], "low_bead_count")
        self.assertTrue(np.isnan(result.loc[1, "raw_median_apc"]))
        self.assertEqual(result.loc[1, "qc_flags"], "missing_population")

    def test_overlapping_events_never_assigned_first_match(self):
        data, setup, panel = fixture()
        setup["populations"][0]["vertices"] = rectangle(0, 0, 40, 40)
        result, qc = summarize(data, setup, panel)
        self.assertEqual(qc["ambiguous_events"], 2)
        self.assertEqual(qc["assigned_events"], 2)
        self.assertEqual(result.loc[0, "n_events"], 0)
        self.assertFalse(result.measurement_valid.any())

    def test_compensation_matrix_convention_channel_order_and_no_mutation(self):
        data, setup, _ = fixture()
        truth = np.array([[10., 20., -2.], [8., 1., 3.]])
        matrix = np.array([[1., .2, .05], [.1, 1., .1], [.03, .05, 1.]])
        matrix_channels = ["FITC-A", "PE-A", "APC-A"]
        cols = [data.channels.index(c) for c in matrix_channels]
        data.events = data.events[:2].copy()
        data.events[:, cols] = truth @ matrix
        before = data.events.copy()
        setup["compensation"] = {"mode": "explicit", "input_state": "uncompensated", "channels": matrix_channels, "matrix": matrix.tolist()}
        compensated = prepare_events(data, setup)
        np.testing.assert_allclose(compensated[:, cols], truth, atol=1e-12)
        np.testing.assert_array_equal(data.events, before)
        setup["compensation"]["input_state"] = "compensated"
        with self.assertRaisesRegex(GateError, "second time"):
            prepare_events(data, setup)

    def test_embedded_matrix_missing_matrix_singular_and_unknown_policy(self):
        data, setup, _ = fixture()
        data.metadata["spillover"] = "2,FITC-A,PE-A,1,0.2,0,1"
        setup["compensation"] = {"mode": "embedded", "input_state": "uncompensated"}
        actual = prepare_events(data, setup)
        np.testing.assert_allclose(actual[:, data.channels.index("PE-A")], data.events[:, 1] - .2 * data.events[:, 3])
        del data.metadata["spillover"]
        with self.assertRaisesRegex(GateError, "no SPILLOVER"):
            prepare_events(data, setup)
        setup["compensation"] = {"mode": "explicit", "input_state": "uncompensated", "channels": ["FITC-A", "PE-A"], "matrix": [[1, 1], [1, 1]]}
        with self.assertRaisesRegex(GateError, "singular"):
            prepare_events(data, setup)
        setup["compensation"] = {"mode": "auto", "input_state": "uncompensated"}
        with self.assertRaises(GateError):
            prepare_events(data, setup)

    def test_review_count_threshold_panel_ids_and_apc_gates_are_required(self):
        data, setup, panel = fixture()
        variants = []
        changed = copy.deepcopy(setup); changed["reviewed"] = False; variants.append(changed)
        changed = copy.deepcopy(setup); changed.pop("min_events_per_population"); variants.append(changed)
        changed = copy.deepcopy(setup); changed["pre_gates"][0]["x"] = "APC-A"; variants.append(changed)
        changed = copy.deepcopy(setup); changed["populations"].pop(); variants.append(changed)
        changed = copy.deepcopy(setup); changed["pre_gates"].pop(); variants.append(changed)
        for variant in variants:
            with self.subTest(variant=variant), self.assertRaises(GateError):
                summarize(data, variant, panel)

    def test_saturation_invalidates_population_without_erasing_raw_value(self):
        data, setup, panel = fixture()
        data.events[0, 0] = 1000
        result, qc = summarize(data, setup, panel)
        self.assertFalse(result.loc[0, "measurement_valid"])
        self.assertIn("apc_saturation", result.loc[0, "qc_flags"])
        self.assertEqual(result.loc[0, "raw_median_apc"], 499.5)
        self.assertEqual(qc["apc_saturated_events"], 1)

    def test_float_pnr_minus_one_and_explicit_saturation_threshold(self):
        data, setup, panel = fixture()
        data.events[0, 0] = 999
        result, qc = summarize(data, setup, panel)
        self.assertIn("apc_saturation", result.loc[0, "qc_flags"])
        self.assertEqual(qc["apc_saturation_threshold"], 999)
        setup["apc_saturation_threshold"] = 1001
        result, qc = summarize(data, setup, panel)
        self.assertTrue(result.measurement_valid.all())
        self.assertEqual(qc["apc_saturated_events"], 0)
        for threshold in [0, -1, True, float("inf"), "1000"]:
            with self.subTest(threshold=threshold), self.assertRaises(GateError):
                summarize(data, {**setup, "apc_saturation_threshold": threshold}, panel)

    def test_saturation_is_checked_before_computational_compensation(self):
        data, setup, panel = fixture()
        data.events[0, 0] = 999
        data.events[0, 3] = 10000
        # FITC spillover into APC makes corrected APC safely below the ADC
        # boundary despite its original saturation; original clipping must win.
        setup["compensation"] = {"mode": "explicit", "input_state": "uncompensated", "channels": ["FITC-A", "APC-A"], "matrix": [[1, .1], [0, 1]]}
        setup["populations"][1]["vertices"] = rectangle(0, 0, 11000, 15)
        result, qc = summarize(data, setup, panel)
        self.assertIn("apc_saturation", result.loc[0, "qc_flags"])
        self.assertEqual(qc["apc_saturated_events"], 1)

    def test_compensated_exports_require_original_saturation_review(self):
        data, setup, panel = fixture()
        setup["compensation"] = {"mode": "none", "input_state": "compensated", "reason": "Export processing documented"}
        result, qc = summarize(data, setup, panel)
        self.assertFalse(result.measurement_valid.any())
        self.assertTrue(result.qc_flags.str.contains("acquisition_saturation_unverified").all())
        self.assertIsNone(qc["apc_saturated_events"])
        self.assertEqual(qc["apc_saturation_source"], "unverified")
        setup["saturation_review_note"] = "Original acquisition review found no clipped APC events in these gates"
        result, qc = summarize(data, setup, panel)
        self.assertTrue(result.measurement_valid.all())
        self.assertIsNone(qc["apc_saturated_events"])
        self.assertEqual(qc["apc_saturation_source"], "manual_acquisition_review")

    def test_pre_gates_require_real_scatter_and_pulse_geometry(self):
        data, setup, panel = fixture()
        # A channel alias that avoids the APC string must not bypass endpoint
        # independence by masquerading as an arbitrary detector.
        data.channels = (*data.channels, "Red-670-A")
        data.events = np.column_stack([data.events, data.events[:, 0]])
        changed = copy.deepcopy(setup)
        changed["pre_gates"][0]["y"] = "Red-670-A"
        with self.assertRaisesRegex(GateError, "conventional"):
            summarize(data, changed, panel)
        changed = copy.deepcopy(setup)
        changed["pre_gates"][0]["y"] = "FSC-H"
        with self.assertRaisesRegex(GateError, "FSC versus SSC"):
            summarize(data, changed, panel)
        changed = copy.deepcopy(setup)
        changed["pre_gates"][1]["y"] = "SSC-A"
        with self.assertRaisesRegex(GateError, "area versus height or width"):
            summarize(data, changed, panel)

    def test_polygon_boundary_and_concave_polygon(self):
        points = np.array([[0, 0], [1, 1], [2, 2], [1, 2], [-1, 0]])
        np.testing.assert_array_equal(polygon_mask(points, rectangle(0, 0, 2, 2)), [True, True, True, True, False])
        concave = [[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]]
        np.testing.assert_array_equal(polygon_mask([[.5, 1.5], [1.5, 1.5]], concave), [True, False])
        with self.assertRaises(GateError):
            polygon_mask(points, [[0, 0], [2, 0], [0, 2], [2, 2]])

    def test_asinh_transform_does_not_mutate_or_clip_negatives(self):
        data = np.array([-10., 0., 10.])
        before = data.copy()
        actual = transform(data, {"kind": "asinh", "cofactor": 5})
        np.testing.assert_allclose(actual, np.arcsinh(data / 5))
        np.testing.assert_array_equal(data, before)
        with self.assertRaises(GateError):
            transform(data, {"kind": "asinh", "cofactor": 0})

    def test_review_fields_must_be_nonempty_strings(self):
        data, setup, panel = fixture()
        for field in ("gate_version", "review_note"):
            for value in (None, False, 1, [], {}, "", "   "):
                changed = copy.deepcopy(setup)
                changed[field] = value
                with self.subTest(field=field, value=value), self.assertRaises(GateError):
                    summarize(data, changed, panel)

    def test_compensation_cannot_modify_scatter_or_time(self):
        data, setup, _ = fixture()
        for channel in ("FSC-A", "SSC-A"):
            changed = copy.deepcopy(setup)
            changed["compensation"] = {"mode": "explicit", "input_state": "uncompensated", "channels": [channel], "matrix": [[1]]}
            with self.subTest(channel=channel), self.assertRaisesRegex(GateError, "cannot modify scatter"):
                prepare_events(data, changed)
        data.channels = (*data.channels, "Time")
        data.events = np.column_stack([data.events, np.arange(len(data.events))])
        setup["compensation"] = {"mode": "explicit", "input_state": "uncompensated", "channels": ["Time"], "matrix": [[1]]}
        with self.assertRaisesRegex(GateError, "cannot modify scatter"):
            prepare_events(data, setup)

    def test_stages_record_cumulative_pregate_event_counts_without_api_change(self):
        data, setup, panel = fixture()
        stages = []
        self.assertEqual(len(gate_events(data, setup, panel, stages=stages)), 3)
        self.assertEqual(stages, [{"gate": "beads", "events": 5}, {"gate": "singlets", "events": 4}])
        _, qc = summarize(data, setup, panel)
        self.assertEqual(qc["gate_counts"], stages)

    def test_malformed_json_objects_fail_as_gate_errors(self):
        data, setup, panel = fixture()
        malformed = [
            ("channels", []), ("compensation", []), ("bead_space", []),
            ("pre_gates", ["beads"]), ("populations", ["22"]),
        ]
        for field, value in malformed:
            with self.subTest(field=field), self.assertRaises(GateError):
                summarize(data, {**setup, field: value}, panel)
        changed = copy.deepcopy(setup)
        changed["pre_gates"][0]["name"] = []
        with self.assertRaises(GateError):
            summarize(data, changed, panel)
        changed = copy.deepcopy(setup)
        changed["pre_gates"][0]["x"] = []
        with self.assertRaises(GateError):
            summarize(data, changed, panel)
        for spec in ([], "asinh", 5):
            with self.subTest(spec=spec), self.assertRaises(GateError):
                transform([0, 1], spec)


if __name__ == "__main__":
    unittest.main()
