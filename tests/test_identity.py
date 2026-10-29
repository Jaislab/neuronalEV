"""Figure 2 consistency checks use synthetic data, not assay calibration."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from facs_neuro.config import EXPECTED_MARKERS
from facs_neuro.gating import GateError
from facs_neuro.identity import FIGURE2_LEVELS, audit_identity, write_identity_reference


def rectangle(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]


def fixture(*, mode="normal", omit=None):
    beads = sorted(EXPECTED_MARKERS)
    # Exactly 5000 total, deliberately uneven. Each cell is distinct, and all
    # counts are sufficiently positive without imposing a biological cutoff.
    counts = np.arange(109, 148)
    counts[:8] += 1
    assert counts.sum() == 5000
    channels = ("APC-A", "FSC-H", "FITC-A", "SSC-A", "PE-A", "FSC-A")
    setup = {"reviewed": True, "gate_version": "synthetic-identity", "review_note": "Known software test only", "channels": {"fitc": "FITC-A", "pe": "PE-A", "apc": "APC-A"}, "compensation": {"input_state": "uncompensated", "mode": "none", "reason": "Independent synthetic channels"}, "pre_gates": [{"name": "beads", "x": "FSC-A", "y": "SSC-A", "vertices": rectangle(50, 50, 150, 150)}, {"name": "singlets", "x": "FSC-A", "y": "FSC-H", "vertices": rectangle(50, 50, 150, 150)}], "bead_space": {}, "min_events_per_population": 1, "populations": []}
    panel, blocks = [], []
    for bead, n in zip(beads, counts):
        x, y = FIGURE2_LEVELS[bead]
        if mode == "transpose":
            x, y = y, x
        elif mode == "reverse":
            x, y = 12 - x, 12 - y
        elif mode == "reverse_fitc":
            x = 12 - x
        elif mode == "reverse_pe":
            y = 12 - y
        elif mode == "warp":
            # Strong monotone, nonuniform spacing plus row/column stagger.
            x, y = x ** 2 + .3 * y, y ** 3 + .2 * x
        x, y = x * 100, y * 100
        setup["populations"].append({"bead_id": bead, "vertices": rectangle(x - 10, y - 10, x + 10, y + 10)})
        role = "isotype_control" if bead in {"97", "99"} else "capture_marker"
        panel.append({"bead_id": bead, "marker": EXPECTED_MARKERS[bead], "role": role})
        if bead == omit:
            continue
        values = np.zeros((n, len(channels)), dtype=float)
        values[:, 0] = np.linspace(-100, 100, n)
        values[:, 1] = values[:, 3] = values[:, 5] = 100
        values[:, 2] = x + np.linspace(-1, 1, n)
        values[:, 4] = y + np.linspace(-1, 1, n)
        blocks.append(values)
    events = np.vstack(blocks)
    np.random.default_rng(1234).shuffle(events, axis=0)
    data = SimpleNamespace(events=events, channels=channels, ranges={"APC-A": 262144}, metadata={"datatype": "F"})
    return data, setup, pd.DataFrame(panel)


class IdentityTests(unittest.TestCase):
    def test_full39_uneven_5000_event_sample_preserves_actual_counts(self):
        data, setup, panel = fixture()
        result = audit_identity(data, setup, panel)
        self.assertEqual(len(data.events), 5000)
        self.assertEqual(len(result), 39)
        self.assertEqual(result.n_events.sum(), 5000)
        self.assertGreater(result.n_events.nunique(), 1)
        self.assertAlmostEqual(result.n_events.mean(), 5000 / 39)
        self.assertTrue(result.layout_status.eq("consistent").all())
        self.assertEqual(result.role.value_counts().to_dict(), {"capture_marker": 37, "isotype_control": 2})
        indexed = result.set_index("bead_id")
        self.assertEqual(indexed.loc["66", "marker"], "CD9")
        self.assertEqual(indexed.loc["99", "marker"], "mIgG1 control")
        self.assertEqual(indexed.loc["22", "figure2_fitc_level"], 2)
        self.assertEqual(indexed.loc["22", "figure2_pe_level"], 2)

    def test_swapped_gate_labels_with_canonical_marker_names_are_inconsistent(self):
        data, setup, panel = fixture()
        lower = next(g for g in setup["populations"] if g["bead_id"] == "22")
        upper = next(g for g in setup["populations"] if g["bead_id"] == "32")
        lower["vertices"], upper["vertices"] = upper["vertices"], lower["vertices"]
        result = audit_identity(data, setup, panel).set_index("bead_id")
        self.assertEqual(result.loc["22", "layout_status"], "inconsistent")
        self.assertIn("figure2_fitc_order_vs_32", result.loc["22", "layout_flags"])
        self.assertIn("figure2_fitc_order_vs_22", result.loc["32", "layout_flags"])

    def test_axis_transposition_and_reflection_fail_ordinal_check(self):
        for mode in ("transpose", "reverse", "reverse_fitc", "reverse_pe"):
            with self.subTest(mode=mode):
                result = audit_identity(*fixture(mode=mode))
                self.assertTrue(result.layout_status.eq("inconsistent").any())

    def test_apc_values_do_not_influence_identity_check(self):
        data, setup, panel = fixture()
        original = audit_identity(data, setup, panel)
        data.events[:, 0] = np.linspace(-1e8, 1e8, len(data.events))
        changed = audit_identity(data, setup, panel)
        pd.testing.assert_frame_equal(original, changed)

    def test_monotone_warp_and_stagger_need_no_schematic_spacing(self):
        result = audit_identity(*fixture(mode="warp"))
        self.assertTrue(result.layout_status.eq("consistent").all())

    def test_missing_population_is_not_recovered_or_filled(self):
        result = audit_identity(*fixture(omit="66")).set_index("bead_id")
        self.assertEqual(result.loc["66", "n_events"], 0)
        self.assertTrue(np.isnan(result.loc["66", "observed_fitc_median"]))
        self.assertEqual(result.loc["66", "layout_status"], "missing")
        self.assertEqual(result.loc["66", "layout_flags"], "missing_population")
        self.assertEqual(len(result), 39)

    def test_panel_validation_and_order_independence(self):
        data, setup, panel = fixture()
        with self.assertRaises(GateError):
            audit_identity(data, setup, panel.iloc[:-1])
        bad = panel.copy()
        bad.loc[0, "marker"] = "incorrect label"
        with self.assertRaises(GateError):
            audit_identity(data, setup, bad)
        expected = audit_identity(data, setup, panel).set_index("bead_id").sort_index()
        setup["populations"].reverse()
        actual = audit_identity(data, setup, panel.iloc[::-1]).set_index("bead_id").sort_index()
        pd.testing.assert_frame_equal(expected, actual)

    def test_reference_files_contain_all39_ordinals_and_refuse_overwrite(self):
        _, _, panel = fixture()
        with tempfile.TemporaryDirectory() as tmp:
            write_identity_reference(tmp, panel)
            table = pd.read_csv(Path(tmp) / "bead-identity-reference.csv", dtype={"bead_id": str})
            self.assertEqual(len(table), 39)
            self.assertEqual(table.set_index("bead_id").loc["98", "marker"], "Ganglioside GD2")
            self.assertGreater((Path(tmp) / "bead-identity-reference.png").stat().st_size, 1000)
            with self.assertRaises(FileExistsError):
                write_identity_reference(tmp, panel)


if __name__ == "__main__":
    unittest.main()
