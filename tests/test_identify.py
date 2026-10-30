"""The setup-bead identification script is checked on synthetic beads only."""

from __future__ import annotations

from contextlib import redirect_stdout
import importlib.util
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from facs_neuro.config import EXPECTED_IDS, EXPECTED_MARKERS
from facs_neuro.demo import generate_demo
from facs_neuro.fcs import read_fcs
from facs_neuro.gating import summarize
from facs_neuro.identity import audit_identity

ROOT = Path(__file__).resolve().parents[1]


def load_script():
    spec = importlib.util.spec_from_file_location("identify_setup_beads", ROOT / "scripts" / "identify_setup_beads.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class IdentifyTests(unittest.TestCase):
    def test_synthetic_setup_beads_are_decoded_named_and_left_unreviewed(self):
        script = load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = generate_demo(root / "fixture", donors=1, run=False)
            with redirect_stdout(StringIO()):
                code = script.main([str(fixture / "fcs" / "setup.fcs"), "--output", str(root / "out"),
                                    "--setup-id", "synthetic_setup", "--transform", "linear", "--singlet-axis", "FSC-H"])
            self.assertEqual(code, 0)
            summary = json.loads((root / "out" / "summary.json").read_text())
            self.assertEqual(summary["status"], "identified_unreviewed")
            self.assertEqual(summary["populations_resolved"], 39)
            self.assertEqual(summary["problems"], [])
            table = pd.read_csv(root / "out" / "populations.csv", dtype={"bead_id": str})
            self.assertEqual(set(table.bead_id), EXPECTED_IDS)
            self.assertTrue(all(m == EXPECTED_MARKERS[b] for b, m in zip(table.bead_id, table.marker, strict=True)))
            self.assertEqual(table.set_index("bead_id").loc["66", "marker"], "CD9")
            self.assertTrue((table.n_events >= 70).all())
            setup = json.loads((root / "out" / "setup-template.json").read_text())["synthetic_setup"]
            self.assertFalse(setup["reviewed"])
            self.assertEqual(setup["gate_version"], "")
            self.assertIsNone(setup["min_events_per_population"])
            self.assertEqual({p["bead_id"] for p in setup["populations"]}, EXPECTED_IDS)
            self.assertEqual([g["name"] for g in setup["pre_gates"]], ["beads", "singlets"])
            # The template reproduces identification through the pipeline's own gating once reviewed.
            setup.update(reviewed=True, gate_version="test", review_note="synthetic test", min_events_per_population=50)
            setup["compensation"]["reason"] = "synthetic independent channels"
            panel = pd.read_csv(fixture / "panel.csv", dtype=str, keep_default_na=False)
            data = read_fcs(fixture / "fcs" / "setup.fcs")
            measured, qc = summarize(data, setup, panel)
            self.assertTrue(measured.measurement_valid.all())
            self.assertEqual(qc["ambiguous_events"], 0)
            self.assertTrue((measured.n_events >= 70).all())
            audit = audit_identity(data, setup, panel)
            self.assertTrue(audit.layout_status.eq("consistent").all())
            self.assertTrue((root / "out" / "identified.png").stat().st_size > 1000)


if __name__ == "__main__":
    unittest.main()
