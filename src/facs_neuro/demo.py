"""Generate explicitly synthetic FCS3.1 files for an end-to-end software demo."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .pipeline import analyze, new_output, save_json


CHANNELS = ("FITC-A", "PE-A", "APC-A", "FSC-A", "SSC-A", "FSC-H", "Time")
BEADS = "22 23 24 32 33 34 35 42 43 44 45 46 52 53 54 55 56 57 63 64 65 66 67 68 74 75 76 77 78 79 85 86 87 88 89 96 97 98 99".split()
NAMES = "ErbB-2,CD49a,CD140a,CD171,CD56,CD13,CD31,CD222,CD11b,CD24,CD45,CD133/1,O4,GLAST,CD90,CX3CR1,EGFR,Podoplanin,A2B5,CD29,CSPG4,CD9,CD63,CD81,CD47,PSA-NCAM,CD44,CD38,CD49f,CD106,CD49e,CD54,CD36,CD107a,CD68,VGLUT2,REA control,Ganglioside GD2,mIgG1 control".split(",")


def write_fcs(path: Path, events: np.ndarray, *, label: str) -> None:
    """Tiny synthetic IEEE float writer, not a general FCS exporter."""
    payload = np.asarray(events, dtype="<f4").tobytes()
    meta = {"$BEGINANALYSIS": "0", "$ENDANALYSIS": "0", "$BEGINSTEXT": "0", "$ENDSTEXT": "0",
            "$BYTEORD": "1,2,3,4", "$DATATYPE": "F", "$MODE": "L", "$NEXTDATA": "0",
            "$PAR": str(len(CHANNELS)), "$TOT": str(len(events)), "$FIL": label,
            "$BEGINDATA": "0000000000", "$ENDDATA": "0000000000", "$SRC": "SYNTHETIC SOFTWARE TEST"}
    for i, name in enumerate(CHANNELS, 1):
        meta.update({f"$P{i}N": name, f"$P{i}B": "32", f"$P{i}R": "262144", f"$P{i}E": "0,0", f"$P{i}G": "1"})
    def encode():
        return ("|" + "|".join(str(v).replace("|", "||") for pair in meta.items() for v in pair) + "|").encode("ascii")
    text = encode()
    start = 58 + len(text)
    end = start + len(payload) - 1
    meta["$BEGINDATA"] = f"{start:010d}"
    meta["$ENDDATA"] = f"{end:010d}"
    text = encode()
    assert 58 + len(text) == start
    header = b"FCS3.1" + b"    " + b"".join(f"{n:8d}".encode() for n in (58, start - 1, start, end, 0, 0))
    path.write_bytes(header + text + payload)


def _manifest(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def generate_demo(output: str | Path, *, donors: int = 24, run: bool = True) -> Path:
    if isinstance(donors, bool) or not isinstance(donors, int) or not 1 <= donors <= 24:
        raise ValueError("Demo donors must be an integer from 1 to 24")
    root = new_output(output)
    (root / "fcs").mkdir()
    rng = np.random.default_rng(20261029)
    config = {
        "schema_version": 1, "dataset_kind": "synthetic", "samples": "samples.csv",
        "controls": "controls.csv", "panel": "panel.csv",
        "normalization": {"enabled": True, "markers": ["CD9", "CD63", "CD81"], "minimum_factor": 1},
        "setups": {"synthetic_setup": {
            "reviewed": True, "gate_version": "SYNTHETIC-v1", "review_note": "Known generated centers; SOFTWARE TEST ONLY",
            "channels": {"fitc": "FITC-A", "pe": "PE-A", "apc": "APC-A"},
            "compensation": {"input_state": "uncompensated", "mode": "none", "reason": "Synthetic independent channels, no spillover"},
            "bead_space": {"x_transform": {"kind": "linear"}, "y_transform": {"kind": "linear"}},
            "pre_gates": [
                {"name": "beads", "x": "FSC-A", "y": "SSC-A", "vertices": [[50, 20], [150, 20], [150, 100], [50, 100]]},
                {"name": "singlets", "x": "FSC-A", "y": "FSC-H", "vertices": [[50, 40], [150, 140], [150, 160], [50, 60]]}],
            "min_events_per_population": 50,
            "populations": [],
        }},
    }
    panel = []
    for bead, name in zip(BEADS, NAMES, strict=True):
        x, y = int(bead[0]) * 1000, int(bead[1]) * 1000
        config["setups"]["synthetic_setup"]["populations"].append(
            {"bead_id": bead, "vertices": [[x - 250, y - 250], [x + 250, y - 250], [x + 250, y + 250], [x - 250, y + 250]]})
        panel.append({"bead_id": bead, "marker": name, "role": "isotype_control" if bead in {"97", "99"} else "capture_marker"})
    _manifest(root / "panel.csv", panel)
    shared = {"batch_id": "SYNTHETIC", "setup_id": "synthetic_setup", "instrument_settings_id": "synthetic_settings",
              "kit_lot": "SYNTHETIC", "detection_reagent": "CD9+CD63+CD81", "incubation_protocol": "synthetic"}
    controls = [
        {"control_id": "setup", "control_type": "setup_beads", "fcs_path": "fcs/setup.fcs", **shared},
        {"control_id": "blank", "control_type": "assay_blank", "fcs_path": "fcs/blank.fcs", **shared}]
    samples = [{"sample_id": str(i), "acquisition_id": str(i), "subject_id": str(i), "fcs_path": f"fcs/{i}.fcs",
                "blank_id": "blank", "setup_control_id": "setup", "source_material": "SYNTHETIC", "enrichment_method": "SYNTHETIC", **shared}
               for i in range(1, donors + 1)]
    for row in controls + samples:
        blocks = []
        donor = int(row["sample_id"]) if "sample_id" in row else 0
        for idx, bead in enumerate(BEADS):
            n = 80
            events = np.empty((n, len(CHANNELS)), dtype=float)
            events[:, 0] = int(bead[0]) * 1000 + rng.normal(0, 20, n)
            events[:, 1] = int(bead[1]) * 1000 + rng.normal(0, 20, n)
            signal = (idx + 1) * 15 * (1 + donor / 24) if donor else 0
            if bead in {"97", "99"}:
                signal = 2 if donor else 0
            if donor == 1 and bead == "22":
                signal = -10  # A negative corrected result must survive processing.
            events[:, 2] = 100 + (idx % 4) + signal + rng.normal(0, .5, n)
            events[:, 3] = rng.normal(100, 5, n)
            events[:, 4] = rng.normal(60, 5, n)
            events[:, 5] = events[:, 3] + rng.normal(0, 1, n)
            events[:, 6] = np.arange(n) + idx * n
            blocks.append(events)
        values = np.vstack(blocks)
        # Shuffle events so identity cannot depend on population order.
        rng.shuffle(values, axis=0)
        write_fcs(root / row["fcs_path"], values, label="SYNTHETIC-" + row["fcs_path"].split("/")[-1])
    _manifest(root / "samples.csv", samples)
    _manifest(root / "controls.csv", controls)
    save_json(root / "analysis.json", config)
    (root / "SYNTHETIC_DATA.txt").write_text("All FCS values and donor IDs here are synthetic software-test data, not experimental results.\n")
    if run:
        analyze(root / "analysis.json", root / "results")
    return root
