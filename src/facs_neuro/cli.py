"""Command-line interface; incomplete real-data configuration never produces results."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile


def main(argv: list[str] | None = None) -> int:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "facs-neuro-mpl"))
    os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "facs-neuro-cache"))
    parser = argparse.ArgumentParser(description="MACSPlex EV Neuro FCS analysis (assay validation required)")
    sub = parser.add_subparsers(dest="command", required=True)
    inv = sub.add_parser("inspect", help="Inspect manifests and original FCS metadata without gating")
    inv.add_argument("--config", required=True)
    inv.add_argument("--output", help="Optional JSON inventory path; must not exist")
    for command, help_text in (("calibrate", "Create unreviewed setup-bead overviews and empty gate templates"),
                               ("analyze", "Analyze with documented compensation and reviewed gates")):
        cmd = sub.add_parser(command, help=help_text)
        cmd.add_argument("--config", required=True)
        cmd.add_argument("--output", required=True, help="New output directory (never overwritten)")
    rep = sub.add_parser("report", help="Regenerate a self-contained report for a completed run")
    rep.add_argument("--run", required=True)
    demo = sub.add_parser("demo", help="Generate synthetic FCS files and analyze them")
    demo.add_argument("--output", required=True, help="New directory")
    demo.add_argument("--donors", type=int, default=24)
    args = parser.parse_args(argv)
    try:
        from . import pipeline
        if args.command == "inspect":
            result = pipeline.inspect(args.config)
            if args.output:
                with Path(args.output).open("x") as handle:
                    json.dump(result, handle, indent=2, allow_nan=False)
                    handle.write("\n")
            else:
                print(json.dumps(result, indent=2, allow_nan=False))
            return 0 if result["inputs_ready"] else 2
        if args.command == "calibrate":
            print(pipeline.calibrate(args.config, args.output))
        elif args.command == "analyze":
            out = pipeline.analyze(args.config, args.output)
            print(out / "report.html")
            status = json.loads((out / "provenance.json").read_text())["status"]
            return 3 if status == "completed_with_qc_failures" else 0
        elif args.command == "report":
            from .report import render_report
            run = Path(args.run)
            if json.loads((run / "provenance.json").read_text()).get("status") not in {"completed", "completed_with_qc_failures"}:
                raise ValueError("Cannot report an incomplete run")
            print(render_report(run))
        elif args.command == "demo":
            from .demo import generate_demo
            root = generate_demo(args.output, donors=args.donors)
            print(f"SYNTHETIC SOFTWARE DEMO: {root / 'results' / 'report.html'}")
        return 0
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"facs-neuro: {exc}", file=sys.stderr)
        return 2
