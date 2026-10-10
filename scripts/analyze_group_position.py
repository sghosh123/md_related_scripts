#!/usr/bin/env python3
"""Export center-of-geometry or center-of-mass coordinates per frame."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
from copper_trafficking.distance import require_nonempty_selection  # noqa: E402
from copper_trafficking.trajectory import scaled_time  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--selection", required=True)
    p.add_argument("--method", choices=["cog", "com"], default="cog")
    p.add_argument("--output-csv", required=True)
    p.add_argument("--time-scale", type=float, default=1.0)
    p.add_argument("--time-offset", type=float, default=0.0)
    p.add_argument("--time-unit", default="ps")
    p.add_argument("--assume-group-whole", action="store_true", help="Confirm selection was made whole/unwrapped for correct center coordinates")
    args = p.parse_args()
    if not args.assume_group_whole:
        raise SystemExit("Refusing to calculate a group center without --assume-group-whole: groups split across PBC boundaries can have misleading centers.")
    if not np.isfinite(args.time_scale) or not np.isfinite(args.time_offset) or not args.time_unit.strip():
        raise SystemExit("Time scale/offset must be finite and time unit non-empty.")
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit("MDAnalysis is required. Install with: python -m pip install -e .") from exc
    u = mda.Universe(args.topology, args.trajectory)
    group = u.select_atoms(args.selection)
    try:
        require_nonempty_selection(group, "--selection")
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if args.method == "com":
        try:
            masses = group.masses
        except Exception as exc:
            raise SystemExit("COM requires mass information in the topology.") from exc
        if not np.all(np.isfinite(masses)) or float(np.sum(masses)) <= 0:
            raise SystemExit("COM requires valid, positive atomic masses for every selected atom.")
    out = Path(args.output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["frame", f"time_{args.time_unit}", "x_A", "y_A", "z_A", "method"])
        for ts in u.trajectory:
            center = group.center_of_geometry() if args.method == "cog" else group.center_of_mass()
            writer.writerow([int(ts.frame), scaled_time(ts, args.time_scale, args.time_offset), *[float(x) for x in center], args.method])
    print(f"Wrote group centers to {out}; method={args.method}; selected atoms={len(group)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
