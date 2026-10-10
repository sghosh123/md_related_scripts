#!/usr/bin/env python3
"""Calculate fitted RMSD to a selected reference frame."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
from copper_trafficking.trajectory import scaled_time  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--selection", required=True, help="Atom selection used for fit and RMSD")
    p.add_argument("--reference-frame", type=int, default=0)
    p.add_argument("--output-csv", required=True)
    p.add_argument("--output-plot")
    p.add_argument("--time-scale", type=float, default=1.0)
    p.add_argument("--time-offset", type=float, default=0.0)
    p.add_argument("--time-unit", default="ps")
    p.add_argument("--assume-group-whole", action="store_true", help="Confirm molecules/groups are whole and PBC imaging has been handled")
    args = p.parse_args()
    if not args.assume_group_whole:
        raise SystemExit("Refusing to calculate fitted RMSD without --assume-group-whole: PBC-split coordinates can produce misleading fits.")
    if not np.isfinite(args.time_scale) or not np.isfinite(args.time_offset) or not args.time_unit.strip():
        raise SystemExit("Time scale/offset must be finite and time unit non-empty.")
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis.rms import RMSD
    except ImportError as exc:
        raise SystemExit("MDAnalysis is required. Install with: python -m pip install -e .") from exc
    try:
        u = mda.Universe(args.topology, args.trajectory)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Could not load topology/trajectory: {exc}") from exc
    if args.reference_frame < 0 or args.reference_frame >= len(u.trajectory):
        raise SystemExit(f"--reference-frame must be in range 0..{len(u.trajectory)-1}.")
    selection = u.select_atoms(args.selection)
    if len(selection) < 3:
        raise SystemExit("At least three selected atoms are required for a stable fitted RMSD.")
    try:
        analysis = RMSD(u, select=args.selection, ref_frame=args.reference_frame).run()
    except Exception as exc:
        raise SystemExit(f"RMSD calculation failed: {exc}") from exc
    values = analysis.results.rmsd
    out = Path(args.output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["frame", f"time_{args.time_unit}", "fitted_RMSD_A"])
        for frame_value, _time_value, rmsd_value in values:
            frame = int(frame_value)
            u.trajectory[frame]
            time = scaled_time(u.trajectory.ts, args.time_scale, args.time_offset)
            writer.writerow([frame, time, float(rmsd_value)])
    if args.output_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
        ax.plot([scaled_time(u.trajectory[frame], args.time_scale, args.time_offset) for frame in range(len(values))], values[:, 2], linewidth=0.8)
        ax.set_xlabel(f"Time ({args.time_unit})")
        ax.set_ylabel("Fitted RMSD (Å)")
        ax.set_title(f"RMSD relative to frame {args.reference_frame}")
        ax.grid(alpha=0.25)
        plot = Path(args.output_plot)
        plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot, dpi=300)
        plt.close(fig)
    print(f"Analysed {len(values)} frames using {len(selection)} selected atoms; output: {out}")
    print("RMSD includes least-squares fitting to the selected reference atoms; confirm selection and PBC preparation.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
