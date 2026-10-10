#!/usr/bin/env python3
"""Calculate and plot a framewise distance between two MDAnalysis selections."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys
import warnings
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
from copper_trafficking.distance import require_single_atom  # noqa: E402
from copper_trafficking.trajectory import scaled_time, valid_box_dimensions  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Framewise single-atom distance analysis; output distances in Å.")
    p.add_argument("--topology", required=True, help="Topology file, e.g. PSF or PDB")
    p.add_argument("--trajectory", required=True, help="Trajectory file, e.g. DCD or XTC")
    p.add_argument("--selection-a", required=True, help="MDAnalysis selection for atom A")
    p.add_argument("--selection-b", required=True, help="MDAnalysis selection for atom B")
    p.add_argument("--output-csv", required=True, help="Output CSV path")
    p.add_argument("--output-plot", help="Optional output plot path (PNG/PDF/SVG)")
    p.add_argument("--time-scale", type=float, default=1.0, help="Multiply MDAnalysis frame time by this factor (default: 1)")
    p.add_argument("--time-offset", type=float, default=0.0, help="Add offset to scaled frame time (default: 0)")
    p.add_argument("--time-unit", default="ps", help="Label for reported time after conversion (default: ps)")
    p.add_argument("--no-pbc", action="store_true", help="Disable periodic boundary correction")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not np.isfinite(args.time_scale) or not np.isfinite(args.time_offset):
        raise SystemExit("--time-scale and --time-offset must be finite numbers")
    if not args.time_unit.strip():
        raise SystemExit("--time-unit must not be empty")
    try:
        import MDAnalysis as mda
        from MDAnalysis.lib.distances import calc_bonds
    except ImportError as exc:
        raise SystemExit("MDAnalysis is required. Install with: python -m pip install -e .") from exc

    universe = mda.Universe(args.topology, args.trajectory)
    atom_a = universe.select_atoms(args.selection_a)
    atom_b = universe.select_atoms(args.selection_b)
    try:
        require_single_atom(atom_a, "--selection-a")
        require_single_atom(atom_b, "--selection-b")
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    rows = []
    invalid_box_frames = 0
    for ts in universe.trajectory:
        dims = ts.dimensions
        has_box = valid_box_dimensions(dims)
        use_pbc = bool(has_box and not args.no_pbc)
        if not has_box:
            invalid_box_frames += 1
        distance = float(calc_bonds(atom_a.positions, atom_b.positions, box=dims if use_pbc else None)[0])
        time = scaled_time(ts, args.time_scale, args.time_offset)
        rows.append((int(ts.frame), time, distance, use_pbc))
    if not rows:
        raise SystemExit("The trajectory contains no frames.")
    if not args.no_pbc and invalid_box_frames:
        warnings.warn(f"{invalid_box_frames}/{len(rows)} frames lack a valid unit cell; Cartesian distances were used in those frames.")

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["frame", f"time_{args.time_unit}", "distance_A", "pbc_applied"])
        writer.writerows(rows)

    if args.output_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        out_plot = Path(args.output_plot)
        out_plot.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
        ax.plot([r[1] for r in rows], [r[2] for r in rows], linewidth=0.8)
        ax.set_xlabel(f"Time ({args.time_unit})")
        ax.set_ylabel("Distance (Å)")
        ax.set_title("Selected-atom distance")
        ax.grid(alpha=0.25)
        fig.savefig(out_plot, dpi=300)
        plt.close(fig)

    print(f"Frames analysed: {len(rows)}")
    print(f"PBC applied: {sum(r[3] for r in rows)}/{len(rows)} frames")
    print(f"CSV written: {out_csv}")
    if args.output_plot:
        print(f"Plot written: {args.output_plot}")
    print(f"Selection A: {args.selection_a}\nSelection B: {args.selection_b}")
    print("Verify selection identities and trajectory time metadata before interpretation.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
