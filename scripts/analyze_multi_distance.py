#!/usr/bin/env python3
"""Calculate many named atom-pair distances from a JSON configuration."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import warnings
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
from copper_trafficking.distance import require_single_atom  # noqa: E402
from copper_trafficking.trajectory import scaled_time, valid_box_dimensions  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--config", required=True, help="JSON file containing named atom pairs")
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--output-plot", help="Optional plot of all distance series")
    parser.add_argument("--no-pbc", action="store_true")
    parser.add_argument("--presence-cutoff", type=float, help="Optional Å cutoff; add a binary <label>_present column")
    parser.add_argument("--trajectory-id", help="Stable trajectory identifier; defaults to the trajectory filename stem")
    parser.add_argument("--segment-id", default="0", help="Trajectory segment identifier (default: 0)")
    args = parser.parse_args()
    try:
        cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"Could not read config JSON: {exc}") from exc
    pairs = cfg.get("pairs")
    if not isinstance(pairs, list) or not pairs:
        raise SystemExit("Config must contain a non-empty 'pairs' list.")
    labels = [str(p.get("label", "")).strip() for p in pairs if isinstance(p, dict)]
    if len(labels) != len(pairs) or any(not x for x in labels):
        raise SystemExit("Every pair must be an object with a non-empty 'label'.")
    if len(set(labels)) != len(labels):
        raise SystemExit("Pair labels must be unique.")
    for pair in pairs:
        if not isinstance(pair.get("selection_a"), str) or not isinstance(pair.get("selection_b"), str):
            raise SystemExit("Every pair needs string fields 'selection_a' and 'selection_b'.")
    if args.presence_cutoff is not None and (not np.isfinite(args.presence_cutoff) or args.presence_cutoff <= 0):
        raise SystemExit("--presence-cutoff must be a finite positive distance in Å.")

    scale = float(cfg.get("time_scale", 1.0))
    offset = float(cfg.get("time_offset", 0.0))
    time_unit = str(cfg.get("time_unit", "ps"))
    if not np.isfinite(scale) or not np.isfinite(offset) or not time_unit.strip():
        raise SystemExit("time_scale/time_offset must be finite and time_unit non-empty.")

    try:
        import MDAnalysis as mda
        from MDAnalysis.lib.distances import calc_bonds
    except ImportError as exc:
        raise SystemExit("MDAnalysis is required. Install with: python -m pip install -e .") from exc
    universe = mda.Universe(args.topology, args.trajectory)
    selected = []
    for pair in pairs:
        atoms_a = universe.select_atoms(pair["selection_a"])
        atoms_b = universe.select_atoms(pair["selection_b"])
        try:
            require_single_atom(atoms_a, f"{pair['label']} selection_a")
            require_single_atom(atoms_b, f"{pair['label']} selection_b")
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        selected.append((pair["label"], atoms_a, atoms_b))

    rows = []
    no_box_frames = 0
    trajectory_id = args.trajectory_id or Path(args.trajectory).stem
    for ts in universe.trajectory:
        has_box = valid_box_dimensions(ts.dimensions)
        pbc = bool(has_box and not args.no_pbc)
        if not has_box:
            no_box_frames += 1
        row = {
            "trajectory_id": trajectory_id,
            "segment_id": str(args.segment_id),
            "frame": int(ts.frame),
            f"time_{time_unit}": scaled_time(ts, scale, offset),
            "pbc_applied": pbc,
        }
        for label, atom_a, atom_b in selected:
            distance = float(calc_bonds(atom_a.positions, atom_b.positions, box=ts.dimensions if pbc else None)[0])
            row[label + "_distance_A"] = distance
            if args.presence_cutoff is not None:
                row[label + "_present"] = int(distance <= args.presence_cutoff)
        rows.append(row)
    if not rows:
        raise SystemExit("The trajectory contains no frames.")
    if not args.no_pbc and no_box_frames:
        warnings.warn(f"{no_box_frames}/{len(rows)} frames lack valid unit-cell dimensions; no PBC correction in those frames.")

    out = Path(args.output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    presence_note = f"; presence cutoff={args.presence_cutoff:g} Å" if args.presence_cutoff is not None else ""
    print(f"Analysed {len(rows)} frames across {len(selected)} distance pair(s){presence_note}; wrote {out}")

    if args.output_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
        x = [r[f"time_{time_unit}"] for r in rows]
        for label, _, _ in selected:
            ax.plot(x, [r[label + "_distance_A"] for r in rows], linewidth=0.8, label=label)
        ax.set_xlabel(f"Time ({time_unit})")
        ax.set_ylabel("Distance (Å)")
        ax.legend(fontsize="small", ncol=1)
        ax.grid(alpha=0.25)
        plot = Path(args.output_plot)
        plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot, dpi=300)
        plt.close(fig)
        print(f"Plot written: {plot}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
