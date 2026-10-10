#!/usr/bin/env python3
"""Analyze per-frame neighbors around an atom or group selection."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
from copper_trafficking.contacts import neighbor_contacts, residue_contact_frequency, residue_pair_contact_frequency  # noqa: E402


def write_csv(path: str | Path, rows: list[dict], empty_fields: list[str] | None = None) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        fields = list(rows[0])
    else:
        # Header-only outputs are still useful when there are zero contacts.
        fields = empty_fields or []
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--center-selection", required=True, help="MDAnalysis selection for center atom(s)")
    p.add_argument("--neighbor-selection", required=True, help="Candidate neighbors; for water analysis, e.g. 'resname TIP3 or resname HOH'")
    p.add_argument("--cutoff", required=True, type=float, help="Contact cutoff in Å")
    p.add_argument("--output-contacts", required=True, help="Long-form atom-pair contact CSV")
    p.add_argument("--output-summary", required=True, help="Per-frame neighbor count CSV")
    p.add_argument("--output-residue-frequency", help="Optional per-neighbor-residue contact occupancy CSV")
    p.add_argument("--output-pair-frequency", help="Optional center-residue/neighbor-residue interaction occupancy CSV")
    p.add_argument("--output-plot", help="Optional per-frame unique-neighbor count plot")
    p.add_argument("--time-scale", type=float, default=1.0)
    p.add_argument("--time-offset", type=float, default=0.0)
    p.add_argument("--time-unit", default="ps")
    p.add_argument("--no-pbc", action="store_true")
    p.add_argument("--exclude-center-residues", action="store_true", help="Exclude all candidate atoms in center residues")
    args = p.parse_args()
    if not args.time_unit.strip():
        raise SystemExit("--time-unit must not be empty")
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit("MDAnalysis is required. Install with: python -m pip install -e .") from exc
    try:
        u = mda.Universe(args.topology, args.trajectory)
        contacts, summaries = neighbor_contacts(
            u, args.center_selection, args.neighbor_selection, args.cutoff,
            time_scale=args.time_scale, time_offset=args.time_offset,
            use_pbc=not args.no_pbc, exclude_center_residues=args.exclude_center_residues,
        )
    except (ValueError, OSError) as exc:
        raise SystemExit(str(exc)) from exc
    # Rename generic time column according to explicit reporting unit.
    for row in contacts + summaries:
        row[f"time_{args.time_unit}"] = row.pop("time")
    contact_fields = ["frame", f"time_{args.time_unit}", "center_atom_index", "center_segid", "center_resindex", "center_resid", "center_resname", "center_atom_name", "neighbor_atom_index", "neighbor_segid", "neighbor_resindex", "neighbor_resid", "neighbor_resname", "neighbor_atom_name", "distance_A", "pbc_applied"]
    summary_fields = ["frame", "neighbor_atom_count", "neighbor_residue_count", "contact_pair_count", "pbc_applied", f"time_{args.time_unit}"]
    write_csv(args.output_contacts, contacts, contact_fields)
    write_csv(args.output_summary, summaries, summary_fields)
    print(f"Frames analysed: {len(summaries)}; atom-pair contact rows: {len(contacts)}")
    print(f"Contacts: {args.output_contacts}\nPer-frame summary: {args.output_summary}")
    if args.output_residue_frequency:
        frequency = residue_contact_frequency(contacts, n_frames=len(summaries))
        frequency_fields = ["neighbor_segid", "neighbor_resindex", "neighbor_resid", "neighbor_resname", "frames_in_contact", "contact_frame_fraction", "minimum_distance_A"]
        write_csv(args.output_residue_frequency, frequency, frequency_fields)
        print(f"Per-residue contact occupancy: {args.output_residue_frequency}")
    if args.output_pair_frequency:
        pair_frequency = residue_pair_contact_frequency(contacts, n_frames=len(summaries))
        pair_fields = ["center_segid", "center_resindex", "center_resid", "center_resname", "neighbor_segid", "neighbor_resindex", "neighbor_resid", "neighbor_resname", "frames_in_contact", "contact_frame_fraction", "minimum_distance_A"]
        write_csv(args.output_pair_frequency, pair_frequency, pair_fields)
        print(f"Residue-pair interaction occupancy: {args.output_pair_frequency}")
    if args.output_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
        x = [row[f"time_{args.time_unit}"] for row in summaries]
        y = [row["neighbor_atom_count"] for row in summaries]
        ax.plot(x, y, linewidth=0.9)
        ax.set_xlabel(f"Time ({args.time_unit})")
        ax.set_ylabel("Unique neighbor atoms within cutoff")
        ax.set_title(f"Neighbor count (cutoff {args.cutoff:g} Å)")
        ax.grid(alpha=0.25)
        out = Path(args.output_plot)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300)
        plt.close(fig)
        print(f"Plot written: {out}")
    if not contacts:
        print("Note: no atom pairs met the cutoff; summary still contains one row per frame.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
