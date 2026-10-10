#!/usr/bin/env python3
"""Detect topology-bonded D-H pairs geometrically interacting with acceptors."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
from copper_trafficking.hydrogen_bonds import hydrogen_bond_contacts  # noqa: E402


def write_csv(path: str, rows: list[dict], empty_fields: list[str] | None = None) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else (empty_fields or [])
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--donor-selection", required=True)
    p.add_argument("--hydrogen-selection", required=True)
    p.add_argument("--acceptor-selection", required=True)
    p.add_argument("--distance-cutoff", type=float, default=3.5, help="D-A distance cutoff in Å")
    p.add_argument("--angle-cutoff", type=float, default=120.0, help="D-H-A angle threshold in degrees")
    p.add_argument("--output-contacts", required=True)
    p.add_argument("--output-summary", required=True)
    p.add_argument("--time-scale", type=float, default=1.0)
    p.add_argument("--time-offset", type=float, default=0.0)
    p.add_argument("--time-unit", default="ps")
    p.add_argument("--no-pbc", action="store_true")
    args = p.parse_args()
    if not args.time_unit.strip():
        raise SystemExit("--time-unit must not be empty")
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit("MDAnalysis is required. Install with: python -m pip install -e .") from exc
    try:
        universe = mda.Universe(args.topology, args.trajectory)
        contacts, summaries = hydrogen_bond_contacts(
            universe, args.donor_selection, args.hydrogen_selection, args.acceptor_selection,
            distance_cutoff_angstrom=args.distance_cutoff,
            angle_cutoff_degrees=args.angle_cutoff,
            time_scale=args.time_scale, time_offset=args.time_offset,
            use_pbc=not args.no_pbc,
        )
    except (ValueError, OSError) as exc:
        raise SystemExit(str(exc)) from exc
    for row in contacts + summaries:
        row[f"time_{args.time_unit}"] = row.pop("time")
    contact_fields = ["frame", f"time_{args.time_unit}", "donor_atom_index", "donor_segid", "donor_resid", "donor_resname", "donor_atom_name", "hydrogen_atom_index", "hydrogen_segid", "hydrogen_resid", "hydrogen_resname", "hydrogen_atom_name", "acceptor_atom_index", "acceptor_segid", "acceptor_resid", "acceptor_resname", "acceptor_atom_name", "DA_distance_A", "DHA_angle_deg", "pbc_applied"]
    summary_fields = ["frame", "hydrogen_bond_geometry_count", "pbc_applied", f"time_{args.time_unit}"]
    write_csv(args.output_contacts, contacts, contact_fields)
    write_csv(args.output_summary, summaries, summary_fields)
    print(f"Frames analysed: {len(summaries)}; geometry-qualified D-H···A contacts: {len(contacts)}")
    print(f"Contacts: {args.output_contacts}\nSummary: {args.output_summary}")
    print("Interpret output as a geometric hydrogen-bond screen, not proof of bond formation.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
