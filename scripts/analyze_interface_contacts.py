#!/usr/bin/env python3
"""Stream inter-protein heavy-atom contacts to a CSV file."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

import MDAnalysis as mda  # noqa: E402
from copper_trafficking.contacts import iter_interfacial_contacts  # noqa: E402

FIELDS = [
    "trajectory_id", "segment_id", "local_frame", "time", "time_unit", "pbc_applied",
    "atom_a_index", "atom_a_segid", "atom_a_resid", "atom_a_resname", "atom_a_name",
    "atom_b_index", "atom_b_segid", "atom_b_resid", "atom_b_resname", "atom_b_name", "distance_A",
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--selection-a", default="protein and segid A and not name H*")
    parser.add_argument("--selection-b", default="protein and segid B and not name H*")
    parser.add_argument("--cutoff", type=float, default=4.5, help="contact cutoff in Å")
    parser.add_argument("--trajectory-id", help="stable trajectory identifier; defaults to filename stem")
    parser.add_argument("--segment-id", default="0")
    parser.add_argument("--time-scale", type=float, default=1.0)
    parser.add_argument("--time-offset", type=float, default=0.0)
    parser.add_argument("--time-unit", default="ps", help="label for scaled time values")
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args(argv)

    universe = mda.Universe(args.topology, args.trajectory)
    try:
        rows = iter_interfacial_contacts(
            universe,
            args.selection_a,
            args.selection_b,
            args.cutoff,
            trajectory_id=args.trajectory_id or Path(args.trajectory).stem,
            segment_id=args.segment_id,
            time_scale=args.time_scale,
            time_offset=args.time_offset,
            time_unit=args.time_unit,
        )
        # Force selection/cutoff validation before creating output files.
        first = next(rows, None)
    except ValueError as exc:
        parser.error(str(exc))

    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        count = 0
        if first is not None:
            writer.writerow(first)
            count = 1
        for row in rows:
            writer.writerow(row)
            count += 1
    print(f"Wrote {count} atom-pair contact(s) to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
