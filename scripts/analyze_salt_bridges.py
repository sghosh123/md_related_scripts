#!/usr/bin/env python3
"""Compute configured salt-bridge-like charged-group distances for every frame."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

import MDAnalysis as mda  # noqa: E402
from copper_trafficking.salt_bridges import iter_salt_bridge_geometry  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--config", required=True, help="JSON list of explicit atom-group pair definitions")
    parser.add_argument("--cutoff", type=float, default=4.0, help="distance-only cutoff in Å")
    parser.add_argument("--trajectory-id", help="stable independent-trajectory ID; defaults to trajectory filename stem")
    parser.add_argument("--segment-id", default="0", help="file segment ID (not an independent replicate by itself)")
    parser.add_argument("--time-scale", type=float, default=1.0)
    parser.add_argument("--time-offset", type=float, default=0.0)
    parser.add_argument("--time-unit", default="ps")
    parser.add_argument("--no-pbc", action="store_true", help="disable minimum-image distances; only use after reviewing periodicity")
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args(argv)
    try:
        config = json.loads(Path(args.config).read_text(encoding="utf-8"))
        definitions = config.get("salt_bridges") if isinstance(config, dict) else config
        if not isinstance(definitions, list):
            raise ValueError("Config must be a JSON list or an object containing a 'salt_bridges' list.")
        universe = mda.Universe(args.topology, args.trajectory)
        rows = iter_salt_bridge_geometry(
            universe, definitions, cutoff_angstrom=args.cutoff,
            trajectory_id=args.trajectory_id or Path(args.trajectory).stem,
            segment_id=args.segment_id, time_scale=args.time_scale,
            time_offset=args.time_offset, time_unit=args.time_unit, use_pbc=not args.no_pbc,
        )
        # Force selection and configuration validation before creating the output.
        first = next(rows, None)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        parser.error(str(exc))
    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(first) if first is not None else ["trajectory_id", "segment_id", "local_frame", "time", "time_unit", "pbc_applied", "bridge_name", "selection_a", "selection_b", "minimum_distance_A", "within_cutoff", "cutoff_A"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        count = 0
        if first is not None:
            writer.writerow(first)
            count = 1
        for row in rows:
            writer.writerow(row)
            count += 1
    print(f"Wrote {count} frame/definition rows to {output}")
    print("Note: within_cutoff is a geometric proximity flag, not proof of salt-bridge energetics.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
