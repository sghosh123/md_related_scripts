#!/usr/bin/env python3
"""Compute the eight-distance/four-angle interface feature vector for each frame."""
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
from copper_trafficking.interface_features import (  # noqa: E402
    DEFAULT_INTERFACE_SELECTIONS,
    iter_interface_features,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, help="topology or structure file")
    parser.add_argument("--trajectory", required=True, help="trajectory file")
    parser.add_argument("--output-csv", required=True, help="output CSV path")
    parser.add_argument("--selections-json", help="JSON file with a 'selections' mapping")
    parser.add_argument("--trajectory-id", help="stable trajectory identifier (defaults to filename stem)")
    parser.add_argument("--segment-id", default="0", help="trajectory-segment identifier (default: 0)")
    parser.add_argument("--time-scale", type=float, default=1.0, help="multiply reader time by this factor")
    parser.add_argument("--time-offset", type=float, default=0.0, help="add this offset after scaling time")
    parser.add_argument("--time-unit", default="ps", help="label for scaled time values (default: ps)")
    args = parser.parse_args(argv)

    selections = DEFAULT_INTERFACE_SELECTIONS
    if args.selections_json:
        try:
            payload = json.loads(Path(args.selections_json).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            parser.error(f"could not read selection JSON: {exc}")
        if not isinstance(payload, dict) or not isinstance(payload.get("selections"), dict):
            parser.error("selection JSON must contain an object-valued 'selections' mapping")
        selections = payload["selections"]

    universe = mda.Universe(args.topology, args.trajectory)
    trajectory_id = args.trajectory_id or Path(args.trajectory).stem
    try:
        rows = iter_interface_features(
            universe,
            selections=selections,
            trajectory_id=trajectory_id,
            segment_id=args.segment_id,
            time_scale=args.time_scale,
            time_offset=args.time_offset,
            time_unit=args.time_unit,
        )
        first = next(rows, None)
    except (ValueError, KeyError) as exc:
        parser.error(str(exc))
    if first is None:
        parser.error("trajectory contains no frames")

    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(first)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(first)
        count = 1
        for row in rows:
            writer.writerow(row)
            count += 1
    print(f"Analysed {count} frame(s); wrote interface features to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
