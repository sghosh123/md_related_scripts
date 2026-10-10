#!/usr/bin/env python3
"""Calculate named torsions from a JSON file containing four atom selections each."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
from copper_trafficking.torsions import dihedral_series  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--config", required=True, help="JSON with a 'dihedrals' list")
    p.add_argument("--output-csv", required=True)
    p.add_argument("--time-scale", type=float, default=1.0)
    p.add_argument("--time-offset", type=float, default=0.0)
    p.add_argument("--time-unit", default="ps")
    p.add_argument("--no-pbc", action="store_true")
    args = p.parse_args()
    if not args.time_unit.strip():
        raise SystemExit("--time-unit must not be empty")
    try:
        cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"Could not read JSON config: {exc}") from exc
    definitions = cfg.get("dihedrals")
    if not isinstance(definitions, list) or not definitions:
        raise SystemExit("Config must contain a non-empty 'dihedrals' list.")
    try:
        import MDAnalysis as mda
        u = mda.Universe(args.topology, args.trajectory)
        rows = dihedral_series(
            u, definitions, time_scale=args.time_scale, time_offset=args.time_offset,
            use_pbc=not args.no_pbc,
        )
    except (ValueError, OSError) as exc:
        raise SystemExit(str(exc)) from exc
    for row in rows:
        row[f"time_{args.time_unit}"] = row.pop("time")
    out = Path(args.output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["frame", f"time_{args.time_unit}", "pbc_applied"]
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Calculated {len(definitions)} dihedral(s) over {len(rows)} frame(s); output: {out}")
    print("Verify atom ordering and chemical interpretation before describing an angle as phi, psi, or a coordination torsion.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
