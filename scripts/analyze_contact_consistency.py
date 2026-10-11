#!/usr/bin/env python3
"""Estimate residue-pair contact occupancy consistency across trajectories."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from copper_trafficking.contact_consistency import summarize_contact_consistency  # noqa: E402


def read_many(paths, sep=None):
    frames = []
    for path in paths:
        delimiter = sep or ("\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ",")
        frames.append(pd.read_csv(path, sep=delimiter))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contacts", nargs="+", required=True, help="One or more interface contact CSV/TSV files")
    parser.add_argument("--frames", nargs="+", required=True, help="Matching frame manifests created by analyze_interface_contacts.py")
    parser.add_argument("--outdir", default="contact_consistency")
    parser.add_argument("--sep", choices=[",", "\\t"], default=None)
    parser.add_argument("--minimum-trajectories", type=int, default=1)
    args = parser.parse_args(argv)
    sep = "\t" if args.sep == "\\t" else args.sep
    try:
        contacts = read_many(args.contacts, sep)
        frame_manifest = read_many(args.frames, sep)
        per_traj, summary = summarize_contact_consistency(
            contacts, frame_manifest, minimum_trajectories=args.minimum_trajectories,
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        parser.error(str(exc))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    per_traj.to_csv(outdir / "contact_occupancy_by_trajectory.tsv", sep="\t", index=False)
    summary.to_csv(outdir / "contact_consistency_summary.tsv", sep="\t", index=False)
    print(f"Wrote {len(per_traj)} pair/trajectory and {len(summary)} summary rows to {outdir.resolve()}")
    print("Trajectory IDs are treated as independent replicates; segment IDs are pooled within each trajectory.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
