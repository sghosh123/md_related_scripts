#!/usr/bin/env python3
"""Write selected trajectory frames as individual PDB files."""
from __future__ import annotations
import argparse
from pathlib import Path


def parse_frames(value: str, n_frames: int) -> list[int]:
    """Parse frame indices or inclusive ranges: '0,10,20-25'."""
    values: set[int] = set()
    try:
        for token in value.split(","):
            token = token.strip()
            if not token:
                raise ValueError("Empty frame token (check for a leading, trailing, or repeated comma).")
            if "-" in token:
                start_text, end_text = token.split("-", 1)
                start, end = int(start_text), int(end_text)
                if end < start:
                    raise ValueError(f"Range end before start: {token}")
                values.update(range(start, end + 1))
            else:
                values.add(int(token))
    except ValueError as exc:
        raise ValueError(f"Invalid frame specification {value!r}: {exc}") from exc
    if not values:
        raise ValueError("No frame indices were requested.")
    invalid = sorted(x for x in values if x < 0 or x >= n_frames)
    if invalid:
        raise ValueError(f"Frame indices outside valid range 0..{n_frames-1}: {invalid}")
    return sorted(values)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--frames", required=True, help="Examples: '0,10,50' or '100-120,500'")
    p.add_argument("--selection", default="all", help="MDAnalysis selection to write (default: all)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--prefix", default="frame")
    args = p.parse_args()
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit("MDAnalysis is required. Install with: python -m pip install -e .") from exc
    u = mda.Universe(args.topology, args.trajectory)
    group = u.select_atoms(args.selection)
    if len(group) == 0:
        raise SystemExit("--selection returned no atoms.")
    try:
        frame_indices = parse_frames(args.frames, len(u.trajectory))
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    for frame in frame_indices:
        u.trajectory[frame]
        output = outdir / f"{args.prefix}_{frame:06d}.pdb"
        with mda.Writer(str(output), n_atoms=len(group)) as writer:
            writer.write(group)
    print(f"Wrote {len(frame_indices)} frame(s) to {outdir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
