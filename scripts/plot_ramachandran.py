#!/usr/bin/env python3
"""Plot a phi-versus-psi scatter plot from an analysis CSV."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np


def read_pair(path: Path, phi_column: str, psi_column: str) -> tuple[np.ndarray, np.ndarray]:
    phi, psi = [], []
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or phi_column not in reader.fieldnames or psi_column not in reader.fieldnames:
            raise ValueError(f"Missing phi/psi columns. Available: {reader.fieldnames}")
        for row_num, row in enumerate(reader, start=2):
            try:
                x, y = float(row[phi_column]), float(row[psi_column])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid numeric phi/psi value at CSV row {row_num}.") from exc
            if np.isfinite(x) and np.isfinite(y):
                phi.append(x)
                psi.append(y)
    if not phi:
        raise ValueError("No finite phi/psi pairs were found.")
    return np.asarray(phi), np.asarray(psi)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--phi-column", required=True)
    p.add_argument("--psi-column", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--title", default="Ramachandran-style phi/psi plot")
    p.add_argument("--point-size", type=float, default=5.0)
    args = p.parse_args()
    if args.point_size <= 0 or not np.isfinite(args.point_size):
        raise SystemExit("--point-size must be finite and positive.")
    try:
        phi, psi = read_pair(Path(args.input), args.phi_column, args.psi_column)
    except (OSError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 6.0), constrained_layout=True)
    ax.scatter(phi, psi, s=args.point_size, alpha=0.35, rasterized=True)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-180, 181, 60))
    ax.set_xlabel("Phi (°)")
    ax.set_ylabel("Psi (°)")
    ax.set_title(args.title)
    ax.grid(alpha=0.2)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Plotted {len(phi)} phi/psi pairs to {out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
