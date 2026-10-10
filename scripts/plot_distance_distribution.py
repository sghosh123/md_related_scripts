#!/usr/bin/env python3
"""Plot a distribution from a named CSV column or a whitespace-delimited file."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np


def positive_int(value: str) -> int:
    try:
        result = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("bins must be a positive integer") from exc
    if result < 1:
        raise argparse.ArgumentTypeError("bins must be a positive integer")
    return result


def read_values(path: Path, column: str | None, column_index: int) -> np.ndarray:
    values = []
    if path.suffix.lower() == ".csv":
        if not column:
            raise ValueError("CSV input requires --column with the exact column name.")
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or column not in reader.fieldnames:
                raise ValueError(f"Column {column!r} not found. Available columns: {reader.fieldnames}")
            for row_no, row in enumerate(reader, start=2):
                value = row.get(column, "")
                if value in (None, ""):
                    continue
                try:
                    values.append(float(value))
                except ValueError as exc:
                    raise ValueError(f"Non-numeric value in {path}:{row_no}, column {column!r}: {value!r}") from exc
    else:
        if column:
            raise ValueError("--column is for CSV only; use --column-index for whitespace-delimited files.")
        with path.open(encoding="utf-8") as handle:
            for row_no, line in enumerate(handle, start=1):
                tokens = line.split()
                if not tokens or line.lstrip().startswith("#"):
                    continue
                if len(tokens) <= column_index:
                    raise ValueError(f"Line {row_no} has only {len(tokens)} column(s), not index {column_index}.")
                try:
                    values.append(float(tokens[column_index]))
                except ValueError:
                    # A leading header line is common in legacy .dat files; skip it
                    # only before any numeric observations have been collected.
                    if not values:
                        continue
                    raise ValueError(f"Non-numeric value in {path}:{row_no}, column index {column_index}.")
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        raise ValueError("No finite numeric values were found in the selected column.")
    return array


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--column", help="Column header in a CSV file")
    p.add_argument("--column-index", type=int, default=1, help="Zero-based column index for whitespace-delimited input (default 1)")
    p.add_argument("--bins", type=positive_int, default=50)
    p.add_argument("--output", required=True)
    p.add_argument("--title", default="Distance distribution")
    p.add_argument("--xlabel", default="Distance (Å)")
    args = p.parse_args()
    source = Path(args.input)
    if args.column_index < 0:
        raise SystemExit("--column-index must be zero or greater.")
    try:
        values = read_values(source, args.column, args.column_index)
    except (OSError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    ax.hist(values, bins=args.bins)
    ax.set_title(args.title)
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel("Count")
    ax.grid(axis="y", alpha=0.25)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Read {values.size} finite values; histogram written to {out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
