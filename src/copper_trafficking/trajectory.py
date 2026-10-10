"""Shared trajectory helpers with explicit periodic-cell validation."""
from __future__ import annotations

from typing import Any
import numpy as np


def valid_box_dimensions(dimensions: Any) -> bool:
    """Return True when MDAnalysis unit-cell dimensions are usable.

    MDAnalysis uses [a, b, c, alpha, beta, gamma] in Å and degrees.
    This check is deliberately conservative: all six values must be finite,
    positive cell lengths are required, and angles must be strictly between
    0 and 180 degrees.
    """
    if dimensions is None:
        return False
    try:
        dims = np.asarray(dimensions, dtype=float)
    except (TypeError, ValueError):
        return False
    if dims.shape != (6,) or not np.all(np.isfinite(dims)):
        return False
    if np.any(dims[:3] <= 0):
        return False
    if np.any(dims[3:] <= 0) or np.any(dims[3:] >= 180):
        return False
    return True


def atom_identity(atom: Any) -> dict[str, Any]:
    """Return a stable, compact identity record for an MDAnalysis atom."""
    residue = atom.residue
    return {
        "atom_index": int(atom.index),
        "segid": str(getattr(atom, "segid", "")),
        "resid": int(atom.resid),
        "resname": str(atom.resname),
        "atom_name": str(atom.name),
        "resindex": int(residue.ix),
    }


def scaled_time(ts: Any, scale: float = 1.0, offset: float = 0.0) -> float:
    """Scale MDAnalysis frame time and apply an explicit offset."""
    value = float(getattr(ts, "time", ts.frame))
    return value * scale + offset
