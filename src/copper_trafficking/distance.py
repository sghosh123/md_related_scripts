"""Distance calculations and input checks used by trajectory workflows."""
from __future__ import annotations
import numpy as np


def distance_angstrom(position_a: np.ndarray, position_b: np.ndarray) -> float:
    """Return Cartesian distance in Å. This function does not apply PBC."""
    a = np.asarray(position_a, dtype=float)
    b = np.asarray(position_b, dtype=float)
    if a.shape != (3,) or b.shape != (3,):
        raise ValueError("Each position must be a 3-element vector.")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("Positions must contain only finite coordinates.")
    return float(np.linalg.norm(a - b))


def require_single_atom(selection, label: str) -> None:
    """Raise a useful error unless an MDAnalysis AtomGroup has one atom."""
    count = len(selection)
    if count != 1:
        raise ValueError(
            f"{label} must select exactly one atom; it selected {count}. "
            "Check the selection string, segment/chain IDs, and residue numbering."
        )


def require_nonempty_selection(selection, label: str) -> None:
    """Raise a clear error when an MDAnalysis selection is empty."""
    count = len(selection)
    if count == 0:
        raise ValueError(f"{label} selected no atoms. Check the selection syntax and topology.")
