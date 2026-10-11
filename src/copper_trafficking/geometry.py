"""Geometry helpers for interfacial features with explicit Å and PBC conventions."""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from MDAnalysis.lib.distances import calc_bonds, minimize_vectors

from .trajectory import valid_box_dimensions


def valid_box(box: Sequence[float] | np.ndarray | None) -> np.ndarray | None:
    """Return a valid MDAnalysis unit cell, or ``None`` if dimensions are unusable.

    MDAnalysis represents boxes as ``[a, b, c, alpha, beta, gamma]`` with lengths in Å
    and angles in degrees. This helper intentionally avoids attempting PBC correction on
    missing, malformed, or nonphysical box dimensions.
    """
    if not valid_box_dimensions(box):
        return None
    return np.asarray(box, dtype=float)


def minimum_image_vector(
    vector: Sequence[float], box: Sequence[float] | np.ndarray | None = None
) -> np.ndarray:
    """Return a directed displacement vector in Å using the minimum-image convention.

    ``vector`` is ``position_end - position_start``. If ``box`` is absent or invalid, the
    vector is returned unchanged; callers should record whether PBC correction was applied.
    """
    vec = np.asarray(vector, dtype=float)
    if vec.shape != (3,) or not np.all(np.isfinite(vec)):
        raise ValueError("vector must contain exactly three finite coordinates")
    pbc_box = valid_box(box)
    if pbc_box is None:
        return vec.copy()
    return np.asarray(minimize_vectors(vec.reshape(1, 3), pbc_box)[0], dtype=float)


def distance_angstrom_pbc(
    position_a: Sequence[float],
    position_b: Sequence[float],
    box: Sequence[float] | np.ndarray | None = None,
) -> float:
    """Calculate an interatomic distance in Å, applying PBC when ``box`` is valid."""
    a = np.asarray(position_a, dtype=float)
    b = np.asarray(position_b, dtype=float)
    if a.shape != (3,) or b.shape != (3,):
        raise ValueError("each position must contain exactly three coordinates")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("positions must contain only finite coordinates")
    pbc_box = valid_box(box)
    return float(calc_bonds(a.reshape(1, 3), b.reshape(1, 3), box=pbc_box)[0])


def angle_degrees(vector_a: Sequence[float], vector_b: Sequence[float]) -> float:
    """Return the angle between two non-zero vectors in degrees, bounded to [0, 180]."""
    a = np.asarray(vector_a, dtype=float)
    b = np.asarray(vector_b, dtype=float)
    if a.shape != (3,) or b.shape != (3,):
        raise ValueError("each vector must contain exactly three components")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("vectors must contain only finite components")
    norm_a = float(np.linalg.norm(a))
    norm_b = float(np.linalg.norm(b))
    if norm_a == 0.0 or norm_b == 0.0:
        raise ValueError("angle is undefined for a zero-length vector")
    cosine = float(np.dot(a, b) / (norm_a * norm_b))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
