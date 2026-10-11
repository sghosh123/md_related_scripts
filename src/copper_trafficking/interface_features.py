"""Inter-protein Cα distance and orientation features for copper-trafficking systems.

The default selections reproduce the current yeast MBD/chaperone feature definition.
They are system-specific starting points, not universal selections. Coordinates and
reported distances use MDAnalysis' Å convention; angles are reported in degrees.
"""
from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import Any

import MDAnalysis as mda
import numpy as np

from .geometry import angle_degrees, distance_angstrom_pbc, minimum_image_vector, valid_box
from .trajectory import scaled_time, valid_box_dimensions

DEFAULT_INTERFACE_SELECTIONS: dict[str, str] = {
    "A1": "segid A and resid 56 and name CA",
    "A2": "segid A and resid 63 and name CA",
    "B1": "segid A and resid 18 and name CA",
    "B2": "segid A and resid 26 and name CA",
    "X1": "segid B and resid 55 and name CA",
    "X2": "segid B and resid 62 and name CA",
    "Y1": "segid B and resid 18 and name CA",
    "Y2": "segid B and resid 25 and name CA",
}

REQUIRED_INTERFACE_ATOMS = tuple(DEFAULT_INTERFACE_SELECTIONS)
INTERFACE_DISTANCE_PAIRS: tuple[tuple[str, str, str], ...] = (
    ("d_A1X1_A", "A1", "X1"),
    ("d_A1Y1_A", "A1", "Y1"),
    ("d_B1X1_A", "B1", "X1"),
    ("d_B1Y1_A", "B1", "Y1"),
    ("d_A2X2_A", "A2", "X2"),
    ("d_A2Y2_A", "A2", "Y2"),
    ("d_B2X2_A", "B2", "X2"),
    ("d_B2Y2_A", "B2", "Y2"),
)
ANGLE_COLUMNS = ("angle_AX_deg", "angle_BX_deg", "angle_AY_deg", "angle_BY_deg")
FEATURE_COLUMNS = ANGLE_COLUMNS + tuple(pair[0] for pair in INTERFACE_DISTANCE_PAIRS)


def resolve_interface_selections(
    universe: mda.Universe,
    selections: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Resolve all required feature selections and fail on missing/ambiguous atoms."""
    selection_map = DEFAULT_INTERFACE_SELECTIONS if selections is None else selections
    missing = set(REQUIRED_INTERFACE_ATOMS) - set(selection_map)
    if missing:
        raise ValueError(f"selection mapping is missing keys: {sorted(missing)}")
    resolved: dict[str, Any] = {}
    for key in REQUIRED_INTERFACE_ATOMS:
        text = selection_map[key]
        if not isinstance(text, str) or not text.strip():
            raise ValueError(f"selection {key!r} must be a non-empty string")
        atoms = universe.select_atoms(text)
        if len(atoms) != 1:
            raise ValueError(
                f"selection {key!r} ({text!r}) resolved to {len(atoms)} atoms; expected exactly 1"
            )
        resolved[key] = atoms
    return resolved


def compute_frame_interface_features(
    universe: mda.Universe,
    selections: Mapping[str, str] | None = None,
) -> dict[str, float]:
    """Compute four orientation angles (degrees) and eight distances (Å) for current frame.

    The four orientation vectors are defined by A1→A2, B1→B2, X1→X2 and Y1→Y2.
    Each vector receives minimum-image correction when the current frame has valid box
    dimensions. Inter-atom distances use the same box convention.
    """
    atoms = resolve_interface_selections(universe, selections)
    return _features_from_resolved_atoms(universe, atoms)


def _features_from_resolved_atoms(universe: mda.Universe, atoms: Mapping[str, Any]) -> dict[str, float]:
    box = valid_box(universe.trajectory.ts.dimensions)
    positions = {key: np.asarray(group.positions[0], dtype=float).copy() for key, group in atoms.items()}
    vectors = {
        "A": minimum_image_vector(positions["A2"] - positions["A1"], box),
        "B": minimum_image_vector(positions["B2"] - positions["B1"], box),
        "X": minimum_image_vector(positions["X2"] - positions["X1"], box),
        "Y": minimum_image_vector(positions["Y2"] - positions["Y1"], box),
    }
    result: dict[str, float] = {
        "angle_AX_deg": angle_degrees(vectors["A"], vectors["X"]),
        "angle_BX_deg": angle_degrees(vectors["B"], vectors["X"]),
        "angle_AY_deg": angle_degrees(vectors["A"], vectors["Y"]),
        "angle_BY_deg": angle_degrees(vectors["B"], vectors["Y"]),
    }
    for column, atom_a, atom_b in INTERFACE_DISTANCE_PAIRS:
        result[column] = distance_angstrom_pbc(positions[atom_a], positions[atom_b], box)
    return result


def iter_interface_features(
    universe: mda.Universe,
    *,
    selections: Mapping[str, str] | None = None,
    trajectory_id: str = "trajectory",
    segment_id: str | int = 0,
    time_scale: float = 1.0,
    time_offset: float = 0.0,
    time_unit: str = "ps",
) -> Iterator[dict[str, object]]:
    """Yield one feature row per frame, preserving trajectory/segment provenance.

    ``time_unit`` labels the scaled time output; verify reader time metadata and provide
    ``time_scale``/``time_offset`` explicitly when the native reader values need conversion.
    Selection errors are raised before the returned iterator begins yielding rows.
    """
    if not np.isfinite(time_scale) or not np.isfinite(time_offset):
        raise ValueError("time_scale and time_offset must be finite numbers")
    if not isinstance(time_unit, str) or not time_unit.strip():
        raise ValueError("time_unit must be a non-empty string")
    resolved = resolve_interface_selections(universe, selections)

    def generate() -> Iterator[dict[str, object]]:
        for ts in universe.trajectory:
            box = valid_box(ts.dimensions)
            row: dict[str, object] = {
                "trajectory_id": str(trajectory_id),
                "segment_id": str(segment_id),
                "local_frame": int(ts.frame),
                "time": scaled_time(ts, time_scale, time_offset),
                "time_unit": time_unit,
                "time_scale": float(time_scale),
                "time_offset": float(time_offset),
                "pbc_applied": bool(valid_box_dimensions(ts.dimensions)),
            }
            row.update(_features_from_resolved_atoms(universe, resolved))
            yield row

    return generate()
