"""Dihedral-angle analysis for explicitly defined four-atom groups."""
from __future__ import annotations
from typing import Any
import numpy as np

from .distance import require_single_atom
from .trajectory import scaled_time, valid_box_dimensions


def dihedral_series(
    universe: Any,
    definitions: list[dict[str, Any]],
    *,
    time_scale: float = 1.0,
    time_offset: float = 0.0,
    use_pbc: bool = True,
) -> list[dict[str, Any]]:
    """Calculate named dihedral angles in degrees for each trajectory frame.

    Each definition must have a unique non-empty ``label`` and exactly four
    MDAnalysis ``selections``. Each selection must resolve to one atom.
    Angles are returned in the conventional wrapped range [-180°, 180°].
    """
    from MDAnalysis.lib.distances import calc_dihedrals

    if not np.isfinite(time_scale) or not np.isfinite(time_offset):
        raise ValueError("time_scale and time_offset must be finite.")
    if not definitions:
        raise ValueError("At least one dihedral definition is required.")
    labels = []
    selections = []
    for index, definition in enumerate(definitions):
        label = str(definition.get("label", "")).strip()
        atoms = definition.get("selections")
        if not label:
            raise ValueError(f"Dihedral #{index + 1} has no label.")
        if label in labels:
            raise ValueError(f"Dihedral labels must be unique; duplicate: {label!r}.")
        if not isinstance(atoms, list) or len(atoms) != 4 or any(not isinstance(x, str) or not x.strip() for x in atoms):
            raise ValueError(f"Dihedral {label!r} needs exactly four non-empty selection strings.")
        resolved = []
        for atom_no, selection_string in enumerate(atoms, start=1):
            group = universe.select_atoms(selection_string)
            require_single_atom(group, f"{label} atom {atom_no} selection")
            resolved.append(group)
        labels.append(label)
        selections.append((label, resolved))

    rows = []
    for ts in universe.trajectory:
        pbc = bool(use_pbc and valid_box_dimensions(ts.dimensions))
        row: dict[str, Any] = {
            "frame": int(ts.frame),
            "time": scaled_time(ts, time_scale, time_offset),
            "pbc_applied": pbc,
        }
        for label, atoms in selections:
            coords = [atom.positions for atom in atoms]
            angle_rad = float(calc_dihedrals(*coords, box=ts.dimensions if pbc else None)[0])
            row[label + "_angle_deg"] = float(np.degrees(angle_rad))
        rows.append(row)
    return rows
