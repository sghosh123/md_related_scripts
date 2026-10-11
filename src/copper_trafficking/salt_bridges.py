"""Geometric salt-bridge/contact analysis with explicit selections and provenance.

The output is a distance-based contact screen. A cutoff crossing is not, by itself,
proof of a salt bridge, protonation state, or electrostatic stabilization.
"""
from __future__ import annotations

from typing import Any, Iterable, Iterator
import numpy as np

from .trajectory import scaled_time, valid_box_dimensions


def iter_salt_bridge_geometry(
    universe: Any,
    definitions: Iterable[dict[str, Any]],
    *,
    cutoff_angstrom: float = 4.0,
    trajectory_id: str = "trajectory",
    segment_id: str | int = 0,
    time_scale: float = 1.0,
    time_offset: float = 0.0,
    time_unit: str = "ps",
    use_pbc: bool = True,
) -> Iterator[dict[str, Any]]:
    """Yield one row per configured atom-group pair and trajectory frame.

    Each definition has ``name``, ``selection_a`` and ``selection_b`` fields. The
    minimum heavy/charged-group atom-pair distance is recorded alongside the atom
    identities. Selection validity is checked before iterating so a typo cannot be
    silently reported as zero contacts.
    """
    from MDAnalysis.lib.distances import distance_array

    if not np.isfinite(cutoff_angstrom) or cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be a finite positive number in Å.")
    if not np.isfinite(time_scale) or not np.isfinite(time_offset):
        raise ValueError("time_scale and time_offset must be finite.")
    if not isinstance(time_unit, str) or not time_unit.strip():
        raise ValueError("time_unit must be a non-empty string.")
    defs = list(definitions)
    if not defs:
        raise ValueError("At least one salt-bridge/contact definition is required.")

    prepared = []
    names = set()
    for i, item in enumerate(defs):
        if not isinstance(item, dict):
            raise ValueError(f"Definition {i} must be an object.")
        name = str(item.get("name", "")).strip()
        sel_a = str(item.get("selection_a", "")).strip()
        sel_b = str(item.get("selection_b", "")).strip()
        if not name or name in names:
            raise ValueError(f"Definition {i} has an empty or duplicate name: {name!r}.")
        names.add(name)
        if not sel_a or not sel_b:
            raise ValueError(f"Definition {name!r} requires selection_a and selection_b.")
        group_a = universe.select_atoms(sel_a)
        group_b = universe.select_atoms(sel_b)
        if len(group_a) == 0:
            raise ValueError(f"Definition {name!r}: selection_a matched no atoms: {sel_a!r}")
        if len(group_b) == 0:
            raise ValueError(f"Definition {name!r}: selection_b matched no atoms: {sel_b!r}")
        prepared.append((name, sel_a, sel_b, group_a, group_b))

    for ts in universe.trajectory:
        dimensions = ts.dimensions if use_pbc and valid_box_dimensions(ts.dimensions) else None
        for name, sel_a, sel_b, group_a, group_b in prepared:
            distances = distance_array(group_a.positions, group_b.positions, box=dimensions)
            ai, bi = np.unravel_index(int(np.argmin(distances)), distances.shape)
            atom_a, atom_b = group_a[int(ai)], group_b[int(bi)]
            distance = float(distances[ai, bi])
            yield {
                "trajectory_id": str(trajectory_id), "segment_id": str(segment_id),
                "local_frame": int(ts.frame), "time": scaled_time(ts, time_scale, time_offset),
                "time_unit": time_unit, "pbc_applied": dimensions is not None,
                "bridge_name": name, "selection_a": sel_a, "selection_b": sel_b,
                "atom_a_index": int(atom_a.index), "atom_a_segid": str(getattr(atom_a, "segid", "")),
                "atom_a_resindex": int(atom_a.residue.ix), "atom_a_resid": int(atom_a.resid),
                "atom_a_resname": str(atom_a.resname), "atom_a_name": str(atom_a.name),
                "atom_b_index": int(atom_b.index), "atom_b_segid": str(getattr(atom_b, "segid", "")),
                "atom_b_resindex": int(atom_b.residue.ix), "atom_b_resid": int(atom_b.resid),
                "atom_b_resname": str(atom_b.resname), "atom_b_name": str(atom_b.name),
                "minimum_distance_A": distance,
                "within_cutoff": bool(distance <= cutoff_angstrom),
                "cutoff_A": float(cutoff_angstrom),
            }


def salt_bridge_geometry(*args: Any, **kwargs: Any) -> list[dict[str, Any]]:
    """Return all salt-bridge geometry rows as a list (convenient for small tests)."""
    return list(iter_salt_bridge_geometry(*args, **kwargs))
