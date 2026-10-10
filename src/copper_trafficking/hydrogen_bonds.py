"""Geometry-based hydrogen-bond detection requiring topology bonds."""
from __future__ import annotations
from typing import Any
import numpy as np

from .trajectory import atom_identity, scaled_time, valid_box_dimensions


def _donor_hydrogen_pairs(universe: Any, donor_indices: set[int], hydrogen_indices: set[int]) -> list[tuple[int, int]]:
    """Match selected hydrogens to selected donors using topology bonds.

    This intentionally refuses to infer bonds from a distance-only heuristic.
    Supply a topology containing bond information (e.g. a suitable PSF or
    topology format with bonds) so donor identity is not guessed.
    """
    try:
        bonds = np.asarray(universe.bonds.indices, dtype=int)
    except Exception as exc:
        raise ValueError("Topology bond data are required to pair donor and hydrogen selections.") from exc
    if bonds.size == 0:
        raise ValueError("Topology contains no bonds; cannot validate donor–hydrogen pairs.")
    if bonds.ndim != 2 or bonds.shape[1] != 2:
        raise ValueError("Unexpected topology bond-index array shape.")
    pairs = set()
    for a, b in bonds:
        a, b = int(a), int(b)
        if a in donor_indices and b in hydrogen_indices:
            pairs.add((a, b))
        elif b in donor_indices and a in hydrogen_indices:
            pairs.add((b, a))
    if not pairs:
        raise ValueError("No donor–hydrogen bonds connect the chosen selections. Verify selections/topology.")
    return sorted(pairs)


def hydrogen_bond_contacts(
    universe: Any,
    donor_selection: str,
    hydrogen_selection: str,
    acceptor_selection: str,
    *,
    distance_cutoff_angstrom: float = 3.5,
    angle_cutoff_degrees: float = 120.0,
    time_scale: float = 1.0,
    time_offset: float = 0.0,
    use_pbc: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Find D–H···A contacts with D–A distance and D–H–A angle criteria.

    The angle is measured at H and must be greater than or equal to the
    specified threshold (e.g. 120°). Distances are in Å, angles in degrees.
    This is a geometric screen, not a quantum-mechanical bond characterization.
    """
    from MDAnalysis.lib.distances import calc_angles, calc_bonds

    if not (np.isfinite(distance_cutoff_angstrom) and distance_cutoff_angstrom > 0):
        raise ValueError("distance_cutoff_angstrom must be finite and positive.")
    if not (np.isfinite(angle_cutoff_degrees) and 0 < angle_cutoff_degrees <= 180):
        raise ValueError("angle_cutoff_degrees must be in (0, 180].")
    if not np.isfinite(time_scale) or not np.isfinite(time_offset):
        raise ValueError("time_scale and time_offset must be finite.")

    donors = universe.select_atoms(donor_selection)
    hydrogens = universe.select_atoms(hydrogen_selection)
    acceptors = universe.select_atoms(acceptor_selection)
    if len(donors) == 0 or len(hydrogens) == 0 or len(acceptors) == 0:
        raise ValueError("Donor, hydrogen, and acceptor selections must all be non-empty.")
    bonded_pairs = _donor_hydrogen_pairs(
        universe, set(int(i) for i in donors.indices), set(int(i) for i in hydrogens.indices)
    )
    donor_by_index = {int(atom.index): atom for atom in donors}
    hydrogen_by_index = {int(atom.index): atom for atom in hydrogens}

    contacts: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for ts in universe.trajectory:
        dims = ts.dimensions
        pbc = bool(use_pbc and valid_box_dimensions(dims))
        box = dims if pbc else None
        frame_contacts = []
        for donor_idx, hydrogen_idx in bonded_pairs:
            donor = donor_by_index[donor_idx]
            hydrogen = hydrogen_by_index[hydrogen_idx]
            donors_xyz = np.repeat(donor.position[np.newaxis, :], len(acceptors), axis=0)
            hydrogens_xyz = np.repeat(hydrogen.position[np.newaxis, :], len(acceptors), axis=0)
            distances = calc_bonds(donors_xyz, acceptors.positions, box=box)
            angles = np.degrees(calc_angles(donors_xyz, hydrogens_xyz, acceptors.positions, box=box))
            for i, acceptor in enumerate(acceptors):
                if float(distances[i]) <= distance_cutoff_angstrom and float(angles[i]) >= angle_cutoff_degrees:
                    d = atom_identity(donor)
                    h = atom_identity(hydrogen)
                    a = atom_identity(acceptor)
                    frame_contacts.append({
                        "frame": int(ts.frame),
                        "time": scaled_time(ts, time_scale, time_offset),
                        "donor_atom_index": d["atom_index"], "donor_segid": d["segid"],
                        "donor_resindex": d["resindex"], "donor_resid": d["resid"],
                        "donor_resname": d["resname"], "donor_atom_name": d["atom_name"],
                        "hydrogen_atom_index": h["atom_index"], "hydrogen_segid": h["segid"],
                        "hydrogen_resindex": h["resindex"], "hydrogen_resid": h["resid"],
                        "hydrogen_resname": h["resname"], "hydrogen_atom_name": h["atom_name"],
                        "acceptor_atom_index": a["atom_index"], "acceptor_segid": a["segid"],
                        "acceptor_resindex": a["resindex"], "acceptor_resid": a["resid"],
                        "acceptor_resname": a["resname"], "acceptor_atom_name": a["atom_name"],
                        "DA_distance_A": float(distances[i]), "DHA_angle_deg": float(angles[i]),
                        "pbc_applied": pbc,
                    })
        contacts.extend(frame_contacts)
        summaries.append({
            "frame": int(ts.frame), "time": scaled_time(ts, time_scale, time_offset),
            "hydrogen_bond_geometry_count": len(frame_contacts), "pbc_applied": pbc,
        })
    return contacts, summaries
