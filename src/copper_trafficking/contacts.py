"""Periodic-boundary-aware atom-neighbor/contact analysis."""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator
from typing import Any
import numpy as np

from .trajectory import atom_identity, scaled_time, valid_box_dimensions


def neighbor_contacts(
    universe: Any,
    center_selection: str,
    neighbor_selection: str,
    cutoff_angstrom: float,
    *,
    time_scale: float = 1.0,
    time_offset: float = 0.0,
    use_pbc: bool = True,
    exclude_center_residues: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return (contact rows, per-frame summary rows).

    Every atom-pair contact at or below ``cutoff_angstrom`` is emitted once.
    Summary counts distinguish unique neighbor atoms from unique residues.
    Distances are in Å; trajectory time inherits the reader's time unit and
    then applies the caller's explicit scale/offset.
    """
    from MDAnalysis.lib.distances import distance_array

    if not np.isfinite(cutoff_angstrom) or cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be a finite positive number.")
    if not np.isfinite(time_scale) or not np.isfinite(time_offset):
        raise ValueError("time_scale and time_offset must be finite numbers.")

    centers = universe.select_atoms(center_selection)
    candidates = universe.select_atoms(neighbor_selection)
    if len(centers) == 0:
        raise ValueError(f"Center selection selected no atoms: {center_selection!r}")
    if len(candidates) == 0:
        raise ValueError(f"Neighbor selection selected no atoms: {neighbor_selection!r}")

    center_indices = set(int(i) for i in centers.indices)
    center_resindices = set(int(atom.resindex) for atom in centers)
    candidate_mask = np.array([int(i) not in center_indices for i in candidates.indices], dtype=bool)
    if exclude_center_residues:
        candidate_mask &= np.array([int(atom.resindex) not in center_resindices for atom in candidates], dtype=bool)
    candidates = candidates[candidate_mask]
    if len(candidates) == 0:
        raise ValueError("No candidate atoms remain after excluding center atoms/residues.")

    contact_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for ts in universe.trajectory:
        dims = ts.dimensions
        pbc = bool(use_pbc and valid_box_dimensions(dims))
        box = dims if pbc else None
        dmat = distance_array(centers.positions, candidates.positions, box=box)
        pair_indices = np.argwhere(dmat <= cutoff_angstrom)
        frame_neighbor_indices: set[int] = set()
        frame_residues: set[tuple[Any, ...]] = set()
        current_contacts: list[dict[str, Any]] = []
        for ci, ni in pair_indices:
            center_atom = centers[int(ci)]
            neighbor_atom = candidates[int(ni)]
            cinfo = atom_identity(center_atom)
            ninfo = atom_identity(neighbor_atom)
            distance = float(dmat[int(ci), int(ni)])
            frame_neighbor_indices.add(int(neighbor_atom.index))
            frame_residues.add((ninfo["resindex"], ninfo["segid"], ninfo["resid"], ninfo["resname"]))
            current_contacts.append({
                "frame": int(ts.frame),
                "time": scaled_time(ts, time_scale, time_offset),
                "center_atom_index": cinfo["atom_index"],
                "center_segid": cinfo["segid"],
                "center_resindex": cinfo["resindex"],
                "center_resid": cinfo["resid"],
                "center_resname": cinfo["resname"],
                "center_atom_name": cinfo["atom_name"],
                "neighbor_atom_index": ninfo["atom_index"],
                "neighbor_segid": ninfo["segid"],
                "neighbor_resindex": ninfo["resindex"],
                "neighbor_resid": ninfo["resid"],
                "neighbor_resname": ninfo["resname"],
                "neighbor_atom_name": ninfo["atom_name"],
                "distance_A": distance,
                "pbc_applied": pbc,
            })
        contact_rows.extend(current_contacts)
        summary_rows.append({
            "frame": int(ts.frame),
            "time": scaled_time(ts, time_scale, time_offset),
            "neighbor_atom_count": len(frame_neighbor_indices),
            "neighbor_residue_count": len(frame_residues),
            "contact_pair_count": len(current_contacts),
            "pbc_applied": pbc,
        })
    return contact_rows, summary_rows


def residue_contact_frequency(contact_rows: list[dict[str, Any]], n_frames: int | None = None) -> list[dict[str, Any]]:
    """Aggregate contacts into per-residue frame occupancy.

    A residue is counted at most once per frame even if many of its atoms
    contact the center selection. Occupancy is the number of distinct frames
    with one or more contacts, not a time-weighted probability.
    """
    per_frame: dict[tuple[int, tuple[Any, ...]], float] = {}
    for row in contact_rows:
        # Include topology residue index to distinguish duplicated residue
        # numbers/names in the same segment (possible in concatenated systems).
        key_res = (
            row["neighbor_segid"],
            int(row.get("neighbor_resindex", row["neighbor_resid"])),
            int(row["neighbor_resid"]),
            row["neighbor_resname"],
        )
        per_frame[(int(row["frame"]), key_res)] = min(
            float(row["distance_A"]), per_frame.get((int(row["frame"]), key_res), float("inf"))
        )
    frames_by_residue: dict[tuple[Any, ...], set[int]] = defaultdict(set)
    mindist_by_residue: dict[tuple[Any, ...], float] = defaultdict(lambda: float("inf"))
    for (frame, residue), minimum in per_frame.items():
        frames_by_residue[residue].add(frame)
        mindist_by_residue[residue] = min(mindist_by_residue[residue], minimum)
    rows = []
    for residue, frames in sorted(frames_by_residue.items(), key=lambda item: (item[0][0], item[0][1], item[0][2])):
        rows.append({
            "neighbor_segid": residue[0],
            "neighbor_resindex": residue[1],
            "neighbor_resid": residue[2],
            "neighbor_resname": residue[3],
            "frames_in_contact": len(frames),
            "contact_frame_fraction": (len(frames) / n_frames if n_frames else float("nan")),
            "minimum_distance_A": mindist_by_residue[residue],
        })
    return rows


def residue_pair_contact_frequency(contact_rows: list[dict[str, Any]], n_frames: int | None = None) -> list[dict[str, Any]]:
    """Aggregate atom contacts into residue-pair contact frequencies.

    Each residue pair contributes at most once per frame, even if multiple
    atom pairs are within the cutoff. The frame fraction denominator includes
    zero-contact frames when ``n_frames`` is supplied.
    """
    by_frame_pair: dict[tuple[int, tuple[Any, ...]], float] = {}
    for row in contact_rows:
        pair = (
            row["center_segid"],
            int(row.get("center_resindex", row["center_resid"])),
            int(row["center_resid"]), row["center_resname"],
            row["neighbor_segid"],
            int(row.get("neighbor_resindex", row["neighbor_resid"])),
            int(row["neighbor_resid"]), row["neighbor_resname"],
        )
        key = (int(row["frame"]), pair)
        by_frame_pair[key] = min(float(row["distance_A"]), by_frame_pair.get(key, float("inf")))
    frames_by_pair: dict[tuple[Any, ...], set[int]] = defaultdict(set)
    min_by_pair: dict[tuple[Any, ...], float] = defaultdict(lambda: float("inf"))
    for (frame, pair), minimum in by_frame_pair.items():
        frames_by_pair[pair].add(frame)
        min_by_pair[pair] = min(min_by_pair[pair], minimum)
    rows = []
    for pair, frames in sorted(frames_by_pair.items(), key=lambda item: item[0]):
        rows.append({
            "center_segid": pair[0], "center_resindex": pair[1],
            "center_resid": pair[2], "center_resname": pair[3],
            "neighbor_segid": pair[4], "neighbor_resindex": pair[5],
            "neighbor_resid": pair[6], "neighbor_resname": pair[7],
            "frames_in_contact": len(frames),
            "contact_frame_fraction": (len(frames) / n_frames if n_frames else float("nan")),
            "minimum_distance_A": min_by_pair[pair],
        })
    return rows



def iter_interfacial_contacts(
    universe: Any,
    selection_a: str = "protein and segid A and not name H*",
    selection_b: str = "protein and segid B and not name H*",
    cutoff_angstrom: float = 4.5,
    *,
    trajectory_id: str = "trajectory",
    segment_id: str | int = 0,
    time_scale: float = 1.0,
    time_offset: float = 0.0,
    time_unit: str = "ps",
) -> Iterator[dict[str, Any]]:
    """Yield interfacial heavy-atom contacts without storing a full trajectory's output.

    Unlike :func:`neighbor_contacts`, which returns all rows and per-frame summaries,
    this iterator yields only atom pairs at or below the specified cutoff. Distances and
    cutoffs are in Å. Each row carries segment/frame provenance and whether PBC was applied.
    Both selections must be non-empty; missing selections raise instead of returning no rows.
    """
    from MDAnalysis.lib.distances import capped_distance

    if not np.isfinite(cutoff_angstrom) or cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be a finite positive number")
    if not np.isfinite(time_scale) or not np.isfinite(time_offset):
        raise ValueError("time_scale and time_offset must be finite numbers")
    if not isinstance(time_unit, str) or not time_unit.strip():
        raise ValueError("time_unit must be a non-empty string")

    group_a = universe.select_atoms(selection_a)
    group_b = universe.select_atoms(selection_b)
    if len(group_a) == 0:
        raise ValueError(f"selection_a matched no atoms: {selection_a!r}")
    if len(group_b) == 0:
        raise ValueError(f"selection_b matched no atoms: {selection_b!r}")

    def generate() -> Iterator[dict[str, Any]]:
        for ts in universe.trajectory:
            dimensions = ts.dimensions if valid_box_dimensions(ts.dimensions) else None
            pairs, distances = capped_distance(
                group_a.positions,
                group_b.positions,
                max_cutoff=float(cutoff_angstrom),
                box=dimensions,
                return_distances=True,
            )
            for (idx_a, idx_b), distance in zip(pairs, distances, strict=True):
                atom_a = group_a[int(idx_a)]
                atom_b = group_b[int(idx_b)]
                yield {
                    "trajectory_id": str(trajectory_id),
                    "segment_id": str(segment_id),
                    "local_frame": int(ts.frame),
                    "time": scaled_time(ts, time_scale, time_offset),
                    "time_unit": time_unit,
                    "pbc_applied": dimensions is not None,
                    "atom_a_index": int(atom_a.index),
                    "atom_a_segid": str(getattr(atom_a, "segid", "")),
                    "atom_a_resindex": int(atom_a.residue.ix),
                    "atom_a_resid": int(atom_a.resid),
                    "atom_a_resname": str(atom_a.resname),
                    "atom_a_name": str(atom_a.name),
                    "atom_b_index": int(atom_b.index),
                    "atom_b_segid": str(getattr(atom_b, "segid", "")),
                    "atom_b_resindex": int(atom_b.residue.ix),
                    "atom_b_resid": int(atom_b.resid),
                    "atom_b_resname": str(atom_b.resname),
                    "atom_b_name": str(atom_b.name),
                    "distance_A": float(distance),
                }

    return generate()
