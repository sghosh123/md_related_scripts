"""Replicate-aware interfacial contact consistency statistics.

Requires a frame manifest that includes every analyzed frame, including frames with no
contacts. Independent trajectory IDs are the replication unit; segment IDs only split
continuous coordinate files and are not treated as independent replicates.
"""
from __future__ import annotations

from typing import Iterable
import warnings
import numpy as np
import pandas as pd

FRAME_KEYS = ["trajectory_id", "segment_id", "local_frame"]


def _required(df: pd.DataFrame, fields: list[str], label: str) -> None:
    missing = [x for x in fields if x not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def summarize_contact_consistency(
    contacts: pd.DataFrame,
    frame_manifest: pd.DataFrame,
    *,
    minimum_trajectories: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (per-trajectory occupancy, cross-trajectory summary) tables.

    The residue-pair universe is the union of residue pairs observed at least once in
    the contact table. For every such pair, trajectories in which it never occurs are
    explicitly assigned zero occupancy. Frames with no contacts are represented by the
    required manifest and remain in each trajectory's denominator.
    """
    if minimum_trajectories < 1:
        raise ValueError("minimum_trajectories must be >= 1")
    _required(frame_manifest, FRAME_KEYS, "frame manifest")
    pair_cols = [
        "atom_a_segid", "atom_a_resindex", "atom_a_resid", "atom_a_resname",
        "atom_b_segid", "atom_b_resindex", "atom_b_resid", "atom_b_resname",
    ]
    # Older contact exports lacked topology-resindex. Permit them with a warning, but
    # label the limitations; new exports now carry residue indices.
    if all(c in contacts.columns for c in pair_cols):
        pass
    else:
        legacy = ["atom_a_segid", "atom_a_resid", "atom_a_resname", "atom_b_segid", "atom_b_resid", "atom_b_resname"]
        _required(contacts, legacy + FRAME_KEYS, "contact table")
        warnings.warn(
            "Contact file lacks atom_a_resindex/atom_b_resindex; duplicate residue labels may be conflated. Re-export contacts using the updated extractor.",
            RuntimeWarning, stacklevel=2,
        )
        contacts = contacts.copy()
        contacts["atom_a_resindex"] = pd.to_numeric(contacts["atom_a_resid"], errors="raise")
        contacts["atom_b_resindex"] = pd.to_numeric(contacts["atom_b_resid"], errors="raise")
        pair_cols = [
            "atom_a_segid", "atom_a_resindex", "atom_a_resid", "atom_a_resname",
            "atom_b_segid", "atom_b_resindex", "atom_b_resid", "atom_b_resname",
        ]
    _required(contacts, FRAME_KEYS + pair_cols, "contact table")
    if not contacts.empty and contacts[pair_cols].isna().any().any():
        raise ValueError("Contact table contains missing residue-pair identity values.")
    fm = frame_manifest[FRAME_KEYS].copy()
    ct = contacts[FRAME_KEYS + pair_cols].copy()
    for name, table in (("frame manifest", fm), ("contact table", ct)):
        if table[FRAME_KEYS].isna().any().any():
            raise ValueError(f"{name} contains missing frame identity values.")
        table["trajectory_id"] = table["trajectory_id"].astype(str)
        table["segment_id"] = table["segment_id"].astype(str)
        table["local_frame"] = pd.to_numeric(table["local_frame"], errors="raise")
        frame_values = table["local_frame"].to_numpy(dtype=float)
        if not np.isfinite(frame_values).all() or (frame_values < 0).any() or (frame_values != np.floor(frame_values)).any():
            raise ValueError(f"{name} local_frame values must be finite non-negative integers.")
        table["local_frame"] = table["local_frame"].astype(int)
    if fm.duplicated(FRAME_KEYS).any():
        raise ValueError("Frame manifest has duplicate trajectory/segment/frame keys.")
    if fm.empty:
        raise ValueError("Frame manifest is empty; zero-contact frames cannot be inferred without it.")
    if not set(map(tuple, ct[FRAME_KEYS].drop_duplicates().itertuples(index=False, name=None))).issubset(
        set(map(tuple, fm[FRAME_KEYS].itertuples(index=False, name=None)))
    ):
        raise ValueError("Contact table contains frame identities absent from the frame manifest.")
    trajectories = sorted(fm["trajectory_id"].unique().tolist())
    frame_counts = fm.groupby("trajectory_id", sort=True).size().to_dict()
    if len(trajectories) < minimum_trajectories:
        raise ValueError(f"Only {len(trajectories)} independent trajectory IDs; need {minimum_trajectories}.")
    if ct.empty:
        # No pair universe can be inferred from zero contacts alone.
        empty = pd.DataFrame(columns=pair_cols + ["trajectory_id", "n_frames", "frames_in_contact", "occupancy"])
        return empty, pd.DataFrame(columns=pair_cols + ["n_trajectories", "mean_occupancy", "std_occupancy", "sem_occupancy"])

    for c in ("atom_a_resindex", "atom_a_resid", "atom_b_resindex", "atom_b_resid"):
        ct[c] = pd.to_numeric(ct[c], errors="raise").astype(int)
    pairs = ct[pair_cols].drop_duplicates().sort_values(pair_cols, kind="stable").reset_index(drop=True)
    # Binary contact per residue pair/frame; multiple atom pairs within one pair count once.
    hit = ct[FRAME_KEYS + pair_cols].drop_duplicates()
    hit_counts = hit.groupby(["trajectory_id"] + pair_cols, sort=False).size().to_dict()
    pair_records = pairs.to_dict("records")
    occupancy_rows = []
    for pair in pair_records:
        pkey = tuple(pair[c] for c in pair_cols)
        for trajectory in trajectories:
            nframes = int(frame_counts[trajectory])
            ncontact = int(hit_counts.get((trajectory,) + pkey, 0))
            occupancy_rows.append({**pair, "trajectory_id": trajectory, "n_frames": nframes,
                                   "frames_in_contact": ncontact, "occupancy": ncontact / nframes})
    per_trajectory = pd.DataFrame(occupancy_rows)
    summary_rows = []
    for pair in pair_records:
        pkey = tuple(pair[c] for c in pair_cols)
        matching = [row for row in occupancy_rows if tuple(row[c] for c in pair_cols) == pkey]
        values = np.asarray([row["occupancy"] for row in matching], dtype=float)
        if len(values) < minimum_trajectories:
            continue
        sd = float(np.std(values, ddof=1)) if len(values) > 1 else float("nan")
        summary_rows.append({**pair, "n_trajectories": int(len(values)),
                             "n_trajectories_with_contact": int(np.count_nonzero(values > 0)),
                             "mean_occupancy": float(np.mean(values)), "std_occupancy": sd,
                             "sem_occupancy": sd / np.sqrt(len(values)) if len(values) > 1 else float("nan"),
                             "min_occupancy": float(np.min(values)), "max_occupancy": float(np.max(values))})
    return per_trajectory, pd.DataFrame(summary_rows)
