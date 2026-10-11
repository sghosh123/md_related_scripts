import numpy as np
import pandas as pd
import pytest
import MDAnalysis as mda

from copper_trafficking.salt_bridges import salt_bridge_geometry
from copper_trafficking.contact_consistency import summarize_contact_consistency


def _charged_pair_universe():
    u = mda.Universe.empty(2, n_residues=2, n_segments=2,
                           atom_resindex=np.array([0, 1]), residue_segindex=np.array([0, 1]), trajectory=True)
    u.add_TopologyAttr("names", ["NZ", "OD1"])
    u.add_TopologyAttr("types", ["N", "O"])
    u.add_TopologyAttr("resids", [24, 61])
    u.add_TopologyAttr("resnames", ["LYS", "ASP"])
    u.add_TopologyAttr("segids", ["A", "B"])
    u.atoms.positions = np.array([[0, 0, 0], [3, 0, 0]], dtype=float)
    u.dimensions = np.array([10, 10, 10, 90, 90, 90], dtype=float)
    return u


def test_salt_bridge_geometry_reports_all_frames_and_minimum_distance():
    u = _charged_pair_universe()
    # Universe.empty with a memory trajectory permits moving the atoms by frame only
    # when coordinates are installed in a MemoryReader; a single frame is sufficient here.
    rows = salt_bridge_geometry(u, [{"name": "Lys24_Asp61", "selection_a": "segid A and name NZ", "selection_b": "segid B and name OD1"}], cutoff_angstrom=3.5, trajectory_id="T1")
    assert len(rows) == 1
    assert rows[0]["trajectory_id"] == "T1"
    assert rows[0]["minimum_distance_A"] == pytest.approx(3.0)
    assert rows[0]["within_cutoff"] is True


def test_salt_bridge_fails_if_selection_is_empty():
    with pytest.raises(ValueError, match="matched no atoms"):
        salt_bridge_geometry(_charged_pair_universe(), [{"name": "bad", "selection_a": "segid Z", "selection_b": "segid B"}])


def test_contact_consistency_counts_zero_contact_frames_and_zero_replicates():
    contacts = pd.DataFrame([{
        "trajectory_id": "T1", "segment_id": "0", "local_frame": 0,
        "atom_a_segid": "A", "atom_a_resindex": 0, "atom_a_resid": 24, "atom_a_resname": "LYS",
        "atom_b_segid": "B", "atom_b_resindex": 1, "atom_b_resid": 61, "atom_b_resname": "ASP",
    }])
    manifest = pd.DataFrame([
        {"trajectory_id": "T1", "segment_id": "0", "local_frame": i} for i in [0, 1]
    ] + [
        {"trajectory_id": "T2", "segment_id": "0", "local_frame": i} for i in [0, 1, 2]
    ])
    by_traj, summary = summarize_contact_consistency(contacts, manifest)
    assert len(by_traj) == 2
    fractions = by_traj.set_index("trajectory_id")["occupancy"].to_dict()
    assert fractions == {"T1": 0.5, "T2": 0.0}
    assert summary.iloc[0]["mean_occupancy"] == pytest.approx(0.25)
    assert summary.iloc[0]["n_trajectories_with_contact"] == 1


def test_contact_consistency_requires_frame_manifest_and_unique_keys():
    contacts = pd.DataFrame(columns=["trajectory_id", "segment_id", "local_frame"])
    with pytest.raises(ValueError, match="required columns"):
        summarize_contact_consistency(contacts, pd.DataFrame())
    manifest = pd.DataFrame({"trajectory_id": ["T1", "T1"], "segment_id": ["0", "0"], "local_frame": [0, 0]})
    contacts = pd.DataFrame(columns=["trajectory_id", "segment_id", "local_frame", "atom_a_segid", "atom_a_resindex", "atom_a_resid", "atom_a_resname", "atom_b_segid", "atom_b_resindex", "atom_b_resid", "atom_b_resname"])
    with pytest.raises(ValueError, match="duplicate"):
        summarize_contact_consistency(contacts, manifest)
