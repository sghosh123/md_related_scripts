import numpy as np
import pytest
import MDAnalysis as mda

from copper_trafficking.contacts import neighbor_contacts, residue_contact_frequency, residue_pair_contact_frequency
from copper_trafficking.hydrogen_bonds import hydrogen_bond_contacts
from copper_trafficking.trajectory import valid_box_dimensions


def _two_atom_neighbor_universe():
    u = mda.Universe.empty(
        2, n_residues=2, n_segments=2,
        atom_resindex=np.array([0, 1]), residue_segindex=np.array([0, 1]), trajectory=True,
    )
    u.add_TopologyAttr("names", ["SG", "O"])
    u.add_TopologyAttr("types", ["S", "O"])
    u.add_TopologyAttr("resids", [13, 101])
    u.add_TopologyAttr("resnames", ["CYS", "HOH"])
    u.add_TopologyAttr("segids", ["A", "W"])
    u.atoms.positions = np.array([[0.5, 0.5, 0.5], [9.5, 0.5, 0.5]], dtype=float)
    u.dimensions = np.array([10, 10, 10, 90, 90, 90], dtype=float)
    return u


def _hbond_universe(angle_coordinates=None):
    u = mda.Universe.empty(
        3, n_residues=3, n_segments=3,
        atom_resindex=np.array([0, 1, 2]), residue_segindex=np.array([0, 1, 2]), trajectory=True,
    )
    u.add_TopologyAttr("names", ["N", "H", "O"])
    u.add_TopologyAttr("types", ["N", "H", "O"])
    u.add_TopologyAttr("resids", [1, 1, 2])
    u.add_TopologyAttr("resnames", ["DON", "DON", "ACC"])
    u.add_TopologyAttr("segids", ["A", "A", "B"])
    u.add_TopologyAttr("bonds", [(0, 1)])
    coordinates = angle_coordinates or [[0, 0, 0], [1, 0, 0], [2.8, 0, 0]]
    u.atoms.positions = np.asarray(coordinates, dtype=float)
    u.dimensions = np.array([10, 10, 10, 90, 90, 90], dtype=float)
    return u


def test_valid_box_dimensions_checks_lengths_and_angles():
    assert valid_box_dimensions([10, 10, 10, 90, 90, 90])
    assert not valid_box_dimensions([10, 10, 0, 90, 90, 90])
    assert not valid_box_dimensions([10, 10, 10, 0, 90, 90])
    assert not valid_box_dimensions([10, 10, 10, 90, np.nan, 90])


def test_neighbor_contacts_apply_pbc_and_write_frame_summary():
    u = _two_atom_neighbor_universe()
    contacts, summary = neighbor_contacts(
        u, "segid A and name SG", "segid W", 1.1, use_pbc=True,
    )
    assert len(contacts) == 1
    assert contacts[0]["distance_A"] == pytest.approx(1.0, abs=1e-4)
    assert contacts[0]["pbc_applied"] is True
    assert summary == [{
        "frame": 0, "time": 0.0, "neighbor_atom_count": 1,
        "neighbor_residue_count": 1, "contact_pair_count": 1, "pbc_applied": True,
    }]


def test_neighbor_contacts_without_pbc_do_not_wrap():
    contacts, summary = neighbor_contacts(
        _two_atom_neighbor_universe(), "segid A and name SG", "segid W", 1.1, use_pbc=False,
    )
    assert contacts == []
    assert summary[0]["neighbor_atom_count"] == 0
    assert summary[0]["pbc_applied"] is False


def test_residue_frequency_counts_a_residue_once_per_frame():
    contacts, _ = neighbor_contacts(
        _two_atom_neighbor_universe(), "segid A and name SG", "segid W", 1.1,
    )
    result = residue_contact_frequency(contacts, n_frames=4)
    assert result[0]["frames_in_contact"] == 1
    assert result[0]["contact_frame_fraction"] == pytest.approx(0.25)


def test_hbond_requires_donor_hydrogen_topology_bond_and_angle():
    u = _hbond_universe()
    contacts, summary = hydrogen_bond_contacts(u, "name N", "name H", "name O")
    assert len(contacts) == 1
    assert contacts[0]["DA_distance_A"] == pytest.approx(2.8, abs=1e-4)
    assert contacts[0]["DHA_angle_deg"] == pytest.approx(180.0, abs=1e-4)
    assert summary[0]["hydrogen_bond_geometry_count"] == 1


def test_hbond_angle_cutoff_excludes_bent_geometry():
    u = _hbond_universe([[0, 0, 0], [1, 0, 0], [1, 2.5, 0]])
    contacts, summary = hydrogen_bond_contacts(
        u, "name N", "name H", "name O", distance_cutoff_angstrom=3.5, angle_cutoff_degrees=150,
    )
    assert contacts == []
    assert summary[0]["hydrogen_bond_geometry_count"] == 0


def test_hbond_fails_clearly_when_bond_info_is_missing():
    u = mda.Universe.empty(
        3, n_residues=3, n_segments=3,
        atom_resindex=np.array([0, 1, 2]), residue_segindex=np.array([0, 1, 2]), trajectory=True,
    )
    u.add_TopologyAttr("names", ["N", "H", "O"])
    u.add_TopologyAttr("resids", [1, 1, 2])
    u.add_TopologyAttr("resnames", ["DON", "DON", "ACC"])
    u.add_TopologyAttr("segids", ["A", "A", "B"])
    u.atoms.positions = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=float)
    with pytest.raises(ValueError, match="bond data are required"):
        hydrogen_bond_contacts(u, "name N", "name H", "name O")


def test_residue_pair_frequency_reports_contact_fraction():
    contacts, _ = neighbor_contacts(
        _two_atom_neighbor_universe(), "segid A and name SG", "segid W", 1.1,
    )
    result = residue_pair_contact_frequency(contacts, n_frames=2)
    assert len(result) == 1
    assert result[0]["center_resid"] == 13
    assert result[0]["neighbor_resid"] == 101
    assert result[0]["contact_frame_fraction"] == pytest.approx(0.5)


def test_dihedral_series_returns_angles_in_degrees():
    from copper_trafficking.torsions import dihedral_series
    u = mda.Universe.empty(
        4, n_residues=4, n_segments=4,
        atom_resindex=np.arange(4), residue_segindex=np.arange(4), trajectory=True,
    )
    u.add_TopologyAttr("names", ["A", "B", "C", "D"])
    u.add_TopologyAttr("resids", [1, 2, 3, 4])
    u.add_TopologyAttr("resnames", ["X", "X", "X", "X"])
    u.add_TopologyAttr("segids", ["A", "A", "A", "A"])
    u.atoms.positions = np.array([[1, 0, 0], [0, 0, 0], [0, 1, 0], [0, 1, 1]], dtype=float)
    u.dimensions = np.array([10, 10, 10, 90, 90, 90], dtype=float)
    rows = dihedral_series(u, [{"label": "torsion", "selections": ["index 0", "index 1", "index 2", "index 3"]}])
    assert len(rows) == 1
    assert -180 <= rows[0]["torsion_angle_deg"] <= 180


def test_residue_frequency_keeps_duplicate_residue_labels_distinct():
    rows = [
        {"frame": 0, "neighbor_segid": "A", "neighbor_resindex": 3,
         "neighbor_resid": 10, "neighbor_resname": "CYS", "distance_A": 2.1},
        {"frame": 0, "neighbor_segid": "A", "neighbor_resindex": 9,
         "neighbor_resid": 10, "neighbor_resname": "CYS", "distance_A": 2.2},
    ]
    result = residue_contact_frequency(rows, n_frames=1)
    assert len(result) == 2
    assert {row["neighbor_resindex"] for row in result} == {3, 9}


def test_residue_pair_frequency_keeps_duplicate_residue_labels_distinct():
    rows = [
        {"frame": 0, "center_segid": "A", "center_resindex": 1, "center_resid": 10,
         "center_resname": "CYS", "neighbor_segid": "B", "neighbor_resindex": 3,
         "neighbor_resid": 20, "neighbor_resname": "ASP", "distance_A": 2.1},
        {"frame": 0, "center_segid": "A", "center_resindex": 2, "center_resid": 10,
         "center_resname": "CYS", "neighbor_segid": "B", "neighbor_resindex": 3,
         "neighbor_resid": 20, "neighbor_resname": "ASP", "distance_A": 2.2},
    ]
    result = residue_pair_contact_frequency(rows, n_frames=1)
    assert len(result) == 2
    assert {row["center_resindex"] for row in result} == {1, 2}
