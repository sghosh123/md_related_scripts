"""Regression tests for inter-protein geometry and output provenance."""
import numpy as np
import pytest
import MDAnalysis as mda

from copper_trafficking.geometry import angle_degrees, distance_angstrom_pbc, minimum_image_vector
from copper_trafficking.interface_features import (
    DEFAULT_INTERFACE_SELECTIONS,
    FEATURE_COLUMNS,
    compute_frame_interface_features,
    iter_interface_features,
)


def make_feature_universe():
    universe = mda.Universe.empty(
        n_atoms=8,
        n_residues=8,
        n_segments=2,
        atom_resindex=np.arange(8),
        residue_segindex=np.asarray([0, 0, 0, 0, 1, 1, 1, 1]),
        trajectory=True,
    )
    universe.add_TopologyAttr("names", ["CA"] * 8)
    universe.add_TopologyAttr("types", ["C"] * 8)
    universe.add_TopologyAttr("resids", [56, 63, 18, 26, 55, 62, 18, 25])
    universe.add_TopologyAttr("resnames", ["ALA"] * 8)
    universe.add_TopologyAttr("segids", ["A", "B"])
    universe.atoms.positions = np.array(
        [
            [1.0, 1.0, 1.0],  # A1
            [2.0, 1.0, 1.0],  # A2
            [1.0, 2.0, 1.0],  # B1
            [2.0, 2.0, 1.0],  # B2
            [1.4, 1.0, 1.0],  # X1
            [2.4, 1.0, 1.0],  # X2
            [1.0, 2.4, 1.0],  # Y1
            [2.0, 2.4, 1.0],  # Y2
        ],
        dtype=np.float32,
    )
    universe.dimensions = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    return universe


def test_pbc_distance_and_displacement_are_minimum_image_corrected():
    box = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
    assert distance_angstrom_pbc([9.8, 1, 1], [0.2, 1, 1], box) == pytest.approx(0.4, abs=1e-5)
    assert minimum_image_vector([-9.0, 0.0, 0.0], box) == pytest.approx([1.0, 0.0, 0.0])


def test_angle_degrees_has_expected_geometry_and_rejects_zero_vectors():
    assert angle_degrees([1, 0, 0], [0, 1, 0]) == pytest.approx(90.0)
    with pytest.raises(ValueError, match="zero-length"):
        angle_degrees([0, 0, 0], [1, 0, 0])


def test_feature_vector_has_8_distances_and_4_angles_with_explicit_units():
    universe = make_feature_universe()
    features = compute_frame_interface_features(universe)
    assert len(features) == 12
    assert set(features) == set(FEATURE_COLUMNS)
    assert features["angle_AX_deg"] == pytest.approx(0.0)
    assert features["d_A1X1_A"] == pytest.approx(0.4, abs=1e-5)
    assert features["d_B2Y2_A"] == pytest.approx(0.4, abs=1e-5)
    assert all(value >= 0 for value in features.values())


def test_feature_selection_fails_loudly_if_atom_missing_or_ambiguous():
    universe = make_feature_universe()
    bad = dict(DEFAULT_INTERFACE_SELECTIONS)
    bad["A1"] = "segid A and resid 999 and name CA"
    with pytest.raises(ValueError, match="resolved to 0 atoms"):
        compute_frame_interface_features(universe, bad)

    ambiguous = dict(DEFAULT_INTERFACE_SELECTIONS)
    ambiguous["A1"] = "segid A and name CA"
    with pytest.raises(ValueError, match="expected exactly 1"):
        compute_frame_interface_features(universe, ambiguous)


def test_feature_iterator_records_trajectory_and_segment_provenance():
    universe = make_feature_universe()
    rows = list(iter_interface_features(
        universe,
        trajectory_id="demo-T7",
        segment_id="part2",
        time_scale=2.0,
        time_offset=-10.0,
        time_unit="ps",
    ))
    assert len(rows) == 1
    row = rows[0]
    assert row["trajectory_id"] == "demo-T7"
    assert row["segment_id"] == "part2"
    assert row["local_frame"] == 0
    assert row["time"] == pytest.approx(-10.0)
    assert row["time_unit"] == "ps"
    assert row["time_scale"] == pytest.approx(2.0)
    assert row["pbc_applied"] is True


def test_time_conversion_rejects_nonfinite_values():
    with pytest.raises(ValueError, match="finite"):
        iter_interface_features(make_feature_universe(), time_scale=float("nan"))
