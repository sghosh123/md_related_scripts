import numpy as np
import pytest
from copper_trafficking.distance import distance_angstrom, require_single_atom


def test_distance_angstrom_basic():
    assert distance_angstrom([0, 0, 0], [3, 4, 0]) == pytest.approx(5.0)


def test_distance_angstrom_rejects_wrong_shape():
    with pytest.raises(ValueError, match="3-element"):
        distance_angstrom([0, 1], [0, 1, 2])


def test_distance_angstrom_rejects_nonfinite_coordinates():
    with pytest.raises(ValueError, match="finite"):
        distance_angstrom([0, float("nan"), 0], [0, 0, 0])


class FakeSelection:
    def __init__(self, count): self.count = count
    def __len__(self): return self.count


def test_require_single_atom_accepts_one():
    require_single_atom(FakeSelection(1), "selection")


@pytest.mark.parametrize("count", [0, 2, 5])
def test_require_single_atom_rejects_ambiguous_selection(count):
    with pytest.raises(ValueError, match="exactly one atom"):
        require_single_atom(FakeSelection(count), "selection")
