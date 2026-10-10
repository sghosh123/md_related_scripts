import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_all_example_configs_are_valid_json():
    files = sorted((ROOT / "configs").glob("*.json"))
    assert files
    for path in files:
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(payload, dict)


def test_multi_distance_and_dihedral_configs_have_consistent_shapes():
    multi = json.loads((ROOT / "configs" / "multi-distance.example.json").read_text(encoding="utf-8"))
    assert multi["pairs"]
    assert len({item["label"] for item in multi["pairs"]}) == len(multi["pairs"])
    for pair in multi["pairs"]:
        assert pair["selection_a"] and pair["selection_b"]
    dihedrals = json.loads((ROOT / "configs" / "dihedral-analysis.example.json").read_text(encoding="utf-8"))
    assert dihedrals["dihedrals"]
    for item in dihedrals["dihedrals"]:
        assert len(item["selections"]) == 4
