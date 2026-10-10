import importlib.util
from pathlib import Path
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "extract_frames.py"
SPEC = importlib.util.spec_from_file_location("extract_frames", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_parse_frames_supports_indices_and_inclusive_ranges():
    assert MODULE.parse_frames("0,3,5-7,3", 10) == [0, 3, 5, 6, 7]


@pytest.mark.parametrize("spec", ["", "5-2", "1,", "abc"])
def test_parse_frames_rejects_invalid_specs(spec):
    with pytest.raises(ValueError):
        MODULE.parse_frames(spec, 10)


def test_parse_frames_rejects_out_of_range_frames():
    with pytest.raises(ValueError, match="outside valid range"):
        MODULE.parse_frames("0,10", 10)
