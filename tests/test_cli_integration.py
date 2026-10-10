"""End-to-end smoke test with a tiny synthetic two-atom trajectory."""
import csv
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter
from MDAnalysis.coordinates.PDB import PDBWriter

ROOT = Path(__file__).resolve().parents[1]


def test_distance_cli_uses_pbc_and_writes_csv(tmp_path):
    universe = mda.Universe.empty(
        2,
        n_residues=2,
        n_segments=2,
        atom_resindex=np.array([0, 1]),
        residue_segindex=np.array([0, 1]),
        trajectory=True,
    )
    universe.add_TopologyAttr("names", ["SG", "SG"])
    universe.add_TopologyAttr("types", ["S", "S"])
    universe.add_TopologyAttr("resids", [13, 55])
    universe.add_TopologyAttr("resnames", ["CYS", "CYS"])
    universe.add_TopologyAttr("segids", ["A", "B"])
    universe.add_TopologyAttr("chainIDs", ["A", "B"])
    universe.add_TopologyAttr("altLocs", ["", ""])
    universe.add_TopologyAttr("icodes", ["", ""])
    universe.add_TopologyAttr("occupancies", [1.0, 1.0])
    universe.add_TopologyAttr("tempfactors", [0.0, 0.0])
    universe.add_TopologyAttr("elements", ["S", "S"])
    universe.add_TopologyAttr("record_types", ["ATOM", "ATOM"])
    universe.add_TopologyAttr("formalcharges", [0, 0])
    universe.atoms.positions = np.array([[1, 1, 1], [9, 1, 1]], dtype=float)
    universe.dimensions = np.array([10, 10, 10, 90, 90, 90], dtype=float)

    topology = tmp_path / "tiny.pdb"
    trajectory = tmp_path / "tiny.dcd"
    with PDBWriter(str(topology)) as writer:
        writer.write(universe.atoms)
    with DCDWriter(str(trajectory), n_atoms=2) as writer:
        for positions in (
            [[1, 1, 1], [9, 1, 1]],
            [[1, 1, 1], [8, 1, 1]],
        ):
            universe.atoms.positions = np.asarray(positions, dtype=float)
            universe.dimensions = np.array([10, 10, 10, 90, 90, 90], dtype=float)
            writer.write(universe.atoms)

    output = tmp_path / "distances.csv"
    command = [
        sys.executable, str(ROOT / "scripts" / "analyze_distance.py"),
        "--topology", str(topology), "--trajectory", str(trajectory),
        "--selection-a", "segid A and resid 13 and name SG",
        "--selection-b", "segid B and resid 55 and name SG",
        "--output-csv", str(output),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert float(rows[0]["distance_A"]) == pytest.approx(2.0, abs=1e-4)
    assert float(rows[1]["distance_A"]) == pytest.approx(3.0, abs=1e-4)
    assert [row["pbc_applied"] for row in rows] == ["True", "True"]
