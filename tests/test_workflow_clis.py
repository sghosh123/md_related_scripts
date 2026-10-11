"""Synthetic end-to-end tests for the migrated analysis commands."""
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter
from MDAnalysis.coordinates.PDB import PDBWriter

ROOT = Path(__file__).resolve().parents[1]


def make_system(tmp_path, names, resids, resnames, segids, positions, bonds=None):
    n = len(names)
    u = mda.Universe.empty(
        n,
        n_residues=n,
        n_segments=n,
        atom_resindex=np.arange(n),
        residue_segindex=np.arange(n),
        trajectory=True,
    )
    u.add_TopologyAttr("names", list(names))
    u.add_TopologyAttr("types", list(names))
    u.add_TopologyAttr("resids", list(resids))
    u.add_TopologyAttr("resnames", list(resnames))
    u.add_TopologyAttr("segids", list(segids))
    u.add_TopologyAttr("chainIDs", [str(s)[:1] or "A" for s in segids])
    u.add_TopologyAttr("altLocs", [""] * n)
    u.add_TopologyAttr("icodes", [""] * n)
    u.add_TopologyAttr("occupancies", [1.0] * n)
    u.add_TopologyAttr("tempfactors", [0.0] * n)
    element_names = {"SG": "S", "S": "S", "O": "O", "N": "N", "H": "H", "CA": "C", "CB": "C", "CU1": "Cu"}
    u.add_TopologyAttr("elements", [element_names.get(name.upper(), name[0].upper()) for name in names])
    u.add_TopologyAttr("record_types", ["ATOM"] * n)
    u.add_TopologyAttr("formalcharges", [0] * n)
    if bonds:
        u.add_TopologyAttr("bonds", bonds)
    u.atoms.positions = np.asarray(positions[0], dtype=float)
    u.dimensions = np.array([10, 10, 10, 90, 90, 90], dtype=float)
    topology = tmp_path / "system.pdb"
    trajectory = tmp_path / "trajectory.dcd"
    with PDBWriter(str(topology)) as writer:
        writer.write(u.atoms)
    with DCDWriter(str(trajectory), n_atoms=n) as writer:
        for frame_positions in positions:
            u.atoms.positions = np.asarray(frame_positions, dtype=float)
            u.dimensions = np.array([10, 10, 10, 90, 90, 90], dtype=float)
            writer.write(u.atoms)
    return topology, trajectory


def run_cli(script, *args):
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
        check=True, capture_output=True, text=True,
    )


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_interface_feature_cli_writes_twelve_features_and_provenance(tmp_path):
    topology, trajectory = make_system(
        tmp_path,
        names=["CA"] * 8,
        resids=[56, 63, 18, 26, 55, 62, 18, 25],
        resnames=["ALA"] * 8,
        segids=["A", "A", "A", "A", "B", "B", "B", "B"],
        positions=[[[1, 1, 1], [2, 1, 1], [1, 2, 1], [2, 2, 1],
                    [1.4, 1, 1], [2.4, 1, 1], [1, 2.4, 1], [2, 2.4, 1]]],
    )
    output = tmp_path / "features.csv"
    run_cli(
        "analyze_interface_features.py", "--topology", topology, "--trajectory", trajectory,
        "--trajectory-id", "T7", "--segment-id", "segment-02", "--output-csv", output,
    )
    rows = read_csv(output)
    assert len(rows) == 1
    assert rows[0]["trajectory_id"] == "T7"
    assert rows[0]["segment_id"] == "segment-02"
    assert rows[0]["local_frame"] == "0"
    assert rows[0]["pbc_applied"] == "True"
    assert float(rows[0]["d_A1X1_A"]) == pytest.approx(0.4, abs=1e-4)
    assert float(rows[0]["d_B2Y2_A"]) == pytest.approx(0.4, abs=1e-4)
    feature_names = [name for name in rows[0] if name.endswith("_A") or name.endswith("_deg")]
    assert len(feature_names) == 12


def test_streaming_interface_contact_cli_writes_contacts(tmp_path):
    topology, trajectory = make_system(
        tmp_path,
        names=["CA", "CA"], resids=[1, 2], resnames=["ALA", "ALA"], segids=["A", "B"],
        positions=[[[0.5, 0.5, 0.5], [9.5, 0.5, 0.5]], [[0.5, 0.5, 0.5], [5.5, 0.5, 0.5]]],
    )
    output = tmp_path / "interface_contacts.csv"
    manifest = tmp_path / "frames.csv"
    run_cli(
        "analyze_interface_contacts.py", "--topology", topology, "--trajectory", trajectory,
        "--selection-a", "segid A and name CA", "--selection-b", "segid B and name CA",
        "--cutoff", "1.1", "--trajectory-id", "T7", "--segment-id", "part1",
        "--output-csv", output, "--output-frame-manifest", manifest,
    )
    rows = read_csv(output)
    assert len(rows) == 1
    assert int(rows[0]["local_frame"]) == 0
    assert float(rows[0]["distance_A"]) == pytest.approx(1.0, abs=1e-4)
    assert rows[0]["pbc_applied"] == "True"
    assert rows[0]["atom_a_resindex"] == "0"
    manifest_rows = read_csv(manifest)
    assert len(manifest_rows) == 2
    assert [row["local_frame"] for row in manifest_rows] == ["0", "1"]
    assert all(row["trajectory_id"] == "T7" and row["segment_id"] == "part1" for row in manifest_rows)


def test_neighbor_cli_exports_contacts_and_pair_occupancy(tmp_path):
    topology, trajectory = make_system(
        tmp_path,
        names=["SG", "O"], resids=[13, 101], resnames=["CYS", "HOH"], segids=["A", "W"],
        positions=[[[1, 1, 1], [9, 1, 1]], [[1, 1, 1], [5, 1, 1]]],
    )
    contacts, summary, pairs = tmp_path / "contacts.csv", tmp_path / "summary.csv", tmp_path / "pairs.csv"
    run_cli(
        "analyze_neighbors.py", "--topology", topology, "--trajectory", trajectory,
        "--center-selection", "segid A and name SG", "--neighbor-selection", "segid W",
        "--cutoff", "3", "--output-contacts", contacts, "--output-summary", summary,
        "--output-pair-frequency", pairs,
    )
    assert len(read_csv(contacts)) == 1
    summary_rows = read_csv(summary)
    assert len(summary_rows) == 2
    assert [int(row["neighbor_atom_count"]) for row in summary_rows] == [1, 0]
    pair_rows = read_csv(pairs)
    assert len(pair_rows) == 1
    assert float(pair_rows[0]["contact_frame_fraction"]) == 0.5


def test_multi_distance_cli_reads_configuration_and_writes_named_columns(tmp_path):
    topology, trajectory = make_system(
        tmp_path,
        names=["SG", "CU1"], resids=[16, 500], resnames=["CYS", "CU1"], segids=["A", "M"],
        positions=[[[1, 1, 1], [9, 1, 1]], [[1, 1, 1], [8, 1, 1]]],
    )
    config = tmp_path / "pairs.json"
    config.write_text(json.dumps({"pairs": [{
        "label": "Cu_Cys16", "selection_a": "segid A and resid 16 and name SG",
        "selection_b": "segid M and resname CU1",
    }]}), encoding="utf-8")
    output = tmp_path / "distances.csv"
    run_cli("analyze_multi_distance.py", "--topology", topology, "--trajectory", trajectory,
            "--config", config, "--output-csv", output, "--presence-cutoff", "2.5",
            "--trajectory-id", "T26", "--segment-id", "part2")
    rows = read_csv(output)
    assert len(rows) == 2
    assert "Cu_Cys16_distance_A" in rows[0]
    assert [row["Cu_Cys16_present"] for row in rows] == ["1", "0"]
    assert np.isclose(float(rows[0]["Cu_Cys16_distance_A"]), 2.0, atol=1e-4)
    assert np.isclose(float(rows[1]["Cu_Cys16_distance_A"]), 3.0, atol=1e-4)
    assert [row["trajectory_id"] for row in rows] == ["T26", "T26"]
    assert [row["segment_id"] for row in rows] == ["part2", "part2"]


def test_hbond_cli_uses_pdb_bond_records_and_angle_criterion(tmp_path):
    topology, trajectory = make_system(
        tmp_path,
        names=["N", "H", "O"], resids=[1, 1, 2], resnames=["DON", "DON", "ACC"], segids=["A", "A", "B"],
        positions=[[[0, 0, 0], [1, 0, 0], [2.8, 0, 0]]], bonds=[(0, 1)],
    )
    contacts, summary = tmp_path / "hbonds.csv", tmp_path / "hbond_summary.csv"
    run_cli(
        "analyze_hbond_geometry.py", "--topology", topology, "--trajectory", trajectory,
        "--donor-selection", "segid A and name N", "--hydrogen-selection", "segid A and name H",
        "--acceptor-selection", "segid B and name O", "--output-contacts", contacts,
        "--output-summary", summary,
    )
    rows = read_csv(contacts)
    assert len(rows) == 1
    assert float(rows[0]["DHA_angle_deg"]) > 179
    assert int(read_csv(summary)[0]["hydrogen_bond_geometry_count"]) == 1


def test_frame_extraction_and_group_position_cli(tmp_path):
    topology, trajectory = make_system(
        tmp_path,
        names=["CA", "CB"], resids=[13, 13], resnames=["CYS", "CYS"], segids=["A", "A"],
        positions=[[[1, 1, 1], [2, 1, 1]], [[2, 1, 1], [3, 1, 1]]],
    )
    outdir = tmp_path / "frames"
    run_cli("extract_frames.py", "--topology", topology, "--trajectory", trajectory,
            "--frames", "0-1", "--selection", "all", "--output-dir", outdir)
    assert (outdir / "frame_000000.pdb").exists()
    assert (outdir / "frame_000001.pdb").exists()
    output = tmp_path / "position.csv"
    run_cli("analyze_group_position.py", "--topology", topology, "--trajectory", trajectory,
            "--selection", "segid A", "--method", "cog", "--assume-group-whole", "--output-csv", output)
    rows = read_csv(output)
    assert len(rows) == 2
    assert float(rows[0]["x_A"]) == 1.5
    assert float(rows[1]["x_A"]) == 2.5


def test_histogram_cli_processes_all_csv_rows(tmp_path):
    data = tmp_path / "data.csv"
    data.write_text("frame,distance_A\n0,1.0\n1,2.0\n2,3.0\n", encoding="utf-8")
    output = tmp_path / "hist.png"
    run_cli("plot_distance_distribution.py", "--input", data, "--column", "distance_A", "--bins", "3", "--output", output)
    assert output.exists() and output.stat().st_size > 0


def test_rmsd_cli_aligns_frames_and_writes_output(tmp_path):
    topology, trajectory = make_system(
        tmp_path,
        names=["CA", "CA", "CA", "CA"], resids=[1, 2, 3, 4],
        resnames=["ALA"] * 4, segids=["A"] * 4,
        positions=[
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[3, 2, 1], [4, 2, 1], [3, 3, 1], [3, 2, 2]],
        ],
    )
    output = tmp_path / "rmsd.csv"
    run_cli("analyze_rmsd.py", "--topology", topology, "--trajectory", trajectory,
            "--selection", "all", "--assume-group-whole", "--output-csv", output)
    rows = read_csv(output)
    assert len(rows) == 2
    assert max(float(row["fitted_RMSD_A"]) for row in rows) < 1e-3


def test_salt_bridge_cli_writes_explicit_geometry(tmp_path):
    topology, trajectory = make_system(
        tmp_path, names=["NZ", "OD1"], resids=[24, 61], resnames=["LYS", "ASP"], segids=["A", "B"],
        positions=[[[1, 1, 1], [4, 1, 1]]],
    )
    config = tmp_path / "salt.json"
    config.write_text(json.dumps({"salt_bridges": [{"name": "Lys24_Asp61", "selection_a": "segid A and name NZ", "selection_b": "segid B and name OD1"}]}))
    output = tmp_path / "salt.csv"
    run_cli("analyze_salt_bridges.py", "--topology", topology, "--trajectory", trajectory,
            "--config", config, "--cutoff", "3.5", "--output-csv", output)
    rows = read_csv(output)
    assert len(rows) == 1
    assert rows[0]["bridge_name"] == "Lys24_Asp61"
    assert float(rows[0]["minimum_distance_A"]) == pytest.approx(3.0, abs=1e-3)
    assert rows[0]["within_cutoff"] == "True"


def test_contact_consistency_cli_writes_replicate_summary(tmp_path):
    contacts = tmp_path / "contacts.csv"
    frames = tmp_path / "frames.csv"
    pd.DataFrame([{
        "trajectory_id": "T1", "segment_id": "0", "local_frame": 0,
        "atom_a_segid": "A", "atom_a_resindex": 0, "atom_a_resid": 24, "atom_a_resname": "LYS",
        "atom_b_segid": "B", "atom_b_resindex": 1, "atom_b_resid": 61, "atom_b_resname": "ASP",
    }]).to_csv(contacts, index=False)
    pd.DataFrame([{"trajectory_id": "T1", "segment_id": "0", "local_frame": 0},
                  {"trajectory_id": "T1", "segment_id": "0", "local_frame": 1},
                  {"trajectory_id": "T2", "segment_id": "0", "local_frame": 0},
                  {"trajectory_id": "T2", "segment_id": "0", "local_frame": 1}]).to_csv(frames, index=False)
    outdir = tmp_path / "consistency"
    run_cli("analyze_contact_consistency.py", "--contacts", contacts, "--frames", frames, "--outdir", outdir)
    summary = pd.read_csv(outdir / "contact_consistency_summary.tsv", sep="\t").to_dict("records")
    assert len(summary) == 1
    assert float(summary[0]["mean_occupancy"]) == pytest.approx(0.25)
    assert int(summary[0]["n_trajectories_with_contact"]) == 1
