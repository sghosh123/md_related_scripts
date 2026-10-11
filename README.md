# Copper Trafficking Mechanism: Molecular Simulation Workflows

Python toolkit for analyzing Cu(I) coordination and protein dynamics in copper-trafficking systems, including chaperone–metal-binding-domain (MBD) association, cysteine-centered environments, contact networks, and structural visualization.


## Included workflows

- **Atom-pair distances:** single or multiple named distances through trajectories, with per-frame periodic-boundary reporting and optional binary proximity columns.
- **Neighbor/contact analysis:** PBC-aware atom-pair contacts, unique neighbor atom/residue counts per frame, and residue-level contact occupancy.
- **Streaming interfacial contacts:** long-form heavy-atom contacts between two explicit selections, written incrementally to avoid retaining all contacts in memory.
- **Interface feature vectors:** eight inter-protein Cα distances and four orientation angles, with validated single-atom selections, PBC-aware geometry, units in column names, and trajectory/segment provenance.
- **Hydrogen-bond geometry screen:** topology-bonded donor–hydrogen pairs, D–A distance, and D–H–A angle criteria.
- **Distribution plots:** histograms of named CSV columns or whitespace-delimited legacy data, without hardcoded row counts.
- **Torsions and binding-loop conformations:** explicitly defined four-atom dihedrals plus a phi/psi scatter plot helper.
- **Fitted RMSD:** aligned RMSD versus a chosen reference frame, guarded by a whole-molecule/PBC acknowledgment.
- **Group position / colvar support:** center-of-geometry or center-of-mass coordinates, guarded by an explicit group-whole acknowledgment.
- **Frame export and visualization:** extract selected trajectory frames with Python, with VMD style/export helpers for structural inspection.

## Repository structure

```text
.
├── configs/                       # Example JSON configurations
├── docs/                          # Scientific assumptions, archive audit and migration notes
├── examples/                      # Commands and selection guidance
├── scripts/                       # User-facing CLI tools
├── simulation_templates/openmm/   # Config-driven bond/angle/torsion-restraint runner
├── src/copper_trafficking/        # Reusable analysis functions
│   ├── geometry.py                # PBC-aware distances, displacements and angles
│   └── interface_features.py      # 8 distances + 4 orientation angles
├── tests/                         # Unit and synthetic trajectory integration tests
├── simulation_templates/namd/     # NVT and NPT templates with review gates
└── visualization/vmd/             # Generic VMD representation/frame-export helpers
```

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
pytest
```

## Quick examples

### Cu–S distance time series

```bash
python scripts/analyze_distance.py \
  --topology system.psf --trajectory trajectory.dcd \
  --selection-a 'segid PROA and resid 16 and name SG' \
  --selection-b 'resname CU1' \
  --output-csv results/cu_cys16.csv \
  --output-plot results/cu_cys16.png
```

Selections must resolve to exactly one atom. The tool writes distances in Å. MDAnalysis trajectory time is used by default; supply `--time-scale`, `--time-offset`, and `--time-unit` only when the conversion has been independently verified. The multi-distance tool can additionally write `<label>_present` columns with `--presence-cutoff 3.0`, for example; that column means distance at or below the supplied cutoff, not chemical bond formation.

### Inter-protein interface features

```bash
python scripts/analyze_interface_features.py \
  --topology system.psf --trajectory trajectory.dcd \
  --trajectory-id T7 --segment-id part1 \
  --selections-json configs/interface-features.example.json \
  --output-csv results/interface_features_T7_part1.csv
```

The bundled selections reproduce a specific yeast MBD/chaperone feature definition and must be checked against your topology. The output has four angle columns (`*_deg`) and eight distance columns (`*_A`), plus trajectory ID, segment ID, local frame, time metadata, and whether valid unit-cell dimensions enabled PBC correction. Override selections with a reviewed JSON file for other systems.

### Stream interfacial heavy-atom contacts

```bash
python scripts/analyze_interface_contacts.py \
  --topology system.psf --trajectory trajectory.dcd \
  --selection-a 'protein and segid A and not name H*' \
  --selection-b 'protein and segid B and not name H*' \
  --cutoff 4.5 --output-csv results/interfacial_contacts.csv
```

This writes one row per atom-pair contact within the cutoff and streams rows to disk. It complements `analyze_neighbors.py`, which retains per-frame summaries and residue/contact occupancy statistics.

### Cysteine-centered water contacts

```bash
python scripts/analyze_neighbors.py \
  --topology system.psf --trajectory trajectory.dcd \
  --center-selection 'segid PROT and resid 13 and name SG' \
  --neighbor-selection 'resname HOH TIP3 WAT' --cutoff 3.5 \
  --output-contacts results/water_contacts.csv \
  --output-summary results/water_neighbor_summary.csv \
  --output-residue-frequency results/water_residue_occupancy.csv
```

The contact file has one row per atom-pair contact. The per-frame summary counts unique neighbor atoms and residues and includes frames with zero contacts. See [`examples/README.md`](examples/README.md) for multi-pair distances, contact-network occupancy, hydrogen-bond geometry, group positions, dihedrals, RMSD, frame extraction, and histogram commands.

## Simulation templates (OpenMM and NAMD)

The repository now includes independent, editable simulation templates:

- [`simulation_templates/openmm/run_restraint_simulation.py`](simulation_templates/openmm/run_restraint_simulation.py) reads a JSON config and can add named-atom harmonic bond restraints, harmonic angle restraints, and periodic torsion terms. Atom selectors are resolved against the loaded topology and a CSV maps selectors to zero-based OpenMM particle indices. Distance/force/angle units are explicit in the config keys.
- [`simulation_templates/namd/equilibration_nvt.conf.example`](simulation_templates/namd/equilibration_nvt.conf.example) is a fresh-start minimization/NVT template with optional positional restraints.
- [`simulation_templates/namd/continue_npt.conf.example`](simulation_templates/namd/continue_npt.conf.example) continues from matching NAMD restart coordinates, velocities and box state.

Read each template README and [`docs/simulation-setups-audit.md`](docs/simulation-setups-audit.md) before use. The numeric targets in the OpenMM example and the NAMD cell lengths are illustrative, not validated for your system. OpenMM restraint forces are additive to the selected force field, so review whether they duplicate any existing bonded terms. NAMD relative paths are interpreted from the process working directory, not automatically from the configuration file directory. Templates intentionally block production until path/settings review acknowledgements are changed; validate the actual force field, atom mapping, geometry, units, and restart continuity before a long run. The archived legacy production files are not copied verbatim because they contain inconsistent model settings and hardcoded atom indices.

## Scientific safeguards

1. Validate topology/trajectory provenance and checksums.
2. Confirm residue numbering, segment IDs, atom names, and selection cardinality.
3. Check that trajectory cell dimensions are valid before using PBC; output reports whether PBC was applied per frame.
4. Distinguish distance-only proximity screens from geometric hydrogen-bond criteria.
5. Treat time scaling/offsets, restraint constants/targets, and cutoffs as scientific choices that must be documented.
6. Make selected groups whole before interpreting centers of geometry/mass.
7. The independent simulation templates are starting points, not chemically validated copper models; the audit records specific legacy discrepancies and unresolved provenance issues.
8. Interface-feature selections and orientation-vector definitions are system-specific. Confirm topology IDs, vector endpoints, PBC behavior, units and results on manually inspected frames before using the features for clustering or enhanced sampling.

### OpenMM copper-restraint template

```bash
python simulation_templates/openmm/run_restraint_simulation.py \
  --config simulation_templates/openmm/copper_restraints.example.json \
  --check-config
```

After copying/editing the config and installing OpenMM in the environment that will run it, use `--validate-model` first. This builds a System and reports the atom mapping and initial potential energy without integrating dynamics. Only after reviewing all model-specific terms, set `template_status` to `user_reviewed`, set every `review_checklist` flag to `true`, and call `--run --confirm-model-reviewed`. The template runner requires Python packages in the main analysis environment only when performing analysis; OpenMM is an optional dependency (`pip install -e '.[openmm-simulation]'`).

### NAMD templates

Read [`simulation_templates/namd/README.md`](simulation_templates/namd/README.md), edit the `.conf.example` files, validate them with the NAMD build installed on your system, and start with a short smoke test. NAMD is an external executable and is not installed by this Python package.

## Data and licensing

Large trajectories, restart files, unpublished datasets, and unreviewed third-party force-field materials are excluded. The repository license remains unresolved; do not make a public release until the owner chooses a license and confirms redistribution rights for any additional material.

## Citation

Until a release and citation record are created, cite the manuscript/preprint associated with the analysis and identify the repository commit used. Do not cite this toolkit alone as evidence for a mechanistic conclusion.
