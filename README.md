# Copper Trafficking Mechanism: Molecular Simulation Workflows

A reproducible Python toolkit for analyzing Cu(I) coordination and protein dynamics in copper-trafficking systems, including chaperone–metal-binding-domain (MBD) association, cysteine-centered environments, contact networks, and structural visualization.

> **Status:** development snapshot v0.3.0. The analysis utilities are tested on synthetic trajectories. Independent, configurable OpenMM restraint and NAMD equilibration/continuation templates are now included; each is explicitly marked as a template and requires system-specific input/parameter review before running.

## Included workflows

- **Atom-pair distances:** single or multiple named distances through trajectories, with per-frame periodic-boundary reporting.
- **Neighbor/contact analysis:** PBC-aware atom-pair contacts, unique neighbor atom/residue counts per frame, and residue-level contact occupancy.
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

Selections must resolve to exactly one atom. The tool writes distances in Å. MDAnalysis trajectory time is used by default; supply `--time-scale`, `--time-offset`, and `--time-unit` only when the conversion has been independently verified.

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

## Source migration

The reusable workflows were identified from reviewed archives including `holo_sim.zip`, `holo_sim_3_coord.zip`, `neighbor_test.zip`, `colvar.zip`, `restrain_test.zip`, the human/yeast chaperone/MBD snapshots, and the yeast complex archive. No research trajectory/output collection is copied here. The full migration table and known limitations are in [`docs/workflow-catalog.md`](docs/workflow-catalog.md).

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
