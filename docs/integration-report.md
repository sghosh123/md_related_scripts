# Repository integration report — v0.3.0 development snapshot

## Outcome

The initial distance-analysis scaffold has been expanded into a modular trajectory-analysis toolkit and now includes generic simulation templates. The migration focuses on reusing analytical behavior rather than bulk-copying historical source/output trees. Historical scripts often hardcoded paths, selections, frame counts, atom indices, and time conversions; maintained tools expose these as explicit inputs instead.

## Integrated analysis commands

| Command | Purpose | Main outputs |
|---|---|---|
| `scripts/analyze_distance.py` | One exact atom pair over a trajectory | CSV time series; optional plot |
| `scripts/analyze_multi_distance.py` | Several named exact atom pairs in one pass | Named CSV columns; optional multi-series plot |
| `scripts/analyze_neighbors.py` | PBC-aware center/candidate atom contacts | Long-form contact CSV, per-frame summary, optional residue and residue-pair occupancy tables, optional plot |
| `scripts/analyze_hbond_geometry.py` | Topology-bonded D–H pair to acceptor geometric screen | Contact CSV with D–A distances and D–H–A angles; per-frame counts |
| `scripts/analyze_dihedrals.py` | Named four-atom torsions from explicit selections | CSV angle series in degrees |
| `scripts/plot_ramachandran.py` | Phi/psi scatter from torsion CSV | Static plot |
| `scripts/analyze_rmsd.py` | Fitted RMSD to a reference frame | CSV RMSD series; optional plot |
| `scripts/analyze_group_position.py` | Center of geometry or mass | CSV x/y/z series |
| `scripts/extract_frames.py` | Export selected frames to individual PDB files | Frame-indexed PDB files |
| `scripts/plot_distance_distribution.py` | Histogram from a named CSV column or numeric whitespace column | Static histogram |

Generic VMD representation and selected-frame export helpers are in `visualization/vmd/`.

## Validation completed

- **50 automated tests pass**, including the original synthetic analysis/integration tests and new OpenMM-config/NAMD-template checks.
- OpenMM example JSON passes an independent schema/numeric check without importing OpenMM. Tests verify negative force constants, ignored switch distances, unbounded minimization settings and an untouched-template production gate are rejected.
- NAMD `.conf.example` files pass Tcl structural-completeness checks. The test suite also runs them with a no-op Tcl shim and temporary placeholder inputs, covers optional restraints on/off and rejects invalid restraint flags, and checks that no Tcl errors are emitted. This catches variable-expansion/control-flow problems without pretending that NAMD dynamics were executed.
- All analysis CLI `--help` checks and Python syntax compilation pass; the existing VMD helper Tcl files also pass structural-completeness checks.
- Package metadata is bumped to `0.3.0`; the wheel build and ZIP integrity are checked as part of release packaging.

These checks do not establish scientific validity for a research system. The original PSF/PDB/trajectory/restart files and custom copper force-field XMLs are not included, so real atom selections, force-field compatibility and production outputs remain unverified.

## Simulation templates and remaining boundary

Two independent OpenMM/NAMD template families were added. The OpenMM template resolves atoms by identity, supports harmonic bond/angle and periodic torsion restraints, writes a particle-index mapping and separates schema checking, model preflight and gated production. The NAMD examples cover fresh NVT equilibration and NPT continuation with explicit parameters and restart contracts. Both families include review gates and detailed usage notes.

The original production scripts are still not copied verbatim. Their specific inconsistencies are documented in `simulation-setups-audit.md`. OpenMM and NAMD are unavailable in this integration environment, and the source copper force-field/structure files are missing; therefore the simulation templates have not been run against an actual model.

## Release boundary

This is a development snapshot. Before public publication, pick a license, add a citation/metadata record, verify that all migrated source patterns can be redistributed, and run representative commands on the actual research systems with documented input checksums and expected observables.
