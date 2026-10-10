# Migration log

## v0.3.0 — independent simulation templates and audit corrections

| New/maintained component | Migration summary | Important scientific correction |
|---|---|---|
| `simulation_templates/openmm/run_restraint_simulation.py` | Adds a config-driven OpenMM template for harmonic bond/angle and periodic torsion restraints | Resolves atoms by explicit topology selectors instead of copying old particle indices; defines units per parameter; adds `--check-config`, `--validate-model`, and gated `--run` modes |
| `simulation_templates/openmm/copper_restraints.example.json` | Makes topology, force-field XMLs, integrator, box, restraints, outputs and review checklist editable | Example values are marked illustrative and production is blocked until the user confirms model review |
| `simulation_templates/namd/equilibration_nvt.conf.example` | Adds independent NVT minimization/equilibration template | Paths, cell geometry, parameters, outputs and optional positional restraints are explicit; a review sentinel prevents accidental submission of the untouched template |
| `simulation_templates/namd/continue_npt.conf.example` | Adds NPT continuation template based on prior NAMD restart coor/vel/xsc files | Requires explicit restart base and global timestep, avoiding implicit cross-system restart reuse |
| `docs/simulation-setups-audit.md` | Documents inconsistencies in the source OpenMM and NAMD runs | Records conflicting 1FD8 angle targets and bond coefficient conventions, excessive minimization, GPU property-name changes, mixed-system NAMD paths, reused NPT config and restraint-PDB assumptions |
| Existing trajectory-analysis tools | Retained without API changes | Continue to use configurable selections and synthetic integration tests; real research trajectories still need validation |

## v0.2.0 — consolidated workflow expansion

The initial `v0.1.0` scaffold exposed only one exact atom-pair distance command. This version added generalized replacements for the reusable analytical patterns reviewed across copper-trafficking/chaperone/MBD archives.

| New/maintained component | Migration summary | Important scientific correction |
|---|---|---|
| `analyze_distance.py` | Generalized repeated Cu–S and SG–NZ distance-versus-time scripts from four yeast-intermediate orientations | Atom selections must match one atom; historical time scaling/offset remains explicit rather than assumed |
| `analyze_multi_distance.py` | Supports several named distances in a single pass | Selection definitions are config-driven; duplicate labels and ambiguous selections are rejected |
| `analyze_neighbors.py` / `contacts.py` | Replaces hardcoded neighbor/water loops and consolidates atom, residue, and residue-pair occupancy | Center and candidate selections are caller-controlled; unique counts are explicit; valid PBC is applied and recorded |
| `analyze_hbond_geometry.py` / `hydrogen_bonds.py` | Replaces distance-only “H-bond” helpers with an operational D–H···A geometry screen | Donor–hydrogen links must be in topology bond records; both D–A distance and D–H–A angle cutoffs are applied |
| `plot_distance_distribution.py` | Replaces many fixed-filename, fixed-row-count histograms | Reads all numeric observations and validates the chosen column |
| `analyze_dihedrals.py` / `plot_ramachandran.py` | Generalizes binding-loop phi/psi/Tcl torsion analysis to explicit four-atom selections | Outputs degrees and requires atom order/residue continuity to be verified |
| `analyze_rmsd.py` | Consolidates alignment/RMSD-style trajectory comparison | Requires an explicit whole-group/PBC acknowledgment and a defined fit selection/reference frame |
| `analyze_group_position.py` | Replaces hardcoded VMD atom-coordinate/centroid exports for colvar design | Requires an explicit acknowledgment that the selected group is whole/unwrapped |
| `extract_frames.py` and VMD helpers | Generalizes PDB-per-frame extraction and common styling operations | Frame ranges and selections are checked; topology-specific visualization selections remain editable |

## Validation status

For v0.3.0, validation includes schema checks for the OpenMM config, a gate test that production cannot run from the untouched example, static Tcl completeness and mocked-command control-flow checks on the NAMD templates (including optional restraints), full unit/integration tests for the existing analysis package, source syntax compilation, and package-build verification. The final suite currently contains 50 passing tests. No actual copper force field/trajectory was available; OpenMM/NAMD simulation dynamics have not been run in the integration environment.
