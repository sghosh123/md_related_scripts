# NAMD configuration templates

These are **editable, independent templates**, generalized from the historical NAMD equilibration workflows. They are not configurations for a specific protein and are deliberately blocked by a `template_acknowledgement` sentinel until the user reviews and edits the files. NAMD interprets relative input/output paths from its process working directory; either launch NAMD from the directory you intend to use as the path base or supply explicit absolute paths.

- `equilibration_nvt.conf.example`: fresh structure + coordinate PDB, minimization, velocity initialization, NVT run, optional positional restraints.
- `continue_npt.conf.example`: continue from matching NAMD restart coordinates, velocities and extended-system box (`.restart.coor`, `.restart.vel`, `.restart.xsc`) into NPT.

## Before running

1. Choose one consistent system. Verify that the PSF, starting PDB, parameter/stream files, copper/metal terms and restart files all belong together.
2. Replace all `CHANGE_ME` paths. Complete the parameter list; the shown list is illustrative, not a complete CHARMM36 distribution.
3. Set the actual unit-cell vectors/origin (NVT template) or verify the input `.xsc` file (NPT template). The NVT example's 60 Å box lengths are placeholders, not defaults for research.
4. Do not copy numerical restraint constants directly between NAMD and OpenMM: the energy expressions and coefficient conventions differ. NAMD positional restraint files store per-atom values in the selected PDB field, and the NAMD manual expresses the default quadratic potential as `k × displacement²`; OpenMM harmonic bond/angle forces use `1/2 × k × displacement²`.
5. Verify NAMD version compatibility. The templates use `oneFourScaling`, the current spelling in NAMD 3 documentation. If an older NAMD build rejects it, use that build’s documented `1-4scaling` spelling and validate the full config with the exact executable you will run.
6. For positional restraints, `consref` provides reference positions; `conskfile` provides the per-atom force constants; `conskcol B` selects the PDB B-factor field. Atom order/count must match the PSF. A zero force constant leaves an atom unrestrained. Treat `constraintScaling` as a multiplier, not as a substitute for checking the PDB constants.
7. Set `template_acknowledgement` to `MODEL_REVIEWED` only after auditing the model, parameters, geometry, PBC, timestep, output cadence, and restraint selection.
8. Run an energy/minimization smoke test and a short stability run before long production. Review energies, temperature, pressure/volume, restart continuity and restraint forces.

## Important limitations

- The NVT template is orthorhombic. Replace the three basis vectors explicitly for triclinic systems.
- The continuation template assumes the prior run wrote all three restart files with the same basename and that `restart_step` is the correct global timestep.
- These templates cannot prove a metal force field is chemically valid. Cu(I) coordination parameters must come from a documented, independently reviewed model.
- NAMD is not bundled. These configuration files have been statically reviewed but not executed in NAMD during repository integration.
