# Simulation setup audit and migration boundary

The archive review did not justify copying the original production scripts unmodified. Instead, the repository now supplies **independent, parameterized templates** that can be adapted to another system after the user validates the topology, force field, atom map, restraint geometry, units, box and protocol.

- OpenMM: [`../simulation_templates/openmm/README.md`](../simulation_templates/openmm/README.md), `run_restraint_simulation.py`, and `copper_restraints.example.json`.
- NAMD: [`../simulation_templates/namd/README.md`](../simulation_templates/namd/README.md), `equilibration_nvt.conf.example`, and `continue_npt.conf.example`.

The actual historical production files are not copied verbatim as runnable scripts. Their repeated logic is represented by the configurable replacement. This preserves useful practices while making system-specific assumptions explicit and preventing guessed particle indices or hidden legacy paths from becoming new defaults.

## OpenMM copper-restraint scripts

Reviewed source archives include `holo_sim.zip` and `restrain_test.zip`, with related output and Tcl helper folders in the human/yeast chaperone and MBD snapshots.

### Included replacement template

The new script:

- Uses modern `openmm` imports rather than the legacy `simtk.openmm` imports and Python 2 `commands` module.
- Defines inputs, force-field XML files, integrator, thermostat/barostat, nonbonded settings, output cadence and restraint lists in JSON.
- Resolves restraint atoms by chain ID, residue ID/name, atom name and optional insertion code/element; selections must resolve to exactly one topology atom.
- Writes a mapping CSV linking each selector to the actual zero-based particle index generated from the same topology.
- Supports harmonic bond restraints (`r0_nm`, `k_kj_mol_nm2`), harmonic angle restraints (`theta0_degrees`, `k_kj_mol_rad2`) and periodic torsion restraints (`periodicity`, `phase_degrees`, `k_kj_mol`).
- Provides `--check-config` without OpenMM, `--validate-model` to build a system/check mapping/evaluate initial energy without dynamics, and `--run` only after the config review checklist and explicit command-line acknowledgement are satisfied.
- Uses a finite, editable minimization iteration cap instead of copying historical minimization limits of tens of billions of iterations.

The example's numeric target values and selectors are demonstration values only. They are not endorsed for 1FD8, 1TL4, or Atx1–Ccc2. A force-field XML set must parameterize every atom/residue, including copper and its coordination environment; the restraint script does not infer or generate missing copper parameters.

### Specific inconsistencies in the legacy source

1. **1FD8 angle target conflict.** Active code contains `HarmonicAngleForce.addAngle(..., 2.28638, 460.24)`, where 2.28638 rad is approximately 131.0°. The nearby commented custom-angle block uses 2.0944 rad (120°) and a different coefficient convention. These are not equivalent targets and must not be merged without checking the intended Cu coordination geometry.
2. **Inconsistent bond-energy convention in comments.** The active `HarmonicBondForce.addBond()` calls use target lengths 0.243/0.232 nm (2.43/2.32 Å) and bare constants 2550 in OpenMM native units (kJ mol⁻¹ nm⁻²). The commented custom energy expression is `k*(r-r0)^2`, without the 1/2 factor used by `HarmonicBondForce`'s `0.5*k*(r-r0)^2`; its commented `100 kcal mol⁻¹ Å⁻²` coefficient is also numerically different after unit conversion. The coefficients cannot be interchanged directly. NAMD position restraints also use a separate per-atom PDB-field coefficient and its documented quadratic energy convention; the numerical constants are not interchangeable with OpenMM harmonic force constants.
3. **Torsion needs chemical review.** A `PeriodicTorsionForce` term in the 1FD8 script uses atom indices 84–11618–56–53, periodicity 3, phase π and k=300. The replacement can express torsion terms, but a user must validate the four-atom ordering, phase, periodicity and coefficient against the intended geometry.
4. **Huge minimization count.** The 1FD8 script calls `minimizeEnergy(20000000000)`, which is not copied. A finite, measured convergence protocol is required.
5. **GPU property names are version-dependent.** Older OpenMM docs/examples use names such as `CudaDeviceIndex`/`CudaPrecision`; modern OpenMM 8 documentation uses `DeviceIndex`/`Precision`. Configure platform properties for the installed OpenMM version rather than assuming a name from the historic scripts.
6. **Reporter cadence mismatch.** The 1FD8 source has trajectory output at 10,000,000 steps but checkpoint/log cadence of 100 steps and force output at 1,000,000 steps. The replacement exposes these independently; they need to be chosen deliberately for the intended trajectory and storage budget.
7. **System-specific files and box are absent from the source snapshots.** Force-field XML, PDB files, custom copper parameters, output restart/checkpoint files and box dimensions are referenced by names/paths not included as reusable, self-contained simulation inputs.

### Required validation before a production run

1. Record the source and checksum for PDB/topology and each force-field XML; confirm the license allows sharing.
2. Check residue templates and charges for copper, ligands, water, ions and termini. Do not silence unmatched-residue errors without resolving the parameterization.
3. Review the CSV atom map from `--validate-model`. Confirm unique selections match the chemically intended atoms and their topology order.
4. Verify every target geometry against experimental/validated model evidence and report units. Examine initial potential energy and forces; run an energy/minimization smoke test before production.
5. Verify box dimensions/cutoff compatibility, constraints, integrator timestep, thermodynamic ensemble, output intervals, checkpoint portability and the NPT barostat protocol.
6. Use a short run to inspect potential/kinetic energy, temperature, pressure/volume if relevant, restraint energy contribution if instrumented, trajectory continuity and whether restraints overconstrain the model.

OpenMM was not installed in the environment used for this integration pass, and the required copper force-field/input files were not available; the template's schema can be tested, but the model-specific preflight and simulation have not been executed here.

## NAMD equilibration configurations

Reviewed source archives include `y-complex-apo-md.zip`, `yeast_chaperone_apo.zip`, `yeast_MBD_apo.zip`, and `Wilson_MBD1_apo.zip` through `Wilson_MBD3_apo.zip`.

### Included replacement templates

- `equilibration_nvt.conf.example` supports a fresh PSF/PDB, a user-maintained CHARMM parameter list, orthorhombic cell lengths/origin, minimization, velocity initialization, an NVT run and optional positional restraints.
- `continue_npt.conf.example` expects consistent `.restart.coor`, `.restart.vel` and `.restart.xsc` from the preceding stage, restores the extended-system cell and uses a Langevin piston to run NPT.
- Both require paths and an explicit `template_acknowledgement` edit. Paths containing `CHANGE_ME` deliberately fail early. Restraint inputs remain optional and are checked only when enabled.

These templates are generic examples, not a byte-for-byte copy of any specific archive's setup. The NVT cell lengths (60 Å) and numerical protocol choices are placeholders that must be set to the actual system. The `.example` suffix is intentional: validate using the NAMD version installed on the user's HPC system.

### Specific inconsistencies in legacy inputs

1. **Mixed system provenance.** Root-level configs in `y-complex-apo-md.zip` refer to `apo_bacillus_chaperone_equilibriation.psf/.pdb`, while `mmz_nvt` configs target a yeast Atx1–Cu–Ccc2 complex. These must remain separate until proven to belong to one workflow.
2. **Cross-system duplicate config.** `step_5_NPT.conf` is byte-identical across `y-complex-apo-md.zip` and the yeast MBD snapshot despite different intended systems. The new NPT template uses a clear restart contract rather than carrying the duplicated system-specific values forward.
3. **Implicit/unresolved paths and Tcl variables.** Legacy configs depend on local `toppar/` installations, `$inputname`, `$temperature`, pre-existing output folders and restart filenames. The generic templates make these variables explicit, fail on unresolved file paths, and document the working-directory assumption.
4. **Restraint PDB semantics need verification.** Legacy files enable `constraints on` and use `conskcol B` with a restraint PDB. The reference coordinates and per-atom force constants are different concepts even when stored in separate files; atom order/count must match the PSF. `constraintScaling` multiplies the per-atom restraint strength and should not be used to disguise missing/incorrect B-factor values.
5. **PBC must match the input.** The old configs contain literal box vectors and origins. Copying them to a different solvated box or system is unsafe. The new NVT template marks all box values as user settings; the NPT template loads the box from the matching restart `.xsc` file.
6. **Version-sensitive keyword.** The templates use `1-4scaling`, matching many NAMD 2.x files. NAMD 3 documentation deprecates it in favor of `oneFourScaling`; users should change it to the keyword supported by the actual installed version and validate the configuration.
7. **Restart continuity is not inferred.** Set `restart_step` to the actual global step at the end of the source run and ensure the NVT base name resolves to the matching restart coordinate, velocity and extended-system files.

NAMD was not installed in the integration environment. The templates received static syntax-completeness checks, but were not executed by NAMD; the user's system files, parameters, restraint PDBs, NAMD build and HPC launcher must be validated locally.

## Third-party and generic test material

`restrain_test.zip` includes unrelated di-phenyl, solvated-methanol and CHARMM-GUI-generated files alongside copper-related tests. These were not copied. Third-party CHARMM-GUI files or force-field materials should not be bundled without version/license checks and upstream attribution.


## References for version-sensitive settings

- OpenMM platform properties and running simulations: [OpenMM User Guide — Running Simulations](https://docs.openmm.org/latest/userguide/application/02_running_sims.html). Platform property names vary by platform and OpenMM version; inspect `platform.getPropertyNames()` as the runner does.
- OpenMM force-field construction: [OpenMM `ForceField` API](https://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html).
- OpenMM custom force expressions and the harmonic energy convention: [OpenMM User Guide — Custom Forces](https://docs.openmm.org/latest/userguide/theory/03_custom_forces.html).
- NAMD positional restraint semantics (`consref`, `conskfile`, `conskcol`, `constraintScaling`): [NAMD 3.0 User Guide — Constraints and Restraints](https://www.ks.uiuc.edu/Research/namd/3.0b3/ug/node29.html).
- NAMD configuration and version-specific options: [NAMD 3.0 User Guide — Configuration File](https://www.ks.uiuc.edu/Research/namd/3.0/ug/node9.html).

These references inform the template conventions; users must still validate exact behavior against the manual for the NAMD/OpenMM version installed on their system.
