# Configurable OpenMM copper-restraint template

`run_restraint_simulation.py` is a new, independent replacement for the repeated historical OpenMM scripts. It supports optional harmonic bond, harmonic angle and periodic torsion restraints defined by explicit atom selectors in JSON. It does not assume the original atom indices or force-field XML filenames.

## Quick start

```bash
python simulation_templates/openmm/run_restraint_simulation.py \
  --config simulation_templates/openmm/copper_restraints.example.json \
  --check-config
```

Copy the example config, replace all input paths and selectors, and install OpenMM in the simulation environment (`pip install -e '.[openmm-simulation]'` for this repository, or install a compatible OpenMM build independently). Then run a preflight:

```bash
python simulation_templates/openmm/run_restraint_simulation.py \
  --config path/to/your-config.json --validate-model
```

The preflight creates the system on the Reference platform, resolves every restraint atom against the supplied topology, writes `restraint_atom_mapping.csv`, and reports the initial potential energy. It does **not** minimize or integrate dynamics. It can fail if the force field does not parameterize Cu/other non-standard residues; this is expected and must be resolved in the force-field model, not by bypassing the check.

After reviewing the preflight mapping/energy and force-field/model provenance, set `template_status` to `user_reviewed`, change all four `review_checklist` flags to `true`, then run:

```bash
python simulation_templates/openmm/run_restraint_simulation.py \
  --config path/to/your-config.json --run --confirm-model-reviewed
```

The runner refuses production dynamics without both the config review checklist and the explicit command-line confirmation.

## Restraint unit convention

| JSON key | Unit / convention |
|---|---|
| `r0_nm` | Equilibrium distance in nm |
| `k_kj_mol_nm2` | Harmonic-bond force constant in kJ mol⁻¹ nm⁻²; OpenMM `HarmonicBondForce` uses `E = 1/2 k (r-r0)^2` |
| `theta0_degrees` | Equilibrium angle in degrees (converted to radians for OpenMM) |
| `k_kj_mol_rad2` | Harmonic-angle force constant in kJ mol⁻¹ rad⁻²; `E = 1/2 k (theta-theta0)^2` |
| `periodicity`, `phase_degrees`, `k_kj_mol` | Periodic-torsion periodicity, phase angle and force constant in radians after conversion / kJ mol⁻¹, following OpenMM's periodic torsion convention |

Do not blindly copy restraint constants from other engines or the legacy scripts. Conversion depends on the exact energy expression: a custom expression `k*(r-r0)^2` does not use the same numerical convention as `0.5*k*(r-r0)^2` unless the coefficient is adjusted.

## Atom mapping and topology

Each atom selector can constrain `chain_id`, `residue_id`, `residue_name`, `atom_name`, `insertion_code`, and/or `element`; each selector must be specific enough to resolve exactly one atom. The generated mapping report records the actual zero-based OpenMM particle index for traceability. Check chain IDs and residue IDs in the actual PDB because renumbering, insertion codes, blank chain IDs, alternate locations and metal-residue naming vary among systems.

The runner reads the supplied PDB as-is: it does not add missing hydrogens, rebuild loops, solvate, add ions, or repair coordinates. The force-field XML list must correctly cover all atoms, including Cu(I), ligands, termini, water and ions. This script does not invent missing metal charges or Lennard-Jones parameters, and the restraint terms are not a replacement for a validated force field. Box vectors must match the coordinates. In periodic simulations, check that the nonbonded cutoff is compatible with the cell dimensions.

The configured restraint terms are **additional forces** layered on top of the force-field energy; they can double-count or conflict with native bonded terms if the same interaction is already parameterized. Decide deliberately whether each term is intended as a physical force-field term or an extra restraint, and compare the resulting energy/forces with a no-restraint baseline. The runner starts from a PDB coordinate structure and initializes velocities from the configured temperature; it is not a restart-continuation workflow. The example settings are illustrative only. They have not been validated against 1FD8, 1TL4 or the Atx1–Ccc2 system. OpenMM is not available in every runtime; the template's JSON schema and safety gate can be tested without it, but a model preflight requires OpenMM plus the actual force-field and structure inputs.
