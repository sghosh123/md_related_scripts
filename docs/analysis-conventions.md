# Analysis conventions and scientific guardrails

## Distances and time
- Distances are reported in Å, consistent with MDAnalysis coordinates.
- Reported time is `ts.time * time_scale + time_offset`.
- Verify timestep, trajectory write interval, and frame origin before applying time conversions.
- Distance-pair and dihedral atom selections must each resolve to exactly one atom. Neighbor/contact selections intentionally accept groups with multiple atoms and report all pair contacts within the chosen cutoff.

## Periodic boundaries
The CLI uses MDAnalysis unit-cell dimensions for minimum-image distance calculations when valid dimensions are present. `--no-pbc` disables this. Missing/invalid dimensions lead to non-PBC distances, recorded as `pbc_applied=False`. Inspect the CSV and trajectory metadata rather than assuming PBC was always applied.

## Metal coordination
A distance trace is an observable, not proof of a bond or transfer mechanism. Cu(I) coordination interpretation may require donor identity, geometry, protonation state, force-field treatment, and sampling. Any Cu–S cutoff must be motivated and reported, not silently embedded in analysis.

## Hydrogen bonds and proximity
Distance-only sulfur/hydrogen/acceptor screens should be described as proximity analysis unless donor–hydrogen–acceptor angles and chemical identities are checked. The geometry tool maps donors to hydrogens only through topology bond records; it does not guess a bond from proximity. State whether counts represent atoms, residues, or unique contacts.

## Provenance manifest
For each real system, record structure provenance and edits; topology/trajectory/restart filenames and checksums; segment/chain and residue numbering; Cu(I) coordination atoms and parameter source; engine/version, integrator, timestep, temperature, and restraints; trajectory write interval and units; PBC assumptions; analysis commands, package versions, outputs, and commit hash. See `system-manifest-template.md`.

## RMSD, dihedrals, and group positions
RMSD requires at least three selected atoms and a valid reference frame. The selected group must be consistently imaged/made whole before analysis; the tool does not reconstruct molecules across PBC boundaries. Dihedral calculations use four explicit single-atom selections and report degrees in a wrapped range. Group centers require the selected group to be whole/unwrapped first; a center calculated across opposite sides of a periodic box is not physically meaningful.
