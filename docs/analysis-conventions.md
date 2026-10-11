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


## Inter-protein interface feature vectors

`interface_features.py` defines eight named Cα pair distances in Å and four orientation angles in degrees. The default key-to-selection mapping is intentionally stored in one constant and also shipped as `configs/interface-features.example.json`. Confirm all eight selections resolve to exactly one atom for the topology used. Distances and within-partner displacement vectors use minimum-image correction when a valid six-component MDAnalysis box is available; every feature row records `pbc_applied`. If box metadata is absent or invalid, geometry uses the raw coordinate differences and records `pbc_applied=False`.

The orientation vectors are A1→A2, B1→B2, X1→X2, and Y1→Y2. The four returned angles are A-vs-X, B-vs-X, A-vs-Y, and B-vs-Y. This convention is part of the feature definition; change it only with a documented scientific rationale and updated reference tests. Output also includes `trajectory_id`, `segment_id`, and `local_frame`; local frame numbers alone are not globally unique across trajectory segments. The `time` field is scaled using explicit scale/offset arguments and is labelled by `time_unit`; verify reader-native time units before analysis.

The streaming interfacial-contact iterator is designed for long contact tables. It writes only pairs within the requested cutoff and does not provide zero-contact frame summaries or residue occupancies. Use `analyze_neighbors.py` when those summaries are needed.
