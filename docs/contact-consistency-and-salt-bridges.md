# Salt-bridge geometry and replicate-aware contact consistency

## What this adds

- `copper_trafficking.salt_bridges.salt_bridge_geometry`: a per-frame, per-definition minimum distance between two explicitly selected atom groups; includes atom identities, PBC use, local frame, trajectory and segment IDs.
- `copper_trafficking.contact_consistency.summarize_contact_consistency`: residue-pair occupancy per independent trajectory plus mean, standard deviation, standard error, extrema, and number of trajectories with contact.
- `analyze_interface_contacts.py --output-frame-manifest`: writes every analyzed frame, including frames with no contacts.
- `scripts/analyze_salt_bridges.py` and `scripts/analyze_contact_consistency.py`: command-line tools.

## Scientific interpretation

A distance below a configured threshold means only that the configured charged-group selections are close in that frame. It does not determine protonation, dielectric response, interaction free energy, or transfer reactivity. Verify protonation states and selected atoms against the actual parameterized topology.

Contact occupancy is computed as frames with any atom-pair contact for a residue pair divided by all analyzed frames for the independent trajectory. Segments are pooled within a trajectory ID; they must not be counted as independent replicates. Frames with no contacts are included only if the complete frame manifest was written and supplied. The cross-trajectory summary assigns zero occupancy to a pair in trajectories where that pair was not observed, for all pairs in the union of observed pairs.

The current legacy-format fallback for contact tables without topology residue indices emits a warning because duplicate residue numbering can be conflated. Re-export through the updated contact tool for robust topology residue identity.

## Example

See `configs/salt-bridges.example.json` and the command examples in the main README. The example residue numbers and segment IDs are not a validated universal definition; replace them after inspecting your topology.
