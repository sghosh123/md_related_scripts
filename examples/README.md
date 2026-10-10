# Examples

All configurations use illustrative selections. Check the topology interactively before production analysis, especially where two protein segments reuse residue numbers.

## 1. One atom pair

```bash
python scripts/analyze_distance.py \
  --topology system.psf --trajectory trajectory.dcd \
  --selection-a 'segid PROA and resid 16 and name SG' \
  --selection-b 'resname CU1' \
  --time-scale 1 --time-unit ps \
  --output-csv results/cu_cys16.csv \
  --output-plot results/cu_cys16.png
```

## 2. Several named distances in one pass

Edit `configs/multi-distance.example.json` so every selection identifies exactly one atom, then run:

```bash
python scripts/analyze_multi_distance.py \
  --topology system.psf --trajectory trajectory.dcd \
  --config configs/multi-distance.example.json \
  --output-csv results/distances.csv \
  --output-plot results/distances.png
```

## 3. Water or protein neighbors

```bash
python scripts/analyze_neighbors.py \
  --topology system.psf --trajectory trajectory.dcd \
  --center-selection 'segid PROT and resid 13 and name SG' \
  --neighbor-selection 'resname HOH TIP3 WAT' \
  --cutoff 3.5 \
  --output-contacts results/water_contacts.csv \
  --output-summary results/water_neighbor_summary.csv \
  --output-residue-frequency results/water_occupancy.csv \
  --output-pair-frequency results/residue_pair_occupancy.csv \
  --output-plot results/water_neighbor_count.png
```

The long-form contact file contains one row per center-atom/neighbor-atom pair within the cutoff. Summary counts are unique neighbor atoms/residues per frame. The default cutoff is not inferred from legacy filenames; set it deliberately.

## Inter-protein contact network / acidic–basic contact screening

The neighbor command can also describe a residue-pair contact network. This example mirrors the *idea* of the historical PROA/PROB analysis but intentionally broadens acidic oxygen atom names rather than silently checking OD1 alone:

```bash
python scripts/analyze_neighbors.py \
  --topology complex.psf --trajectory complex.dcd \
  --center-selection 'segid PROA and (resid 32 37 46 54) and (name OD1 OD2 OE1 OE2)' \
  --neighbor-selection 'segid PROB and (resid 28 32 56) and (name NZ NH1 NH2 NE)' \
  --cutoff 4.0 \
  --output-contacts results/acid_base_atom_contacts.csv \
  --output-summary results/acid_base_frame_summary.csv \
  --output-pair-frequency results/acid_base_residue_pair_frequency.csv
```

These residue numbers/segment IDs come from one historical workspace and are **illustrative only**. Check residue names and protonation states carefully. The output is an atom-distance contact network; do not call every edge a salt bridge solely because of the atom names.

## 4. Geometry-based hydrogen-bond screen

This requires topology bonds, not only coordinates:

```bash
python scripts/analyze_hbond_geometry.py \
  --topology system.psf --trajectory trajectory.dcd \
  --donor-selection 'protein and name N NZ NE NH1 NH2 OG OG1 OH' \
  --hydrogen-selection 'protein and name H HN HZ1 HZ2 HZ3 HG HG1 HE' \
  --acceptor-selection 'segid PROT and resid 13 and name SG' \
  --distance-cutoff 3.5 --angle-cutoff 120 \
  --output-contacts results/hbond_geometry.csv \
  --output-summary results/hbond_geometry_summary.csv
```

This detects the operational criterion `D–A ≤ 3.5 Å` and `D–H–A ≥ 120°`; scientifically justify and report cutoffs. The selection examples are not a universal chemical typing system.

## 5. Group position / colvar support

```bash
python scripts/analyze_group_position.py \
  --topology system.psf --trajectory whole_trajectory.dcd \
  --selection 'segid PROA and resid 16 and name CA CB SG' \
  --method cog --assume-group-whole \
  --output-csv results/group_position.csv
```

The group must first be made whole/consistently unwrapped; an arithmetic center of atoms split across periodic boundaries is misleading. `--assume-group-whole` is an explicit acknowledgment of this requirement.

## 6. Extract frames and plot a distribution

```bash
python scripts/extract_frames.py \
  --topology system.psf --trajectory trajectory.dcd \
  --frames '0,10,100-105' --output-dir results/frames

python scripts/plot_distance_distribution.py \
  --input results/distances.csv --column Cu_to_Cys16_SG_distance_A \
  --bins 60 --output results/distance_histogram.png
```

For whitespace-delimited historical `.dat` files, use `--column-index` (zero-based, default 1); the script reads all available rows and does not assume a fixed sample count.


## 7. Dihedrals and binding-loop conformations

Edit `configs/dihedral-analysis.example.json` so the four atom selections for each torsion are unambiguous, then run:

```bash
python scripts/analyze_dihedrals.py \
  --topology system.psf --trajectory trajectory.dcd \
  --config configs/dihedral-analysis.example.json \
  --output-csv results/phi_psi.csv

python scripts/plot_ramachandran.py \
  --input results/phi_psi.csv --phi-column resid13_phi_angle_deg \
  --psi-column resid13_psi_angle_deg --output results/resid13_ramachandran.png
```

For amino-acid backbone angles, verify neighboring-residue identities and peptide continuity; these are explicit four-atom dihedrals, not automatic residue-chemistry recognition.

## 8. Fitted RMSD

```bash
python scripts/analyze_rmsd.py \
  --topology system.psf --trajectory trajectory.dcd \
  --selection 'protein and name CA' --reference-frame 0 \
  --assume-group-whole --output-csv results/rmsd.csv \
  --output-plot results/rmsd.png
```

The tool performs a least-squares fit on the selected atoms before reporting RMSD in Å. Handle imaging/whole-molecule reconstruction before the analysis and state whether RMSD is calculated for a whole complex, a domain, or a loop.
