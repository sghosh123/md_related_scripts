# Source archive provenance and migration status

The archive names below identify the historical workspaces reviewed. Archive data and trajectories are not bundled with the maintained code; only generalized implementations and useful documentation have been integrated.

| Source archive | Useful workflow | Repository disposition |
|---|---|---|
| `holo_sim_3_coord.zip` | Four-orientation yeast intermediate Cu–S and SG–NZ distance time series | Consolidated into `analyze_distance.py` and `analyze_multi_distance.py`; old time scale/offset left explicit because the reason is undocumented |
| `neighbor_test.zip` | Cysteine-centered water/protein/electronegative neighbor lists, contact counts and residue populations | Refactored into configurable PBC-aware neighbor-contact workflow; erroneous hardcoded selections and repeated wrong phase dictionaries not carried over |
| `colvar.zip` | Atom positions and group centroid/center-of-geometry workflows | Consolidated into `analyze_group_position.py`; requires an explicit “group made whole” acknowledgment |
| `holo_sim.zip` | OpenMM metal-coordination restraint scans, distance histograms, VMD output | Histogram analysis consolidated; independent config-driven OpenMM bond/angle/torsion-restraint template added; original model-specific constants/indices not copied |
| `restrain_test.zip` | Positional restraint examples, topology inspection, bond/angle restraint experiments, NAMD/OpenMM support | General OpenMM restraint runner and NAMD templates added; unrelated di-phenyl, methanol and third-party CHARMM-GUI files excluded |
| `Human_chaperone_apo.zip` | Apo human chaperone cysteine and binding-loop observables | Analysis concept merged into the configurable distance/contact/hydrogen-bond tools; generated `.dat` sets and duplicate archive omitted |
| `yeast_chaperone_apo.zip` | Apo chaperone cysteine proximity, histogram, VMD/movie helpers, NAMD setup | Geometry-screen functionality consolidated; VMD helpers generalized; generic NAMD templates added instead of copying system-specific configs |
| `yeast_MBD_apo.zip` | Apo MBD Cys13/Cys16 distances, cysteine neighborhood and visualization | Geometry/contact functionality consolidated; generic NAMD templates added instead of copying system-specific configs |
| `y-complex-apo-md.zip` | Atx1–Ccc2/MBD/chaperone acid-base distances, cysteine distances, frame extraction, visualization | Distance/export workflows generalized; generic NAMD templates added; mixed Bacillus/yeast configs not reused |
| `Wilson_MBD1_apo.zip`, `Wilson_MBD2_apo.zip`, `Wilson_MBD3_apo.zip` | MBD-specific NAMD/OpenMM setup and distance histogram patterns | Histogram pattern consolidated; independent NAMD templates added; original system-specific configs not copied pending provenance/input recovery |
| `Human_chaperone_apo(1).zip` | Duplicate of the human chaperone archive | Same size as the other archive and likely a duplicate, but not byte-hash-verified; no separate copy included pending checksum confirmation |
| `propKa_cys.zip` | Cysteine protonation/pKa calculation and analysis | Reserved for the separate `cysteine-protonation-analysis` project; not merged here |
| `qmmm_revision_calc.zip`, `paper_1.zip`, `paper_1a.zip`, `paper_1a(1).zip` | QM/QM-MM setup and model-system calculations | Reserved for existing solvent-shell/QM-MM and metal-ligand projects; not duplicated here |

## Data policy

Do not commit raw trajectories, restart/checkpoint files, or custom force-field files without a clear sharing plan and license/provenance review. The original neighbor ZIP contained large generated JSON/PNG collections; they have not been copied. New CSV output formats document atom identity, distances, frame/time, and PBC status directly.
