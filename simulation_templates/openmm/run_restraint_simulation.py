#!/usr/bin/env python3
"""Config-driven OpenMM simulation template with optional metal-geometry restraints.

This script is intentionally independent of any one copper-trafficking system.
It uses atom identities (chain/residue/atom names), not copied particle indices.
Do not run a scientific production simulation until force-field parameters,
atom mappings, restraint targets/units and the protocol have been reviewed.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any


ALLOWED_NONBONDED_METHODS = {
    "NoCutoff",
    "CutoffNonPeriodic",
    "CutoffPeriodic",
    "Ewald",
    "PME",
    "LJPME",
}
ALLOWED_CONSTRAINTS = {None, "HBonds", "AllBonds", "HAngles"}
REVIEW_FLAGS = (
    "force_field_reviewed",
    "atom_mapping_reviewed",
    "restraint_geometry_and_units_reviewed",
    "simulation_protocol_reviewed",
)


class ConfigError(ValueError):
    """Raised when a configuration is incomplete or internally inconsistent."""


def _require_mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigError(f"{label} must be a JSON object")
    return value


def _require_positive_number(value: Any, label: str, *, allow_zero: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigError(f"{label} must be a finite number")
    if value < 0 or (value == 0 and not allow_zero):
        comparator = ">= 0" if allow_zero else "> 0"
        raise ConfigError(f"{label} must be {comparator}")
    return float(value)


def validate_config(config: Any) -> dict[str, Any]:
    """Validate the config schema without importing OpenMM or opening input files."""
    config = _require_mapping(config, "top-level config")
    if config.get("schema_version") != 1:
        raise ConfigError("schema_version must be 1")
    if config.get("template_status") not in ("template", "user_reviewed"):
        raise ConfigError("template_status must be 'template' or 'user_reviewed'")

    inputs = _require_mapping(config.get("inputs"), "inputs")
    if not isinstance(inputs.get("pdb"), str) or not inputs["pdb"].strip():
        raise ConfigError("inputs.pdb must be a non-empty path string")
    xmls = inputs.get("forcefield_xml")
    if not isinstance(xmls, list) or not xmls or any(not isinstance(x, str) or not x.strip() for x in xmls):
        raise ConfigError("inputs.forcefield_xml must be a non-empty list of path strings")

    simulation = _require_mapping(config.get("simulation"), "simulation")
    _require_positive_number(simulation.get("temperature_K"), "simulation.temperature_K")
    _require_positive_number(simulation.get("friction_per_ps"), "simulation.friction_per_ps")
    _require_positive_number(simulation.get("timestep_fs"), "simulation.timestep_fs")
    for key in ("steps", "minimization_max_iterations", "trajectory_interval_steps", "log_interval_steps", "checkpoint_interval_steps"):
        val = simulation.get(key)
        min_ok_zero = key == "minimization_max_iterations"
        _require_positive_number(val, f"simulation.{key}", allow_zero=min_ok_zero)
        if int(val) != val:
            raise ConfigError(f"simulation.{key} must be an integer")
    if int(simulation["steps"]) < 1:
        raise ConfigError("simulation.steps must be >= 1")
    if int(simulation["trajectory_interval_steps"]) < 1 or int(simulation["log_interval_steps"]) < 1 or int(simulation["checkpoint_interval_steps"]) < 1:
        raise ConfigError("output intervals must be >= 1")
    if not isinstance(simulation.get("minimize_before_run"), bool):
        raise ConfigError("simulation.minimize_before_run must be true or false")
    if simulation["minimize_before_run"] and int(simulation["minimization_max_iterations"]) < 1:
        raise ConfigError("simulation.minimization_max_iterations must be >= 1 when minimize_before_run is true")
    if simulation.get("constraints") not in (None, "HBonds", "AllBonds", "HAngles"):
        raise ConfigError("simulation.constraints must be null, 'HBonds', 'AllBonds', or 'HAngles'")
    if simulation.get("nonbonded_method") not in tuple(ALLOWED_NONBONDED_METHODS):
        raise ConfigError(f"simulation.nonbonded_method must be one of {sorted(ALLOWED_NONBONDED_METHODS)}")
    if simulation["nonbonded_method"] != "NoCutoff":
        _require_positive_number(simulation.get("nonbonded_cutoff_nm"), "simulation.nonbonded_cutoff_nm")
    if simulation.get("switch_distance_nm") is not None:
        _require_positive_number(simulation["switch_distance_nm"], "simulation.switch_distance_nm")
        if simulation["nonbonded_method"] == "NoCutoff":
            raise ConfigError("switch_distance_nm has no use with NoCutoff; set it to null or choose a cutoff method")
        if simulation.get("nonbonded_cutoff_nm") is not None and simulation["switch_distance_nm"] >= simulation["nonbonded_cutoff_nm"]:
            raise ConfigError("simulation.switch_distance_nm must be less than nonbonded_cutoff_nm")
    if not isinstance(simulation.get("rigid_water"), bool):
        raise ConfigError("simulation.rigid_water must be true or false")
    if not isinstance(simulation.get("remove_center_of_mass_motion"), bool):
        raise ConfigError("simulation.remove_center_of_mass_motion must be true or false")
    if not isinstance(simulation.get("platform"), str) or not simulation["platform"].strip():
        raise ConfigError("simulation.platform must be a platform name such as 'CPU' or 'CUDA'")
    props = simulation.get("platform_properties", {})
    if not isinstance(props, dict) or any(not isinstance(k, str) or not isinstance(v, (str, int, float)) for k, v in props.items()):
        raise ConfigError("simulation.platform_properties must be an object of scalar property values")
    pressure = simulation.get("pressure_bar")
    if pressure is not None:
        _require_positive_number(pressure, "simulation.pressure_bar", allow_zero=True)
        if simulation["nonbonded_method"] not in {"CutoffPeriodic", "Ewald", "PME", "LJPME"}:
            raise ConfigError("pressure_bar requires a periodic nonbonded method")
    if pressure is not None:
        frequency = simulation.get("barostat_frequency_steps")
        _require_positive_number(frequency, "simulation.barostat_frequency_steps")
        if isinstance(frequency, bool) or not isinstance(frequency, int):
            raise ConfigError("simulation.barostat_frequency_steps must be a positive integer")

    seed = simulation.get("random_seed")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int) or seed < 1):
        raise ConfigError("simulation.random_seed must be null or an integer >= 1 (0 is reserved by OpenMM for a non-fixed seed)")

    box = inputs.get("box_lengths_nm")
    if box is not None:
        if not isinstance(box, list) or len(box) != 3:
            raise ConfigError("inputs.box_lengths_nm must be null or a three-number orthorhombic box")
        for i, v in enumerate(box):
            _require_positive_number(v, f"inputs.box_lengths_nm[{i}]")

    outputs = _require_mapping(config.get("outputs"), "outputs")
    for key in ("directory", "trajectory", "log", "final_pdb", "checkpoint", "atom_mapping_csv", "run_summary_json"):
        if not isinstance(outputs.get(key), str) or not outputs[key].strip():
            raise ConfigError(f"outputs.{key} must be a non-empty path string")
    if not isinstance(outputs.get("write_pdb_trajectory"), bool):
        raise ConfigError("outputs.write_pdb_trajectory must be true or false")
    for key in ("pdb_trajectory", "final_checkpoint"):
        value = outputs.get(key)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ConfigError(f"outputs.{key} must be null or a non-empty path string")
    periodic_output = outputs.get("enforce_periodic_box")
    if periodic_output is not None and not isinstance(periodic_output, bool):
        raise ConfigError("outputs.enforce_periodic_box must be true, false, or null")

    restraints = _require_mapping(config.get("restraints"), "restraints")
    for category, n_atoms in (("bonds", 2), ("angles", 3), ("torsions", 4)):
        entries = restraints.get(category, [])
        if not isinstance(entries, list):
            raise ConfigError(f"restraints.{category} must be a list")
        for i, entry_any in enumerate(entries):
            entry = _require_mapping(entry_any, f"restraints.{category}[{i}]")
            if not isinstance(entry.get("name"), str) or not entry["name"].strip():
                raise ConfigError(f"restraints.{category}[{i}].name must be a non-empty string")
            atoms = entry.get("atoms")
            if not isinstance(atoms, list) or len(atoms) != n_atoms:
                raise ConfigError(f"restraints.{category}[{i}].atoms must contain exactly {n_atoms} atom selectors")
            for j, selector_any in enumerate(atoms):
                selector = _require_mapping(selector_any, f"restraints.{category}[{i}].atoms[{j}]")
                if "atom_name" not in selector or not isinstance(selector["atom_name"], str) or not selector["atom_name"].strip():
                    raise ConfigError("each atom selector must include a non-empty atom_name")
                supported = {"chain_id", "residue_id", "residue_name", "atom_name", "insertion_code", "element"}
                if set(selector) - supported:
                    raise ConfigError(f"unsupported atom selector keys: {sorted(set(selector) - supported)}")
                if "residue_id" not in selector and not ("chain_id" in selector and "residue_name" in selector):
                    raise ConfigError("each atom selector must include residue_id or both chain_id and residue_name")
            periodic = entry.get("periodic", True)
            if not isinstance(periodic, bool):
                raise ConfigError(f"restraints.{category}[{i}].periodic must be true or false")
            if category == "bonds":
                _require_positive_number(entry.get("r0_nm"), f"restraints.bonds[{i}].r0_nm")
                _require_positive_number(entry.get("k_kj_mol_nm2"), f"restraints.bonds[{i}].k_kj_mol_nm2", allow_zero=True)
            elif category == "angles":
                theta = _require_positive_number(entry.get("theta0_degrees"), f"restraints.angles[{i}].theta0_degrees")
                if theta >= 180:
                    raise ConfigError(f"restraints.angles[{i}].theta0_degrees must be < 180")
                _require_positive_number(entry.get("k_kj_mol_rad2"), f"restraints.angles[{i}].k_kj_mol_rad2", allow_zero=True)
            else:
                periodicity = entry.get("periodicity")
                if isinstance(periodicity, bool) or not isinstance(periodicity, int) or periodicity < 1:
                    raise ConfigError(f"restraints.torsions[{i}].periodicity must be an integer >= 1")
                phase = entry.get("phase_degrees")
                if isinstance(phase, bool) or not isinstance(phase, (int, float)) or not math.isfinite(phase):
                    raise ConfigError(f"restraints.torsions[{i}].phase_degrees must be finite")
                _require_positive_number(entry.get("k_kj_mol"), f"restraints.torsions[{i}].k_kj_mol", allow_zero=True)

    review = _require_mapping(config.get("review_checklist"), "review_checklist")
    if any(not isinstance(review.get(key), bool) for key in REVIEW_FLAGS):
        raise ConfigError(f"review_checklist must contain boolean fields: {', '.join(REVIEW_FLAGS)}")
    return config


def load_config(path: Path) -> dict[str, Any]:
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ConfigError(f"Cannot read config file {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise ConfigError(f"Invalid JSON at {path}:{exc.lineno}:{exc.colno}: {exc.msg}") from exc
    return validate_config(loaded)


def _resolve(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _resolve_output(base: Path, output_dir: Path, value: str) -> Path:
    """Resolve output names relative to outputs.directory unless absolute."""
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (output_dir / path).resolve()


def _selector_text(selector: dict[str, Any]) -> str:
    return ", ".join(f"{k}={v}" for k, v in selector.items())


def _resolve_atom(selector: dict[str, Any], topology: Any) -> Any:
    matches = []
    for atom in topology.atoms():
        residue = atom.residue
        chain_id = str(residue.chain.id)
        if "chain_id" in selector and chain_id != str(selector["chain_id"]):
            continue
        if "residue_id" in selector and str(residue.id) != str(selector["residue_id"]):
            continue
        if "residue_name" in selector and str(residue.name) != str(selector["residue_name"]):
            continue
        if "atom_name" in selector and str(atom.name) != str(selector["atom_name"]):
            continue
        if "insertion_code" in selector and str(getattr(residue, "insertionCode", "")) != str(selector["insertion_code"]):
            continue
        if "element" in selector:
            element_symbol = getattr(getattr(atom, "element", None), "symbol", None)
            if str(element_symbol) != str(selector["element"]):
                continue
        matches.append(atom)
    if len(matches) != 1:
        details = f"selector {{{_selector_text(selector)}}} matched {len(matches)} atoms; exactly one is required"
        if matches:
            details += "; matches: " + ", ".join(
                f"index={a.index} {a.residue.chain.id}/{a.residue.id}/{a.residue.name}/{a.name}" for a in matches[:8]
            )
        raise ConfigError(details)
    return matches[0]


def _openmm_imports() -> tuple[Any, Any, Any]:
    try:
        import openmm
        from openmm import app, unit
    except ImportError as exc:
        raise RuntimeError("OpenMM is required for --validate-model and --run. Install a compatible OpenMM build in the target environment.") from exc
    return openmm, app, unit


def _prepare(config: dict[str, Any], config_path: Path, platform_override: str | None = None) -> dict[str, Any]:
    openmm, app, unit = _openmm_imports()
    base = config_path.parent.resolve()
    pdb_path = _resolve(base, config["inputs"]["pdb"])
    xml_paths = [_resolve(base, x) for x in config["inputs"]["forcefield_xml"]]
    for path in [pdb_path, *xml_paths]:
        if not path.is_file():
            raise ConfigError(f"Required input file does not exist: {path}. Edit the template paths first.")

    pdb = app.PDBFile(str(pdb_path))
    box = config["inputs"].get("box_lengths_nm")
    if box is not None:
        pdb.topology.setUnitCellDimensions(openmm.Vec3(*[float(x) for x in box]) * unit.nanometer)

    sim_cfg = config["simulation"]
    periodic_methods = {"CutoffPeriodic", "Ewald", "PME", "LJPME"}
    has_periodic_restraints = any(
        item.get("periodic", True)
        for category in ("bonds", "angles", "torsions")
        for item in config["restraints"].get(category, [])
    )
    if (sim_cfg["nonbonded_method"] in periodic_methods or sim_cfg.get("pressure_bar") is not None or has_periodic_restraints):
        if pdb.topology.getPeriodicBoxVectors() is None:
            raise ConfigError(
                "This setup uses periodic nonbonded interactions, a barostat, or periodic restraints, "
                "but no periodic cell was found. Provide a valid CRYST1 cell in the PDB or set inputs.box_lengths_nm."
            )
    method = getattr(app, sim_cfg["nonbonded_method"])
    constraints = {
        None: None,
        "HBonds": app.HBonds,
        "AllBonds": app.AllBonds,
        "HAngles": app.HAngles,
    }[sim_cfg["constraints"]]
    create_kwargs: dict[str, Any] = {
        "nonbondedMethod": method,
        "constraints": constraints,
        "rigidWater": sim_cfg["rigid_water"],
        "removeCMMotion": sim_cfg["remove_center_of_mass_motion"],
        "ewaldErrorTolerance": float(sim_cfg.get("ewald_error_tolerance", 0.0005)),
    }
    if sim_cfg["nonbonded_method"] != "NoCutoff":
        create_kwargs["nonbondedCutoff"] = float(sim_cfg["nonbonded_cutoff_nm"]) * unit.nanometer
    if sim_cfg.get("switch_distance_nm") is not None:
        create_kwargs["switchDistance"] = float(sim_cfg["switch_distance_nm"]) * unit.nanometer
    forcefield = app.ForceField(*[str(p) for p in xml_paths])
    try:
        system = forcefield.createSystem(pdb.topology, **create_kwargs)
    except Exception as exc:
        raise RuntimeError(
            "OpenMM could not build the system from this topology and force-field XML set. "
            "Check residue templates, ion/copper atom types, charges, and input/force-field provenance. "
            f"Original error: {exc}"
        ) from exc

    topology_atoms = list(pdb.topology.atoms())
    if system.getNumParticles() < len(topology_atoms):
        raise RuntimeError(f"System has {system.getNumParticles()} particles but topology has {len(topology_atoms)} atoms")

    mapping_rows: list[dict[str, Any]] = []
    bond_groups: dict[bool, Any] = {}
    angle_groups: dict[bool, Any] = {}
    torsion_groups: dict[bool, Any] = {}

    def atom_index(selector: dict[str, Any], restraint_name: str, atom_slot: int) -> int:
        atom = _resolve_atom(selector, pdb.topology)
        if atom.index >= system.getNumParticles():
            raise RuntimeError(f"Topology atom index {atom.index} exceeds System particle count")
        mapping_rows.append({
            "restraint": restraint_name,
            "atom_slot": atom_slot,
            "particle_index_zero_based": atom.index,
            "chain_id": str(atom.residue.chain.id),
            "residue_id": str(atom.residue.id),
            "insertion_code": str(getattr(atom.residue, "insertionCode", "")),
            "residue_name": str(atom.residue.name),
            "atom_name": str(atom.name),
            "element": str(getattr(getattr(atom, "element", None), "symbol", "")),
            "selector": _selector_text(selector),
        })
        return atom.index

    for item in config["restraints"].get("bonds", []):
        indices = [atom_index(sel, item["name"], j + 1) for j, sel in enumerate(item["atoms"])]
        pbc = item.get("periodic", True)
        if pbc not in bond_groups:
            bond_groups[pbc] = openmm.HarmonicBondForce()
        force = bond_groups[pbc]
        force.addBond(
            indices[0], indices[1],
            float(item["r0_nm"]) * unit.nanometer,
            float(item["k_kj_mol_nm2"]) * unit.kilojoule_per_mole / unit.nanometer**2,
        )
    for pbc, force in bond_groups.items():
        force.setUsesPeriodicBoundaryConditions(pbc)
        force.setName(f"User harmonic bond restraints (PBC={pbc})")
        system.addForce(force)

    for item in config["restraints"].get("angles", []):
        indices = [atom_index(sel, item["name"], j + 1) for j, sel in enumerate(item["atoms"])]
        pbc = item.get("periodic", True)
        if pbc not in angle_groups:
            angle_groups[pbc] = openmm.HarmonicAngleForce()
        force = angle_groups[pbc]
        force.addAngle(
            indices[0], indices[1], indices[2],
            math.radians(float(item["theta0_degrees"])) * unit.radian,
            float(item["k_kj_mol_rad2"]) * unit.kilojoule_per_mole / unit.radian**2,
        )
    for pbc, force in angle_groups.items():
        force.setUsesPeriodicBoundaryConditions(pbc)
        force.setName(f"User harmonic angle restraints (PBC={pbc})")
        system.addForce(force)

    for item in config["restraints"].get("torsions", []):
        indices = [atom_index(sel, item["name"], j + 1) for j, sel in enumerate(item["atoms"])]
        pbc = item.get("periodic", True)
        if pbc not in torsion_groups:
            torsion_groups[pbc] = openmm.PeriodicTorsionForce()
        force = torsion_groups[pbc]
        force.addTorsion(
            indices[0], indices[1], indices[2], indices[3],
            int(item["periodicity"]), math.radians(float(item["phase_degrees"])) * unit.radian,
            float(item["k_kj_mol"]) * unit.kilojoule_per_mole,
        )
    for pbc, force in torsion_groups.items():
        force.setUsesPeriodicBoundaryConditions(pbc)
        force.setName(f"User periodic torsion restraints (PBC={pbc})")
        system.addForce(force)

    pressure_bar = sim_cfg.get("pressure_bar")
    if pressure_bar is not None:
        if pdb.topology.getPeriodicBoxVectors() is None:
            raise ConfigError("pressure_bar is set but no periodic box is available; provide a valid CRYST1 box or box_lengths_nm")
        barostat = openmm.MonteCarloBarostat(
            float(pressure_bar) * unit.bar,
            float(sim_cfg["temperature_K"]) * unit.kelvin,
            int(sim_cfg["barostat_frequency_steps"]),
        )
        system.addForce(barostat)

    output_base = _resolve(base, config["outputs"]["directory"])
    output_base.mkdir(parents=True, exist_ok=True)
    mapping_path = _resolve_output(base, output_base, config["outputs"]["atom_mapping_csv"])
    mapping_path.parent.mkdir(parents=True, exist_ok=True)
    with mapping_path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["restraint", "atom_slot", "particle_index_zero_based", "chain_id", "residue_id", "insertion_code", "residue_name", "atom_name", "element", "selector"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(mapping_rows)

    platform_name = platform_override or sim_cfg["platform"]
    try:
        platform = openmm.Platform.getPlatformByName(platform_name)
    except Exception as exc:
        raise RuntimeError(f"OpenMM platform '{platform_name}' is unavailable: {exc}") from exc
    properties = {} if platform_override else {str(k): str(v) for k, v in sim_cfg.get("platform_properties", {}).items()}
    if properties:
        supported = set(platform.getPropertyNames())
        unsupported = set(properties) - supported
        if unsupported:
            raise ConfigError(f"Platform {platform_name} does not support property/properties {sorted(unsupported)}; supported: {sorted(supported)}")
    integrator = openmm.LangevinMiddleIntegrator(
        float(sim_cfg["temperature_K"]) * unit.kelvin,
        float(sim_cfg["friction_per_ps"]) / unit.picosecond,
        float(sim_cfg["timestep_fs"]) * unit.femtosecond,
    )
    seed = sim_cfg.get("random_seed")
    if seed is not None:
        integrator.setRandomNumberSeed(int(seed))
    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)
    return {
        "openmm": openmm,
        "app": app,
        "unit": unit,
        "pdb": pdb,
        "system": system,
        "simulation": simulation,
        "platform_name": platform_name,
        "mapping_path": mapping_path,
        "output_base": output_base,
        "config_base": base,
    }


def _energy_kj_mol(ctx: Any, unit: Any) -> float:
    return float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))


def run_validate_model(config: dict[str, Any], config_path: Path) -> int:
    prepared = _prepare(config, config_path, platform_override="Reference")
    simulation = prepared["simulation"]
    initial_energy = _energy_kj_mol(simulation.context, prepared["unit"])
    if not math.isfinite(initial_energy):
        raise RuntimeError("Initial potential energy is not finite; do not run dynamics")
    if simulation.context.getState(getForces=True).getForces(asNumpy=True).shape[0] != prepared["system"].getNumParticles():
        raise RuntimeError("Force evaluation returned an unexpected particle count")
    print("MODEL PREFLIGHT PASSED (no dynamics were run)")
    print(f"Particles: {prepared['system'].getNumParticles()}")
    print(f"Forces: {prepared['system'].getNumForces()}")
    print(f"Platform: {prepared['platform_name']}")
    print(f"Initial potential energy: {initial_energy:.6f} kJ/mol")
    print(f"Atom mapping report: {prepared['mapping_path']}")
    print("Review the reported mapping and potential energy before considering a production run.")
    return 0


def run_simulation(config: dict[str, Any], config_path: Path) -> int:
    if config.get("template_status") != "user_reviewed":
        raise ConfigError("Refusing to run: set template_status to 'user_reviewed' only after editing the paths and model parameters")
    metadata = config.get("metadata", {})
    model_system = metadata.get("model_system") if isinstance(metadata, dict) else None
    if not isinstance(model_system, str) or not model_system.strip() or "CHANGE_ME" in model_system or model_system == "SET_THIS_TO_YOUR_SYSTEM":
        raise ConfigError("Refusing to run: set metadata.model_system to a meaningful system identifier")
    failed = [key for key in REVIEW_FLAGS if config["review_checklist"].get(key) is not True]
    if failed:
        raise ConfigError("Refusing to run: review checklist is incomplete: " + ", ".join(failed))

    prepared = _prepare(config, config_path)
    openmm, app, unit = prepared["openmm"], prepared["app"], prepared["unit"]
    simulation = prepared["simulation"]
    sim_cfg, out_cfg = config["simulation"], config["outputs"]
    base, output_base = prepared["config_base"], prepared["output_base"]

    def output_path(key: str) -> Path:
        value = out_cfg[key]
        return _resolve_output(base, output_base, value)

    log_path = output_path("log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    dcd_path = output_path("trajectory")
    dcd_path.parent.mkdir(parents=True, exist_ok=True)
    final_pdb_path = output_path("final_pdb")
    final_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output_path("checkpoint")
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

    initial_energy = _energy_kj_mol(simulation.context, unit)
    if not math.isfinite(initial_energy):
        raise RuntimeError("Initial potential energy is not finite; aborting before dynamics")
    if sim_cfg["minimize_before_run"]:
        simulation.minimizeEnergy(maxIterations=int(sim_cfg["minimization_max_iterations"]))
    # The input is a coordinate file, not a restart with velocities; initialize
    # Maxwell-Boltzmann velocities explicitly after minimization.
    velocity_seed = int(sim_cfg["random_seed"]) + 1 if sim_cfg.get("random_seed") is not None else 0
    simulation.context.setVelocitiesToTemperature(
        float(sim_cfg["temperature_K"]) * unit.kelvin,
        velocity_seed,
    )
    post_min_energy = _energy_kj_mol(simulation.context, unit)
    if not math.isfinite(post_min_energy):
        raise RuntimeError("Potential energy is not finite after minimization; aborting before dynamics")

    simulation.reporters.append(app.DCDReporter(
        str(dcd_path), int(sim_cfg["trajectory_interval_steps"]),
        enforcePeriodicBox=out_cfg.get("enforce_periodic_box"),
    ))
    if out_cfg["write_pdb_trajectory"]:
        pdb_traj = output_path("pdb_trajectory") if out_cfg.get("pdb_trajectory") else dcd_path.with_suffix(".pdb")
        pdb_traj.parent.mkdir(parents=True, exist_ok=True)
        simulation.reporters.append(app.PDBReporter(str(pdb_traj), int(sim_cfg["trajectory_interval_steps"])))
    simulation.reporters.append(app.StateDataReporter(
        str(log_path), int(sim_cfg["log_interval_steps"]),
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, speed=True, separator="\t",
    ))
    simulation.reporters.append(app.CheckpointReporter(str(checkpoint_path), int(sim_cfg["checkpoint_interval_steps"])))

    print(f"Running {int(sim_cfg['steps'])} steps on {prepared['platform_name']}")
    print(f"Initial potential energy: {initial_energy:.6f} kJ/mol")
    print(f"After minimization: {post_min_energy:.6f} kJ/mol")
    simulation.step(int(sim_cfg["steps"]))
    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    final_energy = float(final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    if not math.isfinite(final_energy):
        raise RuntimeError("Final potential energy is not finite")
    with final_pdb_path.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), handle, keepIds=True)
    final_checkpoint = output_path("final_checkpoint") if out_cfg.get("final_checkpoint") else checkpoint_path.with_name(checkpoint_path.stem + "_final.chk")
    final_checkpoint.parent.mkdir(parents=True, exist_ok=True)
    simulation.saveCheckpoint(str(final_checkpoint))

    summary = {
        "template_status": config["template_status"],
        "model_system": config.get("metadata", {}).get("model_system", "unspecified"),
        "openmm_version": getattr(openmm, "__version__", "unknown"),
        "platform": prepared["platform_name"],
        "particles": prepared["system"].getNumParticles(),
        "forces": prepared["system"].getNumForces(),
        "steps_requested": int(sim_cfg["steps"]),
        "temperature_K": float(sim_cfg["temperature_K"]),
        "timestep_fs": float(sim_cfg["timestep_fs"]),
        "initial_potential_energy_kj_mol": initial_energy,
        "post_minimization_potential_energy_kj_mol": post_min_energy,
        "final_potential_energy_kj_mol": final_energy,
        "input_pdb": str(_resolve(base, config["inputs"]["pdb"])),
        "forcefield_xml": [str(_resolve(base, p)) for p in config["inputs"]["forcefield_xml"]],
        "atom_mapping_csv": str(prepared["mapping_path"]),
        "trajectory": str(dcd_path),
        "final_pdb": str(final_pdb_path),
        "checkpoint": str(checkpoint_path),
        "final_checkpoint": str(final_checkpoint),
    }
    summary_path = output_path("run_summary_json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Finished. Summary: {summary_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path, help="JSON configuration; paths resolve relative to this file")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check-config", action="store_true", help="Check schema and numerical consistency without importing OpenMM")
    mode.add_argument("--validate-model", action="store_true", help="Build the OpenMM System and check atom mapping/initial energy; runs no dynamics")
    mode.add_argument("--run", action="store_true", help="Run the configured simulation (requires all review acknowledgements)")
    parser.add_argument("--confirm-model-reviewed", action="store_true", help="Explicit confirmation required for --run")
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config.resolve())
        if args.check_config:
            print("CONFIG SCHEMA PASSED")
            print(f"Template status: {config['template_status']}")
            if config["template_status"] == "template":
                print("NOTE: This is an example template; edit paths, selections, force constants and protocol before use.")
            return 0
        if args.validate_model:
            return run_validate_model(config, args.config.resolve())
        if args.run:
            if not args.confirm_model_reviewed:
                raise ConfigError("--run requires --confirm-model-reviewed as an explicit final confirmation")
            return run_simulation(config, args.config.resolve())
    except (ConfigError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
