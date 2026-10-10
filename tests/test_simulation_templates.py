from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
OPENMM_SCRIPT = ROOT / "simulation_templates/openmm/run_restraint_simulation.py"
OPENMM_CONFIG = ROOT / "simulation_templates/openmm/copper_restraints.example.json"

spec = importlib.util.spec_from_file_location("openmm_restraint_template", OPENMM_SCRIPT)
assert spec is not None and spec.loader is not None
openmm_template = importlib.util.module_from_spec(spec)
spec.loader.exec_module(openmm_template)


def example_config():
    return json.loads(OPENMM_CONFIG.read_text(encoding="utf-8"))


def test_openmm_template_config_schema_passes_without_openmm():
    config = openmm_template.validate_config(example_config())
    assert config["template_status"] == "template"
    assert len(config["restraints"]["bonds"]) == 2
    assert config["restraints"]["angles"][0]["theta0_degrees"] == 120.0


def test_openmm_rejects_negative_restraint_constants():
    config = example_config()
    config["restraints"]["bonds"][0]["k_kj_mol_nm2"] = -10.0
    with pytest.raises(openmm_template.ConfigError, match="must be >= 0"):
        openmm_template.validate_config(config)


def test_openmm_rejects_ignored_switch_distance_with_nocutoff():
    config = example_config()
    config["simulation"]["nonbonded_method"] = "NoCutoff"
    with pytest.raises(openmm_template.ConfigError, match="has no use with NoCutoff"):
        openmm_template.validate_config(config)


def test_openmm_rejects_unbounded_minimization_by_default():
    config = example_config()
    config["simulation"]["minimization_max_iterations"] = 0
    with pytest.raises(openmm_template.ConfigError, match="must be >= 1"):
        openmm_template.validate_config(config)


def test_openmm_check_config_does_not_require_openmm_installation():
    result = subprocess.run(
        [sys.executable, str(OPENMM_SCRIPT), "--config", str(OPENMM_CONFIG), "--check-config"],
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "CONFIG SCHEMA PASSED" in result.stdout
    assert "Template status: template" in result.stdout


def test_openmm_run_is_blocked_for_untouched_example_config():
    result = subprocess.run(
        [sys.executable, str(OPENMM_SCRIPT), "--config", str(OPENMM_CONFIG), "--run", "--confirm-model-reviewed"],
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == 2
    assert "Refusing to run" in result.stderr
    assert "template_status" in result.stderr


def test_namd_templates_are_structurally_complete_tcl():
    tclsh = shutil.which("tclsh")
    if not tclsh:
        pytest.skip("tclsh not installed; NAMD templates still require engine validation")
    files = [
        ROOT / "simulation_templates/namd/equilibration_nvt.conf.example",
        ROOT / "simulation_templates/namd/continue_npt.conf.example",
    ]
    for path in files:
        tcl_script = f"set f [open {{{path}}} r]\nset contents [read $f]\nclose $f\nputs [info complete $contents]\n"
        result = subprocess.run([tclsh], input=tcl_script, text=True, capture_output=True, check=False)
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "1"


def test_namd_template_paths_include_early_stop_and_model_acknowledgement():
    for relative in (
        "simulation_templates/namd/equilibration_nvt.conf.example",
        "simulation_templates/namd/continue_npt.conf.example",
    ):
        content = (ROOT / relative).read_text(encoding="utf-8")
        assert "template_acknowledgement" in content
        assert "MODEL_REVIEWED" in content
        assert "CHANGE_ME" in content
        assert "proc require_file" in content


def _run_namd_template_with_tcl_shim(source: Path, tmp_path: Path, substitutions: dict[str, str], files: list[str]):
    tclsh = shutil.which("tclsh")
    if not tclsh:
        pytest.skip("tclsh not installed; NAMD templates still require engine validation")
    work = tmp_path / source.stem
    work.mkdir()
    for name in files:
        path = work / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("placeholder\n", encoding="utf-8")
    content = source.read_text(encoding="utf-8")
    for old, new in substitutions.items():
        assert old in content, f"Expected placeholder {old!r} in {source.name}"
        content = content.replace(old, new)
    dest = work / source.name
    dest.write_text(content, encoding="utf-8")
    tcl_driver = f"proc unknown {{command args}} {{ return }}\nsource {{{dest}}}\nputs TEMPLATE_EXECUTED\n"
    result = subprocess.run([tclsh], input=tcl_driver, text=True, capture_output=True, cwd=work, check=False)
    assert result.returncode == 0, result.stderr
    assert not result.stderr.strip(), f"Unexpected Tcl error/warning output: {result.stderr}"
    assert "TEMPLATE_EXECUTED" in result.stdout


def test_namd_nvt_template_control_flow_with_mocked_namd_commands(tmp_path):
    source = ROOT / "simulation_templates/namd/equilibration_nvt.conf.example"
    _run_namd_template_with_tcl_shim(
        source,
        tmp_path,
        {
            'template_acknowledgement "EDIT_ME"': 'template_acknowledgement "MODEL_REVIEWED"',
            '"CHANGE_ME_system_nvt"': '"smoke_nvt"',
            '"CHANGE_ME/system.psf"': '"system.psf"',
            '"CHANGE_ME/system.pdb"': '"system.pdb"',
            '"CHANGE_ME/par_all36m_prot.prm"': '"par.prm"',
            '"CHANGE_ME/toppar_water_ions.str"': '"water.str"',
        },
        ["system.psf", "system.pdb", "par.prm", "water.str"],
    )


def test_namd_npt_template_control_flow_with_mocked_namd_commands(tmp_path):
    source = ROOT / "simulation_templates/namd/continue_npt.conf.example"
    _run_namd_template_with_tcl_shim(
        source,
        tmp_path,
        {
            'template_acknowledgement "EDIT_ME"': 'template_acknowledgement "MODEL_REVIEWED"',
            '"CHANGE_ME_system_npt"': '"smoke_npt"',
            '"CHANGE_ME/system.psf"': '"system.psf"',
            '"CHANGE_ME/par_all36m_prot.prm"': '"par.prm"',
            '"CHANGE_ME/toppar_water_ions.str"': '"water.str"',
            '"CHANGE_ME/output_nvt/CHANGE_ME_system_nvt"': '"prev/system_nvt"',
        },
        [
            "system.psf", "par.prm", "water.str",
            "prev/system_nvt.restart.coor", "prev/system_nvt.restart.vel", "prev/system_nvt.restart.xsc",
        ],
    )


@pytest.mark.parametrize("template_name", [
    "equilibration_nvt.conf.example",
    "continue_npt.conf.example",
])
def test_namd_position_restraints_can_be_enabled_with_matching_files(template_name, tmp_path):
    source = ROOT / "simulation_templates/namd" / template_name
    if template_name.startswith("equilibration"):
        substitutions = {
            'template_acknowledgement "EDIT_ME"': 'template_acknowledgement "MODEL_REVIEWED"',
            '"CHANGE_ME_system_nvt"': '"smoke_nvt"',
            '"CHANGE_ME/system.psf"': '"system.psf"',
            '"CHANGE_ME/system.pdb"': '"system.pdb"',
            '"CHANGE_ME/par_all36m_prot.prm"': '"par.prm"',
            '"CHANGE_ME/toppar_water_ions.str"': '"water.str"',
            'set use_position_restraints 0': 'set use_position_restraints 1',
            '"CHANGE_ME/restraint_reference.pdb"': '"reference.pdb"',
            '"CHANGE_ME/restraint_force_constants.pdb"': '"forces.pdb"',
        }
        files = ["system.psf", "system.pdb", "par.prm", "water.str", "reference.pdb", "forces.pdb"]
    else:
        substitutions = {
            'template_acknowledgement "EDIT_ME"': 'template_acknowledgement "MODEL_REVIEWED"',
            '"CHANGE_ME_system_npt"': '"smoke_npt"',
            '"CHANGE_ME/system.psf"': '"system.psf"',
            '"CHANGE_ME/par_all36m_prot.prm"': '"par.prm"',
            '"CHANGE_ME/toppar_water_ions.str"': '"water.str"',
            '"CHANGE_ME/output_nvt/CHANGE_ME_system_nvt"': '"prev/system_nvt"',
            'set use_position_restraints 0': 'set use_position_restraints 1',
            '"CHANGE_ME/restraint_reference.pdb"': '"reference.pdb"',
            '"CHANGE_ME/restraint_force_constants.pdb"': '"forces.pdb"',
        }
        files = ["system.psf", "par.prm", "water.str", "prev/system_nvt.restart.coor",
                 "prev/system_nvt.restart.vel", "prev/system_nvt.restart.xsc", "reference.pdb", "forces.pdb"]
    _run_namd_template_with_tcl_shim(source, tmp_path, substitutions, files)


@pytest.mark.parametrize("template_name", [
    "equilibration_nvt.conf.example",
    "continue_npt.conf.example",
])
def test_namd_position_restraint_flag_rejects_invalid_boolean(template_name, tmp_path):
    source = ROOT / "simulation_templates/namd" / template_name
    work = tmp_path / source.stem
    work.mkdir()
    content = source.read_text(encoding="utf-8")
    content = content.replace('template_acknowledgement "EDIT_ME"', 'template_acknowledgement "MODEL_REVIEWED"')
    content = content.replace('set use_position_restraints 0', 'set use_position_restraints maybe')
    if template_name.startswith("equilibration"):
        for old, new in {
            '"CHANGE_ME_system_nvt"': '"smoke_nvt"',
            '"CHANGE_ME/system.psf"': '"system.psf"',
            '"CHANGE_ME/system.pdb"': '"system.pdb"',
            '"CHANGE_ME/par_all36m_prot.prm"': '"par.prm"',
            '"CHANGE_ME/toppar_water_ions.str"': '"water.str"',
        }.items():
            content = content.replace(old, new)
        names = ["system.psf", "system.pdb", "par.prm", "water.str"]
    else:
        for old, new in {
            '"CHANGE_ME_system_npt"': '"smoke_npt"',
            '"CHANGE_ME/system.psf"': '"system.psf"',
            '"CHANGE_ME/par_all36m_prot.prm"': '"par.prm"',
            '"CHANGE_ME/toppar_water_ions.str"': '"water.str"',
            '"CHANGE_ME/output_nvt/CHANGE_ME_system_nvt"': '"prev/system_nvt"',
        }.items():
            content = content.replace(old, new)
        names = ["system.psf", "par.prm", "water.str", "prev/system_nvt.restart.coor",
                 "prev/system_nvt.restart.vel", "prev/system_nvt.restart.xsc"]
    for name in names:
        f = work / name
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text("placeholder\n", encoding="utf-8")
    dest = work / source.name
    dest.write_text(content, encoding="utf-8")
    driver = f"proc unknown {{command args}} {{ return }}\nsource {{{dest}}}\n"
    result = subprocess.run([shutil.which("tclsh") or "tclsh"], input=driver, text=True,
                            capture_output=True, cwd=work, check=False)
    assert result.returncode != 0
    assert "use_position_restraints must be 0/1" in result.stderr


@pytest.mark.parametrize("mutator, message", [
    (lambda c: c["simulation"].update(random_seed=0), "random_seed must be null or an integer >= 1"),
    (lambda c: (c["simulation"].update(pressure_bar=1.0), c["simulation"].update(barostat_frequency_steps=2.5)), "barostat_frequency_steps must be a positive integer"),
    (lambda c: c["outputs"].update(enforce_periodic_box=1), "enforce_periodic_box must be true, false, or null"),
])
def test_openmm_rejects_ambiguous_seed_frequency_or_boolean(mutator, message):
    config = example_config()
    mutator(config)
    with pytest.raises(openmm_template.ConfigError, match=message):
        openmm_template.validate_config(config)
