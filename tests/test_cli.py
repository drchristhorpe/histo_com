from pathlib import Path

from click.testing import CliRunner

from histo_com.cli import main

FIXTURES = Path(__file__).parent / "fixtures"
ABD = str(FIXTURES / "1hhk_1_abd.cif")
PEPTIDE = str(FIXTURES / "1hhk_1_peptide.cif")
COMPLEX = str(FIXTURES / "8gvi_1_aligned.cif")


def run(*args):
    return CliRunner().invoke(main, list(args))


def test_cli_all_mode_default():
    result = run(ABD)
    assert result.exit_code == 0, result.output
    assert result.output.startswith("all:")


def test_cli_all_mode_explicit():
    result = run(PEPTIDE, "--mode", "all")
    assert result.exit_code == 0, result.output
    assert result.output.startswith("all:")


def test_cli_residues_mode():
    result = run(PEPTIDE, "--mode", "residues", "--residues", "1-9")
    assert result.exit_code == 0, result.output
    lines = result.output.strip().splitlines()
    assert len(lines) == 9
    assert lines[0].startswith("residue 1:")
    assert lines[-1].startswith("residue 9:")


def test_cli_domains_mode_single_chain():
    result = run(COMPLEX, "--mode", "domains", "--domains", "P")
    assert result.exit_code == 0, result.output
    lines = result.output.strip().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("P:")


def test_cli_domains_mode_chain_range():
    result = run(COMPLEX, "--mode", "domains", "--domains", "L:1-180")
    assert result.exit_code == 0, result.output
    lines = result.output.strip().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("L:1-180:")


def test_cli_domains_mode_chain_list():
    result = run(COMPLEX, "--mode", "domains", "--domains", "A,B")
    assert result.exit_code == 0, result.output
    lines = result.output.strip().splitlines()
    assert len(lines) == 2
    assert lines[0].startswith("A:")
    assert lines[1].startswith("B:")


def test_cli_domains_mode_requires_domains_option():
    result = run(COMPLEX, "--mode", "domains")
    assert result.exit_code != 0
    assert "--domains is required" in result.output


def test_cli_residues_mode_requires_residues_option():
    result = run(COMPLEX, "--mode", "residues")
    assert result.exit_code != 0
    assert "--residues is required" in result.output


def test_cli_unknown_chain_is_a_clean_error():
    result = run(COMPLEX, "--mode", "domains", "--domains", "Z")
    assert result.exit_code != 0
    assert "No such chain" in result.output


def test_cli_missing_file():
    result = run("does_not_exist.cif")
    assert result.exit_code != 0


def test_cli_output_writes_marker_pdb(tmp_path):
    out_path = tmp_path / "markers.pdb"
    result = run(
        COMPLEX, "--mode", "domains", "--domains", "P,L:1-180,A,B", "--output", str(out_path)
    )
    assert result.exit_code == 0, result.output
    assert f"Wrote centre-of-mass marker(s) to {out_path}" in result.output
    assert out_path.is_file()

    from histo_com.core import load_structure

    markers = load_structure(out_path)
    atoms = list(markers[0].get_atoms())
    assert len(atoms) == 4
    assert all(a.get_parent().resname == "COM" for a in atoms)


def test_cli_output_all_mode(tmp_path):
    out_path = tmp_path / "markers.pdb"
    result = run(ABD, "-o", str(out_path))
    assert result.exit_code == 0, result.output
    assert out_path.is_file()
