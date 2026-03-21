"""Tests for CLI auto-detection of input format."""

import os
import pytest

from pals2cosy.cli import _detect_format

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")


def test_detect_official_pals():
    """PALS: root key → official-pals."""
    fodo = os.path.join(EXAMPLES_DIR, "fodo.pals.yaml")
    assert _detect_format(fodo) == "official-pals"


def test_detect_official_pals_excerpt():
    """PALS: root key for UH FEL excerpt."""
    excerpt = os.path.join(EXAMPLES_DIR, "uhfel_excerpt.pals.yaml")
    assert _detect_format(excerpt) == "official-pals"


def test_detect_felsim_v2():
    """beamline: root key → felsim."""
    uhfel = os.path.join(os.path.dirname(__file__), "fixtures", "uhfel_beamline.yaml")
    assert _detect_format(uhfel) == "felsim"


def test_cli_auto_fodo(tmp_path):
    """End-to-end: auto-detect FODO → FOX output."""
    import subprocess
    fodo = os.path.join(EXAMPLES_DIR, "fodo.pals.yaml")
    out = tmp_path / "fodo.fox"
    result = subprocess.run(
        ["python", "-m", "pals2cosy", fodo, "-o", str(out),
         "--ke", "1000", "--particle", "proton"],
        capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    fox = out.read_text()
    assert "RPP 1000" in fox
    assert "MQ" in fox


def test_cli_auto_felsim(tmp_path):
    """End-to-end: auto-detect flat beamline → FOX output."""
    import subprocess
    uhfel = os.path.join(os.path.dirname(__file__), "fixtures", "uhfel_beamline.yaml")
    out = tmp_path / "uhfel.fox"
    result = subprocess.run(
        ["python", "-m", "pals2cosy", uhfel, "-o", str(out), "--ke", "40"],
        capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    fox = out.read_text()
    assert "RPE 40" in fox


def test_cli_explicit_mode_official_pals(tmp_path):
    """Explicit --mode official-pals works."""
    import subprocess
    fodo = os.path.join(EXAMPLES_DIR, "fodo.pals.yaml")
    out = tmp_path / "fodo.fox"
    result = subprocess.run(
        ["python", "-m", "pals2cosy", fodo, "-o", str(out),
         "--mode", "official-pals", "--ke", "1000", "--particle", "proton"],
        capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


def test_cli_beamline_flag(tmp_path):
    """--beamline selects a specific BeamLine from the facility."""
    import subprocess
    fodo = os.path.join(EXAMPLES_DIR, "fodo.pals.yaml")
    out = tmp_path / "cell.fox"
    result = subprocess.run(
        ["python", "-m", "pals2cosy", fodo, "-o", str(out),
         "--beamline", "fodo_cell", "--ke", "1000", "--particle", "proton"],
        capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    fox = out.read_text()
    mq_lines = [l for l in fox.splitlines() if "MQ" in l and "VARIABLE" not in l]
    # fodo_cell has 2 quads (not repeated), so 2 MQ commands
    assert len(mq_lines) == 2


def test_mode_pals_is_strict_felsim():
    """--mode pals routes to the FELsim parser in strict mode (not official PALS parser).

    This documents the current (confusing) behavior: --mode pals rejects
    FELsim-native types but uses the FELsim parser, not the official PALS parser.
    The --mode felsim-strict alias (P3-4) will make this clearer.
    """
    # Official PALS file should fail with --mode pals (FELsim parser, no PALS: root)
    from pals2cosy.parser import parse_lattice
    fodo = os.path.join(EXAMPLES_DIR, "fodo.pals.yaml")
    # fodo.pals.yaml has 'PALS' root, not 'beamline', so FELsim parser will fail
    with pytest.raises(KeyError):
        parse_lattice(fodo, mode="pals")


def test_beamline_flag_non_pals_warns(tmp_path):
    """--beamline with felsim mode emits a warning."""
    import subprocess
    uhfel = os.path.join(os.path.dirname(__file__), "fixtures", "uhfel_beamline.yaml")
    out = tmp_path / "out.fox"
    result = subprocess.run(
        ["python", "-m", "pals2cosy", uhfel, "-o", str(out),
         "--ke", "40", "--beamline", "ignored"],
        capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert "--beamline is only used with official-pals" in result.stderr


def test_detect_format_toml(tmp_path):
    """TOML file with PALS root key is detected as official-pals."""
    p = str(tmp_path / "test.toml")
    with open(p, "w") as f:
        f.write('[PALS]\nname = "test"\n')
    assert _detect_format(p) == "official-pals"


def test_detect_format_toml_felsim(tmp_path):
    """TOML file with beamline root key is detected as felsim."""
    p = str(tmp_path / "test.toml")
    with open(p, "w") as f:
        f.write('[beamline]\nname = "test"\n')
    assert _detect_format(p) == "felsim"


def test_detect_format_malformed_json(tmp_path):
    """Malformed JSON produces a clear error via _detect_format."""
    p = str(tmp_path / "bad.json")
    with open(p, "w") as f:
        f.write("{broken")
    with pytest.raises(ValueError, match="Cannot parse"):
        _detect_format(p)


def test_detect_format_malformed_yaml(tmp_path):
    """Malformed YAML produces a clear error via _detect_format."""
    p = str(tmp_path / "bad.yaml")
    with open(p, "w") as f:
        f.write("key: [unmatched\n  bad: indent")
    with pytest.raises(ValueError, match="Cannot parse"):
        _detect_format(p)


def test_detect_format_non_dict(tmp_path):
    """Non-dict YAML produces a clear error via _detect_format."""
    p = str(tmp_path / "list.yaml")
    with open(p, "w") as f:
        f.write("- item1\n- item2\n")
    with pytest.raises(ValueError, match="expected a mapping"):
        _detect_format(p)


def test_mode_felsim_strict_rejects_felsim_types(tmp_path):
    """--mode felsim-strict rejects FELsim-native type names."""
    import subprocess
    uhfel = os.path.join(os.path.dirname(__file__), "fixtures", "uhfel_beamline.yaml")
    out = tmp_path / "out.fox"
    result = subprocess.run(
        ["python", "-m", "pals2cosy", uhfel, "-o", str(out),
         "--ke", "40", "--mode", "felsim-strict"],
        capture_output=True, text=True
    )
    assert result.returncode != 0
    assert "not allowed in strict PALS mode" in result.stderr


def test_mode_pals_deprecation_warning(tmp_path):
    """--mode pals emits a DeprecationWarning."""
    import subprocess
    uhfel = os.path.join(os.path.dirname(__file__), "fixtures", "uhfel_beamline.yaml")
    out = tmp_path / "out.fox"
    result = subprocess.run(
        ["python", "-W", "all", "-m", "pals2cosy", uhfel, "-o", str(out),
         "--ke", "40", "--mode", "pals"],
        capture_output=True, text=True
    )
    assert "deprecated" in result.stderr.lower()
    assert "felsim-strict" in result.stderr
