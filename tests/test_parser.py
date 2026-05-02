"""Tests for pals2cosy lattice parser."""

import math
import os
import tempfile
import warnings
import pytest

from pals2cosy.parser import parse_lattice

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
UHFEL_YAML = os.path.join(FIXTURES_DIR, "uhfel_beamline.yaml")


@pytest.fixture
def uhfel():
    return parse_lattice(UHFEL_YAML)


def test_beam_params(uhfel):
    bp, _ = uhfel
    assert bp["kinetic_energy_mev"] == 45.0
    assert bp["particle_type"] == "electron"
    assert bp["mass_mev"] == pytest.approx(0.51099895)
    assert bp["rf_frequency_hz"] == pytest.approx(2856e6)
    assert bp["quad_gradient"] == pytest.approx(2.694)


def test_element_count(uhfel):
    _, elements = uhfel
    # Count physical element types (excluding auto-inserted drifts)
    types = [e["type"] for e in elements if e["type"] != "DRIFT"]
    assert len(types) > 80  # UH FEL has ~90 non-drift elements


def test_type_resolution(uhfel):
    _, elements = uhfel
    non_drift = [e for e in elements if e["type"] != "DRIFT"]
    types = {e["type"] for e in non_drift}

    assert "QPF" in types
    assert "QPD" in types
    assert "DPH" in types
    assert "DPW" in types
    assert "BPM" in types
    assert "OTR" in types
    assert "UND" in types


def test_first_quad(uhfel):
    _, elements = uhfel
    quads = [e for e in elements if e["type"] in ("QPF", "QPD")]
    q1 = quads[0]
    assert q1["name"] == "LIN_QPF_004"
    assert q1["type"] == "QPF"
    assert q1["current"] == pytest.approx(0.885719309299156)
    assert q1["length"] == pytest.approx(0.0889, rel=1e-3)


def test_first_dipole_triplet(uhfel):
    _, elements = uhfel
    dpw_elements = [e for e in elements if e["type"] == "DPW"]
    dph_elements = [e for e in elements if e["type"] == "DPH"]

    # First DPW is DC1_DPW_017
    assert dpw_elements[0]["name"] == "DC1_DPW_017"
    assert dpw_elements[0]["wedge_angle"] == pytest.approx(0.0)

    # First DPH is DPHa
    assert dph_elements[0]["name"] == "DPHa"
    assert dph_elements[0]["angle"] == pytest.approx(1.5)


def test_chicane_enge_coefficients(uhfel):
    _, elements = uhfel
    # FC1 chicane dipoles have Enge coefficients on their entrance DPW
    fc1_dpws_with_enge = [
        e for e in elements
        if e["type"] == "DPW" and e.get("enge_coeffs") is not None
    ]
    assert len(fc1_dpws_with_enge) >= 8  # 4 FC1 + 4 FC2 entrance DPWs
    assert fc1_dpws_with_enge[0]["enge_coeffs"][0] == pytest.approx(56.49)


def test_drift_insertion(uhfel):
    _, elements = uhfel
    drifts = [e for e in elements if e["type"] == "DRIFT"]
    assert len(drifts) > 0

    # First drift should cover the gap before the first quad
    assert drifts[0]["length"] == pytest.approx(0.358775, rel=1e-6)


def test_undulator(uhfel):
    _, elements = uhfel
    undulators = [e for e in elements if e["type"] == "UND"]
    assert len(undulators) == 2
    assert undulators[0]["length"] == pytest.approx(0.5405, rel=1e-3)


def test_negative_angle_dipoles(uhfel):
    _, elements = uhfel
    neg_dipoles = [e for e in elements if e["type"] == "DPH" and e.get("angle", 0) < 0]
    assert len(neg_dipoles) > 0
    # DC2 and DC4 dipoles bend negative
    assert neg_dipoles[0]["angle"] == pytest.approx(-4.0)


def test_pals_mode_rejects_felsim_types():
    """Strict PALS mode should reject FELsim-native types like DIPOLE_WEDGE."""
    # The UH FEL YAML uses FELsim-native types (DIPOLE_WEDGE), so pals mode
    # should raise ValueError
    with pytest.raises(ValueError, match="not allowed in strict PALS mode"):
        parse_lattice(UHFEL_YAML, mode="pals")


def test_felsim_mode_accepts_all(uhfel):
    """Default felsim mode should accept all type names."""
    # uhfel fixture uses default mode="felsim", should parse without error
    _, elements = uhfel
    assert len(elements) > 0


# ---------------------------------------------------------------------------
# v3 format support (P2-2)
# ---------------------------------------------------------------------------

def _write_v3_file(data, suffix=".yaml"):
    """Write data to a temp file."""
    import yaml
    f = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    yaml.dump(data, f)
    f.close()
    return f.name


def _v3_lattice(elements, format_version=3, global_settings=None):
    """Build a minimal FELsim v3 lattice dict."""
    d = {
        "beamline": {
            "metadata": {
                "name": "test", "version": "1.0",
                "format_version": format_version,
                "reference_energy_mev": 45.0,
                "particle_type": "electron",
            },
            "beam_parameters": {
                "particle": {
                    "type": "electron",
                    "kinetic_energy_mev": 45.0,
                    "mass_mev": 0.51099895,
                    "charge_e": -1,
                },
                "rf_frequency_hz": 2856e6,
            },
            "elements": elements,
        }
    }
    if global_settings:
        d["beamline"]["global_settings"] = global_settings
    return d


def test_format_version_3_accepted():
    """FELsim v3 files parse without error."""
    data = _v3_lattice([{
        "name": "D1", "type": "DRIFT",
        "s_start_m": 0.0, "s_end_m": 0.5, "length_m": 0.5,
        "parameters": {},
    }])
    path = _write_v3_file(data)
    try:
        _, elements = parse_lattice(path)
        assert len(elements) == 1
    finally:
        os.unlink(path)


def test_v3_bn1_quad_extracted():
    """Bn1 from MagneticMultipoleP is extracted for quads."""
    data = _v3_lattice([{
        "name": "Q1", "type": "Quadrupole", "polarity": "focusing",
        "s_start_m": 0.0, "s_end_m": 0.1, "length_m": 0.1,
        "parameters": {"current_a": 1.0},
        "MagneticMultipoleP": {"Bn1": -0.05},
    }])
    path = _write_v3_file(data)
    try:
        _, elements = parse_lattice(path)
        q = [e for e in elements if e["type"] == "QPF"][0]
        assert q["bn1"] == pytest.approx(-0.05)
    finally:
        os.unlink(path)


def test_v3_bendp_angle():
    """BendP.g_ref → angle computed from g_ref × L."""
    g_ref = 0.3
    L = 0.1
    data = _v3_lattice([{
        "name": "B1", "type": "SBend",
        "s_start_m": 0.0, "s_end_m": L, "length_m": L,
        "parameters": {"bending_angle_deg": 0, "dipole_length_m": L},
        "BendP": {"g_ref": g_ref, "e1": 0.1, "e2": 0.2},
    }])
    path = _write_v3_file(data)
    try:
        _, elements = parse_lattice(path)
        d = [e for e in elements if e["type"] == "DPH"][0]
        assert d["angle"] == pytest.approx(math.degrees(g_ref * L), rel=1e-10)
        assert d["entrance_edge_angle"] == pytest.approx(math.degrees(0.1), rel=1e-10)
        assert d["exit_edge_angle"] == pytest.approx(math.degrees(0.2), rel=1e-10)
    finally:
        os.unlink(path)


def test_format_version_99_rejected():
    """Unsupported format_version raises ValueError."""
    data = _v3_lattice([], format_version=99)
    path = _write_v3_file(data)
    try:
        with pytest.raises(ValueError, match="Unsupported format_version"):
            parse_lattice(path)
    finally:
        os.unlink(path)


# ---------------------------------------------------------------------------
# Error path tests (P2-3)
# ---------------------------------------------------------------------------

def test_empty_file_raises(tmp_path):
    """Empty YAML file raises an error."""
    p = str(tmp_path / "empty.yaml")
    with open(p, "w") as f:
        f.write("")
    with pytest.raises((ValueError, TypeError, KeyError)):
        parse_lattice(p)


def test_non_dict_yaml_root_raises(tmp_path):
    """YAML that parses to a list raises a clear ValueError."""
    import yaml
    p = str(tmp_path / "list.yaml")
    with open(p, "w") as f:
        yaml.dump([1, 2, 3], f)
    with pytest.raises(ValueError, match="mapping"):
        parse_lattice(p)


def test_missing_beamline_key_raises(tmp_path):
    """YAML without 'beamline' root key raises a clear ValueError."""
    import yaml
    p = str(tmp_path / "nobeamline.yaml")
    with open(p, "w") as f:
        yaml.dump({"other": {}}, f)
    with pytest.raises(ValueError, match="beamline"):
        parse_lattice(p)


def test_inverted_positions_rejected():
    """Elements with s_end < s_start are skipped with a warning."""
    data = _v3_lattice([{
        "name": "BAD", "type": "DRIFT",
        "s_start_m": 1.0, "s_end_m": 0.5, "length_m": 0.5,
        "parameters": {},
    }])
    path = _write_v3_file(data)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, elements = parse_lattice(path)
        # Element should be skipped
        non_drift = [e for e in elements if e["type"] != "DRIFT"]
        assert len(non_drift) == 0
        assert any("s_end" in str(warning.message) for warning in w)
    finally:
        os.unlink(path)


# --- Phase 2C: Additional edge case tests ---

def test_zero_current_quad():
    """Quadrupole with current_a=0 produces MQ with b_pole=0."""
    from pals2cosy.converter import convert
    data = _v3_lattice([{
        "name": "Q0", "type": "Quadrupole", "polarity": "focusing",
        "s_start_m": 0.0, "s_end_m": 0.1, "length_m": 0.1,
        "parameters": {"current_a": 0},
    }])
    path = _write_v3_file(data)
    try:
        bp, elements = parse_lattice(path)
        fox = convert(bp, elements)
        mq_lines = [l.strip() for l in fox.splitlines() if l.strip().startswith("MQ")]
        assert len(mq_lines) == 1
        parts = mq_lines[0].split()
        assert float(parts[2]) == pytest.approx(0.0)
    finally:
        os.unlink(path)


def test_octupole_resolves_to_oct():
    """Octupole CamelCase type resolves to OCT and produces a DL drift."""
    from pals2cosy.converter import convert
    data = _v3_lattice([{
        "name": "O1", "type": "Octupole",
        "s_start_m": 0.0, "s_end_m": 0.2, "length_m": 0.2,
        "parameters": {},
    }])
    path = _write_v3_file(data)
    try:
        _, elements = parse_lattice(path)
        assert elements[0]["type"] == "OCT"
        bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
              "quad_gradient": 2.694, "source_file": "test"}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fox = convert(bp, elements)
        assert "DL 0.2" in fox
        assert not any("Unknown element type" in str(x.message) for x in w)
    finally:
        os.unlink(path)


def test_multipole_resolves_to_mult():
    """Multipole CamelCase type resolves to MULT and produces a DL drift."""
    from pals2cosy.converter import convert
    data = _v3_lattice([{
        "name": "M1", "type": "Multipole",
        "s_start_m": 0.0, "s_end_m": 0.15, "length_m": 0.15,
        "parameters": {},
    }])
    path = _write_v3_file(data)
    try:
        _, elements = parse_lattice(path)
        assert elements[0]["type"] == "MULT"
        bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
              "quad_gradient": 2.694, "source_file": "test"}
        fox = convert(bp, elements)
        assert "DL 0.15" in fox
    finally:
        os.unlink(path)


def test_overlapping_elements_raises():
    """Overlapping elements raise ValueError (fail-fast on physically invalid input)."""
    data = _v3_lattice([
        {"name": "D1", "type": "DRIFT",
         "s_start_m": 0.0, "s_end_m": 0.5, "length_m": 0.5,
         "parameters": {}},
        {"name": "D2", "type": "DRIFT",
         "s_start_m": 0.3, "s_end_m": 0.8, "length_m": 0.5,
         "parameters": {}},
    ])
    path = _write_v3_file(data)
    try:
        with pytest.raises(ValueError, match="overlaps"):
            parse_lattice(path)
    finally:
        os.unlink(path)


def test_missing_type_warns():
    """Element without type or kind emits a warning and is skipped."""
    data = _v3_lattice([{
        "name": "MYSTERY",
        "s_start_m": 0.0, "s_end_m": 0.5, "length_m": 0.5,
        "parameters": {},
    }])
    path = _write_v3_file(data)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, elements = parse_lattice(path)
        non_drift = [e for e in elements if e["type"] != "DRIFT"]
        assert len(non_drift) == 0
        assert any("no 'type' or 'kind'" in str(x.message) for x in w)
    finally:
        os.unlink(path)


def test_missing_positions_warns():
    """Element without s_start_m/s_end_m emits a warning and is skipped."""
    data = _v3_lattice([{
        "name": "NOLOC", "type": "DRIFT",
        "length_m": 0.5,
        "parameters": {},
    }])
    path = _write_v3_file(data)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, elements = parse_lattice(path)
        assert len(elements) == 0
        assert any("missing s_start_m" in str(x.message) for x in w)
    finally:
        os.unlink(path)


def test_format_version_1_accepted():
    """format_version 1 is accepted (elements still need s_start_m/s_end_m)."""
    data = _v3_lattice([{
        "name": "D1", "type": "DRIFT",
        "s_start_m": 0.0, "s_end_m": 0.5, "length_m": 0.5,
        "parameters": {},
    }], format_version=1)
    path = _write_v3_file(data)
    try:
        _, elements = parse_lattice(path)
        assert len(elements) == 1
    finally:
        os.unlink(path)
