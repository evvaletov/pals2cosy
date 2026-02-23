"""Tests for pals2cosy lattice parser."""

import os
import pytest

from pals2cosy.parser import parse_lattice

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
UHFEL_YAML = os.path.join(EXAMPLES_DIR, "uhfel_beamline.yaml")


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
