"""Tests for pals2cosy FOX converter.

Validates generated FOX code against the reference input.fox from FELsim's
COSY cross-validation study. The reference was generated from the Excel
beamline at KE=40 MeV, FR=0, order=3, dim=3.
"""

import os
import re
import pytest

from pals2cosy.parser import parse_lattice
from pals2cosy.converter import convert, _consolidate_dipoles

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
UHFEL_YAML = os.path.join(EXAMPLES_DIR, "uhfel_beamline.yaml")
REFERENCE_FOX = os.path.join(
    os.path.dirname(__file__), "..", "..", "..",
    "UH", "fel-merge-workspace", "FELsim", "backend", "test", "results", "input.fox"
)


@pytest.fixture
def uhfel():
    return parse_lattice(UHFEL_YAML)


@pytest.fixture
def fox_output(uhfel):
    bp, elems = uhfel
    return convert(bp, elems, ke_override=40, fringe_field_order=0,
                   computation_order=3, dimensions=3)


def _extract_commands(fox_text, cmd_pattern):
    """Extract FOX commands matching a regex pattern."""
    return [line.strip() for line in fox_text.splitlines()
            if re.match(r'\s+' + cmd_pattern, line)]


# --- Structure tests ---

def test_fox_header(fox_output):
    assert "INCLUDE 'COSY' ;" in fox_output
    assert "PROCEDURE RUN ;" in fox_output
    assert "PROCEDURE LATTICE ;" in fox_output


def test_fox_variables(fox_output):
    assert "VARIABLE A0 100 3" in fox_output
    assert "VARIABLE F0 100 6" in fox_output


def test_fox_settings(fox_output):
    assert "OV 3 3 0 ;" in fox_output
    assert "RPE 40 ;" in fox_output
    assert "FR 0 ;" in fox_output


def test_fox_twiss_output(fox_output):
    assert "GT MAP F0 MU0 A0 B0 G0 R0 ;" in fox_output
    assert '"beta_x"' in fox_output
    assert '"alpha_y"' in fox_output


def test_fox_footer(fox_output):
    assert "RUN ;" in fox_output
    assert "END ;" in fox_output


# --- Consolidation tests ---

def test_dpw_dph_dpw_consolidation(uhfel):
    _, elements = uhfel
    consolidated = _consolidate_dipoles(elements)

    # All DPW-DPH-DPW triplets should be merged
    remaining_dpw = [e for e in consolidated if e["type"] == "DPW"]
    assert len(remaining_dpw) == 0, f"Unconsolidated DPW elements: {[e['name'] for e in remaining_dpw]}"

    dipoles = [e for e in consolidated if e["type"] == "DIPOLE_CONSOLIDATED"]
    assert len(dipoles) > 0


def test_consolidation_preserves_angles(uhfel):
    _, elements = uhfel
    consolidated = _consolidate_dipoles(elements)
    dipoles = [e for e in consolidated if e["type"] == "DIPOLE_CONSOLIDATED"]

    # First dipole: DC1 DPHa, angle=1.5, entrance=0.0, exit=1.5
    d1 = dipoles[0]
    assert d1["angle"] == pytest.approx(1.5)
    assert d1["entrance_angle"] == pytest.approx(0.0)
    assert d1["exit_angle"] == pytest.approx(1.5)


def test_consolidation_preserves_enge(uhfel):
    _, elements = uhfel
    consolidated = _consolidate_dipoles(elements)
    dipoles = [e for e in consolidated if e["type"] == "DIPOLE_CONSOLIDATED"]

    # FC1 chicane dipoles have Enge on entrance
    fc1_with_enge = [d for d in dipoles if d.get("entrance_enge")]
    assert len(fc1_with_enge) == 8  # 4 FC1 + 4 FC2


def test_negative_angle_consolidation(uhfel):
    _, elements = uhfel
    consolidated = _consolidate_dipoles(elements)
    dipoles = [e for e in consolidated if e["type"] == "DIPOLE_CONSOLIDATED"]

    neg_dipoles = [d for d in dipoles if d["angle"] < 0]
    assert len(neg_dipoles) > 0
    assert neg_dipoles[0]["angle"] == pytest.approx(-4.0)


# --- FOX command generation tests ---

def test_quad_formula(fox_output):
    """Verify MQ b_pole = sign * G * I * r."""
    mq_lines = _extract_commands(fox_output, r'MQ')
    assert len(mq_lines) > 0

    # First quad: QPF, I=0.885719309299156, G=2.694, r=0.0135
    # MQ L b_pole r ; → parts = ["MQ", L, b_pole, r, ";"]
    # b_pole = -2.694 * 0.885719309299156 * 0.0135
    expected_bpole = -2.694 * 0.885719309299156 * 0.0135
    parts = mq_lines[0].split()
    actual_bpole = float(parts[2])
    assert actual_bpole == pytest.approx(expected_bpole, rel=1e-12)


def test_dil_commands_match_reference(fox_output):
    """DIL commands should match the reference exactly."""
    dil_lines = _extract_commands(fox_output, r'DIL')

    # Reference DIL commands (from input.fox)
    ref_dils = [
        "DIL 0.08889999999999998 1.5 0.007239 0.0 0 1.5 0 ;",
        "DIL 0.0889000000000002 1.5 0.007239 0.75 0 0.75 0 ;",
        "DIL 0.04063999999999979 4.0 0.007239 2.018 0 2.018 0 ;",
        "DIL 0.04063999999999979 4.0 0.007239 2.018 0 2.018 0 ;",
        "DIL 0.04063999999999979 5.0 0.007239 2.536 0 2.536 0 ;",
        "DIL 0.04063999999999979 5.0 0.007239 2.536 0 2.536 0 ;",
        "DIL 0.04063999999999979 4.0 0.007239 2.018 0 2.018 0 ;",
        "DIL 0.04063999999999979 4.0 0.007239 2.018 0 2.018 0 ;",
        "DIL 0.08890000000000065 1.5 0.007239 0.75 0 0.75 0 ;",
        "DIL 0.08890000000000065 1.5 0.007239 1.5 0 0.0 0 ;",
        "DIL 0.03738899999999923 11.25 0.00635 0.0 0 11.25 0 ;",
        "DIL 0.03738899999999923 11.25 0.00635 11.25 0 0.0 0 ;",
        "DIL 0.03738899999999923 11.25 0.00635 0.0 0 11.25 0 ;",
        "DIL 0.03738899999999923 11.25 0.00635 11.25 0 0.0 0 ;",
    ]

    for i, ref in enumerate(ref_dils):
        assert dil_lines[i] == ref, f"DIL mismatch at index {i}: {dil_lines[i]} != {ref}"


def test_fc_commands(fox_output):
    """FC commands should match reference Enge coefficients."""
    fc_lines = _extract_commands(fox_output, r'FC')
    # All FC commands should have the same chicane Enge coefficients
    expected = "FC 1 1 1 56.49 -50.79 19.32 -3.621 0.3315 -0.01193 ;"
    for fc in fc_lines:
        assert fc == expected


def test_cb_wrapping(fox_output):
    """Negative-angle dipoles should be wrapped with CB."""
    lines = [l.strip() for l in fox_output.splitlines()]
    for i, line in enumerate(lines):
        if "DIL" in line:
            # Extract angle from DIL command
            parts = line.split()
            dil_angle = float(parts[2])
            # Check preceding lines for CB
            if i > 0 and lines[i - 1] == "CB ;":
                # Negative angle dipole — verify there's a trailing CB too
                assert lines[i + 1] == "CB ;", f"Missing trailing CB after DIL at line {i}"


def test_undulators_as_drifts(fox_output):
    """Undulators should appear as DL drifts."""
    dl_lines = _extract_commands(fox_output, r'DL')
    und_drifts = [l for l in dl_lines if "0.5404999999999998" in l]
    assert len(und_drifts) == 2


def test_override_ke():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, ke_override=35)
    assert "RPE 35 ;" in fox


def test_override_dimensions():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, dimensions=2)
    assert "OV 3 2 0 ;" in fox
    assert "VARIABLE F0 100 4" in fox


def test_override_order():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, computation_order=5)
    assert "OV 5 3 0 ;" in fox


@pytest.mark.skipif(not os.path.exists(REFERENCE_FOX),
                     reason="Reference input.fox not available")
def test_reference_drifts(fox_output):
    """Compare drift values against the FELsim reference (overlapping region)."""
    with open(REFERENCE_FOX) as f:
        ref_text = f.read()

    ref_dls = _extract_commands(ref_text, r'DL')
    out_dls = _extract_commands(fox_output, r'DL')

    # Compare the first N drifts that should match (before the YAML/Excel divergence)
    # First 42 drifts match (ref has one split drift at position 43)
    for i in range(min(42, len(ref_dls), len(out_dls))):
        ref_val = float(ref_dls[i].split()[1])
        out_val = float(out_dls[i].split()[1])
        assert out_val == pytest.approx(ref_val, rel=1e-12), \
            f"Drift mismatch at index {i}: {out_val} != {ref_val}"
