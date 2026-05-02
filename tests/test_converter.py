"""Tests for pals2cosy COSYScript converter.

Tests for pals2cosy COSYScript converter.
"""

import os
import re
import warnings
import pytest

from pals2cosy.parser import parse_lattice
from pals2cosy.converter import convert, _consolidate_dipoles, _format_enge

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
UHFEL_YAML = os.path.join(FIXTURES_DIR, "uhfel_beamline.yaml")


@pytest.fixture
def uhfel():
    return parse_lattice(UHFEL_YAML)


@pytest.fixture
def fox_output(uhfel):
    bp, elems = uhfel
    return convert(bp, elems, ke_override=40, fringe_field_order=0,
                   computation_order=3, dimensions=3)


def _extract_commands(fox_text, cmd_pattern):
    """Extract COSYScript commands matching a regex pattern, stripping inline comments."""
    results = []
    for line in fox_text.splitlines():
        if re.match(r'\s+' + cmd_pattern, line):
            # Strip trailing COSYScript comments: { ... }
            cmd = re.sub(r'\s*\{[^}]*\}\s*$', '', line.strip())
            results.append(cmd)
    return results


# --- Structure tests ---

def test_fox_header(fox_output):
    assert "INCLUDE 'COSY' ;" in fox_output
    assert "PROCEDURE RUN ;" in fox_output
    assert "PROCEDURE LATTICE ;" in fox_output


def test_fox_variables_no_twiss(fox_output):
    """Without --twiss, Twiss variables should not appear."""
    assert "VARIABLE A0" not in fox_output
    assert "GT MAP" not in fox_output


def test_fox_variables_with_twiss():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, ke_override=40, twiss=True)
    assert "VARIABLE A0 100 3" in fox
    assert "VARIABLE F0 100 6" in fox
    assert "GT MAP F0 MU0 A0 B0 G0 R0 ;" in fox
    assert '"beta_x"' in fox


def test_fox_settings(fox_output):
    assert "OV 3 3 0 ;" in fox_output
    assert "RPE 40 ;" in fox_output
    assert "FR 0 ;" in fox_output


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


def test_standalone_dph_consolidation():
    """Standalone DPH with edge angles becomes DIPOLE_CONSOLIDATED."""
    elements = [{
        "type": "DPH", "name": "B1", "length": 0.5,
        "angle": 15.0, "entrance_edge_angle": 3.0, "exit_edge_angle": 5.0,
        "pole_gap": 0.02, "dipole_length": 0.45,
        "s_start": 0.0, "s_end": 0.5,
        "current": None, "wedge_angle": None, "enge_coeffs": None, "bn1": None,
    }]
    consolidated = _consolidate_dipoles(elements)
    assert len(consolidated) == 1
    d = consolidated[0]
    assert d["type"] == "DIPOLE_CONSOLIDATED"
    assert d["length"] == pytest.approx(0.45)  # uses dipole_length
    assert d["angle"] == pytest.approx(15.0)
    assert d["entrance_angle"] == pytest.approx(3.0)
    assert d["exit_angle"] == pytest.approx(5.0)
    assert d["pole_gap"] == pytest.approx(0.02)


def test_standalone_dph_enge_symmetric():
    """Standalone DPH with Enge gets same coefficients on both faces."""
    elements = [{
        "type": "DPH", "name": "B1", "length": 0.5,
        "angle": 15.0, "entrance_edge_angle": 3.0, "exit_edge_angle": 5.0,
        "pole_gap": 0.02, "dipole_length": 0.45,
        "s_start": 0.0, "s_end": 0.5,
        "current": None, "wedge_angle": None,
        "enge_coeffs": [1.0, 2.0, 3.0], "bn1": None,
    }]
    consolidated = _consolidate_dipoles(elements)
    d = consolidated[0]
    assert d["entrance_enge"] == [1.0, 2.0, 3.0]
    assert d["exit_enge"] == [1.0, 2.0, 3.0]


def test_standalone_dph_no_enge():
    """Standalone DPH without Enge has None on both faces."""
    elements = [{
        "type": "DPH", "name": "B1", "length": 0.5,
        "angle": 15.0, "entrance_edge_angle": 0.0, "exit_edge_angle": 0.0,
        "pole_gap": 0.02, "dipole_length": 0.45,
        "s_start": 0.0, "s_end": 0.5,
        "current": None, "wedge_angle": None,
        "enge_coeffs": None, "bn1": None,
    }]
    consolidated = _consolidate_dipoles(elements)
    d = consolidated[0]
    assert d["entrance_enge"] is None
    assert d["exit_enge"] is None


# --- COSYScript command generation tests ---

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


def test_fc_commands_fr0(fox_output):
    """FC commands should NOT be emitted when FR=0."""
    fc_lines = _extract_commands(fox_output, r'FC')
    assert len(fc_lines) == 0


def test_fc_commands_fr1():
    """FC commands should match reference Enge coefficients when FR=1."""
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, ke_override=40, fringe_field_order=1)
    fc_lines = _extract_commands(fox, r'FC')
    expected = "FC 1 1 1 56.49 -50.79 19.32 -3.621 0.3315 -0.01193 ;"
    assert len(fc_lines) > 0
    for fc in fc_lines:
        assert fc == expected


def test_fc_commands_fr2():
    """FC commands should NOT be emitted for FR=2 (symplectic scaling uses SYSCA.DAT)."""
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, ke_override=40, fringe_field_order=2)
    fc_lines = _extract_commands(fox, r'FC')
    assert len(fc_lines) == 0
    assert "FR 2 ;" in fox


def test_fc_commands_fr3():
    """FC commands SHOULD be emitted for FR=3 (high-precision ODE uses Enge function)."""
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, ke_override=40, fringe_field_order=3)
    fc_lines = _extract_commands(fox, r'FC')
    assert len(fc_lines) > 0
    assert "FR 3 ;" in fox


def test_no_explicit_fd_in_setup(fox_output):
    """No explicit FD in setup — OV calls FD internally."""
    lines = fox_output.splitlines()
    for i, l in enumerate(lines):
        if l.strip() == "LATTICE ;":
            # Check lines between "FR" and "LATTICE" — no FD should appear
            setup_lines = [x.strip() for x in lines[:i]]
            assert "FD ;" not in setup_lines
            break


def test_cb_wrapping(fox_output):
    """Negative-angle dipoles should be wrapped with CB, positive should not."""
    lines = [l.strip() for l in fox_output.splitlines()]

    def _strip_comment(s):
        return re.sub(r'\s*\{[^}]*\}\s*$', '', s)

    # The reference beamline has negative-angle dipoles (DC2, DC4 at -4°).
    # After consolidation, all angles in DIL are abs(), so we track which
    # DIL commands appear between CB pairs.
    cb_wrapped_count = 0
    non_wrapped_count = 0
    for i, line in enumerate(lines):
        if not line.startswith("DIL"):
            continue
        has_leading_cb = i > 0 and _strip_comment(lines[i - 1]) == "CB ;"
        if has_leading_cb:
            assert _strip_comment(lines[i + 1]) == "CB ;", \
                f"Missing trailing CB after DIL at line {i}"
            cb_wrapped_count += 1
        else:
            non_wrapped_count += 1

    # UH FEL has 8 negative-angle dipoles and 10 positive (18 total)
    assert cb_wrapped_count == 8, f"Expected 8 CB-wrapped DIL, got {cb_wrapped_count}"
    assert non_wrapped_count == 10, f"Expected 10 non-wrapped DIL, got {non_wrapped_count}"


def test_undulators_as_drifts(fox_output):
    """Undulators should appear as DL drifts."""
    dl_lines = _extract_commands(fox_output, r'DL')
    und_drifts = [l for l in dl_lines if "0.5404999999999998" in l]
    assert len(und_drifts) == 2


def test_override_ke():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, ke_override=35)
    assert "RPE 35 ;" in fox


def test_particle_electron_default():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, ke_override=40)
    assert "RPE 40 ;" in fox


def test_particle_proton_override():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, ke_override=250, particle_override="proton")
    assert "RPP 250 ;" in fox


def test_comments_enabled(fox_output):
    """Element name comments should appear when comments=True (default)."""
    assert "{ LIN_QPF_004 }" in fox_output
    assert "{ DPHa }" in fox_output


def test_comments_disabled():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, comments=False)
    assert "{ LIN_QPF_004 }" not in fox
    assert "{ DPHa }" not in fox


def test_override_dimensions():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, dimensions=2)
    assert "OV 3 2 0 ;" in fox

    fox_twiss = convert(bp, elems, dimensions=2, twiss=True)
    assert "VARIABLE F0 100 4" in fox_twiss


def test_override_order():
    bp, elems = parse_lattice(UHFEL_YAML)
    fox = convert(bp, elems, computation_order=5)
    assert "OV 5 3 0 ;" in fox


def test_bn1_quad_fox():
    """MQ command uses Bn1 directly when present (no current_a computation)."""
    from pals2cosy.pals_parser import parse_lattice as parse_pals
    pals_file = os.path.join(EXAMPLES_DIR, "uhfel_excerpt.pals.yaml")
    bp, elems = parse_pals(pals_file, ke_override=40, particle_override="electron")
    fox = convert(bp, elems, ke_override=40, particle_override="electron")
    mq_lines = _extract_commands(fox, r'MQ')
    assert len(mq_lines) == 3  # LQ1, LQ2, DPHQ
    # LQ1: Bn1=-0.03221273, length=0.0889
    parts = mq_lines[0].split()
    assert float(parts[2]) == pytest.approx(-0.03221273, rel=1e-6)
    # LQ2: Bn1=0.03829256
    parts = mq_lines[1].split()
    assert float(parts[2]) == pytest.approx(0.03829256, rel=1e-6)


def test_fodo_fox():
    """End-to-end: official PALS FODO → COSYScript with Bn1 quads."""
    from pals2cosy.pals_parser import parse_lattice as parse_pals
    fodo_file = os.path.join(EXAMPLES_DIR, "fodo.pals.yaml")
    bp, elems = parse_pals(fodo_file, ke_override=1000, particle_override="proton")
    fox = convert(bp, elems, ke_override=1000, particle_override="proton")
    assert "RPP 1000 ;" in fox
    mq_lines = _extract_commands(fox, r'MQ')
    # 6 quads (2 per cell × 3 cells)
    assert len(mq_lines) == 6
    # First quad: Bn1=1.0
    parts = mq_lines[0].split()
    assert float(parts[2]) == pytest.approx(1.0)
    # Second quad: Bn1=-1.0
    parts = mq_lines[1].split()
    assert float(parts[2]) == pytest.approx(-1.0)


def test_uhfel_excerpt_dipole():
    """UH FEL excerpt SBend produces correct DIL command."""
    from pals2cosy.pals_parser import parse_lattice as parse_pals
    pals_file = os.path.join(EXAMPLES_DIR, "uhfel_excerpt.pals.yaml")
    bp, elems = parse_pals(pals_file, ke_override=40, particle_override="electron")
    fox = convert(bp, elems, ke_override=40, fringe_field_order=0)
    dil_lines = _extract_commands(fox, r'DIL')
    assert len(dil_lines) == 1
    # DC1_B1: angle ≈ 1.5°
    parts = dil_lines[0].split()
    angle = float(parts[2])
    assert angle == pytest.approx(1.5, rel=0.02)  # g_ref×L ≈ 1.501°


def test_beamline_total_length(fox_output):
    """Total beamline length from COSYScript (sum of DL+MQ+DIL) should match YAML."""
    dl_lengths = [float(c.split()[1]) for c in _extract_commands(fox_output, r'DL')]
    mq_lengths = [float(c.split()[1]) for c in _extract_commands(fox_output, r'MQ')]
    dil_lengths = [float(c.split()[1]) for c in _extract_commands(fox_output, r'DIL')]
    total = sum(dl_lengths) + sum(mq_lengths) + sum(dil_lengths)
    # UH FEL beamline total: 14.400 m (consolidated, excluding DPW wedge lengths)
    assert total == pytest.approx(14.4, rel=1e-3)


def test_element_counts(fox_output):
    """Verify expected element counts in generated COSYScript."""
    dl_count = len(_extract_commands(fox_output, r'DL'))
    mq_count = len(_extract_commands(fox_output, r'MQ'))
    dil_count = len(_extract_commands(fox_output, r'DIL'))
    assert mq_count == 26   # 13 QPF + 13 QPD
    assert dil_count == 18   # 10 chicane + 8 spectrometer
    assert dl_count > 40     # drifts + undulators


# --- P2-4: Converter edge cases ---

def test_lone_dpw_not_consolidated():
    """A lone DPW (not part of DPW-DPH-DPW triplet) passes through unconsolidated."""
    elements = [
        {"type": "DPW", "name": "W1", "length": 0.01, "wedge_angle": 5.0,
         "angle": 15.0, "pole_gap": 0.01, "dipole_length": 0.2,
         "enge_coeffs": None, "s_start": 0.0, "s_end": 0.01},
        {"type": "DRIFT", "name": "", "length": 0.1,
         "s_start": 0.01, "s_end": 0.11},
    ]
    consolidated = _consolidate_dipoles(elements)
    types = [e["type"] for e in consolidated]
    assert "DPW" in types
    assert "DIPOLE_CONSOLIDATED" not in types


def test_comment_sanitization():
    """Element name containing } doesn't break COSYScript comment delimiters."""
    elements = [
        {"type": "DRIFT", "name": "bad}name{here", "length": 0.5,
         "s_start": 0.0, "s_end": 0.5},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    # Name should appear sanitized: no nested braces
    assert "{ badnamehere }" in fox
    # Verify balanced braces in comments
    for line in fox.splitlines():
        if "{" in line and "WRITE" not in line:
            assert line.count("{") == line.count("}"), f"Unbalanced braces: {line}"


def test_lone_dpw_fox_output():
    """A lone DPW emits a drift + warning in COSYScript output."""
    elements = [
        {"type": "DPW", "name": "W1", "length": 0.01, "wedge_angle": 5.0,
         "angle": 15.0, "pole_gap": 0.01, "dipole_length": 0.2,
         "enge_coeffs": None, "s_start": 0.0, "s_end": 0.01},
        {"type": "DRIFT", "name": "", "length": 0.5,
         "s_start": 0.01, "s_end": 0.51},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        fox = convert(bp, elements)
    assert any("Standalone DPW" in str(warning.message) for warning in w)
    dl_lines = _extract_commands(fox, r'DL')
    # Both the DPW-as-drift and the actual drift should appear
    assert len(dl_lines) == 2
    assert "0.01" in dl_lines[0]


def test_enge_truncation_warning():
    """Enge coefficients longer than 6 emit a warning."""
    coeffs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = _format_enge(coeffs)
    assert any("truncated" in str(warning.message) for warning in w)
    # Still produces 6 values
    assert len(result.split()) == 6


def test_enge_padding():
    """Fewer than 6 Enge coefficients are zero-padded."""
    result = _format_enge([1.0, 2.0])
    vals = result.split()
    assert len(vals) == 6
    assert vals[0] == "1.0"
    assert vals[1] == "2.0"
    assert all(v == "0.0" for v in vals[2:])


def test_enge_empty():
    """Empty Enge coefficient list produces all-zero output."""
    result = _format_enge([])
    vals = result.split()
    assert len(vals) == 6
    assert all(v == "0.0" for v in vals)


# --- Phase 2A: Robustness edge case tests ---

def test_muon_particle():
    """Muon particle type generates RPMU command."""
    bp = {"kinetic_energy_mev": 200, "particle_type": "muon",
          "quad_gradient": 2.694, "source_file": "test"}
    elements = [{"type": "DRIFT", "name": "", "length": 0.5,
                 "s_start": 0.0, "s_end": 0.5}]
    fox = convert(bp, elements, particle_override="muon")
    assert "RPMU 200 ;" in fox


def test_pion_particle():
    """Pion particle type generates RPPI command with charge."""
    bp = {"kinetic_energy_mev": 300, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    elements = [{"type": "DRIFT", "name": "", "length": 0.5,
                 "s_start": 0.0, "s_end": 0.5}]
    fox = convert(bp, elements, particle_override="pion+")
    assert "RPPI 300 1 ;" in fox
    fox2 = convert(bp, elements, particle_override="pi-")
    assert "RPPI 300 -1 ;" in fox2


def test_antimuon_uses_rp_not_rpmu():
    """μ⁺ must use generic RP with charge +1, not RPMU which hardcodes charge -1."""
    bp = {"kinetic_energy_mev": 200, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    elements = [{"type": "DRIFT", "name": "", "length": 0.5,
                 "s_start": 0.0, "s_end": 0.5}]
    for name in ("mu+", "antimuon"):
        fox = convert(bp, elements, particle_override=name)
        assert "RP 200 0.1134289168 1 ;" in fox
        assert "RPMU" not in fox


def test_generic_particle():
    """Generic particles (positron, deuteron, etc.) use RP command."""
    bp = {"kinetic_energy_mev": 500, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    elements = [{"type": "DRIFT", "name": "", "length": 0.5,
                 "s_start": 0.0, "s_end": 0.5}]
    fox = convert(bp, elements, particle_override="positron")
    assert "RP 500 0.000548579911 1 ;" in fox
    fox2 = convert(bp, elements, particle_override="deuteron")
    assert "RP 500 2.01355321271 1 ;" in fox2


def test_particle_case_insensitive():
    """Particle names are case-insensitive."""
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    elements = [{"type": "DRIFT", "name": "", "length": 0.5,
                 "s_start": 0.0, "s_end": 0.5}]
    fox1 = convert(bp, elements, particle_override="Electron")
    fox2 = convert(bp, elements, particle_override="PROTON")
    assert "RPE 40 ;" in fox1
    assert "RPP 40 ;" in fox2


def test_unknown_particle_raises():
    """convert() with an unrecognized particle type raises ValueError."""
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    elements = [{"type": "DRIFT", "name": "", "length": 0.5,
                 "s_start": 0.0, "s_end": 0.5}]
    with pytest.raises(ValueError, match="Unknown particle type"):
        convert(bp, elements, particle_override="tachyon")


def test_sextupole_with_bn2():
    """Sextupole with Bn2 generates MH command."""
    elements = [
        {"type": "SXT", "name": "S1", "length": 0.2,
         "s_start": 0.0, "s_end": 0.2, "bn2": 0.5,
         "bn1": None, "bz": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "MH 0.2 0.5 0.0135 ;" in fox


def test_sextupole_without_bn2():
    """Sextupole without Bn2 falls back to drift."""
    elements = [
        {"type": "SXT", "name": "S1", "length": 0.2,
         "s_start": 0.0, "s_end": 0.2, "bn2": None,
         "bn1": None, "bz": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "MH" not in fox
    assert "DL 0.2" in fox


def test_solenoid_with_bz():
    """Solenoid with Bz generates CMS command."""
    elements = [
        {"type": "SOL", "name": "SOL1", "length": 0.5,
         "s_start": 0.0, "s_end": 0.5, "bz": 1.5,
         "bn1": None, "bn2": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "CMS 1.5 0.0135 0.5 ;" in fox


def test_octupole_with_bn3():
    """Octupole with Bn3 generates MO command."""
    elements = [
        {"type": "OCT", "name": "O1", "length": 0.15,
         "s_start": 0.0, "s_end": 0.15, "bn3": 0.3,
         "bn1": None, "bn2": None, "bz": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None, "rf_voltage_kv": None, "rf_frequency_hz": None,
         "rf_phase_deg": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "MO 0.15 0.3 0.0135 ;" in fox


def test_octupole_without_bn3():
    """Octupole without Bn3 falls back to drift."""
    elements = [
        {"type": "OCT", "name": "O1", "length": 0.15,
         "s_start": 0.0, "s_end": 0.15, "bn3": None,
         "bn1": None, "bn2": None, "bz": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None, "rf_voltage_kv": None, "rf_frequency_hz": None,
         "rf_phase_deg": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "MO" not in fox
    assert "DL 0.15" in fox


def test_rf_cavity_with_params():
    """RF cavity with voltage and frequency generates RF command."""
    elements = [
        {"type": "RFC", "name": "RFC1", "length": 0.5,
         "s_start": 0.0, "s_end": 0.5,
         "rf_voltage_kv": 100.0, "rf_frequency_hz": 2856e6,
         "rf_phase_deg": 90.0,
         "bn1": None, "bn2": None, "bn3": None, "bz": None,
         "current": None, "angle": None, "wedge_angle": None,
         "entrance_edge_angle": None, "exit_edge_angle": None,
         "pole_gap": None, "dipole_length": None, "enge_coeffs": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "RF 100.0 0 2856000000.0 90.0 0.0135 ;" in fox


def test_rf_cavity_without_params():
    """RF cavity without voltage/frequency falls back to drift."""
    elements = [
        {"type": "RFC", "name": "RFC1", "length": 0.5,
         "s_start": 0.0, "s_end": 0.5,
         "rf_voltage_kv": None, "rf_frequency_hz": None,
         "rf_phase_deg": None,
         "bn1": None, "bn2": None, "bn3": None, "bz": None,
         "current": None, "angle": None, "wedge_angle": None,
         "entrance_edge_angle": None, "exit_edge_angle": None,
         "pole_gap": None, "dipole_length": None, "enge_coeffs": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "RF " not in fox
    assert "DL 0.5" in fox


def test_multipole_with_harmonics():
    """Multipole with Bn fields generates M5 command."""
    elements = [
        {"type": "MULT", "name": "M1", "length": 0.2,
         "s_start": 0.0, "s_end": 0.2,
         "bn1": 0.5, "bn2": 0.1, "bn3": 0.05, "bn4": None, "bn5": None,
         "bz": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None, "rf_voltage_kv": None, "rf_frequency_hz": None,
         "rf_phase_deg": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "M5 0.2 0.5 0.1 0.05 0 0 0.0135 ;" in fox


def test_multipole_without_harmonics():
    """Multipole without any Bn fields falls back to drift."""
    elements = [
        {"type": "MULT", "name": "M1", "length": 0.2,
         "s_start": 0.0, "s_end": 0.2,
         "bn1": None, "bn2": None, "bn3": None, "bn4": None, "bn5": None,
         "bz": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None, "rf_voltage_kv": None, "rf_frequency_hz": None,
         "rf_phase_deg": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "M5" not in fox
    assert "DL 0.2" in fox


def test_wiggler_with_field():
    """Wiggler with peak field and period generates WI command."""
    import math
    elements = [
        {"type": "UND", "name": "U1", "length": 2.0,
         "s_start": 0.0, "s_end": 2.0,
         "wiggler_field": 0.8, "wiggler_period": 0.05,
         "bn1": None, "bn2": None, "bn3": None, "bn4": None, "bn5": None,
         "bz": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None, "rf_voltage_kv": None, "rf_frequency_hz": None,
         "rf_phase_deg": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    expected_k = 2 * math.pi / 0.05
    assert f"WI 0.8 {expected_k} 2.0 0.0135 0 0 0 ;" in fox


def test_wiggler_without_field():
    """Wiggler without field parameters falls back to drift."""
    elements = [
        {"type": "UND", "name": "U1", "length": 2.0,
         "s_start": 0.0, "s_end": 2.0,
         "wiggler_field": None, "wiggler_period": None,
         "bn1": None, "bn2": None, "bn3": None, "bn4": None, "bn5": None,
         "bz": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None, "rf_voltage_kv": None, "rf_frequency_hz": None,
         "rf_phase_deg": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "WI" not in fox
    assert "DL 2.0" in fox


def test_solenoid_without_bz():
    """Solenoid without Bz falls back to drift."""
    elements = [
        {"type": "SOL", "name": "SOL1", "length": 0.5,
         "s_start": 0.0, "s_end": 0.5, "bz": None,
         "bn1": None, "bn2": None, "current": None, "angle": None,
         "wedge_angle": None, "entrance_edge_angle": None,
         "exit_edge_angle": None, "pole_gap": None, "dipole_length": None,
         "enge_coeffs": None},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "CMS" not in fox
    assert "DL 0.5" in fox


def test_exit_only_enge():
    """DIPOLE_CONSOLIDATED with exit Enge only emits FC 1 2 1 but not FC 1 1 1."""
    elements = [{
        "type": "DIPOLE_CONSOLIDATED",
        "name": "D1", "length": 0.1, "angle": 5.0,
        "entrance_angle": 0.0, "exit_angle": 0.0,
        "pole_gap": 0.01,
        "entrance_enge": None,
        "exit_enge": [1.0, 2.0, 3.0],
    }]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements, fringe_field_order=1)
    assert "FC 1 2 1" in fox
    assert "FC 1 1 1" not in fox


def test_negative_angle_with_enge():
    """Negative-angle dipole with Enge: FC → CB → DIL → CB → FD ordering."""
    elements = [{
        "type": "DIPOLE_CONSOLIDATED",
        "name": "D1", "length": 0.1, "angle": -5.0,
        "entrance_angle": 2.0, "exit_angle": 2.0,
        "pole_gap": 0.01,
        "entrance_enge": [1.0, 2.0, 3.0],
        "exit_enge": None,
    }]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements, fringe_field_order=1)
    lines = [l.strip() for l in fox.splitlines() if l.strip()]
    # Find FC, CB, DIL, CB, FD in order
    fc_idx = next(i for i, l in enumerate(lines) if l.startswith("FC 1 1 1"))
    cb1_idx = next(i for i, l in enumerate(lines) if i > fc_idx and l.startswith("CB"))
    dil_idx = next(i for i, l in enumerate(lines) if i > cb1_idx and l.startswith("DIL"))
    cb2_idx = next(i for i, l in enumerate(lines) if i > dil_idx and l.startswith("CB"))
    fd_idx = next(i for i, l in enumerate(lines) if i > cb2_idx and l.startswith("FD"))
    assert fc_idx < cb1_idx < dil_idx < cb2_idx < fd_idx


def test_unknown_element_type_warns():
    """Unrecognized element type emits a warning."""
    elements = [{"type": "FOOBAR", "name": "X1", "length": 0.5,
                 "s_start": 0.0, "s_end": 0.5}]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        convert(bp, elements)
    assert any("Unknown element type 'FOOBAR'" in str(warning.message) for warning in w)


def test_zero_length_quad_skipped():
    """A zero-length quadrupole produces no MQ command."""
    elements = [
        {"type": "QPF", "name": "Q0", "length": 0.0,
         "s_start": 0.0, "s_end": 0.0, "bn1": 0.5,
         "current": None, "angle": None, "wedge_angle": None,
         "entrance_edge_angle": None, "exit_edge_angle": None,
         "pole_gap": None, "dipole_length": None, "enge_coeffs": None},
        {"type": "DRIFT", "name": "", "length": 0.5,
         "s_start": 0.0, "s_end": 0.5},
    ]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    mq_lines = [l.strip() for l in fox.splitlines() if l.strip().startswith("MQ")]
    assert len(mq_lines) == 0


def test_zero_length_dipole_skipped():
    """A zero-length dipole produces no DIL command."""
    elements = [{
        "type": "DIPOLE_CONSOLIDATED",
        "name": "D0", "length": 0.0, "angle": 5.0,
        "entrance_angle": 0.0, "exit_angle": 0.0,
        "pole_gap": 0.01,
        "entrance_enge": None, "exit_enge": None,
    }]
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": 2.694, "source_file": "test"}
    fox = convert(bp, elements)
    assert "DIL" not in fox


def test_quad_gradient_none():
    """convert() handles quad_gradient=None (from pals_parser) without crashing."""
    bp = {"kinetic_energy_mev": 40, "particle_type": "electron",
          "quad_gradient": None, "source_file": "test"}
    elements = [
        {"type": "QPF", "name": "Q1", "length": 0.1,
         "s_start": 0.0, "s_end": 0.1, "bn1": None,
         "current": 1.5, "angle": None, "wedge_angle": None,
         "entrance_edge_angle": None, "exit_edge_angle": None,
         "pole_gap": None, "dipole_length": None, "enge_coeffs": None},
    ]
    fox = convert(bp, elements)
    mq_lines = [l.strip() for l in fox.splitlines() if l.strip().startswith("MQ")]
    assert len(mq_lines) == 1
    # Should use G_QUAD_DEFAULT=2.694 as fallback
    parts = mq_lines[0].split()
    expected = -2.694 * 1.5 * 0.0135
    assert float(parts[2]) == pytest.approx(expected)
