"""Physics of the PALS → COSYScript mapping, checked without running COSY."""

import math
import re

import pytest

from pals2cosy.constants import C_LIGHT, E0_PROTON
from pals2cosy.converter import THIN_SLICE, convert, pole_tip_field
from pals2cosy.pals_elements import (
    cosy_rf_phase, resolve_bend_geometry, resolve_edge_angles,
)
from pals2cosy.pals_parser import parse_lattice

R = 0.0135  # converter default bore radius


def _lattice(elements, species="proton", energy="pc_ref: 1.0e9"):
    body = "\n".join("          " + line for line in elements.strip().splitlines())
    return f"""\
PALS:
  facility:
    - main:
        kind: BeamLine
        line:
          - begin:
              kind: BeginningEle
              ReferenceP:
                species_ref: {species}
                {energy}
{body}
"""


def _commands(fox, cmd):
    out = []
    for line in fox.splitlines():
        line = re.sub(r"\s*\{[^}]*\}\s*$", "", line.strip())
        if line.split(" ", 1)[0] == cmd:
            out.append(line.split()[1:-1])
    return out


def _one(fox, cmd):
    found = _commands(fox, cmd)
    assert len(found) == 1, fox
    return [float(x) for x in found[0]]


@pytest.fixture
def convert_text(pals_file):
    def run(text, **kwargs):
        bp, elements = parse_lattice(pals_file(text))
        return bp, elements, convert(bp, elements, **kwargs)
    return run


# --- Reference particle and rigidity ---------------------------------------

def test_brho_one_gev_proton(convert_text):
    bp, _, _ = convert_text(_lattice("- d: {kind: Drift, length: 1.0}"))
    assert bp["pc_mev"] == pytest.approx(1000.0)
    assert bp["brho_tm"] == pytest.approx(1e9 / C_LIGHT, rel=1e-12)
    assert bp["brho_tm"] == pytest.approx(3.3356, abs=5e-5)
    assert bp["kinetic_energy_mev"] == pytest.approx(
        math.hypot(1000.0, E0_PROTON) - E0_PROTON)


def test_brho_is_signed_by_charge(convert_text):
    bp, _, _ = convert_text(_lattice("- d: {kind: Drift, length: 1.0}", "electron"))
    assert bp["brho_tm"] == pytest.approx(-1e9 / C_LIGHT, rel=1e-12)


def test_total_energy_reference(convert_text):
    bp, _, fox = convert_text(_lattice("- d: {kind: Drift, length: 1.0}", "electron",
                                       "E_tot_ref: 50.51099895e6"))
    assert bp["kinetic_energy_mev"] == pytest.approx(50.0)
    assert "RPE 50" in fox


def test_inconsistent_energy_and_momentum(pals_file):
    text = _lattice("- d: {kind: Drift, length: 1.0}", "proton",
                    "E_tot_ref: 2.0e9\n                pc_ref: 1.0e9")
    with pytest.raises(ValueError, match="inconsistent"):
        parse_lattice(pals_file(text))


def test_cli_overrides_win(pals_file):
    bp, _ = parse_lattice(pals_file(_lattice("- d: {kind: Drift, length: 1.0}")),
                          ke_override=40, particle_override="electron")
    assert bp["particle_type"] == "electron"
    assert bp["kinetic_energy_mev"] == 40


# --- Multipole field scaling ------------------------------------------------

@pytest.mark.parametrize("order", [1, 2, 3, 4, 5])
def test_pole_tip_formula(order):
    assert pole_tip_field(7.0, order, 0.02) == pytest.approx(
        7.0 * 0.02 ** order / math.factorial(order))


@pytest.mark.parametrize("kind, key, cmd, index", [
    ("Quadrupole", "Bn1", "MQ", 1),
    ("Sextupole", "Bn2", "MH", 1),
    ("Octupole", "Bn3", "MO", 1),
])
def test_single_order_pole_tip(convert_text, kind, key, cmd, index):
    order = int(key[-1])
    _, _, fox = convert_text(_lattice(
        f"- m: {{kind: {kind}, length: 0.2, MagneticMultipoleP: {{{key}: 3.0}}}}"))
    args = _one(fox, cmd)
    assert args[0] == pytest.approx(0.2)
    assert args[index] == pytest.approx(3.0 * R ** order / math.factorial(order))
    assert args[-1] == pytest.approx(R)


def test_multipole_pole_tips_through_dodecapole(convert_text):
    _, _, fox = convert_text(_lattice(
        "- m: {kind: Multipole, length: 0.3, MagneticMultipoleP: "
        "{Bn1: 1.0, Bn2: 2.0, Bn3: 3.0, Bn4: 4.0, Bn5: 5.0}}"))
    args = _one(fox, "M5")
    expected = [n * R ** n / math.factorial(n) for n in range(1, 6)]
    assert args[1:6] == pytest.approx(expected)


def test_quadrupole_aperture_does_not_change_gradient(convert_text):
    _, _, fox = convert_text(_lattice(
        "- q: {kind: Quadrupole, length: 0.2, MagneticMultipoleP: {Bn1: 3.0}}"),
        quad_aperture=0.1)
    _, b, r = _one(fox, "MQ")
    assert b / r == pytest.approx(3.0)


def test_normalized_strength_uses_rigidity(convert_text):
    bp, elements, fox = convert_text(_lattice(
        "- q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {Kn1: 0.5}}"))
    assert elements[0]["bn"][1] == pytest.approx(0.5 * bp["brho_tm"])
    assert _one(fox, "MQ")[1] == pytest.approx(0.5 * bp["brho_tm"] * R)


def test_normalized_strength_sign_follows_charge(convert_text):
    bp, elements, _ = convert_text(_lattice(
        "- q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {Kn1: 0.5}}",
        "electron"))
    assert elements[0]["bn"][1] == pytest.approx(-0.5 * 1e9 / C_LIGHT)
    assert elements[0]["type"] == "QPD"


def test_integrated_strength_divided_by_length(convert_text):
    bp, elements, _ = convert_text(_lattice(
        "- q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {Kn1L: 0.25, Bn2L: 0.4}}"))
    assert elements[0]["bn"][1] == pytest.approx(0.5 * bp["brho_tm"])
    assert elements[0]["bn"][2] == pytest.approx(0.8)


def test_taper_adds_to_field(convert_text):
    _, elements, _ = convert_text(_lattice(
        "- q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {Bn1: 2.0, Bn1_taper: 0.1}}"))
    assert elements[0]["bn"][1] == pytest.approx(2.1)


def test_thin_multipole_keeps_length_and_kick(convert_text):
    bp, elements, fox = convert_text(_lattice(
        "- q: {kind: Quadrupole, MagneticMultipoleP: {Kn1L: 0.1}}\n"
        "- d: {kind: Drift, length: 1.0}"), fringe_field_order=1)
    assert elements[0]["thin"]
    lines = [l.strip() for l in fox.splitlines()]
    i = next(i for i, l in enumerate(lines) if l.startswith("MQ"))
    assert lines[i - 2:i] == ["FR 0 ;", f"DL {-THIN_SLICE / 2} ;"]
    assert lines[i + 1:i + 3] == [f"DL {-THIN_SLICE / 2} ;", "FR 1 ;"]
    length, b_pole, r = _one(fox, "MQ")
    assert length == THIN_SLICE
    assert b_pole / r * length == pytest.approx(0.1 * bp["brho_tm"])


def test_skew_component_uses_mms(convert_text):
    _, elements, fox = convert_text(_lattice(
        "- q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {Bs1: 2.0}}"))
    assert elements[0]["bs"] == {1: pytest.approx(2.0)}
    assert "VARIABLE PALSBN 1 10 ; VARIABLE PALSBS 1 10 ;" in fox
    assert f"PALSBN(1) := 0.0 ; PALSBS(1) := {2.0 * R} ;" in fox
    assert _commands(fox, "MMS") == [["0.5", "PALSBN", "PALSBS", "1", str(R)]]


def test_tilt_is_folded_into_skew(convert_text):
    _, elements, _ = convert_text(_lattice(
        "- q: {kind: Quadrupole, length: 0.5, "
        "MagneticMultipoleP: {Bn1: 2.0, tilt1: 0.7853981633974483}}"))
    # (Bn + i Bs) exp(-i (N+1) T) with T = pi/4: a normal quad becomes Bs = -Bn
    assert 1 not in elements[0]["bn"]
    assert elements[0]["bs"][1] == pytest.approx(-2.0)


def test_mixed_strength_forms_rejected(pals_file):
    for group in ("{Bn1: 1.0, Kn1: 0.3}", "{Kn1: 0.3, Bs1: 0.2}", "{Kn1: 0.3, Ks1L: 0.2}"):
        text = _lattice(f"- q: {{kind: Quadrupole, length: 0.5, MagneticMultipoleP: {group}}}")
        with pytest.raises(ValueError):
            parse_lattice(pals_file(text))


def test_multipole_order_above_five_uses_mms(convert_text):
    _, _, fox = convert_text(_lattice(
        "- m: {kind: Multipole, length: 0.2, MagneticMultipoleP: {Bn6: 720.0}}"))
    assert f"PALSBN(6) := {720.0 * R ** 6 / 720} ;" in fox
    assert _commands(fox, "MMS")[0][3] == "6"


def test_missing_strength_warns(pals_file):
    text = _lattice("- q: {kind: Quadrupole, length: 0.5}")
    with pytest.warns(UserWarning, match="no readable strength"):
        _, elements = parse_lattice(pals_file(text))
    assert elements[0]["bn"] == {}


# --- Bends --------------------------------------------------------------------

RHO, THETA = 5.0, 0.2


@pytest.mark.parametrize("given", [
    {"g": 1 / RHO, "length": RHO * THETA},
    {"g": 1 / RHO, "chord": 2 * RHO * math.sin(THETA / 2)},
    {"g": 1 / RHO, "rect": RHO * math.sin(THETA)},
    {"angle": THETA, "length": RHO * THETA},
    {"angle": THETA, "chord": 2 * RHO * math.sin(THETA / 2)},
    {"angle": THETA, "rect": RHO * math.sin(THETA)},
    {"g": 1 / RHO, "angle": THETA},
    {"g": 1 / RHO, "angle": THETA, "length": RHO * THETA},
])
def test_bend_two_of_three(given):
    L, g, angle = resolve_bend_geometry("b", **given)
    assert L == pytest.approx(RHO * THETA)
    assert g == pytest.approx(1 / RHO)
    assert angle == pytest.approx(THETA)


@pytest.mark.parametrize("given", [
    {"g": 0.2, "angle": 0.3, "length": 1.0},       # over-specified, inconsistent
    {"length": 1.0, "chord": 1.0},                  # two from the same set
    {"angle": 0.2},                                 # angle alone
    {"g": 1.0, "chord": 3.0},                       # chord longer than diameter
])
def test_bend_geometry_errors(given):
    with pytest.raises(ValueError):
        resolve_bend_geometry("b", **given)


def test_bend_from_radius_and_angle(convert_text):
    _, elements, fox = convert_text(_lattice(
        "- b: {kind: Bend, BendP: {angle_ref: 0.2, radius_ref: 5.0}}"))
    assert elements[0]["length"] == pytest.approx(1.0)
    assert _one(fox, "DIL")[:2] == pytest.approx([1.0, math.degrees(0.2)])


def test_bend_from_reference_field(convert_text):
    bp, elements, _ = convert_text(_lattice(
        "- b: {kind: Bend, length: 1.0, BendP: {Bn0_ref: 0.5}}"))
    assert math.radians(elements[0]["angle"]) == pytest.approx(0.5 / bp["brho_tm"])


@pytest.mark.parametrize("geometry, e1_rect, e2_rect, expected", [
    ("ARC", 0.01, 0.02, (0.01 + THETA / 2, 0.02 + THETA / 2)),
    ("CHORD", None, None, (THETA / 2, THETA / 2)),
    ("ENTRANCE_COORDS", None, None, (0.0, THETA)),
    ("EXIT_COORDS", 0.0, 0.0, (THETA, 0.0)),
    ("ARC", None, None, (0.0, 0.0)),
])
def test_edge_angles_by_ref_geometry(geometry, e1_rect, e2_rect, expected):
    assert resolve_edge_angles("b", THETA, geometry, None, None,
                               e1_rect, e2_rect) == pytest.approx(expected)


def test_edge_angle_consistency_check():
    with pytest.raises(ValueError, match="inconsistent"):
        resolve_edge_angles("b", THETA, "ARC", 0.5, None, 0.0, None)


def test_negative_bend_edges_in_cb_frame(convert_text):
    # A rectangular bend with a negative angle has e1 = e2 = angle/2 < 0; under
    # COSY's CB the mirrored magnet needs positive face angles.
    _, _, fox = convert_text(_lattice(
        "- b: {kind: Bend, length: 1.0, BendP: {angle_ref: -0.2, ref_geometry: CHORD}}"))
    lines = [l.strip() for l in fox.splitlines()]
    i = next(i for i, l in enumerate(lines) if l.startswith("DIL"))
    assert lines[i - 1] == "CB ;" and lines[i + 1] == "CB ;"
    L, angle, d, e1, h1, e2, h2 = _one(fox, "DIL")
    assert angle == pytest.approx(math.degrees(0.2))
    assert e1 == e2 == pytest.approx(math.degrees(0.1))


def test_bend_tilt_and_face_curvature(convert_text):
    _, _, fox = convert_text(_lattice(
        "- b: {kind: Bend, length: 1.0, BendP: {angle_ref: 0.2, tilt_ref: 1.5707963267948966,"
        " h1: 0.5, h2: -0.25}}"))
    lines = [l.strip() for l in fox.splitlines()]
    i = next(i for i, l in enumerate(lines) if l.startswith("DIL"))
    assert lines[i - 1] == "RA 90.0 ;" and lines[i + 1] == "RA -90.0 ;"
    _, _, _, _, h1, _, h2 = _one(fox, "DIL")
    assert (h1, h2) == (0.5, -0.25)


def test_bend_pole_gap_from_aperture(convert_text):
    _, _, fox = convert_text(_lattice(
        "- b: {kind: Bend, length: 1.0, BendP: {g_ref: 0.1}, "
        "ApertureP: {y_min: -0.02, y_max: 0.02}}"))
    assert _one(fox, "DIL")[2] == pytest.approx(0.02)


def test_bend_dipole_error_warns(pals_file):
    text = _lattice("- b: {kind: Bend, length: 1.0, BendP: {g_ref: 0.1}, "
                    "MagneticMultipoleP: {Kn0: 0.12}}")
    with pytest.warns(UserWarning, match="differs from the reference field"):
        parse_lattice(pals_file(text))


# --- RF cavities --------------------------------------------------------------

def test_rf_preserves_length(convert_text):
    _, _, fox = convert_text(_lattice(
        "- c: {kind: RFCavity, length: 0.4, RFP: {voltage: 8.0e6, frequency: 1.3e9}}"))
    lines = [re.sub(r"\s*\{[^}]*\}$", "", l.strip()) for l in fox.splitlines()]
    i = lines.index("PALSRFV(1,1) := 8000.0 ;")
    assert lines[i - 1] == "DL 0.2 ;"
    assert lines[i + 1] == "RF PALSRFV 0 1300000000.0 90.0 0.0135 ;"
    assert lines[i + 2] == "DL 0.2 ;"
    drifts = sum(float(a[0]) for a in _commands(fox, "DL"))
    assert drifts == pytest.approx(0.4)


def test_rf_gradient_and_active_length(convert_text):
    _, elements, _ = convert_text(_lattice(
        "- c: {kind: RFCavity, length: 0.4, RFP: {gradient: 1.0e8, L_active: 0.25,"
        " frequency: 1.3e9, dE_ref: 2.5e7}}"))
    assert elements[0]["rf_voltage_kv"] == pytest.approx(2.5e4)


def test_rf_inconsistent_voltage_and_gradient(pals_file):
    text = _lattice("- c: {kind: RFCavity, length: 0.4, RFP: {gradient: 1.0e7,"
                    " voltage: 8.0e6, frequency: 1.3e9}}")
    with pytest.raises(ValueError, match="disagree"):
        parse_lattice(pals_file(text))


@pytest.mark.parametrize("zero_phase, charge, phase, expected", [
    ("ACCELERATING", +1, 0.0, 90.0),
    ("ACCELERATING", -1, 0.0, -90.0),
    ("ACCELERATING", +1, -0.5, 90.0 - math.degrees(0.5)),
    ("BELOW_TRANSITION", +1, 0.0, 0.0),
    ("ABOVE_TRANSITION", +1, 0.0, 180.0),
    ("BELOW_TRANSITION", -1, 0.0, 180.0),
])
def test_rf_phase_convention(zero_phase, charge, phase, expected):
    assert cosy_rf_phase(phase, zero_phase, charge) == pytest.approx(expected)


def test_rf_on_crest_gain_matches_dE_ref(pals_file):
    """COSY's reference energy gain Z V sin(PHI) equals the PALS crest gain."""
    text = _lattice("- c: {kind: RFCavity, length: 0.4, RFP: {voltage: 8.0e6,"
                    " frequency: 1.3e9, dE_ref: 8.0e6}}", "electron")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        parse_lattice(pals_file(text))


def test_rf_reference_mismatch_warns(pals_file):
    text = _lattice("- c: {kind: RFCavity, length: 0.4, RFP: {voltage: 8.0e6,"
                    " frequency: 1.3e9}}")
    with pytest.warns(UserWarning, match="re-references"):
        parse_lattice(pals_file(text))


def test_rf_harmonic_number(convert_text):
    bp, elements, _ = convert_text(_lattice(
        "- d: {kind: Drift, length: 99.0}\n"
        "- c: {kind: RFCavity, length: 1.0, RFP: {voltage: 0.0, harmonic_number: 4}}\n"
        "- c2: {kind: RFCavity, length: 0.0, RFP: {voltage: 1.0e3, harmonic_number: 4,"
        " zero_phase: BELOW_TRANSITION}}"))
    beta = bp["pc_mev"] / math.hypot(bp["pc_mev"], E0_PROTON)
    assert elements[2]["rf_frequency_hz"] == pytest.approx(4 * beta * C_LIGHT / 100.0)


def test_rf_needs_frequency(pals_file):
    text = _lattice("- c: {kind: RFCavity, length: 0.4, RFP: {voltage: 8.0e6}}")
    with pytest.raises(ValueError, match="frequency or harmonic_number"):
        parse_lattice(pals_file(text))


# --- Solenoids ----------------------------------------------------------------

def test_solenoid_field_and_normalized(convert_text):
    bp, _, fox = convert_text(_lattice(
        "- s1: {kind: Solenoid, length: 0.5, SolenoidP: {Bsol: 1.5}}\n"
        "- s2: {kind: Solenoid, length: 0.5, SolenoidP: {Ksol: 0.3}}"))
    (b1, d1, l1), (b2, d2, l2) = [[float(x) for x in c] for c in _commands(fox, "CMS")]
    assert (b1, d1, l1) == (1.5, R, 0.5)
    assert b2 == pytest.approx(0.3 * bp["brho_tm"])
