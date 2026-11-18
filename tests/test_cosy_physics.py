"""End-to-end checks: PALS lattice → COSYScript → COSY INFINITY transfer map.

Each test compares map coefficients computed by COSY with the value implied
by the PALS definition of the element. Skipped without COSY INFINITY.
"""

import math

import pytest

from pals2cosy.constants import C_LIGHT
from pals2cosy.converter import convert
from pals2cosy.pals_parser import parse_lattice
from tests.cosy_harness import requires_cosy, run_cosy

pytestmark = requires_cosy

PC = 1.0e9           # reference momentum [eV/c]
BRHO = PC / C_LIGHT  # proton rigidity [T m]


def _map(pals_file, cosy_dir, elements, species="proton", order=3, **convert_kwargs):
    body = "\n".join("          " + line for line in elements.strip().splitlines())
    path = pals_file(f"""\
PALS:
  facility:
    - main:
        kind: BeamLine
        line:
          - begin:
              kind: BeginningEle
              ReferenceP: {{species_ref: {species}, pc_ref: {PC}}}
{body}
""")
    bp, elems = parse_lattice(path)
    result = run_cosy(convert(bp, elems, computation_order=order, **convert_kwargs), cosy_dir)
    assert not result.errors, result.errors
    assert result.scalars["SPOS"] == pytest.approx(sum(e["length"] for e in elems))
    return bp, result


def test_rigidity_matches_cosy(pals_file, cosy_dir):
    _, m = _map(pals_file, cosy_dir, "- d: {kind: Drift, length: 1.0}")
    assert m.scalars["CHIM"] == pytest.approx(BRHO, rel=1e-6)


def test_thick_quadrupole_from_gradient(pals_file, cosy_dir):
    _, m = _map(pals_file, cosy_dir,
                "- q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {Bn1: 2.0}}")
    k = math.sqrt(2.0 / BRHO)
    assert m.coeff(0, 1) == pytest.approx(math.cos(k * 0.5), rel=1e-6)
    assert m.coeff(1, 1) == pytest.approx(-k * math.sin(k * 0.5), rel=1e-6)
    assert m.coeff(2, 0, 0, 1) == pytest.approx(math.cosh(k * 0.5), rel=1e-6)


def test_normalized_quadrupole_focuses_electrons(pals_file, cosy_dir):
    _, m = _map(pals_file, cosy_dir,
                "- q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {Kn1: 0.8}}",
                species="electron")
    k = math.sqrt(0.8)
    assert m.coeff(1, 1) == pytest.approx(-k * math.sin(k * 0.5), rel=1e-6)


@pytest.mark.parametrize("order, kind", [(2, "Sextupole"), (3, "Octupole"),
                                         (4, "Multipole"), (5, "Multipole")])
def test_multipole_kick_matches_pals_definition(pals_file, cosy_dir, order, kind):
    # x'' = -(K_N / N!) x^N on the midplane; the order-N coefficient of a
    # pure multipole is exactly -(K_N / N!) L. The tolerance covers the
    # ~4e-8 difference between COSY's proton preset mass and CODATA 2018.
    k_n, length = 3.0, 0.2
    _, m = _map(pals_file, cosy_dir,
                f"- m: {{kind: {kind}, length: {length}, MagneticMultipoleP: {{Kn{order}: {k_n}}}}}",
                order=order)
    assert m.coeff(1, order) == pytest.approx(-k_n / math.factorial(order) * length, rel=1e-6)


def test_skew_quadrupole_sign(pals_file, cosy_dir):
    # Bs1: By = -Bs1 y, Bx = Bs1 x, so a' = +Ks1 y and b' = +Ks1 x (thin lens).
    _, m = _map(pals_file, cosy_dir,
                "- q: {kind: Quadrupole, length: 0.01, MagneticMultipoleP: {Ks1: 2.0}}")
    assert m.coeff(1, 0, 0, 1) == pytest.approx(2.0 * 0.01, rel=1e-3)
    assert m.coeff(3, 1) == pytest.approx(2.0 * 0.01, rel=1e-3)


def test_tilted_quadrupole_equals_skew(pals_file, cosy_dir):
    _, tilted = _map(pals_file, cosy_dir,
                     "- q: {kind: Quadrupole, length: 0.3, MagneticMultipoleP: "
                     "{Kn1: 1.5, tilt1: 0.7853981633974483}}")
    _, skew = _map(pals_file, cosy_dir,
                   "- q: {kind: Quadrupole, length: 0.3, MagneticMultipoleP: {Ks1: -1.5}}")
    for row, exps in ((1, (0, 0, 1)), (3, (1,)), (0, (1,)), (2, (0, 1, 1))):
        assert tilted.coeff(row, *exps) == pytest.approx(skew.coeff(row, *exps),
                                                         rel=1e-9, abs=1e-12)


def test_thin_quadrupole(pals_file, cosy_dir):
    _, m = _map(pals_file, cosy_dir,
                "- q: {kind: Quadrupole, MagneticMultipoleP: {Kn1L: 0.1}}\n"
                "- d: {kind: Drift, length: 1.0}")
    assert m.coeff(1, 1) == pytest.approx(-0.1, rel=1e-6)
    assert m.coeff(0, 1) == pytest.approx(1 - 0.1, rel=1e-6)


@pytest.mark.parametrize("angle", [0.3, -0.3])
def test_rectangular_bend_has_no_horizontal_focusing(pals_file, cosy_dir, angle):
    _, m = _map(pals_file, cosy_dir,
                f"- b: {{kind: Bend, length: 1.0, BendP: {{angle_ref: {angle}, "
                f"e1_rect: 0.0, e2_rect: 0.0}}}}")
    rho = 1.0 / abs(angle)
    assert m.coeff(0, 1) == pytest.approx(1.0, abs=1e-12)
    assert m.coeff(0, 0, 1) == pytest.approx(rho * math.sin(abs(angle)), rel=1e-9)
    # vertical: entrance-edge kick -tan(e)/rho, then the 1 m body
    assert m.coeff(2, 0, 0, 1) == pytest.approx(1 - math.tan(abs(angle) / 2) / rho,
                                                rel=1e-9)


@pytest.mark.parametrize("angle", [0.3, -0.3])
def test_sector_bend(pals_file, cosy_dir, angle):
    _, m = _map(pals_file, cosy_dir,
                f"- b: {{kind: Bend, length: 1.0, BendP: {{g_ref: {angle}}}}}")
    assert m.coeff(0, 1) == pytest.approx(math.cos(angle), rel=1e-9)
    assert m.coeff(2, 0, 0, 1) == pytest.approx(1.0, abs=1e-12)


def test_bend_from_chord(pals_file, cosy_dir):
    _, m = _map(pals_file, cosy_dir,
                "- b: {kind: Bend, BendP: {g_ref: 0.25, L_chord: 1.9}}")
    angle = 2 * math.asin(0.25 * 1.9 / 2)
    assert m.scalars["SPOS"] == pytest.approx(angle / 0.25)
    assert m.coeff(0, 1) == pytest.approx(math.cos(angle), rel=1e-9)


def test_vertical_bend_from_tilt(pals_file, cosy_dir):
    # tilt_ref = +pi/2 rotates the bend direction -x into -y.
    _, m = _map(pals_file, cosy_dir,
                "- b: {kind: Bend, length: 1.0, BendP: {g_ref: 0.3, "
                "tilt_ref: 1.5707963267948966}}")
    dispersion = (1 - math.cos(0.3)) / 0.3
    assert m.coeff(2, 0, 0, 0, 0, 0, 1) > 0
    assert abs(m.coeff(0, 0, 0, 0, 0, 0, 1)) < 1e-12
    assert m.coeff(2, 0, 0, 1) == pytest.approx(math.cos(0.3), rel=1e-9)
    assert m.coeff(0, 1) == pytest.approx(1.0, abs=1e-12)
    assert m.coeff(2, 0, 0, 0, 0, 0, 1) == pytest.approx(
        dispersion * _dp_over_p_per_delta(), rel=1e-6)


def _dp_over_p_per_delta():
    """COSY's delta is the relative kinetic-energy deviation: dp/p = delta E K / (pc)^2."""
    m = 938.27208816
    pc = PC / 1e6
    e_tot = math.hypot(pc, m)
    return (e_tot - m) * e_tot / pc ** 2


@pytest.mark.parametrize("species", ["proton", "electron"])
def test_rf_on_crest_accelerates(pals_file, cosy_dir, species):
    bp, m = _map(pals_file, cosy_dir,
                 "- c: {kind: RFCavity, length: 0.6, RFP: {voltage: 2.0e6, "
                 "frequency: 1.3e9, dE_ref: 2.0e6}}", species=species)
    assert m.scalars["E0"] == pytest.approx(bp["kinetic_energy_mev"] + 2.0, rel=1e-12)


def test_rf_zero_crossing_below_transition(pals_file, cosy_dir):
    # Below transition the stable phase gives late particles (l < 0) more energy.
    bp, m = _map(pals_file, cosy_dir,
                 "- c: {kind: RFCavity, length: 0.6, RFP: {voltage: 2.0e6, "
                 "frequency: 1.3e9, zero_phase: BELOW_TRANSITION}}")
    assert m.scalars["E0"] == pytest.approx(bp["kinetic_energy_mev"], rel=1e-12)
    assert m.coeff(5, 0, 0, 0, 0, 1) < 0


def test_solenoid_rotation(pals_file, cosy_dir):
    _, m = _map(pals_file, cosy_dir,
                "- s: {kind: Solenoid, length: 1.0, SolenoidP: {Bsol: 2.0}}")
    kl = 2.0 / (2 * BRHO)
    # Hard-edge value; COSY's CMS has tanh fringes over the bore radius.
    assert m.coeff(0, 0, 0, 1) == pytest.approx(math.sin(kl) * math.cos(kl), rel=2e-3)
    assert m.coeff(0, 1) == pytest.approx(math.cos(kl) ** 2, rel=2e-3)
