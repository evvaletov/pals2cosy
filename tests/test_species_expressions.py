"""PALS species names and the expression evaluator."""

import pytest

from pals2cosy.expressions import ExpressionError, Scope
from pals2cosy.species import charge_of, mass_of, parse_species


# --- Species ------------------------------------------------------------------

@pytest.mark.parametrize("name, mass_mev, charge", [
    ("electron", 0.51099895, -1),
    ("positron", 0.51099895, 1),
    ("proton", 938.27208816, 1),
    ("anti-proton", 938.27208816, -1),
    ("antiproton", 938.27208816, -1),
    ("muon", 105.6583755, -1),
    ("anti-muon", 105.6583755, 1),
    ("mu+", 105.6583755, 1),
    ("pion+", 139.57039, 1),
    ("pi-", 139.57039, -1),
    ("deuteron", 1875.61294257, 1),
    ("alpha", 3727.3794066, 2),
    ("#4He+2", 3727.3794066, 2),
    ("#3He+2", 2808.39160743, 2),
    ("#2H+1", 1875.61294257, 1),
    ("#12C+6", 12 * 931.49410242 - 6 * 0.51099895, 6),
    ("Au+79", 196.96656879 * 931.49410242 - 79 * 0.51099895, 79),
    ("#238U+92", 238.0507884 * 931.49410242 - 92 * 0.51099895, 92),
    ("Ar+", 39.9623831237 * 931.49410242 - 0.51099895, 1),
    ("H-", 1.00782503223 * 931.49410242 + 0.51099895, -1),
])
def test_species(name, mass_mev, charge):
    sp = parse_species(name)
    # Ion masses neglect electron binding energies (< 1e-5 relative).
    assert sp.mass_mev == pytest.approx(mass_mev, rel=1e-5)
    assert sp.charge == charge


def test_species_presets():
    assert parse_species("electron").cosy_command(40) == "RPE 40"
    assert parse_species("proton").cosy_command(40) == "RPP 40"
    assert parse_species("#1H+1").cosy_command(40) == "RPP 40"
    assert parse_species("pion-").cosy_command(40) == "RPPI 40 -1"
    assert parse_species("muon").cosy_command(40) == "RPMU 40"
    assert parse_species("anti-muon").cosy_command(40).startswith("RP 40 ")


@pytest.mark.parametrize("name", ["tachyon", "#3Xx", "#99He+2", "He+3", "au+79", ""])
def test_species_errors(name):
    with pytest.raises(ValueError):
        parse_species(name)


def test_species_functions():
    assert mass_of("proton") == pytest.approx(938.27208816e6)
    assert charge_of("#3He+2") == 2.0
    assert charge_of("#3He") == 0.0


# --- Expressions ----------------------------------------------------------------

@pytest.mark.parametrize("text, value", [
    ("1 + 2 * 3", 7.0),
    ("(1 + 2) * 3", 9.0),
    ("-2^2", -4.0),
    ("2^-1", 0.5),
    ("2**3**2", 512.0),
    ("1.0e9 / 1e3", 1e6),
    (".5 * 4", 2.0),
    ("sqrt(16) + abs(-1) + exp(0) + log(1)", 6.0),
    ("atan2(1, 1) * 4", 3.141592653589793),
    ("twopi / pi", 2.0),
    ("modulo(7, 3) + sign(-2) + nint(2.5) + int(-2.7)", 1 - 1 + 3 - 2),
    ("expr(3 * 2)", 6.0),
    ("factorial(4) / 24", 1.0),
    ("mass_of(\"electron\") / 1e6", 0.51099895),
])
def test_expression_values(text, value):
    assert Scope().evaluate(text) == pytest.approx(value)


@pytest.mark.parametrize("text", [
    "__import__('os').system('true')", "1 +", "unknown_name", "2 // 3",
    "sin(", "1 / 0", "[1, 2]",
])
def test_expression_errors(text):
    with pytest.raises(ExpressionError):
        Scope().evaluate(text)


def test_expression_scope_and_cycles():
    scope = Scope({"a": "2 * b", "b": 3, "s": "#3He", "m": "mass_of(s)"})
    assert scope.lookup("a") == 6.0
    assert scope.lookup("s") == "#3He"
    assert scope.lookup("m") == pytest.approx(mass_of("#3He"))
    with pytest.raises(ExpressionError, match="circular"):
        Scope({"x": "y", "y": "x + 1"}).lookup("x")
