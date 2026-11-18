"""PALS document features: expressions, versions, unsupported constructs."""

import math
import subprocess
import sys
import warnings

import pytest

from pals2cosy import PALSDeprecationWarning, UnsupportedFeatureError, __pals_commit__
from pals2cosy.constants import C_LIGHT
from pals2cosy.pals_parser import parse_lattice

BEGIN = """\
          - begin:
              kind: BeginningEle
              ReferenceP: {species_ref: electron, E_tot_ref: 1.0e9}
"""


def _doc(facility, root="", line="          - d: {kind: Drift, length: 1.0}\n"):
    return f"""\
PALS:
{root}  facility:
{facility}    - main:
        kind: BeamLine
        line:
{BEGIN}{line}"""


# --- Constants, variables and expressions in lattices ------------------------

def test_constants_variables_and_references(pals_file):
    text = _doc("""\
    - constants:
        - L0: 0.25
        - L1: 2 * L0
    - variables:
        - k: 0.5 / L1
    - g_length:
        kind: constant
        value: L0 + L1
    - q1:
        kind: Quadrupole
        length: L1
        MagneticMultipoleP: {Kn1: k}
""", line="""\
          - q1
          - d1: {kind: Drift, length: g_length + q1>length}
""")
    bp, elements = parse_lattice(pals_file(text))
    assert elements[0]["length"] == pytest.approx(0.5)
    assert elements[0]["bn"][1] == pytest.approx(1.0 * bp["brho_tm"])
    assert elements[1]["length"] == pytest.approx(1.25)


def test_redefinition_with_new_value_rejected(pals_file):
    text = _doc("    - variables:\n        - a: 1\n    - variables:\n        - a: 2\n")
    with pytest.raises(ValueError, match="redefined"):
        parse_lattice(pals_file(text))


def test_species_constant(pals_file):
    text = """\
PALS:
  facility:
    - constants:
        - ion: "#12C+6"
    - main:
        kind: BeamLine
        line:
          - begin:
              kind: BeginningEle
              ReferenceP: {species_ref: ion, pc_ref: 6.0e9}
          - d: {kind: Drift, length: 1.0}
"""
    bp, _ = parse_lattice(pals_file(text))
    assert bp["particle_type"] == "#12C+6"
    assert bp["brho_tm"] == pytest.approx(6.0e9 / 299792458.0 / 6)


# --- Version --------------------------------------------------------------------

@pytest.mark.parametrize("version", ["null", "1", '"1.0"', '"0.1"', '"1.0, 1.2"'])
def test_known_versions(pals_file, version):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        bp, _ = parse_lattice(pals_file(_doc("", root=f"  version: {version}\n")))
    assert bp["pals_version"] == (None if version == "null" else version.strip('"'))


def test_unknown_major_version_warns(pals_file):
    with pytest.warns(UserWarning, match=f"major version 2.*{__pals_commit__}"):
        parse_lattice(pals_file(_doc("", root='  version: "2.0.1"\n')))


def test_spec_commit_recorded():
    import pals2cosy
    assert pals2cosy.__pals_commit__ == "b9d5b355"
    out = subprocess.run([sys.executable, "-m", "pals2cosy", "--version"],
                         capture_output=True, text=True, check=True).stdout
    assert pals2cosy.__version__ in out and "b9d5b355" in out


# --- Unsupported and partially supported constructs ---------------------------

@pytest.mark.parametrize("facility, match", [
    ("    - set: {parameter: d>length, value: 2}\n", "set"),
    ("    - sets: [{d>length: 2}]\n", "sets"),
    ("    - superimpose: {place: d}\n", "superposition"),
    ("    - include: other.subpals.yaml\n", "include"),
    ("    - ps:\n        kind: Controller\n        controls: []\n", "Controller"),
])
def test_unsupported_facility_entries(pals_file, facility, match):
    with pytest.raises(UnsupportedFeatureError, match=match):
        parse_lattice(pals_file(_doc(facility)))


@pytest.mark.parametrize("root, match", [
    ("  load: [a.pals.yaml]\n", "load"),
    ("  include: a.subpals.yaml\n", "include"),
    ("  post_expansion:\n    - set: {parameter: d>length, value: 2}\n", "post_expansion"),
])
def test_unsupported_root_keys(pals_file, root, match):
    with pytest.raises(UnsupportedFeatureError, match=match):
        parse_lattice(pals_file(_doc("", root=root)))


@pytest.mark.parametrize("line, match", [
    ("          - d: {kind: Drift, length: 1.0, direction: -1}\n", "direction"),
    ("          - d: {kind: Drift, length: 1.0, placement: {offset: 1.0}}\n", "placement"),
    ("          - p: {kind: Patch, PatchP: {x_offset: 0.01}}\n", "Patch"),
    ("          - r: {kind: ReferenceChange, ReferenceChangeP: {dE_ref: 1.0e6}}\n",
     "ReferenceChange"),
])
def test_unsupported_line_constructs(pals_file, line, match):
    with pytest.raises(UnsupportedFeatureError, match=match):
        parse_lattice(pals_file(_doc("", line=line)))


def test_unsupported_error_is_value_and_not_implemented_error():
    assert issubclass(UnsupportedFeatureError, ValueError)
    assert issubclass(UnsupportedFeatureError, NotImplementedError)


@pytest.mark.parametrize("facility, root, line, match", [
    ("    - knob:\n        kind: Controller\n        control_type: RELATIVE\n",
     "", None, "RELATIVE Controller"),
    ("", "  adjustments:\n    - optimize: []\n", None, "adjustments"),
    ("", "  girders:\n    g1: {kind: Girder}\n", None, "girders"),
    ("", "", "          - f: {kind: Fork, ForkP: {to_line: main}}\n", "not followed"),
    ("", "", "          - k: {kind: Kicker, length: 0.1, MagneticMultipoleP: {Kn0: 0.001}}\n",
     "kick is not represented"),
    ("", "", "          - q: {kind: Quadrupole, length: 0.1, MagneticMultipoleP: {Bn1: 1},"
             " BodyShiftP: {x_offset: 0.001}}\n", "misalignment"),
    ("", "", "          - t: {kind: Taylor}\n", "no COSYScript mapping"),
    ("", "", "          - q: {kind: Quadrupole, length: 0.1, MagneticMultiPoleP: {Bn1: 1}}\n",
     "unknown key"),
])
def test_warned_constructs(pals_file, facility, root, line, match):
    kwargs = {"root": root} if line is None else {"root": root, "line": line}
    with pytest.warns(UserWarning, match=match):
        parse_lattice(pals_file(_doc(facility, **kwargs)))


def test_multipass_and_extra_branches_warn(pals_file):
    text = """\
PALS:
  facility:
    - begin:
        kind: BeginningEle
        ReferenceP: {species_ref: proton, E_tot_ref: 2.0e9}
    - linac: {kind: BeamLine, multipass: true, line: [d]}
    - d: {kind: Drift, length: 1.0}
    - a: {kind: BeamLine, line: [begin, linac, linac]}
    - b: {kind: BeamLine, line: [begin, d]}
    - lat: {kind: Lattice, branches: [a, b]}
"""
    with pytest.warns(UserWarning) as record:
        bp, elements = parse_lattice(pals_file(text))
    messages = " ".join(str(w.message) for w in record)
    assert "multipass" in messages and "2 branches" in messages
    assert bp["branch"] == "a" and len(elements) == 2
    bp, elements = parse_lattice(pals_file(text), beamline_name="b")
    assert bp["branch"] == "b" and len(elements) == 1


def test_zero_kicker_and_markers_silent(pals_file):
    line = ("          - k: {kind: Kicker, length: 0.1, MagneticMultipoleP: {Kn0: 0.0}}\n"
            "          - m: {kind: Marker}\n"
            "          - f: {kind: Fiducial}\n")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _, elements = parse_lattice(pals_file(_doc("", line=line)))
    assert [e["type"] for e in elements] == ["KICKER", "MARKER", "MARKER"]


def test_marker_with_length_rejected(pals_file):
    with pytest.raises(ValueError, match="zero length"):
        parse_lattice(pals_file(_doc("", line="          - m: {kind: Marker, length: 0.1}\n")))


def test_extension_labels_silence_unknown_keys(pals_file):
    root = "  extension_labels:\n    names: {Rotator: new kind}\n    prefixes: {X_: tool data}\n"
    line = ("          - q: {kind: Quadrupole, length: 0.1, MagneticMultipoleP: {Bn1: 1},"
            " X_data: {a: 1}}\n")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        parse_lattice(pals_file(_doc("", root=root, line=line)))
    with pytest.warns(UserWarning, match="extension kind 'Rotator'"):
        parse_lattice(pals_file(_doc("", root=root,
                                     line="          - r: {kind: Rotator, length: 0.5}\n")))


def test_reminders_are_reported(pals_file):
    from pals2cosy import PALSReminder
    with pytest.warns(PALSReminder, match="phase the RF"):
        parse_lattice(pals_file(_doc("", root="  reminders: [phase the RF]\n")))


def test_lattice_in_place_branch_and_inherit(pals_file):
    text = """\
PALS:
  facility:
    - begin:
        kind: BeginningEle
        ReferenceP: {species_ref: electron, E_tot_ref: 1.0e9}
    - q0: {kind: Quadrupole, length: 0.2, MagneticMultipoleP: {Bn1: 1.0}}
    - q1: {inherit: q0, MagneticMultipoleP: {Bn1: 2.0}}
    - lat:
        kind: Lattice
        branches:
          - ring:
              kind: BeamLine
              line: [begin, q0, q1, {q0: {length: 0.4}}]
    - use: lat
"""
    bp, elements = parse_lattice(pals_file(text))
    assert bp["branch"] == "ring"
    assert [(e["length"], e["bn"][1]) for e in elements] == [(0.2, 1.0), (0.2, 2.0), (0.4, 1.0)]


def test_parameter_group_inherit(pals_file):
    text = _doc("""\
    - ap1: {kind: ApertureP, y_width: 0.03}
    - b0: {kind: Bend, length: 1.0, BendP: {g_ref: 0.1}, ApertureP: {inherit: ap1}}
""", line="          - b0\n          - b1: {inherit: b0, ApertureP: {inherit: b0>ApertureP}}\n")
    _, elements = parse_lattice(pals_file(text))
    assert [e["pole_gap"] for e in elements] == [0.03, 0.03]


# --- Deprecated forms -----------------------------------------------------------

@pytest.mark.parametrize("line", [
    "          - b: {kind: SBend, length: 1.0, BendP: {g_ref: 0.1}}\n",
    "          - b: {kind: RBend, length: 1.0, BendP: {g_ref: 0.1}}\n",
    "          - s: {kind: Solenoid, length: 1.0, SolenoidP: {Bz: 0.5}}\n",
    "          - c: {kind: RFCavity, length: 1.0, RFCavityP: {voltage_kv: 10,"
    " frequency_hz: 1.3e9, phase_deg: 90}}\n",
])
def test_deprecated_forms_still_convert(pals_file, line):
    with pytest.warns(PALSDeprecationWarning):
        _, elements = parse_lattice(pals_file(_doc("", line=line)))
    assert elements


def test_legacy_lattice_particle(pals_file):
    text = """\
PALS:
  facility:
    - d: {kind: Drift, length: 1.0}
    - main: {kind: BeamLine, line: [d]}
    - lat:
        kind: Lattice
        particle: {species: proton, kinetic_energy: 2.5e8}
        branches: [main]
"""
    with pytest.warns(PALSDeprecationWarning, match="particle"):
        bp, _ = parse_lattice(pals_file(text))
    assert (bp["particle_type"], bp["kinetic_energy_mev"]) == ("proton", 250.0)


# --- Command line ---------------------------------------------------------------

def _cli(*args):
    return subprocess.run([sys.executable, "-m", "pals2cosy", *args],
                          capture_output=True, text=True)


def test_cli_strict_turns_warnings_into_errors(pals_file):
    path = pals_file(_doc("", line="          - q: {kind: Quadrupole, length: 0.1}\n"))
    assert _cli(path).returncode == 0
    result = _cli(path, "--strict")
    assert result.returncode == 1
    assert "no readable strength" in result.stderr


def test_cli_reports_unsupported_without_traceback(pals_file):
    path = pals_file(_doc("    - set: {parameter: d>length, value: 2}\n"))
    result = _cli(path)
    assert result.returncode == 1
    assert result.stderr.startswith("pals2cosy: error:")
    assert "Traceback" not in result.stderr


# --- Parameter-group keys, extensions and include -----------------------------

@pytest.mark.parametrize("line, match", [
    ("          - b: {kind: Bend, length: 1.0, BendP: {angle: 0.3}}\n",
     "unknown BendP key 'angle'"),
    ("          - c: {kind: RFCavity, length: 1.0, RFP: {volts: 1.0e6, frequency: 1.0e9}}\n",
     "unknown RFP key 'volts'"),
    ("          - s: {kind: Solenoid, length: 1.0, SolenoidP: {Bsol: 1, Kso: 2}}\n",
     "unknown SolenoidP key 'Kso'"),
    ("          - q: {kind: Quadrupole, length: 0.1, MagneticMultipoleP: {Bn1: 1, Bnn2: 3}}\n",
     "unknown MagneticMultipoleP key 'Bnn2'"),
    ("          - q: {kind: Quadrupole, length: 0.1, MagneticMultipoleP: {Bn1: 1},"
     " ApertureP: {y_widht: 0.02}}\n", "unknown ApertureP key 'y_widht'"),
])
def test_unknown_group_keys_warn(pals_file, line, match):
    with pytest.warns(UserWarning, match=match):
        parse_lattice(pals_file(_doc("", line=line)))


def test_misspelled_reference_key_warns(pals_file):
    text = _doc("", line="").replace("E_tot_ref: 1.0e9", "E_tot_ref: 1.0e9, E_tot: 5")
    with pytest.warns(UserWarning, match="unknown ReferenceP key 'E_tot'"):
        parse_lattice(pals_file(text + "          - d: {kind: Drift, length: 1.0}\n"))


def test_extension_blocks_at_any_level_are_silent(pals_file):
    root = "  my_tool: {extension: FooTool, data: 1}\n"
    line = ("          - q: {kind: Quadrupole, length: 0.1, my_ext: {extension: Foo},\n"
            "              MagneticMultipoleP: {Bn1: 1, fit: {extension: Fit, weight: 2}}}\n"
            "          - aside: {extension: Bar, stuff: 3}\n")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _, elements = parse_lattice(pals_file(_doc("", root=root, line=line)))
    assert [e["name"] for e in elements] == ["q"]


@pytest.mark.parametrize("line", [
    "          - q: {kind: Quadrupole, length: 0.5, include: q.subpals.yaml}\n",
    "          - q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {include: q.subpals.yaml}}\n",
    "          - include: more.subpals.yaml\n",
])
def test_include_at_any_level_is_unsupported(pals_file, line):
    with pytest.raises(UnsupportedFeatureError, match="'include' is not supported"):
        parse_lattice(pals_file(_doc("", line=line)))


@pytest.mark.parametrize("kind", ["Instrument", "Drift", "Wiggler", "Marker"])
def test_unconverted_multipoles_warn(pals_file, kind):
    length = "" if kind == "Marker" else "length: 0.5, "
    strength = "Kn1L: 0.1" if kind == "Marker" else "Kn1: 1.0"
    line = f"          - e: {{kind: {kind}, {length}MagneticMultipoleP: {{{strength}}}}}\n"
    with pytest.warns(UserWarning, match="MagneticMultipoleP fields are not represented"):
        parse_lattice(pals_file(_doc("", line=line)))


# --- BeamLine length ------------------------------------------------------------

def _line_doc(main, sub="    - sub: {kind: BeamLine, line: [d, q]}\n"):
    return f"""\
PALS:
  facility:
    - begin:
        kind: BeginningEle
        ReferenceP: {{species_ref: electron, E_tot_ref: 1.0e9}}
    - d: {{kind: Drift, length: 1.0}}
    - q: {{kind: Quadrupole, length: 0.5, MagneticMultipoleP: {{Bn1: 1.0}}}}
{sub}    - main:
        kind: BeamLine
{main}"""


def _layout(elements):
    return [(e["name"], e["length"], e["s_start"]) for e in elements]


def test_beamline_length_leaves_a_gap_at_its_end(pals_file):
    bp, elements = parse_lattice(pals_file(_line_doc("""\
        length: 6.0
        line: [begin, sub, {sub: {length: 2.0}}, d]
""", sub="    - sub: {kind: BeamLine, length: 1.8, line: [d, q]}\n")))
    assert _layout(elements) == [
        ("d", 1.0, 0.0), ("q", 0.5, 1.0), ("sub_end", pytest.approx(0.3), 1.5),
        ("d", 1.0, 1.8), ("q", 0.5, 2.8), ("sub_end", 0.5, 3.3),
        ("d", 1.0, 3.8), ("main_end", pytest.approx(1.2), 4.8)]
    assert elements[-1]["s_end"] == pytest.approx(6.0)
    assert bp["branch_time_s"] == pytest.approx(6.0 / C_LIGHT, rel=1e-6)


def test_reversed_beamline_keeps_its_gap_at_the_front(pals_file):
    _, elements = parse_lattice(pals_file(_line_doc(
        "        line: [begin, {sub: {repeat: -1}}]\n",
        sub="    - sub: {kind: BeamLine, length: 2.0, line: [d, q]}\n")))
    assert _layout(elements) == [("sub_end", 0.5, 0.0), ("q", 0.5, 0.5), ("d", 1.0, 1.0)]


def test_beamline_length_equal_to_contents_adds_nothing(pals_file):
    _, elements = parse_lattice(pals_file(_line_doc(
        "        length: 1.5\n        line: [begin, d, q]\n")))
    assert [e["name"] for e in elements] == ["d", "q"]


def test_beamline_shorter_than_contents_is_unsupported(pals_file):
    with pytest.raises(UnsupportedFeatureError, match="shorter than its line items"):
        parse_lattice(pals_file(_line_doc("        length: 1.0\n        line: [begin, d, q]\n")))


# --- Reference particle, time and energy ----------------------------------------

def test_beginning_element_must_come_first(pals_file):
    text = _line_doc("        line: [d, sub2]\n",
                     sub="    - sub2:\n        kind: BeamLine\n        line: [begin, q]\n")
    with pytest.warns(UserWarning, match="does not start with a BeginningEle"):
        bp, elements = parse_lattice(pals_file(text))
    assert bp["particle_type"] == "electron"
    assert [e["name"] for e in elements] == ["d", "q"]


@pytest.mark.parametrize("change, period", [
    ("dtime_ref: 1.0e-6", 1.0e-6 + 2.0 / C_LIGHT),
    ("time_ref: 4.0e-6", 4.0e-6 + 1.0 / C_LIGHT),
])
def test_reference_time_changes_enter_harmonic_frequency(pals_file, change, period):
    line = ("          - d: {kind: Drift, length: 1.0}\n"
            f"          - r: {{kind: ReferenceChange, ReferenceChangeP: {{{change}}}}}\n"
            "          - c: {kind: RFCavity, length: 1.0, RFP: {voltage: 1.0e3,"
            " harmonic_number: 2, zero_phase: BELOW_TRANSITION}}\n")
    bp, elements = parse_lattice(pals_file(_doc("", line=line)))
    assert bp["branch_time_s"] == pytest.approx(period, rel=1e-6)
    assert elements[-1]["rf_frequency_hz"] == pytest.approx(2 / period, rel=1e-6)


def test_reference_change_rejects_two_time_settings(pals_file):
    line = ("          - r: {kind: ReferenceChange, ReferenceChangeP:"
            " {dtime_ref: 1.0e-9, time_ref: 0.0}}\n")
    with pytest.raises(ValueError, match="only one of dtime_ref and time_ref"):
        parse_lattice(pals_file(_doc("", line=line)))


def test_floor_shift_with_length_rejected(pals_file):
    with pytest.raises(ValueError, match="zero length"):
        parse_lattice(pals_file(_doc("", line="          - f: {kind: FloorShift, length: 2.0}\n")))


@pytest.mark.parametrize("rfp", ["voltage: 0, frequency: 1.0e9, dE_ref: 5.0e8",
                                 "frequency: 1.0e9, dE_ref: 5.0e8"])
def test_rf_reference_energy_change_without_voltage(pals_file, rfp):
    # dE_ref always moves the PALS reference; COSY, which gets no RF kick, does
    # not follow, and the mismatch is reported.
    line = (f"          - c: {{kind: RFCavity, length: 1.0, RFP: {{{rfp}}}}}\n"
            "          - q: {kind: Quadrupole, length: 0.5, MagneticMultipoleP: {Kn1: 1}}\n")
    with pytest.warns(UserWarning, match="re-references"):
        _, elements = parse_lattice(pals_file(_doc("", line=line)))
    pc = math.sqrt(1.5e9 ** 2 - (0.51099895e6) ** 2)
    assert elements[1]["bn"][1] == pytest.approx(-pc / C_LIGHT, rel=1e-12)


# --- Expressions ----------------------------------------------------------------

@pytest.mark.parametrize("radius", ["Inf", "-Inf"])
def test_inf_radius_is_a_straight_bend(pals_file, radius):
    _, elements = parse_lattice(pals_file(_doc("", line=(
        f"          - b: {{kind: Bend, length: 1.0, BendP: {{radius_ref: {radius}}}}}\n"))))
    assert elements[0]["angle"] == 0.0


def test_random_values_are_shared_by_instances(pals_file):
    line = ("          - q: {kind: Quadrupole, length: 0.5, repeat: 3,"
            " MagneticMultipoleP: {Bn1: '1 + 0.1*random()'}}\n"
            "          - q2: {kind: Quadrupole, length: 0.5,"
            " MagneticMultipoleP: {Bn1: '1 + 0.1*random()'}}\n")
    with pytest.warns(UserWarning, match="fixed-seed"):
        _, elements = parse_lattice(pals_file(_doc("", line=line)))
    draws = [e["bn"][1] for e in elements]
    assert draws[0] == draws[1] == draws[2] != draws[3]


def test_redefinition_with_equal_value_accepted(pals_file):
    text = _doc("    - constants:\n        - a: 2*pi\n    - constants:\n        - a: twopi\n",
                line="          - d: {kind: Drift, length: a}\n")
    _, elements = parse_lattice(pals_file(text))
    assert elements[0]["length"] == pytest.approx(2 * math.pi)


@pytest.mark.parametrize("length", ["'(-8)^(1/3)'", "'10^400'", "'0^-1'"])
def test_failed_powers_are_reported_in_one_line(pals_file, length):
    path = pals_file(_doc("", line=f"          - d: {{kind: Drift, length: {length}}}\n"))
    with pytest.raises(ValueError, match="field 'length' must be numeric"):
        parse_lattice(path)
    result = _cli(path)
    assert result.returncode == 1 and "Traceback" not in result.stderr


def test_cli_twiss_defaults_to_two_dimensions(pals_file, tmp_path):
    out = tmp_path / "t.fox"
    path = pals_file(_doc("", line="          - d: {kind: Drift, length: 1.0}\n"))
    assert _cli(path, "--twiss", "-o", str(out)).returncode == 0
    assert "OV 3 2 0 ;" in out.read_text()
    result = _cli(path, "--twiss", "--dim", "3")
    assert result.returncode == 1 and "needs RF" in result.stderr


def test_cli_solenoid_fringe(pals_file, tmp_path):
    out = tmp_path / "s.fox"
    path = pals_file(_doc("", line="          - s: {kind: Solenoid, length: 0.5,"
                                   " SolenoidP: {Bsol: 1.0}}\n"))
    assert _cli(path, "--strict").returncode == 1
    assert _cli(path, "--strict", "--solenoid-fringe", "1e-4", "-o", str(out)).returncode == 0
    assert "CMS 1.0 0.0001 0.5 ;" in out.read_text()
