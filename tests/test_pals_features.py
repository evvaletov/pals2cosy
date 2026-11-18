"""PALS document features: expressions, versions, unsupported constructs."""

import subprocess
import sys
import warnings

import pytest

from pals2cosy import PALSDeprecationWarning, UnsupportedFeatureError, __pals_commit__
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
