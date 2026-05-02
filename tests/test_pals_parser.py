"""Tests for official PALS format parser (pals_parser.py)."""

import os
import warnings
import pytest

from pals2cosy.pals_parser import parse_lattice, _build_catalog, _expand_beamline

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
FODO_YAML = os.path.join(EXAMPLES_DIR, "fodo.pals.yaml")


@pytest.fixture
def fodo():
    return parse_lattice(FODO_YAML, ke_override=1000, particle_override="proton")


def test_fodo_parse(fodo):
    """Parse fodo.pals.yaml and verify element sequence."""
    bp, elements = fodo
    assert bp["kinetic_energy_mev"] == 1000
    assert bp["particle_type"] == "proton"
    assert len(elements) > 0


def test_fodo_element_sequence(fodo):
    """Verify the expanded FODO cell sequence: drift-quad-drift-quad-drift × 3."""
    _, elements = fodo
    kinds = [e["type"] for e in elements]
    # 3 repetitions of fodo_cell: drift1, quad1, drift2, quad2, drift1
    # = 15 elements total
    assert len(elements) == 15
    # Each cell: DRIFT, QPF, DRIFT, QPD, DRIFT
    cell = kinds[:5]
    assert cell == ["DRIFT", "QPF", "DRIFT", "QPD", "DRIFT"]
    assert kinds[5:10] == cell
    assert kinds[10:15] == cell


def test_beamline_expansion():
    """Verify line: reference resolution from the FODO example."""
    import yaml
    with open(FODO_YAML) as f:
        data = yaml.safe_load(f)
    catalog = _build_catalog(data["PALS"]["facility"])
    flat = _expand_beamline("fodo_cell", catalog)
    names = [e.get("name", "") for e in flat]
    assert names == ["drift1", "quad1", "drift2", "quad2", "drift1"]


def test_repeat():
    """Verify repeat=N produces N copies."""
    import yaml
    with open(FODO_YAML) as f:
        data = yaml.safe_load(f)
    catalog = _build_catalog(data["PALS"]["facility"])
    flat = _expand_beamline("fodo_channel", catalog)
    # fodo_channel repeats fodo_cell 3 times → 15 elements
    assert len(flat) == 15


def test_inherit():
    """Verify element override via inherit (quad2 inherits quad1, changes Bn1)."""
    _, elements = parse_lattice(FODO_YAML, ke_override=1000, particle_override="proton")
    quads = [e for e in elements if e["type"] in ("QPF", "QPD")]
    # quad1: Bn1=1.0, quad2: Bn1=-1.0 (inherits quad1, overrides Bn1)
    q1 = quads[0]
    q2 = quads[1]
    assert q1["bn1"] == pytest.approx(1.0)
    assert q2["bn1"] == pytest.approx(-1.0)


def test_quadrupole_bn1(fodo):
    """Verify Bn1 extraction into normalized dict."""
    _, elements = fodo
    quads = [e for e in elements if e["type"] in ("QPF", "QPD")]
    assert len(quads) == 6  # 2 quads per cell × 3 cells
    # All QPF quads have Bn1=1.0
    qpf = [q for q in quads if q["type"] == "QPF"]
    assert all(q["bn1"] == pytest.approx(1.0) for q in qpf)
    # All QPD quads have Bn1=-1.0
    qpd = [q for q in quads if q["type"] == "QPD"]
    assert all(q["bn1"] == pytest.approx(-1.0) for q in qpd)


def test_cumulative_positions(fodo):
    """Verify s_start/s_end are computed from cumulative lengths."""
    _, elements = fodo
    # drift1(0.25) + quad1(1.0) + drift2(0.5) + quad2(1.0) + drift1(0.25) = 3.0 per cell
    cell_length = 3.0
    assert elements[0]["s_start"] == pytest.approx(0.0)
    assert elements[0]["s_end"] == pytest.approx(0.25)
    # Last element of first cell
    assert elements[4]["s_end"] == pytest.approx(cell_length)
    # First element of second cell
    assert elements[5]["s_start"] == pytest.approx(cell_length)
    # Total beamline
    assert elements[-1]["s_end"] == pytest.approx(3 * cell_length)


def test_element_lengths(fodo):
    """Verify element lengths from the FODO example."""
    _, elements = fodo
    assert elements[0]["length"] == pytest.approx(0.25)  # drift1
    assert elements[1]["length"] == pytest.approx(1.0)    # quad1
    assert elements[2]["length"] == pytest.approx(0.5)    # drift2
    assert elements[3]["length"] == pytest.approx(1.0)    # quad2
    assert elements[4]["length"] == pytest.approx(0.25)   # drift1


def test_sbend_angle_conversion():
    """Verify g_ref × l → degrees conversion for SBend elements."""
    import math
    import yaml
    # Create a minimal PALS file with an SBend
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"b1": {
                    "kind": "SBend",
                    "length": 2.0,
                    "BendP": {"g_ref": 0.5, "e1": 0.1, "e2": 0.2},
                }},
                {"my_line": {
                    "kind": "BeamLine",
                    "line": ["b1"],
                }},
            ]
        }
    }
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        _, elements = parse_lattice(tmp_path, ke_override=1000, particle_override="proton")
        bend = [e for e in elements if e["type"] == "DPH"][0]
        expected_angle = math.degrees(0.5 * 2.0)
        assert bend["angle"] == pytest.approx(expected_angle)
    finally:
        os.unlink(tmp_path)


def test_sbend_edge_angles():
    """Verify e1/e2 radian→degree conversion."""
    import math
    import yaml
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"b1": {
                    "kind": "SBend",
                    "length": 1.0,
                    "BendP": {"g_ref": 0.3, "e1": 0.15, "e2": 0.25},
                }},
                {"line1": {
                    "kind": "BeamLine",
                    "line": ["b1"],
                }},
            ]
        }
    }
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        _, elements = parse_lattice(tmp_path, ke_override=1000, particle_override="proton")
        bend = [e for e in elements if e["type"] == "DPH"][0]
        assert bend["entrance_edge_angle"] == pytest.approx(math.degrees(0.15))
        assert bend["exit_edge_angle"] == pytest.approx(math.degrees(0.25))
    finally:
        os.unlink(tmp_path)


def test_circular_reference_detection():
    """Circular BeamLine references should raise ValueError."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"line_a": {"kind": "BeamLine", "line": ["line_b"]}},
                {"line_b": {"kind": "BeamLine", "line": ["line_a"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with pytest.raises(ValueError, match="Circular reference"):
            parse_lattice(tmp_path, ke_override=1000, particle_override="proton",
                         beamline_name="line_a")
    finally:
        os.unlink(tmp_path)


def test_repeat_zero_produces_no_elements():
    """repeat: 0 on a BeamLine produces no elements and warns."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": 1.0}},
                {"cell": {"kind": "BeamLine", "line": ["d1"]}},
                {"main": {
                    "kind": "BeamLine",
                    "line": [{"cell": {"repeat": 0}}],
                }},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, elements = parse_lattice(
                tmp_path, ke_override=1000, particle_override="proton",
                beamline_name="main"
            )
        assert len(elements) == 0
        assert any("repeat=0" in str(x.message) for x in w)
    finally:
        os.unlink(tmp_path)


# --- Phase 2B: Edge case tests ---

def test_missing_kind_warns():
    """Element without 'kind' emits a warning and defaults to Drift."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"mystery": {"length": 1.0}},  # no 'kind' key
                {"line1": {"kind": "BeamLine", "line": ["mystery"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, elements = parse_lattice(tmp_path, ke_override=1000,
                                        particle_override="proton")
        assert any("no 'kind'" in str(warning.message) for warning in w)
        assert elements[0]["type"] == "DRIFT"
    finally:
        os.unlink(tmp_path)


def test_unresolved_reference_raises():
    """BeamLine referencing a nonexistent element raises ValueError."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"line1": {"kind": "BeamLine", "line": ["no_such_element"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with pytest.raises(ValueError, match="Unresolved reference"):
            parse_lattice(tmp_path, ke_override=1000, particle_override="proton",
                         beamline_name="line1")
    finally:
        os.unlink(tmp_path)


def test_octupole_emitted_as_drift():
    """Octupole element is emitted as a passive drift in COSYScript output."""
    import yaml, tempfile
    from pals2cosy.converter import convert
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"oct1": {"kind": "Octupole", "length": 0.3}},
                {"line1": {"kind": "BeamLine", "line": ["oct1"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        bp, elems = parse_lattice(tmp_path, ke_override=1000,
                                   particle_override="proton")
        fox = convert(bp, elems, ke_override=1000, particle_override="proton")
        assert "DL 0.3" in fox
        # Should not trigger unknown element warning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            convert(bp, elems, ke_override=1000, particle_override="proton")
        assert not any("Unknown element type" in str(x.message) for x in w)
    finally:
        os.unlink(tmp_path)


def test_solenoid_emitted_as_drift():
    """Solenoid element is emitted as a passive drift in COSYScript output."""
    import yaml, tempfile
    from pals2cosy.converter import convert
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"sol1": {"kind": "Solenoid", "length": 0.5}},
                {"line1": {"kind": "BeamLine", "line": ["sol1"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        bp, elems = parse_lattice(tmp_path, ke_override=1000,
                                   particle_override="proton")
        fox = convert(bp, elems, ke_override=1000, particle_override="proton")
        assert "DL 0.5" in fox
    finally:
        os.unlink(tmp_path)


# --- QA4: Warning path coverage ---

def test_negative_repeat_raises():
    """Negative repeat is unsupported and raises ValueError."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": 1.0}},
                {"cell": {"kind": "BeamLine", "line": ["d1"]}},
                {"main": {
                    "kind": "BeamLine",
                    "line": [{"cell": {"repeat": -3}}],
                }},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with pytest.raises(ValueError, match="negative repeat"):
            parse_lattice(tmp_path, ke_override=1000,
                          particle_override="proton",
                          beamline_name="main")
    finally:
        os.unlink(tmp_path)


def test_extra_keys_on_beamline_ref_warns():
    """Extra keys on a BeamLine reference emit a warning."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": 1.0}},
                {"cell": {"kind": "BeamLine", "line": ["d1"]}},
                {"main": {
                    "kind": "BeamLine",
                    "line": [{"cell": {"repeat": 2, "extra_key": "ignored"}}],
                }},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, elements = parse_lattice(tmp_path, ke_override=1000,
                                        particle_override="proton",
                                        beamline_name="main")
        assert any("ignored" in str(x.message) for x in w)
        assert len(elements) == 2  # repeat=2 still works
    finally:
        os.unlink(tmp_path)


def test_direction_modifier_raises():
    """Direction modifier is unsupported and raises ValueError."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": 1.0}},
                {"main": {
                    "kind": "BeamLine",
                    "line": [{"d1": {"direction": "reversed"}}],
                }},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with pytest.raises(ValueError, match="direction"):
            parse_lattice(tmp_path, ke_override=1000,
                          particle_override="proton",
                          beamline_name="main")
    finally:
        os.unlink(tmp_path)


def test_default_particle_warns():
    """Missing particle type in lattice and CLI emits a warning and defaults to proton."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": 1.0}},
                {"main": {"kind": "BeamLine", "line": ["d1"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bp, _ = parse_lattice(tmp_path, ke_override=1000)
        assert any("defaulting to 'proton'" in str(x.message) for x in w)
        assert bp["particle_type"] == "proton"
    finally:
        os.unlink(tmp_path)


def test_unknown_species_mass_warns():
    """Unknown particle species from lattice file emits a mass fallback warning."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": 1.0}},
                {"main": {"kind": "BeamLine", "line": ["d1"]}},
                {"lat": {
                    "kind": "Lattice",
                    "particle": {"species": "tachyon", "kinetic_energy": 500e6},
                    "branches": ["main"],
                }},
                {"use": "lat"},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bp, _ = parse_lattice(tmp_path)
        assert any("Unknown particle species" in str(x.message) for x in w)
        assert bp["particle_type"] == "tachyon"
    finally:
        os.unlink(tmp_path)


def test_extended_species_mass_recognized():
    """Extended particle species (e.g. positron, muon) get correct mass without warning."""
    import yaml, tempfile
    from pals2cosy.constants import E0_MUON
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": 1.0}},
                {"main": {"kind": "BeamLine", "line": ["d1"]}},
                {"lat": {
                    "kind": "Lattice",
                    "particle": {"species": "muon", "kinetic_energy": 500e6},
                    "branches": ["main"],
                }},
                {"use": "lat"},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bp, _ = parse_lattice(tmp_path)
        assert not any("Unknown particle species" in str(x.message) for x in w)
        assert bp["particle_type"] == "muon"
        assert abs(bp["mass_mev"] - E0_MUON) < 1e-6
    finally:
        os.unlink(tmp_path)


def test_nan_length_rejected():
    """A NaN length value is rejected with a clear error."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": float("nan")}},
                {"main": {"kind": "BeamLine", "line": ["d1"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with pytest.raises(ValueError, match="finite"):
            parse_lattice(tmp_path, ke_override=1000,
                          particle_override="proton",
                          beamline_name="main")
    finally:
        os.unlink(tmp_path)


def test_nonnumeric_bn1_rejected():
    """Quadrupole Bn1 with non-numeric value is rejected with a clear error."""
    import yaml, tempfile
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"q1": {"kind": "Quadrupole", "length": 0.1,
                        "MagneticMultipoleP": {"Bn1": "bad"}}},
                {"main": {"kind": "BeamLine", "line": ["q1"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with pytest.raises(ValueError, match="MagneticMultipoleP.Bn1"):
            parse_lattice(tmp_path, ke_override=1000,
                          particle_override="proton",
                          beamline_name="main")
    finally:
        os.unlink(tmp_path)


def test_particle_override_case_insensitive():
    """--particle Positron (mixed case) resolves to positron mass, not proton fallback."""
    import yaml, tempfile
    from pals2cosy.constants import E0_POSITRON
    pals_data = {
        "PALS": {
            "version": None,
            "facility": [
                {"d1": {"kind": "Drift", "length": 1.0}},
                {"main": {"kind": "BeamLine", "line": ["d1"]}},
            ]
        }
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pals.yaml', delete=False) as f:
        yaml.dump(pals_data, f)
        tmp_path = f.name
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bp, _ = parse_lattice(tmp_path, ke_override=1000,
                                  particle_override="Positron",
                                  beamline_name="main")
        assert not any("Unknown particle species" in str(x.message) for x in w)
        assert abs(bp["mass_mev"] - E0_POSITRON) < 1e-6
    finally:
        os.unlink(tmp_path)
