"""Parse v2 PALS JSON/YAML lattice files into normalized element lists.

Self-contained parser with no FELsim dependency. Handles both PALS-native
SBend/RBend (with edge angle attributes) and FELsim DIPOLE_WEDGE elements.

Author: Eremey Valetov
"""

import json
import os

from .constants import E0_ELECTRON, G_QUAD_DEFAULT, F_RF_DEFAULT

# Type resolution table: raw type -> internal short name.
# None means the type needs further resolution (polarity, plane, etc.)
_TYPE_ALIASES = {
    # FELsim-native
    "QUADRUPOLE": None, "QPF": "QPF", "QPD": "QPD",
    "DIPOLE": "DPH", "DPH": "DPH",
    "DIPOLE_WEDGE": "DPW", "DPW": "DPW",
    "UNDULATOR": "UND", "UND": "UND",
    "BPM": "BPM", "OTR": "OTR",
    "CORRECTOR_V": "STV", "STV": "STV",
    "CORRECTOR_H": "STH", "STH": "STH",
    "SPECTROMETER": "SPC", "SPC": "SPC",
    "XRS": "XRS", "BSW": "BSW",
    "DRIFT": "DRIFT",
    # PALS CamelCase
    "Drift": "DRIFT", "Quadrupole": None,
    "SBend": "DPH", "RBend": "DPH",
    "Wiggler": "UND", "Solenoid": "SOL", "RFCavity": "RFC",
    "Sextupole": "SXT",
    "Kicker": None, "Instrument": None, "Marker": "DRIFT",
}

# Types treated as passive (zero-length or drift-like)
_PASSIVE_TYPES = {"BPM", "OTR", "STV", "STH", "SPC", "XRS", "BSW", "UND", "SOL", "RFC", "SXT"}


def parse_lattice(path):
    """Parse a v2 JSON/YAML lattice file.

    Returns:
        (beam_params, elements) where beam_params is a dict with keys like
        kinetic_energy_mev, particle_type, rf_frequency_hz, quad_gradient,
        and elements is a list of normalized dicts with drifts inserted.
    """
    data = _load_file(path)
    beamline = data["beamline"]

    meta = beamline["metadata"]
    fv = meta.get("format_version", 2)
    if fv not in (1, 2):
        raise ValueError(f"Unsupported format_version {fv}")

    bp = beamline["beam_parameters"]
    particle = bp["particle"]
    global_settings = beamline.get("global_settings", {})

    beam_params = {
        "kinetic_energy_mev": particle["kinetic_energy_mev"],
        "particle_type": particle["type"],
        "mass_mev": particle.get("mass_mev", E0_ELECTRON),
        "rf_frequency_hz": bp.get("rf_frequency_hz", F_RF_DEFAULT),
        "quad_gradient": global_settings.get(
            "quadrupole_gradient_coefficient_t_per_a_per_m", G_QUAD_DEFAULT
        ),
        "source_file": os.path.basename(path),
    }

    raw_elements = beamline["elements"]
    positioned = []
    for raw in raw_elements:
        elem = _normalize_element(raw)
        if elem is not None:
            positioned.append(elem)

    positioned.sort(key=lambda e: e["s_start"])

    # Insert drifts between elements
    result = []
    prev_end = 0.0
    for elem in positioned:
        gap = elem["s_start"] - prev_end
        if gap > 1e-9:
            result.append(_make_drift(gap))
        result.append(elem)
        prev_end = elem["s_end"]

    return beam_params, result


def _load_file(path):
    ext = os.path.splitext(path)[1].lower()
    with open(path, "r") as f:
        if ext in (".yaml", ".yml"):
            import yaml
            return yaml.safe_load(f)
        else:
            return json.load(f)


def _resolve_type(raw_type, raw_elem):
    """Resolve element type string to internal short name."""
    if raw_type in ("QUADRUPOLE", "Quadrupole"):
        polarity = raw_elem.get("polarity", "")
        return "QPF" if polarity == "focusing" else "QPD"
    if raw_type == "Kicker":
        plane = raw_elem.get("plane", "vertical")
        return "STV" if plane == "vertical" else "STH"
    if raw_type == "Instrument":
        return raw_elem.get("instrument_type", "BPM")

    short = _TYPE_ALIASES.get(raw_type, raw_type)
    return short if short is not None else raw_type


def _normalize_element(raw):
    """Convert a raw element dict into a normalized internal dict."""
    raw_type = raw.get("type") or raw.get("kind")
    if raw_type is None:
        return None

    internal_type = _resolve_type(raw_type, raw)
    s_start = raw.get("s_start_m")
    s_end = raw.get("s_end_m")
    if s_start is None or s_end is None:
        return None

    length = raw.get("length_m", s_end - s_start)
    params = raw.get("parameters", {})
    name = raw.get("name", "")

    elem = {
        "type": internal_type,
        "name": name,
        "length": length,
        "s_start": s_start,
        "s_end": s_end,
        "current": params.get("current_a"),
        "angle": None,
        "wedge_angle": None,
        "entrance_edge_angle": None,
        "exit_edge_angle": None,
        "pole_gap": None,
        "dipole_length": None,
        "enge_coeffs": None,
    }

    if internal_type == "DPH":
        elem["angle"] = params.get("bending_angle_deg", 0)
        elem["dipole_length"] = params.get("dipole_length_m", length)
        elem["pole_gap"] = params.get("pole_gap_m", 0)
        # PALS-native edge angles (on the SBend element itself)
        elem["entrance_edge_angle"] = params.get("entrance_edge_angle_deg")
        elem["exit_edge_angle"] = params.get("exit_edge_angle_deg")
        # Enge coefficients on the SBend element itself
        ff = raw.get("fringe_fields", {})
        coeffs = ff.get("enge_coefficients")
        if coeffs:
            elem["enge_coeffs"] = list(coeffs)

    elif internal_type == "DPW":
        elem["wedge_angle"] = params.get("wedge_angle_deg", 0)
        elem["angle"] = params.get("dipole_angle_deg", 0)
        elem["dipole_length"] = params.get("dipole_length_m", 0)
        elem["pole_gap"] = params.get("pole_gap_m", 0)
        ff = raw.get("fringe_fields", {})
        coeffs = ff.get("enge_coefficients")
        if coeffs:
            elem["enge_coeffs"] = list(coeffs)

    return elem


def _make_drift(length):
    return {
        "type": "DRIFT",
        "name": "",
        "length": length,
        "s_start": None,
        "s_end": None,
        "current": None,
        "angle": None,
        "wedge_angle": None,
        "entrance_edge_angle": None,
        "exit_edge_angle": None,
        "pole_gap": None,
        "dipole_length": None,
        "enge_coeffs": None,
    }
