"""Parse flat beamline JSON/YAML lattice files into normalized element lists.

Handles PALS CamelCase types (SBend, Quadrupole, etc.) and additional
type aliases. Supports MagneticMultipoleP, BendP, SolenoidP, RFCavityP,
and WigglerP parameter blocks.

Author: Eremey Valetov
"""

import math
import os
import warnings

from ._io import load_file as _load_file
from .constants import E0_ELECTRON, G_QUAD_DEFAULT, F_RF_DEFAULT

# Type resolution table: raw type -> internal short name.
# None means the type needs further resolution (polarity, plane, etc.)
_TYPE_ALIASES = {
    # Uppercase aliases
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
    "Sextupole": "SXT", "Octupole": "OCT", "Multipole": "MULT",
    "Kicker": None, "Instrument": None, "Marker": "DRIFT",
}

# Uppercase alias types (rejected in strict PALS mode)
_FELSIM_ONLY_TYPES = frozenset({
    "QUADRUPOLE", "DIPOLE", "DIPOLE_WEDGE", "DPW", "DPH", "QPF", "QPD",
    "UNDULATOR", "UND", "BPM", "OTR", "SPC", "XRS", "BSW",
    "CORRECTOR_V", "STV", "CORRECTOR_H", "STH", "SPECTROMETER",
})


def parse_lattice(path, mode="felsim"):
    """Parse a v2 JSON/YAML lattice file.

    Parameters
    ----------
    path : str
        Path to JSON or YAML lattice file.
    mode : str
        ``"felsim"`` accepts all type names (uppercase aliases + PALS CamelCase).
        ``"pals"`` rejects uppercase alias types (DIPOLE_WEDGE, etc.).

    Returns
    -------
    tuple
        (beam_params, elements) where beam_params is a dict with keys like
        kinetic_energy_mev, particle_type, rf_frequency_hz, quad_gradient,
        and elements is a list of normalized dicts with drifts inserted.
    """
    if mode not in ("pals", "felsim"):
        raise ValueError(f"Unknown mode '{mode}', expected 'pals' or 'felsim'")

    data = _load_file(path)
    beamline = data["beamline"]

    meta = beamline["metadata"]
    fv = meta.get("format_version", 2)
    if fv not in (1, 2, 3):
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
        elem = _normalize_element(raw, mode=mode)
        if elem is not None:
            positioned.append(elem)

    positioned.sort(key=lambda e: e["s_start"])

    # Insert drifts between elements
    result = []
    prev_end = 0.0
    for elem in positioned:
        gap = elem["s_start"] - prev_end
        if gap < -1e-9:
            warnings.warn(
                f"Element '{elem.get('name', '<unnamed>')}' overlaps previous element "
                f"by {-gap:.6g} m (s_start={elem['s_start']}, prev_end={prev_end})"
            )
        elif gap > 1e-9:
            result.append(_make_drift(gap))
        result.append(elem)
        prev_end = elem["s_end"]

    return beam_params, result


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


def _normalize_element(raw, mode="felsim"):
    """Convert a raw element dict into a normalized internal dict."""
    raw_type = raw.get("type") or raw.get("kind")
    if raw_type is None:
        warnings.warn(
            f"Element '{raw.get('name', '<unnamed>')}' has no 'type' or 'kind'; skipping"
        )
        return None

    if mode == "pals" and raw_type in _FELSIM_ONLY_TYPES:
        name = raw.get("name", "<unnamed>")
        raise ValueError(
            f"Element '{name}' uses type '{raw_type}', "
            f"which is not allowed in strict PALS mode. "
            f"Use PALS CamelCase types (SBend, Quadrupole, etc.) instead."
        )

    internal_type = _resolve_type(raw_type, raw)
    s_start = raw.get("s_start_m")
    s_end = raw.get("s_end_m")
    if s_start is None or s_end is None:
        warnings.warn(
            f"Element '{raw.get('name', '<unnamed>')}' missing s_start_m or s_end_m; skipping"
        )
        return None
    if s_end < s_start:
        warnings.warn(
            f"Element '{raw.get('name', '<unnamed>')}' has s_end ({s_end}) < s_start ({s_start}); skipping"
        )
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
        "bn1": None,
        "bn2": None,
        "bn3": None,
        "bn4": None,
        "bn5": None,
        "bz": None,
        "rf_voltage_kv": None,
        "rf_frequency_hz": None,
        "rf_phase_deg": None,
        "wiggler_field": None,
        "wiggler_period": None,
    }

    # MagneticMultipoleP for quadrupoles, sextupoles, octupoles
    mmp = raw.get("MagneticMultipoleP")
    if mmp is not None:
        bn1 = mmp.get("Bn1")
        if bn1 is not None:
            elem["bn1"] = float(bn1)
        bn2 = mmp.get("Bn2")
        if bn2 is not None:
            elem["bn2"] = float(bn2)
        bn3 = mmp.get("Bn3")
        if bn3 is not None:
            elem["bn3"] = float(bn3)
        bn4 = mmp.get("Bn4")
        if bn4 is not None:
            elem["bn4"] = float(bn4)
        bn5 = mmp.get("Bn5")
        if bn5 is not None:
            elem["bn5"] = float(bn5)

    # Wiggler/undulator parameters
    wig_p = raw.get("WigglerP")
    if wig_p is not None:
        bw = wig_p.get("peak_field")
        if bw is not None:
            elem["wiggler_field"] = float(bw)
        period = wig_p.get("period")
        if period is not None:
            elem["wiggler_period"] = float(period)

    # Solenoid on-axis field
    sol_p = raw.get("SolenoidP")
    if sol_p is not None:
        bz = sol_p.get("Bz")
        if bz is not None:
            elem["bz"] = float(bz)

    # RF cavity parameters
    rf_p = raw.get("RFCavityP")
    if rf_p is not None:
        v = rf_p.get("voltage_kv")
        if v is not None:
            elem["rf_voltage_kv"] = float(v)
        f = rf_p.get("frequency_hz")
        if f is not None:
            elem["rf_frequency_hz"] = float(f)
        phi = rf_p.get("phase_deg")
        if phi is not None:
            elem["rf_phase_deg"] = float(phi)

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

        # v3: BendP overrides
        bend_p = raw.get("BendP")
        if bend_p is not None:
            g_ref = bend_p.get("g_ref")
            if g_ref is not None:
                elem["angle"] = math.degrees(g_ref * length)
            e1 = bend_p.get("e1")
            if e1 is not None:
                elem["entrance_edge_angle"] = math.degrees(e1)
            e2 = bend_p.get("e2")
            if e2 is not None:
                elem["exit_edge_angle"] = math.degrees(e2)

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
        "bn1": None,
        "bn2": None,
        "bn3": None,
        "bn4": None,
        "bn5": None,
        "bz": None,
        "rf_voltage_kv": None,
        "rf_frequency_hz": None,
        "rf_phase_deg": None,
        "wiggler_field": None,
        "wiggler_period": None,
    }
