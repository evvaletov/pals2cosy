"""Parse official PALS format (``PALS:`` root) into normalized element lists.

Handles ``facility:`` element definitions and ``BeamLine`` composition with
``line:`` references, ``inherit:`` overrides, and ``repeat:`` expansion.
Produces the same (beam_params, elements) tuple as ``parser.parse_lattice()``.

Author: Eremey Valetov
"""

import math
import os
import warnings

from ._io import load_file as _load_file
from .constants import E0_ELECTRON, E0_PROTON, F_RF_DEFAULT, PARTICLE_MASSES_MEV


def parse_lattice(path, beamline_name=None, ke_override=None, particle_override=None):
    """Parse an official PALS file (``PALS:`` root key).

    Parameters
    ----------
    path : str
        Path to ``.pals.yaml`` or ``.pals.json`` file.
    beamline_name : str or None
        Which BeamLine to expand. Default: last ``BeamLine`` in facility.
    ke_override : float or None
        Kinetic energy in MeV (required if not in lattice file).
    particle_override : str or None
        Particle type override (``"electron"`` or ``"proton"``).

    Returns
    -------
    tuple
        (beam_params, elements) matching the ``parser.parse_lattice()`` contract.
    """
    data = _load_file(path)
    if not isinstance(data, dict):
        raise ValueError(
            f"PALS file '{path}' must contain a YAML/JSON mapping at the root, "
            f"got {type(data).__name__}"
        )
    if "PALS" not in data:
        raise ValueError(
            f"PALS file '{path}' is missing the required 'PALS' root key"
        )
    pals = data["PALS"]
    if not isinstance(pals, dict):
        raise ValueError(
            f"PALS file '{path}': 'PALS' must be a mapping, "
            f"got {type(pals).__name__}"
        )
    facility_list = pals.get("facility", [])

    # Build name→definition index from the facility list
    catalog = _build_catalog(facility_list)

    # Select beamline to expand
    bl_name = beamline_name or _find_default_beamline(catalog, facility_list)
    if bl_name not in catalog:
        raise ValueError(f"BeamLine '{bl_name}' not found in facility")

    # Expand beamline into a flat element sequence
    flat = _expand_beamline(bl_name, catalog)

    # Extract beam parameters
    beam_params = _extract_beam_params(catalog, facility_list, path,
                                       ke_override, particle_override)

    # Compute cumulative positions and normalize
    elements = _position_and_normalize(flat)

    return beam_params, elements


def _build_catalog(facility_list):
    """Index facility items by name → definition dict."""
    catalog = {}
    for item in facility_list:
        if not isinstance(item, dict):
            continue
        # Skip 'use' directives
        if "use" in item:
            continue
        for name, defn in item.items():
            if isinstance(defn, dict):
                catalog[name] = defn
    return catalog


def _find_default_beamline(catalog, facility_list):
    """Find the last BeamLine defined in the facility (not Lattice)."""
    # First check for a 'use' directive pointing to a Lattice
    for item in reversed(facility_list):
        if isinstance(item, dict) and "use" in item:
            lattice_name = item["use"]
            if lattice_name in catalog:
                lat = catalog[lattice_name]
                branches = lat.get("branches", [])
                if branches:
                    # First branch
                    br = branches[0]
                    if isinstance(br, str):
                        return br
                    elif isinstance(br, dict) and br:
                        return next(iter(br))

    # Fall back to last BeamLine in catalog
    last_bl = None
    for item in facility_list:
        if not isinstance(item, dict):
            continue
        for name, defn in item.items():
            if isinstance(defn, dict) and defn.get("kind") == "BeamLine":
                last_bl = name
    if last_bl:
        return last_bl
    raise ValueError("No BeamLine found in facility")


def _expand_beamline(name, catalog, _stack=None):
    """Recursively expand a BeamLine into a flat list of element dicts."""
    if _stack is None:
        _stack = set()
    if name in _stack:
        raise ValueError(f"Circular reference in BeamLine expansion: {name}")
    _stack = _stack | {name}

    defn = catalog[name]
    line = defn.get("line", [])
    result = []

    for item in line:
        if isinstance(item, str):
            # Reference to a named element
            result.extend(_resolve_reference(item, catalog, _stack))
        elif isinstance(item, dict):
            for elem_name, elem_defn in item.items():
                if elem_defn is None:
                    # Simple reference (name with null value, unusual but handle it)
                    result.extend(_resolve_reference(elem_name, catalog, _stack))
                elif isinstance(elem_defn, dict):
                    result.extend(
                        _resolve_inline(elem_name, elem_defn, catalog, _stack)
                    )
                else:
                    result.extend(_resolve_reference(elem_name, catalog, _stack))

    return result


def _resolve_reference(name, catalog, stack):
    """Resolve a plain string reference to an element or sub-beamline."""
    if name not in catalog:
        raise ValueError(f"Unresolved reference: '{name}'")

    defn = catalog[name]
    if defn.get("kind") == "BeamLine":
        return _expand_beamline(name, catalog, stack)
    else:
        return [{"name": name, **defn}]


def _resolve_inline(name, defn, catalog, stack):
    """Resolve an inline element definition (possibly with inherit/repeat)."""
    repeat = defn.get("repeat", 1)
    inherit_from = defn.get("inherit")

    if repeat == 0:
        warnings.warn(f"Element '{name}': repeat=0 produces no elements")
    elif repeat < 0:
        raise ValueError(
            f"Element '{name}': negative repeat ({repeat}) is not supported "
            f"(reversed expansion would require flipping bend angles, edge "
            f"angles, and element order, which is not implemented)"
        )
    if "direction" in defn:
        raise ValueError(
            f"Element '{name}': 'direction' modifier is not supported "
            f"(would require flipping bend signs and edge angles, "
            f"which is not implemented)"
        )

    # If the inline dict only has modifiers (repeat, inherit, etc.) and the name
    # exists in the catalog, treat it as a modified reference to the catalog entry.
    modifier_keys = {"repeat", "inherit", "direction"}
    if not inherit_from and set(defn.keys()) <= modifier_keys and name in catalog:
        cat_entry = catalog[name]
        if cat_entry.get("kind") == "BeamLine":
            expanded = _expand_beamline(name, catalog, stack)
            return [dict(e) for _ in range(abs(repeat)) for e in expanded]
        else:
            return [{"name": name, **cat_entry} for _ in range(abs(repeat))]

    if inherit_from:
        if inherit_from not in catalog:
            raise ValueError(f"Unresolved inherit reference: '{inherit_from}'")
        merged = {**catalog[inherit_from], **defn}
        # Merge nested dicts one level deep (e.g. MagneticMultipoleP, BendP).
        # Deeper nesting is not merged — no current PALS element uses it.
        for key in defn:
            if isinstance(defn[key], dict) and isinstance(catalog[inherit_from].get(key), dict):
                merged[key] = {**catalog[inherit_from][key], **defn[key]}
        merged.pop("inherit", None)
        merged.pop("repeat", None)
    else:
        merged = {k: v for k, v in defn.items() if k != "repeat"}

    kind = merged.get("kind")
    if kind == "BeamLine":
        catalog[name] = merged
        expanded = _expand_beamline(name, catalog, stack)
        return [dict(e) for _ in range(abs(repeat)) for e in expanded]
    elif name in catalog and catalog[name].get("kind") == "BeamLine":
        # Name references a catalog BeamLine but inline has extra keys
        extra = set(merged.keys()) - {"kind", "repeat", "inherit", "direction"}
        if extra:
            warnings.warn(
                f"Element '{name}': keys {extra} ignored on BeamLine reference"
            )
        expanded = _expand_beamline(name, catalog, stack)
        return [dict(e) for _ in range(abs(repeat)) for e in expanded]
    else:
        return [{"name": name, **merged} for _ in range(abs(repeat))]


# -----------------------------------------------------------------------
# Normalization: PALS elements → pals2cosy internal dict format
# -----------------------------------------------------------------------

# Official PALS kinds → internal short names
# Quadrupole, Kicker, Instrument are handled explicitly in _normalize_kind
_PALS_KIND_MAP = {
    "Drift": "DRIFT",
    "SBend": "DPH",
    "RBend": "DPH",
    "Wiggler": "UND",
    "Solenoid": "SOL",
    "RFCavity": "RFC",
    "Sextupole": "SXT",
    "Octupole": "OCT",
    "Multipole": "MULT",
    "Marker": "DRIFT",
}


def _to_float(value, name, field):
    try:
        f = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Element '{name}': field '{field}' must be numeric, got {value!r}"
        ) from exc
    if not math.isfinite(f):
        raise ValueError(
            f"Element '{name}': field '{field}' must be a finite number, got {f}"
        )
    return f


def _to_nonneg_float(value, name, field):
    f = _to_float(value, name, field)
    if f < 0:
        raise ValueError(
            f"Element '{name}': field '{field}' must be non-negative, got {f}"
        )
    return f


def _to_positive_float(value, name, field):
    f = _to_float(value, name, field)
    if f <= 0:
        raise ValueError(
            f"Element '{name}': field '{field}' must be positive, got {f}"
        )
    return f


def _position_and_normalize(flat_elements):
    """Assign cumulative s_start/s_end and normalize to internal dict format."""
    result = []
    s = 0.0

    for raw in flat_elements:
        name = raw.get("name", "<unnamed>")
        length = _to_nonneg_float(raw.get("length", 0.0), name, "length")
        kind = raw.get("kind")
        if kind is None:
            warnings.warn(
                f"Element '{name}' has no 'kind'; defaulting to Drift"
            )
            kind = "Drift"
        internal_type = _normalize_kind(kind, raw)

        elem = {
            "type": internal_type,
            "name": raw.get("name", ""),
            "length": length,
            "s_start": s,
            "s_end": s + length,
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

        if kind == "Quadrupole":
            mmp = raw.get("MagneticMultipoleP", {})
            bn1 = mmp.get("Bn1")
            if bn1 is not None:
                elem["bn1"] = _to_float(bn1, name, "MagneticMultipoleP.Bn1")

        elif kind == "Sextupole":
            mmp = raw.get("MagneticMultipoleP", {})
            bn2 = mmp.get("Bn2")
            if bn2 is not None:
                elem["bn2"] = _to_float(bn2, name, "MagneticMultipoleP.Bn2")

        elif kind == "Octupole":
            mmp = raw.get("MagneticMultipoleP", {})
            bn3 = mmp.get("Bn3")
            if bn3 is not None:
                elem["bn3"] = _to_float(bn3, name, "MagneticMultipoleP.Bn3")

        elif kind == "Multipole":
            mmp = raw.get("MagneticMultipoleP", {})
            for attr, key in [("bn1", "Bn1"), ("bn2", "Bn2"), ("bn3", "Bn3"),
                              ("bn4", "Bn4"), ("bn5", "Bn5")]:
                val = mmp.get(key)
                if val is not None:
                    elem[attr] = _to_float(val, name, f"MagneticMultipoleP.{key}")

        elif kind == "Wiggler":
            wig_p = raw.get("WigglerP", {})
            bw = wig_p.get("peak_field")
            if bw is not None:
                elem["wiggler_field"] = _to_float(bw, name, "WigglerP.peak_field")
            period = wig_p.get("period")
            if period is not None:
                elem["wiggler_period"] = _to_positive_float(period, name, "WigglerP.period")

        elif kind == "Solenoid":
            sol_p = raw.get("SolenoidP", {})
            bz = sol_p.get("Bz")
            if bz is not None:
                elem["bz"] = _to_float(bz, name, "SolenoidP.Bz")

        elif kind == "RFCavity":
            rf_p = raw.get("RFCavityP", {})
            v = rf_p.get("voltage_kv")
            if v is not None:
                elem["rf_voltage_kv"] = _to_float(v, name, "RFCavityP.voltage_kv")
            f = rf_p.get("frequency_hz")
            if f is not None:
                elem["rf_frequency_hz"] = _to_positive_float(f, name, "RFCavityP.frequency_hz")
            phi = rf_p.get("phase_deg")
            if phi is not None:
                elem["rf_phase_deg"] = _to_float(phi, name, "RFCavityP.phase_deg")

        elif kind in ("SBend", "RBend"):
            bend_p = raw.get("BendP", {})
            g_ref = _to_float(bend_p.get("g_ref", 0.0), name, "BendP.g_ref")
            angle_rad = g_ref * length
            elem["angle"] = math.degrees(angle_rad)
            elem["dipole_length"] = length

            e1 = _to_float(bend_p.get("e1", 0.0), name, "BendP.e1")
            e2 = _to_float(bend_p.get("e2", 0.0), name, "BendP.e2")
            elem["entrance_edge_angle"] = math.degrees(e1)
            elem["exit_edge_angle"] = math.degrees(e2)

            # Pole gap from aperture or explicit parameter
            aperture_p = raw.get("ApertureP", {})
            y_width = aperture_p.get("y_width", 0.0)
            elem["pole_gap"] = _to_nonneg_float(y_width, name, "ApertureP.y_width")

        result.append(elem)
        s += length

    return result


def _normalize_kind(kind, raw):
    """Map PALS kind to internal short name with polarity resolution.

    Note on QPF/QPD labeling: The sign convention here (Bn1 >= 0 → QPF)
    follows the proton convention. For electrons, the labels are inverted
    (positive Bn1 is actually defocusing). This is cosmetic only — the
    converter uses Bn1 directly for MQ commands, bypassing the current/sign
    path entirely.
    """
    if kind == "Quadrupole":
        mmp = raw.get("MagneticMultipoleP", {})
        raw_bn1 = mmp.get("Bn1", 0.0)
        bn1 = _to_float(raw_bn1, raw.get("name", "<unnamed>"),
                        "MagneticMultipoleP.Bn1")
        return "QPF" if bn1 >= 0 else "QPD"
    if kind == "Kicker":
        plane = raw.get("plane", "vertical")
        return "STV" if plane == "vertical" else "STH"
    if kind == "Instrument":
        return raw.get("instrument_type", "BPM")
    return _PALS_KIND_MAP.get(kind, kind)


def _extract_beam_params(catalog, facility_list, path, ke_override, particle_override):
    """Extract beam parameters from Lattice definition or CLI overrides."""
    # Try to find a Lattice with particle info
    particle_type = particle_override.lower() if particle_override else None
    ke = ke_override
    mass_mev = E0_PROTON

    for item in facility_list:
        if not isinstance(item, dict):
            continue
        for name, defn in item.items():
            if not isinstance(defn, dict):
                continue
            if defn.get("kind") == "Lattice":
                particle = defn.get("particle", {})
                if isinstance(particle, dict):
                    species = particle.get("species", particle.get("species_ref"))
                    if species and not particle_override:
                        particle_type = species.lower()
                    energy = particle.get("kinetic_energy")
                    if energy is not None and ke is None:
                        # PALS uses eV; convert to MeV
                        ke = float(energy) / 1e6

    if particle_type is None:
        particle_type = "proton"
        warnings.warn(
            "Particle type not specified in lattice file or CLI; defaulting to 'proton'. "
            "Use --particle to specify explicitly."
        )

    mass_mev = PARTICLE_MASSES_MEV.get(particle_type)
    if mass_mev is None:
        warnings.warn(
            f"Unknown particle species '{particle_type}'; using proton mass as fallback"
        )
        mass_mev = E0_PROTON

    if ke is None:
        raise ValueError(
            "Kinetic energy not found in lattice file. "
            "Use --ke to specify kinetic energy in MeV."
        )

    return {
        "kinetic_energy_mev": ke,
        "particle_type": particle_type,
        "mass_mev": mass_mev,
        "rf_frequency_hz": F_RF_DEFAULT,
        "quad_gradient": None,  # Not needed for Bn1 quads
        "source_file": os.path.basename(path),
    }
