"""Convert a parsed PALS lattice to COSY INFINITY COSYScript.

Author: Eremey Valetov
"""

import math
import re
import warnings

from .constants import G_QUAD_DEFAULT

# Element types that produce no COSYScript optics commands (emitted as DL if length > 0)
_PASSIVE_TYPES = frozenset({
    "BPM", "OTR", "STV", "STH", "SPC", "XRS", "BSW",
})


def convert(beam_params, elements,
            fringe_field_order=0,
            computation_order=3,
            dimensions=3,
            quad_aperture=0.027,
            ke_override=None,
            particle_override=None,
            comments=True,
            twiss=False):
    """Convert parsed lattice to COSY INFINITY COSYScript.

    Parameters
    ----------
    beam_params : dict
        From parser.parse_lattice().
    elements : list[dict]
        Normalized element list from parser.parse_lattice().
    fringe_field_order : int
        FR command value (0, 1, 2, or 3).
    computation_order : int
        DA computation order (OV).
    dimensions : int
        Phase-space dimensions (2 or 3).
    quad_aperture : float
        Quadrupole bore diameter in meters.
    ke_override : float or None
        Override kinetic energy from lattice file.
    particle_override : str or None
        Override particle type.
    comments : bool
        Include element name comments in COSYScript.
    twiss : bool
        Append GT Twiss extraction and JSON output (requires periodic lattice).

    Returns
    -------
    str
        Complete COSYScript.
    """
    ke = ke_override if ke_override is not None else beam_params["kinetic_energy_mev"]
    particle = particle_override or beam_params.get("particle_type", "electron")
    G = beam_params.get("quad_gradient") or G_QUAD_DEFAULT
    r = quad_aperture / 2
    source = beam_params.get("source_file", "")

    consolidated = _consolidate_dipoles(elements)
    # FC (Enge coefficients) are used by FR 1 (soft-edge) and FR 3 (high-precision),
    # but NOT by FR 2 (symplectic scaling, which uses SYSCA.DAT data instead).
    body = _generate_lattice_body(consolidated, G, r, comments,
                                  emit_fc=fringe_field_order in (1, 3))

    return _fox_template(
        source_file=source,
        dimensions=dimensions,
        order=computation_order,
        ke=ke,
        fr=fringe_field_order,
        lattice_body=body,
        particle=particle,
        twiss=twiss,
    )


# ---------------------------------------------------------------------------
# DPW-DPH-DPW consolidation
# ---------------------------------------------------------------------------

def _float_or(val, default=0.0):
    """Return val as float, or default if None."""
    return float(val) if val is not None else default


def _consolidate_dipoles(elements):
    """Detect DPW-DPH-DPW triplets and merge into DIPOLE_CONSOLIDATED."""
    result = []
    i = 0
    n = len(elements)

    while i < n:
        elem = elements[i]

        if (elem["type"] == "DPW"
                and i + 2 < n
                and elements[i + 1]["type"] == "DPH"
                and elements[i + 2]["type"] == "DPW"):
            entrance = elements[i]
            main = elements[i + 1]
            exit_wedge = elements[i + 2]

            result.append({
                "type": "DIPOLE_CONSOLIDATED",
                "name": main.get("name", ""),
                "length": main["length"],
                "angle": _float_or(main.get("angle"), 0.0),
                "entrance_angle": _float_or(entrance.get("wedge_angle"), 0.0),
                "exit_angle": _float_or(exit_wedge.get("wedge_angle"), 0.0),
                "pole_gap": _float_or(main.get("pole_gap"), 0.0),
                # Position-based Enge assignment (more robust than sign detection)
                "entrance_enge": entrance.get("enge_coeffs"),
                "exit_enge": exit_wedge.get("enge_coeffs"),
            })
            i += 3

        elif elem["type"] == "DPH":
            # PALS-native SBend with edge angles on the element itself
            result.append({
                "type": "DIPOLE_CONSOLIDATED",
                "name": elem.get("name", ""),
                "length": elem["dipole_length"] if elem.get("dipole_length") is not None else elem["length"],
                "angle": _float_or(elem.get("angle"), 0.0),
                "entrance_angle": _float_or(elem.get("entrance_edge_angle"), 0.0),
                "exit_angle": _float_or(elem.get("exit_edge_angle"), 0.0),
                "pole_gap": _float_or(elem.get("pole_gap"), 0.0),
                "entrance_enge": elem.get("enge_coeffs"),
                "exit_enge": elem.get("enge_coeffs"),  # symmetric: same Enge for both faces
            })
            i += 1

        else:
            result.append(elem)
            i += 1

    return result


# ---------------------------------------------------------------------------
# COSYScript generation
# ---------------------------------------------------------------------------

def _generate_lattice_body(elements, G, r, comments, emit_fc=True):
    """Generate the body of the LATTICE procedure."""
    lines = []

    def _comment(elem):
        """Return inline comment string if comments enabled and element has a name."""
        if comments and elem.get("name"):
            safe = re.sub(r'[^\x20-\x7e]|[{}]', '', elem["name"])
            return f" {{ {safe} }}"
        return ""

    for elem in elements:
        etype = elem["type"]

        if etype == "DRIFT":
            if elem["length"] > 0:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        elif etype in ("QPF", "QPD"):
            if elem["length"] <= 0:
                continue
            bn1 = elem.get("bn1")
            if bn1 is not None:
                b_pole = bn1
            else:
                sign = -1 if etype == "QPF" else 1
                current = elem.get("current") or 0
                b_pole = sign * G * current * r
            lines.append(f"    MQ {elem['length']} {b_pole} {r} ;{_comment(elem)}")

        elif etype == "DIPOLE_CONSOLIDATED":
            if elem["length"] <= 0:
                continue
            comment = _comment(elem)
            lines.extend(_dipole_fox(elem, comment, emit_fc=emit_fc))

        elif etype == "SXT":
            if elem["length"] <= 0:
                continue
            bn2 = elem.get("bn2")
            if bn2 is not None:
                lines.append(f"    MH {elem['length']} {bn2} {r} ;{_comment(elem)}")
            else:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        elif etype == "OCT":
            if elem["length"] <= 0:
                continue
            bn3 = elem.get("bn3")
            if bn3 is not None:
                lines.append(f"    MO {elem['length']} {bn3} {r} ;{_comment(elem)}")
            else:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        elif etype == "SOL":
            if elem["length"] <= 0:
                continue
            bz = elem.get("bz")
            if bz is not None:
                lines.append(f"    CMS {bz} {r} {elem['length']} ;{_comment(elem)}")
            else:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        elif etype == "RFC":
            if elem["length"] <= 0:
                continue
            v_kv = elem.get("rf_voltage_kv")
            f_hz = elem.get("rf_frequency_hz")
            phi = elem.get("rf_phase_deg", 0.0)
            if v_kv is not None and f_hz is not None:
                # RF V N W PHI D: V=voltage (kV), N=polynomial order (0=uniform),
                # W=frequency (Hz), PHI=phase (deg), D=half-gap
                lines.append(f"    RF {v_kv} 0 {f_hz} {phi} {r} ;{_comment(elem)}")
            else:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        elif etype == "MULT":
            if elem["length"] <= 0:
                continue
            bns = [elem.get(f"bn{i}", 0) or 0 for i in range(1, 6)]
            if any(b != 0 for b in bns):
                lines.append(
                    f"    M5 {elem['length']} {bns[0]} {bns[1]} {bns[2]} "
                    f"{bns[3]} {bns[4]} {r} ;{_comment(elem)}")
            else:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        elif etype == "UND":
            if elem["length"] <= 0:
                continue
            bw = elem.get("wiggler_field")
            period = elem.get("wiggler_period")
            if bw is not None and period is not None:
                if period <= 0:
                    raise ValueError(
                        f"Wiggler '{elem.get('name', '')}': period must be positive, "
                        f"got {period}"
                    )
                k = 2 * math.pi / period
                lines.append(
                    f"    WI {bw} {k} {elem['length']} {r} 0 0 0 ;{_comment(elem)}")
            else:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        elif etype == "DPW":
            # Lone DPW not part of a triplet — emit as drift with warning
            warnings.warn(
                f"Standalone DPW '{elem.get('name', '')}' cannot be converted to COSYScript; "
                f"emitting as drift of length {elem['length']}"
            )
            if elem["length"] > 0:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        elif etype in _PASSIVE_TYPES:
            if elem["length"] > 0:
                lines.append(f"    DL {elem['length']} ;{_comment(elem)}")

        else:
            warnings.warn(
                f"Unknown element type '{etype}' (name='{elem.get('name', '')}'); skipping"
            )

    return "\n".join(lines)


def _dipole_fox(elem, comment="", emit_fc=True):
    """Generate COSYScript lines for a consolidated dipole."""
    lines = []
    angle = elem["angle"]
    d_half = elem["pole_gap"] / 2
    e1 = elem["entrance_angle"]
    e2 = elem["exit_angle"]
    entrance_enge = elem.get("entrance_enge")
    exit_enge = elem.get("exit_enge")
    has_fc = False

    # Fringe field coefficients (only when FR > 0)
    if emit_fc:
        if entrance_enge:
            coeffs = _format_enge(entrance_enge)
            lines.append(f"    FC 1 1 1 {coeffs} ;")
            has_fc = True
        if exit_enge:
            coeffs = _format_enge(exit_enge)
            lines.append(f"    FC 1 2 1 {coeffs} ;")
            has_fc = True

    # CB wrap for negative bending angle
    use_cb = angle < 0
    if use_cb:
        lines.append("    CB ;")

    lines.append(f"    DIL {elem['length']} {abs(angle)} {d_half} {e1} 0 {e2} 0 ;{comment}")

    if use_cb:
        lines.append("    CB ;")

    if has_fc:
        lines.append("    FD ;")

    return lines


def _format_enge(coeffs):
    """Pad/truncate Enge coefficients to exactly 6 values."""
    if len(coeffs) > 6:
        warnings.warn(
            f"Enge coefficients truncated from {len(coeffs)} to 6 "
            f"(COSY INFINITY FC accepts max 6)"
        )
    padded = (list(coeffs) + [0.0] * 6)[:6]
    return " ".join(str(float(c)) for c in padded)


# ---------------------------------------------------------------------------
# COSYScript template
# ---------------------------------------------------------------------------

# Particles with dedicated COSY preset procedures.
# COSY signatures (cosy.fox lines 11508–11514):
#   RPP E       → RP E 1.00727646688  +1  (proton)
#   RPE E       → RP E 5.485799110e-4 -1  (electron)
#   RPMU E      → RP E 0.1134289168   -1  (μ⁻, charge hardcoded to -1)
#   RPPI E Z    → RP E 0.149834758     Z  (pion, charge argument)
_PARTICLE_PRESETS = {
    "electron":  ("RPE",  None),   # RPE E
    "proton":    ("RPP",  None),   # RPP E
    "muon":      ("RPMU", None),   # RPMU E (μ⁻)
    "mu-":       ("RPMU", None),   # alias
    "pion+":     ("RPPI", +1),     # RPPI E 1
    "pion-":     ("RPPI", -1),     # RPPI E -1
    "pi+":       ("RPPI", +1),     # alias
    "pi-":       ("RPPI", -1),     # alias
}

# Generic RP fallback: particle name → (mass_amu, charge).
# Used for particles without a dedicated COSY preset, or where the preset
# hardcodes the wrong charge sign (e.g. μ⁺ needs RP, not RPMU which is μ⁻).
_PARTICLE_GENERIC = {
    "mu+":       (0.1134289168,    +1),   # antimuon — RPMU hardcodes charge -1
    "antimuon":  (0.1134289168,    +1),
    "positron":  (5.485799110e-4,  +1),
    "antiproton": (1.00727646688,  -1),
    "deuteron":  (2.01355321271,   +1),
    "alpha":     (4.00150617912,   +2),
    "carbon12":  (12.0,            +6),
}


def _particle_command(particle, ke):
    """Return the COSYScript command string to set the reference particle."""
    name = particle.lower()
    if name in _PARTICLE_PRESETS:
        cmd, charge = _PARTICLE_PRESETS[name]
        if charge is not None:
            return f"{cmd} {ke} {charge}"
        return f"{cmd} {ke}"
    if name in _PARTICLE_GENERIC:
        mass_amu, charge = _PARTICLE_GENERIC[name]
        return f"RP {ke} {mass_amu} {charge}"
    supported = sorted(set(list(_PARTICLE_PRESETS) + list(_PARTICLE_GENERIC)))
    raise ValueError(
        f"Unknown particle type '{particle}'. "
        f"Supported: {', '.join(supported)}"
    )


def _fox_template(source_file, dimensions, order, ke, fr, lattice_body,
                  particle="electron", twiss=False):
    dim = dimensions
    rp_line = _particle_command(particle, ke)
    safe_source = re.sub(r'[^\x20-\x7e]|[{}]', '', source_file)
    # No explicit FD needed — OV calls FD internally, loading default
    # Enge coefficients (PEP/SLAC) for all multipoles.

    twiss_vars = f"""\
    VARIABLE A0 100 {dim} ; VARIABLE B0 100 {dim} ; VARIABLE G0 100 {dim} ;
    VARIABLE R0 100 {dim} ; VARIABLE MU0 100 {dim} ; VARIABLE F0 100 {2*dim} ;
""" if twiss else ""

    twiss_block = """
    CO 1 ; PM 99 ;

    GT MAP F0 MU0 A0 B0 G0 R0 ;

    OPENF 51 'result.txt' 'UNKNOWN' ;
        WRITE 51 '{' ;
        WRITE 51 '"twiss": {' ;
           WRITE 51 '  "beta_x": "'&S(CONS(B0(1)))&'",' ;
           WRITE 51 '  "beta_y": "'&S(CONS(B0(2)))&'",' ;
           WRITE 51 '  "alpha_x": "'&S(CONS(A0(1)))&'",' ;
           WRITE 51 '  "alpha_y": "'&S(CONS(A0(2)))&'",' ;
           WRITE 51 '  "gamma_x": "'&S(CONS(G0(1)))&'",' ;
           WRITE 51 '  "gamma_y": "'&S(CONS(G0(2)))&'",' ;
           WRITE 51 '  "mu_x": "'&S(CONS(MU0(1)))&'",' ;
           WRITE 51 '  "mu_y": "'&S(CONS(MU0(2)))&'"' ;
        WRITE 51 '}' ;
        WRITE 51 '}' ;
    CLOSEF 51 ;""" if twiss else ""

    return f"""\
{{ Generated by pals2cosy from {safe_source} }}
INCLUDE 'COSY' ;
PROCEDURE RUN ;
{twiss_vars}\
    PROCEDURE LATTICE ;
        UM ; CR ;
{lattice_body}
    ENDPROCEDURE ;

    OV {order} {dim} 0 ;
    {rp_line} ;
    FR {fr} ;
    LATTICE ;
{twiss_block}
ENDPROCEDURE ;
RUN ;
END ;
"""
