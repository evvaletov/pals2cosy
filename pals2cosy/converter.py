"""Convert a parsed PALS lattice to COSY INFINITY FOX code.

Handles both PALS-native SBend/RBend (with edge angle attributes) and
FELsim DIPOLE_WEDGE (DPW-DPH-DPW triplet) representations. Both produce
identical FOX output.

Author: Eremey Valetov
"""

from .constants import G_QUAD_DEFAULT

# Element types that produce no FOX optics commands
_PASSIVE_TYPES = frozenset({
    "BPM", "OTR", "STV", "STH", "SPC", "XRS", "BSW",
    "SOL", "RFC", "SXT",  # stubs — no COSY mapping yet
})


def convert(beam_params, elements,
            fringe_field_order=0,
            computation_order=3,
            dimensions=3,
            quad_aperture=0.027,
            ke_override=None,
            particle_override=None,
            comments=True):
    """Convert parsed lattice to COSY INFINITY FOX code.

    Parameters
    ----------
    beam_params : dict
        From parser.parse_lattice().
    elements : list[dict]
        Normalized element list from parser.parse_lattice().
    fringe_field_order : int
        FR command value (0 or 1).
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
        Include element name comments in FOX code.

    Returns
    -------
    str
        Complete FOX code.
    """
    ke = ke_override if ke_override is not None else beam_params["kinetic_energy_mev"]
    G = beam_params.get("quad_gradient", G_QUAD_DEFAULT)
    r = quad_aperture / 2
    source = beam_params.get("source_file", "")

    consolidated = _consolidate_dipoles(elements)
    body = _generate_lattice_body(consolidated, G, r, comments)

    return _fox_template(
        source_file=source,
        dimensions=dimensions,
        order=computation_order,
        ke=ke,
        fr=fringe_field_order,
        lattice_body=body,
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
                "length": elem.get("dipole_length") or elem["length"],
                "angle": _float_or(elem.get("angle"), 0.0),
                "entrance_angle": _float_or(elem.get("entrance_edge_angle"), 0.0),
                "exit_angle": _float_or(elem.get("exit_edge_angle"), 0.0),
                "pole_gap": _float_or(elem.get("pole_gap"), 0.0),
                "entrance_enge": elem.get("enge_coeffs"),
                "exit_enge": None,
            })
            i += 1

        else:
            result.append(elem)
            i += 1

    return result


# ---------------------------------------------------------------------------
# FOX code generation
# ---------------------------------------------------------------------------

def _generate_lattice_body(elements, G, r, comments):
    """Generate the body of the LATTICE procedure."""
    lines = []

    for elem in elements:
        etype = elem["type"]

        if etype == "DRIFT":
            if elem["length"] > 0:
                lines.append(f"    DL {elem['length']} ;")

        elif etype in ("QPF", "QPD"):
            sign = -1 if etype == "QPF" else 1
            current = elem.get("current") or 0
            b_pole = sign * G * current * r
            lines.append(f"    MQ {elem['length']} {b_pole} {r} ;")

        elif etype == "DIPOLE_CONSOLIDATED":
            lines.extend(_dipole_fox(elem))

        elif etype == "UND":
            # Undulator → drift
            if elem["length"] > 0:
                lines.append(f"    DL {elem['length']} ;")

        elif etype in _PASSIVE_TYPES:
            # Zero-length diagnostics/correctors — skip
            if elem["length"] > 0:
                lines.append(f"    DL {elem['length']} ;")

        # else: unknown type, skip silently

    return "\n".join(lines)


def _dipole_fox(elem):
    """Generate FOX lines for a consolidated dipole."""
    lines = []
    angle = elem["angle"]
    d_half = elem["pole_gap"] / 2
    e1 = elem["entrance_angle"]
    e2 = elem["exit_angle"]
    entrance_enge = elem.get("entrance_enge")
    exit_enge = elem.get("exit_enge")
    has_fc = False

    # Fringe field coefficients
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

    lines.append(f"    DIL {elem['length']} {abs(angle)} {d_half} {e1} 0 {e2} 0 ;")

    if use_cb:
        lines.append("    CB ;")

    if has_fc:
        lines.append("    FD ;")

    return lines


def _format_enge(coeffs):
    """Pad/truncate Enge coefficients to exactly 6 values."""
    padded = (list(coeffs) + [0.0] * 6)[:6]
    return " ".join(str(c) for c in padded)


# ---------------------------------------------------------------------------
# FOX template
# ---------------------------------------------------------------------------

def _fox_template(source_file, dimensions, order, ke, fr, lattice_body):
    dim = dimensions
    return f"""\
{{ Generated by pals2cosy from {source_file} }}
INCLUDE 'COSY' ;
PROCEDURE RUN ;
    VARIABLE A0 100 {dim} ; VARIABLE B0 100 {dim} ; VARIABLE G0 100 {dim} ;
    VARIABLE R0 100 {dim} ; VARIABLE MU0 100 {dim} ; VARIABLE F0 100 {2*dim} ;

    PROCEDURE LATTICE ;
        UM ; CR ;
{lattice_body}
    ENDPROCEDURE ;

    OV {order} {dim} 0 ;
    RPE {ke} ;
    FR {fr} ;
    LATTICE ;

    CO 1 ; PM 99 ;

    GT MAP F0 MU0 A0 B0 G0 R0 ;

    OPENF 51 'result.txt' 'UNKNOWN' ;
        WRITE 51 '{{' ;
        WRITE 51 '"twiss": {{' ;
           WRITE 51 '  "beta_x": "'&S(CONS(B0(1)))&'",' ;
           WRITE 51 '  "beta_y": "'&S(CONS(B0(2)))&'",' ;
           WRITE 51 '  "alpha_x": "'&S(CONS(A0(1)))&'",' ;
           WRITE 51 '  "alpha_y": "'&S(CONS(A0(2)))&'",' ;
           WRITE 51 '  "gamma_x": "'&S(CONS(G0(1)))&'",' ;
           WRITE 51 '  "gamma_y": "'&S(CONS(G0(2)))&'",' ;
           WRITE 51 '  "mu_x": "'&S(CONS(MU0(1)))&'",' ;
           WRITE 51 '  "mu_y": "'&S(CONS(MU0(2)))&'"' ;
        WRITE 51 '}}' ;
        WRITE 51 '}}' ;
    CLOSEF 51 ;
ENDPROCEDURE ;
RUN ;
END ;
"""
