"""Convert a parsed PALS lattice to COSY INFINITY COSYScript.

Author: Eremey Valetov
"""

import math
import re
import warnings

from ._elements import MAX_MULTIPOLE_ORDER
from ._version import __pals_commit__, __version__
from .constants import G_QUAD_DEFAULT
from .diagnostics import PALSConversionWarning
from .species import parse_species

# Element types that produce no COSYScript optics commands (emitted as DL if length > 0)
_PASSIVE_TYPES = frozenset({
    "BPM", "OTR", "STV", "STH", "SPC", "XRS", "BSW",
    "KICKER", "INSTRUMENT", "MARKER",
})

_MULTIPOLE_TYPES = frozenset({"QPF", "QPD", "SXT", "OCT", "MULT"})

# Zero-length (thin) multipoles are modelled as a slice of this length
# between drifts of minus half its length, which keeps the total length and
# reproduces the thin kick to O(K L * THIN_SLICE).
THIN_SLICE = 1e-6

# COSY's CMS solenoid has tanh fringe fields of length D, which reduce its
# focusing relative to the hard-edge PALS solenoid by about D/L. Below this
# ratio the difference is not reported.
SOLENOID_FRINGE_TOLERANCE = 1e-3

# Integral of F (1 - F) over S for COSY's default dipole Enge function F
# (procedure FD), where S is the distance from the effective field boundary
# in units of the full gap 2d. A bend of half-gap d therefore has the PALS
# fringe-field integral edge_int = d * DIPOLE_ENGE_INTEGRAL.
DIPOLE_ENGE_INTEGRAL = 0.3694585

_MM_NORMAL = "PALSBN"
_MM_SKEW = "PALSBS"
_RF_VOLTAGE = "PALSRFV"


def convert(beam_params, elements,
            fringe_field_order=0,
            computation_order=3,
            dimensions=None,
            quad_aperture=0.027,
            ke_override=None,
            particle_override=None,
            comments=True,
            twiss=False,
            solenoid_fringe=None):
    """Convert parsed lattice to COSY INFINITY COSYScript.

    Parameters
    ----------
    beam_params : dict
        From parser.parse_lattice() or pals_parser.parse_lattice().
    elements : list[dict]
        Normalized element list from the same parser.
    fringe_field_order : int
        FR command value (0, 1, 2, or 3).
    computation_order : int
        DA computation order (OV).
    dimensions : int or None
        Phase-space dimensions (2 or 3). Default: 3, or 2 with ``twiss``.
    quad_aperture : float
        Multipole bore diameter in meters. Its radius is the aperture at which
        COSY's pole-tip fields are specified, and the half-gap of bends that
        give no ``ApertureP`` when fringe fields are on.
    ke_override : float or None
        Override kinetic energy from lattice file.
    particle_override : str or None
        Override particle type.
    comments : bool
        Include element name comments in COSYScript.
    twiss : bool
        Append GT Twiss extraction and JSON output (requires periodic lattice).
    solenoid_fringe : float or None
        Fringe-field length [m] of COSY's soft-edge solenoid ``CMS``. Default:
        the bore radius ``quad_aperture / 2``.

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
    if solenoid_fringe is None:
        solenoid_fringe = r
    elif not solenoid_fringe > 0:
        raise ValueError(f"solenoid_fringe must be positive, got {solenoid_fringe}")

    overridden = (
        (ke_override is not None and ke_override != beam_params.get("kinetic_energy_mev"))
        or (particle_override is not None
            and _canonical(particle_override) != _canonical(beam_params.get("particle_type"))))
    if overridden and any(e.get("normalized") for e in elements):
        warnings.warn(
            "convert() overrides the reference particle, but normalized strengths "
            "were turned into fields with the parser's reference; pass the same "
            "overrides to the parser", PALSConversionWarning, stacklevel=2)

    consolidated = _consolidate_dipoles(elements)
    # FC (Enge coefficients) are used by FR 1 (soft-edge) and FR 3 (high-precision),
    # but NOT by FR 2 (symplectic scaling, which uses SYSCA.DAT data instead).
    writer = _LatticeWriter(G, r, comments, fringe_field_order, solenoid_fringe)
    body = writer.write(consolidated)

    if dimensions is None:
        dimensions = 2 if twiss else 3
    if twiss and dimensions == 3 and not writer.needs_rf_array:
        raise ValueError(
            "Twiss extraction with a 3D map needs RF: without longitudinal focusing "
            "COSY's GT cannot compute the Twiss parameters; use 2 dimensions")

    return _fox_template(
        source_file=source,
        dimensions=dimensions,
        order=computation_order,
        ke=ke,
        fr=fringe_field_order,
        lattice_body=body,
        declarations=writer.declarations(),
        particle=particle,
        twiss=twiss,
    )


def _canonical(particle):
    try:
        return parse_species(particle).name
    except ValueError:
        return particle


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
            # PALS-native bend with edge angles on the element itself
            result.append({
                "type": "DIPOLE_CONSOLIDATED",
                "name": elem.get("name", ""),
                "length": elem["dipole_length"] if elem.get("dipole_length") is not None else elem["length"],
                "angle": _float_or(elem.get("angle"), 0.0),
                "entrance_angle": _float_or(elem.get("entrance_edge_angle"), 0.0),
                "exit_angle": _float_or(elem.get("exit_edge_angle"), 0.0),
                "entrance_curvature": _float_or(elem.get("entrance_edge_curvature"), 0.0),
                "exit_curvature": _float_or(elem.get("exit_edge_curvature"), 0.0),
                "tilt": _float_or(elem.get("tilt"), 0.0),
                "pole_gap": _float_or(elem.get("pole_gap"), 0.0),
                "edge_int": elem.get("edge_int"),
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

def pole_tip_field(bn, order, radius):
    """COSY pole-tip field [T] of a PALS multipole coefficient ``B_N`` [T/m^N].

    PALS defines ``By + i Bx = (Bn + i Bs) (x + i y)^N / N!``, so the field
    magnitude at radius r is ``Bn r^N / N!``, which is what COSY's
    MQ/MH/MO/M5/MMS take at aperture r.
    """
    return bn * radius ** order / math.factorial(order)


class _LatticeWriter:
    """Generates the body of the LATTICE procedure."""

    def __init__(self, G, r, comments, fringe_field_order, solenoid_fringe):
        self.G = G
        self.r = r
        self.comments = comments
        self.fr = fringe_field_order
        self.solenoid_fringe = solenoid_fringe
        self.emit_fc = fringe_field_order in (1, 3)
        self.needs_multipole_arrays = False
        self.needs_rf_array = False

    def declarations(self):
        lines = []
        if self.needs_multipole_arrays:
            lines.append(f"VARIABLE {_MM_NORMAL} 1 {MAX_MULTIPOLE_ORDER} ; "
                         f"VARIABLE {_MM_SKEW} 1 {MAX_MULTIPOLE_ORDER} ;")
        if self.needs_rf_array:
            lines.append(f"VARIABLE {_RF_VOLTAGE} 1 1 1 ;")
        return lines

    def _comment(self, elem):
        if self.comments and elem.get("name"):
            safe = re.sub(r'[^\x20-\x7e]|[{}]', '', elem["name"])
            return f" {{ {safe} }}"
        return ""

    def _drift(self, elem, length=None):
        length = elem["length"] if length is None else length
        return [f"    DL {length} ;{self._comment(elem)}"] if length > 0 else []

    def write(self, elements):
        lines = []
        for elem in elements:
            lines.extend(self._element(elem))
        return "\n".join(lines)

    def _element(self, elem):
        etype = elem["type"]
        if etype == "DRIFT" or etype in _PASSIVE_TYPES:
            return self._drift(elem)
        if etype in _MULTIPOLE_TYPES:
            return self._multipole(elem)
        if etype == "DIPOLE_CONSOLIDATED":
            if elem["length"] <= 0:
                return []
            return _dipole_fox(elem, self._half_gap(elem), self._comment(elem),
                               emit_fc=self.emit_fc)
        if etype == "SOL":
            return self._solenoid(elem)
        if etype == "RFC":
            return self._rf(elem)
        if etype == "UND":
            return self._wiggler(elem)
        if etype == "DPW":
            # Lone DPW not part of a triplet — emit as drift with warning
            warnings.warn(
                f"Standalone DPW '{elem.get('name', '')}' cannot be converted to COSYScript; "
                f"emitting as drift of length {elem['length']}",
                PALSConversionWarning,
            )
            return self._drift(elem)
        warnings.warn(
            f"Unknown element type '{etype}' (name='{elem.get('name', '')}'); "
            f"emitting as a drift",
            PALSConversionWarning,
        )
        return self._drift(elem)

    def _half_gap(self, elem):
        """Half-gap of a bend, which sets the extent of COSY's fringe field.

        With fringe fields on, PALS fringe-field integrals take precedence:
        the half-gap is chosen so that COSY's default Enge field has them.
        """
        name = elem.get("name", "")
        edge_int = elem.get("edge_int") or (0.0, 0.0)
        if any(edge_int) and self.fr == 0:
            for face, value in enumerate(edge_int, 1):
                if value:
                    warnings.warn(
                        f"Bend '{name}': BendP.edge{face}_int is not represented at "
                        f"--fr 0, where COSY's bends have hard edges; --fr 1, 2 or 3 "
                        f"includes it", PALSConversionWarning)
        elif any(edge_int):
            if edge_int[0] != edge_int[1]:
                warnings.warn(
                    f"Bend '{name}': COSY's DIL has one fringe-field extent for both "
                    f"faces; it is set from the mean of edge1_int = {edge_int[0]} m and "
                    f"edge2_int = {edge_int[1]} m", PALSConversionWarning)
            return (edge_int[0] + edge_int[1]) / 2 / DIPOLE_ENGE_INTEGRAL

        half_gap = (elem.get("pole_gap") or 0.0) / 2
        if half_gap == 0 and self.fr > 0:
            warnings.warn(
                f"Bend '{name}' gives no pole gap (ApertureP.y_width) and no "
                f"edge1_int/edge2_int; its fringe field uses the bore radius {self.r} m "
                f"as the half-gap", PALSConversionWarning)
            half_gap = self.r
        return half_gap

    def _solenoid(self, elem):
        length = elem["length"]
        bz = elem.get("bz")
        if length <= 0:
            return []
        if bz is None:
            return self._drift(elem)
        fringe = self.solenoid_fringe
        if self.fr == 0 and fringe / length > SOLENOID_FRINGE_TOLERANCE:
            warnings.warn(
                f"Solenoid '{elem.get('name', '')}': COSY's CMS has tanh fringe fields "
                f"of length {fringe} m, so its focusing is about {fringe / length:.2%} "
                f"weaker than the hard-edge solenoid; a shorter --solenoid-fringe "
                f"reduces the difference at the cost of integration time",
                PALSConversionWarning,
            )
        return [f"    CMS {bz} {fringe} {length} ;{self._comment(elem)}"]

    def _multipole(self, elem):
        thin = bool(elem.get("thin"))
        length = THIN_SLICE if thin else elem["length"]
        if length <= 0:
            return []
        scale = 1 / THIN_SLICE if thin else 1.0
        normal, skew = {}, {}
        for order, b in (elem.get("b_pole") or {}).items():
            if b is not None:
                normal[order] = float(b)
        for order, b in (elem.get("bn") or {}).items():
            normal[order] = normal.get(order, 0.0) + pole_tip_field(b * scale, order, self.r)
        for order, b in (elem.get("bs") or {}).items():
            skew[order] = pole_tip_field(b * scale, order, self.r)

        etype = elem["type"]
        if not normal and not skew and etype in ("QPF", "QPD") and elem.get("current") is not None:
            sign = -1 if etype == "QPF" else 1
            normal[1] = sign * self.G * elem["current"] * self.r

        comment = self._comment(elem)
        nonzero = {k for k, v in normal.items() if v} | {k for k, v in skew.items() if v}
        if not nonzero and not (etype in ("QPF", "QPD") and 1 in normal):
            return [] if thin else self._drift(elem)

        r = self.r
        if skew or max(nonzero, default=0) > 5:
            nb = max(nonzero)
            self.needs_multipole_arrays = True
            assign = " ".join(
                f"{_MM_NORMAL}({i}) := {normal.get(i, 0.0)} ; "
                f"{_MM_SKEW}({i}) := {skew.get(i, 0.0)} ;"
                for i in range(1, nb + 1))
            cmd = [f"    {assign}",
                   f"    MMS {length} {_MM_NORMAL} {_MM_SKEW} {nb} {r} ;{comment}"]
        elif etype in ("QPF", "QPD") and nonzero <= {1}:
            cmd = [f"    MQ {length} {normal.get(1, 0.0)} {r} ;{comment}"]
        elif etype == "SXT" and nonzero == {2}:
            cmd = [f"    MH {length} {normal[2]} {r} ;{comment}"]
        elif etype == "OCT" and nonzero == {3}:
            cmd = [f"    MO {length} {normal[3]} {r} ;{comment}"]
        else:
            b = [normal.get(i, 0) for i in range(1, 6)]
            cmd = [f"    M5 {length} {b[0]} {b[1]} {b[2]} {b[3]} {b[4]} {r} ;{comment}"]

        if not thin:
            return cmd
        half = -THIN_SLICE / 2
        lines = []
        if self.fr != 0:
            lines.append("    FR 0 ;")
        lines += [f"    DL {half} ;", *cmd, f"    DL {half} ;"]
        if self.fr != 0:
            lines.append(f"    FR {self.fr} ;")
        return lines

    def _rf(self, elem):
        v_kv = elem.get("rf_voltage_kv")
        f_hz = elem.get("rf_frequency_hz")
        if v_kv is None or f_hz is None or v_kv == 0:
            return self._drift(elem)
        phi = elem.get("rf_phase_deg") or 0.0
        half = elem["length"] / 2
        self.needs_rf_array = True
        # RF V N W PHI D: V = voltage coefficient array (kV), N = polynomial
        # order (0 = uniform), W = frequency (Hz), PHI = phase (deg),
        # D = aperture. COSY's RF is a thin kick, so it sits mid-cavity.
        lines = [f"    DL {half} ;"] if half > 0 else []
        lines += [f"    {_RF_VOLTAGE}(1,1) := {v_kv} ;",
                  f"    RF {_RF_VOLTAGE} 0 {f_hz} {phi} {self.r} ;{self._comment(elem)}"]
        if half > 0:
            lines.append(f"    DL {half} ;")
        return lines

    def _wiggler(self, elem):
        if elem["length"] <= 0:
            return []
        bw = elem.get("wiggler_field")
        period = elem.get("wiggler_period")
        if bw is None or period is None:
            return self._drift(elem)
        if period <= 0:
            raise ValueError(
                f"Wiggler '{elem.get('name', '')}': period must be positive, "
                f"got {period}"
            )
        # WI B0 lambda L d k I A: COSY takes the period lambda, not 2 pi / lambda.
        return [f"    WI {bw} {period} {elem['length']} {self.r} 0 0 0 ;"
                f"{self._comment(elem)}"]


def _dipole_fox(elem, d_half, comment="", emit_fc=True):
    """Generate COSYScript lines for a consolidated dipole of half-gap ``d_half``."""
    lines = []
    angle = elem["angle"]
    e1 = elem["entrance_angle"]
    e2 = elem["exit_angle"]
    h1 = elem.get("entrance_curvature", 0.0) or 0
    h2 = elem.get("exit_curvature", 0.0) or 0
    tilt = elem.get("tilt", 0.0) or 0.0
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

    # A tilted bend is a bend in a frame rotated by tilt_ref (x toward y)
    # about the entrance axis; RA rotates COSY's coordinates by the same angle.
    if tilt:
        lines.append(f"    RA {math.degrees(tilt)} ;")

    # CB wrap for negative bending angle
    use_cb = angle < 0
    if use_cb:
        lines.append("    CB ;")

    lines.append(f"    DIL {elem['length']} {abs(angle)} {d_half} {e1} {h1} {e2} {h2} ;{comment}")

    if use_cb:
        lines.append("    CB ;")

    if tilt:
        lines.append(f"    RA {-math.degrees(tilt)} ;")

    if has_fc:
        lines.append("    FD ;")

    return lines


def _format_enge(coeffs):
    """Pad/truncate Enge coefficients to exactly 6 values."""
    if len(coeffs) > 6:
        warnings.warn(
            f"Enge coefficients truncated from {len(coeffs)} to 6 "
            f"(COSY INFINITY FC accepts max 6)",
            PALSConversionWarning,
        )
    padded = (list(coeffs) + [0.0] * 6)[:6]
    return " ".join(str(float(c)) for c in padded)


# ---------------------------------------------------------------------------
# COSYScript template
# ---------------------------------------------------------------------------

def _particle_commands(particle, ke):
    """COSYScript commands that set the reference particle."""
    try:
        species = parse_species(particle)
    except ValueError as exc:
        raise ValueError(f"Unknown particle type '{particle}': {exc}") from None
    if species.charge == 0:
        raise ValueError(
            f"Particle type '{particle}' is neutral; COSY INFINITY needs a charged "
            f"reference particle")
    commands = [species.cosy_command(ke)]
    if species.anomalous_moment:
        # RP leaves COSY's spin G0 alone; RPS with spin tracking off sets it, as
        # the presets RPP, RPE and RPMU do.
        commands.append(f"RPS 0 {species.anomalous_moment}")
    return commands


def _fox_template(source_file, dimensions, order, ke, fr, lattice_body,
                  particle="electron", twiss=False, declarations=()):
    dim = dimensions
    reference = "".join(f"    {cmd} ;\n" for cmd in _particle_commands(particle, ke))
    safe_source = re.sub(r'[^\x20-\x7e]|[{}]', '', source_file)
    decl = "".join(f"        {line}\n" for line in declarations)
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
{{ Generated by pals2cosy {__version__} (PALS {__pals_commit__}) from {safe_source} }}
INCLUDE 'COSY' ;
PROCEDURE RUN ;
{twiss_vars}\
    PROCEDURE LATTICE ;
{decl}\
        UM ; CR ;
{lattice_body}
    ENDPROCEDURE ;

    OV {order} {dim} 0 ;
{reference}\
    FR {fr} ;
    LATTICE ;
{twiss_block}
ENDPROCEDURE ;
RUN ;
END ;
"""
