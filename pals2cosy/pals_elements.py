"""Translate expanded PALS lattice elements into converter element dicts.

Each PALS element kind is mapped to the internal representation described
in ``_elements``. Field strengths are resolved here: normalized strengths
(``KnN``, ``Ksol``, ``Bn0_ref`` vs ``g_ref``) are turned into fields with the
signed reference rigidity ``Brho = P0 / q``, integrated strengths are divided
by the element length, and multipole tilts are folded into normal and skew
components.

Author: Eremey Valetov
"""

import cmath
import math
import re
import warnings

from ._elements import MAX_MULTIPOLE_ORDER, blank_element
from .constants import C_LIGHT
from .diagnostics import (
    PALSConversionWarning, PALSDeprecationWarning, UnsupportedFeatureError,
)
from .expressions import ExpressionError

_KNOWN_GROUPS = frozenset({
    "ACKickerP", "ApertureP", "BeamBeamP", "BendP", "BodyShiftP", "ConverterP",
    "CoordinateSetP", "ElectricMultipoleP", "FloorP", "FoilP", "ForkFromP",
    "ForkP", "GirderP", "MagneticMultipoleP", "MetaP", "ParticleP", "PatchP",
    "ReferenceP", "ReferenceChangeP", "RFP", "ShiftedFloorP", "SolenoidP",
    "TaylorP", "TrackingP", "TwissP",
})
_NON_GROUP_KEYS = frozenset({
    "kind", "name", "length", "is_on", "s_position", "inherit", "elements",
})
# pals2cosy extensions and keys accepted from pre-0.6 files.
_EXTENSION_KEYS = frozenset({"WigglerP", "RFCavityP", "plane", "instrument_type"})

_MULTIPOLE_KEY = re.compile(r"^(?P<form>[BK])(?P<comp>[ns])(?P<order>\d+)"
                            r"(?P<integrated>L)?(?P<taper>_taper)?$")
_ELECTRIC_KEY = re.compile(r"^E[ns]\d+L?$")
_TILT_KEY = re.compile(r"^tilt(?P<order>\d+)$")

# Components of the parameter groups pals2cosy reads, at the pinned
# specification commit (parameters/*.md), plus the pals2cosy extension groups.
_OFFSETS = ("x_offset", "y_offset", "z_offset", "x_rot", "y_rot", "z_rot")
_GROUP_KEYS = {
    "ApertureP": frozenset({
        "x_min", "x_max", "x_center", "x_width", "y_min", "y_max", "y_center",
        "y_width", "shape", "location", "vertices", "material", "thickness",
        "aperture_shifts_with_body", "aperture_active"}),
    "BendP": frozenset({
        "angle_ref", "Bn0_ref", "e1", "e2", "e1_rect", "e2_rect", "edge1_int",
        "edge2_int", "g_ref", "h1", "h2", "Kn0_from_g_ref", "L_chord", "L_sagitta",
        "L_rectangle", "multipole_geometry", "ref_geometry", "radius_ref", "tilt_ref"}),
    "BodyShiftP": frozenset(_OFFSETS),
    "ForkP": frozenset({
        "to_line", "destination_element", "new_branch", "propagate_reference",
        "forked_to"}),
    "PatchP": frozenset(_OFFSETS + ("flexible", "ref_coords", "user_sets_length")),
    "ReferenceP": frozenset({"species_ref", "pc_ref", "E_tot_ref", "time_ref"}),
    "ReferenceChangeP": frozenset({
        "dtime_ref", "dE_ref", "dpc_ref", "time_ref", "E_tot_ref", "pc_ref",
        "species_ref"}),
    "RFP": frozenset({
        "frequency", "harmonic_number", "voltage", "gradient", "phase",
        "multipass_phase", "cavity_type", "num_cells", "zero_phase", "L_active",
        "dE_ref"}),
    "SolenoidP": frozenset({"Ksol", "Bsol", "Bz"}),
    "WigglerP": frozenset({"peak_field", "period"}),
    "RFCavityP": frozenset({"voltage_kv", "frequency_hz", "phase_deg"}),
}
_GROUP_KEY_PATTERNS = {
    "MagneticMultipoleP": (_MULTIPOLE_KEY, _TILT_KEY),
    "ElectricMultipoleP": (_ELECTRIC_KEY, _TILT_KEY),
}

_ZERO_LENGTH_KINDS = ("Marker", "BeginningEle", "Placeholder", "Fiducial", "Fork",
                      "FloorShift")
_MULTIPOLE_KINDS = frozenset({
    "Quadrupole", "Sextupole", "Octupole", "Multipole", "Kicker", "Bend",
})
_UNMODELLED_KINDS = frozenset({
    "ACKicker", "BeamBeam", "Converter", "CrabCavity", "EGun", "Feedback",
    "Foil", "Girder", "Mask", "Match", "Taylor", "UnionEle",
})

_RF_ZERO_PHASE_DEG = {
    "ACCELERATING": 90.0,
    "BELOW_TRANSITION": 0.0,
    "ABOVE_TRANSITION": 180.0,
}

_REF_GEOMETRIES = ("ARC", "CHORD", "ENTRANCE_COORDS", "EXIT_COORDS")


def _warn(message):
    warnings.warn(message, PALSConversionWarning, stacklevel=3)


def is_extension_block(value):
    """``<name>: {extension: <type>, ...}`` blocks may appear at any level."""
    return isinstance(value, dict) and "extension" in value


class ExtensionLabels:
    """The ``extension_labels`` block of a PALS document."""

    def __init__(self, block=None):
        block = block or {}
        if not isinstance(block, dict):
            raise ValueError("PALS 'extension_labels' must be a mapping")
        self.names = set(block.get("names") or {})
        self.prefixes = tuple(block.get("prefixes") or {})
        self.suffixes = tuple(block.get("suffixes") or {})

    def matches(self, key):
        return (key in self.names
                or any(key.startswith(p) for p in self.prefixes)
                or any(key.endswith(s) for s in self.suffixes))


class ElementContext:
    """State shared while translating the elements of one branch."""

    def __init__(self, scope, catalog, labels, species, e_tot_ev):
        self.scope = scope
        self.catalog = catalog
        self.labels = labels
        self.species = species
        self.e_tot_ev = e_tot_ev

    @property
    def pc_ev(self):
        m = self.species.mass_mev * 1e6
        return math.sqrt(max(self.e_tot_ev ** 2 - m ** 2, 0.0))

    @property
    def brho(self):
        """Signed magnetic rigidity P0/q [T m] at the current element."""
        return self.pc_ev / (C_LIGHT * self.species.charge)

    def number(self, value, name, field, allow_infinite=False):
        try:
            result = self.scope.evaluate(value)
        except ExpressionError as exc:
            raise ValueError(
                f"Element '{name}': field '{field}' must be numeric, "
                f"got {value!r} ({exc})"
            ) from exc
        if isinstance(result, str):
            raise ValueError(
                f"Element '{name}': field '{field}' must be numeric, got {value!r}")
        if math.isnan(result) or (math.isinf(result) and not allow_infinite):
            raise ValueError(
                f"Element '{name}': field '{field}' must be a finite number, got {result}")
        return result

    def optional(self, group, key, name, gname, default=None):
        value = group.get(key) if group else None
        if value is None:
            return default
        return self.number(value, name, f"{gname}.{key}")

    def group(self, raw, gname, name):
        """Return parameter group ``gname`` of ``raw`` with ``inherit`` resolved."""
        group = raw.get(gname)
        if group is None:
            return None
        return self._resolve_group(group, gname, name, ())

    def _resolve_group(self, group, gname, name, seen):
        if not isinstance(group, dict):
            raise ValueError(f"Element '{name}': '{gname}' must be a mapping")
        target = group.get("inherit")
        if target is None:
            return group
        if target in seen:
            raise ValueError(f"Element '{name}': circular {gname} inheritance via '{target}'")
        if ">" in target:
            elem_name, _, sub = target.partition(">")
            base = (self.catalog.get(elem_name) or {}).get(sub)
        else:
            base = self.catalog.get(target)
            if base is not None and base.get("kind") not in (None, gname):
                base = None
        if not isinstance(base, dict):
            raise ValueError(
                f"Element '{name}': cannot resolve {gname} inherit '{target}'")
        base = self._resolve_group(base, gname, name, seen + (target,))
        merged = {k: v for k, v in base.items() if k not in ("kind", "inherit")}
        merged.update((k, v) for k, v in group.items() if k != "inherit")
        return merged


# ---------------------------------------------------------------------------
# Bend geometry
# ---------------------------------------------------------------------------

def resolve_bend_geometry(name, length=None, g=None, chord=None, rect=None,
                          angle=None):
    """Resolve the BendP "two of three" parameter sets.

    The sets are {curvature g}, {arc length, L_chord, L_rectangle} and
    {angle}. Two parameters from different sets fix the geometry; a third,
    redundant parameter is accepted if consistent. Returns
    ``(arc_length, g, angle)``.
    """
    given_lengths = [k for k, v in (("length", length), ("L_chord", chord),
                                    ("L_rectangle", rect)) if v is not None]
    if len(given_lengths) > 1:
        raise ValueError(
            f"Bend '{name}': only one of length, L_chord, L_rectangle may be "
            f"given, got {', '.join(given_lengths)}")

    def from_length(L):
        if g is not None:
            return L, g, g * L
        if angle is not None:
            if L > 0:
                return L, angle / L, angle
            if angle != 0:
                raise ValueError(f"Bend '{name}': nonzero angle_ref with zero length")
        return L, 0.0, 0.0

    if length is not None:
        sol = from_length(length)
    elif chord is not None:
        if g is not None and g != 0:
            x = g * chord / 2
            if abs(x) > 1:
                raise ValueError(f"Bend '{name}': L_chord exceeds the bend diameter")
            ang = 2 * math.asin(x)
            sol = (ang / g, g, ang)
        elif angle is not None and angle != 0:
            if chord <= 0:
                raise ValueError(f"Bend '{name}': nonzero angle_ref with zero L_chord")
            gg = 2 * math.sin(angle / 2) / chord
            sol = (angle / gg, gg, angle)
        else:
            sol = (chord, 0.0, 0.0)
    elif rect is not None:
        if g is not None and g != 0:
            x = g * rect
            if abs(x) > 1:
                raise ValueError(f"Bend '{name}': L_rectangle exceeds the bend radius")
            ang = math.asin(x)
            sol = (ang / g, g, ang)
        elif angle is not None and angle != 0:
            if rect <= 0 or math.sin(angle) == 0:
                raise ValueError(
                    f"Bend '{name}': angle_ref and L_rectangle do not fix the curvature")
            gg = math.sin(angle) / rect
            sol = (angle / gg, gg, angle)
        else:
            sol = (rect, 0.0, 0.0)
    elif g is not None and angle is not None:
        if g != 0:
            sol = (angle / g, g, angle)
        elif angle != 0:
            raise ValueError(f"Bend '{name}': nonzero angle_ref with zero curvature")
        else:
            sol = (0.0, 0.0, 0.0)
    elif angle is not None and angle != 0:
        raise ValueError(
            f"Bend '{name}': angle_ref alone does not fix the bend; give a length "
            f"or a curvature as well")
    else:
        sol = (0.0, 0.0, 0.0)

    L, gg, ang = sol
    if L < 0:
        raise ValueError(f"Bend '{name}': resolved arc length {L} is negative")
    computed = {
        "length": L,
        "g": gg,
        "angle": ang,
        "L_chord": 2 * math.sin(ang / 2) / gg if gg else L,
        "L_rectangle": math.sin(ang) / gg if gg else L,
    }
    for key, value in (("length", length), ("g", g), ("angle", angle),
                       ("L_chord", chord), ("L_rectangle", rect)):
        if value is not None and not math.isclose(value, computed[key],
                                                  rel_tol=1e-9, abs_tol=1e-12):
            raise ValueError(
                f"Bend '{name}': BendP is over-specified and inconsistent "
                f"({key} = {value}, but the other parameters give {computed[key]})")
    return L, gg, ang


def resolve_edge_angles(name, angle, ref_geometry, e1, e2, e1_rect, e2_rect):
    """Sector-referenced face angles (e1, e2) from e1/e1_rect and e2/e2_rect."""
    if ref_geometry not in _REF_GEOMETRIES:
        raise ValueError(
            f"Bend '{name}': unknown ref_geometry {ref_geometry!r}; "
            f"expected one of {', '.join(_REF_GEOMETRIES)}")
    offsets = {
        "ARC": (angle / 2, angle / 2),
        "CHORD": (angle / 2, angle / 2),
        "ENTRANCE_COORDS": (0.0, angle),
        "EXIT_COORDS": (angle, 0.0),
    }[ref_geometry]
    result = []
    for label, e, e_rect, offset in (("1", e1, e1_rect, offsets[0]),
                                     ("2", e2, e2_rect, offsets[1])):
        if e is not None and e_rect is not None:
            if not math.isclose(e, e_rect + offset, rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError(
                    f"Bend '{name}': e{label} = {e} and e{label}_rect = {e_rect} "
                    f"are inconsistent for ref_geometry {ref_geometry}")
        elif e is None:
            if e_rect is None:
                e_rect = None if ref_geometry == "ARC" else 0.0
            e = 0.0 if e_rect is None else e_rect + offset
        result.append(e)
    return tuple(result)


# ---------------------------------------------------------------------------
# Multipoles
# ---------------------------------------------------------------------------

def read_multipoles(ctx, raw, name, length, gname="MagneticMultipoleP"):
    """Read a MagneticMultipoleP group.

    Returns ``(fields, thin, normalized)`` where ``fields`` maps each order N
    to the complex coefficient ``Bn + i Bs`` [T/m^N] with the order's tilt
    folded in. For a zero-length element, ``thin`` is True and the fields are
    length-integrated [T/m^(N-1)]. ``None`` is returned when the group is
    absent.
    """
    group = ctx.group(raw, gname, name)
    if group is None:
        return None
    components = {}
    tilts = {}
    for key, value in group.items():
        if value is None:
            continue
        m = _MULTIPOLE_KEY.match(key)
        if m:
            order = int(m.group("order"))
            form = (m.group("form"), bool(m.group("integrated")))
            slot = components.setdefault((order, m.group("comp")), {})
            slot[form] = slot.get(form, 0.0) + ctx.number(value, name, f"{gname}.{key}")
            continue
        m = _TILT_KEY.match(key)
        if m:
            tilts[int(m.group("order"))] = ctx.number(value, name, f"{gname}.{key}")

    thin = length == 0
    normalized = False
    fields = {}
    forms_by_order = {}
    for (order, comp), slot in components.items():
        if len(slot) > 1:
            spelled = sorted(f"{f}{comp}{order}{'L' if i else ''}" for f, i in slot)
            raise ValueError(
                f"Element '{name}': {gname} gives the same component in several "
                f"forms ({', '.join(spelled)})")
        (form, integrated), value = next(iter(slot.items()))
        other = forms_by_order.setdefault(order, (form, integrated))
        if other != (form, integrated):
            raise ValueError(
                f"Element '{name}': the normal and skew components of order {order} "
                f"must both be fields or both normalized, and both integrated or not")
        if form == "K":
            value *= _brho(ctx, name)
            normalized = True
        if integrated and not thin:
            value /= length
        elif not integrated and thin:
            if value != 0:
                _warn(f"Element '{name}': non-integrated {gname} strength on a "
                      f"zero-length element has no effect")
            value = 0.0
        fields[order] = fields.get(order, 0j) + (value if comp == "n" else 1j * value)

    for order, tilt in tilts.items():
        if order in fields and tilt:
            fields[order] *= cmath.exp(-1j * (order + 1) * tilt)
    return fields, thin and any(fields.values()), normalized


def _brho(ctx, name):
    if ctx.species.charge == 0:
        raise ValueError(
            f"Element '{name}': normalized strength needs a charged reference "
            f"particle, but '{ctx.species.name}' is neutral")
    return ctx.brho


def _store_multipoles(elem, fields, name, kind):
    for order in sorted(fields):
        c = fields[order]
        if c == 0:
            continue
        if order == 0:
            _warn(f"{kind} '{name}': dipole (order 0) multipole component is not "
                  f"represented in COSYScript")
            continue
        if order > MAX_MULTIPOLE_ORDER:
            _warn(f"{kind} '{name}': multipole order {order} exceeds the supported "
                  f"maximum {MAX_MULTIPOLE_ORDER}; dropped")
            continue
        tol = 1e-14 * abs(c)
        if abs(c.real) > tol:
            elem["bn"][order] = c.real
        if abs(c.imag) > tol:
            elem["bs"][order] = c.imag


def _no_strength(kind, name, what):
    _warn(f"{kind} '{name}' has no readable strength ({what}); converted to a drift")


# ---------------------------------------------------------------------------
# Element kinds
# ---------------------------------------------------------------------------

def _multipole(ctx, raw, name, kind, length):
    etype = {"Quadrupole": "QPF", "Sextupole": "SXT", "Octupole": "OCT",
             "Multipole": "MULT"}[kind]
    elem = blank_element(etype, name, length)
    result = read_multipoles(ctx, raw, name, length)
    if result is None:
        _no_strength(kind, name, "no MagneticMultipoleP")
        return elem
    fields, thin, normalized = result
    _store_multipoles(elem, fields, name, kind)
    elem["thin"] = thin
    elem["normalized"] = normalized
    if kind == "Quadrupole" and elem["bn"].get(1, 0.0) < 0:
        elem["type"] = "QPD"
    return elem


def _kicker(ctx, raw, name, kind, length):
    elem = blank_element("KICKER", name, length)
    result = read_multipoles(ctx, raw, name, length)
    if result is None:
        _no_strength(kind, name, "no MagneticMultipoleP")
    elif any(result[0].values()):
        _warn(f"Kicker '{name}': the kick is not represented in COSYScript; "
              f"converted to a drift")
    return elem


def _bend(ctx, raw, name, kind, length):
    bend = ctx.group(raw, "BendP", name) or {}

    def num(key, default=None):
        return ctx.optional(bend, key, name, "BendP", default)

    curvatures = {}
    if bend.get("g_ref") is not None:
        curvatures["g_ref"] = num("g_ref")
    if bend.get("Bn0_ref") is not None:
        curvatures["Bn0_ref"] = num("Bn0_ref") / _brho(ctx, name)
    if bend.get("radius_ref") is not None:
        radius = ctx.number(bend["radius_ref"], name, "BendP.radius_ref",
                            allow_infinite=True)
        if radius == 0:
            raise ValueError(f"Bend '{name}': radius_ref must be nonzero")
        curvatures["radius_ref"] = 0.0 if math.isinf(radius) else 1 / radius
    g = None
    for value in curvatures.values():
        if g is not None and not math.isclose(g, value, rel_tol=1e-9, abs_tol=1e-15):
            raise ValueError(
                f"Bend '{name}': {', '.join(curvatures)} give different curvatures")
        g = value

    L, g, angle = resolve_bend_geometry(
        name, length=length, g=g, chord=num("L_chord"), rect=num("L_rectangle"),
        angle=num("angle_ref"))
    e1, e2 = resolve_edge_angles(
        name, angle, bend.get("ref_geometry", "ARC"),
        num("e1"), num("e2"), num("e1_rect"), num("e2_rect"))

    for key in ("edge1_int", "edge2_int"):
        if num(key, 0.0):
            _warn(f"Bend '{name}': BendP.{key} is not represented; COSY applies its "
                  f"own Enge fringe-field model with --fr 1, 2 or 3")

    mult = read_multipoles(ctx, raw, name, L)
    if mult is not None:
        fields, _, _ = mult
        b0 = fields.pop(0, 0j)
        reference = g * _brho(ctx, name) if g else 0.0
        if b0.imag:
            _warn(f"Bend '{name}': skew dipole field Bs0/Ks0 is not represented")
        if b0.real and not math.isclose(b0.real, reference, rel_tol=1e-9, abs_tol=1e-12):
            _warn(f"Bend '{name}': actual dipole field {b0.real} T differs from the "
                  f"reference field {reference} T; COSY's DIL uses the reference field")
        if any(fields.values()):
            _warn(f"Bend '{name}': multipoles in a bend (combined-function magnets) "
                  f"are not converted")
    elif bend.get("Kn0_from_g_ref") is False and g:
        _warn(f"Bend '{name}': Kn0_from_g_ref is false with no Kn0/Bn0, so the bend "
              f"has no field; COSY's DIL uses the reference field")

    elem = blank_element("DPH", name, L)
    direction = -1.0 if angle < 0 else 1.0
    elem["angle"] = math.degrees(angle)
    elem["dipole_length"] = L
    # COSY bends negative angles through CB, which mirrors x; PALS face angles
    # of a mirrored (negative-angle) bend carry the sign of the bend angle.
    elem["entrance_edge_angle"] = math.degrees(direction * e1)
    elem["exit_edge_angle"] = math.degrees(direction * e2)
    elem["entrance_edge_curvature"] = num("h1", 0.0)
    elem["exit_edge_curvature"] = num("h2", 0.0)
    elem["tilt"] = num("tilt_ref", 0.0)
    elem["pole_gap"] = _pole_gap(ctx, raw, name)
    return elem


def _pole_gap(ctx, raw, name):
    aperture = ctx.group(raw, "ApertureP", name) or {}
    width = ctx.optional(aperture, "y_width", name, "ApertureP")
    if width is None:
        y_min = ctx.optional(aperture, "y_min", name, "ApertureP")
        y_max = ctx.optional(aperture, "y_max", name, "ApertureP")
        width = y_max - y_min if y_min is not None and y_max is not None else 0.0
    if width < 0:
        raise ValueError(
            f"Element '{name}': field 'ApertureP.y_width' must be non-negative, got {width}")
    return width


def _solenoid(ctx, raw, name, kind, length):
    elem = blank_element("SOL", name, length)
    sol = ctx.group(raw, "SolenoidP", name)
    bz = None
    if sol is not None:
        given = [k for k in ("Ksol", "Bsol", "Bz") if sol.get(k) is not None]
        if len(given) > 1:
            raise ValueError(f"Solenoid '{name}': give only one of {', '.join(given)}")
        if "Bsol" in given:
            bz = ctx.number(sol["Bsol"], name, "SolenoidP.Bsol")
        elif "Ksol" in given:
            bz = ctx.number(sol["Ksol"], name, "SolenoidP.Ksol") * _brho(ctx, name)
            elem["normalized"] = True
        elif "Bz" in given:
            warnings.warn(
                f"Solenoid '{name}': SolenoidP.Bz is deprecated; use SolenoidP.Bsol",
                PALSDeprecationWarning, stacklevel=2)
            bz = ctx.number(sol["Bz"], name, "SolenoidP.Bz")
    if bz is None:
        _no_strength(kind, name, "no SolenoidP.Bsol or Ksol")
    elif length == 0 and bz:
        _warn(f"Solenoid '{name}': zero-length solenoid has no effect")
    elem["bz"] = bz
    return elem


def _warn_unconverted_multipoles(ctx, raw, name, kind, length):
    result = read_multipoles(ctx, raw, name, length)
    if result is not None and any(result[0].values()):
        _warn(f"{kind} '{name}': MagneticMultipoleP fields are not represented in COSYScript")


def _warn_unconverted_solenoid(ctx, raw, name, kind):
    sol = ctx.group(raw, "SolenoidP", name) or {}
    if any(ctx.optional(sol, k, name, "SolenoidP", 0.0) for k in ("Ksol", "Bsol")):
        _warn(f"{kind} '{name}': SolenoidP fields are not represented in COSYScript")


def _rf_cavity(ctx, raw, name, kind, length):
    elem = blank_element("RFC", name, length)
    _warn_unconverted_solenoid(ctx, raw, name, kind)
    rfp = ctx.group(raw, "RFP", name)
    if rfp is None and raw.get("RFCavityP") is not None:
        return _legacy_rf(ctx, raw, name, elem)
    if rfp is None:
        _no_strength(kind, name, "no RFP")
        return elem

    def num(key, default=None):
        return ctx.optional(rfp, key, name, "RFP", default)

    l_active = num("L_active", length)
    voltage, gradient = num("voltage"), num("gradient")
    if voltage is not None and gradient is not None:
        if not math.isclose(voltage, gradient * l_active, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(
                f"RFCavity '{name}': voltage {voltage} V and gradient {gradient} V/m "
                f"disagree for L_active = {l_active} m")
    elif voltage is None and gradient is not None:
        voltage = gradient * l_active
    if voltage is None:
        _no_strength(kind, name, "no RFP.voltage or gradient")

    d_e_ref = num("dE_ref", 0.0)
    if not voltage:
        _change_reference_energy(ctx, name, d_e_ref, 0.0)
        return elem

    frequency, harmonic = num("frequency"), num("harmonic_number")
    if frequency is not None and harmonic is not None:
        raise ValueError(
            f"RFCavity '{name}': give either RFP.frequency or harmonic_number, not both")
    if frequency is None and harmonic is None:
        raise ValueError(
            f"RFCavity '{name}': RFP needs frequency or harmonic_number")
    if frequency is not None and frequency <= 0:
        raise ValueError(f"RFCavity '{name}': RFP.frequency must be positive")

    zero_phase = rfp.get("zero_phase", "ACCELERATING")
    if zero_phase not in _RF_ZERO_PHASE_DEG:
        raise ValueError(
            f"RFCavity '{name}': unknown zero_phase {zero_phase!r}; expected one of "
            f"{', '.join(_RF_ZERO_PHASE_DEG)}")
    if num("multipass_phase", 0.0):
        _warn(f"RFCavity '{name}': RFP.multipass_phase is ignored")
    if rfp.get("cavity_type", "STANDING_WAVE") != "STANDING_WAVE":
        _warn(f"RFCavity '{name}': cavity_type {rfp['cavity_type']} is modelled by "
              f"COSY's thin standing-wave RF kick")

    phase_deg = cosy_rf_phase(num("phase", 0.0), zero_phase, ctx.species.charge)
    cosy_gain = ctx.species.charge * voltage * math.sin(math.radians(phase_deg))
    _change_reference_energy(ctx, name, d_e_ref, cosy_gain)

    elem["rf_voltage_kv"] = voltage / 1e3
    elem["rf_frequency_hz"] = frequency
    elem["rf_harmonic"] = harmonic
    elem["rf_phase_deg"] = phase_deg
    elem["rf_de_ref_ev"] = d_e_ref
    return elem


def _change_reference_energy(ctx, name, d_e_ref, cosy_gain):
    """Apply ``RFP.dE_ref`` to the PALS reference; report a COSY mismatch."""
    if not math.isclose(d_e_ref, cosy_gain, rel_tol=1e-6, abs_tol=1e-3):
        _warn(f"RFCavity '{name}': COSY re-references the energy to the reference "
              f"particle's gain of {cosy_gain} eV, while RFP.dE_ref = {d_e_ref} eV; "
              f"downstream maps are relative to the COSY reference")
    ctx.e_tot_ev += d_e_ref


def cosy_rf_phase(phase_rad, zero_phase, charge):
    """COSY ``RF`` phase [deg] for a PALS phase and ``zero_phase`` convention.

    COSY gives the reference particle ``Z * V * sin(PHI)``; PALS measures the
    phase from the crest (``ACCELERATING``) or from the stable zero crossing
    below or above transition, independent of the sign of the charge.
    """
    deg = math.degrees(phase_rad) + _RF_ZERO_PHASE_DEG[zero_phase]
    if charge < 0:
        deg += 180.0
    deg = math.fmod(deg, 360.0)
    if deg > 180.0:
        deg -= 360.0
    elif deg <= -180.0:
        deg += 360.0
    return deg


def _legacy_rf(ctx, raw, name, elem):
    warnings.warn(
        f"RFCavity '{name}': RFCavityP is deprecated; use the PALS RFP group "
        f"(voltage [V], frequency [Hz], phase [rad])",
        PALSDeprecationWarning, stacklevel=3)
    rf = raw["RFCavityP"]
    if not isinstance(rf, dict):
        raise ValueError(f"Element '{name}': 'RFCavityP' must be a mapping")
    v = ctx.optional(rf, "voltage_kv", name, "RFCavityP")
    f = ctx.optional(rf, "frequency_hz", name, "RFCavityP")
    phi = ctx.optional(rf, "phase_deg", name, "RFCavityP")
    if f is not None and f <= 0:
        raise ValueError(
            f"Element '{name}': field 'RFCavityP.frequency_hz' must be positive, got {f}")
    elem["rf_voltage_kv"] = v
    elem["rf_frequency_hz"] = f
    elem["rf_phase_deg"] = phi
    return elem


def _wiggler(ctx, raw, name, kind, length):
    elem = blank_element("UND", name, length)
    wig = ctx.group(raw, "WigglerP", name)
    if wig is None:
        _no_strength(kind, name, "no WigglerP")
        return elem
    elem["wiggler_field"] = ctx.optional(wig, "peak_field", name, "WigglerP")
    period = ctx.optional(wig, "period", name, "WigglerP")
    if period is not None and period <= 0:
        raise ValueError(
            f"Element '{name}': field 'WigglerP.period' must be positive, got {period}")
    elem["wiggler_period"] = period
    return elem


def _drift(ctx, raw, name, kind, length):
    return blank_element("DRIFT", name, length)


def _instrument(ctx, raw, name, kind, length):
    return blank_element("INSTRUMENT", name, length)


def _zero_length(ctx, raw, name, kind, length):
    if length != 0:
        raise ValueError(f"{kind} '{name}' must have zero length, got {length}")
    if kind == "Fork":
        fork = raw.get("ForkP") or {}
        _warn(f"Fork '{name}' to '{fork.get('to_line')}' is not followed; only the "
              f"selected branch is converted")
    return blank_element("MARKER", name, 0.0)


def _patch(ctx, raw, name, kind, length):
    patch = ctx.group(raw, "PatchP", name) or {}
    shifts = [k for k in _OFFSETS if ctx.optional(patch, k, name, "PatchP", 0.0)]
    if shifts or length or patch.get("flexible"):
        raise UnsupportedFeatureError(
            f"Patch '{name}' shifts the reference orbit; Patch elements are not supported")
    return blank_element("MARKER", name, 0.0)


def _reference_change(ctx, raw, name, kind, length):
    rc = ctx.group(raw, "ReferenceChangeP", name) or {}
    changed = [k for k in ("dE_ref", "dpc_ref", "E_tot_ref", "pc_ref", "species_ref")
               if rc.get(k) is not None]
    if changed:
        raise UnsupportedFeatureError(
            f"ReferenceChange '{name}' changes {', '.join(changed)}; changing the "
            f"reference energy or species mid-branch is not supported")
    timing = [k for k in ("dtime_ref", "time_ref") if rc.get(k) is not None]
    if len(timing) > 1:
        raise ValueError(
            f"ReferenceChange '{name}': give only one of dtime_ref and time_ref")
    # Field-free, so a length is kept as a drift. Time adjustments enter the
    # branch reference time that RFP.harmonic_number refers to.
    elem = blank_element("MARKER", name, length)
    for key in timing:
        elem[key] = ctx.number(rc[key], name, f"ReferenceChangeP.{key}")
    return elem


def _unmodelled(ctx, raw, name, kind, length):
    what = (f"extension kind '{kind}'" if ctx.labels.matches(kind)
            else f"{kind} element")
    _warn(f"{what} '{name}' has no COSYScript mapping; "
          + (f"converted to a drift of length {length}" if length else "skipped"))
    return blank_element("DRIFT", name, length)


_HANDLERS = {
    "Drift": _drift,
    "Quadrupole": _multipole,
    "Sextupole": _multipole,
    "Octupole": _multipole,
    "Multipole": _multipole,
    "Kicker": _kicker,
    "Bend": _bend,
    "Solenoid": _solenoid,
    "RFCavity": _rf_cavity,
    "Wiggler": _wiggler,
    "Instrument": _instrument,
    "Patch": _patch,
    "ReferenceChange": _reference_change,
}
_HANDLERS.update((k, _zero_length) for k in _ZERO_LENGTH_KINDS)
_HANDLERS.update((k, _unmodelled) for k in _UNMODELLED_KINDS)

_FIELD_KINDS = frozenset({
    "Quadrupole", "Sextupole", "Octupole", "Multipole", "Kicker", "Bend",
    "Solenoid", "RFCavity", "Wiggler",
})


def normalize_element(ctx, raw):
    """Translate one expanded PALS element into a converter element dict."""
    name = raw.get("name", "<unnamed>")
    kind = raw.get("kind")
    if kind is None:
        _warn(f"Element '{name}' has no 'kind'; defaulting to Drift")
        kind = "Drift"
    if kind in ("SBend", "RBend"):
        warnings.warn(
            f"Element '{name}': kind '{kind}' is deprecated; use kind 'Bend' "
            f"(treated as a Bend with ref_geometry ARC, as in pals2cosy 0.5)",
            PALSDeprecationWarning, stacklevel=2)
        kind = "Bend"

    raw_length = raw.get("length")
    length = 0.0 if raw_length is None else ctx.number(raw_length, name, "length")
    if length < 0:
        raise ValueError(f"Element '{name}': field 'length' must be non-negative, got {length}")
    length_given = None if raw_length is None else length

    _check_keys(ctx, raw, name)
    _check_group_keys(ctx, raw, name)
    _check_body_shift(ctx, raw, name)
    _check_electric(ctx, raw, name, kind)

    if raw.get("is_on") is False and kind in _FIELD_KINDS:
        if kind == "Bend":
            raise UnsupportedFeatureError(
                f"Bend '{name}' has is_on: false; a bend without field still curves "
                f"the reference orbit, which COSY cannot represent")
        return blank_element("DRIFT", name, length)

    handler = _HANDLERS.get(kind, _unmodelled)
    if kind == "Bend":
        return _bend(ctx, raw, name, kind, length_given)
    if kind not in _MULTIPOLE_KINDS and handler is not _unmodelled:
        _warn_unconverted_multipoles(ctx, raw, name, kind, length)
    return handler(ctx, raw, name, kind, length)


def _check_keys(ctx, raw, name):
    for key, value in raw.items():
        if (key in _KNOWN_GROUPS or key in _NON_GROUP_KEYS or key in _EXTENSION_KEYS
                or ctx.labels.matches(key) or is_extension_block(value)):
            continue
        _warn(f"Element '{name}': unknown key '{key}' ignored (declare it in "
              f"extension_labels if it is extension data)")


def _known_group_key(gname, key):
    if gname in _GROUP_KEY_PATTERNS:
        return any(p.match(key) for p in _GROUP_KEY_PATTERNS[gname])
    return key in _GROUP_KEYS[gname]


def _check_group_keys(ctx, raw, name):
    """Report misspelled or unknown components of the groups pals2cosy reads."""
    for gname in raw:
        if gname not in _GROUP_KEYS and gname not in _GROUP_KEY_PATTERNS:
            continue
        for key, value in (ctx.group(raw, gname, name) or {}).items():
            if not (_known_group_key(gname, key) or ctx.labels.matches(key)
                    or is_extension_block(value)):
                _warn(f"Element '{name}': unknown {gname} key '{key}' ignored")


def _check_body_shift(ctx, raw, name):
    shift = ctx.group(raw, "BodyShiftP", name)
    if not shift:
        return
    nonzero = [k for k in _OFFSETS
               if shift.get(k) is not None and ctx.number(shift[k], name, f"BodyShiftP.{k}")]
    if nonzero:
        _warn(f"Element '{name}': BodyShiftP misalignment ({', '.join(nonzero)}) "
              f"is not represented")


def _check_electric(ctx, raw, name, kind):
    group = ctx.group(raw, "ElectricMultipoleP", name)
    if not group:
        return
    values = [ctx.number(v, name, f"ElectricMultipoleP.{k}")
              for k, v in group.items() if v is not None and _ELECTRIC_KEY.match(k)]
    if any(values):
        _warn(f"{kind} '{name}': ElectricMultipoleP fields are not converted")
