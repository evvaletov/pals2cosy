"""Parse official PALS format (``PALS:`` root) into normalized element lists.

Reads documents written against the PALS specification
(github.com/pals-project/pals, see ``pals2cosy.__pals_commit__``): the
``facility`` list with element, ``BeamLine`` and ``Lattice`` definitions,
``use``, constants and variables, expressions, ``inherit`` and ``repeat``
(including negative repeats, which reverse the element order), and the
reference particle of the branch's ``BeginningEle``. Constructs that change
the lattice in ways a single COSY INFINITY transfer map cannot follow raise
``UnsupportedFeatureError``; those that leave the converted branch intact
(other branches, Forks, multipass bookkeeping, adjustments, RELATIVE
controllers, girders) are reported with ``PALSConversionWarning``.

Produces the same (beam_params, elements) tuple as ``parser.parse_lattice()``.

Author: Eremey Valetov
"""

import math
import os
import re
import warnings
from dataclasses import dataclass

from ._elements import blank_element
from ._io import load_file as _load_file
from ._version import __pals_commit__
from .constants import C_LIGHT, F_RF_DEFAULT
from .diagnostics import (
    PALSConversionWarning, PALSDeprecationWarning, PALSReminder,
    UnsupportedFeatureError,
)
from .expressions import ExpressionError, Scope
from .pals_elements import (
    ElementContext, ExtensionLabels, is_extension_block, normalize_element,
)
from .species import parse_species

# Major versions of the PALS schema this release understands. The standard
# has not fixed a versioning scheme yet; files carry null, 0.x or 1.x.
KNOWN_MAJOR_VERSIONS = (0, 1)

_ROOT_KEYS = frozenset({
    "adjustments", "authors", "extension_labels", "facility", "girders",
    "include", "load", "notes", "phase_space_coordinates", "post_expansion",
    "reminders", "version",
})


# Kinds that never appear in the tracked element sequence: the branch's
# BeginningEle only carries initial values, and expansion removes Placeholders.
_BOOKKEEPING_KINDS = ("BeginningEle", "Placeholder")

# Expansion brackets the items of a BeamLine that declares a `length` with
# these marks, so that the remainder can be filled with a drift once the
# element lengths are known.
_LINE_START = "<BeamLine start>"
_LINE_END = "<BeamLine end>"

_RANDOM_CALL = re.compile(r"\brandom(?:_gauss)?\s*\(")


def _warn(message):
    warnings.warn(message, PALSConversionWarning, stacklevel=3)


def parse_lattice(path, beamline_name=None, ke_override=None, particle_override=None):
    """Parse an official PALS file (``PALS:`` root key).

    Parameters
    ----------
    path : str
        Path to ``.pals.yaml`` or ``.pals.json`` file.
    beamline_name : str or None
        Branch of the selected Lattice, or any BeamLine, to expand. Default:
        the first branch of the Lattice named by the last ``use`` (else the
        last Lattice; else the last BeamLine).
    ke_override : float or None
        Kinetic energy in MeV, overriding the file's reference energy.
    particle_override : str or None
        Reference species, overriding ``ReferenceP.species_ref``.

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

    labels = ExtensionLabels(pals.get("extension_labels"))
    version = _check_root(pals, labels)
    facility = _Facility(pals.get("facility") or [], labels)
    branch, root, periodic = facility.select_branch(beamline_name)
    flat = _expand_beamline(root, facility.catalog)

    ctx = ElementContext(facility.scope, facility.catalog, labels, None, 0.0)
    species, kinetic_energy, time_ref = _reference(flat, facility, ctx, ke_override,
                                                   particle_override)
    ctx.species = species
    ctx.e_tot_ev = (kinetic_energy + species.mass_mev) * 1e6
    brho = ctx.brho
    pc_mev = ctx.pc_ev / 1e6

    elements, transit_time = _position_and_normalize(flat, ctx, time_ref)

    return {
        "kinetic_energy_mev": kinetic_energy,
        "particle_type": species.name,
        "mass_mev": species.mass_mev,
        "charge": species.charge,
        "pc_mev": pc_mev,
        "brho_tm": brho,
        "rf_frequency_hz": F_RF_DEFAULT,
        "quad_gradient": None,  # Not needed for PALS multipoles
        "source_file": os.path.basename(path),
        "pals_version": version,
        "branch": branch,
        "periodic": bool(periodic),
        "branch_time_s": transit_time,
    }, elements


# -----------------------------------------------------------------------
# Document level
# -----------------------------------------------------------------------

def _check_root(pals, labels):
    if "load" in pals:
        raise UnsupportedFeatureError(
            "PALS 'load' is not supported; combine the lattice into a single "
            "self-contained document")
    _reject_include(pals, "PALS")
    for key, value in pals.items():
        if is_extension_block(value):
            continue
        if key == "post_expansion" and value:
            raise UnsupportedFeatureError(
                "PALS 'post_expansion' is not supported (it modifies the expanded "
                "lattice, which pals2cosy converts as defined in 'facility')")
        if key == "adjustments" and value:
            _warn("PALS 'adjustments' are not applied; the lattice is converted with "
                  "the parameter values as written")
        if key == "girders" and value:
            _warn("PALS 'girders' are not represented (girder alignment shifts are "
                  "ignored)")
        if key not in _ROOT_KEYS and not labels.matches(key):
            _warn(f"Unknown PALS root key '{key}' ignored")
    for reminder in pals.get("reminders") or []:
        warnings.warn(f"PALS reminder: {reminder}", PALSReminder, stacklevel=3)
    return _check_version(pals.get("version"))


def _reject_include(node, where):
    """``include`` may appear at any level of the document; none is supported."""
    if isinstance(node, dict):
        if "include" in node:
            raise UnsupportedFeatureError(
                f"PALS 'include' is not supported (in {where}); combine the lattice "
                f"into a single self-contained document")
        for key, value in node.items():
            if not is_extension_block(value):
                _reject_include(value, f"{where}.{key}")
    elif isinstance(node, list):
        for i, item in enumerate(node):
            _reject_include(item, f"{where}[{i}]")


def _check_version(version):
    if version is None:
        return None
    text = str(version)
    for part in text.split(","):
        part = part.strip()
        try:
            major = int(part.split(".")[0])
        except ValueError:
            _warn(f"Unrecognized PALS version '{part}'")
            continue
        if major not in KNOWN_MAJOR_VERSIONS:
            _warn(f"PALS version '{part}' has major version {major}; pals2cosy targets "
                  f"the PALS specification at commit {__pals_commit__} "
                  f"(major versions {', '.join(map(str, KNOWN_MAJOR_VERSIONS))})")
    return text


class _Facility:
    """Definitions from the ``facility`` list, indexed by name."""

    def __init__(self, facility_list, labels=None):
        if not isinstance(facility_list, list):
            raise ValueError("PALS 'facility' must be a list")
        self.labels = labels or ExtensionLabels()
        self.catalog = {}
        self.order = []
        self.uses = []
        self.definitions = {}
        self._redefinitions = []
        for item in facility_list:
            if not isinstance(item, dict):
                raise ValueError(f"Malformed facility entry {item!r}")
            for key, value in item.items():
                self._entry(key, value)
        self.scope = Scope(self.definitions, self._parameter)
        self._check_redefinitions()
        for name in self.order:
            self._draw_random_values(self.catalog[name])

    def _entry(self, key, value):
        if key == "use":
            self.uses.append(value)
        elif key in ("constants", "variables"):
            for entry in [value] if isinstance(value, dict) else value or []:
                if not isinstance(entry, dict):
                    raise ValueError(f"Malformed '{key}' entry {entry!r}")
                for name, raw in entry.items():
                    self._define(name, raw)
        elif key in ("set", "sets"):
            raise UnsupportedFeatureError(
                f"PALS '{key}' commands are not supported; write the parameter "
                f"values into the element definitions")
        elif key == "superimpose":
            raise UnsupportedFeatureError("PALS superposition is not supported")
        elif key == "include":
            raise UnsupportedFeatureError(
                "PALS 'include' is not supported; combine the lattice into a single "
                "self-contained document")
        elif not isinstance(value, dict):
            if not self.labels.matches(key):
                _warn(f"Facility entry '{key}' is not a definition; ignored")
        elif is_extension_block(value):
            return
        elif value.get("kind") in ("constant", "variable"):
            if "value" not in value:
                raise ValueError(f"{value['kind'].capitalize()} '{key}' has no value")
            self._define(key, value["value"])
        elif value.get("kind") == "Controller":
            control = value.get("control_type", "ABSOLUTE")
            if control == "RELATIVE":
                _warn(f"RELATIVE Controller '{key}' ignored (it does not change the "
                      f"lattice as read)")
            else:
                raise UnsupportedFeatureError(
                    f"{control} Controller '{key}' sets lattice parameters; "
                    f"Controllers are not supported")
        else:
            if key in self.catalog:
                raise ValueError(f"'{key}' is defined more than once in the facility")
            self.catalog[key] = value
            self.order.append(key)

    def _define(self, name, raw):
        if name in self.definitions:
            if self.definitions[name] != raw:
                self._redefinitions.append((name, raw))
        else:
            self.definitions[name] = raw

    def _check_redefinitions(self):
        """A constant or variable may be defined again only with the same value."""
        for name, raw in self._redefinitions:
            first, again = self.scope.lookup(name), self.scope.evaluate_definition(raw)
            if isinstance(first, str) or isinstance(again, str):
                same = (isinstance(first, str) and isinstance(again, str)
                        and parse_species(first).name == parse_species(again).name)
            else:
                same = math.isclose(first, again, rel_tol=1e-12, abs_tol=0.0)
            if not same:
                raise ValueError(
                    f"Constant or variable '{name}' is redefined with a different value "
                    f"({self.definitions[name]!r} and {raw!r})")

    def _draw_random_values(self, node):
        """Evaluate random() draws once per facility node, before branch expansion,
        so that every instance of an element shares the same value."""
        if isinstance(node, dict):
            items = list(node.items())
        elif isinstance(node, list):
            items = list(enumerate(node))
        else:
            return
        for key, value in items:
            if isinstance(value, str) and _RANDOM_CALL.search(value):
                try:
                    result = self.scope.evaluate(value)
                except ExpressionError:
                    continue
                if not isinstance(result, str):
                    node[key] = result
            elif not is_extension_block(value):
                self._draw_random_values(value)

    def _parameter(self, name, path):
        if name not in self.catalog:
            raise KeyError(f"no element '{name}'")
        node = _definition(name, self.catalog)
        for part in path.split("."):
            if not isinstance(node, dict) or node.get(part) is None:
                raise KeyError(f"'{name}' has no parameter {path} set")
            node = node[part]
        return node

    def select_branch(self, beamline_name=None):
        """Return (branch name, root BeamLine name, periodic) to convert."""
        lattice = None
        for target in reversed(self.uses):
            if self.catalog.get(target, {}).get("kind") != "Lattice":
                raise ValueError(f"'use: {target}' does not name a Lattice")
            lattice = target
            break
        if lattice is None:
            lattices = [n for n in self.order if self.catalog[n].get("kind") == "Lattice"]
            lattice = lattices[-1] if lattices else None
        branches = self._branches(lattice) if lattice else []

        if beamline_name is not None:
            for branch in branches:
                if branch[0] == beamline_name:
                    return branch
            defn = self.catalog.get(beamline_name)
            if defn is not None and defn.get("kind") == "BeamLine":
                return beamline_name, beamline_name, defn.get("periodic", False)
            raise ValueError(f"BeamLine '{beamline_name}' not found in facility")

        if branches:
            if len(branches) > 1:
                others = ", ".join(b[0] for b in branches[1:])
                _warn(f"Lattice '{lattice}' has {len(branches)} branches; converting "
                      f"'{branches[0][0]}' (other branches: {others}; select one "
                      f"with --beamline)")
            return branches[0]
        beamlines = [n for n in self.order if self.catalog[n].get("kind") == "BeamLine"]
        if beamlines:
            name = beamlines[-1]
            return name, name, self.catalog[name].get("periodic", False)
        raise ValueError("No BeamLine found in facility")

    def _branches(self, lattice):
        result = []
        for entry in self.catalog[lattice].get("branches") or []:
            if isinstance(entry, str):
                bname, spec = entry, {}
            elif isinstance(entry, dict) and len(entry) == 1:
                bname, spec = next(iter(entry.items()))
                spec = spec or {}
            else:
                raise ValueError(f"Lattice '{lattice}': malformed branch entry {entry!r}")
            if spec.get("kind") == "BeamLine" or "line" in spec:
                if bname in self.catalog and self.catalog[bname] is not spec:
                    raise ValueError(f"'{bname}' is defined more than once")
                self.catalog[bname] = spec
                root = bname
            else:
                root = spec.get("inherit", bname)
            if self.catalog.get(root, {}).get("kind") != "BeamLine":
                raise ValueError(
                    f"Lattice '{lattice}': branch '{bname}' has no root BeamLine '{root}'")
            periodic = spec.get("periodic", self.catalog[root].get("periodic", False))
            result.append((bname, root, periodic))
        return result


def _build_catalog(facility_list):
    """Index facility definitions by name."""
    return _Facility(facility_list).catalog


# -----------------------------------------------------------------------
# Branch expansion
# -----------------------------------------------------------------------

def _merge(base, override):
    """Override ``base`` with ``override``, merging parameter groups one level deep."""
    merged = {**base, **override}
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merged[key] = {**base[key], **value}
    return merged


def _definition(name, catalog, seen=()):
    """Facility definition of ``name`` with element-level ``inherit`` applied."""
    defn = catalog[name]
    target = defn.get("inherit")
    if target is None:
        return defn
    if target == name or target in seen:
        raise ValueError(f"Circular inherit: {' -> '.join(seen + (name, target))}")
    if target not in catalog:
        raise ValueError(f"Unresolved inherit reference: '{target}'")
    base = _definition(target, catalog, seen + (name,))
    return _merge(base, {k: v for k, v in defn.items() if k != "inherit"})


def _expand_beamline(name, catalog, _stack=None, length=None):
    """Recursively expand a BeamLine into a flat list of element dicts.

    ``length`` overrides the BeamLine's own ``length``. When either is given,
    the items are bracketed by BeamLine marks (see ``_close_beamlines``).
    """
    stack = tuple(_stack or ())
    if name in stack:
        raise ValueError(
            f"Circular reference in BeamLine expansion: {' -> '.join(stack + (name,))}")
    stack += (name,)
    defn = catalog[name]
    if defn.get("multipass"):
        _warn(f"BeamLine '{name}' is a multipass line; it is expanded as a plain "
              f"sequence (multipass element sharing and multipass_phase are not "
              f"represented)")
    result = []
    for item in defn.get("line") or []:
        if isinstance(item, str):
            result.extend(_resolve_reference(item, catalog, stack))
        elif isinstance(item, dict):
            for elem_name, spec in item.items():
                if is_extension_block(spec):
                    continue
                if isinstance(spec, dict):
                    result.extend(_resolve_inline(elem_name, spec, catalog, stack))
                else:
                    result.extend(_resolve_reference(elem_name, catalog, stack))
        else:
            raise ValueError(f"BeamLine '{name}': malformed line item {item!r}")
    if length is None:
        length = defn.get("length")
    if length is not None:
        result = [{"kind": _LINE_START, "name": name}, *result,
                  {"kind": _LINE_END, "name": name, "length": length}]
    return result


def _resolve_reference(name, catalog, stack):
    """Resolve a plain reference to an element or sub-beamline."""
    if name not in catalog:
        raise ValueError(f"Unresolved reference: '{name}'")
    defn = _definition(name, catalog)
    kind = defn.get("kind")
    if kind == "BeamLine":
        return _expand_beamline(name, catalog, stack)
    if kind == "Lattice":
        raise ValueError(f"Lattice '{name}' cannot be used as a line item")
    return [{**defn, "name": name}]


def _resolve_inline(name, spec, catalog, stack):
    """Resolve a line item with modifiers or an in-place definition."""
    if "placement" in spec:
        raise UnsupportedFeatureError(
            f"Line item '{name}': 'placement' is not supported")
    direction = spec.get("direction", 1)
    if direction not in (1, -1) or isinstance(direction, bool):
        raise ValueError(
            f"Line item '{name}': direction must be +1 or -1, got {direction!r}")
    if direction == -1:
        raise UnsupportedFeatureError(
            f"Line item '{name}': direction reversal (direction: -1) is not supported; "
            f"it needs reversed-element tracking and reflection patches. A negative "
            f"repeat reverses the element order without reversing direction")
    repeat = spec.get("repeat", 1)
    if isinstance(repeat, bool) or not isinstance(repeat, int):
        raise ValueError(f"Line item '{name}': repeat must be an integer, got {repeat!r}")
    if repeat == 0:
        warnings.warn(f"Element '{name}': repeat=0 produces no elements",
                      PALSConversionWarning, stacklevel=2)

    body = {k: v for k, v in spec.items() if k not in ("repeat", "direction")}
    inherit = body.pop("inherit", None)
    if inherit is not None:
        if inherit not in catalog:
            raise ValueError(f"Unresolved inherit reference: '{inherit}'")
        merged = _merge(_definition(inherit, catalog), body)
    elif name in catalog and "kind" not in body:
        merged = _merge(_definition(name, catalog), body)
    else:
        merged = body

    if merged.get("kind") == "BeamLine":
        if inherit is None and catalog.get(name, {}).get("kind") == "BeamLine" \
                and "line" not in body:
            extra = set(body) - {"kind", "length"}
            if extra:
                warnings.warn(
                    f"Element '{name}': keys {sorted(extra)} ignored on BeamLine reference",
                    PALSConversionWarning, stacklevel=2)
            expanded = _expand_beamline(name, catalog, stack, length=body.get("length"))
        else:
            expanded = _expand_beamline(name, {**catalog, name: merged}, stack)
    else:
        expanded = [{**merged, "name": name}]

    if repeat < 0:
        expanded = expanded[::-1]
    return [dict(e) for _ in range(abs(repeat)) for e in expanded]


# -----------------------------------------------------------------------
# Reference particle
# -----------------------------------------------------------------------

def _reference(flat, facility, ctx, ke_override, particle_override):
    """Return (species, kinetic energy [MeV], reference time [s]) at the start
    of the branch."""
    refp, source = None, None
    tracked = [e for e in flat if e.get("kind") not in (_LINE_START, _LINE_END)]
    beginning = next((e for e in tracked if e.get("kind") == "BeginningEle"), None)
    if beginning is not None:
        refp = ctx.group(beginning, "ReferenceP", beginning["name"])
        source = beginning["name"]
        if beginning is not tracked[0]:
            _warn(f"Branch does not start with a BeginningEle; the reference particle "
                  f"is taken from '{source}' further down the branch")
    else:
        for elem in flat:
            if elem.get("ReferenceP"):
                refp = ctx.group(elem, "ReferenceP", elem["name"])
                source = elem["name"]
                _warn(f"Branch has no BeginningEle; using the ReferenceP of "
                      f"'{source}' as the reference particle")
                break
    refp = refp or {}
    legacy = _legacy_particle(facility)

    if particle_override:
        species = parse_species(particle_override)
    elif refp.get("species_ref") is not None:
        raw = refp["species_ref"]
        if isinstance(raw, str) and raw in facility.definitions:
            raw = facility.scope.lookup(raw)
        species = parse_species(raw)
    elif legacy.get("species") is not None:
        species = parse_species(legacy["species"])
    else:
        warnings.warn(
            "Particle type not specified in lattice file or CLI; defaulting to 'proton'. "
            "Use --particle to specify explicitly.",
            PALSConversionWarning, stacklevel=3)
        species = parse_species("proton")
    if species.charge == 0:
        raise ValueError(
            f"Reference species '{species.name}' is neutral; COSY INFINITY needs a "
            f"charged reference particle")

    mass_ev = species.mass_mev * 1e6
    if ke_override is not None:
        ke = float(ke_override)
    else:
        e_tot = ctx.optional(refp, "E_tot_ref", source, "ReferenceP", 0.0)
        pc = ctx.optional(refp, "pc_ref", source, "ReferenceP", 0.0)
        if e_tot and pc and not math.isclose(e_tot ** 2, pc ** 2 + mass_ev ** 2,
                                             rel_tol=1e-9):
            raise ValueError(
                f"ReferenceP of '{source}': E_tot_ref = {e_tot} eV and pc_ref = {pc} eV "
                f"are inconsistent for {species.name}")
        if e_tot:
            if e_tot < mass_ev:
                raise ValueError(
                    f"ReferenceP of '{source}': E_tot_ref = {e_tot} eV is below the "
                    f"{species.name} rest energy {mass_ev} eV")
            ke = (e_tot - mass_ev) / 1e6
        elif pc:
            ke = (math.sqrt(pc ** 2 + mass_ev ** 2) - mass_ev) / 1e6
        elif legacy.get("kinetic_energy") is not None:
            ke = ctx.number(legacy["kinetic_energy"], "Lattice", "particle.kinetic_energy") / 1e6
        else:
            raise ValueError(
                "Kinetic energy not found in lattice file. Set ReferenceP.E_tot_ref or "
                "pc_ref on the BeginningEle, or use --ke to specify kinetic energy in MeV.")
    if not ke > 0:
        raise ValueError(f"Reference kinetic energy must be positive, got {ke} MeV")
    return species, ke, ctx.optional(refp, "time_ref", source, "ReferenceP", 0.0)


def _legacy_particle(facility):
    """Pre-0.6 ``Lattice.particle`` block (species, kinetic_energy in eV)."""
    for name in facility.order:
        defn = facility.catalog[name]
        if defn.get("kind") == "Lattice" and isinstance(defn.get("particle"), dict):
            warnings.warn(
                f"Lattice '{name}': 'particle' is deprecated; set the reference "
                f"particle with ReferenceP (species_ref, E_tot_ref or pc_ref) on the "
                f"branch's BeginningEle",
                PALSDeprecationWarning, stacklevel=4)
            particle = defn["particle"]
            return {
                "species": particle.get("species", particle.get("species_ref")),
                "kinetic_energy": particle.get("kinetic_energy"),
            }
    return {}


# -----------------------------------------------------------------------
# Element translation and positions
# -----------------------------------------------------------------------

def _position_and_normalize(flat, ctx, time_ref=0.0):
    """Translate elements, assign s positions, resolve RF harmonic numbers.

    Returns the element list and the change of reference time along the
    branch (``t1_ref``, including ReferenceChange time adjustments).
    """
    items = []
    for raw in flat:
        kind = raw.get("kind")
        if kind in (_LINE_START, _LINE_END):
            items.append(_Mark(kind, raw["name"], raw.get("length"),
                               ctx.e_tot_ev, ctx.pc_ev))
            continue
        e_up, pc_up = ctx.e_tot_ev, ctx.pc_ev
        elem = normalize_element(ctx, raw)
        if kind not in _BOOKKEEPING_KINDS:
            items.append((elem, e_up, pc_up, ctx.e_tot_ev, ctx.pc_ev))

    elements = []
    s = 0.0
    t = time_ref
    for elem, e_up, pc_up, e_down, pc_down in _close_beamlines(items, ctx):
        length = elem["length"]
        elem["s_start"] = s
        elem["s_end"] = s + length
        s += length
        if length:
            t += length * (e_up + e_down) / (C_LIGHT * (pc_up + pc_down))
        if elem.get("time_ref") is not None:
            t = elem["time_ref"]
        t += elem.get("dtime_ref") or 0.0
        elements.append(elem)
    transit_time = t - time_ref

    for elem in elements:
        harmonic = elem.get("rf_harmonic")
        if harmonic is not None:
            if transit_time <= 0:
                raise ValueError(
                    f"RFCavity '{elem['name']}': harmonic_number needs a positive "
                    f"reference time along the branch to define the revolution "
                    f"frequency, got {transit_time} s")
            elem["rf_frequency_hz"] = harmonic / transit_time
    return elements, transit_time


@dataclass
class _Mark:
    kind: str
    line: str
    length: object
    e_tot_ev: float
    pc_ev: float


def _close_beamlines(items, ctx):
    """Replace BeamLine marks by drifts that make up each declared length.

    A BeamLine with ``length`` ends that far from its start, so a length
    beyond its items leaves a field-free gap at its end (at the front, when a
    negative repeat reversed the line). A shorter length would make the next
    item overlap this one, which a transfer map cannot represent.
    """
    result, stack = [], []
    for item in items:
        if not isinstance(item, _Mark):
            result.append(item)
            if stack:
                stack[-1]["length"] += item[0]["length"]
            continue
        top = stack[-1] if stack else None
        if top is None or top["mark"].line != item.line or top["mark"].kind == item.kind:
            stack.append({"mark": item, "index": len(result), "length": 0.0})
            continue
        stack.pop()
        end = item if item.kind == _LINE_END else top["mark"]
        declared = ctx.number(end.length, end.line, "length")
        actual = top["length"]
        if math.isclose(declared, actual, rel_tol=1e-9, abs_tol=1e-12):
            gap = 0.0
        elif declared < actual:
            raise UnsupportedFeatureError(
                f"BeamLine '{end.line}': length {declared} m is shorter than its line "
                f"items ({actual} m); overlapping line items are not supported")
        else:
            gap = declared - actual
            drift = blank_element("DRIFT", f"{end.line}_end", gap)
            entry = (drift, end.e_tot_ev, end.pc_ev, end.e_tot_ev, end.pc_ev)
            if end is item:
                result.append(entry)
            else:
                result.insert(top["index"], entry)
        if stack:
            stack[-1]["length"] += actual + gap
    return result
