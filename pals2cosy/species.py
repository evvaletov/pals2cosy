"""Particle species names: the PALS (OpenPMD SpeciesType) convention.

PALS names species as OpenPMD does, with an optional charge state appended:
``electron``, ``positron``, ``proton``, ``anti-proton``, ``muon``,
``anti-muon``, ``#3He+2`` (helium-3, fully stripped), ``Au+79``,
``#12C+6``. A bare element or isotope (``#3He``) is a neutral atom.
pals2cosy also accepts the short names of earlier releases (``antiproton``,
``mu+``, ``pi-``, ``alpha``, ``carbon12``, ...) and charged pions.

Ion masses are atomic masses minus the removed electrons; electron binding
energies (below 1e-5 relative even for bare uranium) are neglected.

Author: Eremey Valetov
"""

import re
from dataclasses import dataclass
from typing import Optional

from .constants import (
    AMU_MEV, AMU_MEV_COSY, E0_ALPHA, E0_DEUTERON, E0_ELECTRON, E0_HELION,
    E0_MUON, E0_NEUTRON, E0_PION_CHARGED, E0_PROTON, E0_TAU, E0_TRITON,
)


@dataclass(frozen=True)
class Species:
    name: str
    mass_mev: float
    charge: int
    cosy_preset: Optional[str] = None
    anomalous_moment: Optional[float] = None

    def cosy_command(self, kinetic_energy_mev):
        """COSYScript command setting this species as the reference particle."""
        ke = kinetic_energy_mev
        if self.cosy_preset == "RPPI":
            return f"RPPI {ke} {self.charge}"
        if self.cosy_preset:
            return f"{self.cosy_preset} {ke}"
        return f"RP {ke} {self.mass_mev / AMU_MEV_COSY} {self.charge}"


# Canonical OpenPMD names, then aliases accepted for compatibility.
_ELEMENTARY = {
    "electron": Species("electron", E0_ELECTRON, -1, "RPE", 1.15965218128e-3),
    "positron": Species("positron", E0_ELECTRON, +1, None, 1.15965218128e-3),
    "muon": Species("muon", E0_MUON, -1, "RPMU", 1.16592089e-3),
    "anti-muon": Species("anti-muon", E0_MUON, +1, None, 1.16592089e-3),
    "tau": Species("tau", E0_TAU, -1),
    "anti-tau": Species("anti-tau", E0_TAU, +1),
    "proton": Species("proton", E0_PROTON, +1, "RPP", 1.79284734463),
    "anti-proton": Species("anti-proton", E0_PROTON, -1, None, 1.79284734463),
    "neutron": Species("neutron", E0_NEUTRON, 0),
    "anti-neutron": Species("anti-neutron", E0_NEUTRON, 0),
    "pion+": Species("pion+", E0_PION_CHARGED, +1, "RPPI", 0.0),
    "pion-": Species("pion-", E0_PION_CHARGED, -1, "RPPI", 0.0),
    "deuteron": Species("deuteron", E0_DEUTERON, +1, None, -0.1429873),
    "triton": Species("triton", E0_TRITON, +1),
    "helion": Species("helion", E0_HELION, +2),
    "alpha": Species("alpha", E0_ALPHA, +2),
}

_ALIASES = {
    "mu-": "muon", "mu": "muon",
    "mu+": "anti-muon", "antimuon": "anti-muon",
    "antiproton": "anti-proton", "pbar": "anti-proton",
    "antineutron": "anti-neutron",
    "pi+": "pion+", "pi-": "pion-",
    "e-": "electron", "e+": "positron",
    "carbon12": "#12C+6",
}

# Atomic masses [u] of neutral atoms (NIST / AME2020).
_ISOTOPE_MASS_U = {
    ("H", 1): 1.00782503223, ("H", 2): 2.01410177812, ("H", 3): 3.0160492779,
    ("He", 3): 3.0160293201, ("He", 4): 4.00260325413,
    ("Li", 6): 6.0151228874, ("Li", 7): 7.0160034366,
    ("Be", 9): 9.012183065,
    ("B", 10): 10.01293695, ("B", 11): 11.00930536,
    ("C", 12): 12.0, ("C", 13): 13.00335483507, ("C", 14): 14.0032419884,
    ("N", 14): 14.00307400443, ("N", 15): 15.00010889888,
    ("O", 16): 15.99491461957, ("O", 18): 17.99915961286,
    ("F", 19): 18.99840316273,
    ("Ne", 20): 19.9924401762, ("Ne", 22): 21.991385114,
    ("Na", 23): 22.9897692820,
    ("Mg", 24): 23.985041697,
    ("Al", 27): 26.98153853,
    ("Si", 28): 27.97692653465,
    ("P", 31): 30.97376199842,
    ("S", 32): 31.9720711744,
    ("Cl", 35): 34.968852682,
    ("Ar", 36): 35.967545105, ("Ar", 40): 39.9623831237,
    ("K", 39): 38.9637064864,
    ("Ca", 40): 39.962590863, ("Ca", 48): 47.95252276,
    ("Ti", 48): 47.94794198,
    ("Cr", 52): 51.94050623,
    ("Fe", 56): 55.93493633,
    ("Ni", 58): 57.93534241, ("Ni", 64): 63.92796682,
    ("Cu", 63): 62.92959772,
    ("Zn", 64): 63.92914201,
    ("Kr", 78): 77.92036494, ("Kr", 84): 83.9114977282, ("Kr", 86): 85.9106106269,
    ("Zr", 90): 89.9046977,
    ("Ag", 107): 106.9050916, ("Ag", 109): 108.9047553,
    ("Sn", 120): 119.90220163,
    ("Xe", 124): 123.905892, ("Xe", 129): 128.9047808611,
    ("Xe", 132): 131.9041550856, ("Xe", 136): 135.907214484,
    ("Au", 197): 196.96656879,
    ("Pb", 208): 207.9766525,
    ("Bi", 209): 208.9803991,
    ("Th", 232): 232.0380558,
    ("U", 235): 235.0439301, ("U", 238): 238.0507884,
}

# Atomic number and the isotope used when no mass number is given.
_ELEMENTS = {
    "H": (1, 1), "He": (2, 4), "Li": (3, 7), "Be": (4, 9), "B": (5, 11),
    "C": (6, 12), "N": (7, 14), "O": (8, 16), "F": (9, 19), "Ne": (10, 20),
    "Na": (11, 23), "Mg": (12, 24), "Al": (13, 27), "Si": (14, 28),
    "P": (15, 31), "S": (16, 32), "Cl": (17, 35), "Ar": (18, 40),
    "K": (19, 39), "Ca": (20, 40), "Ti": (22, 48), "Cr": (24, 52),
    "Fe": (26, 56), "Ni": (28, 58), "Cu": (29, 63), "Zn": (30, 64),
    "Kr": (36, 84), "Zr": (40, 90), "Ag": (47, 107), "Sn": (50, 120),
    "Xe": (54, 132), "Au": (79, 197), "Pb": (82, 208), "Bi": (83, 209),
    "Th": (90, 232), "U": (92, 238),
}

_ATOM_RE = re.compile(
    r"^(?:#(?P<A>\d+))?(?P<sym>[A-Z][a-z]?)(?P<chg>[+-]\d+|\++|-+)?$")


def parse_species(name):
    """Return the :class:`Species` named by a PALS species string.

    Raises ``ValueError`` for names that are not recognized or for isotopes
    without tabulated mass data.
    """
    if isinstance(name, Species):
        return name
    if not isinstance(name, str) or not name.strip():
        raise ValueError(f"Unknown particle species {name!r}")
    text = name.strip()
    atom = _ALIASES.get(text.lower(), text)
    if atom.lower() in _ELEMENTARY:
        return _ELEMENTARY[atom.lower()]
    m = _ATOM_RE.match(atom)
    if m is None:
        raise ValueError(
            f"Unknown particle species {name!r}. Use an OpenPMD name such as "
            f"'electron', 'proton', 'anti-muon', or an ion such as '#12C+6', 'Au+79'"
        )
    sym = m.group("sym")
    if sym not in _ELEMENTS:
        raise ValueError(
            f"Particle species {name!r}: no mass data for element '{sym}'. "
            f"Tabulated elements: {', '.join(_ELEMENTS)}"
        )
    z, default_a = _ELEMENTS[sym]
    a = int(m.group("A")) if m.group("A") else default_a
    if (sym, a) not in _ISOTOPE_MASS_U:
        raise ValueError(
            f"Particle species {name!r}: no mass data for isotope #{a}{sym}"
        )
    chg = m.group("chg")
    if not chg:
        charge = 0
    elif chg[1:].isdigit():
        charge = int(chg)
    else:
        charge = len(chg) if chg[0] == "+" else -len(chg)
    if charge > z:
        raise ValueError(
            f"Particle species {name!r}: charge {charge} exceeds atomic number {z}"
        )
    mass = _ISOTOPE_MASS_U[(sym, a)] * AMU_MEV - charge * E0_ELECTRON
    canonical = f"#{a}{sym}" + (f"{charge:+d}" if charge else "")
    if (sym, a, charge) == ("H", 1, 1):
        return _ELEMENTARY["proton"]
    return Species(canonical, mass, charge)


def mass_of(name):
    """Rest mass in eV/c^2 (PALS units)."""
    return parse_species(name).mass_mev * 1e6


def charge_of(name):
    """Charge in units of the elementary charge."""
    return float(parse_species(name).charge)


def anomalous_moment_of(name):
    sp = parse_species(name)
    if sp.anomalous_moment is None:
        raise ValueError(f"Anomalous magnetic moment of {sp.name!r} is not tabulated")
    return sp.anomalous_moment
