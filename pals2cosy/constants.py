"""Physical constants for beam physics (CODATA 2018 values).

Author: Eremey Valetov
"""

import math

C_LIGHT = 299792458.0          # m/s
AMU_MEV = 931.49410242         # MeV/c^2, unified atomic mass unit

E0_ELECTRON = 0.51099895000  # MeV
E0_POSITRON = 0.51099895000  # MeV
E0_PROTON = 938.27208816     # MeV
E0_ANTIPROTON = 938.27208816 # MeV
E0_NEUTRON = 939.56542052    # MeV
E0_MUON = 105.6583755        # MeV (μ⁻ and μ⁺)
E0_TAU = 1776.86             # MeV (PDG)
E0_PION_CHARGED = 139.57039  # MeV (π⁺ and π⁻, PDG)
E0_DEUTERON = 1875.61294257  # MeV
E0_TRITON = 2808.92113298    # MeV
E0_HELION = 2808.39160743    # MeV (³He²⁺)
E0_ALPHA = 3727.3794066      # MeV (⁴He²⁺)
E0_CARBON12 = 12 * AMU_MEV - 6 * E0_ELECTRON  # MeV (¹²C⁶⁺, atomic mass minus electrons)
G_QUAD_DEFAULT = 2.694       # T/A/m
F_RF_DEFAULT = 2856e6        # Hz

# COSY INFINITY's own atomic mass unit (from the CODATA 1998 constants in
# cosy.fox). Masses handed to RP are divided by this so that COSY's rest
# energy reproduces the intended value exactly.
AMU_MEV_COSY = 1.66053873e-27 * 2.99792458e8 ** 2 / 1.602176462e-19 * 1e-6

# Constants predefined by the PALS standard, in PALS units (eV, m, s).
_R_ELECTRON = 2.8179403262e-15
PALS_CONSTANTS = {
    "pi": math.pi,
    "twopi": 2 * math.pi,
    "c_light": C_LIGHT,
    "h_planck": 4.135667696e-15,                 # eV s
    "hbar": 6.582119569e-16,                     # eV s
    "k_boltzmann": 8.617333262e-5,               # eV/K
    "r_electron": _R_ELECTRON,                   # m
    "r_proton": _R_ELECTRON * E0_ELECTRON / E0_PROTON,
    "e_charge": 1.602176634e-19,                 # C
    "mu_0": 1.25663706212e-6 * 1.602176634e-19,  # eV s^2/m (charge in units of e)
    "epsilon_0": 8.8541878128e-12 / 1.602176634e-19,  # 1/(eV m)
    "classical_radius_factor": _R_ELECTRON * E0_ELECTRON * 1e6,  # m eV
    "fine_structure": 7.2973525693e-3,
    "n_avogadro": 6.02214076e23,
}

# Map of recognized particle species (lowercase) to rest mass in MeV.
# Kept for API compatibility; pals2cosy.species is the authoritative parser.
PARTICLE_MASSES_MEV = {
    "electron": E0_ELECTRON,
    "positron": E0_POSITRON,
    "proton": E0_PROTON,
    "antiproton": E0_ANTIPROTON,
    "muon": E0_MUON,
    "mu-": E0_MUON,
    "mu+": E0_MUON,
    "antimuon": E0_MUON,
    "pion+": E0_PION_CHARGED,
    "pi+": E0_PION_CHARGED,
    "pion-": E0_PION_CHARGED,
    "pi-": E0_PION_CHARGED,
    "deuteron": E0_DEUTERON,
    "alpha": E0_ALPHA,
    "carbon12": E0_CARBON12,
}
