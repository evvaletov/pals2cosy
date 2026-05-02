"""Physical constants for beam physics (CODATA 2018 values).

Author: Eremey Valetov
"""

E0_ELECTRON = 0.51099895000  # MeV
E0_POSITRON = 0.51099895000  # MeV
E0_PROTON = 938.27208816     # MeV
E0_ANTIPROTON = 938.27208816 # MeV
E0_MUON = 105.6583755        # MeV (μ⁻ and μ⁺)
E0_PION_CHARGED = 139.57039  # MeV (π⁺ and π⁻)
E0_DEUTERON = 1875.61294257  # MeV
E0_ALPHA = 3727.3794066      # MeV (⁴He²⁺)
E0_CARBON12 = 11177.92922    # MeV (¹²C⁶⁺ rest mass; nuclear)
G_QUAD_DEFAULT = 2.694       # T/A/m
F_RF_DEFAULT = 2856e6        # Hz

# Map of recognized particle species (lowercase) to rest mass in MeV.
# The converter accepts a richer set of names than the official PALS
# specification; this table keeps beam_params consistent with the
# COSYScript reference-particle commands emitted by converter._reference_particle().
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
