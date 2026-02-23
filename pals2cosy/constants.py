"""Physical constants for beam physics (CODATA 2018 values).

Author: Eremey Valetov
"""

import math

E0_ELECTRON = 0.51099895000  # MeV
E0_PROTON = 938.27208816     # MeV
G_QUAD_DEFAULT = 2.694       # T/A/m
F_RF_DEFAULT = 2856e6        # Hz

PARTICLES = {
    "electron": (E0_ELECTRON, -1),
    "proton":   (E0_PROTON,    1),
}


def momentum(ke, e0):
    """Relativistic momentum pc in MeV from kinetic energy and rest energy."""
    return math.sqrt(ke**2 + 2 * ke * e0)
