"""Internal element representation shared by the parsers and the converter.

Multipole fields come in two forms:

- ``bn`` / ``bs``: normal and skew PALS field coefficients ``B_N`` in
  T/m^N, with any multipole tilt already folded in. The field is
  ``By + i Bx = sum_N (bn[N] + i bs[N]) (x + i y)^N / N!``. For a ``thin``
  element (zero length) they hold the length-integrated values in T/m^(N-1).
- ``b_pole``: pole-tip fields in T at the converter's aperture radius, the
  convention of the flat ``beamline:`` (FELsim) format.

Angles are in degrees, lengths in meters, edge curvatures in 1/m.

Author: Eremey Valetov
"""

MAX_MULTIPOLE_ORDER = 10  # COSY MMS holds at most 10 normal + 10 skew orders


def blank_element(etype, name="", length=0.0, s_start=None, s_end=None):
    return {
        "type": etype,
        "name": name,
        "length": length,
        "s_start": s_start,
        "s_end": s_end,
        "current": None,
        "angle": None,
        "wedge_angle": None,
        "entrance_edge_angle": None,
        "exit_edge_angle": None,
        "entrance_edge_curvature": None,
        "exit_edge_curvature": None,
        "tilt": 0.0,
        "pole_gap": None,
        "edge_int": None,
        "dipole_length": None,
        "enge_coeffs": None,
        "bn": {},
        "bs": {},
        "b_pole": {},
        "thin": False,
        "normalized": False,
        "bz": None,
        "rf_voltage_kv": None,
        "rf_frequency_hz": None,
        "rf_phase_deg": None,
        "wiggler_field": None,
        "wiggler_period": None,
    }
