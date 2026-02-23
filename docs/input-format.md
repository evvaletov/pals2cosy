# Input format

pals2cosy reads PALS v2 lattice files in JSON or YAML format. The file
extension (`.json`, `.yaml`, `.yml`) determines the parser used.

## Structure

A v2 lattice file has this top-level structure:

```yaml
beamline:
  metadata:
    format_version: 2
    name: My_Beamline
    version: '1.0'
    description: Example beamline
    reference_energy_mev: 45.0
    particle_type: electron
  beam_parameters:
    particle:
      type: electron
      kinetic_energy_mev: 45.0
      mass_mev: 0.51099895
      charge_e: -1
    rf_frequency_hz: 2856000000.0
  elements:
    - ...
  global_settings:
    quadrupole_gradient_coefficient_t_per_a_per_m: 2.694
```

## Element types

### Quadrupole

```yaml
- name: Q1
  type: Quadrupole
  s_start_m: 0.5
  s_end_m: 0.589
  length_m: 0.089
  polarity: focusing      # or "defocusing"
  parameters:
    current_a: 1.5
  aperture_m: 0.027
```

### Dipole (PALS-native)

Edge angles and Enge coefficients are attributes of the SBend element:

```yaml
- name: B1
  type: SBend
  s_start_m: 1.0
  s_end_m: 1.037
  length_m: 0.037389
  parameters:
    bending_angle_deg: 11.25
    dipole_length_m: 0.037389
    pole_gap_m: 0.0127
    entrance_edge_angle_deg: 0.0
    exit_edge_angle_deg: 11.25
  fringe_fields:
    enge_coefficients: [56.49, -50.79, 19.32, -3.621, 0.3315, -0.01193]
```

### Dipole (FELsim DIPOLE_WEDGE)

Edge kicks are separate elements bracketing the main dipole. The converter
auto-detects consecutive DPW → DPH → DPW patterns and consolidates them:

```yaml
- name: W1_entrance
  type: DIPOLE_WEDGE
  s_start_m: 0.99
  s_end_m: 1.00
  length_m: 0.01
  parameters:
    wedge_angle_deg: 0.0
    dipole_angle_deg: 11.25
    dipole_length_m: 0.037389
    pole_gap_m: 0.0127
  fringe_fields:
    enge_coefficients: [56.49, -50.79, 19.32, -3.621, 0.3315, -0.01193]

- name: B1
  type: SBend
  s_start_m: 1.00
  s_end_m: 1.037
  length_m: 0.037389
  parameters:
    bending_angle_deg: 11.25
    dipole_length_m: 0.037389
    pole_gap_m: 0.0127

- name: W1_exit
  type: DIPOLE_WEDGE
  s_start_m: 1.037
  s_end_m: 1.047
  length_m: 0.01
  parameters:
    wedge_angle_deg: 11.25
    dipole_angle_deg: 11.25
    dipole_length_m: 0.037389
    pole_gap_m: 0.0127
```

Both representations produce identical FOX output.

### Undulator / Wiggler

Treated as a drift of the same length:

```yaml
- name: U1
  type: Wiggler
  s_start_m: 12.0
  s_end_m: 12.54
  length_m: 0.54
  parameters: {}
```

### Diagnostics and correctors

Zero-length elements (BPM, OTR, correctors, spectrometers) are skipped in the
FOX output. The drift space around them is preserved.

```yaml
- name: BPM1
  type: BPM
  s_start_m: 2.457
  s_end_m: 2.457
  length_m: 0.0
  parameters: {}
```

## Drift insertion

Drifts are automatically inserted between elements based on the gap between
each element's `s_end_m` and the next element's `s_start_m`. There is no need
to include explicit drift elements in the lattice file.

## Type aliases

The following type names are recognized (case-sensitive):

| PALS name | FELsim name | Internal | FOX command |
|-----------|-------------|----------|-------------|
| `Quadrupole` | `QUADRUPOLE` | QPF/QPD | `MQ` |
| `SBend`, `RBend` | `DIPOLE`, `DPH` | DPH | `DIL` |
| — | `DIPOLE_WEDGE`, `DPW` | DPW | _(consolidated)_ |
| `Wiggler` | `UNDULATOR`, `UND` | UND | `DL` |
| `Kicker` | `CORRECTOR_V/H` | STV/STH | _(skipped)_ |
| `Instrument` | `BPM`, `OTR`, `SPC` | BPM/OTR/SPC | _(skipped)_ |
| `Marker` | — | DRIFT | `DL` |
