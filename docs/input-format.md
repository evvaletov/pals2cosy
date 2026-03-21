# Input format

pals2cosy reads lattice files in PALS format, auto-detected from the `PALS:`
root key. A flat `beamline:` format with positioned elements is also accepted.

---

## Official PALS format

Files with the `PALS:` root key follow the
[official PALS specification](https://pals-project.readthedocs.io/).

### Structure

```yaml
PALS:
  version: null
  facility:
    - element_name:
        kind: ElementType
        length: 1.0
        ...
    - line_name:
        kind: BeamLine
        line:
          - element_name
          - other_element:
              inherit: element_name
              MagneticMultipoleP:
                Bn1: -1.0
    - lattice_name:
        kind: Lattice
        branches:
          - line_name
    - use: lattice_name
```

### Quadrupole (Bn1)

Quadrupoles use `MagneticMultipoleP.Bn1` — the pole-tip field in Tesla:

```yaml
- q1:
    kind: Quadrupole
    length: 1.0
    MagneticMultipoleP:
      Bn1: 1.0    # Tesla; positive = horizontally focusing for positive charge
```

### Sextupole (Bn2)

Sextupoles use `MagneticMultipoleP.Bn2` — the sextupole pole-tip field in Tesla:

```yaml
- s1:
    kind: Sextupole
    length: 0.2
    MagneticMultipoleP:
      Bn2: 0.5    # Tesla
```

### Solenoid (Bz)

Solenoids use `SolenoidP.Bz` — the on-axis field in Tesla:

```yaml
- sol1:
    kind: Solenoid
    length: 0.5
    SolenoidP:
      Bz: 1.5    # Tesla (on-axis field)
```

### RF Cavity

RF cavities use `RFCavityP` with voltage, frequency, and phase:

```yaml
- cav1:
    kind: RFCavity
    length: 0.5
    RFCavityP:
      voltage_kv: 100.0        # peak voltage (kV)
      frequency_hz: 2856.0e6   # RF frequency (Hz)
      phase_deg: 90.0           # RF phase (degrees)
```

### Octupole (Bn3)

Octupoles use `MagneticMultipoleP.Bn3`:

```yaml
- o1:
    kind: Octupole
    length: 0.15
    MagneticMultipoleP:
      Bn3: 0.3    # Tesla
```

### Multipole (Bn1–Bn5)

Generic multipoles can carry up to 5 harmonic fields (quad through dodecapole):

```yaml
- m1:
    kind: Multipole
    length: 0.2
    MagneticMultipoleP:
      Bn1: 0.5    # quadrupole component
      Bn2: 0.1    # sextupole component
      Bn3: 0.05   # octupole component
```

### Dipole (BendP)

Dipoles use `BendP` with `g_ref` (1/m) and edge angles in radians:

```yaml
- b1:
    kind: SBend
    length: 2.0
    BendP:
      g_ref: 0.5    # 1/m → bending angle = g_ref × length
      e1: 0.1        # entrance edge angle (radians)
      e2: 0.2        # exit edge angle (radians)
```

### BeamLine composition

Elements are composed using `line:` references, `inherit:` overrides, and `repeat:`:

```yaml
- fodo_cell:
    kind: BeamLine
    line:
      - drift1
      - quad1
      - drift2:
          kind: Drift
          length: 0.5
      - quad2:
          inherit: quad1
          MagneticMultipoleP:
            Bn1: -1.0
      - drift1

- fodo_channel:
    kind: BeamLine
    line:
      - fodo_cell:
          repeat: 3
```

### Parameter mapping

| PALS field | Normalized dict | Conversion |
|---|---|---|
| `MagneticMultipoleP.Bn1` (T) | `bn1` | Direct |
| `MagneticMultipoleP.Bn2` (T) | `bn2` | Direct |
| `MagneticMultipoleP.Bn3` (T) | `bn3` | Direct |
| `MagneticMultipoleP.Bn4` (T) | `bn4` | Direct |
| `MagneticMultipoleP.Bn5` (T) | `bn5` | Direct |
| `SolenoidP.Bz` (T) | `bz` | Direct |
| `RFCavityP.voltage_kv` (kV) | `rf_voltage_kv` | Direct |
| `RFCavityP.frequency_hz` (Hz) | `rf_frequency_hz` | Direct |
| `RFCavityP.phase_deg` (deg) | `rf_phase_deg` | Direct |
| `WigglerP.peak_field` (T) | `wiggler_field` | Direct |
| `WigglerP.period` (m) | `wiggler_period` | Direct |
| `BendP.g_ref` (1/m) + `length` (m) | `angle` (deg) | $\theta = g_\mathrm{ref} \times L \times 180/\pi$ |
| `BendP.e1` / `BendP.e2` (rad) | edge angles (deg) | $\times 180/\pi$ |

---

## Supported element types

| PALS kind | FOX command | Notes |
|-----------|-------------|-------|
| `Quadrupole` | `MQ L Bn1 r ;` | Bn1 from `MagneticMultipoleP` |
| `Sextupole` | `MH L Bn2 r ;` | Bn2 from `MagneticMultipoleP`; drift if absent |
| `Octupole` | `MO L Bn3 r ;` | Bn3 from `MagneticMultipoleP`; drift if absent |
| `Multipole` | `M5 L BQ BH BO BD BZ r ;` | Bn1–Bn5; drift if all absent |
| `Solenoid` | `CMS Bz D L ;` | Bz from `SolenoidP`; drift if absent |
| `RFCavity` | `RF V 0 W PHI D ;` | Drift if voltage/frequency absent |
| `Wiggler` | `WI B K L D 0 0 0 ;` | Drift if field/period absent |
| `SBend`, `RBend` | `DIL L θ g/2 e₁ 0 e₂ 0 ;` | With optional FC/CB |
| `Drift` | `DL L ;` | |
| `Kicker`, `Instrument`, `Marker` | _(skipped)_ | Zero-length |

## Known limitations

- **Enge coefficient truncation:** COSY INFINITY's `FC` command accepts at most
  6 Enge coefficients. If more are provided, they are truncated with a warning.

- **Bn1 sign convention and particle type:** In the official PALS parser,
  quadrupole polarity labels (QPF/QPD) are assigned from the sign of Bn1 using
  the proton convention (positive Bn1 = horizontally focusing). For electrons,
  the labels are inverted. This is cosmetic only — the converter uses Bn1
  directly for `MQ` commands.

- **Inherit merge depth:** The `inherit:` mechanism merges nested dicts
  (e.g. `MagneticMultipoleP`) one level deep only. Deeper nesting is not
  currently supported by any PALS element type.
