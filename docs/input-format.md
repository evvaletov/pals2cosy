# Input format

pals2cosy reads lattice files in PALS format, auto-detected from the `PALS:`
root key. A flat `beamline:` format with positioned elements is also accepted.

---

## Official PALS format

Files with the `PALS:` root key follow the
[PALS specification](https://pals-project.readthedocs.io/); pals2cosy 0.6
targets commit `b9d5b355` of
[pals-project/pals](https://github.com/pals-project/pals).

### Structure

```yaml
PALS:
  version: null
  facility:
    - start:
        kind: BeginningEle
        ReferenceP:
          species_ref: electron
          E_tot_ref: 40.51099895e6     # [eV]; or pc_ref [eV]
    - element_name:
        kind: ElementKind
        length: 1.0
        ...
    - line_name:
        kind: BeamLine
        line:
          - start
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

The branch converted is the first branch of the Lattice named by the last
`use` statement (or of the last Lattice); `--beamline` picks another branch or
any BeamLine. A negative `repeat` reverses the order of the repeated elements
without reversing their direction. A BeamLine `length` longer than its items
leaves a field-free gap at its end.

### Reference particle

The `ReferenceP` group of the branch's `BeginningEle` sets the reference
species (`species_ref`) and energy (`E_tot_ref` or `pc_ref`, in eV).
`--ke` (MeV) and `--particle` override them. Species use the PALS (OpenPMD)
spelling: `electron`, `positron`, `proton`, `anti-proton`, `muon`,
`anti-muon`, `deuteron`, `alpha`, `pion+`, ions such as `#3He+2` or `Au+79`.

### Constants, variables and expressions

```yaml
- constants:
    - L_q: 0.25
- variables:
    - k_f: 0.8 / L_q
- qf:
    kind: Quadrupole
    length: 2 * L_q
    MagneticMultipoleP:
      Kn1: k_f
- qd:
    inherit: qf
    MagneticMultipoleP:
      Kn1: -qf>MagneticMultipoleP.Kn1
```

Any numeric field may be an expression built from numbers, `+ - * / ^ **`,
parentheses, the PALS constants and functions, constants and variables, and
element parameter references.

### Multipoles (MagneticMultipoleP)

Strengths are PALS field coefficients of order `N`
(`By + i Bx = Σ (BnN + i BsN)(x + i y)^N / N!`), in any of the forms

| Key | Units | Meaning |
|-----|-------|---------|
| `BnN`, `BsN` | T/m^N | normal / skew field coefficient |
| `KnN`, `KsN` | 1/m^(N+1) | normalized: `BnN = KnN · Bρ` |
| `BnNL`, `KnNL`, ... | T/m^(N−1), 1/m^N | length-integrated |
| `tiltN` | rad | rotation of the order-N multipole |
| `BnN_taper`, ... | | added to the corresponding strength |

```yaml
- q1:
    kind: Quadrupole
    length: 0.5
    MagneticMultipoleP:
      Bn1: 2.0          # T/m; positive focuses positive particles horizontally
- s1:
    kind: Sextupole
    length: 0.2
    MagneticMultipoleP:
      Kn2L: 0.3
- sq:
    kind: Quadrupole
    length: 0.1
    MagneticMultipoleP:
      Kn1: 0.5
      tilt1: 0.7853981634   # a skew quadrupole
```

A quadrupole, sextupole, octupole or multipole may carry any orders from 1 to
10. A zero-length element with integrated strengths is a thin multipole.

### Bends (BendP)

```yaml
- b1:
    kind: Bend
    length: 2.0
    BendP:
      angle_ref: 0.1        # any two of {g_ref, Bn0_ref, radius_ref},
      e1_rect: 0.0          #   {length, L_chord, L_rectangle}, {angle_ref}
      e2_rect: 0.0
      ref_geometry: ARC     # ARC, CHORD, ENTRANCE_COORDS, EXIT_COORDS
      h1: 0.0               # face curvatures [1/m]
      tilt_ref: 0.0         # [rad]
    ApertureP:
      y_width: 0.03         # full gap [m]
```

### Solenoids (SolenoidP)

```yaml
- sol1:
    kind: Solenoid
    length: 0.5
    SolenoidP:
      Bsol: 1.5    # T; or Ksol [1/m] = Bsol / Bρ
```

### RF cavities (RFP)

```yaml
- cav1:
    kind: RFCavity
    length: 0.5
    RFP:
      voltage: 1.0e6           # V; or gradient [V/m] with L_active [m]
      frequency: 1.3e9         # Hz; or harmonic_number
      phase: 0.0               # rad
      zero_phase: ACCELERATING # or BELOW_TRANSITION, ABOVE_TRANSITION
      dE_ref: 1.0e6            # eV, change of the reference energy
```

### Parameter mapping

| PALS field | Normalized dict | Conversion |
|---|---|---|
| `MagneticMultipoleP.BnN`, `BsN` (T/m^N) | `bn[N]`, `bs[N]` | tilt folded in, `KnN·Bρ`, `BnNL/L` |
| `SolenoidP.Bsol` / `Ksol` | `bz` | `Ksol·Bρ` |
| `RFP.voltage` (V) | `rf_voltage_kv` | ÷ 1000 |
| `RFP.frequency` / `harmonic_number` | `rf_frequency_hz` | `h / T_branch` |
| `RFP.phase` (rad), `zero_phase` | `rf_phase_deg` | COSY sine convention |
| `BendP` geometry | `angle` (deg), `length` | two-of-three resolution |
| `BendP.e1`, `e2` (rad) | edge angles (deg) | in the bend's own (`CB`) frame |
| `BendP.h1`, `h2` (1/m) | edge curvatures | direct |
| `BendP.tilt_ref` (rad) | `tilt` | direct |
| `WigglerP.peak_field`, `period` | `wiggler_field`, `wiggler_period` | direct (pals2cosy extension) |

---

## Flat `beamline:` format

The flat format, written by FELsim, lists positioned elements
(`s_start_m`, `s_end_m`) and inserts drifts between them. In this format
`MagneticMultipoleP.Bn1`–`Bn5` are pole-tip fields in tesla, passed to COSY
unchanged, and quadrupoles may instead give a current (`parameters.current_a`)
converted with the gradient coefficient `G`.
