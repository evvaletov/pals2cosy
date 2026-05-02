# Element mapping

This page documents the COSYScript commands generated for each element type.

## Drift

```text
DL {length} ;
```

Zero-length elements (BPM, OTR, correctors) produce no output.

## Quadrupole → MQ

```text
MQ {L} {b_pole} {r} ;
```

- $L$ = quadrupole length (m)
- $r$ = bore radius = `quad_aperture / 2` (default 0.0135 m)
- $b_\mathrm{pole}$ = pole-tip magnetic field (T) from `MagneticMultipoleP.Bn1`

## Sextupole → MH

```text
MH {L} {b_pole} {r} ;
```

- $L$ = sextupole length (m)
- $b_\mathrm{pole}$ = pole-tip field (T) from `MagneticMultipoleP.Bn2`
- $r$ = bore radius (m)

Falls back to `DL` if `Bn2` is not provided.

## Octupole → MO

```text
MO {L} {b_pole} {r} ;
```

- $L$ = octupole length (m)
- $b_\mathrm{pole}$ = pole-tip field (T) from `MagneticMultipoleP.Bn3`
- $r$ = bore radius (m)

Falls back to `DL` if `Bn3` is not provided.

## Multipole → M5

```text
M5 {L} {BQ} {BH} {BO} {BD} {BZ} {r} ;
```

- $L$ = multipole length (m)
- BQ, BH, BO, BD, BZ = pole-tip fields (T) for quad through dodecapole
  (`MagneticMultipoleP.Bn1`–`Bn5`; zero if absent)
- $r$ = bore radius (m)

Falls back to `DL` if all Bn fields are absent.

## Solenoid → CMS

```text
CMS {Bz} {D} {L} ;
```

- $B_z$ = on-axis field (T) from `SolenoidP.Bz`
- $D$ = bore radius (m), used as the Enge falloff parameter
- $L$ = solenoid length (m)

Falls back to `DL` if `Bz` is not provided.

## RF Cavity → RF

```text
RF {V} 0 {W} {PHI} {D} ;
```

- $V$ = peak voltage (kV) from `RFCavityP.voltage_kv`
- 0 = polynomial order (uniform field)
- $W$ = frequency (Hz) from `RFCavityP.frequency_hz`
- $\Phi$ = phase (degrees) from `RFCavityP.phase_deg`
- $D$ = bore radius (m)

Falls back to `DL` if voltage or frequency is absent.

## Dipole → DIL

```text
DIL {L} {θ} {g/2} {e₁} 0 {e₂} 0 ;
```

- $L$ = effective magnetic length (m)
- $\theta$ = absolute bending angle (degrees)
- $g/2$ = half pole gap (m)
- $e_1$, $e_2$ = entrance and exit edge angles (degrees)

For **negative bending angles**, the element is wrapped with `CB`:

```text
CB ;
DIL {L} {|θ|} {g/2} {e₁} 0 {e₂} 0 ;
CB ;
```

### Fringe fields (Enge coefficients)

When `--fr 1` or `--fr 3` is used and the dipole has Enge coefficients:

```text
FC 1 1 1 {c₀} {c₁} {c₂} {c₃} {c₄} {c₅} ;
FC 1 2 1 {c₀} {c₁} {c₂} {c₃} {c₄} {c₅} ;
DIL {L} {θ} {g/2} {e₁} 0 {e₂} 0 ;
FD ;
```

`FC 1 1 1` sets entrance Enge coefficients, `FC 1 2 1` sets exit.
`FD` restores default Enge coefficients after the element.

FC is not emitted for FR 0 (no fringe) or FR 2 (symplectic scaling).

## Summary table

| Element | COSYScript | Notes |
|---------|-----|-------|
| Drift | `DL L ;` | |
| Quadrupole (Bn1) | `MQ L Bn1 r ;` | |
| Sextupole (Bn2) | `MH L Bn2 r ;` | Drift if Bn2 absent |
| Octupole (Bn3) | `MO L Bn3 r ;` | Drift if Bn3 absent |
| Multipole (Bn1–Bn5) | `M5 L BQ BH BO BD BZ r ;` | Up to dodecapole; drift if all absent |
| Solenoid (Bz) | `CMS Bz D L ;` | Drift if Bz absent |
| RFCavity | `RF V 0 W PHI D ;` | Drift if voltage/frequency absent |
| Wiggler/Undulator | `WI B K L D 0 0 0 ;` | Drift if field/period absent |
| SBend/RBend (no Enge) | `DIL L θ g/2 e₁ 0 e₂ 0 ;` | |
| SBend/RBend (with Enge) | `FC ..; DIL ..; FD ;` | |
| SBend (negative angle) | `CB ; DIL ..; CB ;` | |
| Marker | _(skipped)_ | Zero-length |
| BPM/OTR/corrector | _(skipped)_ | Zero-length |
