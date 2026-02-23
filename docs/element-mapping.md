# Element mapping

This page documents the exact FOX commands generated for each element type.

## Drift

```fox
DL {length} ;
```

Passive elements with non-zero length (undulators, solenoids) are also emitted
as `DL`. Zero-length elements (BPM, OTR, correctors) produce no output.

## Quadrupole

```fox
MQ {L} {b_pole} {r} ;
```

Where:
- $L$ = quadrupole length (m)
- $r$ = bore radius = `quad_aperture / 2` (default 0.0135 m)
- $b_\mathrm{pole} = \mathrm{sign} \times G \times I \times r$

The sign convention follows COSY INFINITY: **focusing in $x$** (QPF) uses
$\mathrm{sign} = -1$, defocusing (QPD) uses $\mathrm{sign} = +1$.

$G$ is the quadrupole gradient coefficient (T/A/m), read from the lattice
file's `global_settings` or defaulting to 2.694 T/A/m.

## Dipole (without Enge coefficients)

```fox
DIL {L} {θ} {g/2} {e₁} 0 {e₂} 0 ;
```

Where:
- $L$ = effective magnetic length (m)
- $\theta$ = absolute bending angle (degrees)
- $g/2$ = half pole gap (m)
- $e_1$, $e_2$ = entrance and exit edge angles (degrees)

For **negative bending angles**, the element is wrapped with `CB`:

```fox
CB ;
DIL {L} {|θ|} {g/2} {e₁} 0 {e₂} 0 ;
CB ;
```

## Dipole (with Enge coefficients)

Fringe field coefficients are set before the `DIL` command and cleared after:

```fox
FC 1 1 1 {c₀} {c₁} {c₂} {c₃} {c₄} {c₅} ;
DIL {L} {θ} {g/2} {e₁} 0 {e₂} 0 ;
FD ;
```

`FC 1 1 1` sets entrance Enge coefficients. COSY applies them symmetrically to
both entrance and exit by default. If separate exit coefficients are present,
`FC 1 2 1` is emitted before the `DIL`.

`FD` resets all fringe field coefficients after the element.

A negative-angle dipole with Enge coefficients combines both patterns:

```fox
FC 1 1 1 {c₀} ... {c₅} ;
CB ;
DIL {L} {|θ|} {g/2} {e₁} 0 {e₂} 0 ;
CB ;
FD ;
```

## DPW-DPH-DPW consolidation

When the input uses the FELsim DIPOLE_WEDGE representation, consecutive
DPW → SBend → DPW elements are automatically merged into a single dipole:

| Source | Consolidated field |
|--------|--------------------|
| Entrance DPW `wedge_angle_deg` | $e_1$ |
| Main SBend `bending_angle_deg` | $\theta$ |
| Exit DPW `wedge_angle_deg` | $e_2$ |
| Main SBend `length_m` | $L$ |
| Main SBend `pole_gap_m` | $g$ |
| Entrance DPW `enge_coefficients` | entrance Enge |
| Exit DPW `enge_coefficients` | exit Enge |

The physical lengths of the DPW wedge elements are absorbed (they represent
edge-kick regions, not magnetic length). This matches the COSY convention where
edge angles are instantaneous angular kicks at zero path length.

## Summary table

| Element | FOX | Notes |
|---------|-----|-------|
| Drift | `DL L ;` | |
| QPF | `MQ L {-G·I·r} r ;` | Focusing in $x$ |
| QPD | `MQ L {+G·I·r} r ;` | Defocusing in $x$ |
| DPH (no Enge) | `DIL L θ g/2 e₁ 0 e₂ 0 ;` | |
| DPH (with Enge) | `FC ..; DIL ..; FD ;` | |
| DPH (neg. angle) | `CB ; DIL ..; CB ;` | |
| Wiggler/UND | `DL L ;` | Passive drift |
| BPM/OTR/corrector | _(skipped)_ | Zero-length |
