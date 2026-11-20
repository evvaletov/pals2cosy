# Element mapping

This page documents the COSYScript commands generated for each element kind.
Throughout, `r` is the bore radius `--quad-aperture / 2` (default 0.0135 m)
and `Bρ = P0/q` is the signed reference rigidity.

## Drift

```text
DL {length} ;
```

Zero-length bookkeeping elements (`Marker`, `Fiducial`, `Fork`, `FloorShift`)
produce no output; `BeginningEle` and `Placeholder` are not part of the
element list. `Kicker` and `Instrument` become drifts (a nonzero kick or
multipole field warns), and so does a `ReferenceChange` with a length. A
BeamLine whose declared `length` exceeds its items ends with a drift of the
difference, named `<line>_end`.

## Multipoles → MQ, MH, MO, M5, MMS

PALS gives field coefficients `B_N` in T/m^N; COSY takes the pole-tip field

$$B_\mathrm{pole} = \frac{B_N\, r^N}{N!},$$

so that a quadrupole gets $B_1 r$, a sextupole $B_2 r^2/2$, an octupole
$B_3 r^3/6$ and a decapole $B_4 r^4/24$. Normalized strengths are multiplied
by $B\rho$, integrated strengths divided by the length.

```text
MQ {L} {B1·r} {r} ;                         Quadrupole, order 1 only
MH {L} {B2·r²/2} {r} ;                      Sextupole, order 2 only
MO {L} {B3·r³/6} {r} ;                      Octupole, order 3 only
M5 {L} {BQ} {BH} {BO} {BD} {BZ} {r} ;       several normal orders ≤ 5
```

Skew components, tilts, or orders 6–10 use `MMS` with arrays declared in the
`LATTICE` procedure:

```text
PALSBN(1) := ... ; PALSBS(1) := ... ; ...
MMS {L} PALSBN PALSBS {NB} {r} ;
```

with $(B_N' + i B_{sN}') = (B_N + i B_{sN})\, e^{-i(N+1)\,\mathrm{tilt}_N}$.

A zero-length multipole with integrated strengths is a thin kick, represented
by a 1 µm slice between drifts of −0.5 µm (with `FR 0` around it):

```text
DL -5e-07 ; MQ 1e-06 {B1L/1e-06·r} {r} ; DL -5e-07 ;
```

## Solenoid → CMS

```text
CMS {Bsol} {D} {L} ;
```

`Bsol` is `SolenoidP.Bsol`, or `Ksol · Bρ`. COSY's CMS field is always soft
edged,

$$B(s) = \frac{B_\mathrm{sol}}{2\tanh(L/2D)}\left[\tanh\frac{s}{D} - \tanh\frac{s-L}{D}\right],$$

with the fringe length $D$ given by `--solenoid-fringe` (default: the bore
radius $r$). The rotation matches the hard-edge PALS solenoid, but
$\int B^2\,ds \approx B_\mathrm{sol}^2 (L - D)$, so the focusing is weaker by
about $D/L$; at `--fr 0` a warning reports this when it exceeds 0.1%. A
shorter $D$ approaches the hard-edge map but lengthens COSY's integration.

## RF cavity → RF

```text
DL {L/2} ;
PALSRFV(1,1) := {V} ;
RF PALSRFV 0 {f} {φ} {r} ;
DL {L/2} ;
```

- $V$ = `RFP.voltage` (or `gradient · L_active`) in kV
- $f$ = `RFP.frequency`, or `harmonic_number` divided by the reference
  transit time of the branch
- $\varphi = \mathrm{phase}\,[\mathrm{deg}] + 90°$ (`ACCELERATING`), $+0°$
  (`BELOW_TRANSITION`) or $+180°$ (`ABOVE_TRANSITION`), plus 180° for a
  negative reference charge, because COSY gives the reference particle
  $Z V \sin\varphi$.

COSY's RF element is thin, so it sits mid-cavity and the drifts keep the
length.

## Bend → DIL

```text
[RA {tilt_ref} ;] [CB ;] DIL {L} {|θ|} {d} {e1} {h1} {e2} {h2} ; [CB ;] [RA {-tilt_ref} ;]
```

- $L$, $\theta$: arc length and bend angle from the BendP two-of-three
  resolution
- $d$: half gap, which sets the extent of COSY's fringe field at `--fr` 1–3:
  $d = \mathrm{edge\_int} / 0.3695$, so that COSY's default Enge field has
  the PALS fringe-field integral (the mean of `edge1_int` and `edge2_int`, with
  a warning if they differ); otherwise `ApertureP.y_width / 2` (or
  `(y_max − y_min)/2`), or the bore radius $r$ with a warning
- $e_1$, $e_2$: sector face angles (degrees) from `e1`/`e1_rect` per
  `ref_geometry`, multiplied by $\mathrm{sign}(\theta)$ because a negative
  bend is mirrored through `CB`
- $h_1$, $h_2$: face curvatures (1/m, positive convex)
- `RA` rotates the coordinates by `tilt_ref` (degrees) for a tilted bend

### Fringe fields (Enge coefficients)

When `--fr 1` or `--fr 3` is used and a flat-format dipole has Enge
coefficients:

```text
FC 1 1 1 {c₀} {c₁} {c₂} {c₃} {c₄} {c₅} ;
FC 1 2 1 {c₀} {c₁} {c₂} {c₃} {c₄} {c₅} ;
DIL ... ;
FD ;
```

FC is not emitted for FR 0 (no fringe) or FR 2 (symplectic scaling).
PALS bends use COSY's default Enge coefficients, scaled through $d$ to the
given `edge1_int`/`edge2_int`; 0.3695 is $\int F(1-F)\,dS$ of that Enge
function $F$, with $S$ in units of the full gap $2d$. At FR 0 the bend has
hard edges and a nonzero `edge1_int` or `edge2_int` warns.

## Wiggler → WI

```text
WI {B0} {period} {L} {r} 0 0 0 ;
```

From the pals2cosy `WigglerP` extension (`peak_field`, `period`). COSY's `WI`
takes the period λ itself and builds the field `B0 cos(2πz/λ)`.

## Summary table

| Element | COSYScript | Notes |
|---------|-----|-------|
| Drift | `DL L ;` | |
| Quadrupole | `MQ L B1·r r ;` | |
| Sextupole | `MH L B2·r²/2 r ;` | |
| Octupole | `MO L B3·r³/6 r ;` | |
| Multipole (orders 1–5) | `M5 L BQ BH BO BD BZ r ;` | `B_N r^N / N!` |
| Skew, tilted, orders 6–10 | `MMS L PALSBN PALSBS NB r ;` | |
| Thin multipole | `DL -ε/2 ; M.. ε ... ; DL -ε/2 ;` | ε = 1 µm |
| Solenoid | `CMS Bsol D L ;` | soft edge, `D` = `--solenoid-fringe` |
| RFCavity | `DL L/2 ; RF PALSRFV 0 f φ r ; DL L/2 ;` | |
| Wiggler | `WI B period L r 0 0 0 ;` | `WigglerP` extension |
| Bend | `[RA ;] [CB ;] DIL L θ d e1 h1 e2 h2 ; [CB ;] [RA ;]` | |
| Kicker, Instrument | `DL L ;` | a kick or multipole field warns |
| Marker, Fiducial, Fork, FloorShift | _(none)_ | zero length |
