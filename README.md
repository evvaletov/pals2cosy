# pals2cosy

Convert [PALS](https://pals-project.readthedocs.io/) lattice files to
[COSY INFINITY](https://cosyinfinity.org) input code (COSYScript). The
generated program computes the transfer map of one lattice branch.

pals2cosy 0.6 targets the PALS specification at commit
[`b9d5b355`](https://github.com/pals-project/pals/tree/b9d5b35555a63f7ba015460e1c17ac929d35e215)
(13 November 2026); `pals2cosy --version` and `pals2cosy.__pals_commit__`
report it. Every official example at that commit is part of the test suite
and either converts or fails with a documented reason (see
[Conformance](#conformance)). A flat `beamline:` format with positioned
elements is also accepted.

## Installation

```bash
pip install pals2cosy
```

Or for development:

```bash
pip install -e .
```

Requires Python 3.9+ and PyYAML (for YAML input).

## Usage

```bash
# PALS FODO channel (auto-detected from the PALS: root key)
pals2cosy examples/fodo.pals.yaml

# Write to file with fringe fields
pals2cosy beamline.pals.yaml -o output.fox --fr 1

# Select a Lattice branch or a BeamLine from the facility (a BeamLine without
# a BeginningEle needs the reference particle on the command line)
pals2cosy lattice.pals.yaml --beamline fodo_cell --ke 1000 --particle proton

# Override the reference particle; compute Twiss parameters (periodic lattices)
pals2cosy lattice.pals.yaml --ke 1000 --particle proton --twiss

# Fail instead of warning when something is not represented in COSYScript
pals2cosy lattice.pals.yaml --strict
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o FILE` | stdout | Output COSYScript file |
| `--ke FLOAT` | from lattice | Reference kinetic energy (MeV) |
| `--particle STR` | from lattice | Reference species (see [Reference particle](#reference-particle)) |
| `--fr INT` | 0 | Fringe field order (0=none, 1=Enge, 2=symplectic, 3=high-precision Enge) |
| `--order INT` | 3 | DA computation order |
| `--dim INT` | 3 (2 with `--twiss`) | Phase-space dimensions (2 or 3) |
| `--quad-aperture FLOAT` | 0.027 | Multipole bore diameter (m); pole-tip fields are given at its radius |
| `--solenoid-fringe FLOAT` | bore radius | Fringe-field length (m) of COSY's soft-edge solenoid (see [Solenoids](#element-mapping)) |
| `--twiss` | off | Append GT Twiss extraction (requires periodic lattice; a 3D map also needs RF) |
| `--beamline NAME` | first branch | Lattice branch or BeamLine to expand (official PALS only) |
| `--strict` | off | Turn conversion warnings into errors |
| `--mode STR` | auto | Input format override |
| `--no-comments` | off | Omit element name comments in COSYScript output |
| `--version` | | Print the pals2cosy version and the targeted PALS commit |

Exit status is 1, with a one-line message, when a lattice cannot be converted.

## Input format

```yaml
PALS:
  version: null
  facility:
    - start:
        kind: BeginningEle
        ReferenceP:
          species_ref: proton
          pc_ref: 1.0e9          # [eV]; or E_tot_ref
    - q1:
        kind: Quadrupole
        length: 0.5
        MagneticMultipoleP:
          Kn1: 0.4               # [1/m^2]; or Bn1 [T/m], Kn1L, Bn1L
    - b1:
        kind: Bend
        length: 2.0
        BendP:
          angle_ref: 0.1         # [rad]; any two of {g_ref, Bn0_ref, radius_ref},
          e1_rect: 0.0           #   {length, L_chord, L_rectangle}, {angle_ref}
          e2_rect: 0.0
    - my_line:
        kind: BeamLine
        line:
          - start
          - q1
          - b1
          - q1:
              MagneticMultipoleP:
                Kn1: -0.4
```

See `examples/` for a FODO channel with `inherit:` and `repeat:` and a UH FEL
beamline excerpt.

## Supported PALS features

- **Document:** `facility`, `version` (null, `0.x`, `1.x`; another major
  version warns), `extension_labels` (declared names, prefixes and suffixes are
  passed over silently), `authors`, `notes`, `reminders` (shown as
  `PALSReminder` warnings), `phase_space_coordinates`.
- **Lattices and lines:** `Lattice` with `branches` (by name, with
  `inherit`, or defined in place), `use` (the last one wins; default is the
  last Lattice), `BeamLine` with in-place definitions, element and line
  `inherit`, parameter-group `inherit` (named groups and `element>Group`; an
  inherited parameter may be modified in another form, such as `Bn1` for
  `Kn1`),
  `repeat` (a negative repeat reverses the element order without direction
  reversal, so a Bend keeps `e1` on its upstream face), and BeamLine `length`
  (a drift fills the line up to its declared length; a shorter length, which
  would make line items overlap, raises `UnsupportedFeatureError`).
- **Extensions:** extension blocks (`name: {extension: type, ...}`) at any
  level and the names, prefixes and suffixes declared in `extension_labels`
  are passed over silently.
- **Values:** constants and variables (`constants:`/`variables:` lists and
  `kind: constant`/`variable`; a repeated definition must have the same
  value), arithmetic expressions in any numeric field (`+ - * / ^ **`,
  parentheses, unary signs, `Inf`), the PALS constants (`pi`, `twopi`,
  `c_light`, `r_electron`, ...), the PALS functions (`sqrt`, `exp`, `log`,
  trigonometric and hyperbolic functions, `atan2`, `abs`, `factorial`,
  `int`, `nint`, `floor`, `ceiling`, `modulo`, `sign`, `mass_of`,
  `charge_of`, `anomalous_moment_of`, `random`, `random_gauss`, `expr`), and
  element parameter references (`Q1>MagneticMultipoleP.Kn1`). Random values
  are drawn once per facility definition, so every instance of an element
  shares them. Expressions are parsed by a small recursive-descent
  evaluator; Python `eval` is never used.
- **Reference particle:** `ReferenceP` of the branch's `BeginningEle`
  (`species_ref`, `E_tot_ref` or `pc_ref`, `time_ref`), with
  reference-energy changes through RF `dE_ref` and reference-time changes
  through `ReferenceChange` (`dtime_ref`, `time_ref`).
- **Elements:** see the mapping below.

## Element mapping

PALS multipole coefficients are plain field derivatives,
`By + i Bx = Σ (BnN + i BsN) (x + i y)^N / N!` with `BnN` in T/m^N, while
COSY INFINITY's multipole elements take the field at the pole tip of a bore
of radius `r` (`--quad-aperture / 2`). pals2cosy converts with

    B_pole(N) = BnN · r^N / N!          (quad Bn1·r, sextupole Bn2·r²/2,
                                          octupole Bn3·r³/6, decapole Bn4·r⁴/24)

Normalized strengths use the signed reference rigidity
`Bρ = P0 / q = pc / (c q)` (`BnN = KnN · Bρ`, so the field of a given `KnN`
changes sign with the charge), and length-integrated strengths are divided by
the element length. At `--fr 0`, `--quad-aperture` therefore changes only
the radius at which COSY is told the field, not the gradient; with `--fr` 1–3
it also sets the extent of the multipole fringe fields, and so the optics.

| PALS | COSYScript | Notes |
|------|------------|-------|
| `Drift` | `DL L ;` | |
| `Quadrupole` | `MQ L Bn1·r r ;` | `Kn1`, `Bn1L`, `Kn1L` converted to `Bn1` |
| `Sextupole` | `MH L Bn2·r²/2 r ;` | |
| `Octupole` | `MO L Bn3·r³/6 r ;` | |
| `Multipole`, or several orders in one element | `M5 L B1 B2 B3 B4 B5 r ;` | `B_N = BnN·r^N/N!`, N = 1–5 |
| skew components, tilts, or orders 6–10 | `MMS L PALSBN PALSBS NB r ;` | `(BnN' + i BsN') = (BnN + i BsN)·exp(−i(N+1)·tiltN)` |
| zero-length multipole with `BnNL`/`KnNL` | `DL −ε/2 ; M.. ε … ; DL −ε/2 ;` | thin kick, ε = 1 µm, local `FR 0` |
| `Bend` | `[RA τ ;] [CB ;] DIL L abs(θ) d e1 h1 e2 h2 ; [CB ;] [RA −τ ;]` | see below |
| `Solenoid` | `CMS Bsol D L ;` | `Bsol = Ksol·Bρ`; soft edge, see below |
| `RFCavity` | `DL L/2 ; RF PALSRFV 0 f φ r ; DL L/2 ;` | see below |
| `Wiggler` (pals2cosy `WigglerP`) | `WI B period L r 0 0 0 ;` | |
| `Kicker`, `Instrument`, `Drift` | `DL L ;` | a nonzero kick or multipole field warns |
| `Marker`, `Fiducial`, `Fork`, `FloorShift`, `BeginningEle`, `Placeholder` | — | zero-length bookkeeping |
| `ReferenceChange` | `DL L ;` if it has a length | reference-time changes only |

**Bends.** The arc length `L`, curvature `g` and angle `θ` are resolved from
any two of the three BendP sets {`g_ref`, `Bn0_ref` = `g_ref`·Bρ,
`radius_ref`}, {`length`, `L_chord` = 2 sin(θ/2)/g, `L_rectangle` = sin θ/g},
{`angle_ref`}; a consistent third value is accepted, an inconsistent one is an
error. Face angles follow `ref_geometry`: `e1 = e1_rect + θ/2` (ARC, CHORD),
`e1 = e1_rect` / `e2 = e2_rect + θ` (ENTRANCE_COORDS),
`e1 = e1_rect + θ` / `e2 = e2_rect` (EXIT_COORDS), with the spec's defaults
(sector faces for ARC, rectangular otherwise). A negative angle is bent
through `CB`, which mirrors x, so the face angles of a negative bend are
emitted as `sign(θ)·e` (a rectangular bend has `e1 = e2 = abs(θ)/2` in COSY for
either sign). Face curvatures `h1`, `h2` map to the DIL curvature arguments
(positive is convex in both codes). `tilt_ref` rotates the bend about the
entrance axis and is emitted as `RA` rotations (degrees) around it.

The half-gap `d` sets the extent of COSY's fringe field at `--fr` 1–3. It is
chosen so that COSY's default Enge field has the fringe-field integrals
`edge1_int`/`edge2_int`: `d = edge_int / 0.3695` (from their mean, with a
warning, if they differ). Without them, `d` is half of `ApertureP.y_width`
(or `y_max − y_min`), or else the bore radius (`--quad-aperture / 2`) with a
warning. At `--fr 0` COSY's bends have hard edges and a nonzero `edge1_int`
or `edge2_int` warns. In the official IOTA example
(`tests/pals_examples/iota.pals.yaml`, `edge_int` = 0.0145 m) the bends alone
give the fractional tunes (0.3002, 0.3029) at `--fr` 1–3 against the
design (0.300, 0.300), but (0.300, 0.404) at `--fr 0`.

**Solenoids.** COSY INFINITY's `CMS` solenoid always has soft edges,
`B(s) ∝ tanh(s/D) − tanh((s − L)/D)` with fringe length `D`, while the PALS
solenoid is hard-edged. The rotation is the same, but the focusing integral
`∫B² ds` is smaller by about `D·B²`, so the focusing is weaker by about
`D/L`. `D` is `--solenoid-fringe`, which defaults to the bore radius. At
`--fr 0` a warning gives the deficit whenever it exceeds 0.1%; a shorter
`--solenoid-fringe` approaches the hard-edge map, at the cost of COSY
integration time (a few seconds per solenoid for `D = 10⁻⁴ m`).

**RF cavities.** COSY's `RF` is a thin kick, so it is placed between two
drifts of half the cavity length. The voltage is `RFP.voltage` (or
`gradient · L_active`, `L_active` defaulting to the length) in kV; the
frequency is `RFP.frequency`, or `harmonic_number / T` with `T` the reference
transit time of the converted branch. COSY gives the reference particle
`Z·V·sin(φ)`, so the COSY phase is
`φ = phase [deg] + 90 (ACCELERATING) | 0 (BELOW_TRANSITION) | 180 (ABOVE_TRANSITION)`,
plus 180 for a negatively charged reference. COSY then re-references the
energy to that gain; if it differs from `RFP.dE_ref` (PALS default 0) a
warning says so; the same holds for a cavity with zero or no voltage but a
nonzero `dE_ref`. The reference energy used for normalized strengths
downstream follows `dE_ref`, as in PALS. The RF phase of each cavity is
taken relative to the arrival of the reference particle at that cavity,
which is COSY's convention; this is how pals2cosy reads the
`2π·t·frequency + phase + phi0` of the specification, which does not fix the
origin of `t`. `harmonic_number` is converted with the branch reference time
`t1_ref`: the straight-line transit time of the branch (the specification's
formula) plus `ReferenceChange` time adjustments.

## Reference particle

The reference particle is read from `ReferenceP` on the branch's
`BeginningEle`, which must be the first element of the branch (a later one
is used with a warning): `species_ref` (or the name of a constant holding
one) and `E_tot_ref` or `pc_ref` in eV. `--particle` and `--ke` (MeV)
override them; with only `--particle`, the file's `E_tot_ref`/`pc_ref` is
kept and the kinetic energy follows from the new mass. Normalized strengths
are always resolved with the reference actually used.

Species names follow PALS (OpenPMD) spelling: `electron`, `positron`,
`proton`, `anti-proton`, `muon`, `anti-muon`, `tau`, `deuteron`, `triton`,
`helion`, `alpha`, charged pions `pion+`/`pion-`, and ions written as an
optional isotope and a charge state, such as `#3He+2`, `#12C+6`, `Au+79`
(the most abundant isotope when none is given) or `H-`. Ion masses are atomic
masses minus the removed electrons (binding energies neglected). The short
names of earlier releases (`antiproton`, `mu+`, `antimuon`, `pi-`,
`carbon12`, ...) still work. Every species is set with `RP E M Z`, the
mass in COSY's atomic mass unit, so that COSY's reference rigidity is the one
normalized strengths were converted with; COSY's presets (`RPP`, `RPE`,
`RPMU`, `RPPI`) carry masses that differ by up to 1.5e-6. Where the anomalous
magnetic moment `G` is tabulated (electron, positron, muon, anti-muon,
proton, anti-proton, deuteron), `RPS 0 G` then sets COSY's `G0`, as the presets do, with spin
tracking off. A neutral reference species is an error.

## Unsupported constructs

Constructs that would change the converted lattice raise
`UnsupportedFeatureError` (a subclass of both `ValueError` and
`NotImplementedError`): `set`/`sets`, `post_expansion`, ABSOLUTE
`Controller`s, `include` (at any level of the document), `load`,
`superimpose`, line-item `placement`, direction reversal (`direction: -1`),
overlapping line items (a BeamLine `length` shorter than its items), `Patch`
elements with offsets, and `ReferenceChange` of the energy or species.

Constructs that leave the converted branch intact are converted with a
`PALSConversionWarning`: additional Lattice branches (select one with
`--beamline`) and `Fork`s (not followed), `multipass` lines (expanded as plain
sequences), `adjustments` (not applied), RELATIVE controllers, `girders` and
`BodyShiftP` misalignments, `ElectricMultipoleP`, multipoles inside bends,
bend dipole errors (`Kn0` differing from the reference), `edge1_int`/`edge2_int`
at `--fr 0`, kicks, RF
`multipass_phase` and traveling-wave cavities, element kinds without a COSY
counterpart (emitted as drifts), multipole fields on elements that do not
convert them (drifts, instruments, wigglers, markers, solenoids, RF
cavities), unknown keys, including unknown components of the parameter
groups pals2cosy reads (`BendP`, `RFP`, `SolenoidP`, `ReferenceP`,
`ReferenceChangeP`, `ApertureP`, `PatchP`, `ForkP`, `BodyShiftP`, the
multipole groups), and active elements with no readable strength (emitted as
drifts). `--strict` turns these into errors.

## Leniencies and model differences

Where pals2cosy is more lenient or stricter than the specification, or COSY
INFINITY models an element differently:

- `BendP` accepts a redundant parameter from the same set as long as it is
  consistent (`g_ref` with `radius_ref`, or `e1` with `e1_rect`), although the
  specification asks for two parameters from different sets.
- A branch whose first element is not its `BeginningEle` takes the reference
  particle from a later one, with a warning; a `BeginningEle` that opens a
  subline at the start of the branch is accepted silently.
- Solenoids are soft-edged in COSY (see [Solenoids](#element-mapping)).
- Thin multipoles are 1 µm slices, and RF cavities thin kicks between
  half-length drifts.
- Negative element lengths are rejected; the specification has no such rule.
- An RF cavity with a voltage needs `frequency` or `harmonic_number` to build
  COSY's `RF` element; several official examples give neither.

## Deprecations

Accepted with a `PALSDeprecationWarning` (a `FutureWarning`, shown by
default):

| Deprecated | Use instead |
|------------|-------------|
| `kind: SBend`, `kind: RBend` | `kind: Bend` (both are read as a Bend with `ref_geometry: ARC`, as in 0.5) |
| `Lattice.particle` (`species`, `kinetic_energy` in eV) | `ReferenceP` on the `BeginningEle` |
| `SolenoidP.Bz` | `SolenoidP.Bsol` |
| `RFCavityP` (`voltage_kv`, `frequency_hz`, `phase_deg` in COSY's convention) | `RFP` |

## pals2cosy extensions

`WigglerP` (`peak_field` [T], `period` [m]) on a `Wiggler` has no PALS
counterpart yet. Files that use it can declare it, as the specification's
extension mechanism provides:

```yaml
PALS:
  extension_labels:
    names:
      WigglerP: pals2cosy wiggler parameters (peak_field [T], period [m])
```

The flat `beamline:` format (FELsim) is a pals2cosy format; in it,
`MagneticMultipoleP.Bn1`–`Bn5` are pole-tip fields in tesla, as FELsim writes
them, not PALS field derivatives.

## Conformance

`tests/pals_examples/` holds the `examples/` tree of the PALS repository at
the commit in `tests/pals_examples/SPEC_COMMIT`, with a `manifest.yaml`
listing for every file whether it converts or which documented error it
raises. The suite checks that the manifest covers every example, that each
outcome holds, and, when COSY INFINITY is installed, that every converted
example runs in COSY at `--fr 0`, `1`, `2` and `3` with the expected total
length.
Further tests compare COSY transfer maps with the PALS definitions (reference
rigidity, multipole kicks to fifth order, skew and tilted multipoles against
COSY's own frame rotation, bend direction and tilts, sector and rectangular
bends of both signs, RF phase conventions, solenoids against the hard-edge
map, wiggler focusing).

## Known limitations

- One branch per conversion; forks, multipass bookkeeping and other branches
  are not followed.
- Combined-function bends, electrostatic elements and misalignments are not
  converted.
- Solenoids are soft-edged; their focusing is weaker than the hard-edge PALS
  solenoid by about `--solenoid-fringe / L`.
- `--fr` applies to every element. PALS multipoles have no fringe-field
  parameters, so with `--fr` 1–3 they get COSY's default fringe fields at the
  bore radius; in the IOTA example these move the fractional tunes from
  (0.300, 0.303) to (0.270, 0.286).
- COSY INFINITY 10.2 compiles programs of up to about 5,000 element commands.
  Beyond that it stops with `!!! MEMORY EXHAUSTION, INCREASE LCON` and still
  exits with status 0; larger lattices need a COSY build with a larger
  `LCON`.
- Enge coefficients: COSY's `FC` accepts at most 6 per fringe field; longer
  lists are truncated with a warning.
- Flat-format quadrupoles given by current (`±G·I·r`) use an electron sign
  convention; PALS input specifies fields directly.

## Roadmap

- Combined-function magnets (dipole + multipole fields → COSY `MC`/`MS`)
- Electrostatic elements (`ElectricMultipoleP` → `EQ`/`EM`, deflectors)
- `include` and `load`, `set`, ABSOLUTE controllers
- Misalignments (`BodyShiftP` → `SA`/`RA`)
- Particle tracking (`SR` + `TR` + `PRAY`/`WRAY`) and FIT optimization blocks

## Author

Eremey Valetov

- Department of Physics and Astronomy, Michigan State University, East Lansing, MI, USA
- Department of Physics and Astronomy, University of Hawaiʻi at Mānoa, Honolulu, HI, USA

## License

MIT; see `LICENSE`. The PALS examples in `tests/pals_examples/` are under
CC BY 4.0 (`tests/pals_examples/LICENSE`).

## Acknowledgments

This work was supported in part by DOE EPSCoR Award DE-SC0025583. The converter
draws on the author's work on the University of Hawaiʻi free electron laser.
The official examples in `tests/pals_examples/` come from the PALS project
(github.com/pals-project/pals) under CC BY 4.0.

## Tests

```bash
pip install pytest
pytest tests/ -v
```

Tests that run COSY INFINITY are skipped when `/usr/local/bin/cosy` and
`cosy.fox` are not found (override with `PALS2COSY_COSY_BIN` and
`PALS2COSY_COSY_LIBRARY`). The `--fr 1` and `--fr 2` runs also need COSY's
`SYSCA.DAT`, looked for next to `cosy.fox` or given with
`PALS2COSY_COSY_SYSCA`.
