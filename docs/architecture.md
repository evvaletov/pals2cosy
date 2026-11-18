# Architecture

## Module overview

```
pals2cosy/
├── _version.py       Package version and the targeted PALS commit
├── constants.py      Physical constants (CODATA 2018), PALS named constants
├── species.py        PALS/OpenPMD species names → mass, charge, COSY RP command
├── expressions.py    Safe evaluator for PALS arithmetic expressions
├── diagnostics.py    Warning and error categories
├── _io.py            Shared file loading (JSON, YAML, TOML)
├── _elements.py      Internal element dict shared by parsers and converter
├── parser.py         Flat beamline (beamline: root) → element list
├── pals_parser.py    Official PALS document: facility, lattices, expansion, reference
├── pals_elements.py  PALS element kinds → element dicts (field resolution)
├── converter.py      Element list → COSYScript string
├── cli.py            Command-line interface (auto-detection, --strict)
└── __main__.py       python -m pals2cosy entry point
```

## Data flow

```
Input file
    │
    ├─ PALS: root ──► pals_parser.parse_lattice(path)
    │                   Root keys, version, extension labels
    │                   facility → catalog, constants/variables (expressions.Scope)
    │                   Lattice/branch selection, BeamLine expansion
    │                   Reference particle (BeginningEle ReferenceP, species.py)
    │                   pals_elements.normalize_element() per element:
    │                     KnN·Bρ, BnNL/L, tilt folding, Bend geometry, RFP, SolenoidP
    │
    └─ beamline: root ► parser.parse_lattice(path)
                         Reads flat elements array
                         Inserts drifts from s_start/s_end gaps
    │
    ▼
(beam_params, elements) — same tuple format from both parsers
    │
    ├─► beam_params dict
    │     kinetic_energy_mev, particle_type, mass_mev, charge,
    │     pc_mev, brho_tm, branch, periodic, pals_version, ...
    │
    └─► elements list[dict] (see _elements.blank_element)
          type, name, length, s_start, s_end,
          bn/bs (T/m^N) or b_pole (T), thin,
          angle, edge angles/curvatures, tilt, pole_gap,
          bz, rf_voltage_kv, rf_frequency_hz, rf_phase_deg, ...
    │
    ▼
converter.convert(beam_params, elements, ...)
    │
    ├─ _consolidate_dipoles()
    │    Normalizes dipole representations → DIPOLE_CONSOLIDATED
    │
    ├─ _LatticeWriter
    │    DRIFT/passive → DL L ;
    │    multipoles    → MQ/MH/MO/M5 (normal) or MMS (skew, orders > 5),
    │                    pole tip B_N r^N / N!; thin slice for zero length
    │    DIPOLE_CON    → FC/RA/CB/DIL/CB/RA/FD
    │    SOL           → CMS ;  RFC → DL/RF/DL ;  UND → WI
    │
    └─ _fox_template()
         Wraps the lattice body (with any array declarations) in COSY boilerplate
    │
    ▼
Complete COSYScript string
```

## Parsers

### Official PALS parser (`pals_parser.py`, `pals_elements.py`)

- **Catalog-based:** builds a name→definition index from the `facility:`
  list, applying element-level `inherit` lazily.
- **Branch expansion:** recursively resolves `line:` references, in-place
  definitions, `inherit:` overrides, and `repeat:` (negative counts reverse
  the order).
- **Expressions:** numeric fields go through `expressions.Scope`, which
  resolves constants, variables, PALS constants and element parameter
  references with cycle detection.
- **Field resolution:** strengths are reduced to PALS field coefficients
  `bn`/`bs` in T/m^N with the signed rigidity at each element; the converter
  alone knows the bore radius and forms pole-tip fields.
- **Diagnostics:** unsupported constructs raise `UnsupportedFeatureError`;
  partial representations warn with `PALSConversionWarning`; deprecated input
  warns with `PALSDeprecationWarning`.

### Flat beamline parser (`parser.py`)

Parses files with the `beamline:` root key containing a flat `elements:` array
with explicit `s_start_m`/`s_end_m` positions. Supports both PALS CamelCase
type names (SBend, Quadrupole, etc.) and additional type aliases. Its
`MagneticMultipoleP.Bn1`–`Bn5` are pole-tip fields (`b_pole`).

Drifts are inserted between elements based on gaps between `s_end_m` of one
element and `s_start_m` of the next.

## Converter

`converter.py` performs two passes over the element list:

1. **Consolidation:** Detects dipole representations and normalizes them
   into `DIPOLE_CONSOLIDATED` elements for uniform handling.

2. **COSYScript generation:** Iterates over the consolidated list and emits the
   appropriate COSYScript commands for each element type. Elements that need
   COSY arrays (`MMS`, `RF`) make the writer declare them at the top of the
   `LATTICE` procedure. The complete program is assembled from a template
   that includes the lattice procedure, setup commands (`OV`,
   `RPE`/`RPP`/`RP`/etc., `FR`), and optionally Twiss output (`--twiss`).

## Design decisions

**Field coefficients, not pole tips, in the element list.** The parsers store
what the lattice says (PALS `B_N`), and the converter alone applies the
aperture convention of COSY INFINITY's multipole elements. The flat format,
which already stores pole-tip fields, keeps them in a separate key.

**Signed rigidity.** `Bρ = P0/q` carries the sign of the charge, so a given
normalized strength focuses the same way for electrons and positrons and the
COSY reference particle (whose charge enters COSY's own rigidity) sees the
intended field.

**Thin elements as slices.** COSY has no thin multipole; a 1 µm slice between
−0.5 µm drifts reproduces the kick with an error of order `K L · 1 µm` and
keeps the branch length exact.

**Position-based Enge assignment.** When consolidating dipole representations,
Enge coefficients are assigned by position: entrance coefficients from the
leading element, exit coefficients from the trailing element.

**FR-dependent FC emission.** Enge fringe field coefficients (`FC` commands) are
emitted only for FR modes that use the Enge function via the ODE integrator:
FR 1 (soft-edge, two-step) and FR 3 (high-precision). FR 2 uses symplectic
scaling (pre-computed data from `SYSCA.DAT`), not Enge coefficients, so FC is
omitted. FR 0 has no fringe fields.

**Particle presets.** The converter selects COSY's preset procedures where they
exist (`RPE`, `RPP`, `RPMU`, `RPPI`), falling back to `RP` with the mass in
COSY's own atomic mass unit, so COSY's rest energy equals the intended one.

## Scope

**Included:**
- PALS documents at spec commit `b9d5b355`: facility, lattices, branches,
  `use`, BeamLines, `inherit`, `repeat`, constants, variables, expressions
- Reference particle from `BeginningEle.ReferenceP`; species parser
- Quadrupoles, sextupoles, octupoles, general multipoles (normal, skew,
  tilted, normalized, integrated, thin), bends, solenoids, RF cavities,
  wigglers
- Flat beamline format (positioned elements with `current_a` or pole tips)
- FR 0/1/2/3; transfer map computation; optional Twiss extraction (`--twiss`)
- Conformance suite over the official PALS examples

**Deferred:**
- `include`/`load`, `set`, controllers, superposition, placement, direction
  reversal, forks and multipass bookkeeping
- Combined-function magnets (dipole + multipole fields → COSY `MC`/`MS`)
- Electrostatic elements
- Misalignments (`BodyShiftP`)
- Particle tracking (`SR` + `TR` + `PRAY`/`WRAY`), FIT optimization blocks
