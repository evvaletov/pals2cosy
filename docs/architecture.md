# Architecture

## Module overview

```
pals2cosy/
├── constants.py      Physical constants (E₀, G, particle masses)
├── _io.py            Shared file loading (JSON, YAML, TOML)
├── parser.py         Flat beamline (beamline: root) → normalized element list
├── pals_parser.py    Official PALS (PALS: root) → normalized element list
├── converter.py      Element list → COSYScript string
├── cli.py            Command-line interface (auto-detection)
└── __main__.py       python -m pals2cosy entry point
```

## Data flow

```
Input file
    │
    ├─ PALS: root ──► pals_parser.parse_lattice(path)
    │                   Expands facility → BeamLine → flat elements
    │                   Converts Bn1, BendP (SI) → internal dict
    │
    └─ beamline: root ► parser.parse_lattice(path)
                         Reads flat elements array
                         Inserts drifts from s_start/s_end gaps
    │
    ▼
(beam_params, elements) — same tuple format from both parsers
    │
    ├─► beam_params dict
    │     kinetic_energy_mev, particle_type,
    │     quad_gradient, rf_frequency_hz, ...
    │
    └─► elements list[dict]
          type, name, length, current, bn1,
          angle, wedge_angle, entrance_edge_angle,
          exit_edge_angle, pole_gap, dipole_length,
          enge_coeffs
    │
    ▼
converter.convert(beam_params, elements, ...)
    │
    ├─ _consolidate_dipoles()
    │    Normalizes dipole representations → DIPOLE_CONSOLIDATED
    │
    ├─ _generate_lattice_body()
    │    DRIFT      → DL L ;
    │    QPF/QPD    → MQ L b r ; (Bn1 or G·I·r)
    │    DIPOLE_CON → FC/CB/DIL/CB/FD
    │    UND        → DL L ;
    │    passive    → DL L ; (or skip if zero-length)
    │
    └─ _fox_template()
         Wraps lattice body in COSY boilerplate
    │
    ▼
Complete COSYScript string
```

## Parsers

### Official PALS parser (`pals_parser.py`)

Parses files with the `PALS:` root key following the
[official PALS specification](https://pals-project.readthedocs.io/). Key features:

- **Catalog-based:** Builds a name→definition index from the `facility:` list.
- **BeamLine expansion:** Recursively resolves `line:` references, `inherit:`
  overrides, and `repeat:` repetitions into a flat element sequence.
- **SI unit conversion:** `BendP.g_ref` (1/m) → degrees, `BendP.e1`/`e2` (rad) → degrees.
- **Bn1 pass-through:** `MagneticMultipoleP.Bn1` stored directly in the `bn1` field.
- **Cumulative positioning:** `s_start`/`s_end` computed from element lengths
  (no explicit positions in official PALS).

### Flat beamline parser (`parser.py`)

Parses files with the `beamline:` root key containing a flat `elements:` array
with explicit `s_start_m`/`s_end_m` positions. Supports both PALS CamelCase
type names (SBend, Quadrupole, etc.) and additional type aliases.

Drifts are inserted between elements based on gaps between `s_end_m` of one
element and `s_start_m` of the next.

Both parsers produce the same `(beam_params, elements)` tuple for the converter.

## Converter

`converter.py` performs two passes over the element list:

1. **Consolidation:** Detects dipole representations and normalizes them
   into `DIPOLE_CONSOLIDATED` elements for uniform handling.

2. **COSYScript generation:** Iterates over the consolidated list and emits the
   appropriate COSYScript commands for each element type. The complete program
   is assembled from a template that includes variable declarations, the
   lattice procedure, setup commands (`OV`, `RPE`/`RPP`/etc., `FR`), and
   optionally Twiss output (`--twiss`).

## Design decisions

**Position-based Enge assignment.** When consolidating dipole representations,
Enge coefficients are assigned by position: entrance coefficients from the
leading element, exit coefficients from the trailing element.

**FR-dependent FC emission.** Enge fringe field coefficients (`FC` commands) are
emitted only for FR modes that use the Enge function via the ODE integrator:
FR 1 (soft-edge, two-step) and FR 3 (high-precision). FR 2 uses symplectic
scaling (pre-computed data from `SYSCA.DAT`), not Enge coefficients, so FC is
omitted. FR 0 has no fringe fields. Default Enge coefficients are loaded
automatically by `OV` (which calls `FD` internally). Per-element `FC` commands
override the defaults for dipoles with custom Enge, and the per-element `FD`
after the dipole restores them.

**Particle presets.** The converter auto-selects the most specific COSY preset
procedure when available (`RPE`, `RPP`, `RPMU`, `RPPI`), falling back to the
generic `RP` command with explicit mass (AMU) and charge for particles without
a dedicated preset (positron, antiproton, deuteron, alpha, etc.).

## Scope

**Included:**
- Official PALS format (facility/BeamLine composition, Bn1, BendP)
- Flat beamline format (positioned elements with current_a)
- Auto-detection from root key
- PALS-standard dipoles (SBend/RBend with edge angle attributes)
- JSON, YAML, and TOML input
- Quadrupoles (Bn1 or current), dipoles, undulators, diagnostics
- FR 0 (no fringe), FR 1 (soft-edge Enge), FR 2 (symplectic scaling), FR 3 (high-precision Enge)
- Enge FC coefficients emitted for FR 1 and FR 3 (ODE integrator modes); omitted for FR 2 (uses SYSCA.DAT)
- 15 particle types (electron, proton, muon, pion±, positron, antiproton, deuteron, alpha, carbon12)
- Transfer map computation; optional Twiss extraction (`--twiss`)

**Deferred:**
- Combined-function magnets (dipole + multipole fields → COSY `MC`/`MS`)
- Electrostatic deflectors (`EC`, `ECL`, `ES`)
- Skew multipoles (`MMS`, `EMS`)
- Particle tracking (`SR` + `TR` + `PRAY`/`WRAY`)
- FIT optimization blocks
