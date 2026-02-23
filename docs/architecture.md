# Architecture

## Module overview

```
pals2cosy/
├── constants.py    Physical constants (E₀, G, momentum)
├── parser.py       v2 JSON/YAML → normalized element list
├── converter.py    Element list → FOX code string
├── cli.py          Command-line interface
└── __main__.py     python -m pals2cosy entry point
```

## Data flow

```
YAML/JSON file
    │
    ▼
parser.parse_lattice(path)
    │
    ├─► beam_params dict
    │     kinetic_energy_mev, particle_type,
    │     quad_gradient, rf_frequency_hz, ...
    │
    └─► elements list[dict]
          Each dict has: type, name, length,
          current, angle, wedge_angle,
          entrance_edge_angle, exit_edge_angle,
          pole_gap, dipole_length, enge_coeffs
          (Drifts auto-inserted from s_start/s_end gaps)
    │
    ▼
converter.convert(beam_params, elements, ...)
    │
    ├─ _consolidate_dipoles()
    │    DPW-DPH-DPW → DIPOLE_CONSOLIDATED
    │    Standalone DPH → DIPOLE_CONSOLIDATED
    │
    ├─ _generate_lattice_body()
    │    DRIFT      → DL L ;
    │    QPF/QPD    → MQ L b r ;
    │    DIPOLE_CON → FC/CB/DIL/CB/FD
    │    UND        → DL L ;
    │    passive    → (skip)
    │
    └─ _fox_template()
         Wraps lattice body in COSY boilerplate
    │
    ▼
Complete FOX code string
```

## Parser

`parser.py` is self-contained with no FELsim dependency. It replicates the
type resolution logic from FELsim's `latticeLoaderBase.py`:

- `Quadrupole` + `polarity` → `QPF` or `QPD`
- `SBend` / `RBend` / `DIPOLE` → `DPH`
- `DIPOLE_WEDGE` → `DPW`
- `Wiggler` / `UNDULATOR` → `UND`
- `Kicker` + `plane` → `STV` or `STH`
- `Instrument` + `instrument_type` → `BPM`, `OTR`, `SPC`
- `Marker` → zero-length `DRIFT`

Drifts are inserted between elements based on gaps between `s_end_m` of one
element and `s_start_m` of the next.

## Converter

`converter.py` performs two passes over the element list:

1. **Consolidation:** Scans for consecutive DPW → DPH → DPW patterns and
   merges them into `DIPOLE_CONSOLIDATED` elements. Standalone DPH elements
   (PALS-native SBend with edge angles) are also wrapped as
   `DIPOLE_CONSOLIDATED` for uniform handling.

2. **FOX generation:** Iterates over the consolidated list and emits the
   appropriate FOX commands for each element type. The complete FOX program
   is assembled from a template that includes variable declarations, the
   lattice procedure, setup commands (`OV`, `RPE`, `FR`), and Twiss output.

## Design decisions

**No FELsim dependency.** The converter is intended as a standalone community
tool for the COSY INFINITY website. It copies the minimum necessary logic
(type resolution, drift insertion, DPW consolidation) rather than importing
FELsim modules.

**Position-based Enge assignment.** When consolidating DPW-DPH-DPW triplets,
Enge coefficients from the entrance DPW are assigned as entrance fringe fields,
and coefficients from the exit DPW as exit fringe fields. This is more robust
than the sign-based detection used in FELsim's `cosySimulator`.

**v2 format only.** The converter targets the PALS-aligned v2 format. The
legacy v1 format (FELsim-native type names without `s_start_m`/`s_end_m`)
is not supported.

## Scope

**Included in v0.1:**
- PALS-standard dipoles (SBend/RBend with edge angle attributes)
- FELsim DIPOLE_WEDGE auto-detection and consolidation
- v2 JSON + YAML input
- Quadrupoles, dipoles, undulators, diagnostics
- FR 0 (no fringe) and FR 1 (Enge coefficients)
- Transfer map + Twiss output

**Deferred:**
- Solenoid, RF cavity, sextupole (not in UH beamline; stubs only)
- FR 3 / MGE fieldmap support
- Particle tracking (RRAY/WRAY)
- FIT optimization blocks
- v1 format support
