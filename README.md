# pals2cosy

Convert [PALS](https://pals-project.readthedocs.io/) lattice files to
[COSY INFINITY](https://cosyinfinity.org) FOX code.

Supports the **official PALS format** (`PALS:` root with `facility:`/`line:`
composition) with auto-detection from the root key. Also accepts a flat
`beamline:` format with positioned elements.

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
# PALS FODO cell (auto-detected from PALS: root key)
pals2cosy examples/fodo.pals.yaml --ke 1000 --particle proton

# Write to file with fringe fields
pals2cosy beamline.pals.yaml -o output.fox --ke 40 --fr 1

# Select a specific BeamLine from the facility
pals2cosy lattice.pals.yaml --beamline fodo_cell --ke 1000

# Compute Twiss parameters (periodic lattices only)
pals2cosy lattice.pals.yaml --ke 1000 --particle proton --twiss
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o FILE` | stdout | Output FOX file |
| `--ke FLOAT` | from lattice | Kinetic energy (MeV) |
| `--fr INT` | 0 | Fringe field order (0=none, 1=Enge, 2=symplectic, 3=high-precision Enge) |
| `--order INT` | 3 | DA computation order |
| `--dim INT` | 3 | Phase-space dimensions (2 or 3) |
| `--quad-aperture FLOAT` | 0.027 | Quadrupole bore diameter (m) |
| `--particle STR` | from lattice | Particle type (see [Supported particles](#supported-particles)) |
| `--twiss` | off | Append GT Twiss extraction (requires periodic lattice) |
| `--mode STR` | auto | Input format override |
| `--beamline NAME` | last | BeamLine to expand (official PALS only) |
| `--no-comments` | off | Omit element name comments in FOX output |

## Input format

```yaml
PALS:
  version: null
  facility:
    - q1:
        kind: Quadrupole
        length: 1.0
        MagneticMultipoleP:
          Bn1: 1.0           # Pole-tip field (Tesla)
    - b1:
        kind: SBend
        length: 2.0
        BendP:
          g_ref: 0.5         # Bend strength (1/m)
          e1: 0.1            # Entrance edge angle (radians)
          e2: 0.2            # Exit edge angle (radians)
    - my_line:
        kind: BeamLine
        line:
          - q1
          - b1
```

See `examples/fodo.pals.yaml` for a complete FODO cell with `inherit:` and `repeat:`.

## Element mapping

| Input type | FOX command | Notes |
|-----------|-------------|-------|
| Quadrupole (Bn1) | `MQ L Bn1 r ;` | Direct pole-tip field |
| Quadrupole (current) | `MQ L {±G·I·r} r ;` | Sign from polarity |
| Sextupole (Bn2) | `MH L Bn2 r ;` | Drift if Bn2 absent |
| Octupole (Bn3) | `MO L Bn3 r ;` | Drift if Bn3 absent |
| Multipole (Bn1–Bn5) | `M5 L BQ BH BO BD BZ r ;` | Up to dodecapole; drift if all absent |
| Solenoid (Bz) | `CMS Bz D L ;` | Drift if Bz absent |
| RFCavity | `RF V 0 W PHI D ;` | Drift if voltage/frequency absent |
| Wiggler/Undulator | `WI B K L D 0 0 0 ;` | Drift if field/period absent |
| SBend/RBend (no Enge) | `DIL L θ g/2 e₁ 0 e₂ 0 ;` | |
| SBend/RBend (with Enge) | `FC ..; DIL ..; FD ;` | |
| SBend (negative angle) | `CB ; DIL ..; CB ;` | |
| BPM/OTR/Corrector | _(skipped)_ | Zero-length |

## Supported particles

| Name | FOX command | Notes |
|------|------------|-------|
| `electron` | `RPE` | Preset |
| `proton` | `RPP` | Preset |
| `muon`, `mu-` | `RPMU` | μ⁻ preset |
| `mu+`, `antimuon` | `RP` | μ⁺ via generic RP (RPMU hardcodes charge −1) |
| `pion+`, `pi+` | `RPPI E 1` | Preset with charge |
| `pion-`, `pi-` | `RPPI E -1` | Preset with charge |
| `positron` | `RP` | Generic |
| `antiproton` | `RP` | Generic |
| `deuteron` | `RP` | Generic |
| `alpha` | `RP` | Generic (⁴He²⁺) |
| `carbon12` | `RP` | Generic (¹²C⁶⁺) |

## Known limitations

- **Drift fallback:** All active elements fall back to drifts when their field
  parameters are absent.
- **Enge coefficients:** Maximum 6 coefficients per fringe field (COSY INFINITY FC
  command limit). Longer coefficient lists are truncated with a warning.
- **Quad sign convention:** The current-based quadrupole path (`±G·I·r`) uses an
  electron-specific sign convention (QPF→negative B, QPD→positive B). Proton beamlines
  should use Bn1 (pole-tip field) instead, which bypasses the sign logic.

## Roadmap

pals2cosy implements provisional support for element types and parameters that
are present in COSY INFINITY but not yet part of the official PALS specification,
such as `MagneticMultipoleP.Bn2`–`Bn5`, `SolenoidP.Bz`, and `RFCavityP`. As a
documented PALS use case, these extensions may inform future PALS standardization.
If PALS adopts a different parameterization, pals2cosy can be adjusted accordingly.

Planned extensions:
- Combined-function magnets (dipole + multipole fields → COSY `MC`/`MS`)
- Electrostatic deflectors (`EC`, `ECL`, `ES`)
- Skew multipoles (`MMS`, `EMS`)
- Particle tracking (`SR` + `TR` + `PRAY`/`WRAY`)
- FIT optimization blocks

## Acknowledgments

This work was supported in part by DOE EPSCoR Award DE-SC0025583. The converter
is based in part on research performed on the University of Hawaiʻi free electron
laser by the author.

## Tests

```bash
pip install pytest
pytest tests/ -v
```
