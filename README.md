# pals2cosy

Convert PALS v2 lattice files (JSON/YAML) to COSY INFINITY FOX code.

## Installation

```bash
pip install .
```

Or for development:

```bash
pip install -e .
```

Requires Python 3.9+ and PyYAML (for YAML input).

## Usage

```bash
# Convert to stdout
pals2cosy beamline.yaml

# Write to file with custom settings
pals2cosy beamline.yaml -o output.fox --ke 40 --order 3 --dim 3

# With fringe fields enabled
pals2cosy beamline.yaml --fr 1
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o FILE` | stdout | Output FOX file |
| `--ke FLOAT` | from lattice | Kinetic energy (MeV) |
| `--fr INT` | 0 | Fringe field order |
| `--order INT` | 3 | DA computation order |
| `--dim INT` | 3 | Phase-space dimensions (2 or 3) |
| `--quad-aperture FLOAT` | 0.027 | Quadrupole bore diameter (m) |
| `--particle STR` | from lattice | Particle type (electron/proton) |
| `--no-comments` | off | Omit element name comments |

## Input format

PALS v2 YAML or JSON. See the [PALS specification](https://pals-project.readthedocs.io/)
and `examples/uhfel_beamline.yaml` for the supported format.

### Dipole representation

Two dipole representations are supported and produce identical FOX output:

**PALS-native** — edge angles as SBend attributes:
```yaml
- type: SBend
  parameters:
    bending_angle_deg: 11.25
    dipole_length_m: 0.037389
    pole_gap_m: 0.0127
    entrance_edge_angle_deg: 0.0
    exit_edge_angle_deg: 11.25
  fringe_fields:
    enge_coefficients: [56.49, -50.79, 19.32, -3.621, 0.3315, -0.01193]
```

**FELsim DIPOLE_WEDGE** — edge kicks as separate elements:
```yaml
- type: DIPOLE_WEDGE
  parameters: {wedge_angle_deg: 0.0, ...}
  fringe_fields:
    enge_coefficients: [56.49, -50.79, 19.32, -3.621, 0.3315, -0.01193]
- type: SBend
  parameters: {bending_angle_deg: 11.25, ...}
- type: DIPOLE_WEDGE
  parameters: {wedge_angle_deg: 11.25, ...}
```

DPW-DPH-DPW triplets are auto-detected and consolidated.

## Output

The generated FOX code computes the transfer map and Twiss parameters:

```fox
INCLUDE 'COSY' ;
PROCEDURE RUN ;
    ...
    PROCEDURE LATTICE ;
        UM ; CR ;
        DL 0.358775 ;
        MQ 0.0889 -0.0322 0.0135 ;
        ...
    ENDPROCEDURE ;

    OV 3 3 0 ;
    RPE 40 ;
    FR 0 ;
    LATTICE ;
    CO 1 ; PM 99 ;
    GT MAP F0 MU0 A0 B0 G0 R0 ;
    ...
ENDPROCEDURE ;
RUN ;
END ;
```

### Element mapping

| Input type | FOX command | Notes |
|-----------|-------------|-------|
| Quadrupole (focusing) | `MQ L b r ;` | b = -G·I·r |
| Quadrupole (defocusing) | `MQ L b r ;` | b = +G·I·r |
| SBend/RBend (no Enge) | `DIL L θ g/2 e₁ 0 e₂ 0 ;` | |
| SBend/RBend (with Enge) | `FC ..; DIL ..; FD ;` | |
| SBend (negative angle) | `CB ; DIL ..; CB ;` | |
| Wiggler/Undulator | `DL L ;` | Passive drift |
| BPM/OTR/Corrector | _(skipped)_ | Zero-length |

## Tests

```bash
pip install pytest
pytest tests/ -v
```
