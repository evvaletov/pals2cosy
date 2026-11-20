# Usage

## Command line

```bash
# Auto-detect format and print COSYScript code to stdout
pals2cosy beamline.pals.yaml

# Override the reference particle given in the file
pals2cosy beamline.pals.yaml --ke 1000 --particle proton

# Treat anything not represented in COSYScript as an error
pals2cosy beamline.pals.yaml --strict

# Write to file with fringe fields
pals2cosy lattice.pals.yaml -o output.fox --ke 40 --fr 1

# Select a specific BeamLine from the facility (one without a BeginningEle
# needs the reference particle on the command line)
pals2cosy lattice.pals.yaml --beamline fodo_cell --ke 1000 --particle proton

# Compute Twiss parameters (periodic lattices only)
pals2cosy lattice.pals.yaml --ke 1000 --particle proton --twiss

# 2D phase space, 5th-order computation
pals2cosy beamline.pals.yaml --dim 2 --order 5 --ke 40
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o FILE` | stdout | Output COSYScript file |
| `--ke FLOAT` | from lattice | Reference kinetic energy (MeV) |
| `--fr INT` | 0 | Fringe field order (0=none, 1=Enge, 2=symplectic, 3=high-precision Enge) |
| `--order INT` | 3 | DA computation order |
| `--dim INT` | 3 (2 with `--twiss`) | Phase-space dimensions (2 or 3) |
| `--quad-aperture FLOAT` | 0.027 | Multipole bore diameter (m); pole-tip fields are given at its radius, also the half-gap of bends without `edge1_int`/`edge2_int` or `ApertureP` under `--fr 1`–`3` |
| `--solenoid-fringe FLOAT` | bore radius | Fringe-field length (m) of COSY's soft-edge solenoid `CMS` |
| `--particle STR` | from lattice | Reference species (PALS/OpenPMD name or ion such as `#12C+6`) |
| `--twiss` | off | Append GT Twiss extraction (requires periodic lattice; a 3D map also needs RF) |
| `--mode STR` | auto | Input format override |
| `--beamline NAME` | first branch | Lattice branch or BeamLine to expand (official PALS only) |
| `--strict` | off | Turn conversion warnings into errors |
| `--no-comments` | off | Omit element name comments |
| `--version` | | Print the version and the targeted PALS commit |

A lattice that cannot be converted exits with status 1 and a one-line
message. Conversion warnings (`PALSConversionWarning`) name each element or
construct that is not represented, or only approximately, in the output.

## Python API

```python
from pals2cosy import convert, parse_pals_lattice

beam_params, elements = parse_pals_lattice(
    "lattice.pals.yaml",
    beamline_name="fodo_cell",   # optional
    ke_override=1000,            # MeV; needed when the line has no BeginningEle
    particle_override="proton",
)

fox_code = convert(beam_params, elements, fringe_field_order=1)

with open("output.fox", "w") as f:
    f.write(fox_code)
```

## Running the generated COSYScript code

The output is a complete COSY INFINITY input file. To run it:

```bash
# Compile the COSY library (only needed once)
./cosy cosy.fox

# Run the generated lattice
./cosy output.fox
```

If `--twiss` was specified, results are written to `result.txt` as a JSON
object containing Twiss parameters ($\beta$, $\alpha$, $\gamma$, phase advance
$\mu$) in both horizontal and vertical planes. The map is 2D unless `--dim 3`
is given; COSY's `GT` needs RF to handle the longitudinal plane of a 3D map.

Pass reference-particle overrides to the parser rather than to `convert()`:
normalized strengths (`KnN`, `Ksol`, `Bn0_ref`) are turned into fields with
the parser's reference rigidity, and `convert()` warns when its overrides
would make them inconsistent.
