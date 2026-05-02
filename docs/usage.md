# Usage

## Command line

```bash
# Auto-detect format and print COSYScript code to stdout
pals2cosy beamline.pals.yaml --ke 1000 --particle proton

# Write to file with fringe fields
pals2cosy lattice.pals.yaml -o output.fox --ke 40 --fr 1

# Select a specific BeamLine from the facility
pals2cosy lattice.pals.yaml --beamline fodo_cell --ke 1000

# Compute Twiss parameters (periodic lattices only)
pals2cosy lattice.pals.yaml --ke 1000 --particle proton --twiss

# 2D phase space, 5th-order computation
pals2cosy beamline.pals.yaml --dim 2 --order 5 --ke 40
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o FILE` | stdout | Output COSYScript file |
| `--ke FLOAT` | from lattice | Kinetic energy (MeV) |
| `--fr INT` | 0 | Fringe field order (0=none, 1=Enge, 2=symplectic, 3=high-precision Enge) |
| `--order INT` | 3 | DA computation order |
| `--dim INT` | 3 | Phase-space dimensions (2 or 3) |
| `--quad-aperture FLOAT` | 0.027 | Quadrupole bore diameter (m) |
| `--particle STR` | from lattice | Particle type (see README for full list) |
| `--twiss` | off | Append GT Twiss extraction (requires periodic lattice) |
| `--mode STR` | auto | Input format override |
| `--beamline NAME` | last | BeamLine to expand (official PALS only) |
| `--no-comments` | off | Omit element name comments |

## Python API

```python
from pals2cosy import convert, parse_pals_lattice

beam_params, elements = parse_pals_lattice(
    "lattice.pals.yaml",
    ke_override=1000,
    particle_override="proton",
    beamline_name="fodo_cell",
)

fox_code = convert(beam_params, elements, ke_override=1000,
                   particle_override="proton")

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
$\mu$) in both horizontal and vertical planes.
