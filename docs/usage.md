# Usage

## Command line

```bash
# Print FOX code to stdout
pals2cosy beamline.yaml

# Write to a file
pals2cosy beamline.yaml -o output.fox

# Override kinetic energy and enable fringe fields
pals2cosy beamline.yaml -o output.fox --ke 40 --fr 1

# 2D phase space, 5th-order computation
pals2cosy beamline.yaml --dim 2 --order 5
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o FILE` | stdout | Output FOX file |
| `--ke FLOAT` | from lattice | Kinetic energy (MeV) |
| `--fr INT` | 0 | Fringe field order (0 = off, 1 = Enge) |
| `--order INT` | 3 | DA computation order |
| `--dim INT` | 3 | Phase-space dimensions (2 or 3) |
| `--quad-aperture FLOAT` | 0.027 | Quadrupole bore diameter (m) |
| `--particle STR` | from lattice | Particle type (`electron` or `proton`) |
| `--no-comments` | off | Omit element name comments |

## Python API

```python
from pals2cosy.parser import parse_lattice
from pals2cosy.converter import convert

beam_params, elements = parse_lattice("beamline.yaml")

fox_code = convert(
    beam_params, elements,
    ke_override=40,
    fringe_field_order=0,
    computation_order=3,
    dimensions=3,
)

with open("output.fox", "w") as f:
    f.write(fox_code)
```

## Running the generated FOX code

The output is a complete COSY INFINITY input file. To run it:

```bash
# Compile the COSY library (only needed once)
./cosy cosy.fox

# Run the generated lattice
./cosy output.fox
```

Results are written to `result.txt` as a JSON object containing Twiss
parameters ($\beta$, $\alpha$, $\gamma$, phase advance $\mu$) in both
horizontal and vertical planes.
