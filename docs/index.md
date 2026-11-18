# pals2cosy

**pals2cosy** converts [PALS](https://pals-project.readthedocs.io/) lattice
files into [COSY INFINITY](https://cosyinfinity.org) input code.
The generated program computes the transfer map for the beamline described
in the input file.

pals2cosy 0.6 reads the official PALS format (`PALS:` root) as specified at
PALS commit [`b9d5b355`](https://github.com/pals-project/pals/tree/b9d5b35555a63f7ba015460e1c17ac929d35e215):
`facility:` definitions, `Lattice` branches and `use`, `BeamLine` composition
(`line:`, `inherit:`, `repeat:`), constants, variables and expressions, the
reference particle of the branch's `BeginningEle`, and the magnet, bend,
solenoid and RF parameter groups. A flat `beamline:` format with positioned
elements is also accepted.

```{toctree}
:maxdepth: 2

installation
usage
input-format
element-mapping
architecture
```
