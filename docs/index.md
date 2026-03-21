# pals2cosy

**pals2cosy** converts [PALS](https://pals-project.readthedocs.io/) lattice
files into [COSY INFINITY](https://cosyinfinity.org) FOX input code.
The generated FOX program computes the transfer map for the beamline described
in the input file.

The official PALS format (`PALS:` root) is supported with `facility:` definitions,
`BeamLine` composition (`line:`, `inherit:`, `repeat:`), `MagneticMultipoleP.Bn1`
for quadrupoles, and `BendP` for dipoles. A flat `beamline:` format with positioned
elements is also accepted.

```{toctree}
:maxdepth: 2

installation
usage
input-format
element-mapping
architecture
```
