# pals2cosy

**pals2cosy** converts [PALS](https://pals-project.readthedocs.io/) v2 lattice
files (JSON or YAML) into [COSY INFINITY](https://bt.pa.msu.edu/index_cosy.htm)
FOX input code. The generated FOX program computes the transfer map and Twiss
parameters for the beamline described in the input file.

The converter handles both PALS-native dipole representations (SBend/RBend with
edge-angle attributes) and FELsim's DIPOLE_WEDGE triplet convention, producing
identical FOX output from either.

```{toctree}
:maxdepth: 2

installation
usage
input-format
element-mapping
architecture
```
