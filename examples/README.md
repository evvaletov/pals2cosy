# Examples

| File | Format | Description |
|------|--------|-------------|
| `fodo.pals.yaml` | Official PALS | Standard FODO cell from the [PALS project](https://github.com/pals-project/pals). Demonstrates `facility:`, `line:`, `inherit:`, and `repeat:`. |
| `uhfel_excerpt.pals.yaml` | Official PALS | Beamline excerpt with `MagneticMultipoleP.Bn1` quadrupoles and `BendP` dipole. |

## Usage

```bash
# FODO cell (requires --ke since no Lattice.particle.kinetic_energy)
pals2cosy fodo.pals.yaml --ke 1000 --particle proton

# Beamline excerpt
pals2cosy uhfel_excerpt.pals.yaml --ke 40 --particle electron
```
