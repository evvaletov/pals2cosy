# Examples

| File | Description |
|------|-------------|
| `fodo.pals.yaml` | FODO channel adapted from the [PALS project](https://github.com/pals-project/pals) example, with a `BeginningEle` reference particle. Demonstrates `facility:`, `line:`, `inherit:`, and `repeat:`. |
| `uhfel_excerpt.pals.yaml` | UH FEL beamline excerpt: quadrupoles given by field gradient `Bn1` [T/m] and a `Bend` given by `angle_ref` and length. |

Both files carry their reference particle in `ReferenceP` on the branch's
`BeginningEle`, so they convert without `--ke`/`--particle`.

## Usage

```bash
pals2cosy fodo.pals.yaml
pals2cosy uhfel_excerpt.pals.yaml --fr 1

# Override the reference particle from the command line
pals2cosy fodo.pals.yaml --ke 2000 --particle proton
```

The official PALS examples at the pinned specification commit are vendored
under `tests/pals_examples/` as a conformance suite.
