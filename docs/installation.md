# Installation

## Requirements

- Python 3.9+
- [PyYAML](https://pyyaml.org/) (for YAML input; JSON works without it)

## Install from PyPI

```bash
pip install pals2cosy
```

## Install from source

```bash
git clone https://github.com/evvaletov/pals2cosy.git
cd pals2cosy
pip install .
```

For development (editable install):

```bash
pip install -e ".[dev]"
```

## Verify

```bash
pals2cosy --version   # pals2cosy and the targeted PALS specification commit
pals2cosy --demo
```

## Running tests

```bash
pytest tests/ -v
```

Tests that run COSY INFINITY (transfer-map checks and the conformance suite)
are skipped when `/usr/local/bin/cosy` and `cosy.fox` are not present; set
`PALS2COSY_COSY_BIN` and `PALS2COSY_COSY_LIBRARY` to point elsewhere.
