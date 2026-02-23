# Installation

## Requirements

- Python 3.9+
- [PyYAML](https://pyyaml.org/) (for YAML input; JSON works without it)

## Install from source

```bash
git clone https://github.com/evvaletov/pals2cosy.git
cd pals2cosy
pip install .
```

For development (editable install):

```bash
pip install -e .
```

## Verify

```bash
pals2cosy --help
```

## Running tests

```bash
pip install pytest
pytest tests/ -v
```
