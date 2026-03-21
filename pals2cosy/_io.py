"""Shared file loading utilities for JSON, YAML, and TOML inputs.

Author: Eremey Valetov
"""

import json
import os


def load_file(path):
    """Load a JSON, YAML, or TOML file, auto-detecting from extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".toml":
        try:
            import tomllib
        except ModuleNotFoundError:
            import tomli as tomllib
        with open(path, "rb") as f:
            return tomllib.load(f)
    with open(path, "r") as f:
        if ext in (".yaml", ".yml"):
            import yaml
            return yaml.safe_load(f)
        return json.load(f)
