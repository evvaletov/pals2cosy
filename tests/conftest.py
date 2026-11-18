import os
import textwrap

import pytest


@pytest.fixture(scope="session")
def cosy_dir(tmp_path_factory):
    """A directory with COSY INFINITY compiled, shared by the COSY-backed tests."""
    from tests.cosy_harness import COSY_BIN, COSY_LIBRARY, compile_cosy
    if not (os.path.exists(COSY_BIN) and os.path.exists(COSY_LIBRARY)):
        pytest.skip("COSY INFINITY binary or cosy.fox not found")
    return compile_cosy(str(tmp_path_factory.mktemp("cosy")))


@pytest.fixture
def pals_file(tmp_path):
    """Write a PALS YAML document (given as text) and return its path."""
    counter = iter(range(1000))

    def write(text):
        path = tmp_path / f"lattice{next(counter)}.pals.yaml"
        path.write_text(textwrap.dedent(text))
        return str(path)

    return write
