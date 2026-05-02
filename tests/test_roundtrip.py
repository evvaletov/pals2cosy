"""Round-trip validation: PALS YAML → COSYScript → COSY → Twiss comparison.

Uses the FODO cell example (which is periodic, so GT can find a fixed point)
to validate the full pals2cosy → COSY pipeline.

The UH FEL beamline is a one-pass transport line, not a ring. COSY's GT
command requires a periodic system to find matched Twiss parameters. The
FODO cell is the correct test case for GT validation.

Author: Eremey Valetov
"""

import json
import os
import subprocess
import tempfile

import pytest

from pals2cosy.pals_parser import parse_lattice as parse_pals
from pals2cosy.converter import convert

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")
FODO_YAML = os.path.join(EXAMPLES_DIR, "fodo.pals.yaml")
COSY_BIN = "/usr/local/bin/cosy"
COSY_FOX = "/usr/local/bin/cosy.fox"


def _run_cosy(fox_code, timeout=120):
    """Run COSY INFINITY on COSYScript, return parsed result.txt."""
    with tempfile.TemporaryDirectory(prefix="pals2cosy_test_") as tmpdir:
        fox_path = os.path.join(tmpdir, "input.fox")
        with open(fox_path, "w") as f:
            f.write(fox_code)

        import shutil
        shutil.copy2(COSY_BIN, os.path.join(tmpdir, "cosy"))
        os.chmod(os.path.join(tmpdir, "cosy"), 0o755)

        shutil.copy2(COSY_FOX, os.path.join(tmpdir, "cosy.fox"))
        subprocess.run(["./cosy", "cosy.fox"], cwd=tmpdir,
                       stdin=subprocess.DEVNULL, capture_output=True,
                       timeout=timeout)

        result = subprocess.run(
            ["./cosy", "input.fox"],
            cwd=tmpdir,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            timeout=timeout,
        )

        result_path = os.path.join(tmpdir, "result.txt")
        if not os.path.exists(result_path):
            stdout = result.stdout.decode()
            pytest.fail(f"COSY did not produce result.txt\nstdout: {stdout[-500:]}")

        with open(result_path) as f:
            text = f.read()

        return json.loads(text)


@pytest.mark.skipif(not os.path.exists(COSY_BIN),
                    reason="COSY INFINITY binary not found")
@pytest.mark.skipif(not os.path.exists(COSY_FOX),
                    reason="COSY library (cosy.fox) not found")
class TestFODORoundTrip:
    """FODO cell round-trip: periodic lattice with GT Twiss extraction."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Generate COSYScript from FODO PALS example.

        Uses quad_aperture=1.0 m (bore diameter) to give a stable cell.
        With Bn1=1.0 T and R=0.5 m: G=2 T/m, f≈2.83 m at 1 GeV proton.
        """
        bp, elems = parse_pals(FODO_YAML, ke_override=1000,
                               particle_override="proton")
        self.fox = convert(bp, elems, ke_override=1000,
                           particle_override="proton",
                           quad_aperture=1.0,
                           dimensions=2,  # 4D only; longitudinal degenerate without RF
                           fringe_field_order=0, computation_order=3,
                           twiss=True)
        self.result = _run_cosy(self.fox)

    def test_gt_succeeds(self):
        """GT should find a fixed point for the periodic FODO cell."""
        assert "twiss" in self.result

    def test_beta_x_positive(self):
        beta_x = float(self.result["twiss"]["beta_x"])
        assert beta_x > 0

    def test_beta_y_positive(self):
        beta_y = float(self.result["twiss"]["beta_y"])
        assert beta_y > 0

    def test_beta_x_physical(self):
        """Beta function should be O(cell length) for a FODO cell."""
        beta_x = float(self.result["twiss"]["beta_x"])
        # FODO cell = 3.0 m, beta should be roughly 1–10 m
        assert 0.5 < beta_x < 50

    def test_beta_y_physical(self):
        beta_y = float(self.result["twiss"]["beta_y"])
        assert 0.5 < beta_y < 50
