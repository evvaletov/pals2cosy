"""Command-line interface for pals2cosy.

Author: Eremey Valetov
"""

import argparse
import json
import os
import sys
import warnings

from ._io import load_file
from .converter import convert

_EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "examples")


def _detect_format(path):
    """Auto-detect input format from root key.

    Returns ``"official-pals"`` if root key is ``PALS:``,
    ``"felsim"`` if root key is ``beamline:``.
    """
    try:
        data = load_file(path)
    except (json.JSONDecodeError, ValueError) as e:
        raise ValueError(f"Cannot parse '{path}': {e}") from e
    except Exception as e:
        # Catch YAML parse errors (yaml.YAMLError) and other format-specific
        # exceptions without importing yaml at module level
        if type(e).__module__.startswith("yaml"):
            raise ValueError(f"Cannot parse '{path}': {e}") from e
        raise

    if not isinstance(data, dict):
        raise ValueError(
            f"Cannot auto-detect format for '{path}': "
            f"expected a mapping at top level, got {type(data).__name__}"
        )

    if "PALS" in data:
        return "official-pals"
    if "beamline" in data:
        return "felsim"

    raise ValueError(
        f"Cannot auto-detect format for '{path}'. "
        "Expected root key 'PALS:' (official PALS) or 'beamline:' (flat format)."
    )


def main():
    parser = argparse.ArgumentParser(
        prog="pals2cosy",
        description="Convert a PALS lattice to COSY INFINITY input code.",
    )
    parser.add_argument("input", nargs="?", help="Input lattice file (JSON or YAML)")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    parser.add_argument("--fr", type=int, default=0, help="Fringe field order (default: 0)")
    parser.add_argument("--order", type=int, default=3, help="DA computation order (default: 3)")
    parser.add_argument("--dim", type=int, default=3, choices=[2, 3],
                        help="Phase-space dimensions (default: 3)")
    parser.add_argument("--ke", type=float, default=None,
                        help="Override kinetic energy (MeV)")
    parser.add_argument("--particle", type=str, default=None,
                        help="Override particle type. Presets: electron, proton, muon/mu-, "
                             "pion+/pi+, pion-/pi-. Generic RP: mu+/antimuon, positron, "
                             "antiproton, deuteron, alpha, carbon12")
    parser.add_argument("--quad-aperture", type=float, default=0.027,
                        help="Quadrupole bore diameter in meters (default: 0.027)")
    parser.add_argument("--no-comments", action="store_true",
                        help="Omit element name comments")
    parser.add_argument("--twiss", action="store_true",
                        help="Append GT Twiss extraction and JSON output (requires periodic lattice)")
    parser.add_argument("--mode", type=str, default="auto",
                        choices=["auto", "felsim", "felsim-strict", "pals", "official-pals"],
                        help="Input format override (default: auto)")
    parser.add_argument("--beamline", type=str, default=None,
                        help="BeamLine name to expand (official PALS only; default: last)")
    parser.add_argument("--examples-dir", action="store_true",
                        help="Print path to bundled example files and exit")
    parser.add_argument("--demo", action="store_true",
                        help="Convert the bundled FODO example to stdout")

    args = parser.parse_args()

    if args.examples_dir:
        print(_EXAMPLES_DIR)
        return

    if args.demo:
        fodo = os.path.join(_EXAMPLES_DIR, "fodo.pals.yaml")
        with open(fodo) as f:
            yaml_text = f.read()
        from .pals_parser import parse_lattice as parse_pals
        beam_params, elements = parse_pals(fodo, ke_override=1000,
                                           particle_override="proton")
        fox = convert(beam_params, elements, ke_override=1000,
                      particle_override="proton", twiss=True, dimensions=2)
        print("# === Input: fodo.pals.yaml ===\n")
        print(yaml_text)
        print("# === Output: COSY INFINITY ===\n")
        sys.stdout.write(fox)
        return

    if args.input is None:
        parser.error("the following arguments are required: input (or use --demo / --examples-dir)")

    mode = args.mode
    if mode == "auto":
        mode = _detect_format(args.input)
    if mode == "pals":
        warnings.warn(
            "--mode pals is deprecated; use --mode felsim-strict instead",
            DeprecationWarning,
            stacklevel=2,
        )
    if mode == "felsim-strict":
        mode = "pals"  # internal name for strict mode

    if mode == "official-pals":
        from .pals_parser import parse_lattice as parse_pals
        beam_params, elements = parse_pals(
            args.input,
            beamline_name=args.beamline,
            ke_override=args.ke,
            particle_override=args.particle,
        )
    else:
        if args.beamline is not None:
            warnings.warn(
                f"--beamline is only used with official-pals mode; "
                f"ignored for mode '{mode}'"
            )
        from .parser import parse_lattice
        beam_params, elements = parse_lattice(args.input, mode=mode)

    fox = convert(
        beam_params, elements,
        fringe_field_order=args.fr,
        computation_order=args.order,
        dimensions=args.dim,
        quad_aperture=args.quad_aperture,
        ke_override=args.ke,
        particle_override=args.particle,
        comments=not args.no_comments,
        twiss=args.twiss,
    )

    if args.output:
        with open(args.output, "w") as f:
            f.write(fox)
    else:
        sys.stdout.write(fox)
