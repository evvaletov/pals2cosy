"""Command-line interface for pals2cosy.

Author: Eremey Valetov
"""

import argparse
import json
import os
import sys
import warnings

from ._io import load_file
from ._version import PALS_COMMIT_FULL, __version__
from .converter import convert
from .diagnostics import PALSConversionWarning

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


def _format_warning(message, category, filename, lineno, line=None):
    return f"pals2cosy: {category.__name__}: {message}\n"


def main(argv=None):
    warnings.formatwarning = _format_warning
    parser = argparse.ArgumentParser(
        prog="pals2cosy",
        description="Convert a PALS lattice to COSY INFINITY input code.",
    )
    parser.add_argument("--version", action="version",
                        version=f"pals2cosy {__version__} "
                                f"(PALS specification {PALS_COMMIT_FULL})")
    parser.add_argument("input", nargs="?", help="Input lattice file (JSON or YAML)")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    parser.add_argument("--fr", type=int, default=0, help="Fringe field order (default: 0)")
    parser.add_argument("--order", type=int, default=3, help="DA computation order (default: 3)")
    parser.add_argument("--dim", type=int, default=None, choices=[2, 3],
                        help="Phase-space dimensions (default: 3, or 2 with --twiss)")
    parser.add_argument("--ke", type=float, default=None,
                        help="Override the reference kinetic energy (MeV)")
    parser.add_argument("--particle", type=str, default=None,
                        help="Override the reference species: an OpenPMD/PALS name "
                             "(electron, positron, proton, anti-proton, muon, anti-muon, "
                             "pion+, pion-, deuteron, alpha, ...) or an ion such as "
                             "'#12C+6' or 'Au+79'")
    parser.add_argument("--quad-aperture", type=float, default=0.027,
                        help="Multipole bore diameter in meters; pole-tip fields are "
                             "given at its radius, which is also the half-gap of bends "
                             "without edge1_int/edge2_int or ApertureP under --fr 1-3 "
                             "(default: 0.027)")
    parser.add_argument("--solenoid-fringe", type=float, default=None,
                        help="Fringe-field length in meters of COSY's soft-edge "
                             "solenoid; shorter approaches the hard-edge PALS solenoid "
                             "but integrates more slowly (default: bore radius)")
    parser.add_argument("--no-comments", action="store_true",
                        help="Omit element name comments")
    parser.add_argument("--twiss", action="store_true",
                        help="Append GT Twiss extraction and JSON output (requires periodic lattice)")
    parser.add_argument("--mode", type=str, default="auto",
                        choices=["auto", "felsim", "felsim-strict", "pals", "official-pals"],
                        help="Input format override (default: auto)")
    parser.add_argument("--beamline", type=str, default=None,
                        help="Lattice branch or BeamLine to expand (official PALS only; "
                             "default: first branch of the used Lattice)")
    parser.add_argument("--strict", action="store_true",
                        help="Treat conversion warnings (features not represented in "
                             "COSYScript) as errors")
    parser.add_argument("--examples-dir", action="store_true",
                        help="Print path to bundled example files and exit")
    parser.add_argument("--demo", action="store_true",
                        help="Convert the bundled FODO example to stdout")

    args = parser.parse_args(argv)

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

    with warnings.catch_warnings():
        if args.strict:
            warnings.simplefilter("error", PALSConversionWarning)
        try:
            fox = _convert_file(args)
        except (ValueError, ArithmeticError, OSError, RecursionError,
                PALSConversionWarning) as exc:
            print(f"pals2cosy: error: {exc}", file=sys.stderr)
            return 1

    if args.output:
        with open(args.output, "w") as f:
            f.write(fox)
    else:
        sys.stdout.write(fox)


def _convert_file(args):
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
                f"ignored for mode '{mode}'",
                PALSConversionWarning,
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
        solenoid_fringe=args.solenoid_fringe,
    )

    return fox
