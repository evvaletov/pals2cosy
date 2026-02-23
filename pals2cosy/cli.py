"""Command-line interface for pals2cosy.

Author: Eremey Valetov
"""

import argparse
import sys

from .parser import parse_lattice
from .converter import convert


def main():
    parser = argparse.ArgumentParser(
        prog="pals2cosy",
        description="Convert a PALS v2 lattice (JSON/YAML) to COSY INFINITY FOX code.",
    )
    parser.add_argument("input", help="Input lattice file (JSON or YAML)")
    parser.add_argument("-o", "--output", help="Output FOX file (default: stdout)")
    parser.add_argument("--fr", type=int, default=0, help="Fringe field order (default: 0)")
    parser.add_argument("--order", type=int, default=3, help="DA computation order (default: 3)")
    parser.add_argument("--dim", type=int, default=3, choices=[2, 3],
                        help="Phase-space dimensions (default: 3)")
    parser.add_argument("--ke", type=float, default=None,
                        help="Override kinetic energy (MeV)")
    parser.add_argument("--particle", type=str, default=None,
                        choices=["electron", "proton"],
                        help="Override particle type")
    parser.add_argument("--quad-aperture", type=float, default=0.027,
                        help="Quadrupole bore diameter in meters (default: 0.027)")
    parser.add_argument("--no-comments", action="store_true",
                        help="Omit element name comments")

    args = parser.parse_args()

    beam_params, elements = parse_lattice(args.input)

    fox = convert(
        beam_params, elements,
        fringe_field_order=args.fr,
        computation_order=args.order,
        dimensions=args.dim,
        quad_aperture=args.quad_aperture,
        ke_override=args.ke,
        particle_override=args.particle,
        comments=not args.no_comments,
    )

    if args.output:
        with open(args.output, "w") as f:
            f.write(fox)
    else:
        sys.stdout.write(fox)
