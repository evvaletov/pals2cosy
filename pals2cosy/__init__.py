"""pals2cosy — PALS Lattice to COSY INFINITY (COSYScript) converter."""

__version__ = "0.5.0"

from .converter import convert
from .parser import parse_lattice
from .pals_parser import parse_lattice as parse_pals_lattice
