"""pals2cosy — PALS Lattice to COSY INFINITY FOX converter."""

__version__ = "0.4.1"

from .converter import convert
from .parser import parse_lattice
from .pals_parser import parse_lattice as parse_pals_lattice
