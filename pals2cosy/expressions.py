"""Safe evaluation of PALS arithmetic expressions.

Supports numbers, ``+ - * /``, ``^`` and ``**`` (right associative), unary
signs, parentheses, quoted species strings, ``Inf``, the PALS predefined
constants, lattice constants and variables, element parameter references
(``Q1>MagneticMultipoleP.Kn1``) and the function set of the PALS standard.
No Python ``eval`` is involved: expressions are tokenized and evaluated by a
small recursive-descent parser.

Author: Eremey Valetov
"""

import math
import re
import warnings

from .constants import PALS_CONSTANTS
from .diagnostics import PALSConversionWarning
from .species import anomalous_moment_of, charge_of, mass_of, parse_species

_TOKEN_RE = re.compile(r"""
    (?P<ws>\s+)
  | (?P<num>(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)
  | (?P<ref>[A-Za-z_]\w*>(?!>)[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*)
  | (?P<name>[A-Za-z_]\w*)
  | (?P<str>"[^"]*"|'[^']*')
  | (?P<op>\*\*|[-+*/^(),])
""", re.VERBOSE)


class ExpressionError(ValueError):
    pass


def _tokenize(text):
    pos = 0
    tokens = []
    while pos < len(text):
        m = _TOKEN_RE.match(text, pos)
        if m is None:
            raise ExpressionError(
                f"unexpected character {text[pos]!r} at position {pos} in {text!r}")
        pos = m.end()
        kind = m.lastgroup
        if kind == "ws":
            continue
        value = m.group()
        if kind == "op" and value == "**":
            value = "^"
        elif kind == "str":
            value = value[1:-1]
        tokens.append((kind, value))
    return tokens


def _sign(x):
    return float((x > 0) - (x < 0))


def _nint(x):
    n = math.floor(abs(x) + 0.5)
    return float(n if x >= 0 else -n)


def _modulo(a, p):
    return a - math.floor(a / p) * p


def _sinc(x):
    return 1.0 if x == 0 else math.sin(x) / x


_FUNCTIONS = {
    "sqrt": math.sqrt, "log": math.log, "exp": math.exp,
    "sin": math.sin, "cos": math.cos, "tan": math.tan,
    "cot": lambda x: 1 / math.tan(x), "sinc": _sinc,
    "asin": math.asin, "acos": math.acos, "atan": math.atan, "atan2": math.atan2,
    "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
    "coth": lambda x: 1 / math.tanh(x),
    "asinh": math.asinh, "acosh": math.acosh, "atanh": math.atanh,
    "acoth": lambda x: math.atanh(1 / x),
    "abs": abs, "factorial": lambda n: float(math.factorial(int(n))),
    "int": lambda x: float(math.trunc(x)), "nint": _nint, "sign": _sign,
    "floor": lambda x: float(math.floor(x)), "ceiling": lambda x: float(math.ceil(x)),
    "modulo": _modulo,
}

_SPECIES_FUNCTIONS = {
    "mass_of": mass_of,
    "charge_of": charge_of,
    "anomalous_moment_of": anomalous_moment_of,
}


class _Parser:
    def __init__(self, tokens, scope, text):
        self.tokens = tokens
        self.i = 0
        self.scope = scope
        self.text = text

    def peek(self):
        return self.tokens[self.i] if self.i < len(self.tokens) else (None, None)

    def take(self, value=None):
        tok = self.peek()
        if tok[0] is None or (value is not None and tok[1] != value):
            want = repr(value) if value else "a value"
            raise ExpressionError(f"expected {want} in {self.text!r}")
        self.i += 1
        return tok

    def parse(self):
        value = self.expr()
        if self.i != len(self.tokens):
            raise ExpressionError(
                f"unexpected {self.tokens[self.i][1]!r} in {self.text!r}")
        return value

    def expr(self):
        value = self.term()
        while self.peek()[1] in ("+", "-"):
            op = self.take()[1]
            rhs = self.term()
            value = _num(value, self.text) + _num(rhs, self.text) * (1 if op == "+" else -1)
        return value

    def term(self):
        value = self.unary()
        while self.peek()[1] in ("*", "/"):
            op = self.take()[1]
            rhs = _num(self.unary(), self.text)
            lhs = _num(value, self.text)
            if op == "*":
                value = lhs * rhs
            elif rhs == 0:
                raise ExpressionError(f"division by zero in {self.text!r}")
            else:
                value = lhs / rhs
        return value

    def unary(self):
        if self.peek()[1] in ("+", "-"):
            op = self.take()[1]
            value = _num(self.unary(), self.text)
            return -value if op == "-" else value
        return self.power()

    def power(self):
        base = self.atom()
        if self.peek()[1] == "^":
            self.take()
            exponent = _num(self.unary(), self.text)
            try:
                value = _num(base, self.text) ** exponent
            except (OverflowError, ZeroDivisionError) as exc:
                raise ExpressionError(f"power failed in {self.text!r}: {exc}") from exc
            if isinstance(value, complex):
                raise ExpressionError(
                    f"negative number raised to a fractional power in {self.text!r}")
            return value
        return base

    def atom(self):
        kind, value = self.take()
        if kind == "num":
            return float(value)
        if kind == "str":
            return value
        if kind == "ref":
            name, path = value.split(">", 1)
            return self.scope.lookup_ref(name, path)
        if kind == "name":
            if self.peek()[1] == "(":
                return self.call(value)
            return self.scope.lookup(value)
        if value == "(":
            inner = self.expr()
            self.take(")")
            return inner
        raise ExpressionError(f"unexpected {value!r} in {self.text!r}")

    def call(self, fname):
        self.take("(")
        args = []
        if self.peek()[1] != ")":
            args.append(self.expr())
            while self.peek()[1] == ",":
                self.take()
                args.append(self.expr())
        self.take(")")
        if fname == "expr":
            if len(args) != 1:
                raise ExpressionError(f"expr() takes one argument in {self.text!r}")
            return args[0]
        if fname in ("random", "random_gauss"):
            return self.scope.random(fname, args)
        if fname in _SPECIES_FUNCTIONS:
            if len(args) != 1 or not isinstance(args[0], str):
                raise ExpressionError(
                    f"{fname}() takes one species name in {self.text!r}")
            return _SPECIES_FUNCTIONS[fname](args[0])
        func = _FUNCTIONS.get(fname)
        if func is None:
            raise ExpressionError(f"unknown function {fname}() in {self.text!r}")
        try:
            return float(func(*(_num(a, self.text) for a in args)))
        except (TypeError, ValueError, OverflowError, ZeroDivisionError) as exc:
            raise ExpressionError(f"{fname}() failed in {self.text!r}: {exc}") from exc


def _num(value, text):
    if isinstance(value, str):
        raise ExpressionError(f"string {value!r} used as a number in {text!r}")
    return value


class Scope:
    """Name resolution for expressions: PALS constants, lattice constants and
    variables, and element parameter references.

    ``definitions`` maps names to raw values (numbers, species strings, or
    expression strings), evaluated lazily. ``ref_resolver(name, path)``
    returns the raw value of an element parameter, or raises ``KeyError``.
    """

    def __init__(self, definitions=None, ref_resolver=None, seed=0):
        self.definitions = dict(definitions or {})
        self.ref_resolver = ref_resolver
        self._values = {}
        self._active = []
        self._rng = None
        self._seed = seed

    def evaluate(self, raw):
        if isinstance(raw, bool):
            raise ExpressionError(f"logical {raw!r} used as a number")
        if isinstance(raw, (int, float)):
            return float(raw)
        if not isinstance(raw, str):
            raise ExpressionError(f"cannot evaluate {raw!r}")
        return _Parser(_tokenize(raw), self, raw).parse()

    def evaluate_definition(self, raw):
        """Value of a constant or variable definition: a number or a species name
        such as "#3He", "Au+79" or "electron"."""
        try:
            return self.evaluate(raw)
        except ExpressionError:
            if not _is_species_name(raw):
                raise
            return raw.strip()

    def lookup(self, name):
        if name in self.definitions:
            if name in self._values:
                return self._values[name]
            if name in self._active:
                chain = " -> ".join(self._active + [name])
                raise ExpressionError(f"circular definition: {chain}")
            self._active.append(name)
            try:
                value = self.evaluate_definition(self.definitions[name])
            finally:
                self._active.pop()
            self._values[name] = value
            return value
        if name in PALS_CONSTANTS:
            return PALS_CONSTANTS[name]
        if name == "Inf":
            return math.inf
        if name in ("PARAMETER", "SELF"):
            raise ExpressionError(f"'{name}' is only meaningful inside a set command")
        raise ExpressionError(f"undefined name '{name}'")

    def lookup_ref(self, name, path):
        if self.ref_resolver is None:
            raise ExpressionError(f"cannot resolve parameter reference {name}>{path}")
        key = f"{name}>{path}"
        if key in self._active:
            raise ExpressionError(f"circular definition through {key}")
        try:
            raw = self.ref_resolver(name, path)
        except KeyError as exc:
            raise ExpressionError(
                f"cannot resolve parameter reference {key}: {exc.args[0]}") from None
        self._active.append(key)
        try:
            return self.evaluate(raw)
        finally:
            self._active.pop()

    def random(self, fname, args):
        import random
        if self._rng is None:
            self._rng = random.Random(self._seed)
            warnings.warn(
                f"{fname}() in the lattice is evaluated with a fixed-seed generator; "
                "the COSYScript output is one realization of the random values",
                PALSConversionWarning, stacklevel=2,
            )
        if fname == "random":
            return self._rng.random()
        cut = args[0] if args else None
        while True:
            x = self._rng.gauss(0.0, 1.0)
            if cut is None or abs(x) <= cut:
                return x


def _is_species_name(raw):
    if not isinstance(raw, str):
        return False
    try:
        parse_species(raw)
    except ValueError:
        return False
    return True
