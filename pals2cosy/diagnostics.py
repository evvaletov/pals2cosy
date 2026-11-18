"""Warning and error categories raised while reading and converting lattices.

Author: Eremey Valetov
"""


class PALSConversionWarning(UserWarning):
    """Something in the lattice is not represented, or only approximately,
    in the generated COSYScript. ``--strict`` turns these into errors."""


class PALSDeprecationWarning(FutureWarning):
    """Input uses a pre-0.6 pals2cosy convention that the PALS standard
    has since replaced."""


class PALSReminder(UserWarning):
    """A ``reminders`` entry of the PALS file, shown every time it is read."""


class UnsupportedFeatureError(NotImplementedError, ValueError):
    """The lattice uses a PALS construct that pals2cosy cannot convert."""
