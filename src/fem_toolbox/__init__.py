"""Compatibility import for Femantis, supported throughout v2.x."""

import sys
import warnings

import femantis
from femantis import __version__ as __version__
from femantis._legacy_imports import install_legacy_finder

warnings.warn(
    "fem_toolbox is deprecated; import femantis instead. "
    "Legacy imports remain supported throughout v2.x.",
    DeprecationWarning,
    stacklevel=2,
)
install_legacy_finder()
sys.modules[__name__] = femantis
