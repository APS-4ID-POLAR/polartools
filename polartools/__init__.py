"""
Copyright (c) 2020, UChicago Argonne, LLC.

See LICENSE file for details.
"""

try:
    from ._version import __version__
except ImportError:  # package not installed; running from a source checkout
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
