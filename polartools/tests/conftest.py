"""
Copyright (c) 2020, UChicago Argonne, LLC.

See LICENSE file for details.
"""

import os

# Must be set before any Qt import so the offscreen backend is selected.
# This allows GUI tests to run in headless CI environments (GitHub Actions).
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Tell pytest-qt to use PyQt6 only when it is actually installed.
# Setting this unconditionally causes an INTERNALERROR on Python versions
# where PyQt6 is unavailable (e.g. 3.8) or when libEGL is missing.
try:
    import PyQt6  # noqa: F401

    os.environ.setdefault("PYTEST_QT_API", "pyqt6")
except ImportError:
    pass
