"""
Copyright (c) 2020, UChicago Argonne, LLC.

See LICENSE file for details.
"""

import os

# Must be set before any Qt import so the offscreen backend is selected.
# This allows GUI tests to run in headless CI environments (GitHub Actions).
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
