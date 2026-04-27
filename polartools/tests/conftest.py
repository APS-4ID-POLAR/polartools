"""
Copyright (c) 2020, UChicago Argonne, LLC.

See LICENSE file for details.
"""

import os

# Must be set before any Qt import so the offscreen backend is selected.
# This allows GUI tests to run in headless CI environments (GitHub Actions).
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Force pytest-qt to use PyQt6; pyqtgraph (conda-forge) can pull in PyQt5
# as well, and pytest-qt would otherwise default to PyQt5, causing an
# isinstance mismatch with our PyQt6-based MainWindow.
os.environ.setdefault("PYTEST_QT_API", "pyqt6")
