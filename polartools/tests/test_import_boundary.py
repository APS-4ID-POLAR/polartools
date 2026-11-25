# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

"""The analysis library must import without Qt (issue #135).

Qt is a required dependency, so packaging does not stop a library module from
importing it; this test does. Each import runs in a fresh interpreter, because
by the time this file runs another test may already have loaded PyQt6.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

import polartools

LIBRARY_MODULES = [
    "polartools",
    "polartools._larch",
    "polartools._pyrixs",
    "polartools.absorption",
    "polartools.area_detector_handlers",
    "polartools.dichro_maps",
    "polartools.diffraction",
    "polartools.load_data",
    "polartools.manage_database",
    "polartools.pressure_calibration",
    "polartools.process_flyscan",
    "polartools.process_images",
    "polartools.spectrum_fit",
]

CHECK = """
import importlib, sys
importlib.import_module(sys.argv[1])
bad = sorted(
    name for name in sys.modules
    if name.split(".")[0] in ("PyQt5", "PyQt6", "pyqtgraph")
    or name.startswith("polartools.gui")
    or (name.startswith("matplotlib.backends.") and "_qt" in name)
)
print("\\n".join(bad))
sys.exit(1 if bad else 0)
"""


def test_every_library_module_is_checked():
    # A new top-level module is library code and belongs in the list above.
    package = Path(polartools.__file__).parent
    found = {
        f"polartools.{path.stem}"
        for path in package.glob("*.py")
        if path.stem not in ("__init__", "_version")
    }
    assert found == set(LIBRARY_MODULES) - {"polartools"}


@pytest.mark.parametrize("module", LIBRARY_MODULES)
def test_library_module_imports_without_qt(module):
    # Agg, so that a matplotlibrc asking for QtAgg is not mistaken for a
    # polartools import.
    env = dict(os.environ, MPLBACKEND="Agg")
    result = subprocess.run(
        [sys.executable, "-c", CHECK, module],
        capture_output=True,
        text=True,
        env=env,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"{module} imported GUI modules:\n{result.stdout}{result.stderr}"
    )
