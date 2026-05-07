"""
Copyright (c) 2020, UChicago Argonne, LLC.

See LICENSE file for details.
"""

import os
import shutil
import subprocess
from pathlib import Path

# Must be set before any Qt import so the offscreen backend is selected.
# This allows GUI tests to run in headless CI environments (GitHub Actions).
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _ensure_databroker_catalog_registered():
    """Register the packed databroker test catalog with intake.

    Runs at conftest import time (i.e. before any test module imports
    `from databroker import catalog`) so the global catalog singleton sees
    the entry on first read. Replaces the manual `databroker-unpack`
    prerequisite previously documented in the README and run as a separate
    CI step (closes #34).

    `databroker-unpack inplace` is idempotent — re-running just rewrites the
    user-side intake config file pointing at the same packed data.
    """
    if shutil.which("databroker-unpack") is None:
        return
    data_dir = (
        Path(__file__).parent / "data_for_test" / "databroker"
    ).resolve()
    if not (data_dir / "catalog.yml").exists():
        return  # nothing to register
    subprocess.run(
        ["databroker-unpack", "inplace", str(data_dir), "data_3"],
        check=True,
    )


_ensure_databroker_catalog_registered()
