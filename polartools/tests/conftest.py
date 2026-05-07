"""
Copyright (c) 2020, UChicago Argonne, LLC.

See LICENSE file for details.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

# Must be set before any Qt import so the offscreen backend is selected.
# This allows GUI tests to run in headless CI environments (GitHub Actions).
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="session", autouse=True)
def _unpack_databroker_test_data():
    """Ensure the msgpack databroker test catalog is unpacked before tests run.

    Replaces the manual `databroker-unpack inplace ... data_3` step previously
    required in the README/CI. Idempotent: if the catalog has already been
    registered, this is a no-op.
    """
    data_dir = (
        Path(__file__).parent / "data_for_test" / "databroker"
    ).resolve()
    catalog_yml = data_dir / "catalog.yml"
    if catalog_yml.exists():
        return
    if shutil.which("databroker-unpack") is None:
        pytest.skip(
            "databroker-unpack CLI not available; cannot prepare test data"
        )
    subprocess.run(
        ["databroker-unpack", "inplace", str(data_dir), "data_3"],
        check=True,
    )
