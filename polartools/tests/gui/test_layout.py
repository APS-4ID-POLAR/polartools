# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

"""The gui/ layout: modules import, scripts resolve, old paths are gone."""

import importlib
import importlib.util
from importlib.metadata import entry_points

import pytest

pytest.importorskip("PyQt6", reason="PyQt6 not installed — skipping GUI tests")

GUI_MODULES = [
    "_scan_viewer_base",
    "scan_viewer",
    "spectrum_viewer",
    "image_viewer",
    "flyscan_dichro",
    "xmcd",
    "xanes",
]

SCRIPTS = {
    "scan-viewer": "polartools.gui.scan_viewer:main",
    "spectrum-viewer": "polartools.gui.spectrum_viewer:main",
    "image-viewer": "polartools.gui.image_viewer:main",
    "flyscan-dichro-gui": "polartools.gui.flyscan_dichro:main",
    "xmcd-gui": "polartools.gui.xmcd:main",
    "xanes-gui": "polartools.gui.xanes:main",
}

# Removed in #135 with no shims left behind.
OLD_MODULES = [
    "_scan_viewer_base",
    "scan_viewer",
    "spectrum_viewer",
    "image_viewer",
    "flyscan_gui_dichro",
    "xmcd_gui",
    "xanes_gui",
]


@pytest.mark.parametrize("name", GUI_MODULES)
def test_gui_module_imports(name):
    importlib.import_module(f"polartools.gui.{name}")


@pytest.mark.parametrize("name", OLD_MODULES)
def test_old_module_path_is_gone(name):
    assert importlib.util.find_spec(f"polartools.{name}") is None


def test_console_scripts():
    # Reads the installed metadata: after editing [project.scripts], re-run
    # `pip install -e .` or this fails.
    found = {
        ep.name: ep.value
        for ep in entry_points(group="console_scripts")
        if ep.value.startswith("polartools.")
    }
    assert found == SCRIPTS
    for ep in entry_points(group="console_scripts"):
        if ep.name in SCRIPTS:
            assert callable(ep.load())


def test_scan_viewer_rebinds_the_base():
    from polartools.gui import _scan_viewer_base as base
    from polartools.gui import scan_viewer

    assert base.HdfViewer is scan_viewer.SgHdfViewer
    assert base.list_hdf_files is scan_viewer._scandir_hdf_files
