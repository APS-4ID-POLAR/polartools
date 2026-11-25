# Architecture

Companion to [`CLAUDE.md`](../../CLAUDE.md). Read this before changing how
modules fit together; `CLAUDE.md` has the commands and the house rules.

## Library: data flow

The package is organized around data sources and processing stages.

| Module | Role |
| --- | --- |
| `load_data.py` | Base layer. Reads SPEC, CSV, HDF5, databroker or tiled into pandas DataFrames. Every higher-level module goes through it. |
| `absorption.py` | XAS/XMCD spectra: normalization, background subtraction, dichroism. |
| `diffraction.py` | Bragg peak fitting (Gaussian/Lorentzian/Pseudo-Voigt via lmfit) and HKL conversions. |
| `pressure_calibration.py` | Pressure from Au/Ag diffraction standards; builds on `diffraction.py`. |
| `process_images.py` | Area-detector images: thresholding, curvature extraction, dask for large stacks. |
| `process_flyscan.py` | Flyscan reduction: position stream → per-image x/y/I0, detector file → per-image ROI sum. Addresses files as `{folder}/{stream}/scan_{n:06d}.h5` (`FNAME_FORMAT`). |
| `dichro_maps.py` | Pairs CR/CL flyscans, regrids onto a common grid, differences, aligns, stitches, animates. Builds on `process_flyscan.py`. |
| `manage_database.py` | Import/export databroker data to msgpack or CSV/JSON via `databroker-pack` and `suitcase`. |
| `area_detector_handlers.py` | Databroker asset handlers for Lambda HDF5, Eiger and SPE. |
| `_larch.py` | Hoisted from xraylarch (`finde0`, `index_nearest`) to avoid the heavy larch install. |
| `_pyrixs.py` | Hoisted from pyrixs: 2D image → spectrum (curvature fitting, photon-event binning). |

`_larch.py` and `_pyrixs.py` are upstream re-exports: excluded from ruff,
coverage and autoapi, and not to be reformatted.

## GUIs

Five PyQt6 modules, vendored from standalone beamline scripts and ported from
PyQt5. Three are console scripts (`[project.scripts]`):

| Entry point | Module | Opens |
| --- | --- | --- |
| `scan-viewer` | `scan_viewer.py` | Scan browser for a directory of `*_master.hdf` |
| `flyscan-dichro-gui` | `flyscan_gui_dichro.py` | CR/CL flyscan difference maps |
| `xmcd-gui` | `xmcd_gui.py` | Interactive XMCD normalization |

**`_scan_viewer_base.py`** (~9900 lines) is the scan browser itself — listing,
field table, plots, fits, XMCD panels, colour maps. It is what the beamline
runs every day; treat it as load-bearing.

**`scan_viewer.py`** (~950 lines) is a *rebind* of that base, not a fork. It
imports the base and swaps module globals — `SpectrumViewerWindow` so the
browser opens the softglue spectrum viewer, `list_hdf_files` so directories are
listed with `scandir` (0.05 s against 82 s on a 2714-scan directory). Browser
behaviour therefore belongs in the base, where both editions get it; copying a
466 KB file would mean making every later fix twice.

**`spectrum_viewer.py`** (Vortex MCA cubes) and **`image_viewer.py`**
(area-detector stacks, pyqtgraph) are opened by the browser's `Spc` and `Img`
buttons, and also run as `python -m polartools.<module> <file>`. Neither the
Vortex element count nor the dead-time bookkeeping is constant between run
cycles, so both are probed for rather than assumed.

**`xmcd_gui.py`** is deliberately thin; the processing lives in
`absorption.py`. It is the only GUI module with tests
(`tests/test_xmcd_gui.py`, which builds a session-scoped `QApplication`).

The GUI module docstrings are long and carry the reasoning — usually with the
scan numbers that motivated it — behind every default. Read the docstring
before changing a default.

Known gap: `image_viewer.open_hkl_window` does `from hkl_viewer import
HklWindow`, a module that was never vendored, so the HKL-3D feature cannot
work. Tracked separately; leave it alone unless the issue says otherwise.

## Catalog backends

`load_data.py` supports two, selected automatically by `load_catalog()`:

- **databroker** (legacy, MongoDB): `load_catalog("catalog-name")` is tried
  first. Registers the area-detector handlers (`LambdaHDF5Handler`,
  `EigerHandler`, `SPEHandler`) client-side.
- **tiled** (postgres, via `bluesky-tiled-plugins`): falls back to
  `from_profile("profile-name")[tiled_path]` when the name is not a databroker
  catalog. `tiled_path` defaults to `"/raw"`. Handlers are registered
  server-side, so there is nothing to do client-side.

Detection is `_is_tiled(obj)`, which checks `"tiled" in type(obj).__module__`.
Both backends support `cat[scan_id]` integer indexing and
`run.metadata["start"]`.

`manage_database.py` and `process_images.py` are databroker-only and do not yet
support tiled.
