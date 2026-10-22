# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**polartools** is a Python package for reading and processing data from the APS Polar beamline (4-ID-D at Argonne National Laboratory). It handles X-ray absorption (XAS/XMCD), diffraction, and spectroscopic data analysis, and it ships the beamline's PyQt6 viewers as console scripts.

## Commands

### Install (development)
```bash
pip install -e .[dev,docs]
# or, with the native deps (Qt, HDF5 plugins) from conda:
mamba env create -f environment.yml && mamba activate polartools-dev
```

The build backend is `hatchling` + `hatch-vcs`; there is no `setup.py` and no `requirements.txt`. Everything — dependencies, extras, entry points, ruff, pytest and coverage config — lives in `pyproject.toml`.

### Run tests
```bash
# Run all tests with coverage:
pytest --cov=polartools --cov-report=xml

# Run a single test file:
pytest polartools/tests/test_absorption.py

# Run a single test:
pytest polartools/tests/test_absorption.py::test_function_name
```

`polartools/tests/conftest.py` runs `databroker-unpack inplace` at import time and sets `QT_QPA_PLATFORM=offscreen`, so neither the manual unpack step nor a display is needed (closes #34).

### Lint and format
```bash
ruff check .
ruff format --check .
```

Ruff replaced flake8 and black. Config is in `pyproject.toml`: line length 80, target py310, rules `E`/`F`/`W` with `E203`/`E501`/`E741` ignored. `_version.py`, `_larch.py`, `_pyrixs.py` and `docs/source/conf.py` are excluded. `pre-commit install` wires ruff (with `--fix`) and ruff-format into every commit.

### Build the docs
```bash
sphinx-build -W --keep-going -b html docs/source docs/build/html
```

CI builds with `-W`, so a docstring regression breaks the build. Pages are Markdown (MyST); the API reference is generated statically by sphinx-autoapi from every module except tests, `_version`, `_larch` and `_pyrixs`.

### Test data
Test fixtures live in `polartools/tests/data_for_test/`:
- `absorption.dat`, `fluorescence.dat`, `bluesky_spec.dat`, `pressure_calibration.dat` — SPEC files
- `scan_000025_master.hdf` — HDF5 master file
- `csv/` — CSV exports of scan 1049
- `databroker/` — msgpack catalog (scan 1049), registered automatically by conftest
- `lambda250k/` — databroker catalog with 50 Lambda detector images (scan 276)

**Not yet covered by tests** (need new data files to add coverage):
- `EigerHandler` — needs an Eiger-format HDF5 file (`entry/data/data` dataset)
- `SPEHandler` — needs a Princeton Instruments SPE file
- `process_rxes_mcd` — needs a scan whose image count is divisible by 4 (lambda250k has 50 images, which is not)
- the vendored GUI modules, except `xmcd_gui` — all are in the coverage `omit` list

## Architecture

### Library: data flow
The package is organized around data sources and processing stages:

1. **`load_data.py`** — Base layer: reads raw data from SPEC files, CSV, HDF5, Bluesky databroker, or tiled catalogs into pandas DataFrames. All higher-level modules use these functions.
2. **`absorption.py`** — Loads and processes XAS/XMCD spectra (normalization, background subtraction, dichroism calculations).
3. **`diffraction.py`** — Loads and fits Bragg peaks (Gaussian/Lorentzian/Pseudo-Voigt via lmfit); also handles HKL conversions.
4. **`pressure_calibration.py`** — Pressure calibration using Au/Ag diffraction standards; builds on diffraction.py.
5. **`process_images.py`** — Image processing for area detectors (thresholding, curvature extraction) using dask for large datasets.
6. **`process_flyscan.py`** — Flyscan reduction: position stream → per-image x/y/I0, detector file → per-image ROI sum. Addresses files as `{folder}/{stream}/scan_{n:06d}.h5`.
7. **`dichro_maps.py`** — Pairs CR/CL flyscans, regrids them onto a common grid, differences, aligns, stitches and animates. Builds on `process_flyscan.py`.
8. **`manage_database.py`** — Import/export Bluesky/databroker data to/from msgpack or CSV/JSON via `databroker-pack` and `suitcase`.
9. **`area_detector_handlers.py`** — Custom databroker asset handlers for Lambda HDF5, Eiger, and SPE detectors.
10. **`_larch.py`** — Hoisted utilities from xraylarch (`finde0`, `index_nearest`) to avoid the heavy larch install.
11. **`_pyrixs.py`** — Hoisted utilities from pyrixs for 2D image → spectrum extraction (curvature fitting, photon-event binning).

### GUIs

Five PyQt6 modules, vendored from standalone beamline scripts and ported to PyQt6. Three are console scripts (`[project.scripts]` in `pyproject.toml`):

| Entry point | Module | Opens |
| --- | --- | --- |
| `scan-viewer` | `scan_viewer.py` | Scan browser for a directory of `*_master.hdf` |
| `flyscan-dichro-gui` | `flyscan_gui_dichro.py` | CR/CL flyscan difference maps |
| `xmcd-gui` | `xmcd_gui.py` | Interactive XMCD normalization |

- **`_scan_viewer_base.py`** (~9900 lines) is the scan browser itself — listing, field table, plots, fits, XMCD panels, colour maps. It is what the beamline runs every day; treat it as load-bearing.
- **`scan_viewer.py`** is a thin *rebind* of that base, not a fork: it swaps module globals (`SpectrumViewerWindow`, `list_hdf_files`) so the browser opens the softglue spectrum viewer and lists directories with `scandir`. Changes to browser behaviour belong in the base so both editions get them.
- **`spectrum_viewer.py`** (Vortex MCA cubes) and **`image_viewer.py`** (area-detector stacks, pyqtgraph) are opened by the browser's `Spc` and `Img` buttons, and are also runnable as `python -m polartools.<module> <file>`.
- `xmcd_gui.py` is deliberately thin — the processing lives in `absorption.py`.

The module docstrings are long and carry the reasoning (and the scan numbers) behind each default. Read the docstring before changing a default.

Known gap: `image_viewer.open_hkl_window` does `from hkl_viewer import HklWindow`, a module that was never vendored, so the HKL-3D feature cannot work. Tracked separately — leave it alone unless the issue says otherwise.

### Catalog backends

`load_data.py` supports two catalog backends, selected automatically by `load_catalog()`:

- **databroker** (legacy, MongoDB): `load_catalog("catalog-name")` — tries this first. Registers area detector handlers (`LambdaHDF5Handler`, `EigerHandler`, `SPEHandler`) client-side.
- **tiled** (new, postgres/`bluesky-tiled-plugins`): falls back to `from_profile("profile-name")[tiled_path]` when the name is not a databroker catalog. `tiled_path` defaults to `"/raw"`. No client-side handler registration needed (server-side).

Detection: `_is_tiled(obj)` checks `"tiled" in type(obj).__module__`. Both backends support `cat[scan_id]` integer indexing and `run.metadata["start"]`.

`manage_database.py` and `process_images.py` are databroker-only and do not yet support tiled.

### Key dependencies
- `lmfit` — peak fitting throughout diffraction and absorption modules
- `databroker` / `bluesky` — experimental data catalog (Bluesky ecosystem, legacy)
- `tiled` + `bluesky-tiled-plugins` — new catalog backend (required, postgres-backed)
- `PyQt6`, `pyqtgraph`, `matplotlib` (`backend_qtagg`) — the GUIs
- `h5py` + `hdf5plugin` — HDF5; `hdf5plugin` is imported for its side effect of registering the LZ4/bitshuffle codecs, so it looks unused and is not
- `spec2nexus` — reading SPEC data files
- `dask` — lazy/parallel image loading in `process_images.py`
- `numpy`, `scipy`, `pandas`, `xarray`, `matplotlib` — standard scientific stack

Declare a dependency in `pyproject.toml` whenever a module imports it directly, even if it already arrives as a transitive of databroker. Version floors are only added where an API actually requires one (`matplotlib>=3.5` for `backend_qtagg`, `pyqtgraph>=0.13` for PyQt6 with numpy 2) — don't pin for the sake of pinning.

Note for conda: conda-forge's `pyqt` is still PyQt5. `environment.yml` must ask for `pyqt6`, which also pulls `qt6-main` and with it `xcb-util-cursor` — without that library Qt ≥6.5 cannot load the `xcb` platform plugin and every GUI dies at startup. CI never catches this because it runs offscreen.

### Known issues (open GitHub issues)
Snapshot — check GitHub for the current list.
- **#92** — Bug when loading lockin XMCD with "h5" source
- **#89** — Moving devices from polar_instrument to here
- **#69** — Update the lambda handler
- **#61** — Alternative ways to process images
- **#60** — `manage_database.to_databroker` bug
- **#59** — Speed up image processing
- **#56** — `manage_database.to_csv_json` does not work under pytest (logging conflict)
- **#54** — Cannot pack run with Eiger images
- **#48** — Refactor polartools functions
- **#44** — Interpolation bug in `load_multi` absorption functions
- **#35** — Improve handling of fluorescence scans by `load_absorption`
- **#28** — Slow to read data from databroker
- **#27** — Improve code coverage

### Versioning
Versions come from git tags (prefix `v`) via `hatch-vcs`, which writes `polartools/_version.py` at build time. Do not edit that file.
