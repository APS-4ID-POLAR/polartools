# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**polartools** is a Python package for reading and processing data from the APS Polar beamline (4-ID-D at Argonne National Laboratory). It handles X-ray absorption (XAS/XMCD), diffraction, and spectroscopic data analysis.

## Commands

### Install (development)
```bash
python setup.py develop
# or
pip install -e .
```

### Run tests
```bash
# Unpack databroker test data first (required once):
databroker-unpack inplace polartools/tests/data_for_test/databroker data_3

# Run all tests with coverage:
pytest --cov-config=.coveragerc --cov=./ --cov-report=xml

# Run a single test file:
pytest polartools/tests/test_absorption.py

# Run a single test:
pytest polartools/tests/test_absorption.py::test_function_name
```

### Test data
Test fixtures live in `polartools/tests/data_for_test/`:
- `absorption.dat`, `fluorescence.dat`, `bluesky_spec.dat`, `pressure_calibration.dat` — SPEC files
- `scan_000025_master.hdf` — HDF5 master file
- `csv/` — CSV exports of scan 1049
- `databroker/` — msgpack catalog (scan 1049, must be unpacked before use)
- `lambda250k/` — databroker catalog with 50 Lambda detector images (scan 276)

**Not yet covered by tests** (need new data files to add coverage):
- `EigerHandler` — needs an Eiger-format HDF5 file (`entry/data/data` dataset)
- `SPEHandler` — needs a Princeton Instruments SPE file
- `process_rxes_mcd` — needs a scan whose image count is divisible by 4 (lambda250k has 50 images, which is not)

### Lint and format
```bash
flake8
black . -l 80 --check --diff
```

## Architecture

### Data flow
The package is organized around data sources and processing stages:

1. **`load_data.py`** — Base layer: reads raw data from SPEC files, CSV, HDF5, Bluesky databroker, or tiled catalogs into pandas DataFrames. All higher-level modules use these functions.
2. **`absorption.py`** — Loads and processes XAS/XMCD spectra (normalization, background subtraction, dichroism calculations).
3. **`diffraction.py`** — Loads and fits Bragg peaks (Gaussian/Lorentzian/Pseudo-Voigt via lmfit); also handles HKL conversions.
4. **`pressure_calibration.py`** — Pressure calibration using Au/Ag diffraction standards; builds on diffraction.py.
5. **`process_images.py`** — Image processing for area detectors (thresholding, curvature extraction) using dask for large datasets.
6. **`manage_database.py`** — Import/export Bluesky/databroker data to/from msgpack or CSV/JSON via `databroker-pack` and `suitcase`.
7. **`area_detector_handlers.py`** — Custom databroker asset handlers for Lambda HDF5, Eiger, and SPE detectors.
8. **`_larch.py`** — Hoisted utilities from xraylarch (`finde0`, `index_nearest`) to avoid the heavy larch install.
9. **`_pyrixs.py`** — Hoisted utilities from pyrixs for 2D image → spectrum extraction (curvature fitting, photon-event binning).

### Catalog backends

`load_data.py` supports two catalog backends, selected automatically by `load_catalog()`:

- **databroker** (legacy, MongoDB): `load_catalog("catalog-name")` — tries this first. Registers area detector handlers (`LambdaHDF5Handler`, `EigerHandler`, `SPEHandler`) client-side.
- **tiled** (new, postgres/`bluesky-tiled-plugins`): falls back to `from_profile("profile-name")` when the name is not a databroker catalog. No client-side handler registration needed (server-side). Install `tiled` and `bluesky-tiled-plugins` separately (commented out in `requirements.txt`).

Detection: `_is_tiled(obj)` checks `not hasattr(obj, "v2")`. Both backends support `cat[scan_id]` integer indexing and `run.metadata["start"]`.

`manage_database.py` and `process_images.py` are databroker-only and do not yet support tiled.

### Key dependencies
- `lmfit` — peak fitting throughout diffraction and absorption modules
- `databroker` / `bluesky` — experimental data catalog (Bluesky ecosystem, legacy)
- `tiled` + `bluesky-tiled-plugins` — new catalog backend (optional, postgres-backed)
- `spec2nexus` — reading SPEC data files
- `dask` — lazy/parallel image loading in `process_images.py`
- `numpy`, `scipy`, `pandas`, `matplotlib` — standard scientific stack

### Known issues (open GitHub issues)
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
- **#34** — `databroker-unpack` must be run before pytest
- **#28** — Slow to read data from databroker
- **#27** — Improve code coverage

### Versioning
Uses `versioneer` with git tags (prefix `v`). Do not manually edit `polartools/_version.py`.

### Linting config
- Flake8: max line length 115, ignores E203/E741/W503 (see `.flake8`)
- Black: line length 80 (CI enforces this)
