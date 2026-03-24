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

### Lint and format
```bash
flake8
black . -l 80 --check --diff
```

## Architecture

### Data flow
The package is organized around data sources and processing stages:

1. **`load_data.py`** — Base layer: reads raw data from SPEC files, CSV, HDF5, or Bluesky databroker into pandas DataFrames. All higher-level modules use these functions.
2. **`absorption.py`** — Loads and processes XAS/XMCD spectra (normalization, background subtraction, dichroism calculations).
3. **`diffraction.py`** — Loads and fits Bragg peaks (Gaussian/Lorentzian/Pseudo-Voigt via lmfit); also handles HKL conversions.
4. **`pressure_calibration.py`** — Pressure calibration using Au/Ag diffraction standards; builds on diffraction.py.
5. **`process_images.py`** — Image processing for area detectors (thresholding, curvature extraction) using dask for large datasets.
6. **`manage_database.py`** — Import/export Bluesky/databroker data to/from HDF5.
7. **`area_detector_handlers.py`** — Custom databroker asset handlers for Lambda HDF5, Eiger, and SPE detectors.
8. **`_larch.py`**, **`_pyrixs.py`** — Optional wrappers for Larch and PYRIXS spectroscopy libraries.

### Key dependencies
- `lmfit` — peak fitting throughout diffraction and absorption modules
- `databroker` / `bluesky` — experimental data catalog (Bluesky ecosystem)
- `spec2nexus` — reading SPEC data files
- `dask` — lazy/parallel image loading in `process_images.py`
- `numpy`, `scipy`, `pandas`, `matplotlib` — standard scientific stack

### Versioning
Uses `versioneer` with git tags (prefix `v`). Do not manually edit `polartools/_version.py`.

### Linting config
- Flake8: max line length 115, ignores E203/E741/W503 (see `.flake8`)
- Black: line length 80 (CI enforces this)
