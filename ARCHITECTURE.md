# Architecture

Expanded reference for `polartools` internals. See [CLAUDE.md](CLAUDE.md) for commands.

## Data flow

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
10. **`xmcd_gui.py`** / **`xanes_gui.py`** — PyQt6 + pyqtgraph GUIs for interactive XMCD/XANES normalization. Installed as the `xmcd-gui` / `xanes-gui` console scripts.
11. **`_gui_common.py`** — `GuiCommonMixin`: source-widget builders and small parsing/line-sync helpers shared by both GUIs. Keep genuinely shared, verbatim-identical logic here; leave GUI-specific branching (e.g. XMCD's dichro/lockin or H+/H− handling) in the individual GUI files.

## Catalog backends

`load_data.py` supports two catalog backends, selected automatically by `load_catalog()`:

- **databroker** (legacy, MongoDB): `load_catalog("catalog-name")` — tries this first. Registers area detector handlers (`LambdaHDF5Handler`, `EigerHandler`, `SPEHandler`) client-side.
- **tiled** (new, postgres/`bluesky-tiled-plugins`): falls back to `from_profile("profile-name")[tiled_path]` when the name is not a databroker catalog. `tiled_path` defaults to `"/raw"`. No client-side handler registration needed (server-side). `tiled` and `bluesky-tiled-plugins` are required dependencies.

Detection: `_is_tiled(obj)` checks `"tiled" in type(obj).__module__`. Both backends support `cat[scan_id]` integer indexing and `run.metadata["start"]`.

`manage_database.py` and `process_images.py` are databroker-only and do not yet support tiled.

## Key dependencies
- `lmfit` — peak fitting throughout diffraction and absorption modules
- `databroker` / `bluesky` — experimental data catalog (Bluesky ecosystem, legacy)
- `tiled` + `bluesky-tiled-plugins` — new catalog backend (required, postgres-backed)
- `spec2nexus` — reading SPEC data files
- `dask` — lazy/parallel image loading in `process_images.py`
- `PyQt6` + `pyqtgraph` — the XMCD/XANES GUIs
- `numpy`, `scipy`, `pandas`, `matplotlib` — standard scientific stack

## Test data

Fixtures live in `polartools/tests/data_for_test/`:
- `absorption.dat`, `fluorescence.dat`, `bluesky_spec.dat`, `pressure_calibration.dat` — SPEC files
- `scan_000025_master.hdf` — HDF5 master file
- `csv/` — CSV exports of scan 1049
- `databroker/` — msgpack catalog (scan 1049); `polartools/tests/conftest.py` registers it automatically at test-session start, no manual `databroker-unpack` step needed
- `lambda250k/` — databroker catalog with 50 Lambda detector images (scan 276)

GUI tests (`test_xmcd_gui.py`, `test_xanes_gui.py`) instantiate real `MainWindow` objects under a session-scoped, headless `QApplication` fixture — set `QT_QPA_PLATFORM=offscreen` if running outside conftest's default. No `pytest-qt`/`qtbot` needed.

**Not yet covered by tests** (need new data files to add coverage):
- `EigerHandler` — needs an Eiger-format HDF5 file (`entry/data/data` dataset)
- `SPEHandler` — needs a Princeton Instruments SPE file
- `process_rxes_mcd` — needs a scan whose image count is divisible by 4 (lambda250k has 50 images, which is not)
