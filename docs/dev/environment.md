# Environment, packaging and CI

Companion to [`CLAUDE.md`](../../CLAUDE.md). Read this before touching
dependencies, the conda environment, the test fixtures or the workflows.

## Packaging

`hatchling` + `hatch-vcs`. There is no `setup.py`, no `setup.cfg` and no
`requirements.txt` — dependencies, extras, console scripts, ruff, pytest and
coverage config all live in `pyproject.toml`.

Versions come from git tags (prefix `v`); `hatch-vcs` writes
`polartools/_version.py` at build time, so that file is generated and must not
be edited. A source checkout without tags reports `0.0.0+unknown` through the
`__init__.py` fallback.

`python-publish.yml` builds the sdist and wheel and uploads them when a GitHub
Release is **published** (or on manual dispatch) — pushing a tag alone does not
publish. The docs workflow deploys to `/latest` from `main`, to `/<tag>` from a
`v*` tag, and to `/dev` from any other branch.

## Dependency rules

**Declare anything a module imports directly**, even when it already arrives as
a transitive. `xarray`, `suitcase-utils` and `PyYAML` were undeclared for a
long time and worked only because databroker happened to pull them in.

**Add a version floor only where an API requires one**, with a comment saying
which API:

- `matplotlib>=3.5` — `backends.backend_qtagg` does not exist before it
- `pyqtgraph>=0.13` — first release handling PyQt6 and numpy 2 together
- `PyQt6>=6.2`, `h5py>=3`, `hdf5plugin>=3`

numpy needs no floor: `_scan_viewer_base.py` picks between `np.trapezoid` and
`np.trapz` at import, so 1.x and 2.x both work.

**`import hdf5plugin` is load-bearing.** It registers the LZ4/bitshuffle HDF5
codecs by side effect and is never referenced afterwards; the `# noqa: F401`
beside it is deliberate. Removing it makes compressed detector files
unreadable at runtime, not at import.

### The conda Qt trap

conda-forge's `pyqt` package is **PyQt5** (5.15.x); the viewers are PyQt6.
`environment.yml` must ask for `pyqt6`, which also pulls `qt6-main` and with it
`xcb-util-cursor` — without that library Qt ≥6.5 refuses to start the `xcb`
platform plugin and every GUI dies at launch. The PyPI PyQt6 wheel does not
carry it, so a pip install on Linux needs `libxcb-cursor0` from the system
package manager.

CI never catches this: it runs `QT_QPA_PLATFORM=offscreen`, which does not load
the xcb plugin.

## Tests

`polartools/tests/conftest.py` runs `databroker-unpack inplace` at import time
and sets `QT_QPA_PLATFORM=offscreen`, so neither a manual unpack nor a display
is needed (closes #34).

Fixtures in `polartools/tests/data_for_test/`:

- `absorption.dat`, `fluorescence.dat`, `bluesky_spec.dat`,
  `pressure_calibration.dat` — SPEC files
- `scan_000025_master.hdf` — HDF5 master file
- `csv/` — CSV exports of scan 1049
- `databroker/` — msgpack catalog (scan 1049), registered by conftest
- `lambda250k/` — databroker catalog with 50 Lambda images (scan 276)

Not covered, for want of a data file:

- `EigerHandler` — needs an Eiger-format HDF5 file (`entry/data/data`)
- `SPEHandler` — needs a Princeton Instruments SPE file
- `process_rxes_mcd` — needs a scan whose image count is divisible by 4
  (lambda250k has 50, which is not)
- every GUI module except `xmcd_gui` — all are in the coverage `omit` list

## CI workflows

| Workflow | Does |
| --- | --- |
| `python-package-conda.yml` | ruff check + ruff format --check; pytest on conda (3.10/3.11/3.12, from `environment.yml`) and on pip (3.10/3.12); coverage to codecov |
| `docs.yml` | `sphinx-build -W --keep-going`, then deploy to GitHub Pages |
| `python-publish.yml` | Build and upload sdist + wheel on a published release |
| `codeql.yml` | Security scan |

Both test jobs `apt-get install libegl1 libgl1` and run Qt offscreen. The pip
job is the one that would notice an undeclared dependency, since it installs
nothing from conda.

`ruff format` covers Python inside Markdown code blocks, so a snippet in
`docs/` can fail the lint job.

## Known issues

Snapshot — `gh issue list` for the current set.

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
