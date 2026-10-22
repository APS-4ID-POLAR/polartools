# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working in this repository.

**polartools** reads and processes data from the APS POLAR beamline (4-ID-D at
Argonne): X-ray absorption (XAS/XMCD), diffraction, pressure calibration,
image and flyscan reduction — plus the beamline's PyQt6 viewers, shipped as
console scripts.

Detail lives in two companion files. Read the relevant one before working in
that area; this file is only the index.

- **[docs/dev/architecture.md](docs/dev/architecture.md)** — what each module
  does, how the five GUI modules fit together, the two catalog backends.
- **[docs/dev/environment.md](docs/dev/environment.md)** — packaging,
  dependency rules, the conda/Qt trap, test fixtures and gaps, CI, releases,
  known issues.

User-facing documentation is `docs/source/` (MyST Markdown, published to
GitHub Pages); keep it in step when behaviour changes.

## Commands

```bash
pip install -e .[dev,docs]        # or: mamba env create -f environment.yml

pytest                            # conftest unpacks the databroker fixture
pytest --cov=polartools --cov-report=xml
pytest polartools/tests/test_absorption.py::test_name

ruff check . && ruff format --check .
sphinx-build -W --keep-going -b html docs/source docs/build/html
```

Ruff replaced flake8 and black: line length 80, target py310, rules `E`/`F`/`W`
with `E203`/`E501`/`E741` ignored, all configured in `pyproject.toml`.
`pre-commit install` wires `ruff --fix` and `ruff-format` into every commit.
The docs build is `-W`, so a docstring regression breaks CI.

## Repo map

```text
polartools/
  load_data.py              base layer: SPEC / CSV / HDF5 / databroker / tiled → DataFrame
  absorption.py             XAS, XMCD, normalization
  diffraction.py            Bragg peak fitting (lmfit), HKL
  pressure_calibration.py   Au/Ag standards, on top of diffraction
  process_images.py         area-detector reduction (dask)
  process_flyscan.py        flyscan reduction: position stream + detector ROI
  dichro_maps.py            CR/CL pairing, regridding, difference maps
  manage_database.py        databroker-pack / suitcase import-export
  area_detector_handlers.py Lambda / Eiger / SPE asset handlers
  _larch.py  _pyrixs.py     hoisted upstream helpers — do not reformat
  _scan_viewer_base.py      the scan browser (~9900 lines)
  scan_viewer.py            rebinds the base; `scan-viewer`
  spectrum_viewer.py        Vortex MCA cubes, behind the browser's Spc button
  image_viewer.py           area-detector stacks (pyqtgraph), behind Img
  flyscan_gui_dichro.py     `flyscan-dichro-gui`
  xmcd_gui.py               `xmcd-gui`
```

## House rules

1. **`_scan_viewer_base.py` is load-bearing** — it is what the beamline runs
   every day. `scan_viewer.py` rebinds it rather than forking it, so browser
   changes belong in the base, where both editions get them.
2. **Read a GUI module's docstring before changing one of its defaults.** They
   are long on purpose and record why each default is what it is, usually with
   the scan numbers that settled it.
3. **`import hdf5plugin` looks unused and is not** — it registers the HDF5
   codecs by side effect. The `# noqa: F401` beside it is deliberate.
4. **Declare every directly imported package in `pyproject.toml`**, even when
   databroker already pulls it in, and add a version floor only where an API
   requires one (with a comment saying which).
5. **conda-forge's `pyqt` is PyQt5.** The viewers need `pyqt6`, which also
   brings the `xcb-util-cursor` that Qt ≥6.5 cannot start without.
6. **`polartools/_version.py` is generated** by `hatch-vcs` from git tags. Do
   not edit it.
7. `image_viewer.open_hkl_window` imports a module that was never vendored.
   It is tracked separately — leave it alone.
