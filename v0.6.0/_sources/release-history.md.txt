# Release history

For the current and most recent releases, see the [GitHub releases page][releases].
The entries below are kept for historical context.

[releases]: https://github.com/APS-4ID-POLAR/polartools/releases

## v0.4.0 (2023-10-12)

**New features**

- Create `hkl_utils` module.

**Maintenance**

- Fix bug in image processing for Python ≥ 3.9.
- Update the readthedocs documentation build process.
- Remove defunct LGTM support, add CodeQL.
- Update pytest and flake8 workflow.

**Pull requests:** #68, #70, #71, #72, #73, #74

## v0.3.0 (2022-04-05)

`v0.2.4` didn't work with the pip integration again.

**Maintenance**

- Remove dependence on `larch`.

**Pull requests:** #66

## v0.2.4 (2022-04-05)

`v0.2.2` didn't work with the pip integration.

**New features**

- Create `process_images` module.
- Add XMCD-specific functions to `absorption`.

**Maintenance**

- Fixes to pytest.

**Pull requests:** #58, #62, #65

## v0.2.3 (2021-12-20)

`v0.2.2` didn't work with the pip integration.

**New features**

- Functions to read and search the baseline.
- Implement 3D plotting.
- Integrate automatic Bluesky database loader from `apstools`.
- Add area detector handlers.
- Update the databroker import/export functions to handle external (image) files.

**Maintenance**

- Fix `plot_2d` docs.
- Change how fit classes are handled in `diffraction`.
- Fixes to pytest.
- Fix `absorption` module to always export `numpy.array`.

**Pull requests:** #43, #45, #46, #47, #50, #51, #52, #53, #55, #57

## v0.2.1 (2021-03-24)

**New features**

- Diffraction module:
  - Add the option to `read_mesh` into `plot_2d`.
  - Add `load_axes` in `load_series` and create `fit_series`.
- New `db_tools` module:
  - Move `db_query` from `load_data`.
  - Create `collect_meta` and `show_meta`.

**Maintenance**

- Revert to ubuntu 18.04 for GitHub CI.

**Pull requests:** #37, #38, #40, #41, #42

## v0.2.0 (2021-03-02)

**New features**

- Absorption normalization functions.
- Daniel Haskel's `fluo` code (`polartools.absorption.fluo_corr`) to correct
  for over-absorption in fluorescence data.
- `polartools.diffraction.plot_2d` for plotting 2D data.
- `polartools.diffraction.load_series` to load a series of scans.

**Maintenance**

- All `absorption` functions now return `numpy.array`.
- `polartools.load_data.load_databroker` now uses the databroker "v1" because
  it is currently faster.
- Forces `databroker >= 1.2.2` to fix an issue after the latest `intake` update.

**Pull requests:** #29, #30, #31, #32, #33, #36

## v0.1.4 (2021-02-09)

- Fix bug in `polartools.load_data.is_Bluesky_specfile` and add extra testing.
- Fix bug in building the `manage_database` documentation.

## v0.1.3 (2021-02-08)

Add new modules and rearrange some functions:

- The `process_data` module is removed and its functions transferred to the new
  `absorption` and `pressure_calibration` modules.
- Some absorption-related loading functions were also moved from `load_data`
  to `absorption`.
- New `diffraction` module with functions to fit Bragg peaks.
- New `manage_database` module with functions to import/export data
  from/to databroker.

## v0.1.2 (2021-01-11)

Fixed the PyPI issue (it was a problem with the long_description). Package uploaded.

## v0.1.1 (2021-01-11)

Try to fix problems creating the PyPI package.

## v0.1.0 (2021-01-11)

Initial release. Includes functions to load data and calculate pressure using XRD.
