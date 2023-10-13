===============
Release History
===============

-------------------
v0.4.0 (2023-10-12)
-------------------

**New features**

- Create `hkl_utils` module.

**Maintenance**

- Fix bug in image processing for python >= 3.9.
- Update the readthedocs documentation build process.
- Remove defunct LGTM support, adds CodeQL.
- Update pytest and flake8 workflow.

**Pull requests**

#68, #70, #71, #72, #73, #74

-------------------
v0.3.0 (2022-04-05)
-------------------

v0.2.4 didn't work with the pip integration again!

**Maintenance**

- remove dependence on `larch`.

**Pull requests**

#66.

-------------------
v0.2.4 (2022-04-05)
-------------------

v0.2.2 didn't work with the pip integration.

**New features**

- Create `process_images` module.
  
- Add XMCD specific functions to `absorption`.

**Maintenance**

- fixes to pytest.

**Pull requests**

#58, #62, #65

-------------------
v0.2.3 (2021-12-20)
-------------------

v0.2.2 didn't work with the pip integration.

**New features**

- functions to read and search the baseline.
  
- implement 3D plotting.
  
- integrate automatic Bluesky database loader from apstools.
  
- add area detector handlers.
  
- update the databroker import/export functions to handle external (image) files.

**Maintenance**

- fix `plot_2d` docs.
  
- changed how fit classes are handled in `diffraction`.
  
- fixes to pytest.
  
- fix `absorption` module to always export `numpy.array`.

**Pull requests**

#43, #45, #46, #47, #50, #51, #52, #53, #55, #57.

-------------------
v0.2.1 (2021-03-24)
-------------------

**New features**

- Additions to diffraction module:
  
  - Add the option to read_mesh into plot_2d.

  - Add load_axes in load_series and create fit_series.

- Create db_tools module:

  - Moved db_query from load_data.

  - Create collect_meta and show_meta.

**Maintenance**

Revert to ubuntu 18.04 for Github CI.


**Pull requests**

#37, #38, #40, #41, #42.

-------------------
v0.2.0 (2021-03-02)
-------------------

**New features**

- Absorption normalization functions.
  
- Daniel Haskel's `fluo` code (:func:`polartools.absorption.fluo_corr`) to correct for over-absorption in fluorescence data.

- Function to plot 2D data (:func:`polartools.diffraction.plot_2d`).

- Function to load a series of scans (:func:`polartools.diffraction.load_series`)

**Maintenance**

- All :mod:`absorption` functions now return `numpy.array`.

- :func:`polartools.load_data.load databroker` now uses the databroker "v1" because it is currently faster.

- Forces `databroker` >= 1.2.2, it fixes an issue after the latest `intake` update.

**Pull requests**

#29, 30, 31, 32, 33, 36.

-------------------
v0.1.4 (2021-02-09)
-------------------

- Fixes bug in :func:`polartools.load_data.is_Bluesky_specfile` and adds extra testing for it.

- Fixes bug in building the :mod:`manage_database` documentation.

-------------------
v0.1.3 (2021-02-08)
-------------------

Adds new modules and rearranges some functions:

- The :mod:`process_data` module is removed and its functions transfered to 
  the new :mod:`absorption` and :mod:`pressure_calibration` modules.

- Some absorption-related loading functions were also moved from 
  :mod:`load_data` to :mod:`absorption`.
  
- New :mod:`diffraction` module with functions to fit Bragg peaks.
  
- New :mod:`manage_database` module with functions to import/export data 
  from/to databroker.

-------------------
v0.1.2 (2021-01-11)
-------------------

Fixed the pypi issue (it was a problem with the long_description). Package uploaded.

-------------------
v0.1.1 (2021-01-11)
-------------------

Try to fix problems creating the pypi package.

-------------------
v0.1.0 (2021-01-11)
-------------------

Initial release. Include functions to load data and calculate pressure using XRD.
