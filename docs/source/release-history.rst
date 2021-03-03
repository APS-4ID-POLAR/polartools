===============
Release History
===============

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
