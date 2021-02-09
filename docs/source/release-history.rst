===============
Release History
===============

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
