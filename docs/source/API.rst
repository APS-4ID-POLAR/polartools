===
API
===

---------
Load data
---------

The `load_data` module contains functions to load the various sources of data at 4-ID-D. 

Load a single scans:

.. autofunction:: polartools.load_data.load_scan

.. autofunction:: polartools.load_data.load_absorption

.. autofunction:: polartools.load_data.load_lockin

.. autofunction:: polartools.load_data.load_dichro

Load multiple spectroscopy scans:

.. autofunction:: polartools.load_data.load_xanes

.. autofunction:: polartools.load_data.load_xmcd

The following are some of the core data loaders:

.. autofunction:: polartools.load_data.load_spec

.. autofunction:: polartools.load_data.load_csv

.. autofunction:: polartools.load_data.load_databroker

.. autofunction:: polartools.load_data.run_v2_query

.. autofunction:: polartools.load_data.load_table

------------
Process data
------------

The `process_data` module has general data processing functions.

.. autofunction:: polartools.process_data.normalize_absorption

.. autofunction:: polartools.process_data.fit_bragg_peak

--------------------
Pressure Calibration
--------------------

The `pressure_calibration` module has functions that calibrate the pressure in a diamond anvil cell using Ag or Au x-ray diffraction.

.. autofunction:: polartools.pressure_calibration.load_ag_params

.. autofunction:: polartools.pressure_calibration.load_au_params

.. autofunction:: polartools.pressure_calibration.calculate_pressure

.. autofunction:: polartools.pressure_calibration.xrd_calibrate_pressure
