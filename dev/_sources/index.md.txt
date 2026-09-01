# polartools

**Read and process data from the APS POLAR beamline (4-ID-D).**

`polartools` is a Python package for working with the data products of the POLAR
beamline at the Advanced Photon Source: X-ray absorption (XAS/XMCD), diffraction
(Bragg peak fitting, HKL conversions), pressure calibration from Au/Ag
standards, area-detector image processing, and import/export of Bluesky
[databroker]/[tiled] catalogs.

[databroker]: https://blueskyproject.io/databroker/
[tiled]: https://blueskyproject.io/tiled/

```{toctree}
:maxdepth: 1
:hidden:
:caption: Get started

installation
quickstart
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: User guide

examples/load_data
examples/absorption
examples/diffraction
examples/pressure_calibration
examples/process_images
examples/manage_database
examples/xmcd_gui
```

```{toctree}
:maxdepth: 1
:hidden:
:caption: Reference

api/polartools/index
release-history
```

---

::::{grid} 2 2 3 3
:gutter: 3

:::{grid-item-card} Installation
:link: installation
:link-type: doc

Install from PyPI, conda, or source. Set up a development environment.
:::

:::{grid-item-card} Quickstart
:link: quickstart
:link-type: doc

Load a scan, plot it, fit a peak — in five minutes.
:::

:::{grid-item-card} Loading data
:link: examples/load_data
:link-type: doc

Read SPEC, CSV, HDF5, databroker, and tiled scans into pandas DataFrames.
:::

:::{grid-item-card} XAS / XMCD
:link: examples/absorption
:link-type: doc

Edge-step normalize, compute dichroism, plot ±-field spectra.
:::

:::{grid-item-card} Diffraction
:link: examples/diffraction
:link-type: doc

Fit Bragg peaks (Gaussian / Lorentzian / Voigt) and convert HKL.
:::

:::{grid-item-card} Pressure calibration
:link: examples/pressure_calibration
:link-type: doc

Au and Ag diffraction-standard pressure calibration.
:::

:::{grid-item-card} Image processing
:link: examples/process_images
:link-type: doc

Threshold, curvature extraction, and RIXS image-to-spectrum reduction.
:::

:::{grid-item-card} Database management
:link: examples/manage_database
:link-type: doc

Pack/unpack databroker catalogs, export to CSV/JSON.
:::

:::{grid-item-card} XMCD GUI
:link: examples/xmcd_gui
:link-type: doc

The bundled `xmcd-gui` PyQt application for interactive XMCD processing.
:::

:::{grid-item-card} API reference
:link: api/polartools/index
:link-type: doc

Auto-generated reference for every public function and class.
:::

::::

---

## Getting help

Open a [discussion] or [file an issue] on the GitHub project. Bug reports are
welcome — please include the scan number, the full traceback, and the version
of polartools you're running (`python -c "import polartools; print(polartools.__version__)"`).

[discussion]: https://github.com/APS-4ID-POLAR/polartools/discussions
[file an issue]: https://github.com/APS-4ID-POLAR/polartools/issues
