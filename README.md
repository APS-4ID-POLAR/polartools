# polartools

Python package to read and process data from the APS Polar beamline

[![Python version](https://img.shields.io/pypi/pyversions/polartools.svg)](https://pypi.python.org/pypi/polartools)
![pytest](https://github.com/APS-4ID-POLAR//polartools/actions/workflows/python-package-conda.yml/badge.svg)
[![pages-build-deployment](https://github.com/APS-4ID-POLAR/polartools/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/APS-4ID-POLAR/polartools/actions/workflows/pages/pages-build-deployment)

[![codecov](https://codecov.io/gh/APS-4ID-POLAR/polartools/branch/main/graph/badge.svg?token=pniQQrvLIA)](https://codecov.io/gh/APS-4ID-POLAR/polartools)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b94ee8138eec41d680ece65b77b8f7c2)](https://app.codacy.com/gh/APS-4ID-POLAR/polartools?utm_source=github.com&utm_medium=referral&utm_content=APS-4ID-POLAR/polartools&utm_campaign=Badge_Grade)
[![CodeQL](https://github.com/APS-4ID-POLAR/polartools/actions/workflows/codeql.yml/badge.svg)](https://github.com/APS-4ID-POLAR/polartools/actions/workflows/codeql.yml)

[![PyPi](https://img.shields.io/pypi/v/polartools.svg)](https://pypi.python.org/pypi/polartool)
[![release](https://img.shields.io/github/release/APS-4ID-POLAR/polartools.svg)](https://github.com/APS-4ID-POLAR/polartools/releases)
[![GitHub commits](https://img.shields.io/github/commits-since/APS-4ID-POLAR/polartools/v0.6.0.svg)](https://github.com/APS-4ID-POLAR/polartools/commit/)

-   License: [**ANL OPEN SOURCE LICENSE**](LICENSE.txt)
-   Documentation: [https://aps-4id-polar.github.io/polartools/latest/](https://aps-4id-polar.github.io/polartools/latest/)

## Install

```bash
pip install polartools
```

Or, for a development checkout with the heavier native dependencies (Qt, HDF5
plugins) from conda:

```bash
mamba env create -f environment.yml
mamba activate polartools-dev
```

See the [installation guide](https://aps-4id-polar.github.io/polartools/latest/installation.html)
for the X11 system libraries Qt needs on Linux.

## Graphical tools

Installing the package puts the beamline's viewers on your `PATH`:

| Command | What it opens |
| --- | --- |
| `scan-viewer [directory]` | Scan browser for a directory of `*_master.hdf` files, with the fluorescence-spectrum and image-stack windows behind its `Spc` and `Img` buttons |
| `flyscan-dichro-gui` | Post-scan flyscan dichroism viewer: CR/CL maps and their difference |
| `xmcd-gui` | Interactive XMCD normalization and dichroism |

## Library

```python
from polartools.load_data import load_table

table = load_table(scan_id=10, source="my_spec_file.dat", folder="/data/2024-1")
```

`load_table` reads SPEC files, HDF5, CSV folders, databroker catalogs and
tiled catalogs into a pandas `DataFrame`. See the
[quickstart](https://aps-4id-polar.github.io/polartools/latest/quickstart.html).
