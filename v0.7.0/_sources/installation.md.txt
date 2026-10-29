# Installation

## From PyPI

```bash
pip install polartools
```

## From conda-forge

If you already manage your scientific Python stack with conda or mamba:

```bash
mamba install -c conda-forge polartools
```

## From source

For development, clone and install in editable mode with the dev extras:

```bash
git clone https://github.com/APS-4ID-POLAR/polartools
cd polartools
pip install -e .[dev,docs]
```

Or with conda for the heavier native dependencies (Qt, HDF5 plugins):

```bash
git clone https://github.com/APS-4ID-POLAR/polartools
cd polartools
mamba env create -f environment.yml
mamba activate polartools-dev
```

The conda route is the recommended one on a beamline workstation: the
environment file installs `pyqt6` (and with it `qt6-main`, which carries the
X11 libraries Qt needs) from conda-forge rather than leaving Qt to a PyPI
wheel. Note that conda-forge's `pyqt` package is still PyQt5 — the viewers are
PyQt6, so ask for `pyqt6` if you build an environment by hand.

## Qt system libraries on Linux

The graphical tools need PyQt6, which comes in as a normal dependency either
way. What a `pip install` cannot provide is the handful of X11 libraries Qt
loads at runtime. If a viewer exits immediately with

```text
qt.qpa.plugin: Could not load the Qt platform plugin "xcb"
```

the system libraries are missing. On Debian/Ubuntu:

```bash
sudo apt-get install libxcb-cursor0 libegl1 libgl1
```

On RHEL/Rocky:

```bash
sudo dnf install xcb-util-cursor libglvnd-egl libglvnd-glx
```

`libxcb-cursor0` is the one most often absent: Qt 6.5 and newer refuse to start
the `xcb` platform plugin without it. A conda install of `pyqt6` brings its own
copy (`xcb-util-cursor`) and so sidesteps this entirely.

Over SSH, pass `-X` or `-Y` and check that `$DISPLAY` is set. With no display at
all, set `QT_QPA_PLATFORM=offscreen` — enough to import and test the GUI
modules, but nothing will be drawn.

## Verifying

```bash
python -c "import polartools; print(polartools.__version__)"
```

Should print a version string like `0.6.0` (or `0.6.0.postN+gHASH` for a
development checkout).

To check that the graphical side works too:

```bash
python -c "from PyQt6.QtCore import QT_VERSION_STR; print(QT_VERSION_STR)"
scan-viewer .
```

## Console scripts

Installing the package puts three applications on your `PATH`:

| Command | Opens | Documented in |
| --- | --- | --- |
| `scan-viewer [directory]` | The scan browser, and the spectrum and image viewers behind its buttons | [Scan viewers](examples/viewers.md) |
| `flyscan-dichro-gui [options]` | The post-scan flyscan dichroism viewer | [Flyscan dichroism](examples/flyscan_dichro.md) |
| `xmcd-gui` | Interactive XMCD normalization and dichroism | [XMCD GUI](examples/xmcd_gui.md) |

## Optional: enable pre-commit hooks

If you'll be contributing, install the pre-commit hooks so ruff runs on every
commit:

```bash
pre-commit install
```
