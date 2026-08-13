# XANES GUI

`polartools` ships an interactive PyQt6 application for XANES normalization.
It's installed as a console script when you `pip install polartools`, so you
can launch it from anywhere:

```bash
xanes-gui
```

## What it does

Wraps the same functions documented in the [XAS / XMCD example](absorption.md)
with a graphical workflow, for the simpler single-spectrum case (no ±-field
pairing):

- Load one set of scans from a SPEC file or databroker catalog, or load
  directly from a plain column file by picking which columns are energy and
  mu
- Set the edge energy (E₀), pre-edge, and post-edge windows interactively
- See the normalized and flattened XANES update live as you adjust parameters
- Save the processed XANES to a column file with a metadata header, to a
  folder of your choosing
- Overlay a previously saved XANES file (same format) on the normalized plot
  for comparison

## Running headless (CI / scripted use)

The GUI uses Qt's `offscreen` platform when `QT_QPA_PLATFORM=offscreen` is
set — this is how the test suite exercises it without a display. Useful for
scripted batch processing on a server, but for that workflow you're usually
better off calling `polartools.absorption.normalize_absorption` directly.

## Architecture

The GUI is a single file: [`polartools/xanes_gui.py`][src]. It's intentionally
thin — almost all the actual processing is delegated to functions in
`polartools.absorption`. So if a scriptable equivalent of any GUI feature is
useful, you can usually find it (or build it) in that module.

[src]: https://github.com/APS-4ID-POLAR/polartools/blob/main/polartools/xanes_gui.py

## Reporting GUI bugs

Open an issue with:

- The polartools version (`xanes-gui --version` or `python -c "import polartools; print(polartools.__version__)"`)
- The PyQt6 version (`python -c "from PyQt6.QtCore import QT_VERSION_STR; print(QT_VERSION_STR)"`)
- A screenshot of the failure if it's visual
- The full traceback from the terminal you launched `xanes-gui` from

## See also

- [XAS / XMCD](absorption.md) — the underlying scriptable API
- [XMCD GUI](xmcd_gui.md) — the paired ±-field counterpart to this GUI
