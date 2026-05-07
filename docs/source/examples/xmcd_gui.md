# XMCD GUI

`polartools` ships an interactive PyQt6 application for the common XMCD
processing workflow. It's installed as a console script when you `pip install
polartools`, so you can launch it from anywhere:

```bash
xmcd-gui
```

## What it does

Wraps the same functions documented in the [XAS / XMCD example](absorption.md)
with a graphical workflow:

- Load paired ±-field scans from a SPEC file or databroker catalog
- Set the edge energy (E₀), pre-edge, and post-edge windows interactively
- See the normalized XAS and XMCD update live as you adjust parameters
- Save the processed XMCD to a column file with a metadata header

## Independent ±H normalization

As of v0.5.5, the +H and −H scans can be normalized independently —
useful when the two field directions have different self-absorption or
detector response. Toggle this from the **Settings** menu.

## Running headless (CI / scripted use)

The GUI uses Qt's `offscreen` platform when `QT_QPA_PLATFORM=offscreen` is
set — this is how the test suite exercises it without a display. Useful for
scripted batch processing on a server, but for that workflow you're usually
better off calling `polartools.absorption.process_xmcd` directly.

## Architecture

The GUI is a single file: [`polartools/xmcd_gui.py`][src]. It's intentionally
thin — almost all the actual processing is delegated to functions in
`polartools.absorption`. So if a scriptable equivalent of any GUI feature is
useful, you can usually find it (or build it) in that module.

[src]: https://github.com/APS-4ID-POLAR/polartools/blob/main/polartools/xmcd_gui.py

## Reporting GUI bugs

Open an issue with:

- The polartools version (`xmcd-gui --version` or `python -c "import polartools; print(polartools.__version__)"`)
- The PyQt6 version (`python -c "from PyQt6.QtCore import QT_VERSION_STR; print(QT_VERSION_STR)"`)
- A screenshot of the failure if it's visual
- The full traceback from the terminal you launched `xmcd-gui` from

## See also

- [XAS / XMCD](absorption.md) — the underlying scriptable API
