# RXES Image Processor GUI

`polartools` ships an interactive PyQt6 application that wraps the RIXS/RXES
image-processing workflow described in the [Image processing
example](process_images.md). It's installed as a console script when you `pip
install polartools`:

```bash
process-images-gui
```

## What it does

The window has a shared **Data source** panel (a databroker/tiled catalog
name, or an HDF5 folder — see [Loading data](load_data.md)) used by all three
tabs below it:

- **Curvature** — loads and averages an image (`process_images.load_images`),
  fits its curvature (`process_images.get_curvature`), and displays the image
  with the fitted curve overlaid. The result can be saved to a text file.
- **RXES** — runs `process_images.process_rxes` using a curvature (typed in,
  or copied from the Curvature tab) and plots either a single spectrum or a
  2D RIXS map, depending on whether a positioner is set. Check "No curvature
  (integrate vertically)" to skip the curvature correction and just sum each
  image in the vertical direction.
- **RXES-MCD** — same as RXES but for `process_images.process_rxes_mcd`,
  plotting the RXES and MCD spectra/maps side by side. It has its own
  curvature fields and "No curvature" checkbox.

Each tab has its own "Save result" button — spectra are saved as text
columns, RIXS maps as `.npz` archives.

## Running headless (CI / scripted use)

Like `xmcd-gui`, this GUI uses Qt's `offscreen` platform when
`QT_QPA_PLATFORM=offscreen` is set — this is how the test suite exercises it
without a display. For scripted/batch processing, call
`polartools.process_images.process_rxes`/`process_rxes_mcd` directly instead.

## Architecture

The GUI is a single file: [`polartools/process_images_gui.py`][src]. All
processing is delegated to `polartools.process_images`; matplotlib canvases
(`FigureCanvasQTAgg`) are embedded directly so the existing
`polartools._pyrixs.plot_curvature` overlay can be reused as-is.

[src]: https://github.com/APS-4ID-POLAR/polartools/blob/main/polartools/process_images_gui.py

## See also

- [Image processing](process_images.md) — the underlying scriptable API
- [XMCD GUI](xmcd_gui.md) — the sibling GUI for XAS/XMCD processing
