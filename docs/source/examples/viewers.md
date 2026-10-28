# Scan viewers

Three PyQt6 windows for looking at the files `MyNXWriter` leaves in a data
directory: a scan browser, a fluorescence-spectrum viewer, and an image-stack
viewer. They were developed at the beamline as standalone scripts and are now
part of the package, so they install and update with everything else.

Only the browser has a console script; it opens the other two itself.

```bash
scan-viewer [directory]
```

With no argument it starts in the current directory.

## The scan browser

`scan-viewer` lists every `*_master.hdf` in the directory that carries Bluesky
scan metadata — one row per scan with its number, plan, motor, range, points
and count time, filled in as the listing is read rather than on selection.
Detector image files, whether beside the masters or in a per-detector
subdirectory, are not listed as scans; instead an **Img** tick marks the scans
whose images were found on disk.

Select one scan, or ctrl/shift-click several, and they are overlaid:

- **Stream box** — picks which of `entry/instrument/bluesky/streams` to show.
  Most scans only record `primary`; a dichro scan also writes the XAS and XMCD
  it computed to `dichro_monitor` and opens on that.
- **Field table** — X / Y / Norm ticks drive the plot. X and Norm take one
  field each, Y any number. They start on whatever the scan's own `hints`
  nominate, so a scan plots itself the way it was meant to be read. Norm
  starts clear, because normalizing rescales every number on the plot and in
  the statistics table.
- **log(mon/det)** — plots the absorption rather than a plain ratio, which is
  what makes a transmission XAS scan readable straight off `primary`.
- **Statistics, Integrate, Fit** — max, min, half-max centre and FWHM per
  curve; trapezoidal and background-subtracted areas; and peak fitting
  (Gauss, Lorentzian, Lorentzian², pseudo-Voigt) over a sloping line, with
  centre, FWHM, height, area and R² reported per curve. **Multi** fits several
  peaks at once over one shared line.
- **Show dY/dX** — overlays each curve's derivative and moves the fit onto it,
  which is how a knife edge or an absorption edge becomes measurable.
- **Log X / Log Y** — beside the plot toolbar, put either axis on a log
  scale (the right-hand Y axes too). Points at or below zero are left out.
- **XMCD panels** — reduces the selected scans to XAS, XMCD, XMuCD and XNCD
  against energy, grouping each energy's points into the four
  (polarization, field) quadrants.
- **Remove spikes** — drops points taken while the beam was gone, judged by
  the monitor against its own scan's median. Only lone points go; a stretch of
  dead points is a beam-down run and is left visible.
- **Two-motor grid scans** are drawn as colour maps instead, one per ticked Y
  field, with **Clip map scale** for the odd wild cell and **Center on 0** for
  signed channels in a diverging colormap.
- **Auto refresh** adds scans to the listing as they land, so the browser can
  be left open during a run. A scan whose file is still being written gets a
  row with blank cells rather than being dropped.

:::{note}
Read the module docstrings for the full behaviour — they are long and they
explain *why* each default is what it is, usually with the scan numbers that
motivated it. They are reproduced in the
[API reference](../api/polartools/index.rst).
:::

## The spectrum viewer — `Spc`

The multi-element silicon-drift (Vortex) detector does not write its MCA
spectra into the master file. They go into a file of their own, named in
`entry/instrument/bluesky/metadata/detectors_file_relative_path`:

```text
vortex: vortex/scan_000055.h5
```

holding `entry/data/data` shaped `(frames, elements, bins)`. The browser's
**Spc** button opens that file for the selected scan. It can also be opened
directly:

```bash
python -m polartools.spectrum_viewer path/to/scan_000055_master.hdf
```

Both the master file and the spectrum file itself are accepted — the master is
resolved to its companion.

Neither the element count nor the dead-time bookkeeping is constant between
run cycles, so both are read from the file. Dead-time correction
(`corrected = raw * DTFactor`, or `ICR / OCR` on older files) is **on by
default**, against this package's general habit of showing what was recorded:
on one test scan element 3 runs 8–33% dead while element 6 runs 0.2%, so seven
raw curves in a window meant for comparing elements state a wrong conclusion
about the hardware rather confidently. The checkbox turns it off.

The **Combine** box handles scans that took more than one frame per point, and
`vortex_total_el1` / `vortex_total_el7` are offered as monitors.

## The image viewer — `Img`

Area-detector frames are likewise written outside the master file, by one of
two writers — the pvapy streaming consumer (a sibling file named in
`pva_streaming`) or an areaDetector plugin (a per-detector subdirectory, named
in `detectors_file_relative_path`). The viewer resolves either, and falls back
to deriving `pva_<stem>.h5` when the master records neither.

```bash
python -m polartools.image_viewer path/to/pva_tise2_2_00006.h5 [more.h5 ...]
```

Built on pyqtgraph: pan/zoom, level and colormap control, rectangular ROIs
with live sums, and per-frame metadata (`entry/data/metadata/ca/<pv>` from
pvapy, `entry/instrument/NDAttributes/<name>` from areaDetector) plotted
against frame number. Several files can be opened side by side with ROIs
shared between the panels.

Compressed frames (LZ4, bitshuffle) need `hdf5plugin`, which is a dependency
of the package and imported for its codec registration.

## Troubleshooting

**The window never appears, and the terminal says `Could not load the Qt
platform plugin "xcb"`.** Missing X11 system libraries — see
[Qt system libraries on Linux](../installation.md).

**A scan is missing from the listing.** Only files carrying Bluesky scan
metadata are listed. A detector file, or a master written by something other
than `MyNXWriter`, is skipped on purpose.

**`Spc` or `Img` is greyed out or finds nothing.** The companion file is
resolved relative to the master's own directory, because the paths recorded
inside name the filesystem the scan was taken on — which a copy of the data is
not. Keep the subdirectory layout intact when copying data off the beamline.

## See also

- [XMCD GUI](xmcd_gui.md) — interactive XMCD normalization
- [Flyscan dichroism](flyscan_dichro.md) — the post-scan flyscan viewer
- [Image processing](process_images.md) — the scriptable image reduction API
