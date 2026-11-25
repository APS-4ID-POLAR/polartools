# Flyscan dichroism

A dichro flyscan measurement is a pair of scans taken at the two phase-plate
settings; the map of interest is the difference between them. `polartools`
covers this in three layers:

| Module | Role |
| --- | --- |
| {py:mod}`polartools.process_flyscan` | Reduce one scan: position stream → per-image `x`, `y`, `I0`; detector file → per-image ROI sum |
| {py:mod}`polartools.dichro_maps` | Pair scans up, regrid them, difference, align, stitch and animate |
| {py:mod}`polartools.flyscan_gui_dichro` | The `flyscan-dichro-gui` application over both |

## Expected layout

The reduction addresses files by scan number inside a sample folder:

```text
<folder>/pos_stream/scan_000490.h5    interferometer positions and I0
<folder>/eiger/scan_000490.h5         area-detector frames
<folder>/vortex/scan_000490.h5        fluorescence spectra
```

That pattern is `process_flyscan.FNAME_FORMAT`; pass `fname_format=` to any of
the reduction functions if your directory is laid out differently.

## The GUI

```bash
flyscan-dichro-gui [--folder NdFeB] [--cr 490] [--cl 491] [--config FILE]
```

It loads two scans — circular-right and circular-left — draws their per-image
scalar maps side by side, and computes a `(CL - CR)` difference map
interpolated onto a common `(x, y)` grid. Either polarization can be disabled,
in which case only that map is drawn and the difference panel stays blank, so
the single-scan case is covered by the same window.

Eiger and Vortex each get a tab, with their own ROI and their own copy of the
three maps. `--cr` and `--cl` accept a single scan, a comma-separated list, or
a strided range: `490`, `490,492,494`, `754-764:2`.

### Settings

Settings live in a JSON file rather than in `QSettings`, because the folder,
ROIs, scan numbers and colour limits worth keeping are the ones for the
analysis directory being worked in — a second experiment should not inherit
them.

The file is read at startup (`./flyscan_gui_dichro.config.json` unless
`--config` names another), but **nothing is written until the Autosave button
in the status bar is pressed**, so launching in a directory leaves no trace in
it. The destination can be changed while running; neither it nor the Autosave
state carries over to the next run, by design, so a session never writes
somewhere forgotten.

## Scripting the same reduction

Everything the GUI draws is available without it.

```python
from polartools.process_flyscan import reduce_position_stream, process_images

folder = "/data/2026-3/NdFeB"

# Per-image positions and monitor, grouped by detector trigger.
groups = reduce_position_stream(490, folder)
groups.xs, groups.ys, groups.i0s  # one value per image

# Eiger ROI sum, read in batches so a long scan doesn't land in memory at once.
roi_sum = process_images(490, folder, roi=((678, 888), (333, 480)), batch=100)
```

`process_vortex` does the same for a fluorescence ROI, and
`read_scan_geometry` returns the beam centre, wavelength and detector distance
recorded for the scan.

For the pair-wise work:

```python
from polartools.dichro_maps import (
    build_series,
    find_pairs,
    pair_difference,
    plot_difference_series,
)

roi = ((678, 888), (333, 480))

# Consecutive scans sharing a label value but differing in phase plate,
# returned as [(label, scan_a, scan_b), ...] with scan_a the lower number.
pairs = find_pairs(range(490, 520), "magnet_field", folder)

label, scan_a, scan_b = pairs[0]
XG, YG, difference, mean = pair_difference(scan_a, scan_b, roi, folder)

# Every pair at once, ordered by the label rather than by scan number.
series = build_series(pairs, roi, folder)
plot_difference_series(series, "NdFeB dichroic difference")
```

`dichro_maps` also carries the corrections the raw maps need: `align_series`
undoes the xeryon stage runout row by row, `stitch_x` joins x-adjacent tiles
using each point only once, `realign_pair` re-derives a difference with the
two polarizations aligned, and `write_gif` animates a series.

:::{note}
`process_flyscan` and `dichro_maps` read HDF5 files written with LZ4 and
bitshuffle compression. They import `hdf5plugin` at module level to register
those codecs — don't remove that import, even though nothing references the
name afterwards.
:::

## See also

- [Scan viewers](viewers.md) — the step-scan browser and its image viewer
- [Image processing](process_images.md) — ROI and curvature reduction for
  step scans
