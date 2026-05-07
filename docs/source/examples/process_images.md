# Image processing

The `polartools.process_images` module handles area-detector images
(Lambda, Eiger, Princeton SPE) and reduces them to spectra. It uses
[dask][dask] for lazy/parallel loading so you can work with large image
stacks without loading them all into memory.

[dask]: https://docs.dask.org/

## Loading a stack of images

```python
from polartools.process_images import load_images

images = load_images(
    scan_id=276,
    source=db,                # databroker catalog
    detector="lambda250k",
)
# images is a dask.array — operations are lazy
print(images.shape)
print(images.dtype)
```

For large scans (thousands of frames), all subsequent operations are chunked
and parallelized.

## Cleaning images: thresholding

Most analysis starts with cleaning out cosmic rays, hot pixels, or below-noise
counts:

```python
from polartools.process_images import clean_threshold

clean = clean_threshold(images, threshold=20)  # zero out anything <= 20 counts
```

## Photon-event extraction (RIXS)

For low-flux RIXS, you typically want to centroid each photon event. The
internal `_cleanup_photon_events` and `_cleanup_images` helpers do this; the
public entry points fold it in.

## Curvature fitting

RIXS spectrometers usually project a curved energy-dispersive line onto the
detector. `get_curvature` fits the centroid line as a polynomial of column
position so subsequent integration follows the actual energy axis:

```python
from polartools.process_images import get_curvature

# Fit on a single integrated image
sum_image = images.sum(axis=0).compute()
curvature = get_curvature(sum_image, binx=10, biny=1, plot=True)
# curvature is a numpy.poly1d that maps column → row offset
```

`binx` / `biny` control the binning before the centroid fit; `plot=True`
shows the fit overlaid on the image (useful for tuning).

## Image → spectrum

Once you have a curvature, integrate along the dispersive direction:

```python
from polartools.process_images import get_spectrum, get_spectra

# Single-image reduction
spectrum = get_spectrum(sum_image, curvature=curvature, biny=1)

# Stack reduction (returns one row per input image)
all_spectra = get_spectra(images, curvature=curvature, biny=1)
```

## Resonant inelastic X-ray scattering (RIXS) workflows

For an energy-loss × incident-energy 2D map (RXES):

```python
from polartools.process_images import process_rxes

energy_in, energy_loss, intensity = process_rxes(
    scan_ids=range(100, 150),
    source=db,
    detector="lambda250k",
    curvature=curvature,
    threshold=20,
)
```

For a magnetic-circular-dichroism RXES (RXES-MCD), where each "frame" comes
in groups of 4 polarization states:

```python
from polartools.process_images import process_rxes_mcd

results = process_rxes_mcd(
    scan_ids=range(100, 150),
    source=db,
    detector="lambda250k",
    curvature=curvature,
    threshold=20,
)
```

> **Note:** `process_rxes_mcd` requires the per-scan image count to be a
> multiple of 4. The bundled `lambda250k` test scan has 50 images so this
> path isn't currently covered by tests (issue [#27]).

[#27]: https://github.com/APS-4ID-POLAR/polartools/issues/27

## Detector handlers

For databroker scans, the relevant asset handlers are registered automatically
by `polartools.load_data.load_catalog`. The `polartools.area_detector_handlers`
module defines:

- `LambdaHDF5Handler` — Lambda 250k frames in HDF5
- `EigerHandler` — Eiger frames in HDF5
- `SPEHandler` — Princeton SPE files

You generally don't need to interact with these directly.

## Performance tips

- Use `images = load_images(...)` and operate on the dask array. Only
  `.compute()` when you need the result in memory.
- For curvature fitting, `binx >= 10` is much faster and rarely changes the
  fit by more than the pixel size.
- If you're hitting memory limits, smaller dask chunks (`images.rechunk(...)`)
  often help more than smaller datasets.

## See also

- [Database management](manage_database.md) — exporting image-bearing scans
