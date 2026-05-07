# Diffraction

The `polartools.diffraction` module fits Bragg peaks and helps with reciprocal-
space data. Peak shapes are powered by [`lmfit`][lmfit] (Gaussian, Lorentzian,
Pseudo-Voigt). The module also includes 2D plotting helpers for mesh scans
and a databroker-aware `dbplot` shortcut.

[lmfit]: https://lmfit.github.io/lmfit-py/

## Fitting a single Bragg peak

```python
from polartools.diffraction import fit_peak, Model
from polartools.load_data import load_table

df = load_table(scan_id=42, source="myfile.dat", folder="/data/2024-1")

result = fit_peak(
    xdata=df["th"],
    ydata=df["apd"] / df["I0"],     # normalize by monitor
    model=Model.PseudoVoigt,
)

print(result.fit_report())
result.plot()
```

`Model` is an enum: `Model.Gaussian`, `Model.Lorentzian`, `Model.PseudoVoigt`.
The result is an `lmfit.model.ModelResult`.

## Loading scan metadata

```python
from polartools.diffraction import load_info, get_type

info = load_info(scan_id=42, info=["sample", "T", "B"],
                 source="myfile.dat", folder="/data/2024-1")
scan_type = get_type(scan_id=42, source="myfile.dat", folder="/data/2024-1")
```

## Fitting a series of scans

`fit_series` fits the same peak shape across a range of scans (e.g. a
temperature or field series) and returns a DataFrame of fit parameters:

```python
from polartools.diffraction import fit_series

# scan_series: [start1, end1, step1, start2, end2, step2, ...]
fits = fit_series(
    scan_series=[100, 119, 1],
    source="myfile.dat",
    folder="/data/2024-1",
    x="th", y="apd", monitor="I0",
    model=Model.PseudoVoigt,
)

# fits is a DataFrame indexed by scan_id with columns:
#   center, sigma, amplitude, height, fwhm, ...
fits.plot(y="center", marker="o")
```

`load_series` returns the raw concatenated data without fitting if you want
to plot or process the series differently.

## Mesh scans (2D maps)

For mesh scans (typically motor1 × motor2 with intensity at each grid point):

```python
from polartools.diffraction import load_mesh, plot_2d

mesh = load_mesh(
    scan_id=200,
    source="myfile.dat",
    folder="/data/2024-1",
    x="h", y="k", intensity="apd",
    monitor="I0",
)

plot_2d(mesh)  # heatmap with axis labels
```

For paired ±-field mesh scans:

```python
from polartools.diffraction import load_dichromesh
mesh_plus, mesh_minus, dichro = load_dichromesh(scan_plus=200, scan_minus=201, ...)
```

## Quick plotting

`load_axes` infers sensible defaults for the x/y/monitor axes from the scan
metadata; `plot_data` and `plot_fit` give one-line plots. `dbplot` is the
databroker shortcut:

```python
from polartools.diffraction import plot_data, plot_fit, dbplot

plot_data(scan_id=42, source="myfile.dat", folder="/data/2024-1")
plot_fit(scan_id=42, source="myfile.dat", folder="/data/2024-1",
         model=Model.PseudoVoigt)

# straight from a databroker catalog
dbplot(db, scan_id=42)
```

## Notes

- Peak fitting expects a single peak per scan. For multi-peak fits, build the
  composite model directly with `lmfit`.
- `monitor` should be a column name; the monitor is divided into the
  intensity column before fitting/plotting.
- `fit_series` accepts the SPEC-style `[start, end, step]` triplets so you
  can stitch together discontiguous scan ranges.
