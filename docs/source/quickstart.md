# Quickstart

A short tour of `polartools`: load a scan, look at it, fit a peak.

## Load a SPEC scan into a DataFrame

```python
from polartools.load_data import load_table

table = load_table(scan_id=10, source="my_spec_file.dat", folder="/data/2024-1")
table.head()
```

`load_table` is the universal entry point — pass it a SPEC file, an HDF5 file,
a CSV folder, a databroker catalog, or a tiled catalog and it returns a pandas
`DataFrame` with the scan's data columns.

## Plot it

```python
import matplotlib.pyplot as plt

ax = table.plot(x="Energy", y="IfluorM/I0", marker=".")
ax.set_xlabel("Energy (eV)")
ax.set_ylabel("normalized fluorescence")
```

## Process an XMCD scan

```python
from polartools.absorption import process_xmcd

xmcd = process_xmcd(
    scan_plus=10,  # +H scan
    scan_minus=11, # -H scan
    source="my_spec_file.dat",
    folder="/data/2024-1",
    monitor="I0",
    detector="IfluorM",
    transmission=False,
)
xmcd.plot(x="energy", y="xmcd")
```

## Fit a Bragg peak

```python
from polartools.diffraction import fit_peak

result = fit_peak(
    scan_id=42,
    source="my_spec_file.dat",
    folder="/data/2024-1",
    x="th",
    y="apd",
    monitor="I0",
    model="PseudoVoigt",
)
print(result.fit_report())
result.plot()
```

## Where to next

- [Load data](examples/load_data.md) — every supported source backend, with examples
- [XAS / XMCD](examples/absorption.md) — normalization, dichroism, edge-step
- [Diffraction](examples/diffraction.md) — peak fitting, peak series, HKL conversions
- [API reference](api/polartools/index.rst) — full function-by-function reference
