# XAS / XMCD

The `polartools.absorption` module covers loading and processing of X-ray
absorption spectroscopy (XAS) and X-ray magnetic circular dichroism (XMCD)
scans. It handles edge-step normalization, pre/post-edge background
subtraction, fluorescence over-absorption correction, and the dichroic
combination of ±-field scans.

## Loading single absorption scans

```python
from polartools.absorption import load_absorption

energy, mu = load_absorption(
    scan_id=10,
    source="myfile.dat",
    folder="/data/2024-1",
    monitor="I0",
    detector="It",       # transmission detector
    transmission=True,   # μ = -log(It / I0)
)
```

For fluorescence:

```python
energy, mu = load_absorption(
    scan_id=10,
    source="myfile.dat",
    folder="/data/2024-1",
    monitor="I0",
    detector="IfluorM",  # fluorescence detector
    transmission=False,  # μ = IfluorM / I0
)
```

## Loading lockin / dichroism scans

For lockin-detected XMCD on a single scan (modulated polarization):

```python
from polartools.absorption import load_lockin

energy, xmcd, xas = load_lockin(
    scan_id=12,
    source="myfile.dat",
    folder="/data/2024-1",
    monitor="I0",
    lockin="lockin",
)
```

## Edge-step normalization

```python
from polartools.absorption import normalize_absorption

mu_norm = normalize_absorption(
    energy,
    mu,
    e0=7112,         # Fe K-edge
    pre1=-150, pre2=-30,
    post1=80,  post2=300,
    post_order=2,
)
```

`normalize_absorption` fits a linear pre-edge and a polynomial post-edge,
subtracts the pre-edge, and divides by the edge step. Returns the normalized
μ on the same energy grid.

If you need the backgrounds separately:

```python
from polartools.absorption import (
    pre_edge_background,
    post_edge_background,
    post_edge_flatten,
)

pre = pre_edge_background(energy, mu, e0=7112, pre1=-150, pre2=-30)
post = post_edge_background(energy, mu, e0=7112, post1=80, post2=300, post_order=2)
mu_flat = post_edge_flatten(energy, mu - pre, e0=7112, post1=80, post2=300, post_order=2)
```

## Fluorescence over-absorption correction

For thick concentrated samples in fluorescence, self-absorption distorts the
spectrum. The hoisted `fluo_corr` (from xraylarch's FLUO algorithm) corrects
it given the chemical formula, edge, and incidence/exit angles:

```python
from polartools.absorption import fluo_corr

mu_corrected = fluo_corr(
    norm=mu_norm,         # edge-step-normalized
    formula="Fe2O3",
    elem="Fe",
    edge="K",
    line="Ka",
    anginp=45,            # degrees from sample surface
    angout=45,
)
```

## XMCD: combining +H / −H scans

```python
from polartools.absorption import process_xmcd, plot_xmcd, save_xmcd

# Process a paired +H / -H scan into XAS and XMCD
plus, minus = process_xmcd(
    scan_plus=10,
    scan_minus=11,
    source="myfile.dat",
    folder="/data/2024-1",
    monitor="I0",
    detector="IfluorM",
    transmission=False,
    e0=7112,
    pre1=-150, pre2=-30,
    post1=80, post2=300,
)

# Quick look
plot_xmcd(plus, minus)

# Save to a header-stamped column file
save_xmcd(plus, minus, "Fe2O3_xmcd.dat",
          header="# sample=Fe2O3, T=10K, B=±1T\n")
```

`process_xmcd` returns two `pandas.DataFrame`s (one per field direction) with
both the normalized XAS and the dichroism column, ready to plot or stack.

## Bulk processing many scans

`load_multi_xas`, `load_multi_dichro`, and `load_multi_lockin` accept lists
of scan IDs and return a single concatenated DataFrame interpolated onto a
common energy grid. Useful for averaging:

```python
from polartools.absorption import load_multi_xas

bulk = load_multi_xas(
    scan_ids=[10, 11, 12, 13, 14],
    source="myfile.dat",
    folder="/data/2024-1",
    monitor="I0",
    detector="IfluorM",
    transmission=False,
)

mean_mu = bulk.iloc[:, 1:].mean(axis=1)  # average across scans
```

## See also

- [`xmcd-gui`](xmcd_gui.md) — interactive PyQt GUI for the XMCD workflow above
- [Loading data](load_data.md) — how `source=` is resolved
