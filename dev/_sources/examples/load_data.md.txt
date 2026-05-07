# Loading data

The `polartools.load_data` module is the universal entry point for reading
scans into `pandas.DataFrame`s, regardless of how they were stored on disk.
Five backends are supported:

- **SPEC** files (`.dat`, `.spec`) — text format
- **HDF5** "master" files — Bluesky's HDF5 export
- **CSV** folders — Bluesky's CSV/JSON export (via `suitcase-csv`)
- **databroker** catalogs — MongoDB-backed, the legacy Bluesky catalog
- **tiled** catalogs — postgres-backed, the new Bluesky catalog

The high-level helper {py:func}`polartools.load_data.load_table` auto-routes
to the right backend.

## Auto-routed loading

```python
from polartools.load_data import load_table

# SPEC file
df = load_table(scan_id=10, source="myfile.dat", folder="/data/2024-1")

# HDF5 master file
df = load_table(scan_id=25, source="hdf", folder="/data/run3")

# CSV folder
df = load_table(scan_id=1049, source="csv", folder="/data/csv_export")

# Databroker (legacy)
from databroker import catalog
db = catalog["my_catalog"]
df = load_table(scan_id=1049, source=db)

# Tiled (new)
from tiled.client import from_profile
cat = from_profile("polar")["raw"]
df = load_table(scan_id=2025, source=cat)
```

## Choosing a backend explicitly

If you know your source, the backend-specific functions are slightly faster
and let you pass backend-specific options:

```python
from polartools.load_data import load_spec, load_hdf5_data, load_csv

df = load_spec(scan_id=10, spec_file="myfile.dat", folder="/data/2024-1")
df = load_hdf5_data(scan=25, folder="/data/run3")
df = load_csv(scan_id=1049, folder="/data/csv_export")
```

## Catalog backends side-by-side

`load_catalog` picks the right backend automatically from the catalog name —
it tries `databroker.catalog[name]` first, then falls back to a tiled profile:

```python
from polartools.load_data import load_catalog

cat = load_catalog("polar")          # tries databroker, then tiled
cat = load_catalog("polar",          # tiled-only with explicit subpath
                   tiled_path="/raw/2024-1")
```

For tiled, the default `tiled_path` is `"/raw"`. Pass a deeper path if your
profile is laid out differently.

## Discovering scans

```python
from polartools.load_data import db_query, show_meta, collect_meta

# Query by metadata
runs = db_query(db, query={"sample": "Fe2O3", "scan_type": "ascan"})

# Inspect what metadata is on a scan
show_meta(scan_id=42, db=db, meta_keys=["sample", "user", "uid"])

# Aggregate the same metadata across many scans
table = collect_meta(scan_numbers=range(40, 60), meta_keys=["sample", "energy"], db=db)
```

## Lookup a motor position

```python
from polartools.load_data import lookup_position

th = lookup_position(db, scan=42, search_string="th")
```

Returns the position of any motor whose name contains `"th"` at the start of
the scan.

## Notes

- All loaders return a `pandas.DataFrame` with one row per scan point.
- For SPEC files, the first call parses and caches the file; subsequent calls
  for other scans in the same file are fast.
- For databroker / tiled, the asset handlers (Lambda, Eiger, SPE detectors)
  are registered automatically by `load_catalog`.
- Known issue: `load_multi_*` functions in `absorption` use linear interpolation
  between scans (issue [#44]).

[#44]: https://github.com/APS-4ID-POLAR/polartools/issues/44
