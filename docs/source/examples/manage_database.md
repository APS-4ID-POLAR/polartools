# Database management

The `polartools.manage_database` module wraps [`databroker-pack`][dbpack] and
[`suitcase`][suitcase] to import/export Bluesky catalogs. It's how you
- distribute a beamtime's data as a self-contained directory,
- ingest someone else's exported catalog into your local databroker,
- export to plain CSV/JSON for downstream analysis.

[dbpack]: https://blueskyproject.io/databroker-pack/
[suitcase]: https://blueskyproject.io/suitcase-csv/

> **Scope:** these helpers work with **databroker** catalogs only — tiled is
> not yet supported here (issue [#60]).

[#60]: https://github.com/APS-4ID-POLAR/polartools/issues/60

## Export a catalog to portable msgpack

Pack a subset of a catalog into a folder you can copy to another machine:

```python
from polartools.manage_database import to_databroker

to_databroker(
    db,                           # source databroker catalog
    folder="/exports/run3",       # destination
    query={"scan_id": range(100, 200)},
)
```

The destination will contain a `documents/` directory with msgpack files,
a `catalog.yml`, and a `documents_manifest.txt`.

## Re-import a packed catalog

On the receiving machine:

```python
from databroker import catalog
from polartools.manage_database import from_databroker_inplace

from_databroker_inplace(
    folder="/data/run3",          # the directory created by to_databroker
    name="run3",                  # the catalog name to register
    catalog=catalog,
)

# Now use it like any other catalog
db = catalog["run3"]
df = db[100].primary.read().to_dataframe()
```

`from_databroker_inplace` writes a config file under
`~/.local/share/intake/` that points databroker at the packed data.

## CSV / JSON export

For analysis tools that don't speak Bluesky:

```python
from polartools.manage_database import to_csv_json

to_csv_json(
    db,
    folder="/exports/csv",
    query={"sample": "Fe2O3"},
    overwrite=True,
)
```

Each scan becomes a `.csv` file (per primary stream) plus a `.json`
sidecar with the start/stop documents.

> **Known issue:** `to_csv_json` doesn't currently work under pytest because
> of a logging conflict (issue [#56]). It's fine in interactive use.

[#56]: https://github.com/APS-4ID-POLAR/polartools/issues/56

## Removing a registered catalog

```python
from polartools.manage_database import remove_catalog

remove_catalog("run3", catalog=catalog)
```

This removes the catalog *registration* but leaves the underlying data files
on disk — you can re-add them later with `from_databroker_inplace`.

## Catalog with image data

If your scans reference external image files (Lambda HDF5, Eiger HDF5, SPE),
pass `external=True` to also pack the image files alongside the documents:

```python
to_databroker(db, folder="/exports/run3",
              query={"scan_id": range(100, 200)},
              external=True)
```

Without `external=True`, the packed catalog will reference the original image
file paths — works if the destination machine can see them, breaks otherwise.

## See also

- [Loading data](load_data.md) — how to read what you exported
- [Image processing](process_images.md) — workflows for image-bearing scans
