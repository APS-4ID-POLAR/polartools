"""
Functions to import and export Bluesky data.

.. autosummary::
   ~to_databroker
   ~to_csv_json
   ~from_databroker_inplace
   ~remove_catalog
"""

# Copyright (c) 2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from databroker_pack import (
    export_catalog, write_documents_manifest, write_msgpack_catalog_file,
    unpack_inplace, copy_external_files, write_external_files_manifest
)

from databroker import catalog_search_path
from .load_data import db_query
from os import makedirs, remove
from os.path import exists, join
from itertools import tee
from suitcase.utils import MultiFileManager
from suitcase.csv import export as csv_export
from suitcase.json_metadata import export as json_export
from warnings import warn
from pathlib import Path


def to_databroker(db, folder, query=None, external=False):
    """
    Exports databroker database into msgpack files.
    WARNING: While you can pass a query dictionary here, it is advised to run
    the query and check the results before running this function as you
    may inadvertely export a very large number of scans. See
    :func:`polartools.load_data.db_query`.
    This is a narrow usage of the `databroker-pack` package_. Note that this
    package includes a convenient command line tool.
    .. _package: https://blueskyproject.io/databroker-pack/index.html
    Parameters
    ----------
    db :
        Databroker database.
    folder : str
        Destination directory.
    query : dict, optional
        Search parameters to select a subsection of `db`. See
        :func:`polartools.load_data.db_query` for more details.
    Notes
    ------
    - The scans are saved in msgpack files placed in the `folder/documents` \
    folder.
    - `catalog.yml` and `documents_manifest.txt` are located in `folder`.
    See also
    --------
    :func:`polartools.load_data.db_query`
    :func:`databroker-pack.export_catalog`
    :func:`databroker-pack.write_documents_manifest`
    :func:`databroker-pack.write_msgpack_catalog_file`
    """
    results = db_query(db, query) if query else db.v2

    makedirs(folder, exist_ok=True)
    manager = MultiFileManager(folder)

    artifacts, external_files, _, _ = export_catalog(results, manager)
    write_documents_manifest(manager, folder, artifacts["all"])

    root_map = {}
    if external:
        target_directory = Path(folder, "external_files")
        for ((_, root, unique_id), files) in external_files.items():
            new_root, new_files, _ = copy_external_files(
                target_directory, root, unique_id, files
            )
            # copying_failures.extend(copying_failures_)
            # The root_map value will be the relative path to
            # the data within the pack directory.
            relative_root = new_root.relative_to(folder)
            root_map.update({unique_id: relative_root})
            rel_paths = [Path(f).relative_to(folder) for f in new_files]
            write_external_files_manifest(manager, unique_id, rel_paths)
    write_msgpack_catalog_file(
        manager, folder, ["./documents/*.msgpack"], root_map
    )


def to_csv_json(db, folder, query=None, fname_format='scan_{}_',
                overwrite=False, max_attempts=100):
    """
    Exports scans into .csv and .json files.


    The scans will be labeled by their `scan_id` metadata. If two or more scans
    have the same `scan_id`, it will write the new scan with a `-number` suffix
    where number will be the first available integer starting with 2.

    WARNING: While you can pass a query dictionary here, it is advised to run
    the query and check the results before running this function as you
    may inadvertely export a very large number of scans. See
    :func:`polartools.load_data.db_query`.

    Parameters
    ----------
    db :
        Databroker database.
    folder : str
        Destination directory.
    query : dict, optional
        Search parameters to select a subsection of `db`. See
        :func:`polartools.load_data.db_query` for more details.
    fname_format : str, optional
        Format of the string to be used for the file names. Note that one has
        to be able to add the scan numbers into this string by doing:
        `fname_format.format(scan_number)`.
    overwrite : bool, optional
        Flag to determine if an existing folder should be overwritten.
    max_attempts : int, optional
        Maximum number of times that a new suffix will be added to the file
        name. Once it reaches this maximum, it will overwrite the last file.

    Notes
    -----
    - Each scan has one json "-metadata.json" file, plus ome csv file for \
    each data stream, for instance "-primary.csv".

    See also
    --------
    :func:`polartools.load_data.db_query`
    :func:`suitcase.csv.export`
    :func:`suitcase.json_metadata.export`
    """
    def my_exporter(docs, directory, file_prefix):
        docs1, docs2 = tee(docs, 2)
        csv_export(docs1, directory, file_prefix)
        json_export(docs2, directory, file_prefix)

    def next_available_fname(folder, fname, max_attempts):
        new_fname = fname
        base = join(folder, fname)
        for i in range(2, max_attempts+3):
            if exists(base + 'baseline.csv'):
                new_fname = fname + f'{i}-'
                base = join(folder, new_fname)
            else:
                return new_fname

        warn(f'Reached maximum number of attempts ({max_attempts})'
             f', the files with "{new_fname}" prefix will be overwritten.')
        return new_fname

    if not overwrite and exists(folder):
        raise FileExistsError(f'{folder} already exists. Either select '
                              '`overwrite=True` or enter a another folder.')

    results = db_query(db, query) if query else db.v2
    for uid in list(results):
        scanno = results[uid].metadata['start']['scan_id']
        fname = next_available_fname(folder, fname_format.format(scanno),
                                     max_attempts)
        print('Exporting uid #{}, scan_id #{}'.format(uid[:8], scanno),
              end='... ')
        my_exporter(results[uid].documents(fill='yes'), folder, fname)
        print('Done!')


def from_databroker_inplace(folder, name, catalog, merge=False):
    """
    Load the exported databroker database.

    This is a narrow usage of the `databroker-pack` package_. Note that this
    package includes a convenient command line tool.

    .. _package: https://blueskyproject.io/databroker-pack/index.html

    Parameters
    ----------
    folder : str
        Folder with files exported by
        :func:`polartools.manage_database.to_databroker`.
    name : str
        Unique name that will be set in the databroker catalog.
    merge : bool, optional
        Flag to decide if this data will be merged into an existing catalog.

        Example ::

            from databroker import catalog
            from polartools.manage_database import from_databroker_inplace

            from_databroker_inplace('folder/to/files', 'my_data')
            db = catalog['my_data']
    """
    config_path = unpack_inplace(folder, name, merge=merge)
    catalog.force_reload()
    print(f"Placed configuration file at {config_path!s}")


def remove_catalog(name, catalog=None):
    """
    Removes a catalog created by `to_databroker`.

    Parameters
    ----------
    name : str
        Catalog name

    Notes
    -----
    - This will not remove the data files, only the catalog, which will not be\
     discoverable using `databroker.catalog`.

    - It assumes that the catalog was created by `to_databroker` or \
    `databroker_pack`. The actual name of the file is \
    'databroker_unpack_NAME.yml'.
    """

    found = False
    for path in catalog_search_path():
        filepath = join(path, f'databroker_unpack_{name}.yml')
        if exists(filepath):
            remove(filepath)
            if catalog is not None:
                catalog.force_reload()
            print(f'The {name} catalog was removed.')
            found = True

    if not found:
        print(f'The catalog {name} was not found.')
