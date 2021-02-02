# Copyright (c) 2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from databroker_pack import (export_catalog, write_documents_manifest,
                             write_msgpack_catalog_file)
from .load_data import run_v2_query
from os import makedirs
from os.path import exists, join
from itertools import tee
from suitcase.utils import MultiFileManager
from suitcase.csv import export as csv_export
from suitcase.json_metadata import export as json_export
from warnings import warn


def to_databroker(db, folder, query=None):
    """
    Exports databroker database into msgpack files.

    WARNING: While you can pass a query dictionary here, it is advised to run
    the query and check the results before running this function as you
    may inadvertely export a very large number of scans. See
    :func:`polartools.load_data.run_v2_query`.

    This is a narrow usage of the databroker-pack package. Note that this
    package includes a convenient command line tool.
    https://blueskyproject.io/databroker-pack/index.html

    Parameters
    ----------
    db :
        Databroker database.
    folder : str
        Destination directory.
    query : dict, optional
        Search parameters to select a subsection of `db`. See
        :func:`polartools.load_data.run_v2_query` for more details.

    Notes:
    ------
    - The scans are saved in msgpack files placed in the `folder/documents`
    folder.
    - `catalog.yml` and `documents_manifest.txt` are located in `folder`.

    See also
    --------
    :func:`polartools.load_data.run_v2_query`
    :func:`databroker-pack.export_catalog`
    :func:`databroker-pack.write_documents_manifest`
    :func:`databroker-pack.write_msgpack_catalog_file`
    """
    results = run_v2_query(db, query) if query else db.v2

    makedirs(folder, exist_ok=True)
    manager = MultiFileManager(folder)

    artifacts, _, _, _ = export_catalog(results, manager)
    write_documents_manifest(manager, folder, artifacts["all"])
    write_msgpack_catalog_file(manager, folder, ["./documents/*.msgpack"], {})


def to_csv_json(db, folder, query=None, fname_format='scan-{}-',
                overwrite=True, max_attempts=100):
    """
    Exports scans into *.csv and *.json files.


    The scans will be labeled by their `scan_id` metadata. If two or more scans
    have the same `scan_id`, it will write the new scan with a `-number` suffix
    where number will be the first available integer starting with 2.

    WARNING: While you can pass a query dictionary here, it is advised to run
    the query and check the results before running this function as you
    may inadvertely export a very large number of scans. See
    :func:`polartools.load_data.run_v2_query`.

    Parameters
    ----------
    db :
        Databroker database.
    folder : str
        Destination directory.
    query : dict, optional
        Search parameters to select a subsection of `db`. See
        :func:`polartools.load_data.run_v2_query` for more details.
    fname_format : str, optional
        Format of the string to be used for the file names. Note that one has
        to be able to add the scan numbers into this string by doing:
        `fname_format.format(scan_number)`.
    overwrite : bool, optional
        Flag to determine if an existing folder should be overwritten.
    max_attempts : int, optional
        Maximum number of times that a new suffix will be added to the file
        name. Once it reaches this maximum, it will overwrite the last file.

    Notes:
    ------
    - Each scan has one json "*-metadata.json" file, plus ome csv file for each
    data stream, for instance "*-primary.csv".

    See also
    --------
    :func:`polartools.load_data.run_v2_query`
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

    results = run_v2_query(db, query) if query else db.v2
    for uid in list(results):
        scanno = results[uid].metadata['start']['scan_id']
        fname = next_available_fname(folder, fname_format.format(scanno),
                                     max_attempts)
        my_exporter(results[uid].documents(fill='yes'), folder, fname)
