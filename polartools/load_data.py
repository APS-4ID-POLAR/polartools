"""
Base functions to load data from various sources.

.. autosummary::
    ~load_spec
    ~load_csv
    ~load_databroker
    ~load_table
    ~is_Bluesky_specfile
    ~db_query
    ~show_meta
    ~collect_meta
    ~lookup_position
"""

# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from pandas import read_csv, DataFrame
from os.path import join
from spec2nexus.spec import SpecDataFile
from warnings import warn
import copy
from datetime import datetime
from collections import OrderedDict
from pyRestTable import Table
from databroker.queries import TimeRange
from apstools.utils import getDatabase
from polartools.area_detector_handlers import EigerHandler, LambdaHDF5Handler

# TODO: This should be just temp fix
LambdaHDF5Handler.specs = {'AD_HDF5_lambda'} | LambdaHDF5Handler.specs


def load_catalog(name=None, query=None, handlers=None):
    """
    Loads a databroker catalog and register data handlers.

    Parameters
    ----------
    name : str, optional
        Name of the database. Defaults to 4-ID-D name.
    query : dict, optional
        Dictionary with search parameters for the database.
    handlers : dict, optional
        Dictionary organized as {handler_name: handler_class}. If None,
        defaults to handlers used at 4-ID-D.

    Returns
    -------
    cat : databroker catalog
        Catalog after running the query, and registering the handler.
    """
    cat = getDatabase(catalog_name=name)
    if query is not None:
        cat = db_query(cat, query)
    if handlers is None:
        handlers = dict(
                AD_HDF5_Lambda250k_APSPolar=LambdaHDF5Handler,
                AD_HDF5_lambda=LambdaHDF5Handler,  # Temporary fix
                AD_EIGER_APSPolar=EigerHandler,
            )
    for name, handler in handlers.items():
        cat.register_handler(name, handler, overwrite=True)
    return cat


def load_spec(scan_id, spec_file, folder=""):
    """
    Load data from spec file.

    If `spec_file` is the file name, it will load the spec file internally
    which is time consuming.

    Parameters
    ----------
    scan_id : int
        Scan_id of the scan to be retrieved.
    spec_file : string or spec2nexus.spec.SpecDataFile
        Either the spec file name or a SpecDataFile instance.
    folder : string, optional
        Folder where spec file is located.

    Returns
    -------
    data : pandas.DataFrame
        Table with the data from scan.

    See also
    --------
    :func:`spec2nexus.spec.SpecDataFile`
    """

    if isinstance(spec_file, str):
        path = join(folder, spec_file)
        spec_file = SpecDataFile(path)

    return DataFrame(spec_file.getScan(scan_id).data)


def load_csv(scan_id, folder="", name_format="scan_{}_primary.csv"):
    """
    Load data from the 'primary' stream from exported csv files.

    Parameters
    ----------
    scan_id : int
        Scan_id of the scan to be retrieved.
    folder : string, optional
        Folder where csv files are located.
    name_format : string, optional
        General format of file name. The correct name must be
        retrievable through: `file_name_format.format(scan_id)`

    Returns
    -------
    data : pandas.DataFrame
        Table with the data from the primary stream.

    See also
    --------
    :func:`pandas.read_csv`
    """

    return read_csv(join(folder, name_format.format(scan_id)))


def load_databroker(
    scan_id, db=None, stream="primary", query=None, use_db_v1=True
):
    """
    Load data of the first scan with the provided scan_id.

    Currently defaults to databroker.v1 because it is faster. See issue #28.

    For further details, refer to the `databroker` `documentation`_.

    .. _documentation: https://blueskyproject.io/databroker/

    Parameters
    ----------
    scan_id : int
        Scan_id of the scan to be retrieved
    db :
        `databroker` database
    stream : string, optional
        Selects the stream from which data will be loaded.
    query : dict, optional
        Dictionary with search parameters for the database.
    use_db_v1 : bool, optional
        Chooses databroker API version between 'v1' or 'v2', defaults to 'v1'.

    Returns
    -------
    data : pandas.DataFrame
        Table with the data from the primary stream.
    """
    db = getDatabase(db=db)
    _db = db_query(db, query) if query else db
    if use_db_v1:
        if stream in _db.v1[scan_id].stream_names:
            return _db.v1[scan_id].table(stream_name=stream)
        else:
            raise ValueError(
                f"The stream {stream} does not exist in scan {scan_id}."
            )
    else:
        try:
            return getattr(_db.v2[scan_id], stream).read().to_dataframe()
        except AttributeError:
            raise ValueError(
                f"The stream {stream} does not exist in scan {scan_id}."
            )


def load_table(scan, source=None, **kwargs):
    """
    Automated scan table loader.

    The automation is based on the source argument.

    - if source == 'csv' -> uses `load_csv`.
    - else if source is a string or nexus2spec.spec.SpecDataFile -> uses\
    `load_spec`.
    - else -> uses `load_databroker`.

    Parameters
    ----------
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query, use_db_v1.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

    Returns
    -------
    table : pandas.DataFrame
        Table with the scan data.

    See also
    --------
    :func:`polartools.load_data.load_databroker`
    :func:`polartools.load_data.load_csv`
    :func:`polartools.load_data.load_spec`
    """

    _kwargs = copy.deepcopy(kwargs)
    folder = _kwargs.pop("folder", "")
    if source == "csv":
        name_format = kwargs.pop("name_format", "scan_{}_primary.csv")
        table = load_csv(scan, folder=folder, name_format=name_format)
    elif isinstance(source, str) or isinstance(source, SpecDataFile):
        table = load_spec(scan, source, folder=folder)
    else:
        stream = _kwargs.pop("stream", "primary")
        query = _kwargs.pop("query", None)
        use_db_v1 = _kwargs.pop("use_db_v1", True)
        table = load_databroker(
            scan, db=source, stream=stream, query=query, use_db_v1=use_db_v1
        )

    if len(_kwargs) != 0:
        warn(f"The following kwargs were not used! {list(_kwargs.keys())}")

    return table


def is_Bluesky_specfile(source, folder=""):
    """
    Check if the specfile was created by Bluesky.

    It looks for a "Bluesky" comment in the file header.

    Parameters
    ----------
    source : string or spec2nexus.spec.SpecDataFile
        Either the spec file name or a SpecDataFile instance.
    folder : string, optional
        Folder where spec file is located.

    Returns
    -------
    y : bool
        True if spec_file was generated by Bluesky, False otherwise.

    See also
    --------
    :func:`spec2nexus.spec.SpecDataFile`
    """
    if isinstance(source, str):
        path = join(folder, source)
        source = SpecDataFile(path)

    if len(source.headers[0].comments):
        return source.headers[0].comments[0].startswith("Bluesky")
    else:
        return False


def db_query(db, query):
    """
    Searches the databroker v2 database.

    Parameters
    ----------
    db :
        `databroker` database.
    query: dict
        Search parameters.

    Returns
    -------
    _db :
        Subset of db that satisfy the search parameters. Note that it has the
        same format as db.

    See also
    --------
    :func:`databroker.catalog.search`
    """

    since = query.pop("since", None)
    until = query.pop("until", None)

    if since or until:
        if not since:
            since = "2010"
        if not until:
            until = "2050"

        _db = db.v2.search(TimeRange(since=since, until=until))
    else:
        _db = db

    if len(query) != 0:
        _db = _db.v2.search(query)

    return _db


def show_meta(
    scans,
    scan_to=None,
    db=None,
    query=None,
    meta_keys="short",
    table_fmt="plain",
):
    """
    Print metadata of scans.

    Parameters
    ----------
    scans : int or iterable
        Scan numbers or uids. If an integer is passed, it will process scans
        from `scans` to `scan_to`.
    scan_to : int, optional
        Final scan number to process. Note that this is only meaningful if
        an integer is passed to `scans`.
    db : databroker database (optional)
        Searcheable database
    query : dictionary, optional
        Search parameters.
    meta_keys : string or iterable, optional
        List with metadata keys to read. There are two preset metadata lists
        that can be used with `meta_keys="short"` or `meta_keys="long"`.
    """

    if isinstance(scans, int):
        if scan_to is None:
            scan_to = scans

        if scan_to < scans:
            raise ValueError(
                "scans must be larger than scan_to, but you "
                f"entered: scans = {scans} and scan_to = "
                f"{scan_to}."
            )

        scans = range(scans, scan_to + 1)

    if meta_keys == "short":
        meta_keys = [
            "motors",
            "scan_type",
            "plan_name",
            "plan_pattern_args",
            "num_points",
            "exit_status",
        ]
    elif meta_keys == "long":
        meta_keys = [
            "motors",
            "scan_type",
            "plan_name",
            "plan_pattern_args",
            "num_points",
            "exit_status",
            "time",
            "hints",
        ]

    meta = collect_meta(scans, meta_keys, db=db, query=query)
    table = Table()

    if "plan_pattern_args" in meta_keys:
        index = meta_keys.index("plan_pattern_args")
        meta_keys.remove("plan_pattern_args")
        meta_keys.insert(index, "final pos.")
        meta_keys.insert(index, "init. pos.")

    table.labels = ["Scan # "] + list(meta_keys)
    for scanno, values in meta.items():
        row = [scanno]
        for key, item in values.items():
            # TODO: I don't like this. We need better metadata.
            if key == "plan_pattern_args":
                if item[0] is None:
                    row.append(None)
                    row.append(None)
                elif isinstance(item[0]["args"][1], list):
                    row.append("{:0.4f}".format(item[0]["args"][1][0]))
                    row.append("{:0.4f}".format(item[0]["args"][1][-1]))
                else:
                    row.append(item[0]["args"][-2])
                    row.append(item[0]["args"][-1])

            elif key == "time":
                time = []
                for t in item:
                    time.append(
                        datetime.fromtimestamp(t).strftime("%m/%d/%Y %H:%M:%S")
                    )
                row.append(time)

            else:
                if None in item:
                    item.remove(None)
                if len(item) == 1:
                    item = item[0]
                row.append(item)
        table.rows.append(row)

    print(table.reST(fmt=table_fmt))


def collect_meta(scan_numbers, meta_keys, db=None, query=None):
    """
    Extracts metadata of a list of scans.

    Parameters
    ----------
    scans : iterable
        Scan numbers or uids.
    db : databroker database
        Searcheable database
    scan_to : int, optional
        Final scan number to process. Note that this is only meaningful if
        an integer is passed to `scans`.
    meta_keys : iterable
        List with metadata keys to read.
    query : dictionary, optional
        Search parameters.

    Returns
    -------
    meta : dictionary
        Metadata organized by scan number or uid (whatever is given in
        `scans`).
    """
    db = getDatabase(db=db)
    # print(f"db = {db}")
    # print(f"scan_numbers = {scan_numbers}")
    db_range = db_query(db, query=query) if query else db
    output = OrderedDict()
    for scan in scan_numbers:
        try:
            start = db_range[scan].metadata["start"]
            stop = db_range[scan].metadata.get("stop", None)

            output[scan] = OrderedDict()
            for key in meta_keys:
                output[scan][key] = [start.get(key, None)]
                if stop is None:
                    output[scan][key] += [None]
                else:
                    output[scan][key] += [stop.get(key, None)]
                output[scan][key] = _flatten_list(output[scan][key])

        except KeyError:
            warn(f"The scan number {scan} was not found.")

    return output


def _flatten_list(inp):
    """Only handles lists to second level"""
    output = []
    for item in inp:
        if isinstance(item, list):
            output += [i for i in item]
        else:
            output.append(item)
    return output


def lookup_position(db, scan, search_string="", query=None):
    """
    Lookup positioner values in past scans.


    Parameters
    ----------

    db : databroker database
        Searcheable database
    scan : integer
        Scan numbers or uids.
    search_string : string
        Full or part of positioner name.
    query: dict
        Search parameters.


    Returns
    -------
    output: list
    """

    db_range = db_query(db, query=query) if query else db

    baseline = load_databroker(scan, db_range, "baseline", use_db_v1=True)
    if len(baseline["time"]) == 2:
        date1 = baseline["time"][1].strftime("%m/%d/%y %H:%M:%S")
        date2 = baseline["time"][2].strftime("%m/%d/%y %H:%M:%S")
        print("=".center(100, "="))
        print(f"{'Positioner':>50}{date1:>25}{date2:>25}")
        print("-".center(100, "-"))
        for key in baseline.keys():
            if search_string in key:
                if isinstance(baseline[key][1], list):
                    print(f"{key:>50}{baseline[key][1]}{baseline[key][2]}")
                else:
                    print(
                        f"{key:>50}{baseline[key][1]:>25}{baseline[key][2]:>25}"
                    )

    else:
        date1 = baseline["time"][1].strftime("%m/%d/%y %H:%M:%S")
        print("=".center(100, "="))
        print(f"{'Positioner':>50}{date1:>25}")
        print("-".center(100, "-"))
        for key in baseline.keys():
            if search_string in key:
                if isinstance(baseline[key][1], list):
                    print(f"{key:>50}{baseline[key][1]}")
                else:
                    print(f"{key:>50}{baseline[key][1]:>25}")

    print("-".center(100, "-"))
