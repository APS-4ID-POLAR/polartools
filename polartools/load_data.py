"""
Base functions to load data from various sources.

.. autosummary::
   ~load_spec
   ~load_csv
   ~load_databroker
   ~db_query
   ~load_table
   ~is_Bluesky_specfile
"""

# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from pandas import read_csv, DataFrame
from os.path import join
from spec2nexus.spec import SpecDataFile
from warnings import warn
import copy

try:
    from databroker.queries import TimeRange
except ModuleNotFoundError:
    pass


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


def load_databroker(scan_id, db, stream="primary", query=None, use_db_v1=True):
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

    _db = db_query(db, query) if query else db
    if use_db_v1:
        return _db.v1[scan_id].table(stream_name=stream)
    else:
        return getattr(_db.v2[scan_id], stream).read().to_dataframe()


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


def load_table(scan, source, **kwargs):
    """
    Automated scan table loader.

    The automation is based on the source argument.

    - if source == 'csv' -> uses `load_csv`.
    - else if source is a string or nexus2spec.spec.SpecDataFile -> uses \
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
            scan, source, stream=stream, query=query, use_db_v1=use_db_v1
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
