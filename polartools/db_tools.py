"""
Functions to access database information.

.. autosummary::
    ~db_query
    ~show_meta
    ~collect_meta
"""

# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from datetime import datetime
from collections import OrderedDict
from pyRestTable import Table
from warnings import warn
from databroker.queries import TimeRange


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
    scans, db, scan_to=None, query=None, meta_keys="short", table_fmt="plain"
):
    """
    Print metadata of scans.

    Parameters
    ----------
    scans : int or iterable
        Scan numbers or uids. If an integer is passed, it will process scans
        from `scans` to `scan_to`.
    db : databroker database
        Searcheable database
    scan_to : int, optional
        Final scan number to process. Note that this is only meaningful if
        an integer is passed to `scans`.
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

        scans = range(scans, scan_to+1)

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

    meta = collect_meta(scans, db, meta_keys, query=query)
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


def collect_meta(scan_numbers, db, meta_keys, query=None):
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
