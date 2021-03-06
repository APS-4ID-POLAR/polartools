"""
Functions to access database information.

.. autosummary::
    ~show_meta
"""

from datetime import datetime
from polartools.load_data import db_query
from collections import OrderedDict
from pyRestTable import Table
from warnings import warn


def show_meta(db, scan_from, scan_to=None, query=None, long=False):
    """
    Searches the databroker v2 database.
    Parameters
    ----------
    db :
        `databroker` database.
    from: int
        Scan number start
    to: int
        Scan number end
    query: dict
        Search parameters.
    long: boolean
        long output
    Returns
    -------
    output :
        List of scan metadata
    """

    db_range = db_query(db, query=query) if query else db
    time_old = 0
    if scan_to and scan_from < scan_to:
        scan_to += 1
    else:
        scan_to = scan_from + 1
    for scan_number in range(scan_from, scan_to):
        try:
            meta = db_range[scan_number].metadata["start"]

        except KeyError:
            print(f"Scan {scan_number} not existing!")
            return

        time_new = meta["time"]
        time = datetime.fromtimestamp(time_new)
        if time_new < time_old:
            return
        time_old = time_new
        scan_type = meta["plan_name"] if "plan_name" in meta else None
        motors = meta["motors"] if "motors" in meta else None

        plan_args = (
            meta["plan_pattern_args"] if "plan_pattern_args" in meta else None
        )

        if scan_type == "list_scan":
            scan_from = plan_args["args"][1][0]
            scan_to = plan_args["args"][1][-1]
            scan_type2 = (
                meta["hints"]["scan_type"]
                if "scan_type" in meta["hints"]
                else None
            )
        else:
            scan_from = plan_args["args"][1]
            scan_to = plan_args["args"][2]

        number_of_points = meta["num_points"]
        hint = meta["hints"]
        det = hint["detectors"] if "detectors" in hint else None
        mon = hint["monitor"] if "monitor" in hint else None
        meta_stop = db[scan_number].metadata["stop"]
        status = (
            meta_stop["exit_status"]
            if "exit_status" in meta_stop
            else "running"
        )

        if scan_type == "list_scan":
            print(
                f"#{scan_number} {scan_type}[{scan_type2}] {motors} "
                f"{scan_from} {scan_to} {number_of_points} {status}"
            )
        else:
            print(
                f"#{scan_number} {scan_type} {motors} {scan_from} {scan_to} "
                f"{number_of_points} {status}"
            )

        if long:
            print(f"        detector={det}, monitor={mon}")
            print(f"        time: {time}")


def show_meta_2(scans, db, scan_to=None, query=None, meta_keys='short'):
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
            scan_to = scans + 1

        if scan_to < scans:
            raise ValueError("scans must be larger than scan_to, but you "
                             f"entered: scans = {scans} and scan_to = "
                             f"{scan_to}.")

        scans = range(scans, scan_to)

    if meta_keys == 'short':
        meta_keys = ["motors", "scan_type", "plan_name", "plan_pattern_args",
                     "num_points", "exit_status"]
    elif meta_keys == 'long':
        meta_keys = ["motors", "scan_type", "plan_name", "plan_pattern_args",
                     "num_points", "exit_status", "time", "hints"]

    meta = collect_meta(scans, db, meta_keys, query=query)
    table = Table()

    if "plan_pattern_args" in meta_keys:
        index = meta_keys.index("plan_pattern_args")
        meta_keys.remove("plan_pattern_args")
        meta_keys.insert(index, "final pos.")
        meta_keys.insert(index, "init. pos.")

    table.labels = ['Scan #'] + list(meta_keys)

    for scanno, values in meta.items():

        row = [scanno]
        for key, item in values.items():
            if key == "plan_pattern_args":
                if item[0] is not None:
                    row.append(item[0]['args'][-2])
                    row.append(item[0]['args'][-1])
                else:
                    row.append(None)
                    row.append(None)

            elif key == 'time':
                time = []
                for t in item:
                    time.append(datetime.fromtimestamp(t).strftime(
                        "%m/%d/%Y %H:%M:%S"))
                row.append(time)

            else:
                if None in item:
                    item.remove(None)
                if len(item) == 1:
                    item = item[0]
                row.append(item)
        table.rows.append(row)

    print(table.reST(fmt='grid'))


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
            start = db_range[scan].metadata['start']
            stop = db_range[scan].metadata.get('stop', None)

            output[scan] = OrderedDict()
            for key in meta_keys:
                output[scan][key] = [start.get(key, None), stop.get(key, None)]
                output[scan][key] = _flatten_list(output[scan][key])

        except KeyError:
            warn(f'The scan number {scan} was not found.')

    return output


def _flatten_list(input):
    """Only handles lists to second level"""
    output = []
    for inp in input:
        if isinstance(inp, list):
            output += [item for item in inp]
        else:
            output.append(inp)
    return output
