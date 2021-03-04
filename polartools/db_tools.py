"""
Functions to access database information.

.. autosummary::
    ~show_meta
"""

from datetime import datetime
from .load_data import db_query


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
        except:
            print(f"Scan {scan_number} not existing!")
            return
        try:
            scan_type = meta["plan_name"]
        except:
            scan_type = None
        try:
            motors = meta["motors"]
        except:
            motors = None
        try:
            plan_args = meta["plan_pattern_args"]
        except:
            plan_args = None
        if scan_type == "list_scan":
            scan_from = plan_args["args"][1][0]
            scan_to = plan_args["args"][1][-1]
            try:
                scan_type2 = meta["hints"]["scan_type"]
            except:
                scan_type2 = None
        else:
            scan_from = plan_args["args"][1]
            scan_to = plan_args["args"][2]

        number_of_points = meta["num_points"]
        hint = meta["hints"]
        try:
            det = hint["detectors"]
            mon = hint["monitor"]
        except:
            det = None
            mon = None
        time_new = meta["time"]
        time = datetime.fromtimestamp(time_new)
        if time_new < time_old:
            return
        time_old = time_new
        try:
            status = db[scan_number].metadata["stop"]["exit_status"]
        except:
            status = "running"

        if scan_type == "list_scan":
            print(
                f"#{scan_number} {scan_type}[{scan_type2}] {motors} {scan_from} {scan_to} {number_of_points} {status}"
            )
        else:
            print(
                f"#{scan_number} {scan_type} {motors} {scan_from} {scan_to} {number_of_points} {status}"
            )

        if long:
            print(f"        detector={det}, monitor={mon}")
            print(f"        time: {time}")
