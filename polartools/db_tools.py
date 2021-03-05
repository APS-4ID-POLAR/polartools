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
                f"#{scan_number} {scan_type}[{scan_type2}] {motors} {scan_from} {scan_to} {number_of_points} {status}"
            )
        else:
            print(
                f"#{scan_number} {scan_type} {motors} {scan_from} {scan_to} {number_of_points} {status}"
            )

        if long:
            print(f"        detector={det}, monitor={mon}")
            print(f"        time: {time}")
