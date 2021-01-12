# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

import numpy as np
from pandas import read_csv, DataFrame
from scipy.interpolate import interp1d
from os.path import join
from spec2nexus.spec import SpecDataFile

try:
    from databroker.queries import TimeRange
except ModuleNotFoundError:
    pass

# TODO: For the beamline we can import db directly, and not worry about it
# all over the place here.


def load_spec(scan_id, file_name, folder=''):
    """
    Load data from spec file.

    Parameters
    ----------
    scan_id : int
        Scan_id of the scan to be retrieved.
    file_name : string
        Name of spec file.
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

    path = join(folder, file_name)
    data = SpecDataFile(path).getScan(scan_id).data
    return DataFrame(data)


def load_csv(scan_id, folder='', file_name_format='scan_{}_primary.csv'):
    """
    Load data from the 'primary' stream from exported csv files.

    Parameters
    ----------
    scan_id : int
        Scan_id of the scan to be retrieved.
    folder : string, optional
        Folder where csv files are located.
    file_name_format : string, optional
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

    return read_csv(join(folder, file_name_format.format(scan_id)))


def load_databroker(scan_id, db, version='v2', stream='primary', **kwargs):
    """
    Load data of the first scan with the provided scan_id.

    Parameters
    ----------
    scan_id : int
        Scan_id of the scan to be retrieved
    db :
        `databroker` database
    version : string, optional
        Version of the databroker to be used in retrieving the data. Must be
        'v1' or 'v2'.
    stream : string, optional
        Selects the stream from which data will be loaded.
    kwargs : optional
        These are passed to the database search. Note that if no kwarg is
        passed, it will skip the search to save time.

    Returns
    -------
    data : pandas.DataFrame
        Table with the data from the primary stream.
    """

    if version == 'v1':
        if len(kwargs) == 0:
            data = db.v1[scan_id].table()
        else:
            data = list(db.v1(scan_id=scan_id, **kwargs))[0].table()
    elif version == 'v2':
        if len(kwargs) == 0:
            _db = db
        else:
            _db = run_v2_query(db, kwargs)
        data = getattr(_db.v2[scan_id], stream).read().to_dataframe()
    else:
        raise ValueError(f"version must be 'v1' or 'v2', but you entered \
            {version}")

    return data


def run_v2_query(db, query):
    """
    Searches the V2 database.

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

    since = query.pop('since', None)
    until = query.pop('until', None)

    if since or until:
        if not since:
            since = '2010'
        if not until:
            until = '2050'

        _db = db.v2.search(TimeRange(since=since, until=until))
    else:
        _db = db

    if len(query) != 0:
        _db = _db.v2.search(query)

    return _db


def load_table(scan, db, file_format='spec', **kwargs):
    """
    Load generic pandas dataframe from one scan.

    Parameters
    ----------
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    db : database
        Databroker database. If None, it will attempt to read from spec or csv
        sfiles.
    file_format : string
        If db = None, then this selects the type of file to open. Options are
        'spec' or 'csv'.
    kwargs:
        If db = None these ware passed to `load_csv`, otherwise passed to
        `load_databroker`.

    Returns
    -------
    (x, y1, y2, ..., yn) : tuple
        The size `n` will be `len(detectors)+1`. The first item is the
        positioner values, and the remaining follow the same order as
        `detectors`.

    See also
    --------
    :func:`polartools.load_data.load_bluesky`
    :func:`polartools.load_data.load_csv`
    """

    if not db:
        if file_format == 'csv':
            return load_csv(scan, **kwargs)
        elif file_format == 'spec':
            file_name = kwargs.pop('file_name', None)
            if not file_name:
                raise NameError("file_name kwarg is needed to load spec files")
            return load_spec(scan, file_name, **kwargs)
        else:
            raise ValueError(f"text_format can only be 'csv' or 'spec', but \
                {file_format} was entered.")
    else:
        return load_databroker(scan, db, **kwargs)


def load_scan(db, scan, positioner, detectors, monitor=None, **kwargs):
    """
    Load generic data of multiple detectors from one scan.

    Parameters
    ----------
    db : database
        Databroker database. If None, it will attempt to read from csv or spec
        files.
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    positioner : string
        Name of the positioner, this needs to be the same as defined in
        Bluesky.
    detectors : iterable
        List of detector to be read from this scan, again it needs to be the
        same name as in Bluesky.
    monitor : string, optional
        Name of the monitor detector. The returned scans will be
        detectors/monitor.
    kwargs :
        Passed to `load_table`.

    Returns
    -------
    x : numpy.array
        Positioner values.

    y1, y2, ..., yn : numpy.array
        The size `n` will be `len(detectors)`. It follows the same order as
        `detectors`.

    See also
    --------
    :func:`polartools.load_data.load_table`
    :func:`polartools.load_data.load_bluesky`
    :func:`polartools.load_data.load_csv`
    """

    table = load_table(scan, db=db, **kwargs)
    data = [np.array(table[positioner])]
    if monitor is None:
        for detector in detectors:
            data.append(np.array(table[detector]))
    else:
        for detector in detectors:
            data.append(np.array(table[detector])/np.array(table[monitor]))

    return tuple(data)


def load_absorption(db, scan, positioner='monochromator_energy',
                    detector='Ion Ch 5', monitor='Ion Ch 4',
                    transmission=True, **kwargs):
    """
    Load the x-ray absorption from one scan.

    Parameters
    ----------
    db : database
        Databroker database. If None, it will attempt to read from csv files.
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky. Defauts to the x-ray energy.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. Defaults to the ion chamber 5.
    monitor : string, optional
        Name of the monitor detector. Defaults to the ion chamber 4.
    transmission: boolean, optional
        Flag to select between transmission mode -> ln(monitor/detector)
        or fluorescence mode -> detector/monitor
    kwargs :
        Passed to `load_scan`.

    Returns
    -------
    x : numpy.array
        Positioner data.
    y : numpy.array
        X-ray absorption.

    See also
    --------
    :func:`polartools.load_data.load_scan`
    """

    x, y = load_scan(db, scan, positioner, [detector], monitor=monitor,
                     **kwargs)

    if transmission:
        return x, np.log(1/y)
    else:
        return x, y


def load_lockin(db, scan, positioner='monochromator_energy', dc_col='Lock DC',
                ac_col='Lock AC', acoff_col='Lock AC off', **kwargs):
    """
    Load the x-ray magnetic dichroism from one scan taken in lock-in mode.

    Parameters
    ----------
    db : database
        Databroker database. If None, it will attempt to read from csv files.
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky. Defauts to the x-ray energy.
    dc_col : string, optional
        Name of the DC scaler, again it needs to be the same name as in
        Bluesky. Defaults to 'Lock DC'.
    ac_col : string, optional
        Name of the AC scaler. Defaults to 'Lock AC'.
    acoff_col : string, optional
        Name of the AC offset scaler. Defaults to 'Lock AC off'.
    kwargs :
        Passed to `load_scan`.

    Returns
    -------
    x : numpy.array
        Positioner data.
    dc : numpy.array
        DC response, normally corresponds to the x-ray absorption.
    ac-ac_off : numpy.array
        AC response minus the offset, normally corresponds to the x-ray
        magnetic dichroism.

    See also
    --------
    :func:`polartools.load_data.load_scan`
    """

    x, dc, ac, ac_off = load_scan(db, scan, positioner, [dc_col, ac_col,
                                  acoff_col], **kwargs)
    return x, dc, ac-ac_off


def load_dichro(scan, db=None, positioner='monochromator_energy',
                detector='Ion Ch 5', monitor='Ion Ch 4', **kwargs):

    """
    Load the x-ray magnetic dichroism from one scan taken in non-lock-in mode
    ('dichro').

    Parameters
    ----------
    db : database
        Databroker database. If None, it will attempt to read from csv files.
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky. Defauts to the x-ray energy.
    detector : string, optional
        Name of the DC scaler, again it needs to be the same name as in
        Bluesky. Defaults to ion chamber 5.
    monitor : string, optional
        Name of the AC scaler. Defaults to ion chamber 4.
    kwargs :
        Passed to `load_scan`.

    Returns
    -------
    x : numpy.array
        Positioner data.
    xanes : numpy.array
        X-ray absorption.
    xmcd : numpy.array
        X-ray magnetic dichroism.

    See also
    --------
    :func:`polartools.load_data.load_scan`
    """

    x0, y0 = load_absorption(db, scan, positioner=positioner, monitor=monitor,
                             detectors=[detector], **kwargs)
    size = x0.size//4

    x0 = x0.reshape(size, 4)
    y0 = y0.reshape(size, 4)

    x = x0.mean(axis=1)
    xanes = y0.mean(axis=1)
    xmcd = y0[:, [0, 3]].mean(axis=1) - y0[:, [1, 2]].mean(axis=1)

    return x, xanes, xmcd


def load_xanes(db, scans, return_mean=True, **kwargs):
    """
    Load multiple x-ray absorption energy scans.

    Parameters
    ----------
    db : database
        Databroker database. If None, it will attempt to read from csv files.
    scans : iterable
        Sequence of scan_ids our uids. If scan_id is passed, it will load the
        last scan with that scan_id. Use kwargs for search options.
    return_mean : boolean
        Flag to indicate if the averaging of multiple scans will be performed.
        Note that if True three outputs are generated, otherwise two.
    kwargs :
        Passed to `load_absorption`. Can be used to select the detector/monitor
        to be used. It is also passed to `load_bluesky` or `load_csv`, see
        the respective docstring for details.

    Returns
    -------
    energy : numpy.array
        X-ray energy.
    xanes : numpy.array
        X-ray absorption. If return_mean = True:
        xanes.shape = (len(scans), len(energy)), otherwise:
        xanes.shape = (len(energy))
    xanes_std : numpy.array, optional
        Error of the mean of x-ray absorption. This will only be returned if
        return_mean = True.

    See also
    --------
    :func:`polartools.load_data.load_absorption`
    :func:`polartools.load_data.load_bluesky`
    :func:`polartools.load_data.load_csv`
    """

    for scan in scans:
        if 'energy' not in locals():
            energy, xanes = load_absorption(db, scan, **kwargs)
        else:
            energy_tmp, xanes_tmp = load_absorption(db, scan, **kwargs)
            xanes = np.vstack((xanes, interp1d(energy_tmp, xanes_tmp,
                                               kind='linear',
                                               fill_value=np.nan,
                                               bounds_error=False)(energy)))

    if return_mean:
        if len(xanes.shape) == 2:
            xanes_std = xanes.std(axis=0)/np.sqrt(len(scans))
            xanes = xanes.mean(axis=0)
            return energy, xanes, xanes_std
        else:
            return energy, xanes, np.zeros((xanes.size))

    else:
        return energy, xanes


def load_xmcd(db, scans, return_mean=True, func=load_dichro, **kwargs):
    """
    Load multiple x-ray magnetic dichroism energy scans.

    Parameters
    ----------
    db : database
        Databroker database. If None, it will attempt to read from csv files.
    scans : iterable
        Sequence of scan_ids our uids. If scan_id is passed, it will load the
        last scan with that scan_id. Use kwargs for search options.
    return_mean : boolean
        Flag to indicate if the averaging of multiple scans will be performed.
        Note that if True five outputs are generated, otherwise three.
    func : function, optional
        Function that will load the XMCD signal. It has to have the call:
        func(db, scan, \*\*kwargs). Defaults to `load_dichro`, but see also
        `load_lockin`.
    kwargs :
        Passed to `load_absorption`. Can be used to select the detector/monitor
        to be used. It is also passed to `load_bluesky` or `load_csv`, see
        the respective docstring for details.

    Returns
    -------
    energy : numpy.array
        X-ray energy.
    xanes : numpy.array
        X-ray absorption. If return_mean = True:
        xanes.shape = (len(scans), len(energy)), otherwise:
        xanes.shape = (len(energy))
    xmcd : numpy.array
        X-ray magnetic dichroism. It has the same shape as the xanes.
    xanes_std : numpy.array, optional
        Error of the mean of x-ray absorption. This will only be returned if
        return_mean = True.
    xmcd_std : numpy.array, optional
        Error of the mean of x-ray magnetic dichroism. This will only be
        returned if return_mean = True.

    See also
    --------
    :func:`polartools.load_data.load_dichro`
    :func:`polartools.load_data.load_lockin`
    :func:`polartools.load_data.load_bluesky`
    :func:`polartools.load_data.load_csv`
    """

    for scan in scans:
        if 'energy' not in locals():
            energy, xanes, xmcd = func(db, scan, **kwargs)
        else:
            energy_tmp, xanes_tmp, xmcd_tmp = func(db, scan, **kwargs)
            xanes = np.vstack((xanes, interp1d(energy_tmp, xanes_tmp,
                                               kind='linear',
                                               fill_value=np.nan,
                                               bounds_error=False)(energy)))
            xmcd = np.vstack((xmcd, interp1d(energy_tmp, xmcd_tmp,
                                             kind='linear',
                                             fill_value=np.nan,
                                             bounds_error=False)(energy)))

    if return_mean:
        if len(xanes.shape) == 2:
            xanes_std = xanes.std(axis=0)/np.sqrt(len(scans))
            xanes = xanes.mean(axis=0)
            xmcd_std = xmcd.std(axis=0)/np.sqrt(len(scans))
            xmcd = xmcd.mean(axis=0)
            return energy, xanes, xmcd, xanes_std, xmcd_std
        else:
            return (energy, xanes, xmcd, np.zeros((xanes.size)),
                    np.zeros((xanes.size)))

    else:
        return energy, xanes, xmcd
