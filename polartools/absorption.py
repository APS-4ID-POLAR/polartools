"""
Functions to load and process x-ray absorption data.

.. autosummary::
   ~load_absorption
   ~load_dichro
   ~load_lockin
   ~load_multi_xas
   ~load_multi_dichro
   ~load_multi_lockin
"""

# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from .load_data import load_table, is_Bluesky_specfile
import numpy as np
from scipy.interpolate import interp1d
from spec2nexus.spec import (SpecDataFile, SpecDataFileNotFound,
                             NotASpecDataFile)

_spec_default_cols = dict(
    positioner='Energy',
    detector='IC5',
    monitor='IC4',
    dc_col='Lock DC',
    ac_col='Lock AC',
    acoff_col='Lock ACoff',
    )

_bluesky_default_cols = dict(
    positioner='monochromator_energy',
    detector='Ion Ch 5',
    monitor='Ion Ch 4',
    dc_col='Lock DC',
    ac_col='Lock AC',
    acoff_col='Lock AC off',
    )


def _select_default_names(source, **kwargs):
    # Select default parameters
    if isinstance(source, (str, SpecDataFile)):
        # It is csv or spec.
        try:
            # Checks spec origin.
            check = is_Bluesky_specfile(source, **kwargs)
            _defaults = _bluesky_default_cols if check else _spec_default_cols
        except (NotASpecDataFile, SpecDataFileNotFound):
            # If not a spec file, it must be csv, and use bluesky defaults.
            _defaults = _bluesky_default_cols
    else:
        # It is databroker.
        _defaults = _bluesky_default_cols
    return _defaults


def load_absorption(scan, source, positioner=None, detector=None, monitor=None,
                    transmission=True, **kwargs):
    """
    Load the x-ray absorption from one scan.

    Parameters
    ----------
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to the x-ray energy.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the ion chamber 5.
    monitor : string, optional
        Name of the monitor detector. If None is passed, it defaults to the ion
        chamber 4.
    transmission : bool, optional
        Flag to select between transmission mode -> ln(monitor/detector)
        or fluorescence mode -> detector/monitor
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

    Returns
    -------
    x : numpy.array
        Positioner data.
    y : numpy.array
        X-ray absorption.

    See also
    --------
    :func:`polartools.load_data.load_table`
    """

    _defaults = _select_default_names(source, **kwargs)

    if not positioner:
        positioner = _defaults['positioner']
    if not detector:
        detector = _defaults['detector']
    if not monitor:
        monitor = _defaults['monitor']

    # Load data
    table = load_table(scan, source, **kwargs)

    if transmission:
        return table[positioner], np.log(table[monitor]/table[detector])
    else:
        return table[positioner], table[detector]/table[monitor]


def load_lockin(scan, source, positioner=None, dc_col=None, ac_col=None,
                acoff_col=None, **kwargs):
    """
    Load the x-ray magnetic dichroism from one scan taken in lock-in mode.

    Parameters
    ----------
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defaults to the x-ray energy.
    dc_col : string, optional
        Name of the DC scaler, again it needs to be the same name as in
        Bluesky or SPEC. If None is passed, it defaults to 'Lock DC'.
    ac_col : string, optional
        Name of the AC scaler. If None is passed, it defaults to 'Lock AC'.
    acoff_col : string, optional
        Name of the AC offset scaler. If None is passed, it defaults to
        'Lock AC off'.
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

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
    :func:`polartools.load_data.load_table`
    """

    _defaults = _select_default_names(source, **kwargs)

    if not positioner:
        positioner = _defaults['positioner']
    if not dc_col:
        dc_col = _defaults['dc_col']
    if not ac_col:
        ac_col = _defaults['ac_col']
    if not acoff_col:
        acoff_col = _defaults['acoff_col']

    # Load data
    table = load_table(scan, source, **kwargs)

    return table[positioner], table[dc_col], table[ac_col] - table[acoff_col]


def load_dichro(scan, source, positioner=None, detector=None, monitor=None,
                transmission=True, **kwargs):

    """
    Load the x-ray magnetic dichroism from one scan taken in non-lock-in mode
    ('dichro').

    Parameters
    ----------
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id. See kwargs for search options.
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky. Defauts to the x-ray energy.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. Defaults to the ion chamber 5.
    monitor : string, optional
        Name of the monitor detector. Defaults to the ion chamber 4.
    transmission: bool, optional
        Flag to select between transmission mode -> ln(monitor/detector)
        or fluorescence mode -> detector/monitor.
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

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
    :func:`polartools.load_data.load_absorption`
    :func:`polartools.load_data.load_table`
    """

    # In SPEC the columns are different.
    if (isinstance(source, (str, SpecDataFile)) and source != 'csv' and not
            is_Bluesky_specfile(source, **kwargs)):

        if not positioner:
            positioner = _spec_default_cols['positioner']
        if not monitor:
            monitor = _spec_default_cols['monitor']
        if not detector:
            detector = _spec_default_cols['detector']

        table = load_table(scan, source, **kwargs)
        x = table[positioner]
        monp = table[monitor+'(+)']
        detp = table[detector+'(+)']
        monm = table[monitor+'(-)']
        detm = table[detector+'(-)']

        if transmission:
            xanes = (np.log(monp/detp) + np.log(monm/detm))/2.
            xmcd = np.log(monp/detp) - np.log(monm/detm)
        else:
            xanes = (detp/monp + detm/monm)/2.
            xmcd = detp/monp - detm/monm

    else:
        x0, y0 = load_absorption(scan, source, positioner, detector, monitor,
                                 transmission, **kwargs)
        size = x0.size//4

        x0 = x0.reshape(size, 4)
        y0 = y0.reshape(size, 4)

        x = x0.mean(axis=1)
        xanes = y0.mean(axis=1)
        xmcd = y0[:, [0, 3]].mean(axis=1) - y0[:, [1, 2]].mean(axis=1)

    return x, xanes, xmcd


def load_multi_xas(scans, source, return_mean=True, positioner=None,
                   detector=None, monitor=None, transmission=True, **kwargs):
    """
    Load multiple x-ray absorption energy scans.

    Parameters
    ----------
    scans : iterable
        Sequence of scan_ids our uids. If scan_id is passed, it will load the
        last scan with that scan_id. Use kwargs for search options.
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    return_mean : boolean, optional
        Flag to indicate if the averaging of multiple scans will be performed.
        Note that if True three outputs are generated, otherwise two.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to the x-ray energy.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the ion chamber 5.
    monitor : string, optional
        Name of the monitor detector. If None is passed, it defaults to the ion
        chamber 4.
    transmission: bool, optional
        Flag to select between transmission mode -> ln(monitor/detector)
        or fluorescence mode -> detector/monitor
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

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
    """

    for scan in scans:
        if 'energy' not in locals():
            energy, xanes = load_absorption(
                scan, source, positioner, detector, monitor, transmission,
                **kwargs
                )
        else:
            energy_tmp, xanes_tmp = load_absorption(
                scan, source, positioner, detector, monitor, transmission,
                **kwargs
                )
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


def load_multi_dichro(scans, source, return_mean=True, positioner=None,
                      detector=None, monitor=None, transmission=True,
                      **kwargs):
    """
    Load multiple x-ray magnetic dichroism energy "dichro" scans.

    Parameters
    ----------
    scans : iterable
        Sequence of scan_ids our uids. If scan_id is passed, it will load the
        last scan with that scan_id. Use kwargs for search options.
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    return_mean : boolean, optional
        Flag to indicate if the averaging of multiple scans will be performed.
        Note that if True three outputs are generated, otherwise two.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to the x-ray energy.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the ion chamber 5.
    monitor : string, optional
        Name of the monitor detector. If None is passed, it defaults to the ion
        chamber 4.
    transmission: bool, optional
        Flag to select between transmission mode -> ln(monitor/detector)
        or fluorescence mode -> detector/monitor
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

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
    """

    for scan in scans:
        if 'energy' not in locals():
            energy, xanes, xmcd = load_dichro(
                scan, source, positioner, detector, monitor, transmission,
                **kwargs
                )
        else:
            energy_tmp, xanes_tmp, xmcd_tmp = load_dichro(
                scan, source, positioner, detector, monitor, transmission,
                **kwargs
                )
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


def load_multi_lockin(scans, source, return_mean=True, positioner=None,
                      dc_col=None, ac_col=None, acoff_col=None, **kwargs):
    """
    Load multiple x-ray magnetic dichroism energy "lockin" scans.

    Parameters
    ----------
    scans : iterable
        Sequence of scan_ids our uids. If scan_id is passed, it will load the
        last scan with that scan_id. Use kwargs for search options.
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    return_mean : boolean, optional
        Flag to indicate if the averaging of multiple scans will be performed.
        Note that if True three outputs are generated, otherwise two.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to the x-ray energy.
    dc_col : string, optional
        Name of the DC scaler, again it needs to be the same name as in
        Bluesky or SPEC. If None is passed, it defaults to 'Lock DC'.
    ac_col : string, optional
        Name of the AC scaler. If None is passed, it defaults to 'Lock AC'.
    acoff_col : string, optional
        Name of the AC offset scaler. If None is passed, it defaults to
        'Lock AC off'.
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

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
    :func:`polartools.load_data.load_lockin`
    """

    for scan in scans:
        if 'energy' not in locals():
            energy, xanes, xmcd = load_lockin(
                scan, source, positioner, dc_col, ac_col, acoff_col, **kwargs
                )
        else:
            energy_tmp, xanes_tmp, xmcd_tmp = load_lockin(
                scan, source, positioner, dc_col, ac_col, acoff_col, **kwargs
                )
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


def normalize_absorption(energy, xanes, pre_edge_range, pos_edge_range, e0,
                         pre_edge_order=1, pos_edge_order=1):
    """
    Extract pre- and post-edge normalization curves by fitting polynomials.

    Parameters
    ----------
    energy : list
        Incident energy.
    xanes : list
        X-ray absorption.
    pre_edge_range : list
        List with the energy ranges [initial, final] of the pre-edge region
        **relative** to the absorption edge.
    pos_edge_range : list
        List with the energy ranges [initial, final] of the post-edge region
        **relative** to the absorption edge.
    e0 : float
        Absorption edge energy.
    pre_edge_order : int, optional
        Order of the polynomial to be used in the pre-edge. Defauts to 1.
    pos_edge_order : int, optional
        Order of the polynomial to be used in the post-edge. Defauts to 1.

    Returns
    -------
    pre_edge : numpy.array
        Pre-edge polynomial.
    pos_edge : numpy.array
        Post-edge polynomial.
    jump : float
        Size of the absorption jump.

    See also
    --------
    :func:`numpy.polyfit`
    """

    energy = np.array(energy)
    xanes = np.array(xanes)

    # Process pre-edge
    index = (energy > pre_edge_range[0]) & (energy < pre_edge_range[1])
    pre_edge = np.poly1d(np.polyfit(energy[index], xanes[index],
                                    pre_edge_order))(energy)
    processed_xanes = xanes - pre_edge

    # Process pos-edge
    index = (energy > pos_edge_range[0]) & (energy < pos_edge_range[1])
    pos_edge_func = np.poly1d(np.polyfit(energy[index], processed_xanes[index],
                              pos_edge_order))
    pos_edge = pos_edge_func(energy)
    jump = pos_edge_func(e0)

    return pre_edge, pre_edge + pos_edge, jump
