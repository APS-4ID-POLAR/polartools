"""
Functions to load and process x-ray absorption data.

.. autosummary::
   ~load_absorption
   ~load_dichro
   ~load_lockin
   ~load_multi_xas
   ~load_multi_dichro
   ~load_multi_lockin
   ~normalize_absorption
"""

# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from .load_data import load_table, is_Bluesky_specfile
import numpy as np
from scipy.interpolate import interp1d
from spec2nexus.spec import (SpecDataFile, SpecDataFileNotFound,
                             NotASpecDataFile)
from larch.xafs.pre_edge import _finde0
from larch.math import index_nearest
from lmfit.models import PolynomialModel
from warnings import warn

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


def normalize_absorption(energy, mu, *, e0=None, edge_step=None,
                         pre_range=None, pre_order=1, nvict=0, pre_pars=None,
                         post_range=None, post_order=None, post_pars=None,
                         flat_range=None, flat_order=None, flat_pars=None):
    """
    Extract pre- and post-edge normalization curves by fitting polynomials.

    This is a wrapper of `larch.xafs.preedge` that adds the flattening of the
    post-edge.

    Please see the larch documentation_ for details on how the optional
    parameters are determined.

    .. _documentation:\
    https://xraypy.github.io/xraylarch/xafs/preedge.html#the-pre-edge-function

    Parameters
    ----------
    energy : iterable
        Incident energy. Must be in eV. It will raise an error if it notices
        the maximum of this list is below 100.
    xanes : iterable
        X-ray absorption.
    e0 : float or int, optional
        Absorption edge energy.
    pre_range : list, optional
        List with the energy ranges [initial, final] of the pre-edge region
        **relative** to the absorption edge.
    post_range : list, optional
        List with the energy ranges [initial, final] of the post-edge region
        **relative** to the absorption edge.
    pre_exponent : int, optional
        Energy exponent to use. The pre-edge background is modelled by a line
        that is fit to xanes(energy)*energy**pre_exponent. Defaults to 0.
    post_order : int, optional
        Order of the polynomial to be used in the post-edge. If None, it will
        be determined by `larch.xafs.preedge`:
        - nnorm = 2 if post_range[1]-post_range[0]>350, 1 if
        50 < post_range[1]-post_range[0] < 350, or 0 otherwise.
    fpars : lmfit.Parameters, optional
        Option to input the initial parameters to the polynomial used in the
        data flattening. These will be labelled 'c0', 'c1', ..., depending on
        `post_order`. See lmfit.models.PolynomialModel_ for details.

        .. _lmfit.models.PolynomialModel:\
        https://lmfit.github.io/lmfit-py/builtin_models.html#polynomialmodel

    Returns
    -------
    results : dict
        Dictionary with the results and parameters of the normalization and
        flattening. The most important items are:

        - 'energy' -> incident energy.
        - 'raw' -> raw xanes.
        - 'norm' -> normalized xanes.
        - 'flat' -> flattened xanes.

    See also
    --------
    :func:`larch.xafs.preedge`
    """

    # Start output dictionary
    results = {}
    sort = np.argsort(energy)
    results['energy'] = np.array(energy)[sort]
    results['mu'] = np.array(mu)[sort]

    # Process pre-edge
    pre1, pre2 = (None, None) if not pre_range else tuple(pre_range)

    pre_results = pre_edge_background(
        results['energy'], results['mu'], e0=e0, pre1=pre1, pre2=pre2,
        pre_order=pre_order, nvict=nvict, pre_pars=pre_pars
        )

    results.update(pre_results)

    # Process post-edge
    post1, post2 = (None, None) if not post_range else tuple(post_range)

    post_results = post_edge_background(
        results['energy'], results['mu'], preedge=results['preedge'],
        e0=results['e0'], edge_step=edge_step, post1=post1, post2=post2,
        post_order=post_order, post_pars=post_pars
        )

    results.update(post_results)

    # Flatten post-edge
    use_post = not bool(flat_order or flat_range or flat_pars)
    postedge_results = post_results if use_post else None

    flat1, flat2 = (None, None) if not flat_range else tuple(flat_range)

    flat_results = post_edge_flatten(
        results['energy'], results['norm'], e0=results['e0'],
        flat1=flat1, flat2=flat2, flat_order=flat_order, flat_pars=flat_pars,
        bkg_results=postedge_results
        )

    results.update(flat_results)

    return results


def pre_edge_background(energy, mu, e0=None, pre1=None, pre2=None, pre_order=1,
                        nvict=0, pre_pars=None):

    # Edge energy
    e0, unit = _process_e0(energy, mu, e0)

    # Process pre-edge range
    pre1, pre2 = _process_preedge_params(energy, e0, pre1, pre2)

    # Find pre-edge background
    index = np.logical_and(energy > pre1+e0, energy < pre2+e0)

    mu_fit = mu[index]*energy[index]**nvict
    energy_fit = energy[index]

    fit = _fit_polynomial(energy_fit, mu_fit, pre_order, pars=pre_pars)

    precoefs = fit.best_values
    preedge = fit.eval(x=energy)*energy**(-nvict)

    return dict(preedge=preedge, pre1=pre1, pre2=pre2, pre_order=pre_order,
                nvict=nvict, e0=e0, precoefs=precoefs, energy_unit=unit)


def post_edge_background(energy, mu, preedge=None, e0=None, edge_step=None,
                         post1=None, post2=None, post_order=None,
                         post_pars=None):

    # Edge energy
    e0, unit = _process_e0(energy, mu, e0)

    # Post-edge params
    post1, post2, post_order = _process_postedge_params(
        energy, e0, post1, post2, post_order
        )

    # Pre-edge
    if preedge is None:
        preedge = np.zeros(len(energy))
    else:
        if len(preedge) != len(energy):
            raise TypeError('preedge must have the same dimensions as energy.')

    # Normalization
    index = np.logical_and(energy > post1+e0, energy < post2+e0)
    energy_fit = energy[index]
    mu_fit = (mu-preedge)[index]

    fit = _fit_polynomial(energy_fit, mu_fit, post_order, pars=post_pars)

    postcoefs = fit.best_values
    postedge = preedge + fit.eval(x=energy)

    if not edge_step:
        ie0 = index_nearest(energy, e0)
        edge_step = postedge[ie0] - preedge[ie0]

    norm = (mu - preedge)/edge_step

    return dict(energy=energy, mu=mu, postedge=postedge, preedge=preedge,
                postcoefs=postcoefs, post1=post1, post2=post2,
                post_order=post_order, energy_unit=unit, norm=norm,
                edge_step=edge_step)


def post_edge_flatten(energy, norm, e0=None, flat1=None, flat2=None,
                      flat_order=None, flat_pars=None, bkg_results=None):
    """
    Flattens the normalized absorption.

    This is a slightly modified version from that in the larch_ package.

    .. _larch: https://xraypy.github.io/xraylarch/index.html

    The energies inputs must have the same units.

    Parameters
    ----------
    energy : iterable
        Incident energy in eV.
    norm : iterable
        Normalized x-ray absorption.
    e0 : float or int
        Absorption edge energy.
    norm1 : float or int
        Low energy limit of normalization range with respect to e0.
    norm2 : float or int
        High energy limit of normalization range with respect to e0.
    nnorm : int
        Degree of polynomial to be used.
    fpars : lmfit.Parameters, optional
        Option to input the initial parameters. These will be labelled 'c0',
        'c1', ..., depending on `nnorm`. See lmfit.models.PolynomialModel_ for
        details.

        .. _lmfit.models.PolynomialModel:\
        https://lmfit.github.io/lmfit-py/builtin_models.html#polynomialmodel

    Returns
    -------
    flat : numpy.array
        Flattened x-ray absorption.

    See also
    --------
    :func:`larch.xafs.pre_edge`
    :func:`lmfit.models.PolynomialModel`
    """

    # Edge energy
    e0, unit = _process_e0(energy, norm, e0)
    ie0 = index_nearest(energy, e0)

    # Runs the post edge fit with flat parameters if needed.
    if not bkg_results:
        bkg_results = post_edge_background(
            energy, norm, e0=e0, post1=flat1, post2=flat2,
            post_order=flat_order, post_pars=flat_pars
            )

        flat_function = bkg_results['postedge']
    else:
        flat_function = (bkg_results['postedge'] - bkg_results['preedge'])
        flat_function /= bkg_results['edge_step']

    flat1 = bkg_results['post1']
    flat2 = bkg_results['post2']
    flat_order = bkg_results['post_order']
    flatcoefs = bkg_results['postcoefs']

    flat = norm - (flat_function - flat_function[ie0])
    flat[:ie0] = norm[:ie0]

    return dict(energy=energy, norm=norm, flat=flat, flat1=flat1, flat2=flat2,
                flat_order=flat_order, flatcoefs=flatcoefs,
                flat_function=flat_function, energy_unit=unit)


def _process_e0(energy, mu, e0):

    # Energy unit
    unit = 'keV' if np.nanmax(energy) < 100 else 'eV'

    # Edge energy
    if e0 is None:
        e0 = _finde0(energy, mu)
    elif e0 < energy.min() or e0 > energy.max():
        e0 = _finde0(energy, mu)
        warn("e0 is not contained in the provided energy range, "
             "using e0 = {:0.2f} {}".format(e0, unit))
    else:
        e0 = energy[index_nearest(energy, e0)]

    return e0, unit


def _process_preedge_params(energy, e0, pre1, pre2):

    # Process pre-edge range
    if pre1 is None:
        pre1 = np.nanmin(energy[1:-1])-e0  # avoid first/last points
    if pre2 is None:
        pre2 = 5.0*round(pre1/15.0)  # from larch, don't understand logic.
    if pre1 > pre2:
        pre1, pre2 = pre2, pre1

    return pre1, pre2


def _process_postedge_params(energy, e0, post1, post2, post_order):

    # Post-edge params
    if post2 is None:
        post2 = np.nanmax(energy[1:-1])-e0  # avoid first/last two points
    if post1 is None:
        post1 = min(150, 5.0*round(post2/15.0))
    if post1 > post2:
        post1, post2 = post2, post1
    if post_order is None:
        if post2-post1 < 50:
            post_order = 0
        elif post2-post1 < 350:
            post_order = 1
        else:
            post_order = 2

    return post1, post2, post_order


def _fit_polynomial(x, y, order, pars=None):

    model = PolynomialModel(order)
    if not pars:
        pars = model.guess(y, x=x)
    return model.fit(y, pars, x=x)
