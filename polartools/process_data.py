# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

import numpy as np
from lmfit.models import LinearModel, PseudoVoigtModel


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


def _generate_initial_guess(params, x, y, center, sigma, amplitude, fraction,
                            fit_fraction, slope, intercept, fit_slope):
    """
    Adds initial guess to the lmfit parameters.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    WARNING: The following initial guesses (if None is passed) and constrains
    are applied:
    - center
        value: x position where y is maximum.
        constrain: min = min(x), max = max(x)
    - sigma
        value: (max(x) - min(x))/10
        constrain: min = 0
    - amplitude
        value: max(y)*sigma*sqrt(pi/2)
        constrain: min = 0
    - fraction
        value: 0.5
        constrain: min = 0, max = 1
    - m
        value: (y[-1]-y[0])/(x[-1]-x[0])
        No constrain
    - b
        value: x[0]
        No constrain

    Parameters
    ----------
    params: lmfit Parameters class
        Will hold the initial parameters initial setup.
    x : iterable
        List of x-axis values.
    y : iterable
        List of y-axis values.
    center, sigma, amplitude, fraction : float, optional
        Initial guess parameters of pseudo-voigt function. For more details,
        see: :func: `lmfit.models.PseudoVoigtModel`
    fit_fraction : boolean, optional
        Flag to control if fraction will be varied.
    slope, intercept : float, optional
        Initial guess parameters of linear function. For more details,
        see: :func: `lmfit.models.LinearModel`
    fit_slope : boolean, optional
        Flag to control if slope will be varied.

    Returns
    -------
    params : lmfit Parameters class
        Holds the initial parameters initial setup.

    See also
    --------
    :func:`lmfit.models.PseudoVoigtModel`
    :func:`lmfit.models.LinearModel`
    """

    if not center:
        index = y == np.nanmax(y)
        center = x[index][0]
    params['center'].set(center, min=np.nanmin(x), max=np.nanmax(x))

    if not sigma:
        sigma = (np.nanmax(x)-np.nanmin(x))/10
    params['sigma'].set(sigma, min=0)

    if not amplitude:
        amplitude = np.nanmax(y)*sigma*np.sqrt(np.pi/np.log(2))
    params['amplitude'].set(amplitude, min=0)

    if not fraction:
        fraction = 0.5
    if (fraction < 0) or (fraction > 1):
        print(f'WARNING: fraction must be 0 <= fraction <= 1, but you provided\
         fraction = {fraction}. Using fraction = 0.5.')
        fraction = 0.5
    params['fraction'].set(fraction, min=0, max=1, vary=fit_fraction)

    if not slope:
        slope = (y[-1]-y[0])/(x[-1]-x[0])
    params['slope'].set(slope, vary=fit_slope)

    if not intercept:
        intercept = x[0]
    params['intercept'].set(intercept)

    return params


def fit_bragg_peak(x, y, center=None, sigma=None, amplitude=None, fraction=None,
                   fit_fraction=True, slope=None, intercept=None,
                   fit_slope=True):
    """
    Fit Bragg peak with a pseudo-voigt function.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    WARNING: This imposes constrains to the fit that are described in
    :func: `polartools.process_data._generate_initial_guess`

    Parameters
    ----------
    x : iterable
        List of x-axis values.
    y : iterable
        List of y-axis values.
    center, sigma, amplitude, fraction : float, optional
        Initial guess parameters of pseudo-voigt function. For more details,
        see: :func: `lmfit.models.PseudoVoigtModel`
    fit_fraction : boolean, optional
        Flag to control if fraction will be varied.
    slope, intercept : float, optional
        Initial guess parameters of linear function. For more details,
        see: :func: `lmfit.models.LinearModel`
    fit_slope : boolean, optional
        Flag to control if slope will be varied.

    Returns
    -------
    fit : lmfit ModelResult class
        Contains the fit results. See:
        https://lmfit.github.io/lmfit-py/model.html#the-modelresult-class

    See also
    --------
    :func:`lmfit.models.PseudoVoigtModel`
    :func:`lmfit.models.LinearModel`
    """

    x = np.copy(x)
    y = np.copy(y)

    model = PseudoVoigtModel() + LinearModel()
    params = model.make_params()

    params = _generate_initial_guess(params, x, y, center, sigma, amplitude,
                                     fraction, fit_fraction,  slope, intercept,
                                     fit_slope)

    return model.fit(y, params=params, x=x)
