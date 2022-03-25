"""
Functions to do pressure calibration using Au or Ag diffraction.

.. autosummary::
   ~load_ag_params
   ~load_au_params
   ~calculate_tth
   ~calculate_pressure
   ~xrd_calibrate_pressure
"""

# Copyright (c) 2020, UChicago Argonne, LLC
# See LICENSE file.

from numpy import (sin, exp, log, pi, sqrt, array, linspace, nanmax, nanmin,
                   copy)
from scipy.interpolate import interp1d
from .load_data import load_table
from lmfit.models import LinearModel, PseudoVoigtModel


def load_ag_params(temperature):
    """
    Load the Ag parameters for calculating the pressure.

    These parameters were extracted from Holzapfel et al., J. Phys. Chem. Ref.
    Data 30, 515 (2001). It linearly interpolates the data.

    Parameters
    -----------
    temperature : float
        Measurement temperature in Kelvin.

    Returns
    -----------
    v0_out, k0_out, kp0_out : float
        Volume, K and K' calibration parameters.
    """
    if temperature > 500:
        raise ValueError(f'temperature must be < 500 K, but {temperature} \
        was entered.')

    # Data from paper
    v0 = [16.8439, 16.8439, 16.8512, 16.8815, 16.9210, 16.9644, 17.0099,
          17.057, 17.1055, 17.1553, 17.2063, 17.2585]
    k0 = [110.85, 110.85, 110.31, 108.68, 106.83, 104.91, 102.96, 101.0, 99.03,
          97.05, 95.07, 93.07]
    kp0 = [6, 6, 6.01, 6.03, 6.05, 6.08, 6.12, 6.15, 6.19, 6.22, 6.26, 6.3]
    temp0 = [0, 10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

    v0_out = float(interp1d(temp0, v0, kind='linear')(temperature))
    k0_out = float(interp1d(temp0, k0, kind='linear')(temperature))
    kp0_out = float(interp1d(temp0, kp0, kind='linear')(temperature))
    return v0_out, k0_out, kp0_out


def load_au_params(temperature):
    """
    Load the Au parameters for calculating the pressure.

    These parameters were extracted from Holzapfel et al., J. Phys. Chem. Ref.
    Data 30, 515 (2001). It linearly interpolates the data.

    Parameters
    -----------
    temperature : float
        Measurement temperature in Kelvin.

    Returns
    -----------
    v0_out, k0_out, kp0_out : float
        Volume, K and K' calibration parameters.
    """

    if temperature > 500:
        raise ValueError(f'temperature must be < 500 K, but {temperature} \
        was entered.')

    # Data from paper
    v0 = [16.7905, 16.7906, 16.7984, 16.8238, 16.8550, 16.8885, 16.9232,
          16.959, 16.9956, 17.0329, 17.071, 17.1098]
    k0 = [180.93, 180.93, 179.94, 177.51, 174.86, 172.16, 169.43, 166.7,
          163.96, 161.21, 158.46, 155.7]
    kp0 = [6.08, 6.08, 6.09, 6.11, 6.13, 6.15, 6.17, 6.20, 6.23, 6.25, 6.28,
           6.31]
    temp0 = [0, 10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

    v0_out = float(interp1d(temp0, v0, kind='linear')(temperature))
    k0_out = float(interp1d(temp0, k0, kind='linear')(temperature))
    kp0_out = float(interp1d(temp0, kp0, kind='linear')(temperature))
    return v0_out, k0_out, kp0_out


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
        index = y == nanmax(y)
        center = x[index][0]
    params['center'].set(center, min=nanmin(x), max=nanmax(x))

    if not sigma:
        sigma = (nanmax(x)-nanmin(x))/10
    params['sigma'].set(sigma, min=0)

    if not amplitude:
        amplitude = nanmax(y)*sigma*sqrt(pi/log(2))
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
    `polartools.pressure_calibration._generate_initial_guess`

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

    x = copy(x)
    y = copy(y)

    model = PseudoVoigtModel() + LinearModel()
    params = model.make_params()

    params = _generate_initial_guess(params, x, y, center, sigma, amplitude,
                                     fraction, fit_fraction,  slope, intercept,
                                     fit_slope)

    return model.fit(y, params=params, x=x)


def calculate_tth(pressure, temperature, energy, bragg_peak, calibrant,
                  tth_off=0.0):
    """
    Calculate the two theta of Au or Ag Bragg peak of given a pressure.

    See Holzapfel et al., J. Phys. Chem. Ref. Data 30, 515 (2001) for more
    details.

    Parameters
    -----------
    pressure : float or iterable
        Pressure (or list of) to be converted.
    temperature : float
        Measurement temperature in Kelvin.
    energy : float
        X-ray energy used in keV.
    bragg_peak : iterable
        List containing the Bragg peak indices [H, K, L].
    calibrant : string
        Selects the calibrant used. Options are 'Au' or 'Ag'.
    tth_off : float, optional
        Offset between the reference two theta and the measured value.

    Returns
    -----------
    tth : float or numpy.array
        Calculated tth in degrees. A numpy.array is returned if an iterable is
        passed to pressure.
    """

    tth_ref = linspace(1e-2, 180, 10000)
    pressure_ref = calculate_pressure(tth_ref, temperature, energy, bragg_peak,
                                      calibrant, tth_off=tth_off)

    tth = interp1d(pressure_ref, tth_ref)(pressure)

    return float(tth) if tth.size == 1 else tth


def calculate_pressure(tth, temperature, energy, bragg_peak, calibrant,
                       tth_off=0.0):
    """
    Calculate the pressure using diffraction from Au or Ag.

    See Holzapfel et al., J. Phys. Chem. Ref. Data 30, 515 (2001) for more
    details.

    Parameters
    -----------
    tth : float or iterable
        Two theta of the selected Bragg peak. It can be a list of two theta.
    temperature : float
        Measurement temperature in Kelvin.
    energy : float
        X-ray energy used in keV.
    bragg_peak : iterable
        List containing the Bragg peak indices [H, K, L].
    calibrant : string
        Selects the calibrant used. Options are 'Au' or 'Ag'.
    tth_off : float, optional
        Offset between the reference two theta and the measured value.

    Returns
    -----------
    pressure : float or numpy.array
        Calculated pressure in GPa. A numpy.array is returned if an iterable is
        passed to tth.
    """
    # Constants
    afg = 2337  # GPa.AA^5
    h = 4.135667662E-15  # eV.s
    c = 299792458E10  # AA/s

    # Loading parameters
    if calibrant == 'Au':
        z = 79
        v0, k0, kp0 = load_au_params(temperature)
    elif calibrant == 'Ag':
        z = 47
        v0, k0, kp0 = load_ag_params(temperature)
    else:
        raise ValueError(f'calibrant must be "Au" or "Ag", but {calibrant} was\
        entered')

    # If it's not a number it will turn tth into a numpy.array
    if not isinstance(tth, (int, float)):
        tth = array(tth)

    # Calculate atomic volume
    lamb = h*c/energy/1000.
    d = lamb/2/sin((tth-tth_off)/2.*pi/180.)
    a = d*sqrt(bragg_peak[0]**2+bragg_peak[1]**2+bragg_peak[2]**2)
    v = a**3/4.

    # Calculate pressure
    x = (v/v0)**0.3333
    pfg0 = afg*(z/v0)**1.6666
    c0 = -1*log(3*k0/pfg0)
    c2 = (3/2)*(kp0-3)-c0
    pressure = 3*k0*(1-x)/x**5*exp(c0*(1-x))*(1+c2*x*(1-x))

    return pressure


def xrd_calibrate_pressure(scan, source, bragg_peak=(1, 1, 1),
                           calibrant='Au', temperature=300,
                           energy='monochromator_energy',
                           positioner='huber_tth', detector='APDSector4',
                           monitor=None, center=None, sigma=None,
                           amplitude=None, fraction=None, fit_fraction=True,
                           slope=None, intercept=None, fit_slope=True,
                           tth_offset=0.0, **kwargs):
    """
    Calibrate pressure using x-ray diffraction.

    This is a wrapper of multiple functions in `polartools` meant to make it
    more convenient to extract the pressure in XRD measurements.

    For a given experiment, it is recommended to create a kwarg dictionary that
    holds the input parameters that are different from the standard defaults,
    and make it easier to process multiple scans. For instance:

    .. code-block:: python

        kwargs = {'db': db, 'calibrant': 'Ag', 'temperature': 10}
        pressure = xrd_calibrate_pressure(100, **kwargs)

    Parameters
    -----------
    scan : int or string
        Scan our uid. It will load the last scan with that scan_id. See kwargs
        for search options.
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    bragg_peak : iterable, optional
        List containing the Bragg peak indices [H, K, L].
    calibrant : string, optional
        Selects the calibrant used. Options are 'Au' or 'Ag'.
    temperature : float or string, optional
        A string can only be passed if using databroker. In this case it will
        read the temperature from the database baseline stream. If float is
        passed, then it is the temperature in Kelvin.
    energy : float or string, optional
        A string can only be passed if using databroker. In this case it will
        read the energy from the database baseline stream. If float is passed,
        then it is the energy in keV.
    positioner : string, optional
        2 theta motor name.
    detector : string, optional
        XRD detector name.
    monitor : string, optional
        Monitor detector name.
    center, sigma, amplitude, fraction : float, optional
        Initial guess parameters of pseudo-voigt function. For more details,
        see: `lmfit.models.PseudoVoigtModel`.
    fit_fraction : boolean, optional
        Flag to control if fraction will be varied.
    slope, intercept : float, optional
        Initial guess parameters of linear function. For more details,
        see: `lmfit.models.LinearModel`.
    fit_slope : boolean, optional
        Flag to control if slope will be varied.
    tth_off : float, optional
        Offset between the reference two theta and the measured value.
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query, use_db_v1.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

    Returns
    -----------
    pressure : float
        Calculated pressure in GPa.

    See also
    --------
    :func:`polartools.load_data.load_table`
    :func:`polartools.load_data.load_scan`
    :func:`polartools.load_data.load_csv`
    :func:`polartools.load_data.load_databroker`
    :func:`polartools.process_data.fit_bragg_peak`
    :func:`polartools.pressure_calibration.calculate_pressure`
    :func:`lmfit.models.LinearModel`
    :func:`lmfit.models.PseudoVoigtModel`
    """

    # TODO: Is there a better way to define some of these keyword arguments?
    # for example, positioner and detectors might be retrievable from metadata.

    table = load_table(scan, source, **kwargs)
    x = table[positioner]
    y = table[detector]
    if monitor:
        y /= table[monitor]

    fit = fit_bragg_peak(x, y, center=center, sigma=sigma, amplitude=amplitude,
                         fraction=fraction, fit_fraction=fit_fraction,
                         slope=slope, intercept=intercept, fit_slope=fit_slope)

    tth = fit.best_values['center']

    # TODO: we can probably add the option to search both 'baseline' and
    # 'primary' streams here.
    if isinstance(temperature, str):
        temperature = load_table(
            scan, source, stream='baseline', name_format='scan_{}_primary.csv'
            )[temperature].mean()

    if isinstance(energy, str):
        energy = load_table(
            scan, source, stream='baseline', name_format='scan_{}_primary.csv'
            )[energy].mean()

    return calculate_pressure(tth, temperature, energy, bragg_peak, calibrant,
                              tth_off=tth_offset)
