'''
 Copyright (c) 2020, UChicago Argonne, LLC
 See LICENSE file.
'''

import numpy as np
from scipy.interpolate import interp1d
from lmfit.models import LinearModel, PseudoVoigtModel


def _generate_initial_guess(params, x, y, center, sigma, amplitude, fraction,
                            fit_fraction, m, b, fit_m):
    """
    Adds initial guess to the lmfit parameters.

    Uses lmfit (https://lmfit.github.io/lmfit-py/index.html).

    WARNING: The following initial guess (if None is passed) and constrains are
    applied:
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
    m, b : float, optional
        Initial guess parameters of linear function. For more details,
        see: :func: `lmfit.models.LinearModel`
    fit_n : boolean, optional
        Flag to control if m will be varied.

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
        center = x[index]
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

    if not m:
        m = (y[-1]-y[0])/(x[-1]-x[0])
    params['m'].set(m, vary=fit_m)

    if not b:
        b = x[0]
    params['b'].set(b)

    return params


def fit_bragg_peak(x, y, center=None, sigma=None, amplitude=None, fraction=None,
                   fit_fraction=True, m=None, b=None, fit_m=True):
    """
    Fit Bragg peak with a pseudo-voigt function.

    Uses lmfit (https://lmfit.github.io/lmfit-py/index.html).

    WARNING: This imposes constrains in the fit that are described in
    :func: `polartools.pressure_calibration._generate_initial_guess`

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
    m, b : float, optional
        Initial guess parameters of linear function. For more details,
        see: :func: `lmfit.models.LinearModel`
    fit_n : boolean, optional
        Flag to control if m will be varied.

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
                                     fraction, fit_fraction,  m, b, fit_m)

    return model.fit(y, params=params, x=x)


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


def calculate_pressure(tth, temperature, energy, bragg_peak, calibrant,
                       tth_off=0.0):
    """
    Calculate the pressure using diffraction from Au or Ag.

    See Holzapfel et al., J. Phys. Chem. Ref. Data 30, 515 (2001) for more
    details.

    Parameters
    -----------
    tth : float
        Two theta of the selected Bragg peak.
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
    pressure : float
        Calculated pressure in GPa.
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

    # Calculate atomic volume
    lamb = h*c/energy/1000.
    d = lamb/2/np.sin((tth-tth_off)/2.*np.pi/180.)
    a = d*np.sqrt(bragg_peak[0]**2+bragg_peak[1]**2+bragg_peak[2]**2)
    v = a**3/4.

    # Calculate pressure
    x = (v/v0)**0.3333
    pfg0 = afg*(z/v0)**1.6666
    c0 = -1*np.log(3*k0/pfg0)
    c2 = (3/2)*(kp0-3)-c0
    pressure = 3*k0*(1-x)/x**5*np.exp(c0*(1-x))*(1+c2*x*(1-x))

    return pressure
