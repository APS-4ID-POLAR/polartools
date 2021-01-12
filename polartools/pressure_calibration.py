# Copyright (c) 2020, UChicago Argonne, LLC
# See LICENSE file.

import numpy as np
from scipy.interpolate import interp1d
from .load_data import load_scan, load_databroker
from .process_data import fit_bragg_peak


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


def xrd_calibrate_pressure(scan, *, db=None, bragg_peak=[1, 1, 1],
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
    >> kwargs = {'db': db, 'calibrant': 'Ag', 'temperature': 10}
    >> pressure = xrd_calibrate_pressure(100, \*\*kwargs)

    Parameters
    -----------
    scan : int or string
        Scan our uid. It will load the last scan with that scan_id. See kwargs
        for search options.
    db : databroker database, optional
        If None is passed it will try to load data from \*.csv files using
        `polartools.load_data.load_csv`.
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
        Passed to `polartools.load_data.load_table`.

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

    x, y = load_scan(db, scan, positioner, [detector], monitor=monitor,
                     **kwargs)

    fit = fit_bragg_peak(x, y, center=center, sigma=sigma, amplitude=amplitude,
                         fraction=fraction, fit_fraction=fit_fraction,
                         slope=slope, intercept=intercept, fit_slope=fit_slope)

    tth = fit.best_values['center']

    # TODO: we can probably add the option to search both 'baseline' and
    # 'primary' streams here.
    if isinstance(temperature, str):
        temperature = load_databroker(scan, db,
                                      stream='baseline')[temperature].mean()

    if isinstance(energy, str):
        energy = load_databroker(scan, db,
                                 stream='baseline')[energy].mean()

    return calculate_pressure(tth, temperature, energy, bragg_peak, calibrant,
                              tth_off=tth_offset)
