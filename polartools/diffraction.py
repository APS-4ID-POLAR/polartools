# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

import numpy as np
from pandas import DataFrame
from lmfit.models import (
    GaussianModel,
    LorentzianModel,
    VoigtModel,
    LinearModel,
    PseudoVoigtModel,
)
import matplotlib.pyplot as plt
from .load_data import load_table


def fit_peak(xdata, ydata, model="Gaussian", output=False):
    """
    Fit Bragg peak with model of choice: Gaussian, Lorentzian, Voigt, PseudoVoigt.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    WARNING: This imposes constrains to the fit that are described in
    :func: `polartools.process_data._generate_initial_guess`

    Parameters
    ----------
    xdata : iterable
        List of x-axis values.
    yydata : iterable
        List of y-axis values.
    model:
        fit model: Gaussian, Lorentzian, Voigt, PseudoVoigt
    output:
        Output fit parameters and plot data+fit for each scan.

    Returns
    -------
    fit : lmfit ModelResult class
        Contains the fit results. See:
        https://lmfit.github.io/lmfit-py/model.html#the-modelresult-class
    """

    models = {
        "Gaussian": GaussianModel(),
        "Lorentzian": LorentzianModel(),
        "Voigt": VoigtModel(),
        "PseudoVoigt": PseudoVoigtModel(),
    }
    try:
        peak_mod = models[model]
    except KeyError:
        raise ValueError(
            f"model = {model} is invalid. Only these values are accepted: "
            f"{list(models.keys())}"
        )
    background = LinearModel()
    mod = peak_mod + background
    pars = background.make_params(intercept=ydata.min(), slope=0)
    pars += peak_mod.guess(ydata, x=xdata)
    pars["sigma"].set(min=0)
    pars["amplitude"].set(min=0)

    fit = mod.fit(ydata, pars, x=xdata)
    print(f"Fitting with {model} model")
    if output:
        for key in fit.params:
            print(
                key, "=", fit.params[key].value, "+/-", fit.params[key].stderr
            )
        plt.plot(xdata, ydata)
        plt.plot(xdata, fit.best_fit)
        plt.show()
    return fit


def fit_series(
    db,
    scan_series,
    model="Gaussian",
    output=False,
    var_series=None,
    positioner="4C Theta",
    detector="APD",
    monitor=None,
    **kwargs,
):
    """
    Fit series of scan with reflection profile of choice and provide fit parameters as list.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).


    Parameters
    ----------
    db : database
        Databroker database. If None, it will attempt to read from csv or spec
        files.
    scan_series : int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    model:
        fit model: Gaussian, Lorentian, Voigt, PseidoVoigt
    output:
        Output fit parameters and plot data+fit for each scan.
    var_series:
        Varying variable for scan series, e.g. SampK (sample temperature), optional.
    positioner : string
        Name of the positioner, this needs to be the same as defined in
        Bluesky.
    detector : string
        Detector to be read, needs to be the same name as in Bluesky.
    monitor : string, optional
        Name of the monitor detector. The returned scans will be
        detectors/monitor.
    kwargs :
        Passed to `load_table` and 'fit_peak'.

    Returns
    -------
    fit : lmfit ModelResult class
        Contains the fit results. See:
        https://lmfit.github.io/lmfit-py/model.html#the-modelresult-class
    """

    if len(scan_series) % 3:
        raise ValueError(
            f"expected 3*n={3*(len(scan_series)//3)} arguments, got {len(scan_series)}"
        )

    nbp = 0
    for series in range(1, len(scan_series), 3):
        nbp = (
            nbp
            + int(scan_series[series] - scan_series[series - 1])
            / scan_series[series + 1]
            + 1
        )
    fit_result = [np.zeros(7)]
    for i in range(int(nbp) - 1):
        fit_result.append(np.zeros(7))

    index = 0
    for series in range(1, len(scan_series), 3):
        start = scan_series[series - 1]
        stop = scan_series[series]
        step = scan_series[series + 1]
        print("Intervals: {} to {} with step {}".format(start, stop, step))
        # fitnr=0
        for scan in range(start, stop + 1, step):

            table = load_table(
                scan, db, detector=detector, monitor=monitor, **kwargs
            )
            data = [np.array(table[positioner])]
            if monitor is None:
                data.append(np.array(table[detector]))
            else:
                data.append(
                    np.array(table[detector]) / np.array(table[monitor])
                )
            try:
                data.append(np.array(table[var_series]))
                fit_result[index][0] = data[2].mean()
                print(f"{var_series} = {fit_result[index][0]}")
            except:
                fit_result[index][0] = index

            fit = fit_peak(data[0], data[1], model=model, output=output)

            fit_result[index][1] = fit.params["amplitude"].value
            fit_result[index][2] = fit.params["amplitude"].stderr
            fit_result[index][3] = fit.params["center"].value
            fit_result[index][4] = fit.params["center"].stderr
            fit_result[index][5] = fit.params["fwhm"].value
            fit_result[index][6] = fit.params["fwhm"].stderr
            index += 1

    return DataFrame(
        fit_result,
        columns=[
            "Index",
            "Intensity",
            "sigma I",
            "Position",
            "sigma P",
            "Width",
            "sigma W",
        ],
    )