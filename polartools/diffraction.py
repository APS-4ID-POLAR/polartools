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
from os.path import join
from spec2nexus.spec import SpecDataFile
from .load_data import load_scan


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
    if output:
        for key in fit.params:
            print(f"Fitting with {model} model")
            print(
                key, "=", fit.params[key].value, "+/-", fit.params[key].stderr
            )
        plt.plot(xdata, ydata)
        plt.plot(xdata, fit.best_fit)
        plt.show()
    return fit


def load_info(db, spec_file, scan, folder, info):
    """
    Load metadata variable value.

    Parameters
    ----------
    db : database
        Databroker database. If None, it will attempt to read from csv or spec
        files.
    spec_file : string or spec2nexus.spec.SpecDataFile
        Either the spec file name or a SpecDataFile instance.
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id.
    folder : string, optional
        Folder where spec file is located.
    info: list
        Information on metadata to be read: List starting with #P, #U or #Q for
        motor positions, user values or Q-position:
        #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
        #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
        #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]

    Returns
    -------
    value : number
        Variable value.

    """

    if not db:
        if isinstance(spec_file, str):
            path = join(folder, spec_file)
            spec_file = SpecDataFile(path)
            print(spec_file)
        specscan = spec_file.getScan(scan)
        if isinstance(info, str):
            raise ValueError(
                "expect list [#P, (#Q, #U), Variable (string or line number), element number]"
            )

        if info[0] == "#P":
            data_array = specscan.P
            if isinstance(info[1], int):
                value = data_array[info[1]][info[2]]
            else:
                raise ValueError(
                    "For #P, expect row and column integer numbers"
                )

        elif info[0] == "#U":
            data_array = specscan.U
            if isinstance(info[1], str):
                for item in data_array:
                    ival = item.split(":")
                    if ival[0] == info[1]:
                        value = ival[info[2]].split()[0]
                        break
            else:
                raise ValueError("For #U, expect string and item number")

        elif info[0] == "#Q":
            data_array = specscan.Q
            value = data_array[info[2]]
        else:
            raise ValueError(
                "expect list [#P, (#Q, #U), Variable (string or line number), element number]"
            )

    else:
        # routine to load from db
        # To be implemented
        pass

    return value


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
    fit_result = [np.zeros(7) for i in range(int(nbp))]

    spec_file = kwargs.pop("spec_file", None)
    folder = kwargs.pop("folder", None)
    if spec_file and isinstance(spec_file, str):
        path = join(folder, spec_file)
        spec_file = SpecDataFile(path)

    index = 0
    for series in range(1, len(scan_series), 3):
        start = scan_series[series - 1]
        stop = scan_series[series]
        step = scan_series[series + 1]
        print("Intervals: {} to {} with step {}".format(start, stop, step))
        # fitnr=0
        for scan in range(start, stop + 1, step):
            if var_series and var_series[0][0] == "#":
                fit_result[index][0] = load_info(
                    db, spec_file, scan, folder, info=var_series
                )
                x, y = load_scan(
                    db,
                    scan,
                    positioner,
                    [detector],
                    monitor=[monitor],
                    spec_file=spec_file,
                    **kwargs,
                )
            elif var_series:
                x, y, parameter = load_scan(
                    db,
                    scan,
                    positioner,
                    [detector, var_series],
                    monitor=[monitor, None],
                    spec_file=spec_file,
                    **kwargs,
                )
                fit_result[index][0] = parameter.mean()
            else:
                x, y = load_scan(
                    db,
                    scan,
                    positioner,
                    [detector],
                    monitor=[monitor],
                    **kwargs,
                )
                fit_result[index][0] = index

            fit = fit_peak(x, y, model=model, output=output)

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
