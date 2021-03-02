"""
Functions to load and process x-ray diffraction data.

.. autosummary::
    ~fit_peak
    ~load_info
    ~fit_series
    ~load_series
    ~plot_2d
    ~plot_fit
"""

import numpy as np
from pandas import DataFrame
from lmfit.models import (
    GaussianModel,
    LorentzianModel,
    LinearModel,
    PseudoVoigtModel,
)
import matplotlib.pyplot as plt
from os.path import join
from spec2nexus.spec import SpecDataFile

from .load_data import load_table, load_csv, is_Bluesky_specfile


_spec_default_cols = dict(
    positioner="4C Theta",
    detector="APD",
    monitor="IC3",
)

_bluesky_default_cols = dict(
    positioner="fourc_theta",
    detector="APDSector4",
    monitor="Ion Ch 3",
)


def fit_peak(xdata, ydata, model="Gaussian", scan=None, output=False):
    """
    Fit Bragg peak with model of choice: Gaussian, Lorentzian, PseudoVoigt.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    Parameters
    ----------
    xdata : iterable
        List of x-axis values.
    yydata : iterable
        List of y-axis values.
    model: string
        fit model: Gaussian, Lorentzian, PseudoVoigt
    scan: integer
        scan number of scan fitted
    output: boolean, optional
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
        print(f"Fitting scan #{scan} with {model} model")
        for key in fit.params:
            print(
                key, "=", fit.params[key].value, "+/-", fit.params[key].stderr
            )
        plt.plot(xdata, ydata)
        plt.plot(xdata, fit.best_fit)
        plt.show()
    return fit


def load_info(source, scan_id, info, **kwargs):
    """
    Load metadata variable value.

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id.
    info : list
        If SPEC, information on metadata to be read: List starting with #P, #U
        or #Q for motor positions, user values or Q-position:

        - #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
        - #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
        - #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]

        If CSV, #metadata_name
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
    value : number
        Variable value.

    """
    folder = kwargs.pop("folder", "")
    if source == "csv":
        value = load_csv(
            scan_id, folder=folder, name_format="scan_{}_baseline.csv"
        )[info[1:]].mean()
        # to be implemented for csv
    elif isinstance(source, str) or isinstance(source, SpecDataFile):
        if isinstance(source, str):
            path = join(folder, source)
            source = SpecDataFile(path)
        specscan = source.getScan(scan_id)
        if isinstance(info, str):
            raise ValueError(
                "expect list [#P, (#Q, #U), Variable (string or line number), element number]"
            )

        if info[0] == "#P":
            data_array = specscan.P
            if isinstance(info[1], int):
                value = data_array[info[1]][int(info[2])]
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
                        value = ival[1].split()[int(info[2])]
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
        pass
        # to be implemented for database
    return value


def fit_series(
    source,
    scan_series,
    model="Gaussian",
    output=False,
    var_series=None,
    positioner=None,
    detector=None,
    monitor=None,
    normalize=False,
    **kwargs,
):
    """
    Fit series of scans with chosen functional and returns fit parameters.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    model : string, optional
        fit model: Gaussian, Lorentian, Voigt, PseidoVoigt
    output : boolean, optional
        Output fit parameters and plot data+fit for each scan.
    var_series : string or list
        If string

        - Varying variable for scan series to be read from scan (detector),
                e.g. SampK (sample temperature), optional.
        - String starting with #metadata, reads metadata from CSV baseline

        If list, information on metadata to be read: List starting with #P, #U
        or #Q for motor positions, user values or Q-position, optional:

        - #P ['#P', row, element_number], e.g. ['#P', 2, 0]
        - #U ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
        - #Q ['#Q', None, element_number], e.g. ['#Q', None, 0]

        If None, successive scans will be numbered starting from zero.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to '4C Theta' motor.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the APD detector.
    monitor : string, optional
        Name of the monitor detector. If None is passed, it defaults to the ion
        chamber 3.
    normalize : boolean, optional
        Normalization to selected/default monitor on/off
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
    fit : lmfit ModelResult class
        Contains the fit results. See:
        https://lmfit.github.io/lmfit-py/model.html#the-modelresult-class
    """
    # Select default parameters
    folder = kwargs.pop("folder", "")
    if isinstance(source, (str, SpecDataFile)) and source != "csv":
        if isinstance(source, str):
            path = join(folder, source)
            source = SpecDataFile(path)
        if is_Bluesky_specfile(source):
            _defaults = _bluesky_default_cols
        else:
            _defaults = _spec_default_cols

    else:
        _defaults = _bluesky_default_cols

    if not positioner:
        positioner = _defaults["positioner"]
    if not detector:
        detector = _defaults["detector"]
    if not monitor:
        monitor = _defaults["monitor"]

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
    fit_result = [np.zeros(8) for i in range(int(nbp))]

    index = 0

    for series in range(1, len(scan_series), 3):
        start = scan_series[series - 1]
        stop = scan_series[series]
        step = scan_series[series + 1]
        print("Intervals: {} to {} with step {}".format(start, stop, step))
        for scan in range(start, stop + 1, step):
            if var_series and var_series[0][0] == "#":
                fit_result[index][0] = load_info(
                    source, scan, info=var_series, folder=folder, **kwargs
                )
                fit_result[index][1] = 0
                table = load_table(
                    scan,
                    source,
                    folder=folder,
                    **kwargs,
                )
            elif var_series:
                table = load_table(
                    scan,
                    source,
                    folder=folder,
                    **kwargs,
                )
                fit_result[index][0] = table[var_series].mean()
                fit_result[index][1] = table[var_series].std()
            else:
                table = load_table(
                    scan,
                    source,
                    folder=folder,
                    **kwargs,
                )
                fit_result[index][0] = index
                fit_result[index][1] = 0
            x = table[positioner]
            y = table[detector]
            y0 = table[monitor]
            if normalize:
                y = y / y0

            fit = fit_peak(x, y, model=model, scan=scan, output=output)

            fit_result[index][2] = fit.params["amplitude"].value
            fit_result[index][3] = fit.params["amplitude"].stderr
            fit_result[index][4] = fit.params["center"].value
            fit_result[index][5] = fit.params["center"].stderr
            fit_result[index][6] = fit.params["fwhm"].value
            fit_result[index][7] = fit.params["fwhm"].stderr
            index += 1

    return DataFrame(
        fit_result,
        columns=[
            "Index",
            "Std Index",
            "Intensity",
            "Std I",
            "Position",
            "Std P",
            "Width",
            "Std W",
        ],
    )


def load_series(
    source,
    scan_series,
    log=False,
    var_series=None,
    positioner=None,
    detector=None,
    monitor=None,
    **kwargs,
):
    """
    Load series of scans as function of variable like temperature of field which can be
    spaced unequally. Generates input arrays for plot_2d.

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : list, int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    var_series: string or list
        string:
            - Varying variable for scan series to be read from scan (detector),
                e.g. SampK (sample temperature), optional.
            - String starting with #metadata, reads metadata from CSV baseline
        list: Information on metadata to be read: List starting with #P, #U or #Q for
            motor positions, user values or Q-position, optional:
            #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
            #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
            #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]
        If None, successive scans will be numbered starting from zero.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to '4C Theta' motor.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the APD detector.
    monitor : string, optional
        Name of the monitor detector for normalization. If None is passed, data are not normalized.
    log: boolean
        If True, z-axis plotted in logarithmic scale.
    kwargs:
        The necessary kwargs are passed to the loading and fitting functions defined by the
        `source` argument:
            - csv        -> possible kwargs: folder, name_format, e.g. "scan_{}_primary.csv"
            - spec       -> possible kwargs: folder
            - databroker -> possible kwargs: stream, query
        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

    Returns
    -------
    data : arrays with x, y and z information for 2D plot
    """
    nbp = 0
    for series in range(1, len(scan_series), 3):
        nbp = (
            nbp
            + int(scan_series[series] - scan_series[series - 1])
            / scan_series[series + 1]
            + 1
        )

    folder = kwargs.pop("folder", "")
    if isinstance(source, (str, SpecDataFile)) and source != "csv":
        if isinstance(source, str):
            path = join(folder, source)
            source = SpecDataFile(path)
        if is_Bluesky_specfile(source):
            _defaults = _bluesky_default_cols
        else:
            _defaults = _spec_default_cols

    else:
        _defaults = _bluesky_default_cols

    if not positioner:
        positioner = _defaults["positioner"]
    if not detector:
        detector = _defaults["detector"]
    if not monitor:
        monitor = _defaults["monitor"]
    if len(scan_series) % 3:
        raise ValueError(
            f"expected 3*n={3*(len(scan_series)//3)} arguments, got {len(scan_series)}"
        )

    table = load_table(
        scan_series[1],
        source,
        folder=folder,
        **kwargs,
    )
    data_len = len(table[detector])
    datax = [np.zeros(data_len) for i in range(int(nbp))]
    datay = [np.zeros(data_len) for i in range(int(nbp))]
    dataz = [np.zeros(data_len) for i in range(int(nbp))]

    index = 0
    for series in range(1, len(scan_series), 3):
        start = scan_series[series - 1]
        stop = scan_series[series]
        step = scan_series[series + 1]
        print("Intervals: {} to {} with step {}".format(start, stop, step))
        for scan in range(start, stop + 1, step):
            if var_series and var_series[0][0] == "#":
                table = load_table(
                    scan,
                    source,
                    folder=folder,
                    **kwargs,
                )
                y_value = load_info(
                    source, scan, info=var_series, folder=folder, **kwargs
                )
                tt = np.empty(data_len)
                tt.fill(y_value)
                datay[index] = tt
            elif var_series:
                table = load_table(
                    scan,
                    source,
                    folder=folder,
                    **kwargs,
                )
                datay[index] = table[var_series]
            else:
                table = load_table(
                    scan,
                    source,
                    folder=folder,
                    **kwargs,
                )
                tt = np.empty(data_len)
                tt.fill(index)
                datay[index] = tt
            datax[index] = table[positioner]
            if monitor:
                dataz[index] = table[detector] / table[monitor]
            else:
                dataz[index] = table[detector]

            if log:
                dataz[index][dataz[index] == 0] = 1
                dataz[index] = np.log10(dataz[index])
            index += 1
    return datax, datay, dataz


def plot_2d(
    source,
    scan_series,
    var_series=None,
    positioner=None,
    detector=None,
    monitor=None,
    log=False,
    **kwargs,
):
    """
    Plot 2d: Creates 2D plot of scans as function of variable parameter

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : list, int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
        e.g. [10,14,2,23,27,4] will use scan #10,12,14,23,27
    var_series: string or list
        string:
            - Varying variable for scan series to be read from scan (detector),
                e.g. SampK (sample temperature), optional.
            - String starting with #metadata, reads metadata from CSV baseline
        list: Information on metadata to be read: List starting with #P, #U or #Q for
            motor positions, user values or Q-position, optional:
            #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
            #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
            #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]
        If None, successive scans will be numbered starting from zero.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to '4C Theta' motor.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the APD detector.
    monitor : string, optional
        Name of the monitor detector for normalization. If None is passed, data are not normalized.
    log: boolean
        If True, z-axis plotted in logarithmic scale.
    output: string
        Output file for png file of plot.
    kwargs:
        The necessary kwargs are passed to the loading and fitting functions defined by the
        `source` argument:
            - csv        -> possible kwargs: folder, name_format, e.g. "scan_{}_primary.csv"
            - spec       -> possible kwargs: folder
            - databroker -> possible kwargs: stream, query
        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

    Returns
    -------
    2D plot, png-file

    """
    folder = kwargs.pop("folder", "")
    output = kwargs.pop("output", None)
    if isinstance(source, (str, SpecDataFile)) and source != "csv":
        if isinstance(source, str):
            path = join(folder, source)
            source = SpecDataFile(path)
        if is_Bluesky_specfile(source):
            _defaults = _bluesky_default_cols
        else:
            _defaults = _spec_default_cols

    else:
        _defaults = _bluesky_default_cols

    if not positioner:
        positioner = _defaults["positioner"]
    if not detector:
        detector = _defaults["detector"]
    if not monitor:
        monitor = _defaults["monitor"]

    datax, datay, dataz = load_series(
        source=source,
        scan_series=scan_series,
        log=log,
        var_series=var_series,
        positioner=positioner,
        detector=detector,
        monitor=monitor,
        **kwargs,
    )

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cmap = plt.get_cmap("rainbow")

    c = ax.pcolormesh(datax, datay, dataz, cmap=cmap, shading="auto")
    plt.colorbar(c)
    z_label = detector
    x_label = positioner
    if isinstance(var_series, list):
        y_label = " ".join(map(str, var_series))
    else:
        y_label = var_series
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))

    nlabel = ""
    for series in range(1, len(scan_series), 3):
        start = scan_series[series - 1]
        stop = scan_series[series]
        nlabel = nlabel + (", #{}-{}".format(start, stop))
    areas = len(scan_series) / 3
    if areas < 3:
        ax.set_title("{}{}".format(z_label, nlabel), fontsize=14)
    elif areas < 5:
        ax.set_title("{}{}".format(z_label, nlabel), fontsize=11)
    else:
        ax.set_title("{}{}".format(z_label, nlabel), fontsize=8)

    SIZE = 14
    plt.rc("font", size=SIZE)
    plt.rc("axes", titlesize=SIZE)
    plt.rc("axes", labelsize=SIZE)
    plt.rc("xtick", labelsize=SIZE)
    plt.rc("ytick", labelsize=SIZE)
    plt.rc("legend", fontsize=SIZE)

    if output:
        plt.savefig(output, dpi=600, transparent=True)
    plt.show()


def plot_fit(
    source,
    scan_series,
    model="Gaussian",
    output=False,
    var_series=None,
    positioner=None,
    detector=None,
    monitor=None,
    normalize=False,
    **kwargs,
):
    """
    Fit and plot series of scans with chosen functional and returns fit parameters.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    model : string, optional
        fit model: Gaussian, Lorentian, Voigt, PseidoVoigt
    output : boolean, optional
        Output fit parameters and plot data+fit for each scan.
    var_series : string or list
        If string

        - Varying variable for scan series to be read from scan (detector),
                e.g. SampK (sample temperature), optional.
        - String starting with #metadata, reads metadata from CSV baseline

        If list, information on metadata to be read: List starting with #P, #U
        or #Q for motor positions, user values or Q-position, optional:

        - #P ['#P', row, element_number], e.g. ['#P', 2, 0]
        - #U ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
        - #Q ['#Q', None, element_number], e.g. ['#Q', None, 0]

        If None, successive scans will be numbered starting from zero.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to '4C Theta' motor.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the APD detector.
    monitor : string, optional
        Name of the monitor detector. If None is passed, it defaults to the ion
        chamber 3.
    normalize : boolean, optional
        Normalization to selected/default monitor on/off
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
    fit : lmfit ModelResult class
        Contains the fit results. See:
        https://lmfit.github.io/lmfit-py/model.html#the-modelresult-class

    plot: plot of fitted data as function of variable changing between scans

    """

    data = fit_series(
        source,
        scan_series,
        output=output,
        var_series=var_series,
        positioner=positioner,
        detector=detector,
        monitor=monitor,
        normalize=normalize,
        **kwargs,
    )
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(3, 1, 1)
    ax1.errorbar(
        data["Index"],
        data["Intensity"],
        yerr=data["Std I"],
        xerr=data["Std Index"],
        color="orange",
        marker="o",
        linewidth=2,
        markersize=10,
    )
    ax1.set_ylabel("Intensity")
    ax2 = fig.add_subplot(3, 1, 2)
    ax2.errorbar(
        data["Index"],
        data["Position"],
        yerr=data["Std P"],
        xerr=data["Std Index"],
        color="blue",
        marker="o",
        linewidth=2,
        markersize=10,
    )
    ax2.set_ylabel("Position")
    ax3 = fig.add_subplot(3, 1, 3)
    ax3.errorbar(
        data["Index"],
        data["Width"],
        yerr=data["Std W"],
        xerr=data["Std Index"],
        color="green",
        marker="o",
        linewidth=2,
        markersize=10,
    )
    ax3.set_ylabel("FWHM")
    if isinstance(var_series, list):
        x_label = " ".join(map(str, var_series))
    else:
        x_label = var_series
    ax3.set_xlabel(x_label)

    print(data)
    plt.show()
