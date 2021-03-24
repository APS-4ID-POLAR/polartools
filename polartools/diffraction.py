"""
Functions to load and process x-ray diffraction data.

.. autosummary::
    ~fit_peak
    ~load_info
    ~fit_series
    ~load_series
    ~get_type
    ~load_mesh
    ~plot_2d
    ~plot_fit
    ~load_axes
    ~load_data
    ~dbplot
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
import copy

plt.ion()
from os.path import join
from spec2nexus.spec import SpecDataFile

from .load_data import load_table, load_csv, is_Bluesky_specfile
from .db_tools import collect_meta

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
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(xdata, ydata)
        ax.plot(xdata, fit.best_fit)
        plt.show(block=True)

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
        - databroker -> possible kwargs: stream, query, use_db_v1.

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
                "expect list [#P, (#Q, #U), Variable (string or line number), "
                "element number]"
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
                "expect list [#P, (#Q, #U), Variable (string or line number), "
                "element number]"
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
        fit model: Gaussian, Lorentzian, PseudoVoigt
    output : boolean, optional
        Output fit parameters and plot data+fit for each scan.
    var_series : string or list
        If string:

            - Varying variable for scan series to be read from scan\
                {detector), e.g. SampK (sample temperature), optional.
            - String starting with #metadata, reads metadata from CSV baseline

        If list:
        information on metadata to be read: List starting with #P, #U
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
            - databroker -> possible kwargs: stream, query, use_db_v1.

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

    if len(scan_series) % 3:
        raise ValueError(
            f"expected 3*n={3*(len(scan_series)//3)} arguments, got "
            f"{len(scan_series)}"
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
            if (
                not isinstance(source, SpecDataFile)
                or isinstance(source, str)
                or source == "csv"
            ):
                positioner, detector, monitor = (
                    load_axes(
                        source,
                        scan,
                        positioner=positioner,
                        detector=detector,
                        monitor=monitor,
                        defaults=_defaults,
                        read=False,
                        **kwargs,
                    )
                    if positioner in table.columns
                    else load_axes(
                        source,
                        scan,
                        positioner=positioner,
                        detector=detector,
                        monitor=monitor,
                        defaults=_defaults,
                        read=True,
                        **kwargs,
                    )
                )
            else:
                if not positioner:
                    positioner = _defaults["positioner"]
                if not detector:
                    detector = _defaults["detector"]
                if not monitor:
                    monitor = _defaults["monitor"]

            x = table[positioner].to_numpy()
            y = table[detector].to_numpy()
            if normalize:
                y0 = table[monitor].to_numpy()
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
    normalize=False,
    scale=None,
    **kwargs,
):
    """
    Load series of scans as function of variable like temperature of field
    which can be spaced unequally. Generates input arrays for plot_2d.

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : list, int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    var_series : string or list
        string:

            - Varying variable for scan series to be read from scan\
                (detector), e.g. SampK (sample temperature), optional.
            - String starting with #metadata, reads metadata from CSV baseline

        list:
        Information on metadata to be read: List starting with #P, #U or
        #Q for motor positions, user values or Q-position, optional.

            - #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
            - #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
            - #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]

        If None, successive scans will be numbered starting from zero.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to '4C Theta' motor.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the APD detector.
    monitor : string, optional
        Name of the monitor detector for normalization. If None is passed, data
        are not normalized.
    log : boolean
        If True, z-axis plotted in logarithmic scale.
    scale : list, int
        intensity limits: [z_min,z_max]

    kwargs :
        The necessary kwargs are passed to the loading and fitting functions
        defined by the `source` argument:

            - csv        -> possible kwargs: folder, name_format, e.g.\
                "scan_{}_primary.csv"
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
    if len(scan_series) % 3:
        raise ValueError(
            f"expected 3*n={3*(len(scan_series)//3)} arguments, got "
            f"{len(scan_series)}"
        )
    table = load_table(
        scan_series[1],
        source,
        folder=folder,
        **kwargs,
    )
    if (
        not isinstance(source, SpecDataFile)
        or isinstance(source, str)
        or source == "csv"
    ):
        positioner, detector, monitor = (
            load_axes(
                source,
                scan_series[1],
                positioner=positioner,
                detector=detector,
                monitor=monitor,
                defaults=_defaults,
                read=False,
                **kwargs,
            )
            if positioner in table.columns
            else load_axes(
                source,
                scan_series[1],
                positioner=positioner,
                detector=detector,
                monitor=monitor,
                defaults=_defaults,
                read=True,
                **kwargs,
            )
        )
    else:
        if not positioner:
            positioner = _defaults["positioner"]
        if not detector:
            detector = _defaults["detector"]
        if not monitor:
            monitor = _defaults["monitor"]

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
            if normalize:
                dataz[index] = table[detector] / table[monitor]
            else:
                dataz[index] = table[detector]
            if log:
                dataz[index].replace(0, 1, inplace=True)
                dataz[index] = np.log10(dataz[index])
            if scale:
                if len(scale) > 1:
                    dataz[index].values[
                        dataz[index] < float(scale[0])
                    ] = float(scale[0])
                    dataz[index].values[
                        dataz[index] > float(scale[1])
                    ] = float(scale[1])
                else:
                    dataz[index].values[
                        dataz[index] > float(scale[0])
                    ] = float(scale[0])
            index += 1
    return datax, datay, dataz, detector, positioner


def get_type(source, scan_id, **kwargs):
    """
    get_type returns type of scan and scan parameters.

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_id : int
        scan number

    kwargs :
        The necessary kwargs are passed to the loading and fitting functions
        defined by the `source` argument:
            - csv        -> possible kwargs: folder, name_format, e.g.\
                "scan_{}_primary.csv"
            - spec       -> possible kwargs: folder
            - databroker -> possible kwargs: stream, query

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

    Additional parameters in kwargs:
        folder: location of scan file

    Returns
    -------
    data : type of scan and scan parameters.

    """
    _kwargs = copy.deepcopy(kwargs)
    folder = _kwargs.pop("folder", "")
    scan_info = {
        "scan_no": 0,
        "scan_type": "rel_scan",
        "motor0": None,
        "x0": 0,
        "x1": 0,
        "xint": 0,
        "motor1": None,
        "y0": 0,
        "y1": 0,
        "yint": 0,
        "detector": None,
    }
    if source == "csv":
        pass
        # to be implemented for csv
    elif isinstance(source, str) or isinstance(source, SpecDataFile):
        if isinstance(source, str):
            path = join(folder, source)
            source = SpecDataFile(path)
        specscan = source.getScan(scan_id)
        scan_cmd = specscan.scanCmd.split()
        scan_type = scan_cmd[0]
        scan_info["x0"] = scan_cmd[2]
        scan_info["x1"] = scan_cmd[3]
        scan_info["xint"] = scan_cmd[4]
        if scan_type == "mesh" or scan_type == "hklmesh":
            scan_info["scan_type"] = scan_type
            scan_info["y0"] = scan_cmd[6]
            scan_info["y1"] = scan_cmd[7]
            scan_info["yint"] = scan_cmd[8]
        elif scan_type != "qxscan":
            scan_info["scan_type"] = scan_type
        else:
            pass
    else:
        scan_read = collect_meta(
            [scan_id],
            source,
            ["plan_name", "plan_pattern_args", "num_points", "hints"],
        )
        for scanno, plan in scan_read.items():
            scan_info["scan_no"] = scanno
            for key, item in plan.items():
                if key == "plan_name":
                    scan_info["scan_type"] = item[0]
                if key == "plan_pattern_args":
                    if scan_info["scan_type"] == "grid_scan":
                        scan_info["x0"] = item[0]["args"][1]
                        scan_info["x1"] = item[0]["args"][2]
                        scan_info["xint"] = item[0]["args"][3]
                        scan_info["y0"] = item[0]["args"][5]
                        scan_info["y1"] = item[0]["args"][6]
                        scan_info["yint"] = item[0]["args"][7]
                    else:
                        scan_info["x0"] = item[0]["args"][-2]
                        scan_info["x1"] = item[0]["args"][-1]
                        scan_info["xint"] = item[0]["num"]

                if key == "hints":
                    if scan_info["scan_type"] == "grid_scan":
                        scan_info["motor0"] = item[0]["dimensions"][0][0][0]
                        scan_info["motor1"] = item[0]["dimensions"][1][0][0]
                        scan_info["detector"] = item[0]["detectors"][0]
    return scan_info


def load_mesh(scan, source, scan_range, log=False, scale=None, **kwargs):
    """
    Load mesh generates input array for plot_2d from mesh_scan.

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : list, int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    scan_range : list, int
        scan parameters of mesh scan [x0, x1, xinterval, y0, y1, yinterval]
    log: boolean
        If True, z-axis plotted in logarithmic scale.
    scale : list, int
        intensity limits: [z_min,z_max]

    kwargs :
        The necessary kwargs are passed to the loading and fitting functions
        defined by the `source` argument:
            - csv        -> possible kwargs: folder, name_format, e.g.\
                "scan_{}_primary.csv"
            - spec       -> possible kwargs: folder
            - databroker -> possible kwargs: stream, query

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

    Additional parameters in kwargs:
        scale : list, int
            intensity limits: [z_min,z_max]

    Returns
    -------
    data : arrays with x, y and z information for 2D plot and axes names
    """

    data = load_table(scan=scan, source=source, **kwargs)
    if scan_range["scan_type"] == "grid_scan":
        x_label = scan_range["motor0"]
        y_label = scan_range["motor1"]
        z_label = scan_range["detector"]
        xr = int(scan_range["xint"])
        yr = int(scan_range["yint"])
    else:
        x_label = data.columns[0]
        y_label = data.columns[1]
        z_label = data.columns[-1]
        xr = int(scan_range["xint"]) + 1
        yr = int(scan_range["yint"]) + 1
    x1 = [x1[:] for x1 in [[1.01] * (xr)] * (yr)]
    y1 = [y1[:] for y1 in [[1.01] * (xr)] * (yr)]
    z1 = [z1[:] for z1 in [[1.01] * (xr)] * (yr)]
    r0 = data[x_label]
    r1 = data[y_label]
    r2 = data[z_label]
    for ii in range(0, yr):
        x1[ii] = r0[ii * xr : ii * xr + xr]
        y1[ii] = r1[ii * xr : ii * xr + xr]
        z1[ii] = r2[ii * xr : ii * xr + xr]
        if log:
            z1[ii].replace(0, 1, inplace=True)
            z1[ii] = np.log10(z1[ii])
        if scale:
            if len(scale) > 1:
                z1[ii].values[z1[ii] < float(scale[0])] = float(scale[0])
                z1[ii].values[z1[ii] > float(scale[1])] = float(scale[1])
            else:
                z1[ii].values[z1[ii] > float(scale[0])] = float(scale[0])
    return x1, y1, z1, x_label, y_label, z_label


def plot_2d(
    source,
    scan_series,
    var_series=None,
    positioner=None,
    detector=None,
    monitor=None,
    normalize=False,
    log=False,
    scale=None,
    direction=[1, 1],
    output=False,
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

            - Varying variable for scan series to be read from scan\
                (detector), e.g. SampK (sample temperature), optional.
            - String starting with #metadata, reads metadata from CSV baseline

        list:
        Information on metadata to be read: List starting with #P, #U or #Q for
        motor positions, user values or Q-position, optional:

            - #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
            - #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
            - #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]

        If None, successive scans will be numbered starting from zero.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to '4C Theta' motor.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the APD detector.
    monitor : string, optional
        Name of the monitor detector for normalization. If None is passed,
        data are not normalized.
    log: boolean
        If True, z-axis plotted in logarithmic scale.
    scale : list, int
        intensity limits: [z_min,z_max]
    direction : list, int
        multiply axes for inversion: [1,-1]
    output: string
        Output file for png file of plot.
    kwargs:
        The necessary kwargs are passed to the loading and fitting functions
        defined by the `source` argument:
            - csv        -> possible kwargs: folder, name_format, e.g.\
                "scan_{}_primary.csv"
            - spec       -> possible kwargs: folder
            - databroker -> possible kwargs: stream, query


    Returns
    -------
    2D plot, png-file

    """
    scan_info = get_type(source=source, scan_id=scan_series[0], **kwargs)
    if (
        scan_info["scan_type"] == "mesh"
        or scan_info["scan_type"] == "hklmesh"
        or scan_info["scan_type"] == "grid_scan"
    ):
        datax, datay, dataz, positioner, var_series, detector = load_mesh(
            scan_series[0], source, scan_info, log=log, scale=scale, **kwargs
        )

    else:
        datax, datay, dataz, detector, positioner = load_series(
            source=source,
            scan_series=scan_series,
            log=log,
            var_series=var_series,
            positioner=positioner,
            detector=detector,
            monitor=monitor,
            normalize=normalize,
            scale=scale,
            **kwargs,
        )
    plt.close("all")
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cmap = plt.get_cmap("rainbow")

    datax = np.multiply(datax, direction[0])
    datay = np.multiply(datay, direction[1])
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
    if (
        scan_info["scan_type"] == "mesh"
        or scan_info["scan_type"] == "hklmesh"
        or scan_info["scan_type"] == "grid_scan"
    ):
        nlabel = nlabel + (", #{}".format(scan_series[0]))
    else:
        for series in range(1, len(scan_series), 3):
            start = scan_series[series - 1]
            stop = scan_series[series]
            nlabel = nlabel + (", #{}-{}".format(start, stop))
    areas = len(scan_series) / 3
    if areas < 3:
        ax.set_title("{}{}".format(z_label, nlabel), fontsize=12)
    elif areas < 5:
        ax.set_title("{}{}".format(z_label, nlabel), fontsize=10)
    else:
        ax.set_title("{}{}".format(z_label, nlabel), fontsize=8)

    SIZE = 12
    plt.rc("font", size=SIZE)
    plt.rc("axes", titlesize=SIZE)
    plt.rc("axes", labelsize=SIZE)
    plt.rc("xtick", labelsize=SIZE)
    plt.rc("ytick", labelsize=SIZE)
    plt.rc("legend", fontsize=SIZE)

    if output:
        plt.savefig(output, dpi=600, transparent=True, bbox_inches="tight")


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
        fit model: Gaussian, Lorentzian, PseudoVoigt
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

    plt.close("all")
    data = fit_series(
        source,
        scan_series,
        model=model,
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


def load_axes(
    source,
    scan,
    positioner=None,
    detector=None,
    monitor=None,
    defaults=None,
    read=False,
    **kwargs,
):
    """
    Plot and fit data.

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : int, list
        single scan
        or list [start, stop, step, start2, stop2, step2, ... ,startn, stopn, stepn]
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to '4C Theta' motor.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the APD detector.
    monitor : string, optional
        Name of the monitor detector. If None is passed, it defaults to the ion
        chamber 3.
    defaults : string, optional
        Default values for positioner and detector
    read: boolean, optional
        Determines if positioner is read from metadata

    Output
    -------
    Plot

    """
    _kwargs = copy.deepcopy(kwargs)
    query = _kwargs.pop("query", None)
    meta = collect_meta(
        [scan], source, meta_keys=["motors", "hints"], query=query
    )
    if not positioner or read:
        positioner = meta[scan]["motors"][0]
    det = meta[scan]["hints"] if "hints" in meta[scan] else None
    if not detector:
        detector = (
            det[0]["detectors"][0]
            if "detectors" in det[0]
            else defaults["detector"]
        )
    if not monitor:
        monitor = (
            det[0]["monitor"] if "monitor" in det[0] else defaults["monitor"]
        )
    return positioner, detector, monitor


def plot_data(
    source,
    scan_series,
    positioner=None,
    detector=None,
    monitor=None,
    fit=False,
    normalize=False,
    **kwargs,
):
    """
    Plot and fit data.

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : list
        list [start, stop, step, start2, stop2, step2, ... ,startn, stopn, stepn]
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky.
    monitor : string, optional
        Monitor to be read from this scan, again it needs to be the same name
        as in Bluesky.
    normalize : boolean, optional
        Normalization to selected/default monitor on/off
    fit : boolean, optional
        Fitting of peak using model on/off
    kwargs :
        model : string, optional
            - fit model: Gaussian, Lorentzian, PseudoVoigt
        `source` argument:
            - csv -> possible kwargs: folder, name_format.
            - spec -> possible kwargs: folder.
            - databroker -> possible kwargs: stream, query, use_db_v1.

    Output
    -------
    Plot

    """

    model = kwargs.pop("model", "Gaussian")
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

    plt.close("all")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.clear()
    index = 0
    if isinstance(scan_series, list):
        if len(scan_series) % 3:
            raise ValueError(
                f"expected 3*n={3*(len(scan_series)//3)} arguments, got {len(scan_series)}"
            )
        for series in range(1, len(scan_series), 3):
            start = scan_series[series - 1]
            stop = scan_series[series]
            step = scan_series[series + 1]
            print("Intervals: {} to {} with step {}".format(start, stop, step))
            for scan in range(start, stop + 1, step):
                data = load_table(
                    scan,
                    source,
                    **kwargs,
                )
                positioner, detector, monitor = (
                    load_axes(
                        source,
                        scan,
                        positioner=positioner,
                        detector=detector,
                        monitor=monitor,
                        defaults=_defaults,
                        read=False,
                        **kwargs,
                    )
                    if positioner in data.columns
                    else load_axes(
                        source,
                        scan,
                        positioner=positioner,
                        detector=detector,
                        monitor=monitor,
                        defaults=_defaults,
                        read=True,
                        **kwargs,
                    )
                )
                if normalize:
                    data[detector] = data[detector] / data[monitor]
                if fit:
                    x = data[positioner].to_numpy()
                    y = data[detector].to_numpy()
                    fit_data = fit_peak(
                        x, y, model=model, scan=None, output=False
                    )
                    text1 = f"{fit_data.params['center'].value:.3f}"
                    text2 = f"{fit_data.params['fwhm'].value:.3f}"
                    ax.plot(
                        x,
                        fit_data.best_fit,
                        color="black",
                        linewidth=2,
                    )

                    ax.errorbar(
                        data[positioner],
                        data[detector],
                        color=(f"C{index}"),
                        marker="o",
                        linewidth=2,
                        markersize=10,
                        label=(f"#{scan} [{text1}, {text2}]"),
                    )
                else:
                    ax.errorbar(
                        data[positioner],
                        data[detector],
                        color=(f"C{index}"),
                        marker="o",
                        linewidth=2,
                        markersize=10,
                        label=(f"#{scan}"),
                    )
                index += 1
    else:
        raise ValueError(f"expected list got '{scan_series}'")

    ax.set_xlabel(positioner)
    ax.set_ylabel(detector)
    ax.legend()
    plt.show(block=False)


def dbplot(
    source,
    scan,
    positioner=None,
    detector=None,
    monitor=None,
    normalize=False,
    fit=False,
    **kwargs,
):
    """
    Plot and fit data.

    Parameters
    ----------
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    scan_series : int, list
        single scan
        or list [start, stop, step, start2, stop2, step2, ... ,startn, stopn, stepn]
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
    fit : boolean, optional
        Fitting of peak using model on/off
    kwargs :
        model : string, optional
            - fit model: Gaussian, Lorentzian, PseudoVoigt
        `source` argument:
            - csv -> possible kwargs: folder, name_format.
            - spec -> possible kwargs: folder.
            - databroker -> possible kwargs: stream, query, use_db_v1.

    Output
    -------
    Plot

    """

    if isinstance(scan, int):
        scan_series = [scan, scan, 1]
    elif isinstance(scan, list):
        scan_series = scan
    else:
        raise ValueError(f"expected int or list got '{scan}'")
    plot_data(
        source=source,
        scan_series=scan_series,
        positioner=positioner,
        detector=detector,
        monitor=monitor,
        fit=fit,
        normalize=normalize,
        **kwargs,
    )
