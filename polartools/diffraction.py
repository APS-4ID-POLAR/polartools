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
    ~plot_data
    ~dbplot
"""

from enum import Enum
import numpy as np
from pandas import DataFrame
import lmfit.models
import matplotlib.pyplot as plt
import copy

plt.ion()
from os.path import join
from spec2nexus.spec import SpecDataFile
from xarray import DataArray

from .load_data import (
    load_table,
    load_csv,
    is_Bluesky_specfile,
    collect_meta,
    load_databroker,
)

rng = np.random.default_rng(seed=42)

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


class Model(Enum):
    Gaussian = lmfit.models.GaussianModel
    Lorentzian = lmfit.models.LorentzianModel
    PseudoVoigt = lmfit.models.PseudoVoigtModel


def fit_peak(xdata, ydata, model=Model.Gaussian):
    """
    Fit Bragg peak with model of choice: Gaussian, Lorentzian, PseudoVoigt.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    Parameters
    ----------
    xdata : iterable
        List of x-axis values.
    yydata : iterable
        List of y-axis values.
    model: enumeration
        fit model: model = Model.Gaussian (default), Model.Lorentzian, Model.PseudoVoigt

    Returns
    -------
    fit : lmfit ModelResult class
        Contains the fit results. See:
        https://lmfit.github.io/lmfit-py/model.html#the-modelresult-class
    """

    peak_mod = model.value()
    background = lmfit.models.LinearModel()
    mod = peak_mod + background
    pars = background.make_params(intercept=ydata.min(), slope=0)
    pars += peak_mod.guess(ydata, x=xdata)
    pars["sigma"].set(min=0)
    pars["amplitude"].set(min=0)

    fit = mod.fit(ydata, pars, x=xdata)

    return fit


def load_info(scan_id, info, source=None, **kwargs):
    """
    Load metadata variable value.

    Parameters
    ----------
    scan : int
        Scan_id our uid. If scan_id is passed, it will load the last scan with
        that scan_id.
    info : list
        If SPEC, information on metadata to be read: List starting with #P, #U,
        #Q for motor positions, user values or Q-position or general #xx:

        - #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
        - #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
        - #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]
        - #xx like #UA etc.: ['#UA', row, element_number],

        If CSV, #metadata_name
        If db, #baseline_information (e.g. #lakeshore340_sample)
    source: databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
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
        value = ""
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

        elif info[0][0] == "#":
            data_array = specscan.raw.split("\n")
            index = 0
            for element in data_array:
                if element[0 : len(info[0])] == info[0]:
                    if index == info[1]:
                        value = element.split()[info[2] + 1]
                    index += 1
        else:
            raise ValueError(
                "expect list [#P, (#Q, #U), Variable (string or line number), "
                "element number]"
            )
        if not value:
            raise ValueError(
                "expect list [#P, (#Q, #U, #xx), Variable (string or line number), "
                "element number]"
            )

    else:
        table = load_databroker(scan_id, db=source, stream="baseline")
        value = table[info[1 : len(info)]].mean()

    return value


def fit_series(
    scan_series,
    source=None,
    output=False,
    var_series=None,
    positioner=None,
    detector=None,
    monitor=None,
    normalize=False,
    xrange=None,
    **kwargs,
):
    """
    Fit series of scans with chosen functional and returns fit parameters.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    Parameters
    ----------
    scan_series : int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    output : boolean, optional
        Output fit parameters and plot data+fit for each scan.
    var_series : string or list
        If string:

            - Varying variable for scan series to be read from scan\
                {detector), e.g. SampK (sample temperature), optional.
            - String starting with #metadata, reads metadata from CSV baseline

        If list:
        SPEC: Information on metadata to be read: List starting with #P, #U,
            #Q for motor positions, user values or Q-position or general #xx:
            - #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
            - #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
            - #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]
            - #xx like #UA etc.: ['#UA', row, element_number],
        CSV: #metadata_name
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
    xrange: list
        Set positioner range for fitting
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

            - csv -> possible kwargs: folder, name_format.
            - spec -> possible kwargs: folder.
            - databroker -> possible kwargs: stream, query, use_db_v1.

        Note that a warning will be printed if the an unnecessary kwarg is
        passed.

        model: enumeration
            fit model: model = Model.Gaussian (default), Model.Lorentzian, Model.PseudoVoigt

    Returns
    -------
    fit : lmfit ModelResult class
        Contains the fit results. See:
        https://lmfit.github.io/lmfit-py/model.html#the-modelresult-class
    """
    # Select default parameters
    model = kwargs.pop("model", Model.Gaussian)
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
    fit_result = [np.zeros(9) for i in range(int(nbp))]
    if output:
        fig = plt.figure(num="Fitting", figsize=(6, 4), clear=True)
        ax = fig.add_subplot(1, 1, 1)

    index = 0
    for series in range(1, len(scan_series), 3):
        start = scan_series[series - 1]
        stop = scan_series[series]
        step = scan_series[series + 1]
        print("Intervals: {} to {} with step {}".format(start, stop, step))
        for scan in range(start, stop + 1, step):
            if var_series and var_series[0][0] == "#":
                fit_result[index][1] = load_info(
                    scan,
                    source=source,
                    info=var_series,
                    folder=folder,
                    **kwargs,
                )
                fit_result[index][2] = 0
                table = load_table(
                    scan,
                    source=source,
                    folder=folder,
                    **kwargs,
                )
            elif var_series:
                table = load_table(
                    scan,
                    source=source,
                    folder=folder,
                    **kwargs,
                )
                fit_result[index][1] = table[var_series].mean()
                fit_result[index][2] = table[var_series].std()
            else:
                table = load_table(
                    scan,
                    source=source,
                    folder=folder,
                    **kwargs,
                )
                fit_result[index][1] = index
                fit_result[index][2] = 0
            if (
                not isinstance(source, SpecDataFile)
                or isinstance(source, str)
                or source == "csv"
            ):
                positioner, detector, monitor = (
                    load_axes(
                        scan,
                        source=source,
                        positioner=positioner,
                        detector=detector,
                        monitor=monitor,
                        defaults=_defaults,
                        read=False,
                        **kwargs,
                    )
                    if positioner in table.columns
                    else load_axes(
                        scan,
                        source=source,
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
            table = table.set_index(positioner)
            if xrange:
                table = table.loc[xrange[0] : xrange[1]]
            x = table.index.to_numpy()
            y = table[detector].to_numpy()
            if normalize:
                y0 = table[monitor].to_numpy()
                y = y / y0

            fit = fit_peak(x, y, model=model)
            if output:
                print(f"Fitting scan #{scan} with {model} model")
                for key in fit.params:
                    print(
                        key,
                        "=",
                        fit.params[key].value,
                        "+/-",
                        fit.params[key].stderr,
                    )
                ax.plot(x, y, color=(f"C{index}"), label=f"#{scan}")
                ax.plot(
                    x,
                    fit.best_fit,
                    color=(f"C{index}"),
                    linestyle="dotted",
                )

            fit_result[index][0] = scan
            fit_result[index][3] = fit.params["amplitude"].value
            fit_result[index][4] = fit.params["amplitude"].stderr
            fit_result[index][5] = fit.params["center"].value
            fit_result[index][6] = fit.params["center"].stderr
            fit_result[index][7] = fit.params["fwhm"].value
            fit_result[index][8] = fit.params["fwhm"].stderr
            index += 1
    if output:
        ax.legend(loc=0)
        plt.get_current_fig_manager().show()

    return DataFrame(
        fit_result,
        columns=[
            "Scan #",
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
    scan_series,
    source=None,
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
    scan_series : list, int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    var_series : string or list
        string:

            - Varying variable for scan series to be read from scan\
                (detector), e.g. SampK (sample temperature), optional.
            - String starting with #metadata, reads metadata from CSV baseline

        list:
            SPEC: Information on metadata to be read: List starting with #P, #U,
                #Q for motor positions, user values or Q-position or general #xx:
                - #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
                - #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
                - #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]
                - #xx like #UA etc.: ['#UA', row, element_number],
            CSV: #metadata_name
        None, successive scans will be numbered starting from zero.
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
        source=source,
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
                scan_series[1],
                source=source,
                positioner=positioner,
                detector=detector,
                monitor=monitor,
                defaults=_defaults,
                read=False,
                **kwargs,
            )
            if positioner in table.columns
            else load_axes(
                scan_series[1],
                source=source,
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
                    source=source,
                    folder=folder,
                    **kwargs,
                )
                y_value = load_info(
                    scan,
                    source=source,
                    info=var_series,
                    folder=folder,
                    **kwargs,
                )
                tt = np.empty(data_len)
                tt.fill(y_value)
                datay[index] = tt
            elif var_series:
                table = load_table(
                    scan,
                    source=source,
                    folder=folder,
                    **kwargs,
                )
                datay[index] = table[var_series]
            else:
                table = load_table(
                    scan,
                    source=source,
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


def get_type(scan_id, source=None, **kwargs):
    """
    get_type returns type of scan and scan parameters.

    Parameters
    ----------
    scan_id : int
        scan number
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
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
    detector = _kwargs.pop("detector", "")
    scan_info = {
        "scan_no": 0,
        "plan_name": "rel_scan",
        "scan_type": None,
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
        if specscan is None:
            raise ValueError(f'Filename "{source}" not existing!')
        scan_cmd = specscan.scanCmd.split()
        plan_name = scan_cmd[0]
        scan_info["x0"] = scan_cmd[2]
        scan_info["x1"] = scan_cmd[3]
        scan_info["xint"] = scan_cmd[4]
        if plan_name == "mesh" or plan_name == "hklmesh":
            scan_info["plan_name"] = plan_name
            scan_info["y0"] = scan_cmd[6]
            scan_info["y1"] = scan_cmd[7]
            scan_info["yint"] = scan_cmd[8]
        elif plan_name == "dichromesh":
            scan_info["plan_name"] = plan_name
            scan_info["y0"] = scan_cmd[6]
            scan_info["y1"] = scan_cmd[7]
            scan_info["yint"] = scan_cmd[8]
            scan_info["scan_type"] = "dichro"
        elif plan_name != "qxscan":
            scan_info["plan_name"] = plan_name
        else:
            pass
    else:
        scan_read = collect_meta(
            [scan_id],
            ["plan_name", "plan_pattern_args", "num_points", "hints"],
            db=source,
        )
        for scanno, plan in scan_read.items():
            scan_info["scan_no"] = scanno
            for key, item in plan.items():
                if key == "plan_name":
                    scan_info["plan_name"] = item[0]
                if key == "num_points":
                    if not scan_info["xint"]:
                        scan_info["xint"] = item[0]
                if key == "plan_pattern_args":
                    if (
                        scan_info["plan_name"] == "grid_scan"
                        or scan_info["plan_name"] == "rel_grid_scan"
                    ):
                        scan_info["x0"] = item[0]["args"][1]
                        scan_info["x1"] = item[0]["args"][2]
                        scan_info["xint"] = item[0]["args"][3]
                        scan_info["y0"] = item[0]["args"][5]
                        scan_info["y1"] = item[0]["args"][6]
                        scan_info["yint"] = item[0]["args"][7]
                    else:
                        scan_info["x0"] = item[0]["args"][-2]
                        scan_info["x1"] = item[0]["args"][-1]

                if key == "hints":
                    scan_info["motor0"] = item[0]["dimensions"][0][0][0]
                    if (
                        scan_info["plan_name"] == "grid_scan"
                        or scan_info["plan_name"] == "rel_grid_scan"
                    ):
                        scan_info["motor1"] = item[0]["dimensions"][1][0][0]
                        scan_info["scan_type"] = item[0].get("scan_type", None)
                    scan_info["detector"] = (
                        detector if detector else item[0]["detectors"][0]
                    )
    return scan_info


def load_mesh(
    scan,
    scan_range,
    source=None,
    log=False,
    mrange="reduced",
    detector=None,
    **kwargs,
):
    """
    Load mesh generates input array for plot_2d from mesh scans:
        mesh, dichromesh, hklmesh (SPEC)
        grid_scan, rel_grid_scan (BlueSky)

    Parameters
    ----------
    scan_series : list, int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    scan_range : list, int
        scan parameters of mesh scan [x0, x1, xinterval, y0, y1, yinterval]
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    log: boolean
        If True, z-axis plotted in logarithmic scale.
    mrange: string, list
        reduced, full, [xmin,ymin,xmx,ymax]

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
    data = load_table(scan, source=source, **kwargs)
    if (
        scan_range["plan_name"] == "grid_scan"
        or scan_range["plan_name"] == "rel_grid_scan"
    ):
        x_label = scan_range["motor0"]
        y_label = scan_range["motor1"]
        if x_label == "nanopositioner_nanox":
            x_label = "nanopositioner_nanox_user_setpoint"
        if y_label == "nanopositioner_nanoy":
            y_label = "nanopositioner_nanoy_user_setpoint"
        z_label = scan_range["detector"]
        xr = int(scan_range["xint"])
        yr = int(scan_range["yint"])
    else:
        x_label = data.columns[0]
        y_label = data.columns[1]
        z_label = detector if detector else data.columns[-1]
        yr = int(scan_range["yint"]) + 1
        xr = int(scan_range["xint"]) + 1
    x = data[x_label]
    y = data[y_label]
    zp = data[z_label]
    ya = float(scan_range["y0"])
    yb = float(scan_range["y1"])
    ys = (yb - ya) / (yr - 1)
    if log:
        zp.replace(0, 1, inplace=True)
        zp = np.log10(zp)
    xi = x.unique()
    yi = y.unique()
    if xi.size > xr:
        xi = x[0 : xr * yr : yr].to_numpy()
        yi = y[0:yr:1].to_numpy()
    if yi.size < yr and mrange == "full":
        app = np.arange(yi[-1] + ys, yb, ys)
        yi = np.append(yi, app)
        z = np.zeros((xi.size * yi.size))
        z[: zp.size] = zp
        z[zp.size :] = np.nan
    else:

        data = data.groupby([y_label, x_label]).sum()[z_label]
    zp_left = data.unstack()
    yi = zp_left.index.values
    xi = zp_left.columns.values
    zi = zp_left.values

    return xi, yi, zi, x_label, y_label, z_label


def load_dichromesh(
    scan,
    scan_range,
    source=None,
    detector=None,
    **kwargs,
):
    """
    Load dichromesh generates input array for plot_2d from mesh scans with several polarizations per scan point:
        dichromesh (SPEC)
        dichro_grid_scan (BlueSky)

    Parameters
    ----------
    scan_series : list, int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    scan_range : list, int
        scan parameters of mesh scan [x0, x1, xinterval, y0, y1, yinterval]
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
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
    data : arrays with x, y and z information for 2D plot for left and right polarizatio and axes names
    """
    data = load_table(scan, source=source, **kwargs)
    if (
        scan_range["plan_name"] == "grid_scan"
        or scan_range["plan_name"] == "rel_grid_scan"
    ):
        x_label = scan_range["motor0"]
        y_label = scan_range["motor1"]
        if x_label == "nanopositioner_nanox":
            x_label = "nanopositioner_nanox_user_setpoint"
        if y_label == "nanopositioner_nanoy":
            y_label = "nanopositioner_nanoy_user_setpoint"
        z_label = scan_range["detector"]
    else:
        # needs to be tested for dichromesh (spec)
        x_label = data.columns[0]
        y_label = data.columns[1]
        z_label = detector if detector else data.columns[-1]
    data = (
        data.groupby([x_label, y_label, "pr2_pzt_localDC"])
        .sum()
        .unstack("pr2_pzt_localDC")[z_label]
    )
    left_p = data.columns[0]
    right_p = data.columns[1]
    zp_left = data[left_p].unstack()
    zp_right = data[right_p].unstack()
    yi = zp_left.index.values
    xi = zp_left.columns.values
    zl = zp_left.values
    zr = zp_right.values
    return xi, yi, zl, zr, x_label, y_label, z_label


def plot_2d(
    scans,
    source=None,
    var_series=None,
    positioner=None,
    detector=None,
    monitor=None,
    normalize=False,
    log=False,
    scale=None,
    mrange="reduced",
    direction=[1, 1],
    output=False,
    xcut=None,
    ycut=None,
    plot2=None,
    scale2=None,
    **kwargs,
):
    """
    Plot 2d:
    - Creates 2D plot from individual 1D scans as function of variable parameter
    - OR -
    - Plots a 2D mesh scan
        Supported mesh scans:
            mesh, dichromesh, hklmesh (SPEC)
            grid_scan, rel_grid_scan (BlueSky)


    Parameters
    ----------
    scans : list, int
        1D-scans:   start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
                    e.g. [10,14,2,23,27,4] will use scan #10,12,14,23,27
        mesh-scan:  scan_number, e.g. 10 will read scan #10
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
    var_series: string or list, optional
        string:
            - Varying variable for scan series to be read from scan
                (detector), e.g. SampK (sample temperature), optional.
            - String starting with #metadata, reads metadata from CSV baseline
        list:
        Information on metadata to be read: List starting with #P, #U or #Q for
        motor positions, user values or Q-position, or general #xx, optional:
            - #P: ['#P', row, element_number], e.g. ['#P', 2, 0]
            - #U: ['#U', Variable, element_number], e.g. ['#U', 'KepkoI', 1]
            - #Q: ['#Q', None, element_number], e.g. ['#Q', None, 0]
            - #xx like #UA etc.: ['#UA', row, element_number],
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
    log: boolean, optional
        If True, z-axis plotted in logarithmic scale.
    scale : list, int, optional
        intensity limits: [z_min,z_max]
    scale2 : list, int, optional
        intensity limits: [z_min,z_max] for sum if 2 plots are provided
    mrange: list, optional
        full, reduced, [xmin,ymin,xmax,ymax]
    direction : list, int, optional
        multiply axes for inversion: [1,-1]
    output : string, optional
        Output file for png file of plot.
    xcut : list, optional
        plot 1D cuts for x-axis values
    ycut : list, optional
        plot 1D cuts for y-axis values
    plot2 : list, int
        same as scans. Only z-values are taken. x and y used from scans
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
    dichro = False
    if isinstance(scans, int):
        scan_series = [scans, scans, 1]
    elif isinstance(scans, list):
        scan_series = scans
    else:
        raise ValueError(f"expected int or list got '{scans}'")
    scan_info = get_type(
        scan_id=scan_series[0], source=source, detector=detector, **kwargs
    )
    if (
        scan_info["plan_name"] == "mesh"
        or scan_info["plan_name"] == "dichromesh"
        or scan_info["plan_name"] == "hklmesh"
        or scan_info["plan_name"] == "grid_scan"
        or scan_info["plan_name"] == "rel_grid_scan"
    ):
        if scan_info["scan_type"] == "dichro":
            dichro = True
            (
                datax,
                datay,
                dataz,
                dataz2,
                positioner,
                var_series,
                detector,
            ) = load_dichromesh(
                scan_series[0],
                scan_info,
                source=source,
                log=log,
                mrange=mrange,
                detector=detector,
                **kwargs,
            )
        else:
            datax, datay, dataz, positioner, var_series, detector = load_mesh(
                scan_series[0],
                scan_info,
                source=source,
                log=log,
                mrange=mrange,
                detector=detector,
                **kwargs,
            )
            if plot2:
                _, _, dataz2, _, _, _ = load_mesh(
                    plot2,
                    scan_info,
                    source=source,
                    log=log,
                    mrange=mrange,
                    detector=detector,
                    **kwargs,
                )

    else:
        datax, datay, dataz, detector, positioner = load_series(
            scan_series=scan_series,
            source=source,
            log=log,
            var_series=var_series,
            positioner=positioner,
            detector=detector,
            monitor=monitor,
            normalize=normalize,
            **kwargs,
        )
        if plot2:
            _, _, dataz2, _, _ = load_series(
                scan_series=plot2,
                source=source,
                log=log,
                var_series=var_series,
                positioner=positioner,
                detector=detector,
                monitor=monitor,
                normalize=normalize,
                **kwargs,
            )
    if bool(xcut) ^ bool(ycut):
        fig = (
            plt.figure(
                num="Plot_2d - sum/diff - xcut/ycut",
                clear=True,
                figsize=(15, 4),
                constrained_layout=True,
            )
            if plot2
            else plt.figure(
                num="Plot_2d - xcut/ycut",
                clear=True,
                figsize=(10.3, 4),
                constrained_layout=True,
            )
        )
        ax = (
            fig.subplots(
                ncols=5, gridspec_kw={"width_ratios": [1, 0.05, 1, 0.05, 1]}
            )
            if plot2
            else fig.subplots(
                ncols=3, gridspec_kw={"width_ratios": [1, 0.05, 1]}
            )
        )
    elif xcut and ycut:
        fig = (
            plt.figure(
                num="Plot_2d - sum/diff - xcut and ycut",
                clear=True,
                figsize=(20, 4),
                constrained_layout=True,
            )
            if plot2
            else plt.figure(
                num="Plot_2d - xcut and ycut",
                clear=True,
                figsize=(15, 4),
                constrained_layout=True,
            )
        )
        ax = (
            fig.subplots(
                ncols=6, gridspec_kw={"width_ratios": [1, 0.05, 1, 0.05, 1, 1]}
            )
            if plot2
            else fig.subplots(
                ncols=4, gridspec_kw={"width_ratios": [1, 0.05, 1, 1]}
            )
        )
    elif plot2:
        fig = plt.figure(
            num="Plot_2d - sum/diff",
            clear=True,
            figsize=(11, 4),
            constrained_layout=True,
        )
        ax = fig.subplots(
            ncols=4, gridspec_kw={"width_ratios": [1, 0.05, 1, 0.05]}
        )
    elif dichro:
        fig = plt.figure(
            num="Plot_2d - dichro",
            clear=True,
            figsize=(20, 4),
            constrained_layout=True,
        )
        ax = fig.subplots(
            ncols=8,
            gridspec_kw={"width_ratios": [1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05]},
        )
    else:
        fig = plt.figure(
            num="Plot_2d",
            clear=True,
            figsize=(5.5, 4),
            constrained_layout=True,
        )
        ax = fig.subplots(ncols=2, gridspec_kw={"width_ratios": [1, 0.05]})

    cmap = plt.get_cmap("rainbow")
    datax = np.multiply(datax, direction[0])
    datay = np.multiply(datay, direction[1])
    datazm = np.subtract(dataz, dataz2) if plot2 else dataz

    if scale is None:
        scale = (np.nanpercentile(datazm, 1), np.nanpercentile(datazm, 99))
    vmin = float(scale[0])
    vmax = float(scale[1])
    z_label = detector
    x_label = positioner
    if isinstance(var_series, list):
        y_label = " ".join(map(str, var_series))
    else:
        y_label = var_series
    nlabel = ""
    if (
        scan_info["plan_name"] == "mesh"
        or scan_info["plan_name"] == "dichromesh"
        or scan_info["plan_name"] == "hklmesh"
        or scan_info["plan_name"] == "grid_scan"
        or scan_info["plan_name"] == "rel_grid_scan"
    ):
        nlabel = nlabel + (", #{}".format(scan_series[0]))
    else:
        for series in range(1, len(scan_series), 3):
            start = scan_series[series - 1]
            stop = scan_series[series]
            nlabel = nlabel + (", #{}-{}".format(start, stop))
    c = ax[0].pcolormesh(
        datax,
        datay,
        datazm,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        shading="auto",
    )
    if plot2:
        datazp = np.add(dataz2, dataz)

        if scale2 is None:
            scale2 = (
                np.nanpercentile(datazp, 1),
                np.nanpercentile(datazp, 99),
            )
        vmin = float(scale2[0])
        vmax = float(scale2[1])
        c1 = ax[2].pcolormesh(
            datax,
            datay,
            datazp,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            shading="auto",
        )
        plt.colorbar(c1, cax=ax[3])
        ax[2].set_xlabel(x_label)
        ax[2].set_ylabel(y_label)
        nlabel2 = (
            f", #{scan_series[0]}+#{plot2}"
            if (isinstance(scan_series, list) and isinstance(plot2, int))
            else f", #{scan_series}+#{plot2}"
        )
        nlabel = (
            f", #{scan_series[0]}-#{plot2}"
            if (isinstance(scan_series, list) and isinstance(plot2, int))
            else f", #{scan_series}-#{plot2}"
        )

    if dichro:
        scale = (np.nanpercentile(dataz2, 1), np.nanpercentile(dataz2, 99))
        vmin = float(scale[0])
        vmax = float(scale[1])
        c1 = ax[2].pcolormesh(
            datax,
            datay,
            dataz2,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            shading="auto",
        )
        plt.colorbar(c1, cax=ax[3])
        datazm = np.subtract(dataz, dataz2)
        scale = (np.nanpercentile(datazm, 1), np.nanpercentile(datazm, 99))
        vmin = float(scale[0])
        vmax = float(scale[1])

        c2 = ax[4].pcolormesh(
            datax,
            datay,
            datazm,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            shading="auto",
        )
        plt.colorbar(c2, cax=ax[5])
        datazm = np.add(dataz, dataz2)
        scale = (np.nanpercentile(datazm, 1), np.nanpercentile(datazm, 99))
        vmin = float(scale[0])
        vmax = float(scale[1])

        c3 = ax[6].pcolormesh(
            datax,
            datay,
            datazm,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            shading="auto",
        )
        plt.colorbar(c3, cax=ax[7])
        ax[2].set_xlabel(x_label)
        ax[4].set_xlabel(x_label)
        ax[6].set_xlabel(x_label)
        nlabel = f", #{scan_series[0]}: left circular"
        nlabel2 = f", #{scan_series[0]}: right circular"
        ax[4].set_title("difference", fontsize=12)
        ax[6].set_title("sum", fontsize=12)

    if len(z_label + nlabel) < 35:
        ax[0].set_title("{}{}".format(z_label, nlabel), fontsize=12)
        if plot2 or dichro:
            ax[2].set_title("{}{}".format(z_label, nlabel2), fontsize=12)
    elif len(z_label + nlabel) < 45:
        ax[0].set_title("{}{}".format(z_label, nlabel), fontsize=10)
        if plot2 or dichro:
            ax[2].set_title("{}{}".format(z_label, nlabel2), fontsize=10)
    else:
        ax[0].set_title("{}{}".format(z_label, nlabel), fontsize=8)
        if plot2 or dichro:
            ax[2].set_title("{}{}".format(z_label, nlabel2), fontsize=8)
    ax[0].set_xlabel(x_label)
    ax[0].set_ylabel(y_label)
    plt.colorbar(c, cax=ax[1])

    SIZE = 12
    plt.rc("font", size=SIZE)
    plt.rc("axes", titlesize=SIZE)
    plt.rc("axes", labelsize=SIZE)
    plt.rc("xtick", labelsize=SIZE)
    plt.rc("ytick", labelsize=SIZE)
    plt.rc("legend", fontsize=6)

    if datax.ndim == 2:
        datax = datax[0, :]
    if datay.ndim == 2:
        datay = datay[:, 0]
    data = DataArray(
        datazm,
        dims=("y_label", "x_label"),
        coords={"x_label": datax, "y_label": datay},
    )
    if plot2:
        data2 = DataArray(
            datazp,
            dims=("y_label", "x_label"),
            coords={"x_label": datax, "y_label": datay},
        )

    if xcut:
        num = 4 if plot2 else 2
        for pos in xcut:
            color = "#{:06x}".format(rng.integers(0, 16777215))
            test = data.sel(x_label=pos, method="nearest")
            ax[0].vlines(pos, datay[0], datay[-1], color=color)
            ax[num].plot(
                datay,
                test,
                color=color,
                label=(f"{x_label}={pos}"),
            )
            if plot2:
                color2 = "#{:06x}".format(rng.integers(0, 16777215))
                test = data2.sel(x_label=pos, method="nearest")
                ax[2].vlines(pos, datay[0], datay[-1], color=color2)
                ax[num].plot(
                    datay,
                    test,
                    color=color2,
                    label=(f"{x_label}={pos} (sum)"),
                )

            ax[num].set_xlabel(y_label)
            ax[num].set_ylabel(z_label)
            ax[num].legend(loc=0)
    if ycut:
        num = 4 if plot2 else 2
        if xcut:
            num += 1
        for pos in ycut:
            color = "#{:06x}".format(rng.integers(0, 16777215))
            test = data.sel(y_label=pos, method="nearest")
            ax[0].hlines(pos, datax[0], datax[-1], color=color)
            ax[num].plot(datax, test, color=color, label=(f"{y_label}={pos}"))
            if plot2:
                color2 = "#{:06x}".format(rng.integers(0, 16777215))
                test = data2.sel(y_label=pos, method="nearest")
                ax[2].hlines(pos, datax[0], datax[-1], color=color2)
                ax[num].plot(
                    datax, test, color=color2, label=(f"{y_label}={pos} (sum)")
                )
            ax[num].set_xlabel(x_label)
            ax[num].set_ylabel(z_label)
            ax[num].legend(loc=0)

    plt.get_current_fig_manager().show()
    if output:
        plt.savefig(output, dpi=600, transparent=True, bbox_inches="tight")


def plot_fit(
    scan_series,
    source=None,
    output=False,
    var_series=None,
    positioner=None,
    detector=None,
    monitor=None,
    normalize=False,
    errorbar=True,
    xrange=None,
    **kwargs,
):
    """
    Fit and plot series of scans with chosen functional and returns fit parameters.

    Uses lmfit (https://lmfit.github.io/lmfit-py/).

    Parameters
    ----------
    scan_series : int
        start, stop, step, [start2, stop2, step2, ... ,startn, stopn, stepn]
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
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
        - #xx like #UA etc.: ['#UA', row, element_number]

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
    xrange: list
        Set positioner range for fitting
    noerror : boolean, optional
        Plotting of errorbars on/off
    kwargs :
        The necessary kwargs are passed to the loading functions defined by the
        `source` argument:

        - csv -> possible kwargs: folder, name_format.
        - spec -> possible kwargs: folder.
        - databroker -> possible kwargs: stream, query.

        model : string, optional
            fit model: Gaussian, Lorentzian, PseudoVoigt

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
        scan_series,
        source=source,
        output=output,
        var_series=var_series,
        positioner=positioner,
        detector=detector,
        monitor=monitor,
        normalize=normalize,
        xrange=xrange,
        **kwargs,
    )
    fig = plt.figure(num="Plot_fit", figsize=(8, 8), clear=True)
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)
    if errorbar:
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
    else:
        ax1.plot(
            data["Index"],
            data["Intensity"],
            color="orange",
            marker="o",
            linewidth=2,
            markersize=10,
        )
        ax2.plot(
            data["Index"],
            data["Position"],
            color="blue",
            marker="o",
            linewidth=2,
            markersize=10,
        )
        ax3.plot(
            data["Index"],
            data["Width"],
            color="green",
            marker="o",
            linewidth=2,
            markersize=10,
        )
    ax1.set_ylabel("Intensity")
    ax2.set_ylabel("Position")
    ax3.set_ylabel("FWHM")
    if isinstance(var_series, list):
        x_label = " ".join(map(str, var_series))
    else:
        x_label = var_series
    ax3.set_xlabel(x_label)
    plt.get_current_fig_manager().show()

    return data


def load_axes(
    scan,
    source=None,
    positioner=None,
    detector=None,
    monitor=None,
    defaults=None,
    read=False,
    **kwargs,
):
    """
    Load default positioner and detector for plot

    Parameters
    ----------
    scan : int, list
        single scan
        or list [start, stop, step, start2, stop2, step2, ... ,startn, stopn, stepn]
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
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
        [scan], meta_keys=["motors", "hints"], db=source, query=query
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
    scan_series,
    source=None,
    positioner=None,
    detector=None,
    monitor=None,
    fit=False,
    normalize=False,
    log=False,
    deriv=False,
    direction=[1, 1],
    **kwargs,
):
    """
    Plot and fit data.

    Parameters
    ----------
    scan_series : list
        list [start, stop, step, start2, stop2, step2, ... ,startn, stopn, stepn]
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
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
    log: boolean
        If True, y-axis plotted in logarithmic scale.
    deriv : boolean, optional
        calculated derivative
    fit : boolean, optional
        Fitting of peak using model on/off. In case of deriv=True, derivative is fitted
    direction : list, int
        multiply axes for inversion: [1,-1] inverts y-axis
    kwargs :
        model: enumeration
            fit model: model = Model.Gaussian (default), Model.Lorentzian, Model.PseudoVoigt
        `source` argument:
            - csv -> possible kwargs: folder, name_format.
            - spec -> possible kwargs: folder.
            - databroker -> possible kwargs: stream, query, use_db_v1.

    Output
    -------
    Plot

    """

    model = kwargs.pop("model", Model.Gaussian)
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

    fig = plt.figure(
        num="Plot_data", figsize=(6, 4), constrained_layout=True, clear=True
    )
    ax = fig.add_subplot(1, 1, 1)
    if deriv:
        ax2 = ax.twinx()
    ax.clear()
    index = 0
    if isinstance(scan_series, list):
        if len(scan_series) % 3:
            raise ValueError(
                f"expected 3*n={3*((len(scan_series)+1)//3)} arguments, got {len(scan_series)}"
            )
        for series in range(1, len(scan_series), 3):
            start = scan_series[series - 1]
            stop = scan_series[series]
            step = scan_series[series + 1]
            print("Intervals: {} to {} with step {}".format(start, stop, step))
            for scan in range(start, stop + 1, step):
                data = load_table(
                    scan,
                    source=source,
                    **kwargs,
                )
                if (
                    not isinstance(source, SpecDataFile)
                    or isinstance(source, str)
                    or source == "csv"
                ):
                    positioner, detector, monitor = (
                        load_axes(
                            scan,
                            source=source,
                            positioner=positioner,
                            detector=detector,
                            monitor=monitor,
                            defaults=_defaults,
                            read=False,
                            **kwargs,
                        )
                        if positioner in data.columns
                        else load_axes(
                            scan,
                            source=source,
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

                data[positioner] = np.multiply(data[positioner], direction[0])
                data[detector] = np.multiply(data[detector], direction[1])
                if normalize:
                    data[detector] = data[detector] / data[monitor]
                if log:
                    data[detector].replace(0, 1, inplace=True)
                    data[detector] = np.log10(data[detector])
                x = data[positioner].to_numpy()
                y = data[detector].to_numpy()
                if deriv:
                    y = np.diff(y) / np.diff(x)
                    x = (x[:-1] + x[1:]) / 2
                if fit:
                    fit_data = fit_peak(x, y, model=model)
                    text1 = f"{fit_data.params['center'].value:.3f}"
                    text2 = f"{fit_data.params['fwhm'].value:.3f}"
                if deriv and fit:
                    ax.errorbar(
                        data[positioner],
                        data[detector],
                        color=(f"C{index}"),
                        marker="o",
                        linewidth=2,
                        markersize=8,
                        label=(f"#{scan}"),
                    )
                    ax2.plot(
                        x,
                        fit_data.best_fit,
                        color="black",
                        linewidth=2,
                    )
                    ax2.plot(
                        x,
                        y,
                        color=("#{:06x}".format(rng.integers(0, 16777215))),
                        marker="o",
                        linewidth=2,
                        markersize=8,
                        label=(f"#{scan}_deriv [{text1}, {text2}]"),
                    )
                elif fit:
                    ax.errorbar(
                        data[positioner],
                        data[detector],
                        color=(f"C{index}"),
                        marker="o",
                        linewidth=2,
                        markersize=8,
                        label=(f"#{scan} [{text1}, {text2}]"),
                    )
                    ax.plot(
                        x,
                        fit_data.best_fit,
                        color="black",
                        linewidth=2,
                    )
                elif deriv:
                    ax.errorbar(
                        data[positioner],
                        data[detector],
                        color=(f"C{index}"),
                        marker="o",
                        linewidth=2,
                        markersize=8,
                        label=(f"#{scan}"),
                    )
                    ax2.plot(
                        x,
                        y,
                        color=("#{:06x}".format(rng.integers(0, 16777215))),
                        marker="o",
                        linewidth=2,
                        markersize=8,
                        label=(f"#{scan}_deriv"),
                    )

                else:
                    ax.errorbar(
                        data[positioner],
                        data[detector],
                        color=(f"C{index}"),
                        marker="o",
                        linewidth=2,
                        markersize=8,
                        label=(f"#{scan}"),
                    )
                index += 1
    else:
        raise ValueError(f"expected list got '{scan_series}'")

    ax.set_xlabel(positioner)
    ax.set_ylabel(detector)
    ax.legend(loc=0)
    if deriv:
        ax2.legend(loc=4)
    plt.get_current_fig_manager().show()


def dbplot(
    scan,
    source=None,
    positioner=None,
    detector=None,
    monitor=None,
    normalize=False,
    fit=False,
    deriv=False,
    direction=[1, 1],
    **kwargs,
):
    """
    Plot and fit data.

    Parameters
    ----------
    scan_series : int, list
        single scan
        or list [start, stop, step, start2, stop2, step2, ... ,startn, stopn, stepn]
    source : databroker database, name of the spec file, or 'csv'
        Note that applicable kwargs depend on this selection.
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
    deriv : boolean, optional
        calculated derivative
    fit : boolean, optional
        Fitting of peak using model on/off. In case of deriv=True, derivative is fitted
    direction : list, int
        multiply axes for inversion: [1,-1] inverts y-axis
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
        scan_series=scan_series,
        source=source,
        positioner=positioner,
        detector=detector,
        monitor=monitor,
        fit=fit,
        normalize=normalize,
        deriv=deriv,
        direction=direction,
        **kwargs,
    )
