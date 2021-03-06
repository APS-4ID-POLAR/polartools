
def plot_data(source,
    scan_series,
    positioner=None,
    detector=None,
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
    output : boolean, optional
        Output fit parameters and plot data+fit for each scan.
    positioner : string, optional
        Name of the positioner, this needs to be the same as defined in
        Bluesky or SPEC. If None is passed, it defauts to '4C Theta' motor.
    detector : string, optional
        Detector to be read from this scan, again it needs to be the same name
        as in Bluesky. If None is passed, it defaults to the APD detector.
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
    
    model= kwargs.pop("model","Gaussian")
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

    plt.close('all')
    index=0
    if isinstance(scan_series,int):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1,1,1)
        ax.clear()
        data = load_table(scan_series,
        source,
        **kwargs,
        )
        if len(data.columns)==0:
            raise ValueError(f"No data in scan")
        meta=collect_meta([scan_series],source,meta_keys=['motors','hints'])
        positioner=meta[scan_series]['motors'][0]
        det=meta[scan_series]['hints'] if 'hints' in meta[scan_series] else None
        if det:
            detector=det[0]['detectors'][0] 

        if fit:
            x=data[positioner].to_numpy()
            y=data[detector].to_numpy()
            #print(x,y)
            fit_data=fit_peak(x,y, model=model, scan=None, output=False)
            ax.plot(
            x,
            fit_data.best_fit,
            color="black",
            linewidth=2,
            )
            text1=(f"{fit_data.params['center'].value:.3f}")
            text2=(f"{fit_data.params['fwhm'].value:.3f}")
            ax.errorbar(
            data[positioner],
            data[detector],
            color="orange",
            marker="o",
            linewidth=2,
            markersize=10,
            label=(f'#{scan_series} [{text1}, {text2}]')
            )
        else:
            ax.errorbar(
            data[positioner],
            data[detector],
            color="orange",
            marker="o",
            linewidth=2,
            markersize=10,
            label=(f'#{scan_series}')
            )

    elif isinstance(scan_series,list):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1,1,1)
        ax.clear()
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
                data = load_table(scan,
                source,
                **kwargs,
                )
                meta=collect_meta([scan],source,meta_keys=['motors','hints'])
                positioner=meta[scan]['motors'][0]
                det=meta[scan]['hints'] if 'hints' in meta[scan] else None
                if det:
                    detector=det[0]['detectors'][0] 
                if fit:
                    x=data[positioner].to_numpy()
                    y=data[detector].to_numpy()
                    fit_data=fit_peak(x,y, model=model, scan=None, output=False)
                    text1=(f"{fit_data.params['center'].value:.3f}")
                    text2=(f"{fit_data.params['fwhm'].value:.3f}")
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
                    label=(f'#{scan} [{text1}, {text2}]')
                    )
                else:
                    ax.errorbar(
                    data[positioner],
                    data[detector],
                    color=(f"C{index}"),
                    marker="o",
                    linewidth=2,
                    markersize=10,
                    label=(f'#{scan}')
                )
                index+=1
    else:
        raise ValueError(
            f"expected str or list got '{scan_series}'"
        )

    ax.set_xlabel(positioner)
    ax.set_ylabel(detector)
    ax.legend()
    plt.show(block=False)

