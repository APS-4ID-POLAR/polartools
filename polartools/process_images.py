"""
Functions to process image files.

.. autosummary::
   ~clean_threshold
   ~load_images
   ~get_curvature
   ~get_spectrum
   ~get_spectra
"""

import dask.array as da
import numpy as np
import matplotlib.pyplot as plt
from .load_data import load_hdf5_images, load_table, HDF_DEFAULT_FNAME_FORMAT
from ._pyrixs import (
    image_to_photon_events,
    plot_curvature,
    fit_curvature,
    extract,
)


def clean_threshold(images, threshold):
    """
    Cleans image stack using a fixed threshold value.

    Sets all value above the threshold to nan.

    Parameters
    ----------
    images : stack of images
        Images to be cleaned.
    threshold : float, integer
        Value of the threshold to be applied.

    Returns
    -------
    clean_images : stack of images
        Clean images.
    """
    clean_images = images.copy()
    clean_images[clean_images > threshold] = np.nan
    return clean_images


def _cleanup_images(images, parameters):
    """
    Clean up images.

    Parameters
    ----------
    images : stack of images
        Images to be cleaned.
    parameters : dictionary
        Clean up functions and arguments. Default options:
        - {'threshold': threshold_value}
        For custom functions, use:
        - {'function': (myfunc, (arg1, arg2, ...))}
        where myfunc is a function with call: myfunc(images, arg1, arg2, ...)
        These function can be stacked, but keys have to be different.
        For example, the call:
        {'threshold1': 100,
        'threshold2': 10,
        'function': (myfunc, (arg1, arg2))}
        will run the threshold function twice with 100 and 10 as argument, then
        myfunc with (arg1, arg2).

    Returns
    -------
    images : dask array
        Processed images.
    """

    for function, args in parameters.items():
        if "threshold" in function.lower():
            images = clean_threshold(images, *args)
        elif "function" in function.lower():
            images = args[0](images, *args[1:])
        else:
            raise ValueError(
                "function must be either 'threshold' for preset clean method, "
                "or 'function' for a custom clean method"
            )
    return images


def load_images(
    scans,
    cat,
    detector_key="lamb",
    cleanup=None,
    normalize=None,
    positioner=None,
    **kwargs,
):
    """
    Load scans with 2D images.

    Note that this function forces explicity outputs:

    - If multiple scans are passed, it will return the average of the scans.
    - If positioner = None, it will average all images. The resulting image \
        will have shape = (pixels_horizontal, pixels_vertical).
    - If positioner not None, it will not average all images. The resulting \
        image will have \
            shape = (positioner_size, pixels_horizontal, pixels_vertical).

    PARAMETERS
    ----------
    scans : iterable
        List of scan_id or uids.
    cat : databroker catalog, or "hdf5"/"h5"/"hdf"
        Catalog. If "hdf5", "h5", or "hdf", images and scan data are instead
        read from HDF5 master files on disk via
        :func:`polartools.load_data.load_hdf5_images` and
        :func:`polartools.load_data.load_table`. In that case, the following
        kwargs are used: folder, fname_format, externals_location,
        image_location (see :func:`polartools.load_data.load_hdf5_images`).
    detector_key : string, optional
        Name of item that holds the images. Defaults to "lamb".
    cleanup : dictionary, optional
        Clean up functions and arguments. Available functions:

        - {'threshold': threshold_value}

        For custom functions, use:

        - {'function': (myfunc, (arg1, arg2, ...))} \
        where myfunc is a function with call: myfunc(images, arg1, arg2, ...) \
        These function can be stacked. For example, the call: \
        {'threshold': 100, 'threshold': 10, 'function': (myfunc, (arg1, arg2))} \
        will run the threshold function twice with 100 and 10 as argument, then \
        myfunc with (arg1, arg2).

    normalize : string, optional
        Name of detector that will be used to normalize data. Default is None.
    positioner : string, optional
        Name of positioner to be read. Defaults to None.

    Returns
    -------
    images : dask array
        Processed images.
    positioner_values : numpy ndarray, optional
        Values of the positioner. It is only returned if positioner is not None.

    See also
    --------
    :func:`polartools.load_data.load_hdf5_images`
    :func:`polartools.load_data.load_table`
    """

    from_hdf5 = isinstance(cat, str) and cat in ("hdf5", "h5", "hdf")
    if from_hdf5:
        folder = kwargs.pop("folder", "")
        fname_format = kwargs.pop("fname_format", HDF_DEFAULT_FNAME_FORMAT)
        externals_location = kwargs.pop("externals_location", "entry/externals")
        image_location = kwargs.pop("image_location", "detector/data")

    output = []
    for scan in scans:
        if from_hdf5:
            data = load_table(
                scan, source=cat, folder=folder, fname_format=fname_format
            )
            images = (
                da.array(
                    load_hdf5_images(
                        scan,
                        folder,
                        detector_key,
                        fname_format=fname_format,
                        externals_location=externals_location,
                        image_location=image_location,
                    )
                )
                .compute()
                .astype(np.float64)
            )
            # Match the databroker layout: (n_points, frame=1, height, width)
            images = images.reshape(images.shape[0], 1, *images.shape[1:])
        else:
            data = cat[scan].primary.to_dask()
            images = da.array(data[detector_key]).compute().astype(np.float64)

        if cleanup is not None:
            if not isinstance(cleanup, dict):
                raise TypeError(
                    f"cleanup must be a dictionary, but a {type(cleanup)} "
                    "was entered."
                )

            images = _cleanup_images(images, cleanup)

        if normalize is not None:
            if from_hdf5:
                images[:] /= data[normalize].values[:, None, None, None]
            else:
                images[:] /= data[normalize].broadcast_like(data[detector_key])

        output.append(images)

    if positioner is None:
        return np.nanmean(da.stack(output), axis=(0, 1, 2))
    else:
        return (
            np.nanmean(da.stack(output), axis=(0, 2)),
            data[positioner].values,
        )


def _cleanup_photon_events(photon_events):
    index = np.where(np.isfinite(photon_events[:, 2]))
    return photon_events[index[0], :]


def _get_constant_offset(image, rng=10):
    timg = image.transpose()
    pos = int(timg.shape[1] / 2)
    col = np.nanmean(timg[:, pos - rng : pos + rng], axis=1)
    return np.where(col == np.nanmax(col))[0][0]


def get_curvature(image, binx=10, biny=1, constant_offset=None, plot=False):
    """
    Get the curvature of an image.

    A second order polynomial is fit to the cross correlation of each slice
    using the `pyrixs` algorithm. See pyrixs_ for more details. Ideally used
    with an image with large intensity and a single peak.

    .. _pyrixs: https://github.com/mpmdean/pyrixs/blob/master/pyrixs/process2d.py

    Parameters
    ----------
    image : numpy.array or dask.array
        Image to be used.
    binx : integer, optional
        Bin size in the horizontal direction in pixels. Defaults to 10.
    biny : integer, optional
        Bin size in the vertical direction in pixels. Defaults to 1.
    constant_offset : integer
        Offset for the curvature. This is only relevant for the plot, does not
        affect the use of the curvature to extract the spectra. Defaults to
        None, which triggers an automatic offset calculation.
    plot : boolean
        Flag to plot image and curvature.

    Returns
    -------
    curvature : list
        List of polynomial coefficients [c0, c1, c2] where:
        y = c0 + c1.x + c2.x^2.

    See also
    --------
    :func:`pyrixs.fit_curvature`
    :func:`pyrixs.get_curvature_offsets`
    """
    im = image.compute()
    ph = _cleanup_photon_events(image_to_photon_events(im.transpose()))
    if constant_offset is None:
        constant_offset = _get_constant_offset(im)
    curvature = fit_curvature(
        ph, binx=binx, biny=biny, CONSTANT_OFFSET=constant_offset
    )
    if plot:
        vmax = np.nanpercentile(im, 99.99)
        _, ax = plt.subplots()
        plt.pcolor(
            image.transpose(),
            vmin=0,
            vmax=vmax if vmax > 0 else np.nanmax(im),
            cmap="plasma",
        )
        plt.colorbar()
        plot_curvature(ax, curvature, ph)

    print(f"curvature = {curvature}")

    return curvature


def get_spectrum(image, curvature=None, biny=1):
    """
    Extract the spectrum of a single image.

    Uses the `pyrixs` extract_ function.

    .. _extract: \
        https://github.com/mpmdean/pyrixs/blob/e94531f325478ce18b88d9fb3d42068\
            f269bd118/pyrixs/process2d.py#L325

    Parameters
    ----------
    image : numpy.array or dask.array
        Image to be processed.
    curvature : iterable, optional
        List of polynomial coefficients [c0, c1, c2] where:
        y = c0 + c1.x + c2.x^2. Defaults to None, which skips the curvature
        correction and simply integrates the image in the vertical
        direction.
    biny : integer, optional
        Bin size in the vertical direction in pixels. Defaults to 1.

    Returns
    -------
    spectrum : numpy array
        Extracted 1D spectrum.

    See also
    --------
    :func:`pyrixs.extract`
    """
    if curvature is None:
        curvature = (0.0, 0.0, 0.0)
    if isinstance(image, da.core.Array):
        image = image.compute()
    ph = _cleanup_photon_events(image_to_photon_events(image.transpose()))
    return extract(ph, curvature, biny=biny)


def get_spectra(images, curvature=None, biny=1):
    """
    Extract several spectrum.

    Uses the `pyrixs` extract_ function.

    .. _extract: \
        https://github.com/mpmdean/pyrixs/blob/e94531f325478ce18b88d9fb3d42068\
            f269bd118/pyrixs/process2d.py#L325

    Parameters
    ----------
    images : iterable of numpy.array or dask.array
        Images to be processed.
    curvature : iterable, optional
        List of polynomial coefficients [c0, c1, c2] where:
        y = c0 + c1.x + c2.x^2. Defaults to None, which skips the curvature
        correction and simply integrates each image in the vertical
        direction.
    biny : integer, optional
        Bin size in the vertical direction in pixels. Defaults to 1.

    Returns
    -------
    spectra : numpy.array
        List of 1D spectrum.

    See also
    --------
    :func:`polartools.process_images.get_spectrum`
    :func:`pyrixs.extract`
    """
    spectra = []
    for image in images:
        if isinstance(image, da.core.Array):
            image = image.compute()
        spectra.append(get_spectrum(image, curvature, biny=biny))
    return np.array(spectra)


def process_rxes(
    scans,
    cat,
    detector_key,
    curvature=None,
    cleanup=None,
    normalize=None,
    positioner=None,
    biny=1,
    **kwargs,
):
    """
    Wrapper with typical RXES data processing.

    PARAMETERS
    ----------
    scans : iterable
        List of scan_id or uids.
    cat : databroker catalog
        Catalog.
    detector_key : string
        Name of item that holds the images
    curvature : iterable, optional
        List of polynomial coefficients [c0, c1, c2] where:
        y = c0 + c1.x + c2.x^2. Defaults to None, which skips the curvature
        correction and simply integrates each image in the vertical
        direction.
    cleanup : dictionary, optional
        Clean up functions and arguments. Available functions:

        - {'threshold': threshold_value}

        For custom functions, use:

        - {'function': (myfunc, (arg1, arg2, ...))} \
        where myfunc is a function with call: myfunc(images, arg1, arg2, ...) \
        These function can be stacked. For example, the call: \
        {'threshold': 100, 'threshold': 10, 'function': (myfunc, (arg1, arg2))} \
        will run the threshold function twice with 100 and 10 as argument, then \
        myfunc with (arg1, arg2).

    normalize : string, optional
        Name of detector that will be used to normalize data. Default is None.
    positioner : string, optional
        Name of positioner to be read. Defaults to None.
    biny : integer, optional
        Bin size in the vertical direction in pixels. Defaults to 1.

    Returns
    -------
    spectrum or spectra : array
        Processed spectrum (if positioner is None) or spectra (if positioner
        is not None).
    positioner_values : numpy ndarray, optional
        Values of the positioner. It is only returned if positioner is not None.
    """

    data = load_images(
        scans,
        cat,
        detector_key,
        cleanup=cleanup,
        normalize=normalize,
        positioner=positioner,
        **kwargs,
    )

    if positioner is None:  # "Count" scans
        return get_spectrum(data.compute(), curvature, biny=biny)
    else:  # Scans with a positioner
        return (get_spectra(data[0].compute(), curvature, biny=biny), data[1])


def process_rxes_mcd(
    scans,
    cat,
    detector_key,
    curvature=None,
    cleanup=None,
    normalize=None,
    positioner=None,
    biny=1,
    **kwargs,
):
    """
    Wrapper with typical RXES-MCD data processing.

    PARAMETERS
    ----------
    scans : iterable
        List of scan_id or uids.
    cat : databroker catalog
        Catalog.
    detector_key : string
        Name of item that holds the images
    curvature : iterable, optional
        List of polynomial coefficients [c0, c1, c2] where:
        y = c0 + c1.x + c2.x^2. Defaults to None, which skips the curvature
        correction and simply integrates each image in the vertical
        direction.
    cleanup : dictionary, optional
        Clean up functions and arguments. Available functions:
        - {'threshold': threshold_value}
        For custom functions, use:
        - {'function': (myfunc, (arg1, arg2, ...))} \
        where myfunc is a function with call: myfunc(images, arg1, arg2, ...) \
        These function can be stacked. For example, the call: \
        {'threshold': 100, 'threshold': 10, 'function': (myfunc, (arg1, arg2))} \
        will run the threshold function twice with 100 and 10 as argument, then \
        myfunc with (arg1, arg2).
    normalize : string, optional
        Name of detector that will be used to normalize data. Default is None.
    positioner : string, optional
        Name of positioner to be read. Defaults to None.
    biny : integer, optional
        Bin size in the vertical direction in pixels. Defaults to 1.
    Returns
    -------
    rxes : array
        Processed x-ray emission spectrum(if positioner is None) or spectra
        (if positioner is not None).
    mcd : array
        Processed x-ray emission dichroism spectrum (if positioner is None) or
        spectra (if positioner is not None).
    positioner_values : numpy ndarray, optional
        Values of the positioner. It is only returned if positioner is not None.
    """

    # Need to get all images...
    if positioner is None:
        positioner = "Time"

    images, positions = load_images(
        scans,
        cat,
        detector_key,
        cleanup=cleanup,
        normalize=normalize,
        positioner=positioner,
        **kwargs,
    )

    ims = images.reshape((-1, 4, images.shape[1], images.shape[2]))
    ims_plus = ims[:, [0, 3], :, :].mean(axis=1)
    ims_minus = ims[:, [1, 2], :, :].mean(axis=1)

    specs_plus = get_spectra(ims_plus.compute(), curvature, biny=biny)
    specs_minus = get_spectra(ims_minus.compute(), curvature, biny=biny)

    rxes = (specs_plus + specs_minus) / 2.0
    mcd = specs_plus.copy()
    mcd[:, :, 1] -= specs_minus[:, :, 1]

    if positioner == "Time":
        return rxes.mean(axis=0), mcd.mean(axis=0)
    else:
        return rxes, mcd, positions.reshape(-1, 4).mean(axis=1)
