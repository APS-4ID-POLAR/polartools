import dask.array as da
from numpy import nan, array, float64
import numpy as np
import matplotlib.pyplot as plt
from ._pyrixs import (
    image_to_photon_events, plot_curvature, fit_curvature, extract
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
    clean_images[clean_images > threshold] = nan
    return clean_images


def _cleanup_images(images, parameters):
    """
    Clean up images

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
        These function can be stacked. For example, the call:
        {'threshold': 100, 'threshold': 10, 'function': (myfunc, (arg1, arg2))}
        will run the threshold function twice with 100 and 10 as argument, then
        myfunc with (arg1, arg2).

    Returns
    -------
    images : dask array
        Processed images.
    """

    for function, args in parameters.items():
        if function.lower() == 'threshold':
            images = clean_threshold(images, *args)
        else:
            func = args.pop('function', None)
            if func is None:
                raise ValueError(
                    'The custom cleanup function must be passed in '
                    'the cleanup dictionary using the "fuction" key.'
                )
            images = func[0](images, *func[1])


def load_images(scans, cat, detector_key, cleanup=None, normalize=None,
                positioner=None):
    """
    Load scans with 2D images.

    If multiple scans are passed, it will return the average.

    PARAMETERS
    ----------
    scans : iterable
        List of scan_id or uids.
    cat : databroker catalog
        Catalog.
    detector_key : string
        Name of item that holds the images
    cleanup : dictionary, optional
        Clean up functions and arguments. Available functions:

        - {'threshold': threshold_value}
        For custom functions, use:
        - {'function': (myfunc, (arg1, arg2, ...))}
        where myfunc is a function with call: myfunc(images, arg1, arg2, ...)
        These function can be stacked. For example, the call:
        {'threshold': 100, 'threshold': 10, 'function': (myfunc, (arg1, arg2))}
        will run the threshold function twice with 100 and 10 as argument, then
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
    """

    output = []
    for scan in scans:
        data = cat[scan].primary.to_dask()
        images = da.array(data[detector_key].astype(float64)).compute()

        if cleanup is not None:
            if not isinstance(cleanup, dict):
                raise TypeError(
                    f"cleanup must be a dictionary, but a {type(cleanup)} "
                    "was entered."
                )

            images = _cleanup_images(images, cleanup)

        if normalize is not None:
            images[:] /= data[normalize].broadcast_like(data[detector_key])

        output.append(images)

    if positioner is None:
        return np.nanmean(da.stack(output), axis=0)
    else:
        return np.nanmean(da.stack(output), axis=0), data[positioner].values


def _cleanup_photon_events(photon_events):
    index = np.where(np.isfinite(photon_events[:, 2]))
    return photon_events[index[0], :]


def _get_constant_offset(image, rng=10):
    timg = image.transpose()
    pos = int(timg.shape[1]/2)
    col = np.nanmean(timg[:, pos-rng:pos+rng], axis=1)
    return np.where(col == np.nanmax(col))[0][0]


def get_curvature(image, binx=10, biny=1, constant_offset=None, plot=False):
    ph = _cleanup_photon_events(image_to_photon_events(image.transpose()))
    if constant_offset is None:
        constant_offset = _get_constant_offset(image)
    curv = fit_curvature(
        ph, binx=binx, biny=biny, CONSTANT_OFFSET=constant_offset
    )
    if plot:
        _, ax = plt.subplots()
        plt.pcolor(
            image.transpose(),
            vmin=0,
            vmax=np.nanpercentile(image, 99),
            cmap="plasma")
        plt.colorbar()
        plot_curvature(ax, curv, ph)

    print(f'curvature = {curv}')

    return curv


def get_spectrum(image, curvature, biny=1):
    ph = _cleanup_photon_events(image_to_photon_events(image.transpose()))
    return extract(ph, curvature, biny=biny)


def get_spectra(images, curvature, biny=1):
    spectra = []
    for image in images:
        spectra.append(get_spectrum(image[0].compute(), curvature, biny=biny))
    return array(spectra)
