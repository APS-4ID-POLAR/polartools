# Stuff hoisted from Mark Dean's pyrixs https://github.com/mpmdean/pyrixs.
#
# Not sure I will keep these since may need to operate in images rather
# than photon events.

import numpy as np
import lmfit


def poly(x, p2, p1, p0):
    """Third order polynominal function for fitting curvature.
    Returns p2*x**2 + p1*x + p0
    """
    return p2*x**2 + p1*x + p0


def image_to_photon_events(image):
    """ Convert 2D image into photon_events
    Parameters
    -----------
    image : np.array
        2D image
    Returns
    -----------
    photon_events : np.array
        three column x, y, I photon locations and intensities
    """
    X, Y = np.meshgrid(
        np.arange(image.shape[1]) + 0.5, np.arange(image.shape[0]) + 0.5
    )
    return np.vstack((X.ravel(), Y.ravel(), image.ravel())).transpose()


def bin_edges_centers(minvalue, maxvalue, binsize):
    """Make bin edges and centers for use in histogram
    The rounding of the bins edges is such that all bins are fully populated.
    Parameters
    -----------
    minvalue/maxvalue : array/array
        minimun/ maximum
    binsize : float (usually a whole number)
        difference between subsequnt points in edges and centers array
    Returns
    -----------
    edges : array
        edge of bins for use in np.histrogram
    centers : array
        central value of each bin. One shorter than edges
    """
    edges = binsize * np.arange(minvalue//binsize + 1, maxvalue//binsize)
    centers = (edges[:-1] + edges[1:])/2
    return edges, centers


def get_curvature_offsets(photon_events, binx=64, biny=0.5):
    """ Determine the offests that define the isoenergetic line.
    This is determined as the maximum of the cross correlation function with
    a reference taken from the center of the image.
    Parameters
    ------------
    photon_events : array
        three column x, y, I photon locations and intensities
    binx/biny : float/float
        width of columns/rows binned together prior to computing
        convolution. binx should be increased for noisy data.
    Returns
    -------------
    x_centers : array
        columns positions where offsets were determined
        i.e. binx/2, 3*binx/2, 5*binx/2, ...
    offests : array
        np.array of row offsets defining curvature. This is referenced
        to the center of the image.
    """
    x = photon_events[:, 0]
    y = photon_events[:, 1]
    intensity = photon_events[:, 2]
    x_edges, x_centers = bin_edges_centers(np.nanmin(x), np.nanmax(x), binx)
    y_edges, y_centers = bin_edges_centers(np.nanmin(y), np.nanmax(y), biny)

    H, _, _ = np.histogram2d(x, y, bins=(x_edges, y_edges), weights=intensity)

    ref_column = H[H.shape[0]//2, :]

    offsets = np.array([])
    for column in H:
        cross_correlation = np.correlate(column, ref_column, mode='same')
        offsets = np.append(offsets, y_centers[np.argmax(cross_correlation)])
    return x_centers, offsets - offsets[offsets.shape[0]//2]


def fit_poly(x_centers, offsets):
    """Fit curvature to vaues for curvature offsets.
    Parameters
    ----------
    x_centers, y_centers : float, float
        Shifts of the isoenergetic line as a function of column, x
    Returns
    --------
    result : lmfit result object
        object describing polynominal fit
    """
    poly_model = lmfit.Model(poly)
    params = poly_model.make_params()
    params['p0'].value = offsets[0]
    params['p1'].value = 0.
    params['p2'].value = 0.
    result = poly_model.fit(offsets, x=x_centers, params=params)
    if not result.success:
        print("Fitting failed")
    return result


def fit_curvature(photon_events, binx=32, biny=0.5, CONSTANT_OFFSET=500):
    """Get offsets, fit them and return polynomial that defines the curvature
    Parameters
    -------------
    photon_events : array
        two column x, y photon location coordinates
    binx/biny : float/float
        width of columns/rows binned together prior to computing
        convolution. binx should be increased for noisy data.
    CONSTANT_OFFSET : float
        offset is pass into last value of curvature
    Returns
    -----------
    """
    x_centers, offsets = get_curvature_offsets(
        photon_events, binx=binx, biny=biny
    )
    result = fit_poly(x_centers, offsets)
    curvature = np.array(
        [result.best_values['p2'], result.best_values['p1'], CONSTANT_OFFSET]
    )
    return curvature


def plot_curvature(ax1, curvature, photon_events):
    """ Plot a red line defining curvature on ax1
    Parameters
    ----------
    ax1 : matplotlib axes object
        axes for plotting on
    curvature : array
        n2d order polynominal defining image curvature
        np.array([x^2 coef, x coef, offset])
    photon_events : array
        two column x, y photon location coordinates
    Returns
    ---------
    curvature_artist : matplotlib artist object
        artist from image scatter plot
    """
    x = np.arange(np.nanmax(photon_events[:, 0]))
    y = poly(x, *curvature)
    return ax1.plot(x, y, 'r-')


def extract(photon_events, curvature, biny=0.5):
    """Apply curvature to photon events to create pixel versus intensity spectrum
    Parameters
    ----------
    photon_events : array
        two column x, y photon location coordinates
    curvature : array
        n2d order polynominal defining image curvature
        np.array([x^2 coef, x coef, offset])
    biny : float
        difference between subsequent points in the spectrum
    """
    x = photon_events[:, 0]
    y = photon_events[:, 1]
    intensity = photon_events[:, 2]
    corrected_y = y - poly(x, curvature[0], curvature[1], 0.)
    pix_edges, pix_centers = bin_edges_centers(
        np.nanmin(corrected_y), np.nanmax(corrected_y), biny
    )
    I, _ = np.histogram(corrected_y, bins=pix_edges, weights=intensity)
    spectrum = np.vstack((pix_centers, I)).transpose()
    return spectrum
