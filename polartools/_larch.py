# Functions hoisted from Matt Newville's xraylarch
# https://github.com/xraypy/xraylarch/tree/master/larch.
#
# We've had repeated issues with installing and testing polartools because of
# larch, but we use only a couple of rather simple functions.

from numpy import where, isfinite, gradient, abs


def finde0(energy, mu):
    if len(energy.shape) > 1:
        energy = energy.squeeze()
    if len(mu.shape) > 1:
        mu = mu.squeeze()

    dmu = gradient(mu)/gradient(energy)
    # find points of high derivative
    dmu[where(~isfinite(dmu))] = -1.0
    nmin = max(3, int(len(dmu)*0.05))
    maxdmu = max(dmu[nmin:-nmin])

    high_deriv_pts = where(dmu > maxdmu*0.1)[0]
    idmu_max, dmu_max = 0, 0

    for i in high_deriv_pts:
        if i < nmin or i > len(energy) - nmin:
            continue
        if ((dmu[i] > dmu_max and (i+1 in high_deriv_pts) and
           (i-1 in high_deriv_pts))):
            idmu_max, dmu_max = i, dmu[i]

    return energy[idmu_max]


def index_nearest(array, value):
    """
    return index of array *nearest* to value
    >>> ix = index_nearest(array, value)
    Arguments
    ---------
    array  (ndarray-like):  array to find index in
    value  (float): value to find index of
    Returns
    -------
    integer for index in array nearest value
    """
    return abs(array-value).argmin()
