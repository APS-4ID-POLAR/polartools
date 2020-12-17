import numpy as np


def normalize_absorption(energy, xanes, pre_edge_range, pos_edge_range, e0,
                         pre_edge_order=1, pos_edge_order=1):
    """
    Extract pre- and post-edge normalization curves by fitting polynomials.

    Parameters
    ----------
    energy : list
        Incident energy.
    xanes : list
        X-ray absorption.
    pre_edge_range : list
        List with the energy ranges [initial, final] of the pre-edge region
        **relative** to the absorption edge.
    pos_edge_range : list
        List with the energy ranges [initial, final] of the post-edge region
        **relative** to the absorption edge.
    e0 : float
        Absorption edge energy.
    pre_edge_order : int, optional
        Order of the polynomial to be used in the pre-edge. Defauts to 1.
    pos_edge_order : int, optional
        Order of the polynomial to be used in the post-edge. Defauts to 1.

    Returns
    -------
    pre_edge : numpy.array
        Pre-edge polynomial.
    pos_edge : numpy.array
        Post-edge polynomial.
    jump : float
        Size of the absorption jump.

    See also
    --------
    `numpy.polyfit`
    """

    energy = np.array(energy)
    xanes = np.array(xanes)

    # Process pre-edge
    index = (energy > pre_edge_range[0]) & (energy < pre_edge_range[1])
    pre_edge = np.poly1d(np.polyfit(energy[index], xanes[index],
                                    pre_edge_order))(energy)
    processed_xanes = xanes - pre_edge

    # Process pos-edge
    index = (energy > pos_edge_range[0]) & (energy < pos_edge_range[1])
    pos_edge_func = np.poly1d(np.polyfit(energy[index], processed_xanes[index],
                              pos_edge_order))
    pos_edge = pos_edge_func(energy)
    jump = pos_edge_func(e0)

    return pre_edge, pre_edge + pos_edge, jump
