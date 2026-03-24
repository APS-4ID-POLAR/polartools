# Copyright (c) 2024, UChicago Argonne, LLC.
# See LICENSE file for details.

import numpy as np
from polartools._larch import finde0, index_nearest


def test_index_nearest():
    array = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert index_nearest(array, 3.1) == 2
    assert index_nearest(array, 0.5) == 0
    assert index_nearest(array, 6.0) == 4
    assert index_nearest(array, 2.6) == 2  # equidistant, argmin picks lower


def test_finde0():
    # tanh edge: derivative peaks at e0
    energy = np.linspace(7000, 7100, 300)
    e0 = 7050.0
    mu = np.tanh((energy - e0) / 2.0)
    result = finde0(energy, mu)
    assert abs(result - e0) < 5.0


def test_finde0_2d_input():
    # finde0 should squeeze 2D inputs
    energy_1d = np.linspace(7000, 7100, 300)
    e0 = 7050.0
    mu_1d = np.tanh((energy_1d - e0) / 2.0)
    energy_2d = energy_1d.reshape(1, -1)
    mu_2d = mu_1d.reshape(-1, 1)
    result = finde0(energy_2d, mu_2d)
    assert abs(result - e0) < 5.0
