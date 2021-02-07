# Copyright (c) 2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools import absorption
from numpy import allclose
from os.path import join


def test_load_multi_xas():
    path = join('polartools', 'tests', 'data_for_test')
    scans = [28, 29, 30, 31, 32]
    _, xas, _ = absorption.load_multi_xas(
        scans, 'absorption.dat', detector='IC4', monitor='IC3', folder=path)
    assert allclose(xas.mean(), -7.64847338574727)


def test_load_multi_dichro():
    path = join('polartools', 'tests', 'data_for_test')
    scans = [39, 40, 41]
    _, _, xmcd, _, _ = absorption.load_multi_dichro(
       scans, 'absorption.dat', detector='IC4', monitor='IC3', folder=path)
    assert allclose(xmcd.mean(), -0.09862229544158153)


def test_load_multi_lockin():
    path = join('polartools', 'tests', 'data_for_test')
    scans = [102, 103, 104]
    _, _, xmcd, _, _ = absorption.load_multi_lockin(scans, 'absorption.dat',
                                                    folder=path)
    assert allclose(xmcd.mean(), -2935.7871485943774)
