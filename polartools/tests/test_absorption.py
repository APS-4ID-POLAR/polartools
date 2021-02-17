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


def test_load_absorption():
    # Meant to test loading bluesky vs. SPEC generated spec_file.
    path = join('polartools', 'tests', 'data_for_test')

    _, xas = absorption.load_absorption(
        44, 'bluesky_spec.dat', positioner='fourc_h', folder=path)
    assert allclose(xas.mean(), 0.7141861035115614)

    _, xas = absorption.load_absorption(18, 'absorption.dat', folder=path)
    assert allclose(xas.mean(), 0.7581461103557833)


def test_normalization():
    path = join('polartools', 'tests', 'data_for_test')
    scans = [28, 29, 30, 31, 32]
    energy, xas, _ = absorption.load_multi_xas(
        scans, 'absorption.dat', detector='IC4', monitor='IC3', folder=path)

    print(energy)
    result = absorption.normalize_absorption(energy.to_numpy()*1000., xas)

    assert allclose(result['e0'], 7233.5)
    assert allclose(result['flat'].mean(), -0.25136020653859914)
    assert allclose(result['edge_step'], 0.005179323346776066)
