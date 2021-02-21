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
        scans, 'absorption.dat', detector='IC5', monitor='IC4', folder=path
        )
    result = absorption.normalize_absorption(
        energy*1000., xas, pre_range=[-30, -20],
        post_range=[25, None]
        )

    assert allclose(result['e0'], 7244.49)
    assert allclose(result['flat'].mean(), 0.9401761497318735)
    assert allclose(result['edge_step'], 0.09597937635663531)


def test_fluo():
    path = join('polartools', 'tests', 'data_for_test')
    scans = [199, 200]
    energy, xas, _, _, _ = absorption.load_multi_dichro(
       scans, 'fluorescence.dat', detector='xspRoi1totc', monitor='IC3',
       folder=path, transmission=False
       )

    result = absorption.normalize_absorption(
        energy*1000., xas, pre_range=[-30, -20], pre_order=0,
        post_range=[15, None], post_order=1
        )

    norm_corr = absorption.fluo_corr(
        result['norm'], 'EuFe2As2', 'Eu', 'L3', 'La', 10, 80
        )

    assert allclose(norm_corr.mean(), 0.7611857565650514)
    assert allclose(norm_corr.max(), 2.6009413580639853)
