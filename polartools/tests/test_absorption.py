# Copyright (c) 2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools import absorption
from numpy import allclose
from os.path import join


def test_load_multi_xas():
    path = join("polartools", "tests", "data_for_test")
    scans = [28, 29, 30, 31, 32]
    _, xas, _ = absorption.load_multi_xas(
        scans, "absorption.dat", detector="IC4", monitor="IC3", folder=path
    )
    assert allclose(xas.mean(), -7.64847338574727)


def test_load_multi_dichro():
    path = join("polartools", "tests", "data_for_test")
    scans = [39, 40, 41]
    _, _, xmcd, _, _ = absorption.load_multi_dichro(
        scans, "absorption.dat", detector="IC4", monitor="IC3", folder=path
    )
    assert allclose(xmcd.mean(), -0.09862229544158153)


def test_load_multi_lockin():
    path = join("polartools", "tests", "data_for_test")
    scans = [102, 103, 104]
    _, _, xmcd, _, _ = absorption.load_multi_lockin(
        scans, "absorption.dat", folder=path
    )
    assert allclose(xmcd.mean(), -2935.7871485943774)


def test_load_absorption():
    # Meant to test loading bluesky vs. SPEC generated spec_file.
    path = join("polartools", "tests", "data_for_test")

    _, xas = absorption.load_absorption(
        44, "bluesky_spec.dat", positioner="fourc_h", folder=path
    )
    assert allclose(xas.mean(), 0.7141861035115614)

    _, xas = absorption.load_absorption(18, "absorption.dat", folder=path)
    assert allclose(xas.mean(), 0.7581461103557833)


def test_normalization():
    path = join("polartools", "tests", "data_for_test")
    scans = [28, 29, 30, 31, 32]
    energy, xas, _ = absorption.load_multi_xas(
        scans, "absorption.dat", detector="IC5", monitor="IC4", folder=path
    )
    result = absorption.normalize_absorption(
        energy * 1000.0, xas, pre_range=[-30, -20], post_range=[25, None]
    )

    assert allclose(result["e0"], 7244.49)
    assert allclose(result["flat"].mean(), 0.9401761497318735)
    assert allclose(result["edge_step"], 0.09597937635663531)


def test_fluo():
    path = join("polartools", "tests", "data_for_test")
    scans = [199, 200]
    energy, xas, _, _, _ = absorption.load_multi_dichro(
        scans,
        "fluorescence.dat",
        detector="xspRoi1totc",
        monitor="IC3",
        folder=path,
        transmission=False,
    )

    result = absorption.normalize_absorption(
        energy * 1000.0,
        xas,
        pre_range=[-30, -20],
        pre_order=0,
        post_range=[15, None],
        post_order=1,
    )

    norm_corr = absorption.fluo_corr(
        result["norm"], "EuFe2As2", "Eu", "L3", "La", 10, 80
    )

    assert allclose(norm_corr.mean(), 0.7611857565650514)
    assert allclose(norm_corr.max(), 2.6009413580639853)


def test_process_xmcd():
    path = join("polartools", "tests", "data_for_test")
    scans = [39, 40, 41]
    load_parameters = dict(detector="IC4", monitor="IC3", folder=path)
    normalization_parameters = dict(
        pre_range=[-30, -10],
        pre_order=0,
        post_range=[10, 30],
        post_order=0,
    )
    plus, minus = absorption.process_xmcd(
        scans,
        scans,
        "absorption.dat",
        xmcd_kind="dichro",
        load_parameters=load_parameters,
        normalization_parameters=normalization_parameters,
    )

    assert allclose(plus["xmcd"].mean(), 3.6117537051437343)
    assert list(minus.keys()) == [
        "energy",
        "mu",
        "preedge",
        "pre1",
        "pre2",
        "pre_order",
        "nvict",
        "e0",
        "precoefs",
        "energy_unit",
        "postedge",
        "postcoefs",
        "post1",
        "post2",
        "post_order",
        "norm",
        "edge_step",
        "flat",
        "flat1",
        "flat2",
        "flat_order",
        "flatcoefs",
        "flat_function",
        "xmcd",
    ]


def test_plot_xmcd():
    path = join("polartools", "tests", "data_for_test")
    scans = [39, 40, 41]
    load_parameters = dict(detector="IC4", monitor="IC3", folder=path)
    normalization_parameters = dict(
        pre_range=[-30, -10],
        pre_order=0,
        post_range=[10, 30],
        post_order=0,
    )
    plus, minus = absorption.process_xmcd(
        scans,
        scans,
        "absorption.dat",
        xmcd_kind="dichro",
        load_parameters=load_parameters,
        normalization_parameters=normalization_parameters,
    )

    _, axs = absorption.plot_xmcd(plus, minus)
    assert len(axs) == 6


def test_save_xmcd():
    path = join("polartools", "tests", "data_for_test")
    scans = [39, 40, 41]
    load_parameters = dict(detector="IC4", monitor="IC3", folder=path)
    normalization_parameters = dict(
        pre_range=[-30, -10],
        pre_order=0,
        post_range=[10, 30],
        post_order=0,
    )
    plus, minus = absorption.process_xmcd(
        scans,
        scans,
        "absorption.dat",
        xmcd_kind="dichro",
        load_parameters=load_parameters,
        normalization_parameters=normalization_parameters,
    )
    absorption.save_xmcd(plus, minus, "xmcd_save_test.dat")
