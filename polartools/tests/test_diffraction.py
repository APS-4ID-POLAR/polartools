# Copyright (c) 2024, UChicago Argonne, LLC.
# See LICENSE file for details.

import numpy as np
from numpy import allclose
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from os.path import join
from polartools.diffraction import (
    Model,
    fit_peak,
    get_type,
    fit_series,
    load_series,
    plot_fit,
    plot_data,
)

PATH = join("polartools", "tests", "data_for_test")
SPEC_FILE = "pressure_calibration.dat"

# Scan 419 in pressure_calibration.dat is an ascan on "4C Two Theta"
SCAN = 419
POSITIONER = "4C Two Theta"
DETECTOR = "CyberMag"
MONITOR = "IC3"


def _make_gaussian(center=0.0, sigma=1.0, amplitude=100.0, n=100):
    rng = np.random.default_rng(42)
    x = np.linspace(center - 5 * sigma, center + 5 * sigma, n)
    y = amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)
    y += rng.normal(0, amplitude * 0.01, n)  # 1% noise
    return x, y


# ---------------------------------------------------------------------------
# fit_peak
# ---------------------------------------------------------------------------

def test_fit_peak_gaussian():
    x, y = _make_gaussian()
    result = fit_peak(x, y, model=Model.Gaussian)
    assert result.success
    assert abs(result.best_values["center"]) < 0.3


def test_fit_peak_lorentzian():
    x, y = _make_gaussian()
    result = fit_peak(x, y, model=Model.Lorentzian)
    assert result.success
    assert abs(result.best_values["center"]) < 0.3


def test_fit_peak_pseudovoigt():
    x, y = _make_gaussian()
    result = fit_peak(x, y, model=Model.PseudoVoigt)
    assert result.success
    assert abs(result.best_values["center"]) < 0.3


# ---------------------------------------------------------------------------
# get_type
# ---------------------------------------------------------------------------

def test_get_type_spec():
    info = get_type(SCAN, source=SPEC_FILE, folder=PATH)
    assert isinstance(info, dict)
    assert "plan_name" in info
    assert "scan_no" in info
    assert "motor0" in info


# ---------------------------------------------------------------------------
# fit_series
# ---------------------------------------------------------------------------

def test_fit_series_spec_single_scan():
    result = fit_series(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=DETECTOR,
        folder=PATH,
    )
    assert result.shape == (1, 9)
    assert list(result.columns) == [
        "Scan #", "Index", "Std Index",
        "Intensity", "Std I",
        "Position", "Std P",
        "Width", "Std W",
    ]
    assert result["Position"].iloc[0] > 0  # 2-theta must be positive


def test_fit_series_spec_normalized():
    result = fit_series(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=DETECTOR,
        monitor=MONITOR,
        normalize=True,
        folder=PATH,
    )
    assert result.shape == (1, 9)
    assert result["Intensity"].iloc[0] > 0


def test_fit_series_invalid_length():
    with pytest.raises(ValueError, match="expected 3"):
        fit_series(
            [SCAN, SCAN],  # not a multiple of 3
            source=SPEC_FILE,
            folder=PATH,
        )


# ---------------------------------------------------------------------------
# load_series
# ---------------------------------------------------------------------------

def test_load_series_spec():
    datax, datay, dataz, det, pos = load_series(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=DETECTOR,
        folder=PATH,
    )
    assert len(datax) == 1
    assert len(dataz) == 1
    assert len(dataz[0]) > 0
    assert det == DETECTOR
    assert pos == POSITIONER


def test_load_series_normalize():
    datax, datay, dataz, det, pos = load_series(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=DETECTOR,
        monitor=MONITOR,
        normalize=True,
        folder=PATH,
    )
    assert len(dataz) == 1
    assert dataz[0].max() < dataz[0].max() + 1  # finite values


def test_load_series_invalid_length():
    # load_series checks length after the initial nbp loop, so use length 4
    # (which completes the loop before hitting the check), not length 2
    with pytest.raises(ValueError, match="expected 3"):
        load_series(
            [SCAN, SCAN, SCAN, SCAN],
            source=SPEC_FILE,
            positioner=POSITIONER,
            detector=DETECTOR,
            folder=PATH,
        )


# ---------------------------------------------------------------------------
# plot_fit (wraps fit_series + matplotlib)
# ---------------------------------------------------------------------------

def test_plot_fit_spec():
    result = plot_fit(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=DETECTOR,
        folder=PATH,
    )
    assert result.shape == (1, 9)
    plt.close("all")


# ---------------------------------------------------------------------------
# plot_data
# ---------------------------------------------------------------------------

def test_plot_data_spec():
    # plot_data returns None; just verify it runs without error
    plot_data(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=DETECTOR,
        folder=PATH,
    )
    plt.close("all")


def test_plot_data_invalid_length():
    with pytest.raises(ValueError):
        plot_data(
            [SCAN, SCAN],  # not a multiple of 3
            source=SPEC_FILE,
            folder=PATH,
        )
