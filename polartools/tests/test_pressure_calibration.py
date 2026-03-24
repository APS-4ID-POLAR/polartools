# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools.pressure_calibration import (
    xrd_calibrate_pressure,
    calculate_pressure,
    load_au_params,
    load_ag_params,
    load_pt_params,
    fit_bragg_peak,
    calculate_tth,
)
import numpy as np
from numpy import allclose, linspace
from os.path import join
from spec2nexus.spec import SpecDataFile
from pytest import raises


# TODO: These tests are not great as they rely on fitting the data...
def test_xrd_calibrate_pressure_Au():
    path = join("polartools", "tests", "data_for_test")
    pressure = xrd_calibrate_pressure(
        419,
        "pressure_calibration.dat",
        temperature=9,
        energy=12.0,
        positioner="4C Two Theta",
        detector="CyberMag",
        monitor="IC3",
        folder=path,
    )

    assert allclose(pressure, 0.743387, atol=1e-5)


def test_xrd_calibrate_pressure_Ag():
    path = join("polartools", "tests", "data_for_test")
    spec_file = SpecDataFile(join(path, "pressure_calibration.dat"))
    pressure = xrd_calibrate_pressure(
        419,
        spec_file,
        temperature=9,
        energy=12.0,
        calibrant="Ag",
        positioner="4C Two Theta",
        detector="CyberMag",
        monitor="IC3",
    )

    assert allclose(pressure, 0.81876923, atol=1e-5)


def test_xrd_calculate_pressure_Pt():
    pressure = calculate_pressure(11, 300, 30, [1, 1, 1], "Pt")
    assert allclose(pressure, 60.2033739, atol=1e-5)


def test_wrong_calibrant():
    calibrant = "bla"
    with raises(ValueError) as excinfo:
        _ = calculate_pressure(11, 300, 30, [1, 1, 1], calibrant)
    assert str(excinfo.value) == (
        f'Calibrant must be "Au", "Ag" or "Pt", but {calibrant} was entered.'
    )


def test_load_au_params():
    v0, k0, kp0 = load_au_params(300)
    # At T=300 K, values from Holzapfel table: v0=16.959, k0=166.7, kp0=6.20
    assert allclose(v0, 16.959, atol=1e-3)
    assert allclose(k0, 166.7, atol=0.1)
    assert allclose(kp0, 6.20, atol=1e-3)
    # T=0 K should return the low-temp values
    v0_0, k0_0, _ = load_au_params(0)
    assert allclose(v0_0, 16.7906, atol=1e-3)


def test_load_au_params_out_of_range():
    with raises(ValueError):
        load_au_params(600)


def test_load_ag_params():
    v0, k0, kp0 = load_ag_params(300)
    # At T=300 K: v0=17.057, k0=101.0, kp0=6.15
    assert allclose(v0, 17.057, atol=1e-3)
    assert allclose(k0, 101.0, atol=0.1)
    assert allclose(kp0, 6.15, atol=1e-3)


def test_load_ag_params_out_of_range():
    with raises(ValueError):
        load_ag_params(600)


def test_load_pt_params():
    v0, bt, btp, alphat = load_pt_params()
    assert v0 > 0
    assert bt > 0
    assert btp > 0
    assert alphat > 0


def test_fit_bragg_peak():
    x = linspace(10, 20, 200)
    # Synthetic PseudoVoigt-like peak at center=15
    y = 500.0 * np.exp(-0.5 * ((x - 15.0) / 0.3) ** 2) + 5.0
    result = fit_bragg_peak(x, y)
    assert result.success
    assert allclose(result.best_values["center"], 15.0, atol=0.1)


def test_fit_bragg_peak_with_guesses():
    x = linspace(10, 20, 200)
    y = 500.0 * np.exp(-0.5 * ((x - 15.0) / 0.3) ** 2) + 5.0
    result = fit_bragg_peak(x, y, center=15.0, sigma=0.3, amplitude=500.0)
    assert result.success
    assert allclose(result.best_values["center"], 15.0, atol=0.1)


def test_calculate_tth_round_trip():
    # Round-trip: tth -> pressure -> tth should recover original tth
    tth_original = 26.0  # degrees; above ambient for Au(111) at 12 keV -> positive P
    pressure = calculate_pressure(tth_original, 300, 12.0, [1, 1, 1], "Au")
    assert pressure > 0
    tth_back = calculate_tth(pressure, 300, 12.0, [1, 1, 1], "Au")
    assert allclose(tth_back, tth_original, atol=0.05)


def test_calculate_tth_ag():
    tth_original = 26.0
    pressure = calculate_pressure(tth_original, 300, 12.0, [1, 1, 1], "Ag")
    assert pressure > 0
    tth_back = calculate_tth(pressure, 300, 12.0, [1, 1, 1], "Ag")
    assert allclose(tth_back, tth_original, atol=0.05)
