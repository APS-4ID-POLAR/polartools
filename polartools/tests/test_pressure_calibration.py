# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools.pressure_calibration import (
    xrd_calibrate_pressure,
    calculate_pressure,
)
from numpy import allclose
from os.path import join
from spec2nexus.spec import SpecDataFile


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
