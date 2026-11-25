# Copyright (c) 2024, UChicago Argonne, LLC.
# See LICENSE file for details.

import numpy as np
import pytest
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from os.path import join  # noqa: E402
from polartools.diffraction import (  # noqa: E402
    Model,
    fit_peak,
    get_type,
    fit_series,
    load_mesh,
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
# CyberMag spans 49..879 over scan 419, which is what lets the scale= tests
# below clamp at both ends and still leave the open end untouched.
DETECTOR_MIN = 49.0
DETECTOR_MAX = 879.0

# Iaps is identically zero across scan 419, and that is the point: log=True
# has to turn those zeros into 1 before log10 sees them. The unfixed code
# wrote the replacement through `replace(0, 1, inplace=True)`, which pandas 3
# treats as chained assignment and discards, leaving log10 to produce -inf.
ZERO_DETECTOR = "Iaps"

# Scan 12 is a mesh over cryoz/cryox -- the only two-motor plan in this file.
MESH_SCAN = 12


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
        "Scan #",
        "Index",
        "Std Index",
        "Intensity",
        "Std I",
        "Position",
        "Std P",
        "Width",
        "Std W",
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
    datax, _, dataz, det, pos = load_series(
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
    _, _, dataz, _, _ = load_series(
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


def test_load_series_log_replaces_zeros():
    _, _, dataz, _, _ = load_series(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=ZERO_DETECTOR,
        log=True,
        folder=PATH,
    )
    values = np.asarray(dataz[0], dtype=float)
    assert np.isfinite(values).all()
    assert (values == 0.0).all()  # log10(1), not log10(0)


def test_load_series_scale_two_sided():
    _, _, dataz, _, _ = load_series(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=DETECTOR,
        scale=(100, 500),
        folder=PATH,
    )
    values = np.asarray(dataz[0], dtype=float)
    assert values.min() == 100.0
    assert values.max() == 500.0


def test_load_series_scale_one_sided():
    # A single-element scale clamps from above only, so the minimum has to
    # come through untouched.
    _, _, dataz, _, _ = load_series(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=DETECTOR,
        scale=(500,),
        folder=PATH,
    )
    values = np.asarray(dataz[0], dtype=float)
    assert values.max() == 500.0
    assert values.min() == DETECTOR_MIN


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


def test_plot_data_log_replaces_zeros():
    plot_data(
        [SCAN, SCAN, 1],
        source=SPEC_FILE,
        positioner=POSITIONER,
        detector=ZERO_DETECTOR,
        log=True,
        folder=PATH,
    )
    lines = [line for ax in plt.gcf().axes for line in ax.get_lines()]
    assert lines
    for line in lines:
        y = np.asarray(line.get_ydata(), dtype=float)
        assert np.isfinite(y).all()
    plt.close("all")


def test_plot_data_invalid_length():
    with pytest.raises(ValueError):
        plot_data(
            [SCAN, SCAN],  # not a multiple of 3
            source=SPEC_FILE,
            folder=PATH,
        )


# ---------------------------------------------------------------------------
# load_mesh
# ---------------------------------------------------------------------------


def test_load_mesh_spec():
    # Scan 12 was aborted after 23 of its 30x30 points, so this exercises
    # the partial-mesh path: one y row, 23 x columns. zi is (y, x)-ordered.
    scan_range = get_type(MESH_SCAN, source=SPEC_FILE, folder=PATH)
    assert scan_range["plan_name"] == "mesh"
    xi, yi, zi, x_label, y_label, z_label = load_mesh(
        MESH_SCAN,
        scan_range,
        source=SPEC_FILE,
        folder=PATH,
    )
    assert xi.size > 0
    assert yi.size > 0
    assert zi.shape == (yi.size, xi.size)
    assert (x_label, y_label, z_label) == ("Cryo Z", "Cryo X", "IC4")


def test_load_mesh_log_replaces_zeros():
    scan_range = get_type(MESH_SCAN, source=SPEC_FILE, folder=PATH)
    _, _, zi, _, _, _ = load_mesh(
        MESH_SCAN,
        scan_range,
        source=SPEC_FILE,
        log=True,
        detector=ZERO_DETECTOR,
        folder=PATH,
    )
    measured = np.asarray(zi, dtype=float)
    measured = measured[~np.isnan(measured)]  # padding, where mrange trims
    assert measured.size > 0
    assert np.isfinite(measured).all()
