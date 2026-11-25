# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

import numpy as np
import pytest

from polartools.spectrum_fit import (
    PeakSpec,
    ShapeFit,
    area_weights,
    decompose,
    fit_shape,
    suggest_peaks,
    unit_gaussian,
)

BINS = 400
X = np.arange(BINS, dtype=float)
# Two overlapping lines on a sloped background, as on a softglue scan's
# element 1: a strong line with a weaker one on its tail.
TRUE = [(150.0, 9.0, 4.0e5), (178.0, 9.0, 8.0e4)]


def spectrum(areas=None, slope=40.0, offset=300.0):
    areas = [a for *_, a in TRUE] if areas is None else areas
    y = offset + slope * (X - 150.0) / 100.0
    for (centre, sigma, _), area in zip(TRUE, areas):
        y = y + area * unit_gaussian(X, centre, sigma)
    return y


def starting_peaks():
    return [PeakSpec(147.0, 7.0), PeakSpec(182.0, 7.0)]


def test_fit_recovers_overlapping_peaks():
    rng = np.random.default_rng(1)
    y = rng.poisson(spectrum()).astype(float)
    fit = fit_shape(y, 90, 240, starting_peaks(), "linear")
    assert fit.success, fit.message
    np.testing.assert_allclose(fit.centres, [150.0, 178.0], atol=0.3)
    np.testing.assert_allclose(fit.sigmas, [9.0, 9.0], atol=0.3)
    np.testing.assert_allclose(fit.areas, [4.0e5, 8.0e4], rtol=0.02)
    assert fit.redchi < 1.5


def test_area_weights_give_each_peak_area():
    shape = ShapeFit(
        lo=90,
        hi=240,
        background="linear",
        centres=np.array([150.0, 178.0]),
        sigmas=np.array([9.0, 9.0]),
    )
    window = spectrum()[90:241]
    # Exact on a noiseless spectrum the model describes: the background is
    # solved away, and the second peak's tail is not counted in the first.
    assert area_weights(shape, [0]) @ window == pytest.approx(4.0e5, rel=1e-6)
    assert area_weights(shape, [1]) @ window == pytest.approx(8.0e4, rel=1e-6)
    assert area_weights(shape, [0, 1]) @ window == pytest.approx(
        4.8e5, rel=1e-6
    )
    assert not area_weights(shape, []).any()


def test_area_is_additive_over_frames():
    # The property the viewer relies on: the area of a summed spectrum is
    # the sum of the areas, so per-frame areas can be totalled like counts.
    rng = np.random.default_rng(2)
    frames = rng.poisson(spectrum() / 20.0, size=(20, BINS)).astype(float)
    fit = fit_shape(frames.sum(axis=0), 90, 240, starting_peaks(), "linear")
    weights = area_weights(fit, [0])
    per_frame = frames[:, 90:241] @ weights
    assert per_frame.sum() == pytest.approx(
        frames.sum(axis=0)[90:241] @ weights
    )


def test_background_none_counts_the_background():
    shape = ShapeFit(
        lo=90,
        hi=240,
        background="none",
        centres=np.array([150.0, 178.0]),
        sigmas=np.array([9.0, 9.0]),
    )
    with_bg = area_weights(shape, [0]) @ spectrum()[90:241]
    without = area_weights(shape, [0]) @ spectrum(offset=0.0, slope=0.0)[90:241]
    # The flat 300 and the slope under the peak are soaked up into it when
    # nothing models them: about 2% here.
    assert with_bg > without * 1.01


def test_shared_width_ties_the_sigmas():
    y = spectrum()
    fit = fit_shape(y, 90, 240, starting_peaks(), "linear", shared_width=True)
    assert fit.success, fit.message
    assert fit.sigmas[1] == pytest.approx(fit.sigmas[0])


def test_decompose_matches_area_weights():
    y = spectrum()
    fit = fit_shape(y, 90, 240, starting_peaks(), "constant")
    coefficients, matrix = decompose(y, fit)
    assert matrix.shape == (151, 3)
    assert coefficients[0] == pytest.approx(area_weights(fit, [0]) @ y[90:241])


def test_fit_reports_failure_rather_than_raising():
    assert not fit_shape(np.zeros(BINS), 90, 240, starting_peaks()).success
    assert not fit_shape(spectrum(), 90, 240, []).success
    assert not fit_shape(spectrum(), 100, 104, starting_peaks()).success


def test_suggest_peaks_finds_both_lines():
    found = suggest_peaks(spectrum(), 90, 240)
    centres = sorted(p.centre for p in found)
    assert len(centres) == 2
    np.testing.assert_allclose(centres, [150.0, 178.0], atol=2.0)


def test_centres_stay_near_their_start():
    # A spurious third peak started between the two lines must not walk
    # under either of them and split its area.
    peaks = starting_peaks() + [PeakSpec(164.0, 3.0)]
    fit = fit_shape(spectrum(), 90, 240, peaks, "linear")
    assert fit.success, fit.message
    assert 164.0 - 9.0 <= fit.centres[2] <= 164.0 + 9.0
