# Copyright (c) 2024, UChicago Argonne, LLC.
# See LICENSE file for details.

import numpy as np
from numpy import allclose
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from polartools._pyrixs import (  # noqa: E402
    poly,
    image_to_photon_events,
    bin_edges_centers,
    get_curvature_offsets,
    fit_poly,
    fit_curvature,
    plot_curvature,
    extract,
)


def test_poly():
    assert allclose(poly(0, 1, 2, 3), 3.0)  # 0 + 0 + 3
    assert allclose(poly(1, 1, 2, 3), 6.0)  # 1 + 2 + 3
    assert allclose(poly(2, 1, 2, 3), 11.0)  # 4 + 4 + 3
    x = np.array([0.0, 1.0, 2.0])
    assert allclose(poly(x, 1, 2, 3), [3, 6, 11])


def test_image_to_photon_events():
    image = np.array([[1.0, 2.0], [3.0, 4.0]])
    events = image_to_photon_events(image)
    assert events.shape == (4, 3)
    # x centers at 0.5 and 1.5, y centers at 0.5 and 1.5
    assert allclose(sorted(events[:, 2]), [1.0, 2.0, 3.0, 4.0])
    # All x values should be 0.5 or 1.5
    assert set(events[:, 0]).issubset({0.5, 1.5})


def test_bin_edges_centers():
    edges, centers = bin_edges_centers(0, 10, 2)
    assert len(edges) == len(centers) + 1
    assert edges[0] >= 0
    assert edges[-1] <= 10
    assert allclose(np.diff(edges), 2.0)
    assert allclose(np.diff(centers), 2.0)


def test_get_curvature_offsets():
    # Horizontal stripe -> offsets should all be near 0 (after referencing
    # center)
    rng = np.random.default_rng(42)
    n = 10000
    x = rng.uniform(0, 512, n)
    y = rng.uniform(200, 202, n)  # narrow horizontal band
    intensity = np.ones(n)
    photon_events = np.column_stack([x, y, intensity])
    x_centers, offsets = get_curvature_offsets(photon_events, binx=64, biny=0.5)
    assert len(x_centers) == len(offsets)
    # Center offset is 0 by construction; others should be small
    mid = len(offsets) // 2
    assert abs(offsets[mid]) < 1.0


def test_fit_poly():
    x = np.linspace(0, 10, 50)
    y = 0.05 * x**2 - 0.3 * x + 1.0  # noise-free polynomial
    result = fit_poly(x, y)
    assert result.success
    assert allclose(result.best_values["p2"], 0.05, atol=1e-5)
    assert allclose(result.best_values["p1"], -0.3, atol=1e-5)
    assert allclose(result.best_values["p0"], 1.0, atol=1e-5)


def test_fit_curvature():
    rng = np.random.default_rng(42)
    n = 10000
    x = rng.uniform(0, 512, n)
    y = rng.uniform(200, 202, n)
    intensity = np.ones(n)
    photon_events = np.column_stack([x, y, intensity])
    curvature = fit_curvature(
        photon_events, binx=64, biny=0.5, CONSTANT_OFFSET=201
    )
    assert len(curvature) == 3
    assert abs(curvature[0]) < 0.01  # near-zero curvature for flat stripe
    assert curvature[2] == 201  # CONSTANT_OFFSET passthrough


def test_plot_curvature():
    rng = np.random.default_rng(42)
    n = 500
    x = rng.uniform(0, 256, n)
    y = rng.uniform(100, 200, n)
    intensity = np.ones(n)
    photon_events = np.column_stack([x, y, intensity])
    curvature = np.array([0.0, 0.0, 150.0])
    _, ax = plt.subplots()
    artists = plot_curvature(ax, curvature, photon_events)
    assert len(artists) == 1
    plt.close("all")


def test_extract():
    rng = np.random.default_rng(42)
    n = 5000
    x = rng.uniform(0, 100, n)
    y = rng.uniform(50, 150, n)
    intensity = np.ones(n)
    photon_events = np.column_stack([x, y, intensity])
    curvature = np.array([0.0, 0.0, 0.0])
    spectrum = extract(photon_events, curvature, biny=1)
    assert spectrum.shape[1] == 2  # (pixel, intensity) columns
    assert spectrum[:, 1].sum() > 0  # non-zero intensity
