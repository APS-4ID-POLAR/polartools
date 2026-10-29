# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

import pytest
import numpy as np

pytest.importorskip("PyQt6", reason="PyQt6 not installed — skipping GUI tests")

from polartools._scan_viewer_base import merge_grids  # noqa: E402


def grid(x, y, z):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    return {
        "x_name": "k",
        "y_name": "h",
        "x": x,
        "y": y,
        "z": np.asarray(z, dtype=float),
        "points": int(np.isfinite(z).sum()),
    }


def test_merge_grids_single_is_unchanged():
    one = grid([0, 1], [0, 1], [[1, 2], [3, 4]])
    assert merge_grids([one]) is one


def test_merge_grids_averages_overlapping_row():
    # Rows 0-1 from the first scan, 1-2 from the second; row 1 in both.
    first = grid([0, 1], [0, 1], [[1, 2], [3, 4]])
    second = grid([0, 1], [1, 2], [[5, 6], [7, 8]])
    merged = merge_grids([first, second])
    np.testing.assert_array_equal(merged["y"], [0, 1, 2])
    np.testing.assert_array_equal(merged["z"], [[1, 2], [4, 5], [7, 8]])
    assert merged["averaged"] == 2
    assert merged["scans"] == 2
    assert merged["points"] == 8


def test_merge_grids_keeps_finer_rows():
    coarse = grid([0, 1], [0, 1], [[1, 1], [2, 2]])
    fine = grid([0, 1], [1, 1.5, 2], [[4, 4], [5, 5], [6, 6]])
    merged = merge_grids([coarse, fine])
    np.testing.assert_array_equal(merged["y"], [0, 1, 1.5, 2])
    np.testing.assert_array_equal(merged["z"][:, 0], [1, 3, 5, 6])


def test_merge_grids_joins_levels_a_rounding_apart():
    first = grid([0, 0.5], [0, 0.5], [[1, 1], [1, 1]])
    second = grid([1e-12, 0.5], [0.5 + 1e-12, 1], [[3, 3], [3, 3]])
    merged = merge_grids([first, second])
    assert len(merged["x"]) == 2
    assert len(merged["y"]) == 3
    np.testing.assert_array_equal(merged["z"][1], [2, 2])


def test_merge_grids_ignores_unmeasured_cells():
    first = grid([0, 1], [0, 1], [[1, 2], [3, np.nan]])
    second = grid([0, 1], [1, 2], [[5, 6], [np.nan, np.nan]])
    merged = merge_grids([first, second])
    np.testing.assert_array_equal(merged["z"][1], [4, 6])
    assert merged["averaged"] == 1
    assert np.isnan(merged["z"][2]).all()
