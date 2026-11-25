# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

import h5py
import numpy as np
import pytest

pytest.importorskip("PyQt6", reason="PyQt6 not installed — skipping GUI tests")

# Three separated lines on a flat background, at these bins.
CENTRES = (700, 740, 780)
SIGMA = 8.0


@pytest.fixture(scope="session")
def qapp():
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def fitted(qapp, tmp_path):
    """A viewer on a synthetic 4-element cube, fitted with three peaks."""
    from polartools import spectrum_viewer as sv

    bins = np.arange(4096)
    spectrum = 20.0 + sum(
        height * np.exp(-0.5 * ((bins - centre) / SIGMA) ** 2)
        for centre, height in zip(CENTRES, (2000.0, 800.0, 1200.0))
    )
    rng = np.random.default_rng(0)
    path = tmp_path / "scan_0001.h5"
    with h5py.File(path, "w") as f:
        f[sv.SPECTRUM_DATA_PATH] = rng.poisson(
            np.tile(spectrum, (20, 4, 1))
        ).astype("uint16")

    window = sv.SpectrumViewerWindow([str(path)])
    window._ensure_fit_defaults()
    kev = window._bin_to_kev
    window.fit_lo_kev, window.fit_hi_kev = float(kev(640)), float(kev(840))
    fwhm = SIGMA * sv.FWHM_PER_SIGMA * window.gain_spin.value()
    window.fit_peaks = [window._new_peak(float(kev(c)), fwhm) for c in CENTRES]
    window._open_fit_window()
    window._fit_all()
    yield window
    window.fit_window.close()
    window.close()


def test_removed_peak_does_not_shift_use(fitted):
    # Untick the second peak, then remove the first without refitting: the
    # trend must sum the old third peak, not the old second one.
    assert fitted._fit_used() == [0, 1, 2]
    fitted.fit_peaks[1]["use"] = False
    del fitted.fit_peaks[0]
    fitted.fit_window._changed()

    assert fitted._fitted_rows() == [None, 0, 1]
    assert fitted._fit_used_pairs() == [(2, 1)]
    assert fitted._fit_used() == [2]
    assert fitted.fit_dirty
    assert "peak 2 (7.80 keV)" in fitted._fit_description(short=True)
    # The table's row 1 (the old third peak) shows that peak's fitted centre.
    centre = float(fitted.fit_window.table.item(1, 4).text())
    assert centre == pytest.approx(float(fitted._bin_to_kev(780)), abs=1e-3)


def test_added_peak_is_not_summed_until_refit(fitted):
    fitted.fit_peaks.append(fitted._new_peak(7.6, 150.0))
    fitted.fit_window._changed()
    assert fitted._fit_used() == [0, 1, 2]
    assert fitted.fit_window.table.item(3, 4).text() == ""


def test_energy_axis_change_marks_fit_stale(fitted):
    assert not fitted.fit_dirty
    fitted.gain_spin.setValue(fitted.gain_spin.value() + 1.0)
    assert fitted.fit_dirty
    # The header reports the fit's own bins, not the window's new ones.
    assert "bins 640-840" in fitted._fit_description()


def test_unfit_element_is_left_out_of_sum(fitted):
    from polartools import spectrum_viewer as sv

    fitted.fit_results.pop((0, 3))
    fitted.signal_combo.setCurrentIndex(1)
    settings = fitted._trend_settings()
    _x, curves, _label, notes, _used = fitted._trend_curves(
        [0, 1, 2, 3], True, "all", False, settings
    )
    assert any("element(s) 4, left out" in note for note in notes)
    assert np.isfinite(curves[sv.SUM_KEY]).all()
