# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

import os
import pytest
import numpy as np
from unittest.mock import patch, MagicMock

# Skip the entire module if PyQt6 or pyqtgraph are not installed.
pytest.importorskip("PyQt6", reason="PyQt6 not installed — skipping GUI tests")
pytest.importorskip(
    "pyqtgraph", reason="pyqtgraph not installed — skipping GUI tests"
)

from PyQt6.QtCore import QSignalBlocker  # noqa: E402
from PyQt6.QtWidgets import QMessageBox  # noqa: E402


@pytest.fixture(scope="session")
def qapp():
    """Session-scoped QApplication — created once, reused across all tests."""
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def win(qapp):
    from polartools.xanes_gui import MainWindow

    window = MainWindow()
    yield window
    window.close()


@pytest.fixture
def synthetic_data():
    """Fake Fe-edge arrays: energy in keV, mu a sigmoid step."""
    energy = np.linspace(7.0, 7.3, 200)  # keV (GUI multiplies by 1000)
    mu = 1 / (1 + np.exp(-(energy * 1000 - 7112) / 3))
    return energy, mu


@pytest.fixture
def loaded_win(win, synthetic_data):
    """Window with synthetic data injected — simulates state after a load."""
    e, mu = synthetic_data
    win._energy = e * 1000
    win._mu = mu
    win._e0_val = 7112.0
    return win


@pytest.fixture
def normalized_win(loaded_win):
    """Window with a results dict set — simulates state after normalization."""
    energy = loaded_win._energy
    stub = {
        "energy": energy,
        "mu": loaded_win._mu,
        "norm": np.ones_like(energy),
        "flat": np.ones_like(energy),
        "preedge": np.zeros_like(energy),
        "postedge": np.ones_like(energy),
        "e0": 7112.0,
        "edge_step": 1.0,
    }
    loaded_win._results = stub.copy()
    return loaded_win


def _norm_stub(energy, e0=7112.0):
    return {
        "energy": energy,
        "mu": np.ones_like(energy),
        "norm": np.ones_like(energy),
        "flat": np.ones_like(energy),
        "preedge": np.zeros_like(energy),
        "postedge": np.ones_like(energy),
        "e0": e0,
        "edge_step": 1.0,
    }


# ─── Window creation ──────────────────────────────────────────────────────────


def test_window_title(win):
    assert win.windowTitle() == "XANES Processor"


def test_initial_state_no_data(win):
    assert win._energy is None
    assert win._results is None


def test_save_button_disabled_on_start(win):
    assert not win.btn_save.isEnabled()


# ─── Source stacked widgets ───────────────────────────────────────────────────


def test_source_stack_follows_combo(win):
    for i in range(win.cb_source.count()):
        win.cb_source.setCurrentIndex(i)
        assert win._source_stack.currentIndex() == i


# ─── _parse_scan_list ─────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "text, expected",
    [
        ("1, 2, 3", [1, 2, 3]),
        ("4;5;6", [4, 5, 6]),
        ("10", [10]),
        (" 7 , 8 ", [7, 8]),
    ],
)
def test_parse_scan_list_valid(win, text, expected):
    assert win._parse_scan_list(text) == expected


def test_parse_scan_list_empty_raises(win):
    with pytest.raises(ValueError):
        win._parse_scan_list("")


def test_parse_scan_list_whitespace_raises(win):
    with pytest.raises(ValueError):
        win._parse_scan_list("   ")


def test_parse_scan_list_string_fallback(win):
    assert win._parse_scan_list("myscan") == ["myscan"]


# ─── _build_norm_kwargs ───────────────────────────────────────────────────────


def test_norm_kwargs_defaults(win):
    kw = win._build_norm_kwargs()
    assert kw["e0"] is None
    assert kw["edge_step"] is None
    assert kw["pre_range"] is None
    assert kw["post_range"] is None
    assert kw["pre_order"] == 1
    assert kw["nvict"] == 0
    assert kw["post_order"] is None
    assert kw["flat_range"] is None
    assert kw["flat_order"] is None


def test_norm_kwargs_manual_e0(win):
    win.chk_e0_auto.setChecked(False)
    win.le_e0.setText("7112.5")
    kw = win._build_norm_kwargs()
    assert kw["e0"] == pytest.approx(7112.5)


def test_norm_kwargs_manual_edge_step(win):
    win.chk_es_auto.setChecked(False)
    win.le_edge_step.setText("0.85")
    kw = win._build_norm_kwargs()
    assert kw["edge_step"] == pytest.approx(0.85)


def test_norm_kwargs_pre_range(win):
    win.le_pre1.setText("-150")
    win.le_pre2.setText("-30")
    kw = win._build_norm_kwargs()
    assert kw["pre_range"] == [-150.0, -30.0]


def test_norm_kwargs_post_range(win):
    win.le_post1.setText("50")
    win.le_post2.setText("300")
    kw = win._build_norm_kwargs()
    assert kw["post_range"] == [50.0, 300.0]


def test_norm_kwargs_flat_range(win):
    win.le_flat1.setText("20")
    win.le_flat2.setText("200")
    kw = win._build_norm_kwargs()
    assert kw["flat_range"] == [20.0, 200.0]


def test_norm_kwargs_post_order(win):
    win.cb_post_order.setCurrentText("2")
    kw = win._build_norm_kwargs()
    assert kw["post_order"] == 2


# ─── _resolve_load_kwargs ─────────────────────────────────────────────────────


def test_resolve_load_kwargs_defaults(win):
    kw = win._resolve_load_kwargs({})
    assert kw["transmission"] is True
    assert "positioner" not in kw
    assert "detector" not in kw
    assert "monitor" not in kw


def test_resolve_load_kwargs_custom(win):
    win.le_positioner.setText("4idenergy")
    win.le_detector.setText("IC5")
    win.chk_transmission.setChecked(False)
    kw = win._resolve_load_kwargs({})
    assert kw["positioner"] == "4idenergy"
    assert kw["detector"] == "IC5"
    assert kw["transmission"] is False


def test_resolve_load_kwargs_monitor(win):
    win.le_monitor.setText("IC4")
    kw = win._resolve_load_kwargs({})
    assert kw["monitor"] == "IC4"


def test_resolve_load_kwargs_passes_extra(win):
    kw = win._resolve_load_kwargs({"folder": "/data/folder"})
    assert kw["folder"] == "/data/folder"


# ─── e0 auto / manual toggle ──────────────────────────────────────────────────


def test_e0_entry_disabled_when_auto(win):
    win.chk_e0_auto.setChecked(True)
    assert not win.le_e0.isEnabled()


def test_e0_entry_enabled_when_manual(win):
    win.chk_e0_auto.setChecked(False)
    assert win.le_e0.isEnabled()


def test_edge_step_entry_disabled_when_auto(win):
    win.chk_es_auto.setChecked(True)
    assert not win.le_edge_step.isEnabled()


# ─── _resolve_source ──────────────────────────────────────────────────────────


def test_resolve_source_spec_path_only(win):
    win.cb_source.setCurrentIndex(0)  # SPEC
    win.le_spec_path.setText("/data/scan.dat")
    win.le_spec_folder.setText("")
    src, kw = win._resolve_source()
    assert src == "/data/scan.dat"
    assert kw == {}


def test_resolve_source_spec_with_folder(win):
    win.cb_source.setCurrentIndex(0)
    win.le_spec_path.setText("/data/scan.dat")
    win.le_spec_folder.setText("/data/subdir")
    src, kw = win._resolve_source()
    assert src == "/data/scan.dat"
    assert kw["folder"] == "/data/subdir"


def test_resolve_source_spec_empty_path_raises(win):
    win.cb_source.setCurrentIndex(0)
    win.le_spec_path.setText("")
    with pytest.raises(ValueError, match="SPEC file path"):
        win._resolve_source()


def test_resolve_source_hdf5_defaults(win):
    win.cb_source.setCurrentIndex(1)  # HDF5
    win.le_hdf_folder.setText("/data/hdf5/")
    src, kw = win._resolve_source()
    assert src == "hdf5"
    assert kw["folder"] == "/data/hdf5/"
    assert "fname_format" in kw
    assert "h5_location" in kw


def test_resolve_source_hdf5_custom(win):
    win.cb_source.setCurrentIndex(1)
    win.le_hdf_folder.setText("/data/hdf5/")
    win.le_hdf_format.setText("custom_{:05d}.h5")
    win.le_hdf_loc.setText("entry/data")
    src, kw = win._resolve_source()
    assert src == "hdf5"
    assert kw["fname_format"] == "custom_{:05d}.h5"
    assert kw["h5_location"] == "entry/data"


def test_resolve_source_hdf5_empty_folder_raises(win):
    win.cb_source.setCurrentIndex(1)
    win.le_hdf_folder.setText("")
    with pytest.raises(ValueError, match="HDF5 folder"):
        win._resolve_source()


def test_resolve_source_csv(win):
    win.cb_source.setCurrentIndex(2)  # CSV
    win.le_csv_folder.setText("/data/csv/")
    src, kw = win._resolve_source()
    assert src == "csv"
    assert kw["folder"] == "/data/csv/"


def test_resolve_source_csv_empty_raises(win):
    win.cb_source.setCurrentIndex(2)
    win.le_csv_folder.setText("")
    with pytest.raises(ValueError, match="CSV folder"):
        win._resolve_source()


def test_resolve_source_db_empty_raises(win):
    win.cb_source.setCurrentIndex(3)  # Databroker
    win.le_db_name.setText("")
    with pytest.raises(ValueError, match="catalog name"):
        win._resolve_source()


def test_resolve_source_db_calls_load_catalog(win):
    win.cb_source.setCurrentIndex(3)
    win.le_db_name.setText("my_catalog")
    fake_cat = MagicMock()
    with patch(
        "polartools.load_data.load_catalog", return_value=fake_cat
    ) as mock_lc:
        src, kw = win._resolve_source()
    mock_lc.assert_called_once_with("my_catalog")
    assert src is fake_cat
    assert kw == {}


# ─── _schedule_normalize ─────────────────────────────────────────────────────


def test_schedule_normalize_no_data(win):
    win._energy = None
    win._schedule_normalize()
    assert not win._norm_timer.isActive()


def test_schedule_normalize_with_data(loaded_win):
    loaded_win._schedule_normalize()
    assert loaded_win._norm_timer.isActive()
    loaded_win._norm_timer.stop()


# ─── _run_normalization ───────────────────────────────────────────────────────


def test_run_normalization_no_data_noop(win):
    win._run_normalization()
    assert not win.btn_save.isEnabled()


def test_run_normalization_success(loaded_win):
    energy = loaded_win._energy
    with patch(
        "polartools.xanes_gui.normalize_absorption",
        return_value=_norm_stub(energy),
    ):
        loaded_win._run_normalization()
    assert loaded_win.btn_save.isEnabled()
    assert loaded_win._results is not None


def test_run_normalization_error(loaded_win):
    with patch(
        "polartools.xanes_gui.normalize_absorption",
        side_effect=RuntimeError("bad fit"),
    ):
        loaded_win._run_normalization()
    msg = loaded_win.status_bar.currentMessage()
    assert "Normalization error" in msg


def test_run_normalization_error_clears_stale_results(loaded_win):
    energy = loaded_win._energy
    with patch(
        "polartools.xanes_gui.normalize_absorption",
        return_value=_norm_stub(energy),
    ):
        loaded_win._run_normalization()
    assert loaded_win._results is not None
    assert loaded_win.btn_save.isEnabled()

    with patch(
        "polartools.xanes_gui.normalize_absorption",
        side_effect=RuntimeError("bad fit"),
    ):
        loaded_win._run_normalization()
    assert loaded_win._results is None
    assert not loaded_win.btn_save.isEnabled()


def test_run_normalization_syncs_pre_post_lines_to_new_e0(loaded_win):
    energy = loaded_win._energy
    loaded_win.le_pre1.setText("-30.0")
    loaded_win.le_post1.setText("10.0")

    with patch(
        "polartools.xanes_gui.normalize_absorption",
        return_value=_norm_stub(energy, e0=7150.0),
    ):
        loaded_win._run_normalization()

    assert loaded_win.line_pre1.value() == pytest.approx(7150.0 - 30.0)
    assert loaded_win.line_post1.value() == pytest.approx(7150.0 + 10.0)


# ─── _save_results ────────────────────────────────────────────────────────────


def test_save_results_no_results_noop(win):
    with patch("polartools.xanes_gui.save_xas") as mock_save:
        win._save_results()
    mock_save.assert_not_called()


def test_save_results_absolute_path(normalized_win, tmp_path):
    fname = str(tmp_path / "out.dat")
    normalized_win.le_savename.setText(fname)
    with patch("polartools.xanes_gui.save_xas") as mock_save:
        normalized_win._save_results()
    mock_save.assert_called_once_with(normalized_win._results, fname)
    assert "Saved" in normalized_win.status_bar.currentMessage()


def test_save_results_save_error(normalized_win, tmp_path):
    fname = str(tmp_path / "out.dat")
    normalized_win.le_savename.setText(fname)
    with patch(
        "polartools.xanes_gui.save_xas",
        side_effect=OSError("disk full"),
    ):
        with patch("polartools.xanes_gui.QMessageBox.critical") as mock_err:
            normalized_win._save_results()
    mock_err.assert_called_once()


# ─── _resolve_save_path (folder selection) ───────────────────────────────────


def test_resolve_save_path_absolute_name_ignores_folder(win):
    win.le_savefolder.setText("/some/folder")
    win.le_savename.setText("/abs/out.dat")
    path, resolved = win._resolve_save_path()
    assert path == "/abs/out.dat"
    assert resolved is True


def test_resolve_save_path_folder_plus_relative_name(win):
    win.le_savefolder.setText("/some/folder")
    win.le_savename.setText("out.dat")
    path, resolved = win._resolve_save_path()
    assert path == os.path.join("/some/folder", "out.dat")
    assert resolved is True


def test_resolve_save_path_no_folder_needs_dialog(win):
    win.le_savefolder.setText("")
    win.le_savename.setText("out.dat")
    path, resolved = win._resolve_save_path()
    assert path == "out.dat"
    assert resolved is False


def test_save_results_uses_folder_without_dialog(normalized_win, tmp_path):
    normalized_win.le_savefolder.setText(str(tmp_path))
    normalized_win.le_savename.setText("out.dat")
    expected = str(tmp_path / "out.dat")
    with (
        patch("polartools.xanes_gui.save_xas") as mock_save,
        patch("polartools.xanes_gui.QFileDialog.getSaveFileName") as mock_dlg,
    ):
        normalized_win._save_results()
    mock_dlg.assert_not_called()
    mock_save.assert_called_once_with(normalized_win._results, expected)


# ─── Overwrite confirmation ────────────────────────────────────────────────────


def test_save_results_prompts_when_file_exists(normalized_win, tmp_path):
    fname = tmp_path / "out.dat"
    fname.write_text("existing content")
    normalized_win.le_savename.setText(str(fname))
    with (
        patch("polartools.xanes_gui.save_xas") as mock_save,
        patch(
            "polartools.xanes_gui.QMessageBox.question",
            return_value=QMessageBox.StandardButton.Yes,
        ) as mock_q,
    ):
        normalized_win._save_results()
    mock_q.assert_called_once()
    mock_save.assert_called_once_with(normalized_win._results, str(fname))


def test_save_results_overwrite_declined_skips_save(normalized_win, tmp_path):
    fname = tmp_path / "out.dat"
    fname.write_text("existing content")
    normalized_win.le_savename.setText(str(fname))
    with (
        patch("polartools.xanes_gui.save_xas") as mock_save,
        patch(
            "polartools.xanes_gui.QMessageBox.question",
            return_value=QMessageBox.StandardButton.No,
        ),
    ):
        normalized_win._save_results()
    mock_save.assert_not_called()
    assert "cancelled" in normalized_win.status_bar.currentMessage().lower()


def test_save_results_no_prompt_when_file_missing(normalized_win, tmp_path):
    fname = str(tmp_path / "new_out.dat")
    normalized_win.le_savename.setText(fname)
    with (
        patch("polartools.xanes_gui.save_xas") as mock_save,
        patch("polartools.xanes_gui.QMessageBox.question") as mock_q,
    ):
        normalized_win._save_results()
    mock_q.assert_not_called()
    mock_save.assert_called_once_with(normalized_win._results, fname)


def test_save_results_dialog_path_no_extra_prompt(normalized_win, tmp_path):
    # QFileDialog.getSaveFileName already confirms overwrite natively, so the
    # unresolved (dialog) path must not add a second confirmation.
    fname = tmp_path / "out.dat"
    fname.write_text("existing content")
    normalized_win.le_savefolder.setText("")
    normalized_win.le_savename.setText("out.dat")
    with (
        patch("polartools.xanes_gui.save_xas") as mock_save,
        patch(
            "polartools.xanes_gui.QFileDialog.getSaveFileName",
            return_value=(str(fname), ""),
        ),
        patch("polartools.xanes_gui.QMessageBox.question") as mock_q,
    ):
        normalized_win._save_results()
    mock_q.assert_not_called()
    mock_save.assert_called_once_with(normalized_win._results, str(fname))


# ─── _load_reference / _clear_reference ────────────────────────────────────────


def _write_reference_file(path):
    energy = np.linspace(7000, 7300, 50)
    mu = np.linspace(0, 1, 50)
    norm = np.linspace(0, 1.1, 50)
    flat = np.linspace(0, 1.2, 50)
    data = np.vstack((energy, mu, norm, flat)).transpose()
    np.savetxt(
        path,
        data,
        header="XANES\nEnergy\tXANES\tNormalized\tFlattened",
        fmt="%0.5e",
    )
    return energy, norm


def test_load_reference_no_path_warns(win):
    with patch("polartools.xanes_gui.QMessageBox.warning") as mock_warn:
        win._load_reference()
    mock_warn.assert_called_once()
    assert win._ref_energy is None


def test_load_reference_success(win, tmp_path):
    fname = str(tmp_path / "reference.dat")
    energy, norm = _write_reference_file(fname)
    win.le_reference.setText(fname)
    win._load_reference()
    np.testing.assert_allclose(win._ref_energy, energy, rtol=1e-4)
    np.testing.assert_allclose(win._ref_norm, norm, rtol=1e-4)
    x, y = win.curve_reference.getData()
    np.testing.assert_allclose(x, energy, rtol=1e-4)
    np.testing.assert_allclose(y, norm, rtol=1e-4)
    assert win.btn_clear_reference.isEnabled()
    assert "Reference loaded" in win.status_bar.currentMessage()


def test_load_reference_missing_file_shows_error(win):
    win.le_reference.setText("/no/such/file.dat")
    with patch("polartools.xanes_gui.QMessageBox.critical") as mock_err:
        win._load_reference()
    mock_err.assert_called_once()
    assert win._ref_energy is None


def test_load_reference_bad_shape_shows_error(win, tmp_path):
    fname = str(tmp_path / "bad.dat")
    np.savetxt(fname, np.linspace(0, 1, 10))
    win.le_reference.setText(fname)
    with patch("polartools.xanes_gui.QMessageBox.critical") as mock_err:
        win._load_reference()
    mock_err.assert_called_once()
    assert win._ref_energy is None


def test_clear_reference_resets_state(win, tmp_path):
    fname = str(tmp_path / "reference.dat")
    _write_reference_file(fname)
    win.le_reference.setText(fname)
    win._load_reference()
    win._clear_reference()
    assert win._ref_energy is None
    assert win._ref_norm is None
    assert not win.btn_clear_reference.isEnabled()
    x, y = win.curve_reference.getData()
    assert x is None or len(x) == 0


# ─── _init_markers ───────────────────────────────────────────────────────────


def test_init_markers_lines_visible(loaded_win):
    loaded_win._init_markers()
    for line in (
        loaded_win.line_e0,
        loaded_win.line_pre1,
        loaded_win.line_pre2,
        loaded_win.line_post1,
        loaded_win.line_post2,
    ):
        assert line.isVisible()


def test_init_markers_entries_populated(loaded_win):
    loaded_win._init_markers()
    for le in (
        loaded_win.le_pre1,
        loaded_win.le_pre2,
        loaded_win.le_post1,
        loaded_win.le_post2,
    ):
        assert le.text() != ""


def test_init_markers_e0_val_set(loaded_win):
    loaded_win._e0_val = None
    loaded_win._init_markers()
    assert loaded_win._e0_val is not None


# ─── _line_moved ─────────────────────────────────────────────────────────────


def test_line_moved_e0_none_noop(win):
    win._e0_val = None
    old_text = win.le_pre1.text()
    win._line_moved(win.line_pre1, win.le_pre1)
    assert win.le_pre1.text() == old_text


def test_line_moved_block_flag_noop(loaded_win):
    loaded_win._block_line_update = True
    old_text = loaded_win.le_pre1.text()
    loaded_win._line_moved(loaded_win.line_pre1, loaded_win.le_pre1)
    assert loaded_win.le_pre1.text() == old_text
    loaded_win._block_line_update = False


def test_line_moved_updates_entry(loaded_win):
    loaded_win.line_pre1.setPos(7000.0)
    loaded_win._line_moved(loaded_win.line_pre1, loaded_win.le_pre1)
    expected = f"{7000.0 - loaded_win._e0_val:.1f}"
    assert loaded_win.le_pre1.text() == expected


# ─── _entry_changed ──────────────────────────────────────────────────────────


def test_entry_changed_empty_text_noop(loaded_win):
    loaded_win.le_pre1.setText("")
    old_pos = loaded_win.line_pre1.value()
    loaded_win._entry_changed(loaded_win.le_pre1, loaded_win.line_pre1)
    assert loaded_win.line_pre1.value() == old_pos


def test_entry_changed_invalid_text_noop(loaded_win):
    loaded_win.le_pre1.setText("abc")
    old_pos = loaded_win.line_pre1.value()
    loaded_win._entry_changed(loaded_win.le_pre1, loaded_win.line_pre1)
    assert loaded_win.line_pre1.value() == old_pos


def test_entry_changed_valid_text(loaded_win):
    loaded_win.le_pre1.setText("-100.0")
    loaded_win._entry_changed(loaded_win.le_pre1, loaded_win.line_pre1)
    assert loaded_win.line_pre1.value() == pytest.approx(
        loaded_win._e0_val - 100.0
    )


# ─── _on_load ────────────────────────────────────────────────────────────────


def test_on_load_parse_error_shows_dialog(win):
    win.le_scans.setText("")
    with patch("polartools.xanes_gui.QMessageBox.critical") as mock_msg:
        win._on_load()
    mock_msg.assert_called_once()


def test_on_load_source_error_shows_dialog(win):
    win.le_scans.setText("1")
    win.cb_source.setCurrentIndex(0)  # SPEC
    win.le_spec_path.setText("")  # empty → ValueError in _resolve_source
    with patch("polartools.xanes_gui.QMessageBox.critical") as mock_msg:
        win._on_load()
    mock_msg.assert_called_once()


def test_on_load_load_error_shows_dialog(win):
    win.le_scans.setText("1")
    win.cb_source.setCurrentIndex(0)
    win.le_spec_path.setText("/fake/path.dat")
    with patch(
        "polartools.xanes_gui.load_multi_xas",
        side_effect=OSError("file not found"),
    ):
        with patch("polartools.xanes_gui.QMessageBox.critical") as mock_msg:
            win._on_load()
    mock_msg.assert_called_once()


def test_on_load_success_sets_state(win, synthetic_data):
    e, mu = synthetic_data
    win.le_scans.setText("1")
    win.cb_source.setCurrentIndex(0)
    win.le_spec_path.setText("/fake/path.dat")

    fake_return = (e, mu, np.zeros_like(mu))
    with patch("polartools.xanes_gui.load_multi_xas", return_value=fake_return):
        with patch(
            "polartools.xanes_gui.normalize_absorption",
            side_effect=lambda *a, **k: _norm_stub(e * 1000),
        ):
            win._on_load()

    assert win._energy is not None
    assert win._results is not None


# ─── Guess-only-on-first-load behavior ───────────────────────────────────────


def test_on_load_guesses_only_first_load(win, synthetic_data):
    e, mu = synthetic_data
    win.le_scans.setText("1")
    win.cb_source.setCurrentIndex(0)
    win.le_spec_path.setText("/fake/path.dat")
    fake_return = (e, mu, np.zeros_like(mu))
    with (
        patch("polartools.xanes_gui.load_multi_xas", return_value=fake_return),
        patch(
            "polartools.xanes_gui.normalize_absorption",
            side_effect=lambda *a, **k: _norm_stub(e * 1000),
        ),
        patch.object(win, "_init_markers", wraps=win._init_markers) as spy,
    ):
        win._on_load()
        assert spy.call_count == 1
        assert win._markers_initialized
        win._on_load()
        # Parameters kept — not re-guessed on the second load.
        assert spy.call_count == 1


def test_on_guess_noop_without_data(win):
    with patch.object(win, "_init_markers") as spy:
        win._on_guess()
    spy.assert_not_called()


def test_on_guess_reguesses_with_data(loaded_win):
    energy = loaded_win._energy
    with (
        patch(
            "polartools.xanes_gui.normalize_absorption",
            side_effect=lambda *a, **k: _norm_stub(energy),
        ),
        patch.object(loaded_win, "_init_markers") as spy,
    ):
        loaded_win._on_guess()
    spy.assert_called_once()


# ─── Column File source ───────────────────────────────────────────────────────


def _write_column_file(path, ncols=2, nrows=20):
    energy = np.linspace(7100, 7300, nrows)
    cols = [energy] + [
        np.random.default_rng(0).random(nrows) for _ in range(ncols - 1)
    ]
    data = np.column_stack(cols)
    np.savetxt(path, data, header="a comment line")
    return data


def test_column_source_initial_state(win):
    assert win.cb_energy_col.count() == 1
    assert win.cb_mu_col.count() == 1
    assert win.cb_energy_col.itemText(0) == "0"
    assert win.cb_mu_col.itemText(0) == "1"


def test_source_stack_has_column_file_entry(win):
    assert win.cb_source.itemText(win.cb_source.count() - 1) == "Column File"


def test_preview_column_file_populates_columns(win, tmp_path):
    path = str(tmp_path / "data.dat")
    _write_column_file(path, ncols=4)
    win._preview_column_file(path)
    assert win.cb_energy_col.count() == 4
    assert win.cb_energy_col.currentIndex() == 0
    assert win.cb_mu_col.currentIndex() == 1


def test_preview_column_file_single_column(win, tmp_path):
    path = str(tmp_path / "data.dat")
    _write_column_file(path, ncols=1)
    win._preview_column_file(path)
    assert win.cb_energy_col.count() == 1
    assert win.cb_mu_col.currentIndex() == 0


def test_preview_column_file_missing_shows_error(win):
    with patch("polartools.xanes_gui.QMessageBox.critical") as mock_msg:
        win._preview_column_file("/no/such/file.dat")
    mock_msg.assert_called_once()


def test_browse_column_file_cancelled_noop(win):
    with (
        patch(
            "polartools.xanes_gui.QFileDialog.getOpenFileName",
            return_value=("", ""),
        ),
        patch.object(win, "_preview_column_file") as spy,
    ):
        win._browse_column_file()
    spy.assert_not_called()
    assert win.le_column_path.text() == ""


def test_browse_column_file_sets_path_and_previews(win, tmp_path):
    path = str(tmp_path / "data.dat")
    _write_column_file(path)
    with patch(
        "polartools.xanes_gui.QFileDialog.getOpenFileName",
        return_value=(path, "Data files (*.dat)"),
    ):
        win._browse_column_file()
    assert win.le_column_path.text() == path
    assert win.cb_energy_col.count() == 2


def test_on_column_path_edited_empty_noop(win):
    win.le_column_path.setText("")
    with patch.object(win, "_preview_column_file") as spy:
        win._on_column_path_edited()
    spy.assert_not_called()


def test_on_column_path_edited_triggers_preview(win, tmp_path):
    path = str(tmp_path / "data.dat")
    _write_column_file(path)
    win.le_column_path.setText(path)
    with patch.object(win, "_preview_column_file") as spy:
        win._on_column_path_edited()
    spy.assert_called_once_with(path)


def test_on_source_changed_disables_scan_widgets(win):
    win.cb_source.setCurrentIndex(win.cb_source.count() - 1)  # Column File
    assert not win.le_scans.isEnabled()
    assert not win._xas_params_widget.isEnabled()

    win.cb_source.setCurrentIndex(0)  # SPEC
    assert win.le_scans.isEnabled()
    assert win._xas_params_widget.isEnabled()


def test_on_load_column_file_no_path_shows_error(win):
    win.cb_source.setCurrentIndex(win.cb_source.count() - 1)
    win.le_column_path.setText("")
    with patch("polartools.xanes_gui.QMessageBox.critical") as mock_msg:
        win._on_load()
    mock_msg.assert_called_once()


def test_on_load_column_file_bad_file_shows_error(win):
    win.cb_source.setCurrentIndex(win.cb_source.count() - 1)
    win.le_column_path.setText("/no/such/file.dat")
    with patch("polartools.xanes_gui.QMessageBox.critical") as mock_msg:
        win._on_load()
    mock_msg.assert_called_once()


def test_on_load_column_file_out_of_range_column(win, tmp_path):
    path = str(tmp_path / "data.dat")
    _write_column_file(path, ncols=2)
    win.cb_source.setCurrentIndex(win.cb_source.count() - 1)
    win.le_column_path.setText(path)
    win._preview_column_file(path)
    with QSignalBlocker(win.cb_mu_col):
        win.cb_mu_col.setCurrentIndex(-1)
    with patch("polartools.xanes_gui.QMessageBox.critical") as mock_msg:
        win._on_load()
    mock_msg.assert_called_once()


def test_on_load_column_file_success(win, tmp_path):
    path = str(tmp_path / "data.dat")
    data = _write_column_file(path, ncols=2)
    win.cb_source.setCurrentIndex(win.cb_source.count() - 1)
    win.le_column_path.setText(path)
    win._preview_column_file(path)

    with patch(
        "polartools.xanes_gui.normalize_absorption",
        side_effect=lambda *a, **k: _norm_stub(win._energy),
    ):
        win._on_load()

    assert win._energy is not None
    assert win._mu is not None
    # No keV->eV conversion — energy taken as-is from the file.
    assert np.allclose(sorted(data[:, 0]), win._energy)


def test_on_load_column_file_reuses_preview_cache(win, tmp_path):
    path = str(tmp_path / "data.dat")
    _write_column_file(path, ncols=2)
    win.cb_source.setCurrentIndex(win.cb_source.count() - 1)
    win.le_column_path.setText(path)
    win._preview_column_file(path)

    with (
        patch("polartools.xanes_gui.np.loadtxt") as mock_loadtxt,
        patch(
            "polartools.xanes_gui.normalize_absorption",
            side_effect=lambda *a, **k: _norm_stub(win._energy),
        ),
    ):
        win._on_load()

    mock_loadtxt.assert_not_called()
    assert win._energy is not None


def test_on_load_column_file_reparses_on_cache_miss(win, tmp_path):
    path = str(tmp_path / "data.dat")
    _write_column_file(path, ncols=2)
    win.cb_source.setCurrentIndex(win.cb_source.count() - 1)
    win.le_column_path.setText(path)
    # No preview taken for this path — cache stays empty/mismatched.
    assert win._column_cache is None

    with patch(
        "polartools.xanes_gui.normalize_absorption",
        side_effect=lambda *a, **k: _norm_stub(win._energy),
    ):
        win._on_load()

    assert win._energy is not None
    assert win._results is not None
