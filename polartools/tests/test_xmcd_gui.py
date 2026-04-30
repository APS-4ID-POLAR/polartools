# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

import pytest
import numpy as np
from unittest.mock import patch, MagicMock

# Skip the entire module if PyQt6 or pyqtgraph are not installed.
pytest.importorskip("PyQt6", reason="PyQt6 not installed — skipping GUI tests")
pytest.importorskip(
    "pyqtgraph", reason="pyqtgraph not installed — skipping GUI tests"
)


@pytest.fixture(scope="session")
def qapp():
    """Session-scoped QApplication — created once, reused across all tests."""
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def win(qapp):
    from polartools.xmcd_gui import MainWindow

    window = MainWindow()
    yield window
    window.close()


@pytest.fixture
def synthetic_data():
    """Fake Fe-edge arrays: energy in keV, mu a sigmoid step, xmcd zeros."""
    energy = np.linspace(7.0, 7.3, 200)  # keV (GUI multiplies by 1000)
    mu = 1 / (1 + np.exp(-(energy * 1000 - 7112) / 3))
    xmcd_raw = np.zeros_like(energy)
    return energy, mu, xmcd_raw


@pytest.fixture
def loaded_win(win, synthetic_data):
    """Window with synthetic data injected — simulates state after a load."""
    ep, yp, zp = synthetic_data
    win._plus_energy = ep * 1000
    win._plus_mu = yp
    win._plus_xmcd_raw = zp
    win._minus_energy = ep * 1000
    win._minus_mu = yp
    win._minus_xmcd_raw = zp
    win._e0_val = 7112.0
    return win


@pytest.fixture
def normalized_win(loaded_win):
    """Window with results dicts set — simulates state after normalization."""
    energy = loaded_win._plus_energy
    stub = {
        "energy": energy,
        "mu": loaded_win._plus_mu,
        "norm": np.ones_like(energy),
        "xmcd": np.zeros_like(energy),
        "preedge": np.zeros_like(energy),
        "postedge": np.ones_like(energy),
        "e0": 7112.0,
        "edge_step": 1.0,
        "flat": np.ones_like(energy),
    }
    loaded_win._plus_results = stub.copy()
    loaded_win._minus_results = stub.copy()
    return loaded_win


# ─── Window creation ──────────────────────────────────────────────────────────


def test_window_title(win):
    assert win.windowTitle() == "XMCD Processor"


def test_initial_state_no_data(win):
    assert win._plus_energy is None
    assert win._minus_energy is None
    assert win._plus_results is None
    assert win._minus_results is None


def test_save_button_disabled_on_start(win):
    assert not win.btn_save.isEnabled()


# ─── Source / kind stacked widgets ───────────────────────────────────────────


def test_source_stack_follows_combo(win):
    for i in range(win.cb_source.count()):
        win.cb_source.setCurrentIndex(i)
        assert win._source_stack.currentIndex() == i


def test_kind_stack_follows_combo(win):
    win.cb_kind.setCurrentIndex(0)  # dichro
    assert win._kind_stack.currentIndex() == 0
    win.cb_kind.setCurrentIndex(1)  # lockin
    assert win._kind_stack.currentIndex() == 1


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


# ─── _resolve_load_kwargs (dichro / lockin) ───────────────────────────────────


def test_resolve_load_kwargs_dichro_defaults(win):
    win.cb_kind.setCurrentIndex(0)  # dichro
    kw = win._resolve_load_kwargs({})
    # All optional fields blank → only transmission kwarg present
    assert kw["transmission"] is True
    assert "positioner" not in kw
    assert "detector" not in kw
    assert "monitor" not in kw


def test_resolve_load_kwargs_dichro_custom(win):
    win.cb_kind.setCurrentIndex(0)
    win.le_d_positioner.setText("4idenergy")
    win.le_d_detector.setText("IC5")
    win.chk_transmission.setChecked(False)
    kw = win._resolve_load_kwargs({})
    assert kw["positioner"] == "4idenergy"
    assert kw["detector"] == "IC5"
    assert kw["transmission"] is False


def test_resolve_load_kwargs_dichro_monitor(win):
    win.cb_kind.setCurrentIndex(0)
    win.le_d_monitor.setText("IC4")
    kw = win._resolve_load_kwargs({})
    assert kw["monitor"] == "IC4"


def test_resolve_load_kwargs_lockin_defaults(win):
    win.cb_kind.setCurrentIndex(1)  # lockin
    kw = win._resolve_load_kwargs({})
    assert "positioner" not in kw
    assert "dc_col" not in kw


def test_resolve_load_kwargs_lockin_custom(win):
    win.cb_kind.setCurrentIndex(1)
    win.le_l_dc.setText("Lock DC")
    win.le_l_ac.setText("Lock AC")
    kw = win._resolve_load_kwargs({})
    assert kw["dc_col"] == "Lock DC"
    assert kw["ac_col"] == "Lock AC"


def test_resolve_load_kwargs_lockin_acoff(win):
    win.cb_kind.setCurrentIndex(1)
    win.le_l_acoff.setText("Lock AC off")
    kw = win._resolve_load_kwargs({})
    assert kw["acoff_col"] == "Lock AC off"


def test_resolve_load_kwargs_passes_extra(win):
    win.cb_kind.setCurrentIndex(0)
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
    win._plus_energy = None
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
    energy = loaded_win._plus_energy
    stub = {
        "energy": energy,
        "mu": loaded_win._plus_mu,
        "norm": np.ones_like(energy),
        "xmcd": np.zeros_like(energy),
        "preedge": np.zeros_like(energy),
        "postedge": np.ones_like(energy),
        "e0": 7112.0,
        "edge_step": 1.0,
        "flat": np.ones_like(energy),
    }
    with patch("polartools.xmcd_gui.normalize_absorption", return_value=stub):
        loaded_win._run_normalization()
    assert loaded_win.btn_save.isEnabled()
    assert loaded_win._plus_results is not None
    assert loaded_win._minus_results is not None


def test_run_normalization_error(loaded_win):
    with patch(
        "polartools.xmcd_gui.normalize_absorption",
        side_effect=RuntimeError("bad fit"),
    ):
        loaded_win._run_normalization()
    msg = loaded_win.status_bar.currentMessage()
    assert "Normalization error" in msg


# ─── _save_results ────────────────────────────────────────────────────────────


def test_save_results_no_results_noop(win):
    with patch("polartools.xmcd_gui.save_xmcd") as mock_save:
        win._save_results()
    mock_save.assert_not_called()


def test_save_results_absolute_path(normalized_win, tmp_path):
    fname = str(tmp_path / "out.dat")
    normalized_win.le_savename.setText(fname)
    with patch("polartools.xmcd_gui.save_xmcd") as mock_save:
        normalized_win._save_results()
    mock_save.assert_called_once_with(
        normalized_win._plus_results, normalized_win._minus_results, fname
    )
    assert "Saved" in normalized_win.status_bar.currentMessage()


def test_save_results_save_error(normalized_win, tmp_path):
    fname = str(tmp_path / "out.dat")
    normalized_win.le_savename.setText(fname)
    with patch(
        "polartools.xmcd_gui.save_xmcd",
        side_effect=OSError("disk full"),
    ):
        with patch("polartools.xmcd_gui.QMessageBox.critical") as mock_err:
            normalized_win._save_results()
    mock_err.assert_called_once()


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
    win.le_scans_plus.setText("")
    win.le_scans_minus.setText("2")
    with patch("polartools.xmcd_gui.QMessageBox.critical") as mock_msg:
        win._on_load()
    mock_msg.assert_called_once()


def test_on_load_source_error_shows_dialog(win):
    win.le_scans_plus.setText("1")
    win.le_scans_minus.setText("2")
    win.cb_source.setCurrentIndex(0)  # SPEC
    win.le_spec_path.setText("")  # empty → ValueError in _resolve_source
    with patch("polartools.xmcd_gui.QMessageBox.critical") as mock_msg:
        win._on_load()
    mock_msg.assert_called_once()


def test_on_load_load_error_shows_dialog(win):
    win.le_scans_plus.setText("1")
    win.le_scans_minus.setText("2")
    win.cb_source.setCurrentIndex(0)
    win.le_spec_path.setText("/fake/path.dat")
    with patch(
        "polartools.xmcd_gui.load_multi_dichro",
        side_effect=OSError("file not found"),
    ):
        with patch("polartools.xmcd_gui.QMessageBox.critical") as mock_msg:
            win._on_load()
    mock_msg.assert_called_once()


def test_on_load_success_sets_state(win, synthetic_data):
    ep, yp, zp = synthetic_data
    win.le_scans_plus.setText("1")
    win.le_scans_minus.setText("2")
    win.cb_source.setCurrentIndex(0)
    win.le_spec_path.setText("/fake/path.dat")
    win.cb_kind.setCurrentIndex(0)  # dichro

    fake_return = (ep, yp, zp, None, None)
    with patch(
        "polartools.xmcd_gui.load_multi_dichro", return_value=fake_return
    ):
        with patch("polartools.xmcd_gui.normalize_absorption") as mock_norm:
            mock_norm.return_value = {
                "energy": ep * 1000,
                "mu": yp,
                "norm": np.ones_like(yp),
                "xmcd": zp,
                "preedge": np.zeros_like(yp),
                "postedge": np.ones_like(yp),
                "e0": 7112.0,
                "edge_step": 1.0,
                "flat": np.ones_like(yp),
            }
            win._on_load()

    assert win._plus_energy is not None
    assert win._minus_energy is not None


# ─── Independent H+/H− normalization ─────────────────────────────────────────


def test_independent_panel_hidden_by_default(win):
    assert not win.chk_independent.isChecked()
    assert not win._minus_norm_panel.isVisibleTo(win)


def test_independent_panel_shows_when_toggled(win):
    win.show()
    win.chk_independent.setChecked(True)
    assert win._minus_norm_panel.isVisibleTo(win)
    win.chk_independent.setChecked(False)
    assert not win._minus_norm_panel.isVisibleTo(win)


def test_h_minus_markers_hidden_by_default(win):
    for line in (
        win.m_line_pre1,
        win.m_line_pre2,
        win.m_line_post1,
        win.m_line_post2,
    ):
        assert not line.isVisible()


def test_h_minus_markers_show_when_toggled(loaded_win):
    loaded_win.chk_independent.setChecked(True)
    for line in (
        loaded_win.m_line_pre1,
        loaded_win.m_line_pre2,
        loaded_win.m_line_post1,
        loaded_win.m_line_post2,
    ):
        assert line.isVisible()


def test_run_normalization_independent_uses_separate_kwargs(loaded_win):
    loaded_win.chk_independent.setChecked(True)
    loaded_win.chk_e0_auto.setChecked(False)
    loaded_win.le_e0.setText("7112.0")
    loaded_win.m_chk_e0_auto.setChecked(False)
    loaded_win.m_le_e0.setText("7115.0")

    energy = loaded_win._plus_energy
    stub = {
        "energy": energy,
        "mu": loaded_win._plus_mu,
        "norm": np.ones_like(energy),
        "xmcd": np.zeros_like(energy),
        "preedge": np.zeros_like(energy),
        "postedge": np.ones_like(energy),
        "e0": 7112.0,
        "edge_step": 1.0,
        "flat": np.ones_like(energy),
    }
    with patch(
        "polartools.xmcd_gui.normalize_absorption", return_value=stub
    ) as mock_norm:
        loaded_win._run_normalization()

    assert mock_norm.call_count == 2
    plus_e0 = mock_norm.call_args_list[0].kwargs["e0"]
    minus_e0 = mock_norm.call_args_list[1].kwargs["e0"]
    assert plus_e0 == pytest.approx(7112.0)
    assert minus_e0 == pytest.approx(7115.0)


def test_run_normalization_dependent_uses_same_kwargs(loaded_win):
    loaded_win.chk_independent.setChecked(False)
    loaded_win.chk_e0_auto.setChecked(False)
    loaded_win.le_e0.setText("7112.0")
    loaded_win.m_chk_e0_auto.setChecked(False)
    loaded_win.m_le_e0.setText("7115.0")  # ignored when toggle is off

    energy = loaded_win._plus_energy
    stub = {
        "energy": energy,
        "mu": loaded_win._plus_mu,
        "norm": np.ones_like(energy),
        "xmcd": np.zeros_like(energy),
        "preedge": np.zeros_like(energy),
        "postedge": np.ones_like(energy),
        "e0": 7112.0,
        "edge_step": 1.0,
        "flat": np.ones_like(energy),
    }
    with patch(
        "polartools.xmcd_gui.normalize_absorption", return_value=stub
    ) as mock_norm:
        loaded_win._run_normalization()

    assert mock_norm.call_count == 2
    plus_kw = mock_norm.call_args_list[0].kwargs
    minus_kw = mock_norm.call_args_list[1].kwargs
    assert plus_kw == minus_kw
    assert plus_kw["e0"] == pytest.approx(7112.0)
