# Copyright (c) 2026, UChicago Argonne, LLC.
# See LICENSE file for details.

from os.path import join
from unittest.mock import patch

import numpy as np
import pytest

# Skip the entire module if PyQt6 is not installed.
pytest.importorskip("PyQt6", reason="PyQt6 not installed — skipping GUI tests")

FOLDER = join("polartools", "tests", "data_for_test")


@pytest.fixture(scope="session")
def qapp():
    """Session-scoped QApplication — created once, reused across all tests."""
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def win(qapp):
    from polartools.process_images_gui import MainWindow

    window = MainWindow()
    yield window
    window.close()


@pytest.fixture
def hdf5_win(win):
    win.cb_source.setCurrentIndex(1)
    win.le_folder.setText(FOLDER)
    return win


# ─── Window creation ──────────────────────────────────────────────────────────


def test_window_title(win):
    assert win.windowTitle() == "RXES Image Processor"


def test_three_tabs_present(win):
    assert win.tabs.count() == 3
    titles = [win.tabs.tabText(i) for i in range(win.tabs.count())]
    assert titles == ["Curvature", "RXES", "RXES-MCD"]


def test_curvature_save_disabled_on_start(win):
    assert not win.btn_curv_save.isEnabled()


def test_rxes_save_disabled_on_start(win):
    assert not win.btn_rxes_save.isEnabled()


def test_mcd_save_disabled_on_start(win):
    assert not win.btn_mcd_save.isEnabled()


# ─── Source panel ─────────────────────────────────────────────────────────────


def test_source_stack_follows_combo(win):
    for i in range(win.cb_source.count()):
        win.cb_source.setCurrentIndex(i)
        assert win._source_stack.currentIndex() == i


def test_resolve_cat_kwargs_hdf5(hdf5_win):
    cat, kwargs = hdf5_win._resolve_cat_kwargs()
    assert cat == "hdf5"
    assert kwargs == {"folder": FOLDER}


def test_resolve_cat_kwargs_hdf5_empty_folder_raises(win):
    win.cb_source.setCurrentIndex(1)
    win.le_folder.setText("")
    with pytest.raises(ValueError, match="HDF5 folder"):
        win._resolve_cat_kwargs()


def test_resolve_cat_kwargs_catalog_empty_name_raises(win):
    win.cb_source.setCurrentIndex(0)
    win.le_catalog.setText("")
    with pytest.raises(ValueError, match="Catalog name"):
        win._resolve_cat_kwargs()


def test_resolve_cat_kwargs_catalog_calls_load_catalog(win):
    from unittest.mock import MagicMock

    win.cb_source.setCurrentIndex(0)
    win.le_catalog.setText("my_catalog")
    fake_cat = MagicMock()
    with patch(
        "polartools.load_data.load_catalog", return_value=fake_cat
    ) as mock_lc:
        cat, kwargs = win._resolve_cat_kwargs()
    mock_lc.assert_called_once_with("my_catalog")
    assert cat is fake_cat
    assert kwargs == {}


# ─── Curvature tab ────────────────────────────────────────────────────────────


def test_curvature_tab_hdf5(hdf5_win):
    hdf5_win._curv["scans"].setText("322")
    hdf5_win.sb_binx.setValue(64)
    hdf5_win._on_curvature_run()

    assert hdf5_win._curvature is not None
    assert len(hdf5_win._curvature) == 3
    assert hdf5_win.btn_curv_save.isEnabled()
    assert hdf5_win.le_curv_result.text() != ""


def test_curvature_tab_input_error_shows_dialog(win):
    win._curv["scans"].setText("")
    with patch(
        "polartools.process_images_gui.QMessageBox.critical"
    ) as mock_msg:
        win._on_curvature_run()
    mock_msg.assert_called_once()


def test_curvature_save_noop_without_result(win, tmp_path):
    with patch(
        "polartools.process_images_gui.QFileDialog.getSaveFileName"
    ) as mock_dialog:
        win._on_curvature_save()
    mock_dialog.assert_not_called()


def test_curvature_save_writes_file(win, tmp_path):
    win._curvature = np.array([1.0, 2.0, 3.0])
    path = str(tmp_path / "curvature.txt")
    with patch(
        "polartools.process_images_gui.QFileDialog.getSaveFileName",
        return_value=(path, ""),
    ):
        win._on_curvature_save()
    assert np.allclose(np.loadtxt(path), [1.0, 2.0, 3.0])


# ─── RXES tab ─────────────────────────────────────────────────────────────────


def test_rxes_tab_no_positioner(hdf5_win):
    hdf5_win._rxes["scans"].setText("322")
    hdf5_win._rxes_curv[0].setText("1")
    hdf5_win._rxes_curv[1].setText("2")
    hdf5_win._rxes_curv[2].setText("3")
    hdf5_win._on_rxes_run()

    positioner, result = hdf5_win._rxes_result
    assert positioner is None
    assert result.shape[1] == 2
    assert hdf5_win.btn_rxes_save.isEnabled()


def test_rxes_tab_with_positioner(hdf5_win):
    hdf5_win._rxes["scans"].setText("322")
    hdf5_win._rxes["positioner"].setText("4idgI0")
    hdf5_win._rxes_curv[0].setText("1")
    hdf5_win._rxes_curv[1].setText("2")
    hdf5_win._rxes_curv[2].setText("3")
    hdf5_win._on_rxes_run()

    positioner, (spectra, positioner_values) = hdf5_win._rxes_result
    assert positioner == "4idgI0"
    assert spectra.ndim == 3
    assert positioner_values.shape[0] == spectra.shape[0]


def test_rxes_tab_invalid_curvature_shows_dialog(hdf5_win):
    hdf5_win._rxes["scans"].setText("322")
    hdf5_win._rxes_curv[0].setText("abc")
    with patch(
        "polartools.process_images_gui.QMessageBox.critical"
    ) as mock_msg:
        hdf5_win._on_rxes_run()
    mock_msg.assert_called_once()


def test_rxes_save_1d_writes_file(win, tmp_path):
    win._rxes_result = (None, np.column_stack([np.arange(5), np.arange(5)]))
    path = str(tmp_path / "rxes.txt")
    with patch(
        "polartools.process_images_gui.QFileDialog.getSaveFileName",
        return_value=(path, ""),
    ):
        win._on_rxes_save()
    assert np.loadtxt(path).shape == (5, 2)


def test_rxes_save_2d_writes_file(win, tmp_path):
    spectra = np.random.rand(3, 4, 2)
    positioner_values = np.arange(3)
    win._rxes_result = ("pos", (spectra, positioner_values))
    path = str(tmp_path / "rxes.npz")
    with patch(
        "polartools.process_images_gui.QFileDialog.getSaveFileName",
        return_value=(path, ""),
    ):
        win._on_rxes_save()
    data = np.load(path)
    assert data["intensity"].shape == (3, 4)


# ─── RXES-MCD tab ─────────────────────────────────────────────────────────────


def test_mcd_copy_curvature_button(win):
    win._curvature = [1.0, 2.0, 3.0]
    from PyQt6.QtWidgets import QPushButton

    buttons = win.tabs.widget(2).findChildren(QPushButton)
    copy_btn = [b for b in buttons if b.text() == "Copy from Curvature tab"][0]
    copy_btn.click()
    assert [le.text() for le in win._mcd_curv] == ["1", "2", "3"]


def test_mcd_copy_curvature_without_data_warns(win):
    win._curvature = None
    from PyQt6.QtWidgets import QPushButton

    buttons = win.tabs.widget(2).findChildren(QPushButton)
    copy_btn = [b for b in buttons if b.text() == "Copy from Curvature tab"][0]
    with patch(
        "polartools.process_images_gui.QMessageBox.warning"
    ) as mock_warn:
        copy_btn.click()
    mock_warn.assert_called_once()


def test_mcd_tab_no_positioner_mocked(hdf5_win):
    hdf5_win._mcd["scans"].setText("322")
    hdf5_win._mcd_curv[0].setText("1")
    hdf5_win._mcd_curv[1].setText("2")
    hdf5_win._mcd_curv[2].setText("3")

    fake_rxes = np.column_stack([np.arange(5), np.arange(5) * 1.0])
    fake_mcd = np.column_stack([np.arange(5), np.arange(5) * 0.5])
    with patch(
        "polartools.process_images_gui.process_rxes_mcd",
        return_value=(fake_rxes, fake_mcd),
    ) as mock_proc:
        hdf5_win._on_mcd_run()

    mock_proc.assert_called_once()
    assert mock_proc.call_args.kwargs["folder"] == FOLDER
    is_map, result = hdf5_win._mcd_result
    assert not is_map
    assert hdf5_win.btn_mcd_save.isEnabled()


def test_mcd_tab_with_positioner_mocked(hdf5_win):
    hdf5_win._mcd["scans"].setText("322")
    hdf5_win._mcd["positioner"].setText("4idgI0")
    hdf5_win._mcd_curv[0].setText("1")
    hdf5_win._mcd_curv[1].setText("2")
    hdf5_win._mcd_curv[2].setText("3")

    rxes = np.random.rand(2, 4, 2)
    mcd = np.random.rand(2, 4, 2)
    positioner_values = np.arange(2)
    with patch(
        "polartools.process_images_gui.process_rxes_mcd",
        return_value=(rxes, mcd, positioner_values),
    ):
        hdf5_win._on_mcd_run()

    is_map, result = hdf5_win._mcd_result
    assert is_map


def test_mcd_save_1d_writes_file(win, tmp_path):
    rxes = np.column_stack([np.arange(5), np.arange(5) * 1.0])
    mcd = np.column_stack([np.arange(5), np.arange(5) * 0.5])
    win._mcd_result = (False, (rxes, mcd))
    path = str(tmp_path / "mcd.txt")
    with patch(
        "polartools.process_images_gui.QFileDialog.getSaveFileName",
        return_value=(path, ""),
    ):
        win._on_mcd_save()
    assert np.loadtxt(path).shape == (5, 3)


def test_mcd_save_2d_writes_file(win, tmp_path):
    rxes = np.random.rand(2, 4, 2)
    mcd = np.random.rand(2, 4, 2)
    positioner_values = np.arange(2)
    win._mcd_result = (True, (rxes, mcd, positioner_values))
    path = str(tmp_path / "mcd.npz")
    with patch(
        "polartools.process_images_gui.QFileDialog.getSaveFileName",
        return_value=(path, ""),
    ):
        win._on_mcd_save()
    data = np.load(path)
    assert data["rxes"].shape == (2, 4)
    assert data["mcd"].shape == (2, 4)
