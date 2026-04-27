"""
Copyright (c) 2020, UChicago Argonne, LLC.

See LICENSE file for details.
"""

import pytest

pytest.importorskip("PyQt6", reason="PyQt6 not installed — skipping GUI tests")
pytest.importorskip(
    "pyqtgraph", reason="pyqtgraph not installed — skipping GUI tests"
)


@pytest.fixture
def win(qtbot):
    from polartools.xmcd_gui import MainWindow

    window = MainWindow()
    qtbot.addWidget(window)
    return window


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


def test_norm_kwargs_manual_e0(win, qtbot):
    win.chk_e0_auto.setChecked(False)
    win.le_e0.setText("7112.5")
    kw = win._build_norm_kwargs()
    assert kw["e0"] == pytest.approx(7112.5)


def test_norm_kwargs_pre_range(win):
    win.le_pre1.setText("-150")
    win.le_pre2.setText("-30")
    kw = win._build_norm_kwargs()
    assert kw["pre_range"] == [-150.0, -30.0]


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


def test_resolve_load_kwargs_passes_extra(win):
    win.cb_kind.setCurrentIndex(0)
    kw = win._resolve_load_kwargs({"folder": "/tmp/data"})
    assert kw["folder"] == "/tmp/data"


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
