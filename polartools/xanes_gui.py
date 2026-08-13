#!/usr/bin/env python
"""XANES normalization GUI using PyQt6 + pyqtgraph."""

import sys
import os
import numpy as np

from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QComboBox,
    QSpinBox,
    QCheckBox,
    QFileDialog,
    QGroupBox,
    QStatusBar,
    QSplitter,
    QMessageBox,
    QStackedWidget,
    QScrollArea,
)
from PyQt6.QtCore import Qt, QSignalBlocker, QTimer
from PyQt6.QtGui import QFont

import pyqtgraph as pg
from pyqtgraph import mkPen

from polartools.absorption import (
    load_multi_xas,
    normalize_absorption,
    save_xas,
)

# ─── Color palette ────────────────────────────────────────────────────────────
C_RAW = "#4C72B0"
C_PRE = "#FF8800"
C_POST = "#00AA44"
C_E0 = "#222222"
C_PRE1 = "#FF6600"
C_PRE2 = "#CC4400"
C_POST1 = "#00AA44"
C_POST2 = "#006622"
C_NORM = "#4C72B0"
C_FLAT = "#DD8452"
C_REF = "#888888"

NORM_DELAY_MS = 50


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("XANES Processor")
        self.resize(1300, 750)

        self._energy = None
        self._mu = None
        self._results = None
        self._e0_val = None
        self._block_line_update = False
        self._markers_initialized = False
        self._ref_energy = None
        self._ref_norm = None

        self._norm_timer = QTimer(singleShot=True)
        self._norm_timer.timeout.connect(self._run_normalization)

        self._build_ui()
        self._connect_signals()

    # ─── UI construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        pg.setConfigOption("background", "w")
        pg.setConfigOption("foreground", "k")

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 4)
        root.setSpacing(6)

        root.addLayout(self._build_load_bar())

        body = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(body, 1)

        body.addWidget(self._build_plots())
        body.addWidget(self._build_param_panel())
        body.setSizes([1000, 380])

        root.addLayout(self._build_reference_bar())
        root.addLayout(self._build_save_bar())

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage(
            "Ready — enter scan numbers and click Load & Process"
        )

        self._init_plot_items()

    def _build_load_bar(self):
        layout = QVBoxLayout()
        layout.setSpacing(4)

        # Row 1: source selector
        row1 = QHBoxLayout()

        row1.addWidget(QLabel("Source:"))
        self.cb_source = QComboBox()
        self.cb_source.addItems(
            ["SPEC", "HDF5", "CSV", "Databroker", "Tiled", "Column File"]
        )
        row1.addWidget(self.cb_source)

        # Stacked contextual source widgets (one per source type)
        self._source_stack = QStackedWidget()
        self._source_stack.addWidget(self._build_spec_source())  # 0
        self._source_stack.addWidget(self._build_hdf5_source())  # 1
        self._source_stack.addWidget(self._build_csv_source())  # 2
        self._source_stack.addWidget(self._build_db_source())  # 3
        self._source_stack.addWidget(self._build_tiled_source())  # 4
        self._source_stack.addWidget(self._build_column_source())  # 5
        row1.addWidget(self._source_stack, 1)
        layout.addLayout(row1)

        # Row 2: scan numbers + load params + load button
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Scans:"))
        self.le_scans = QLineEdit()
        self.le_scans.setPlaceholderText("1, 2, 3")
        self.le_scans.setToolTip("Comma-separated scan numbers")
        row2.addWidget(self.le_scans, 1)

        self._xas_params_widget = self._build_xas_params()
        row2.addWidget(self._xas_params_widget)

        self.btn_load = QPushButton("Load & Process")
        self.btn_load.setFont(QFont("", -1, QFont.Weight.Bold))
        self.btn_load.setMinimumWidth(130)
        row2.addWidget(self.btn_load)
        layout.addLayout(row2)

        return layout

    # ── Source parameter sub-widgets ──────────────────────────────────────────

    def _browse_file(self, line_edit, caption, filt):
        path, _ = QFileDialog.getOpenFileName(self, caption, "", filt)
        if path:
            line_edit.setText(path)

    def _browse_dir(self, line_edit, caption):
        path = QFileDialog.getExistingDirectory(self, caption)
        if path:
            line_edit.setText(path)

    def _build_spec_source(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QLabel("File:"))
        self.le_spec_path = QLineEdit()
        self.le_spec_path.setPlaceholderText("/path/to/spec.dat")
        h.addWidget(self.le_spec_path, 1)
        btn = QPushButton("Browse…")
        btn.clicked.connect(
            lambda: self._browse_file(
                self.le_spec_path,
                "Open SPEC file",
                "SPEC files (*.dat *.txt);;All (*)",
            )
        )
        h.addWidget(btn)
        h.addWidget(QLabel("Folder:"))
        self.le_spec_folder = QLineEdit()
        self.le_spec_folder.setPlaceholderText("(optional)")
        self.le_spec_folder.setMaximumWidth(160)
        h.addWidget(self.le_spec_folder)
        return w

    def _build_hdf5_source(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QLabel("Folder:"))
        self.le_hdf_folder = QLineEdit()
        self.le_hdf_folder.setPlaceholderText("/path/to/hdf5/")
        h.addWidget(self.le_hdf_folder, 1)
        btn = QPushButton("Browse…")
        btn.clicked.connect(
            lambda: self._browse_dir(self.le_hdf_folder, "HDF5 folder")
        )
        h.addWidget(btn)
        h.addWidget(QLabel("Format:"))
        self.le_hdf_format = QLineEdit("scan_{:06d}_master.hdf")
        self.le_hdf_format.setMaximumWidth(200)
        h.addWidget(self.le_hdf_format)
        h.addWidget(QLabel("H5 loc:"))
        self.le_hdf_loc = QLineEdit("entry/instrument/bluesky/streams/primary")
        self.le_hdf_loc.setMaximumWidth(280)
        h.addWidget(self.le_hdf_loc)
        return w

    def _build_csv_source(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QLabel("Folder:"))
        self.le_csv_folder = QLineEdit()
        self.le_csv_folder.setPlaceholderText("/path/to/csv/")
        h.addWidget(self.le_csv_folder, 1)
        btn = QPushButton("Browse…")
        btn.clicked.connect(
            lambda: self._browse_dir(self.le_csv_folder, "CSV folder")
        )
        h.addWidget(btn)
        return w

    def _build_db_source(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QLabel("Catalog:"))
        self.le_db_name = QLineEdit()
        self.le_db_name.setPlaceholderText("catalog-name")
        h.addWidget(self.le_db_name, 1)
        return w

    def _build_tiled_source(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QLabel("Profile:"))
        self.le_tiled_profile = QLineEdit()
        self.le_tiled_profile.setPlaceholderText("profile-name")
        h.addWidget(self.le_tiled_profile)
        h.addWidget(QLabel("Path:"))
        self.le_tiled_path = QLineEdit("/raw")
        self.le_tiled_path.setMaximumWidth(120)
        h.addWidget(self.le_tiled_path)
        return w

    def _build_column_source(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QLabel("File:"))
        self.le_column_path = QLineEdit()
        self.le_column_path.setPlaceholderText("/path/to/data.dat")
        self.le_column_path.editingFinished.connect(self._on_column_path_edited)
        h.addWidget(self.le_column_path, 1)
        btn = QPushButton("Browse…")
        btn.clicked.connect(self._browse_column_file)
        h.addWidget(btn)
        h.addWidget(QLabel("Energy col:"))
        self.cb_energy_col = QComboBox()
        self.cb_energy_col.addItem("0")
        h.addWidget(self.cb_energy_col)
        h.addWidget(QLabel("mu col:"))
        self.cb_mu_col = QComboBox()
        self.cb_mu_col.addItem("1")
        h.addWidget(self.cb_mu_col)
        return w

    def _browse_column_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open column file",
            "",
            "Data files (*.dat *.txt *.xy);;All files (*)",
        )
        if not path:
            return
        self.le_column_path.setText(path)
        self._preview_column_file(path)

    def _on_column_path_edited(self):
        path = self.le_column_path.text().strip()
        if path:
            self._preview_column_file(path)

    def _preview_column_file(self, path):
        """Populate the energy/mu column pickers from the file's columns."""
        try:
            data = np.loadtxt(path, comments="#")
        except Exception as exc:
            QMessageBox.critical(self, "Load error", str(exc))
            return
        if data.ndim == 1:
            data = data.reshape(-1, 1)
        ncols = data.shape[1]
        for cb in (self.cb_energy_col, self.cb_mu_col):
            with QSignalBlocker(cb):
                cb.clear()
                cb.addItems([str(i) for i in range(ncols)])
        with QSignalBlocker(self.cb_energy_col):
            self.cb_energy_col.setCurrentIndex(0)
        with QSignalBlocker(self.cb_mu_col):
            self.cb_mu_col.setCurrentIndex(min(1, ncols - 1))

    # ── Load-parameter widgets ────────────────────────────────────────────────

    def _build_xas_params(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        for label, attr, tip in [
            (
                "Positioner:",
                "le_positioner",
                "Energy positioner name (blank = default)",
            ),
            (
                "Detector:",
                "le_detector",
                "Detector name (blank = default IC5)",
            ),
            ("Monitor:", "le_monitor", "Monitor name (blank = default IC4)"),
        ]:
            h.addWidget(QLabel(label))
            le = QLineEdit()
            le.setPlaceholderText("default")
            le.setToolTip(tip)
            le.setMaximumWidth(120)
            setattr(self, attr, le)
            h.addWidget(le)
        self.chk_transmission = QCheckBox("Transmission")
        self.chk_transmission.setChecked(True)
        self.chk_transmission.setToolTip(
            "Transmission mode: ln(monitor/detector)"
        )
        h.addWidget(self.chk_transmission)
        return w

    # ── 2-plot grid ───────────────────────────────────────────────────────────

    def _build_plots(self):
        container = QWidget()
        grid = QGridLayout(container)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(4)

        def make_plot(title, ylabel):
            pw = pg.PlotWidget(title=title)
            pw.setLabel("left", ylabel)
            pw.setLabel("bottom", "Energy (eV)")
            pw.showGrid(x=True, y=True, alpha=0.3)
            pw.addLegend(offset=(10, 10))
            return pw

        self.pw_raw = make_plot("Raw XANES", "μ (a.u.)")
        self.pw_norm = make_plot("Normalized XANES", "μ (norm)")

        self.pw_norm.setXLink(self.pw_raw)

        grid.addWidget(self.pw_raw, 0, 0)
        grid.addWidget(self.pw_norm, 0, 1)

        return container

    def _init_plot_items(self):
        dash = Qt.PenStyle.DashLine
        dot = Qt.PenStyle.DotLine

        # Raw XANES
        self.curve_raw = self.pw_raw.plot(pen=mkPen(C_RAW, width=2), name="μ")
        self.curve_pre = self.pw_raw.plot(
            pen=mkPen(C_PRE, width=1.5, style=dash), name="pre-edge"
        )
        self.curve_post = self.pw_raw.plot(
            pen=mkPen(C_POST, width=1.5, style=dash), name="post-edge"
        )

        # Draggable markers on the raw plot
        self.line_e0 = pg.InfiniteLine(
            angle=90,
            movable=False,
            pen=mkPen(C_E0, width=1.5, style=dot),
            label="e0",
            labelOpts={"position": 0.95},
        )
        self.line_pre1 = pg.InfiniteLine(
            angle=90,
            movable=True,
            pen=mkPen(C_PRE1, width=1.5),
            label="pre1",
            labelOpts={"position": 0.90, "color": C_PRE1},
        )
        self.line_pre2 = pg.InfiniteLine(
            angle=90,
            movable=True,
            pen=mkPen(C_PRE2, width=1.5),
            label="pre2",
            labelOpts={"position": 0.85, "color": C_PRE2},
        )
        self.line_post1 = pg.InfiniteLine(
            angle=90,
            movable=True,
            pen=mkPen(C_POST1, width=1.5),
            label="post1",
            labelOpts={"position": 0.90, "color": C_POST1},
        )
        self.line_post2 = pg.InfiniteLine(
            angle=90,
            movable=True,
            pen=mkPen(C_POST2, width=1.5),
            label="post2",
            labelOpts={"position": 0.85, "color": C_POST2},
        )

        for line in (
            self.line_e0,
            self.line_pre1,
            self.line_pre2,
            self.line_post1,
            self.line_post2,
        ):
            self.pw_raw.addItem(line)
            line.setVisible(False)

        # Normalized XANES
        self.curve_norm = self.pw_norm.plot(
            pen=mkPen(C_NORM, width=2), name="norm"
        )
        self.curve_flat = self.pw_norm.plot(
            pen=mkPen(C_FLAT, width=2), name="flat"
        )
        self.curve_reference = self.pw_norm.plot(
            pen=mkPen(C_REF, width=1.5, style=dot), name="reference"
        )
        self.pw_norm.addLine(y=1.0, pen=mkPen("#cccccc", width=1, style=dash))

    # ── Parameter panel (normalization) ───────────────────────────────────────

    def _build_param_panel(self):
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(8)

        self.btn_guess = QPushButton("Guess parameters")
        self.btn_guess.setToolTip(
            "Re-guess e0 and the pre/post-edge ranges from the loaded data"
        )
        layout.addWidget(self.btn_guess)

        self._norm_panel, self._norm_boxes = self._build_norm_widgets()
        layout.addWidget(self._norm_panel)

        layout.addStretch()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(container)
        return scroll

    def _build_norm_widgets(self):
        """Build the 4 normalization GroupBoxes."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        bold = QFont()
        bold.setBold(True)

        def section(title):
            gb = QGroupBox(title)
            gb.setFont(bold)
            gl = QGridLayout(gb)
            gl.setVerticalSpacing(4)
            return gb, gl

        def entry(tip=""):
            w = QLineEdit()
            w.setPlaceholderText("auto")
            w.setMaximumWidth(90)
            if tip:
                w.setToolTip(tip)
            return w

        boxes = {}

        # Edge
        gb_edge, gl = section("Edge")
        gl.addWidget(QLabel("e0 (eV):"), 0, 0)
        self.le_e0 = entry("Absorption edge energy; blank = auto-detect")
        self.chk_e0_auto = QCheckBox("Auto")
        self.chk_e0_auto.setChecked(True)
        gl.addWidget(self.le_e0, 0, 1)
        gl.addWidget(self.chk_e0_auto, 0, 2)

        gl.addWidget(QLabel("edge_step:"), 1, 0)
        self.le_edge_step = entry("Edge step size; blank = auto")
        self.chk_es_auto = QCheckBox("Auto")
        self.chk_es_auto.setChecked(True)
        gl.addWidget(self.le_edge_step, 1, 1)
        gl.addWidget(self.chk_es_auto, 1, 2)
        layout.addWidget(gb_edge)
        boxes["edge"] = gb_edge

        # Pre-edge
        gb_pre, gl = section("Pre-edge  (drag orange lines)")
        gl.addWidget(QLabel("Start (rel. eV):"), 0, 0)
        self.le_pre1 = entry("Pre-edge start relative to e0")
        gl.addWidget(self.le_pre1, 0, 1)
        gl.addWidget(QLabel("End (rel. eV):"), 1, 0)
        self.le_pre2 = entry("Pre-edge end relative to e0")
        gl.addWidget(self.le_pre2, 1, 1)
        gl.addWidget(QLabel("Order:"), 2, 0)
        self.sp_pre_order = QSpinBox()
        self.sp_pre_order.setRange(0, 5)
        self.sp_pre_order.setValue(1)
        gl.addWidget(self.sp_pre_order, 2, 1)
        gl.addWidget(QLabel("nvict:"), 3, 0)
        self.sp_nvict = QSpinBox()
        self.sp_nvict.setRange(0, 3)
        self.sp_nvict.setValue(0)
        self.sp_nvict.setToolTip("Energy exponent for pre-edge fit")
        gl.addWidget(self.sp_nvict, 3, 1)
        layout.addWidget(gb_pre)
        boxes["pre"] = gb_pre

        # Post-edge
        gb_post, gl = section("Post-edge  (drag green lines)")
        gl.addWidget(QLabel("Start (rel. eV):"), 0, 0)
        self.le_post1 = entry("Post-edge start relative to e0")
        gl.addWidget(self.le_post1, 0, 1)
        gl.addWidget(QLabel("End (rel. eV):"), 1, 0)
        self.le_post2 = entry("Post-edge end relative to e0")
        gl.addWidget(self.le_post2, 1, 1)
        gl.addWidget(QLabel("Order:"), 2, 0)
        self.cb_post_order = QComboBox()
        self.cb_post_order.addItems(["Auto", "0", "1", "2", "3"])
        gl.addWidget(self.cb_post_order, 2, 1)
        layout.addWidget(gb_post)
        boxes["post"] = gb_post

        # Flatten
        gb_flat, gl = section("Flatten  (blank = same as post-edge)")
        gl.addWidget(QLabel("Start (rel. eV):"), 0, 0)
        self.le_flat1 = entry()
        gl.addWidget(self.le_flat1, 0, 1)
        gl.addWidget(QLabel("End (rel. eV):"), 1, 0)
        self.le_flat2 = entry()
        gl.addWidget(self.le_flat2, 1, 1)
        gl.addWidget(QLabel("Order:"), 2, 0)
        self.cb_flat_order = QComboBox()
        self.cb_flat_order.addItems(["Auto", "0", "1", "2", "3"])
        gl.addWidget(self.cb_flat_order, 2, 1)
        layout.addWidget(gb_flat)
        boxes["flat"] = gb_flat

        return panel, boxes

    # ── Reference overlay bar ─────────────────────────────────────────────────

    def _build_reference_bar(self):
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Reference:"))
        self.le_reference = QLineEdit()
        self.le_reference.setPlaceholderText(
            "(optional) previously saved XANES file"
        )
        self.le_reference.setMinimumWidth(280)
        bar.addWidget(self.le_reference, 1)
        btn_browse_ref = QPushButton("Browse…")
        btn_browse_ref.clicked.connect(
            lambda: self._browse_file(
                self.le_reference,
                "Reference file",
                "Data files (*.dat *.txt);;All files (*)",
            )
        )
        bar.addWidget(btn_browse_ref)
        self.btn_load_reference = QPushButton("Load reference")
        bar.addWidget(self.btn_load_reference)
        self.btn_clear_reference = QPushButton("Clear")
        self.btn_clear_reference.setEnabled(False)
        bar.addWidget(self.btn_clear_reference)
        return bar

    # ── Save bar ──────────────────────────────────────────────────────────────

    def _build_save_bar(self):
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Folder:"))
        self.le_savefolder = QLineEdit()
        self.le_savefolder.setPlaceholderText(
            "(current directory: " + os.getcwd() + ")"
        )
        self.le_savefolder.setMinimumWidth(280)
        bar.addWidget(self.le_savefolder, 1)
        btn_browse_folder = QPushButton("Browse…")
        btn_browse_folder.clicked.connect(
            lambda: self._browse_dir(self.le_savefolder, "Save folder")
        )
        bar.addWidget(btn_browse_folder)
        bar.addWidget(QLabel("Save as:"))
        self.le_savename = QLineEdit()
        self.le_savename.setPlaceholderText("xanes_output.dat")
        self.le_savename.setMinimumWidth(200)
        bar.addWidget(self.le_savename)
        self.btn_save = QPushButton("Save")
        self.btn_save.setEnabled(False)
        bar.addWidget(self.btn_save)
        return bar

    # ─── Signal connections ────────────────────────────────────────────────────

    def _connect_signals(self):
        self.btn_load.clicked.connect(self._on_load)
        self.btn_save.clicked.connect(self._save_results)
        self.btn_guess.clicked.connect(self._on_guess)
        self.btn_load_reference.clicked.connect(self._load_reference)
        self.btn_clear_reference.clicked.connect(self._clear_reference)

        self.cb_source.currentIndexChanged.connect(self._on_source_changed)

        self.chk_e0_auto.toggled.connect(lambda c: self.le_e0.setEnabled(not c))
        self.chk_es_auto.toggled.connect(
            lambda c: self.le_edge_step.setEnabled(not c)
        )
        self.le_e0.setEnabled(False)
        self.le_edge_step.setEnabled(False)

        # Draggable lines → entries
        self.line_pre1.sigPositionChanged.connect(
            lambda: self._line_moved(self.line_pre1, self.le_pre1)
        )
        self.line_pre2.sigPositionChanged.connect(
            lambda: self._line_moved(self.line_pre2, self.le_pre2)
        )
        self.line_post1.sigPositionChanged.connect(
            lambda: self._line_moved(self.line_post1, self.le_post1)
        )
        self.line_post2.sigPositionChanged.connect(
            lambda: self._line_moved(self.line_post2, self.le_post2)
        )

        # Entries → draggable lines
        self.le_pre1.editingFinished.connect(
            lambda: self._entry_changed(self.le_pre1, self.line_pre1)
        )
        self.le_pre2.editingFinished.connect(
            lambda: self._entry_changed(self.le_pre2, self.line_pre2)
        )
        self.le_post1.editingFinished.connect(
            lambda: self._entry_changed(self.le_post1, self.line_post1)
        )
        self.le_post2.editingFinished.connect(
            lambda: self._entry_changed(self.le_post2, self.line_post2)
        )

        for widget in (
            self.le_e0,
            self.le_edge_step,
            self.le_flat1,
            self.le_flat2,
            self.chk_e0_auto,
            self.chk_es_auto,
        ):
            if isinstance(widget, QCheckBox):
                widget.toggled.connect(self._schedule_normalize)
            else:
                widget.editingFinished.connect(self._schedule_normalize)

        self.sp_pre_order.valueChanged.connect(self._schedule_normalize)
        self.sp_nvict.valueChanged.connect(self._schedule_normalize)
        self.cb_post_order.currentIndexChanged.connect(self._schedule_normalize)
        self.cb_flat_order.currentIndexChanged.connect(self._schedule_normalize)

    def _schedule_normalize(self):
        if self._energy is not None:
            self._norm_timer.start(NORM_DELAY_MS)

    def _on_source_changed(self, index):
        self._source_stack.setCurrentIndex(index)
        is_column = self.cb_source.currentText() == "Column File"
        self.le_scans.setEnabled(not is_column)
        self._xas_params_widget.setEnabled(not is_column)

    # ─── Source resolution ─────────────────────────────────────────────────────

    def _resolve_source(self):
        """Return (source, extra_kwargs) based on current source widget state."""
        source_name = self.cb_source.currentText()

        if source_name == "SPEC":
            path = self.le_spec_path.text().strip()
            if not path:
                raise ValueError("SPEC file path is required.")
            folder = self.le_spec_folder.text().strip() or ""
            return path, {"folder": folder} if folder else {}

        if source_name == "HDF5":
            folder = self.le_hdf_folder.text().strip()
            if not folder:
                raise ValueError("HDF5 folder is required.")
            kwargs = {"source": "hdf5", "folder": folder}
            fmt = self.le_hdf_format.text().strip()
            if fmt:
                kwargs["fname_format"] = fmt
            loc = self.le_hdf_loc.text().strip()
            if loc:
                kwargs["h5_location"] = loc
            return "hdf5", {k: v for k, v in kwargs.items() if k != "source"}

        if source_name == "CSV":
            folder = self.le_csv_folder.text().strip()
            if not folder:
                raise ValueError("CSV folder is required.")
            return "csv", {"folder": folder}

        if source_name == "Databroker":
            from polartools.load_data import load_catalog

            name = self.le_db_name.text().strip()
            if not name:
                raise ValueError("Databroker catalog name is required.")
            cat = load_catalog(name)
            return cat, {}

        if source_name == "Tiled":
            from tiled.client import from_profile

            profile = self.le_tiled_profile.text().strip()
            if not profile:
                raise ValueError("Tiled profile name is required.")
            path = self.le_tiled_path.text().strip() or "/raw"
            cat = from_profile(profile)[path]
            return cat, {}

        raise ValueError(f"Unknown source: {source_name}")

    def _resolve_load_kwargs(self, extra_kwargs):
        """Return load kwargs dict for load_multi_xas."""
        kwargs = dict(extra_kwargs)
        pos = self.le_positioner.text().strip() or None
        det = self.le_detector.text().strip() or None
        mon = self.le_monitor.text().strip() or None
        if pos is not None:
            kwargs["positioner"] = pos
        if det is not None:
            kwargs["detector"] = det
        if mon is not None:
            kwargs["monitor"] = mon
        kwargs["transmission"] = self.chk_transmission.isChecked()
        return kwargs

    # ─── Load phase ───────────────────────────────────────────────────────────

    def _parse_scan_list(self, text):
        parts = [
            p.strip() for p in text.replace(";", ",").split(",") if p.strip()
        ]
        if not parts:
            raise ValueError("No scan numbers provided.")
        result = []
        for p in parts:
            try:
                result.append(int(p))
            except ValueError:
                result.append(p)
        return result

    def _on_load(self):
        if self.cb_source.currentText() == "Column File":
            self._on_load_column_file()
            return

        try:
            scans = self._parse_scan_list(self.le_scans.text())
        except ValueError as exc:
            QMessageBox.critical(self, "Input error", str(exc))
            return

        try:
            source, extra_kwargs = self._resolve_source()
        except Exception as exc:
            QMessageBox.critical(self, "Source error", str(exc))
            return

        load_kwargs = self._resolve_load_kwargs(extra_kwargs)

        self.status_bar.showMessage("Loading…")
        QApplication.processEvents()

        try:
            energy, mu, _ = load_multi_xas(scans, source, **load_kwargs)
        except Exception as exc:
            QMessageBox.critical(self, "Load error", str(exc))
            self.status_bar.showMessage("Load failed.")
            return

        sort = np.argsort(energy)
        self._energy = energy[sort] * 1000
        self._mu = mu[sort]

        self._finish_load(f"Loaded {len(scans)} scans — normalizing…")

    def _on_load_column_file(self):
        path = self.le_column_path.text().strip()
        if not path:
            QMessageBox.critical(self, "Input error", "Select a file to load.")
            return

        try:
            data = np.loadtxt(path, comments="#")
        except Exception as exc:
            QMessageBox.critical(self, "Load error", str(exc))
            self.status_bar.showMessage("Load failed.")
            return

        if data.ndim == 1:
            data = data.reshape(-1, 1)
        ncols = data.shape[1]
        ecol = self.cb_energy_col.currentIndex()
        mcol = self.cb_mu_col.currentIndex()
        if ecol < 0 or mcol < 0 or ecol >= ncols or mcol >= ncols:
            QMessageBox.critical(
                self, "Column error", "Selected column index is out of range."
            )
            self.status_bar.showMessage("Load failed.")
            return

        sort = np.argsort(data[:, ecol])
        self._energy = data[sort, ecol]
        self._mu = data[sort, mcol]

        self._finish_load(f"Loaded {os.path.basename(path)} — normalizing…")

    def _finish_load(self, status_msg):
        # Show raw data immediately
        self.curve_raw.setData(self._energy, self._mu)

        # Guess normalization range markers only on the first load. On later
        # loads, keep the current parameters (relative pre/post ranges) so
        # tuned settings survive; the marker lines are repositioned onto the
        # new e0 after normalization.
        first_load = not self._markers_initialized
        if first_load:
            self._init_markers()
            self._markers_initialized = True

        self.status_bar.showMessage(status_msg)
        self._run_normalization()

        if not first_load:
            self._sync_lines_to_entries()

    def _init_markers(self):
        energy = self._energy
        deriv = np.gradient(self._mu, energy)
        e0_est = float(energy[np.argmax(deriv)])
        self._e0_val = e0_est

        for line in (
            self.line_e0,
            self.line_pre1,
            self.line_pre2,
            self.line_post1,
            self.line_post2,
        ):
            line.setVisible(True)
        self.line_e0.setPos(e0_est)

        span = energy[-1] - energy[0]
        defaults = {
            self.line_pre1: -0.10 * span,
            self.line_pre2: -0.03 * span,
            self.line_post1: 0.05 * span,
            self.line_post2: 0.40 * span,
        }
        for line, rel in defaults.items():
            self._set_line_silent(line, e0_est + rel)

        for line, entry in [
            (self.line_pre1, self.le_pre1),
            (self.line_pre2, self.le_pre2),
            (self.line_post1, self.le_post1),
            (self.line_post2, self.le_post2),
        ]:
            entry.setText(f"{line.value() - e0_est:.1f}")

    def _on_guess(self):
        """Re-guess e0 and pre/post ranges from the loaded data on demand."""
        if self._energy is None:
            return
        self._init_markers()
        self._run_normalization()

    def _sync_lines_to_entries(self):
        """Move pre/post marker lines onto the new e0 after a reload."""
        for entry, line in [
            (self.le_pre1, self.line_pre1),
            (self.le_pre2, self.line_pre2),
            (self.le_post1, self.line_post1),
            (self.le_post2, self.line_post2),
        ]:
            rel = self._parse_entry(entry)
            if rel is not None and self._e0_val is not None:
                self._set_line_silent(line, self._e0_val + rel)

    # ─── Line ↔ entry synchronization ────────────────────────────────────────

    def _set_line_silent(self, line, pos):
        self._block_line_update = True
        line.setPos(pos)
        self._block_line_update = False

    def _line_moved(self, line, entry):
        if self._block_line_update or self._e0_val is None:
            return
        with QSignalBlocker(entry):
            entry.setText(f"{line.value() - self._e0_val:.1f}")
        self._schedule_normalize()

    def _entry_changed(self, entry, line):
        txt = entry.text().strip()
        if not txt:
            return
        try:
            rel = float(txt)
        except ValueError:
            return
        self._set_line_silent(line, (self._e0_val or 0.0) + rel)
        self._schedule_normalize()

    # ─── Normalize phase ──────────────────────────────────────────────────────

    def _parse_entry(self, entry):
        txt = entry.text().strip()
        if not txt:
            return None
        try:
            return float(txt)
        except ValueError:
            return None

    def _parse_order(self, combo):
        txt = combo.currentText()
        return None if txt == "Auto" else int(txt)

    def _build_norm_kwargs(self):
        e0 = (
            None
            if self.chk_e0_auto.isChecked()
            else self._parse_entry(self.le_e0)
        )
        edge_step = (
            None
            if self.chk_es_auto.isChecked()
            else self._parse_entry(self.le_edge_step)
        )

        pre1, pre2 = (
            self._parse_entry(self.le_pre1),
            self._parse_entry(self.le_pre2),
        )
        pre_range = (
            [pre1, pre2] if (pre1 is not None or pre2 is not None) else None
        )

        post1, post2 = (
            self._parse_entry(self.le_post1),
            self._parse_entry(self.le_post2),
        )
        post_range = (
            [post1, post2] if (post1 is not None or post2 is not None) else None
        )

        flat1, flat2 = (
            self._parse_entry(self.le_flat1),
            self._parse_entry(self.le_flat2),
        )
        flat_range = (
            [flat1, flat2] if (flat1 is not None or flat2 is not None) else None
        )

        return dict(
            e0=e0,
            edge_step=edge_step,
            pre_range=pre_range,
            pre_order=self.sp_pre_order.value(),
            nvict=self.sp_nvict.value(),
            post_range=post_range,
            post_order=self._parse_order(self.cb_post_order),
            flat_range=flat_range,
            flat_order=self._parse_order(self.cb_flat_order),
        )

    def _run_normalization(self):
        if self._energy is None:
            return

        norm_kw = self._build_norm_kwargs()
        try:
            results = normalize_absorption(self._energy, self._mu, **norm_kw)
        except Exception as exc:
            self.status_bar.showMessage(f"Normalization error: {exc}")
            return

        self._results = results

        # Update e0 marker
        self._e0_val = float(results["e0"])
        self.line_e0.setPos(self._e0_val)
        with QSignalBlocker(self.le_e0):
            self.le_e0.setText(f"{self._e0_val:.2f}")

        energy = results["energy"]

        # Raw XANES with pre/post fits
        self.curve_raw.setData(energy, results["mu"])
        self.curve_pre.setData(energy, results["preedge"])
        self.curve_post.setData(energy, results["postedge"])

        # Normalized XANES
        self.curve_norm.setData(energy, results["norm"])
        self.curve_flat.setData(energy, results["flat"])

        self.btn_save.setEnabled(True)
        self.status_bar.showMessage(
            f"e0 = {results['e0']:.2f} eV  |  "
            f"edge_step = {results['edge_step']:.4f}"
        )

    # ─── Reference overlay ──────────────────────────────────────────────────────

    def _load_reference(self):
        path = self.le_reference.text().strip()
        if not path:
            QMessageBox.warning(
                self, "Reference", "Select a reference file first."
            )
            return
        try:
            data = np.loadtxt(path)
        except Exception as exc:
            QMessageBox.critical(self, "Reference load error", str(exc))
            return
        if data.ndim != 2 or data.shape[1] < 3:
            QMessageBox.critical(
                self,
                "Reference load error",
                "Unexpected file format: expected columns "
                "Energy, XANES, Normalized[, Flattened].",
            )
            return
        self._ref_energy = data[:, 0]
        self._ref_norm = data[:, 2]
        self.curve_reference.setData(self._ref_energy, self._ref_norm)
        self.btn_clear_reference.setEnabled(True)
        self.status_bar.showMessage(f"Reference loaded from {path}")

    def _clear_reference(self):
        self._ref_energy = None
        self._ref_norm = None
        self.curve_reference.setData([], [])
        self.btn_clear_reference.setEnabled(False)

    # ─── Save ─────────────────────────────────────────────────────────────────

    def _resolve_save_path(self):
        """Return the target file path for saving, given folder/name entries.

        An explicit folder (typed or picked via Browse…) is always honored
        without prompting a dialog, even if it's a relative path — the user
        already made the choice. The save dialog only appears as a fallback
        when neither an absolute file name nor a folder was given.
        """
        fname = self.le_savename.text().strip() or "xanes_output.dat"
        if os.path.isabs(fname):
            return fname, True
        folder = self.le_savefolder.text().strip()
        if folder:
            return os.path.join(folder, fname), True
        return fname, False

    def _save_results(self):
        if self._results is None:
            return
        fname, resolved = self._resolve_save_path()
        if not resolved:
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Save XANES output",
                fname,
                "Data files (*.dat *.txt);;All files (*)",
            )
            if not path:
                return
            fname = path
        elif os.path.exists(fname):
            reply = QMessageBox.question(
                self,
                "Overwrite file?",
                f"{fname} already exists. Overwrite it?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                self.status_bar.showMessage("Save cancelled.")
                return
        try:
            save_xas(self._results, fname)
            self.status_bar.showMessage(f"Saved → {fname}")
        except Exception as exc:
            QMessageBox.critical(self, "Save error", str(exc))


def main():
    """Console entry point for the ``xanes-gui`` application."""
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
