#!/usr/bin/env python
"""XMCD processing GUI using PyQt6 + pyqtgraph."""

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
    QFrame,
    QMessageBox,
    QStackedWidget,
)
from PyQt6.QtCore import Qt, QSignalBlocker, QTimer
from PyQt6.QtGui import QFont

import pyqtgraph as pg
from pyqtgraph import mkPen

from polartools.absorption import (
    load_multi_dichro,
    load_multi_lockin,
    normalize_absorption,
    save_xmcd,
)

# ─── Color palette ────────────────────────────────────────────────────────────
C_PLUS_RAW = "#4C72B0"
C_MINUS_RAW = "#DD8452"
C_PRE = "#FF8800"
C_POST = "#00AA44"
C_E0 = "#222222"
C_PRE1 = "#FF6600"
C_PRE2 = "#CC4400"
C_POST1 = "#00AA44"
C_POST2 = "#006622"
C_PLUS_NORM = "#4C72B0"
C_MINUS_NORM = "#DD8452"
C_PLUS_XMCD = "#4C72B0"
C_MINUS_XMCD = "#DD8452"
C_MEAN_XMCD = "#5A3E8A"
C_ARTIFACT = "#888888"

NORM_DELAY_MS = 50


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("XMCD Processor")
        self.resize(1500, 900)

        self._plus_energy = None
        self._plus_mu = None
        self._plus_xmcd_raw = None
        self._minus_energy = None
        self._minus_mu = None
        self._minus_xmcd_raw = None
        self._plus_results = None
        self._minus_results = None
        self._e0_val = None
        self._block_line_update = False

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

        root.addLayout(self._build_save_bar())

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready — enter scan numbers and click Load & Process")

        self._init_plot_items()

    def _build_load_bar(self):
        bold = QFont()
        bold.setBold(True)

        layout = QVBoxLayout()
        layout.setSpacing(4)

        # Row 1: source selector
        row1 = QHBoxLayout()

        row1.addWidget(QLabel("Source:"))
        self.cb_source = QComboBox()
        self.cb_source.addItems(["SPEC", "HDF5", "CSV", "Databroker", "Tiled"])
        row1.addWidget(self.cb_source)

        # Stacked contextual source widgets (one per source type)
        self._source_stack = QStackedWidget()
        self._source_stack.addWidget(self._build_spec_source())   # 0
        self._source_stack.addWidget(self._build_hdf5_source())   # 1
        self._source_stack.addWidget(self._build_csv_source())    # 2
        self._source_stack.addWidget(self._build_db_source())     # 3
        self._source_stack.addWidget(self._build_tiled_source())  # 4
        row1.addWidget(self._source_stack, 1)
        layout.addLayout(row1)

        # Row 2: scan numbers + kind + load button
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("H+ scans:"))
        self.le_scans_plus = QLineEdit()
        self.le_scans_plus.setPlaceholderText("1, 2, 3")
        self.le_scans_plus.setToolTip("Comma-separated scan numbers for H+ field")
        row2.addWidget(self.le_scans_plus)

        row2.addWidget(QLabel("H− scans:"))
        self.le_scans_minus = QLineEdit()
        self.le_scans_minus.setPlaceholderText("4, 5, 6")
        self.le_scans_minus.setToolTip("Comma-separated scan numbers for H− field")
        row2.addWidget(self.le_scans_minus)

        row2.addWidget(QLabel("Kind:"))
        self.cb_kind = QComboBox()
        self.cb_kind.addItems(["dichro", "lockin"])
        row2.addWidget(self.cb_kind)

        # Stacked load-param widgets (dichro / lockin)
        self._kind_stack = QStackedWidget()
        self._kind_stack.addWidget(self._build_dichro_params())  # 0
        self._kind_stack.addWidget(self._build_lockin_params())  # 1
        row2.addWidget(self._kind_stack, 1)

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
                self.le_spec_path, "Open SPEC file", "SPEC files (*.dat *.txt);;All (*)"
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
        btn.clicked.connect(lambda: self._browse_dir(self.le_hdf_folder, "HDF5 folder"))
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
        btn.clicked.connect(lambda: self._browse_dir(self.le_csv_folder, "CSV folder"))
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

    # ── Load-parameter sub-widgets (dichro / lockin) ──────────────────────────

    def _build_dichro_params(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        for label, attr, tip in [
            ("Positioner:", "le_d_positioner", "Energy positioner name (blank = default)"),
            ("Detector:", "le_d_detector", "Detector name (blank = default IC5)"),
            ("Monitor:", "le_d_monitor", "Monitor name (blank = default IC4)"),
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
        self.chk_transmission.setToolTip("Transmission mode: ln(monitor/detector)")
        h.addWidget(self.chk_transmission)
        return w

    def _build_lockin_params(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        for label, attr, tip, placeholder in [
            ("Positioner:", "le_l_positioner", "Energy positioner (blank = default)", "default"),
            ("DC col:", "le_l_dc", "DC scaler column (blank = 'Lock DC')", "Lock DC"),
            ("AC col:", "le_l_ac", "AC scaler column (blank = 'Lock AC')", "Lock AC"),
            ("AC off:", "le_l_acoff", "AC offset column (blank = 'Lock AC off')", "Lock AC off"),
        ]:
            h.addWidget(QLabel(label))
            le = QLineEdit()
            le.setPlaceholderText(placeholder)
            le.setToolTip(tip)
            le.setMaximumWidth(120)
            setattr(self, attr, le)
            h.addWidget(le)
        return w

    # ── 6-plot grid ───────────────────────────────────────────────────────────

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

        self.pw_raw_plus = make_plot("Raw XANES  H+", "μ (a.u.)")
        self.pw_raw_minus = make_plot("Raw XANES  H−", "μ (a.u.)")
        self.pw_norm = make_plot("Normalized XANES", "μ (norm)")
        self.pw_xmcd_pm = make_plot("Normalized XMCD  H±", "XMCD (%)")
        self.pw_xmcd_mean = make_plot("Mean XMCD  (H+−H−)/2", "XMCD (%)")
        self.pw_artifact = make_plot("Artifact  (H++H−)/2", "Artifact (%)")

        # Link all X axes to pw_raw_plus
        for pw in (self.pw_raw_minus, self.pw_norm, self.pw_xmcd_pm,
                   self.pw_xmcd_mean, self.pw_artifact):
            pw.setXLink(self.pw_raw_plus)

        grid.addWidget(self.pw_raw_plus, 0, 0)
        grid.addWidget(self.pw_raw_minus, 0, 1)
        grid.addWidget(self.pw_norm, 1, 0)
        grid.addWidget(self.pw_xmcd_pm, 1, 1)
        grid.addWidget(self.pw_xmcd_mean, 2, 0)
        grid.addWidget(self.pw_artifact, 2, 1)

        return container

    def _init_plot_items(self):
        dash = Qt.PenStyle.DashLine
        dot = Qt.PenStyle.DotLine

        # Raw XANES plots
        self.curve_plus_raw = self.pw_raw_plus.plot(
            pen=mkPen(C_PLUS_RAW, width=2), name="H+ μ")
        self.curve_plus_pre = self.pw_raw_plus.plot(
            pen=mkPen(C_PRE, width=1.5, style=dash), name="pre-edge")
        self.curve_plus_post = self.pw_raw_plus.plot(
            pen=mkPen(C_POST, width=1.5, style=dash), name="post-edge")

        self.curve_minus_raw = self.pw_raw_minus.plot(
            pen=mkPen(C_MINUS_RAW, width=2), name="H− μ")
        self.curve_minus_pre = self.pw_raw_minus.plot(
            pen=mkPen(C_PRE, width=1.5, style=dash))
        self.curve_minus_post = self.pw_raw_minus.plot(
            pen=mkPen(C_POST, width=1.5, style=dash))

        # Draggable markers on H+ raw plot
        self.line_e0 = pg.InfiniteLine(
            angle=90, movable=False,
            pen=mkPen(C_E0, width=1.5, style=dot),
            label="e0", labelOpts={"position": 0.95})
        self.line_pre1 = pg.InfiniteLine(
            angle=90, movable=True,
            pen=mkPen(C_PRE1, width=1.5),
            label="pre1", labelOpts={"position": 0.90, "color": C_PRE1})
        self.line_pre2 = pg.InfiniteLine(
            angle=90, movable=True,
            pen=mkPen(C_PRE2, width=1.5),
            label="pre2", labelOpts={"position": 0.85, "color": C_PRE2})
        self.line_post1 = pg.InfiniteLine(
            angle=90, movable=True,
            pen=mkPen(C_POST1, width=1.5),
            label="post1", labelOpts={"position": 0.90, "color": C_POST1})
        self.line_post2 = pg.InfiniteLine(
            angle=90, movable=True,
            pen=mkPen(C_POST2, width=1.5),
            label="post2", labelOpts={"position": 0.85, "color": C_POST2})

        for line in (self.line_e0, self.line_pre1, self.line_pre2,
                     self.line_post1, self.line_post2):
            self.pw_raw_plus.addItem(line)
            line.setVisible(False)

        # Normalized XANES
        self.curve_norm_plus = self.pw_norm.plot(
            pen=mkPen(C_PLUS_NORM, width=2), name="H+")
        self.curve_norm_minus = self.pw_norm.plot(
            pen=mkPen(C_MINUS_NORM, width=2), name="H−")
        self.pw_norm.addLine(
            y=1.0, pen=mkPen("#cccccc", width=1, style=dash))

        # XMCD per polarity
        self.curve_xmcd_plus = self.pw_xmcd_pm.plot(
            pen=mkPen(C_PLUS_XMCD, width=2), name="H+")
        self.curve_xmcd_minus = self.pw_xmcd_pm.plot(
            pen=mkPen(C_MINUS_XMCD, width=2), name="H−")
        self.pw_xmcd_pm.addLine(y=0.0, pen=mkPen("#cccccc", width=1, style=dash))

        # Mean XMCD
        self.curve_xmcd_mean = self.pw_xmcd_mean.plot(
            pen=mkPen(C_MEAN_XMCD, width=2))
        self.pw_xmcd_mean.addLine(y=0.0, pen=mkPen("#cccccc", width=1, style=dash))

        # Artifact
        self.curve_artifact = self.pw_artifact.plot(
            pen=mkPen(C_ARTIFACT, width=2))
        self.pw_artifact.addLine(y=0.0, pen=mkPen("#cccccc", width=1, style=dash))

    # ── Parameter panel (normalization) ───────────────────────────────────────

    def _build_param_panel(self):
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(4, 4, 4, 4)
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

        layout.addStretch()
        return panel

    # ── Save bar ──────────────────────────────────────────────────────────────

    def _build_save_bar(self):
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Save as:"))
        self.le_savename = QLineEdit()
        self.le_savename.setPlaceholderText("xmcd_output.dat")
        self.le_savename.setMinimumWidth(280)
        bar.addWidget(self.le_savename)
        self.btn_save = QPushButton("Save")
        self.btn_save.setEnabled(False)
        bar.addWidget(self.btn_save)
        bar.addStretch()
        return bar

    # ─── Signal connections ────────────────────────────────────────────────────

    def _connect_signals(self):
        self.btn_load.clicked.connect(self._on_load)
        self.btn_save.clicked.connect(self._save_results)

        self.cb_source.currentIndexChanged.connect(self._source_stack.setCurrentIndex)
        self.cb_kind.currentIndexChanged.connect(self._kind_stack.setCurrentIndex)

        self.chk_e0_auto.toggled.connect(lambda c: self.le_e0.setEnabled(not c))
        self.chk_es_auto.toggled.connect(lambda c: self.le_edge_step.setEnabled(not c))
        self.le_e0.setEnabled(False)
        self.le_edge_step.setEnabled(False)

        # Draggable lines → entries
        self.line_pre1.sigPositionChanged.connect(
            lambda: self._line_moved(self.line_pre1, self.le_pre1))
        self.line_pre2.sigPositionChanged.connect(
            lambda: self._line_moved(self.line_pre2, self.le_pre2))
        self.line_post1.sigPositionChanged.connect(
            lambda: self._line_moved(self.line_post1, self.le_post1))
        self.line_post2.sigPositionChanged.connect(
            lambda: self._line_moved(self.line_post2, self.le_post2))

        # Entries → draggable lines
        self.le_pre1.editingFinished.connect(
            lambda: self._entry_changed(self.le_pre1, self.line_pre1))
        self.le_pre2.editingFinished.connect(
            lambda: self._entry_changed(self.le_pre2, self.line_pre2))
        self.le_post1.editingFinished.connect(
            lambda: self._entry_changed(self.le_post1, self.line_post1))
        self.le_post2.editingFinished.connect(
            lambda: self._entry_changed(self.le_post2, self.line_post2))

        for widget in (
            self.le_e0, self.le_edge_step, self.le_flat1, self.le_flat2,
            self.chk_e0_auto, self.chk_es_auto,
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
        if self._plus_energy is not None:
            self._norm_timer.start(NORM_DELAY_MS)

    # ─── Source resolution ─────────────────────────────────────────────────────

    def _resolve_source(self):
        """Return (source, extra_kwargs) based on current source widget state."""
        idx = self.cb_source.currentIndex()
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
            # source kwarg handled specially below
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
        """Return load kwargs dict for the selected xmcd_kind."""
        kwargs = dict(extra_kwargs)
        if self.cb_kind.currentText() == "dichro":
            pos = self.le_d_positioner.text().strip() or None
            det = self.le_d_detector.text().strip() or None
            mon = self.le_d_monitor.text().strip() or None
            if pos is not None:
                kwargs["positioner"] = pos
            if det is not None:
                kwargs["detector"] = det
            if mon is not None:
                kwargs["monitor"] = mon
            kwargs["transmission"] = self.chk_transmission.isChecked()
        else:
            pos = self.le_l_positioner.text().strip() or None
            dc = self.le_l_dc.text().strip() or None
            ac = self.le_l_ac.text().strip() or None
            acoff = self.le_l_acoff.text().strip() or None
            if pos is not None:
                kwargs["positioner"] = pos
            if dc is not None:
                kwargs["dc_col"] = dc
            if ac is not None:
                kwargs["ac_col"] = ac
            if acoff is not None:
                kwargs["acoff_col"] = acoff
        return kwargs

    # ─── Load phase ───────────────────────────────────────────────────────────

    def _parse_scan_list(self, text):
        parts = [p.strip() for p in text.replace(";", ",").split(",") if p.strip()]
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
        try:
            scans_plus = self._parse_scan_list(self.le_scans_plus.text())
            scans_minus = self._parse_scan_list(self.le_scans_minus.text())
        except ValueError as exc:
            QMessageBox.critical(self, "Input error", str(exc))
            return

        try:
            source, extra_kwargs = self._resolve_source()
        except Exception as exc:
            QMessageBox.critical(self, "Source error", str(exc))
            return

        load_kwargs = self._resolve_load_kwargs(extra_kwargs)
        load_func = (
            load_multi_dichro if self.cb_kind.currentText() == "dichro"
            else load_multi_lockin
        )

        self.status_bar.showMessage("Loading…")
        QApplication.processEvents()

        try:
            ep, yp, zp, _, _ = load_func(scans_plus, source, **load_kwargs)
            em, ym, zm, _, _ = load_func(scans_minus, source, **load_kwargs)
        except Exception as exc:
            QMessageBox.critical(self, "Load error", str(exc))
            self.status_bar.showMessage("Load failed.")
            return

        sort_p = np.argsort(ep)
        sort_m = np.argsort(em)
        self._plus_energy = ep[sort_p] * 1000
        self._plus_mu = yp[sort_p]
        self._plus_xmcd_raw = zp[sort_p]
        self._minus_energy = em[sort_m] * 1000
        self._minus_mu = ym[sort_m]
        self._minus_xmcd_raw = zm[sort_m]

        # Show raw data immediately
        self.curve_plus_raw.setData(self._plus_energy, self._plus_mu)
        self.curve_minus_raw.setData(self._minus_energy, self._minus_mu)

        # Set default normalization range markers on first load
        self._init_markers()

        self.status_bar.showMessage(
            f"Loaded {len(scans_plus)} H+ and {len(scans_minus)} H− scans — normalizing…"
        )
        self._run_normalization()

    def _init_markers(self):
        energy = self._plus_energy
        deriv = np.gradient(self._plus_mu, energy)
        e0_est = float(energy[np.argmax(deriv)])
        self._e0_val = e0_est

        for line in (self.line_e0, self.line_pre1, self.line_pre2,
                     self.line_post1, self.line_post2):
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
        e0 = None if self.chk_e0_auto.isChecked() else self._parse_entry(self.le_e0)
        edge_step = None if self.chk_es_auto.isChecked() else self._parse_entry(self.le_edge_step)

        pre1, pre2 = self._parse_entry(self.le_pre1), self._parse_entry(self.le_pre2)
        pre_range = [pre1, pre2] if (pre1 is not None or pre2 is not None) else None

        post1, post2 = self._parse_entry(self.le_post1), self._parse_entry(self.le_post2)
        post_range = [post1, post2] if (post1 is not None or post2 is not None) else None

        flat1, flat2 = self._parse_entry(self.le_flat1), self._parse_entry(self.le_flat2)
        flat_range = [flat1, flat2] if (flat1 is not None or flat2 is not None) else None

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
        if self._plus_energy is None:
            return

        norm_kw = self._build_norm_kwargs()
        try:
            plus = normalize_absorption(self._plus_energy, self._plus_mu, **norm_kw)
            minus = normalize_absorption(self._minus_energy, self._minus_mu, **norm_kw)
        except Exception as exc:
            self.status_bar.showMessage(f"Normalization error: {exc}")
            return

        plus["xmcd"] = self._plus_xmcd_raw / plus["edge_step"]
        minus["xmcd"] = self._minus_xmcd_raw / minus["edge_step"]

        self._plus_results = plus
        self._minus_results = minus

        # Update e0 marker
        self._e0_val = float(plus["e0"])
        self.line_e0.setPos(self._e0_val)
        with QSignalBlocker(self.le_e0):
            self.le_e0.setText(f"{self._e0_val:.2f}")

        ep, em = plus["energy"], minus["energy"]

        # Raw XANES with pre/post fits
        self.curve_plus_raw.setData(ep, plus["mu"])
        self.curve_plus_pre.setData(ep, plus["preedge"])
        self.curve_plus_post.setData(ep, plus["postedge"])
        self.curve_minus_raw.setData(em, minus["mu"])
        self.curve_minus_pre.setData(em, minus["preedge"])
        self.curve_minus_post.setData(em, minus["postedge"])

        # Normalized XANES
        self.curve_norm_plus.setData(ep, plus["norm"])
        self.curve_norm_minus.setData(em, minus["norm"])

        # XMCD per polarity
        self.curve_xmcd_plus.setData(ep, plus["xmcd"] * 100)
        self.curve_xmcd_minus.setData(em, minus["xmcd"] * 100)

        # Mean XMCD and artifact (interpolate minus onto plus energy grid)
        from scipy.interpolate import interp1d
        minus_xmcd_interp = interp1d(
            em, minus["xmcd"], bounds_error=False, fill_value=np.nan
        )(ep)
        self.curve_xmcd_mean.setData(
            ep, (plus["xmcd"] - minus_xmcd_interp) / 2 * 100
        )
        self.curve_artifact.setData(
            ep, (plus["xmcd"] + minus_xmcd_interp) / 2 * 100
        )

        self.btn_save.setEnabled(True)
        self.status_bar.showMessage(
            f"e0 = {plus['e0']:.2f} eV  |  "
            f"edge_step H+ = {plus['edge_step']:.4f}  |  "
            f"edge_step H− = {minus['edge_step']:.4f}"
        )

    # ─── Save ─────────────────────────────────────────────────────────────────

    def _save_results(self):
        if self._plus_results is None:
            return
        fname = self.le_savename.text().strip() or "xmcd_output.dat"
        if not os.path.isabs(fname):
            path, _ = QFileDialog.getSaveFileName(
                self, "Save XMCD output", fname,
                "Data files (*.dat *.txt);;All files (*)"
            )
            if not path:
                return
            fname = path
        try:
            save_xmcd(self._plus_results, self._minus_results, fname)
            self.status_bar.showMessage(f"Saved → {fname}")
        except Exception as exc:
            QMessageBox.critical(self, "Save error", str(exc))


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
