#!/usr/bin/env python
"""RXES image processing GUI using PyQt6 + matplotlib."""

import sys
import numpy as np

from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QComboBox,
    QSpinBox,
    QFileDialog,
    QGroupBox,
    QStatusBar,
    QMessageBox,
    QStackedWidget,
    QTabWidget,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from polartools.process_images import (
    load_images,
    get_curvature,
    process_rxes,
    process_rxes_mcd,
    _cleanup_photon_events,
)
from polartools._pyrixs import image_to_photon_events, plot_curvature


def _parse_scans(text):
    parts = [p.strip() for p in text.replace(";", ",").split(",") if p.strip()]
    if not parts:
        raise ValueError("No scan numbers provided.")
    return [int(p) for p in parts]


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RXES Image Processor")
        self.resize(1200, 800)

        self._curvature = None

        self._build_ui()

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

    # ─── UI construction ────────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)

        root.addWidget(self._build_source_panel())

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_curvature_tab(), "Curvature")
        self.tabs.addTab(self._build_rxes_tab(), "RXES")
        self.tabs.addTab(self._build_rxes_mcd_tab(), "RXES-MCD")
        root.addWidget(self.tabs, 1)

    def _browse_dir(self, line_edit):
        path = QFileDialog.getExistingDirectory(self, "Select folder")
        if path:
            line_edit.setText(path)

    def _build_source_panel(self):
        box = QGroupBox("Data source")
        layout = QHBoxLayout(box)

        layout.addWidget(QLabel("Source:"))
        self.cb_source = QComboBox()
        self.cb_source.addItems(["Databroker/Tiled catalog", "HDF5 folder"])
        layout.addWidget(self.cb_source)

        self._source_stack = QStackedWidget()

        cat_page = QWidget()
        h = QHBoxLayout(cat_page)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QLabel("Catalog:"))
        self.le_catalog = QLineEdit()
        self.le_catalog.setPlaceholderText("catalog-name")
        h.addWidget(self.le_catalog)
        self._source_stack.addWidget(cat_page)  # 0

        hdf5_page = QWidget()
        h = QHBoxLayout(hdf5_page)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QLabel("Folder:"))
        self.le_folder = QLineEdit()
        self.le_folder.setPlaceholderText("/path/to/hdf5/")
        h.addWidget(self.le_folder, 1)
        btn = QPushButton("Browse…")
        btn.clicked.connect(lambda: self._browse_dir(self.le_folder))
        h.addWidget(btn)
        self._source_stack.addWidget(hdf5_page)  # 1

        self.cb_source.currentIndexChanged.connect(
            self._source_stack.setCurrentIndex
        )
        layout.addWidget(self._source_stack, 1)

        layout.addWidget(QLabel("Detector key:"))
        self.le_detector_key = QLineEdit("lamb")
        self.le_detector_key.setMaximumWidth(100)
        layout.addWidget(self.le_detector_key)

        return box

    def _resolve_cat_kwargs(self):
        if self.cb_source.currentIndex() == 0:
            name = self.le_catalog.text().strip()
            if not name:
                raise ValueError("Catalog name is required.")
            from polartools.load_data import load_catalog

            return load_catalog(name), {}

        folder = self.le_folder.text().strip()
        if not folder:
            raise ValueError("HDF5 folder is required.")
        return "hdf5", {"folder": folder}

    def _build_common_inputs(self, layout, with_positioner=True):
        widgets = {}

        row = QHBoxLayout()
        row.addWidget(QLabel("Scans:"))
        le_scans = QLineEdit()
        le_scans.setPlaceholderText("1, 2, 3")
        row.addWidget(le_scans)
        widgets["scans"] = le_scans
        layout.addLayout(row)

        row = QHBoxLayout()
        row.addWidget(QLabel("Threshold (blank = none):"))
        le_threshold = QLineEdit()
        row.addWidget(le_threshold)
        widgets["threshold"] = le_threshold

        row.addWidget(QLabel("Normalize (blank = none):"))
        le_normalize = QLineEdit()
        row.addWidget(le_normalize)
        widgets["normalize"] = le_normalize
        layout.addLayout(row)

        if with_positioner:
            row = QHBoxLayout()
            row.addWidget(QLabel("Positioner (blank = none):"))
            le_positioner = QLineEdit()
            row.addWidget(le_positioner)
            widgets["positioner"] = le_positioner
            layout.addLayout(row)

        return widgets

    def _cleanup_from_text(self, text):
        text = text.strip()
        if not text:
            return None
        return dict(threshold=(float(text),))

    # ─── Curvature tab ──────────────────────────────────────────────────────

    def _build_curvature_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        self._curv = self._build_common_inputs(layout, with_positioner=False)

        row = QHBoxLayout()
        row.addWidget(QLabel("binx:"))
        self.sb_binx = QSpinBox()
        self.sb_binx.setRange(1, 10000)
        self.sb_binx.setValue(10)
        row.addWidget(self.sb_binx)

        row.addWidget(QLabel("biny:"))
        self.sb_curv_biny = QSpinBox()
        self.sb_curv_biny.setRange(1, 10000)
        self.sb_curv_biny.setValue(1)
        row.addWidget(self.sb_curv_biny)

        row.addWidget(QLabel("Constant offset (blank = auto):"))
        self.le_offset = QLineEdit()
        row.addWidget(self.le_offset)
        layout.addLayout(row)

        row = QHBoxLayout()
        self.btn_curv_run = QPushButton("Load && Fit")
        self.btn_curv_run.clicked.connect(self._on_curvature_run)
        row.addWidget(self.btn_curv_run)

        self.le_curv_result = QLineEdit()
        self.le_curv_result.setReadOnly(True)
        row.addWidget(self.le_curv_result, 1)

        self.btn_curv_save = QPushButton("Save curvature")
        self.btn_curv_save.setEnabled(False)
        self.btn_curv_save.clicked.connect(self._on_curvature_save)
        row.addWidget(self.btn_curv_save)
        layout.addLayout(row)

        self.curv_figure = Figure()
        self.curv_canvas = FigureCanvasQTAgg(self.curv_figure)
        layout.addWidget(self.curv_canvas, 1)

        return tab

    def _on_curvature_run(self):
        try:
            scans = _parse_scans(self._curv["scans"].text())
            cat, kwargs = self._resolve_cat_kwargs()
        except Exception as exc:
            QMessageBox.critical(self, "Input error", str(exc))
            return

        detector_key = self.le_detector_key.text().strip() or "lamb"
        cleanup = self._cleanup_from_text(self._curv["threshold"].text())
        normalize = self._curv["normalize"].text().strip() or None
        offset_text = self.le_offset.text().strip()
        constant_offset = int(offset_text) if offset_text else None

        self.status_bar.showMessage("Loading and fitting…")
        QApplication.processEvents()

        try:
            image = load_images(
                scans,
                cat,
                detector_key,
                cleanup=cleanup,
                normalize=normalize,
                **kwargs,
            )
            curvature = get_curvature(
                image,
                binx=self.sb_binx.value(),
                biny=self.sb_curv_biny.value(),
                constant_offset=constant_offset,
                plot=False,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Processing error", str(exc))
            self.status_bar.showMessage("Failed.")
            return

        self._curvature = curvature
        self.le_curv_result.setText(", ".join(f"{c:.6g}" for c in curvature))
        self.btn_curv_save.setEnabled(True)

        im = image.compute() if hasattr(image, "compute") else image
        ph = _cleanup_photon_events(image_to_photon_events(im.transpose()))
        vmax = np.nanpercentile(im, 99.99)

        self.curv_figure.clear()
        ax = self.curv_figure.add_subplot(111)
        mesh = ax.pcolor(
            im.transpose(),
            vmin=0,
            vmax=vmax if vmax > 0 else np.nanmax(im),
            cmap="plasma",
        )
        self.curv_figure.colorbar(mesh, ax=ax)
        plot_curvature(ax, curvature, ph)
        self.curv_canvas.draw()

        self.status_bar.showMessage("Done.")

    def _on_curvature_save(self):
        if self._curvature is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save curvature", "curvature.txt", "Text files (*.txt)"
        )
        if not path:
            return
        try:
            np.savetxt(path, self._curvature)
            self.status_bar.showMessage(f"Saved → {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Save error", str(exc))

    # ─── Shared curvature-row builder for RXES / RXES-MCD tabs ─────────────

    def _build_curvature_row(self, layout):
        row = QHBoxLayout()
        row.addWidget(QLabel("Curvature (c0, c1, c2):"))
        le_c0 = QLineEdit()
        le_c1 = QLineEdit()
        le_c2 = QLineEdit()
        for le in (le_c0, le_c1, le_c2):
            row.addWidget(le)
        btn_copy = QPushButton("Copy from Curvature tab")
        row.addWidget(btn_copy)
        layout.addLayout(row)

        def copy_curvature():
            if self._curvature is None:
                QMessageBox.warning(
                    self,
                    "No curvature",
                    "Run the Curvature tab first.",
                )
                return
            c0, c1, c2 = self._curvature
            le_c0.setText(f"{c0:.6g}")
            le_c1.setText(f"{c1:.6g}")
            le_c2.setText(f"{c2:.6g}")

        btn_copy.clicked.connect(copy_curvature)

        return le_c0, le_c1, le_c2

    def _read_curvature(self, fields):
        try:
            return [float(le.text()) for le in fields]
        except ValueError:
            raise ValueError("Curvature coefficients must be numeric.")

    # ─── RXES tab ───────────────────────────────────────────────────────────

    def _build_rxes_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        self._rxes = self._build_common_inputs(layout)
        self._rxes_curv = self._build_curvature_row(layout)

        row = QHBoxLayout()
        row.addWidget(QLabel("biny:"))
        self.sb_rxes_biny = QSpinBox()
        self.sb_rxes_biny.setRange(1, 10000)
        self.sb_rxes_biny.setValue(1)
        row.addWidget(self.sb_rxes_biny)

        self.btn_rxes_run = QPushButton("Process")
        self.btn_rxes_run.clicked.connect(self._on_rxes_run)
        row.addWidget(self.btn_rxes_run)

        self.btn_rxes_save = QPushButton("Save result")
        self.btn_rxes_save.setEnabled(False)
        self.btn_rxes_save.clicked.connect(self._on_rxes_save)
        row.addWidget(self.btn_rxes_save)
        layout.addLayout(row)

        self.rxes_figure = Figure()
        self.rxes_canvas = FigureCanvasQTAgg(self.rxes_figure)
        layout.addWidget(self.rxes_canvas, 1)

        self._rxes_result = None
        return tab

    def _on_rxes_run(self):
        try:
            scans = _parse_scans(self._rxes["scans"].text())
            cat, kwargs = self._resolve_cat_kwargs()
            curvature = self._read_curvature(self._rxes_curv)
        except Exception as exc:
            QMessageBox.critical(self, "Input error", str(exc))
            return

        detector_key = self.le_detector_key.text().strip() or "lamb"
        cleanup = self._cleanup_from_text(self._rxes["threshold"].text())
        normalize = self._rxes["normalize"].text().strip() or None
        positioner = self._rxes["positioner"].text().strip() or None
        biny = self.sb_rxes_biny.value()

        self.status_bar.showMessage("Processing…")
        QApplication.processEvents()

        try:
            result = process_rxes(
                scans,
                cat,
                detector_key,
                curvature,
                cleanup=cleanup,
                normalize=normalize,
                positioner=positioner,
                biny=biny,
                **kwargs,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Processing error", str(exc))
            self.status_bar.showMessage("Failed.")
            return

        self._rxes_result = (positioner, result)
        self.btn_rxes_save.setEnabled(True)

        self.rxes_figure.clear()
        ax = self.rxes_figure.add_subplot(111)
        if positioner is None:
            spectrum = result
            ax.plot(spectrum[:, 0], spectrum[:, 1])
        else:
            spectra, positioner_values = result
            mesh = ax.pcolormesh(
                positioner_values, spectra[0, :, 0], spectra[:, :, 1].T
            )
            self.rxes_figure.colorbar(mesh, ax=ax)
        self.rxes_canvas.draw()

        self.status_bar.showMessage("Done.")

    def _on_rxes_save(self):
        if self._rxes_result is None:
            return
        positioner, result = self._rxes_result
        if positioner is None:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save spectrum", "rxes.txt", "Text files (*.txt)"
            )
            if not path:
                return
            try:
                np.savetxt(path, result)
            except Exception as exc:
                QMessageBox.critical(self, "Save error", str(exc))
                return
        else:
            spectra, positioner_values = result
            path, _ = QFileDialog.getSaveFileName(
                self, "Save RIXS map", "rxes.npz", "NumPy archive (*.npz)"
            )
            if not path:
                return
            try:
                np.savez(
                    path,
                    pixel=spectra[0, :, 0],
                    positioner=positioner_values,
                    intensity=spectra[:, :, 1],
                )
            except Exception as exc:
                QMessageBox.critical(self, "Save error", str(exc))
                return
        self.status_bar.showMessage(f"Saved → {path}")

    # ─── RXES-MCD tab ───────────────────────────────────────────────────────

    def _build_rxes_mcd_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        self._mcd = self._build_common_inputs(layout)
        self._mcd_curv = self._build_curvature_row(layout)

        row = QHBoxLayout()
        row.addWidget(QLabel("biny:"))
        self.sb_mcd_biny = QSpinBox()
        self.sb_mcd_biny.setRange(1, 10000)
        self.sb_mcd_biny.setValue(1)
        row.addWidget(self.sb_mcd_biny)

        self.btn_mcd_run = QPushButton("Process")
        self.btn_mcd_run.clicked.connect(self._on_mcd_run)
        row.addWidget(self.btn_mcd_run)

        self.btn_mcd_save = QPushButton("Save result")
        self.btn_mcd_save.setEnabled(False)
        self.btn_mcd_save.clicked.connect(self._on_mcd_save)
        row.addWidget(self.btn_mcd_save)
        layout.addLayout(row)

        self.mcd_figure = Figure()
        self.mcd_canvas = FigureCanvasQTAgg(self.mcd_figure)
        layout.addWidget(self.mcd_canvas, 1)

        self._mcd_result = None
        return tab

    def _on_mcd_run(self):
        try:
            scans = _parse_scans(self._mcd["scans"].text())
            cat, kwargs = self._resolve_cat_kwargs()
            curvature = self._read_curvature(self._mcd_curv)
        except Exception as exc:
            QMessageBox.critical(self, "Input error", str(exc))
            return

        detector_key = self.le_detector_key.text().strip() or "lamb"
        cleanup = self._cleanup_from_text(self._mcd["threshold"].text())
        normalize = self._mcd["normalize"].text().strip() or None
        positioner = self._mcd["positioner"].text().strip() or None
        biny = self.sb_mcd_biny.value()

        self.status_bar.showMessage("Processing…")
        QApplication.processEvents()

        try:
            result = process_rxes_mcd(
                scans,
                cat,
                detector_key,
                curvature,
                cleanup=cleanup,
                normalize=normalize,
                positioner=positioner,
                biny=biny,
                **kwargs,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Processing error", str(exc))
            self.status_bar.showMessage("Failed.")
            return

        is_map = len(result) == 3
        self._mcd_result = (is_map, result)
        self.btn_mcd_save.setEnabled(True)

        self.mcd_figure.clear()
        ax_rxes = self.mcd_figure.add_subplot(121)
        ax_mcd = self.mcd_figure.add_subplot(122)
        ax_rxes.set_title("RXES")
        ax_mcd.set_title("MCD")

        if not is_map:
            rxes, mcd = result
            ax_rxes.plot(rxes[:, 0], rxes[:, 1])
            ax_mcd.plot(mcd[:, 0], mcd[:, 1])
        else:
            rxes, mcd, positioner_values = result
            mesh1 = ax_rxes.pcolormesh(
                positioner_values, rxes[0, :, 0], rxes[:, :, 1].T
            )
            self.mcd_figure.colorbar(mesh1, ax=ax_rxes)
            mesh2 = ax_mcd.pcolormesh(
                positioner_values, mcd[0, :, 0], mcd[:, :, 1].T
            )
            self.mcd_figure.colorbar(mesh2, ax=ax_mcd)
        self.mcd_canvas.draw()

        self.status_bar.showMessage("Done.")

    def _on_mcd_save(self):
        if self._mcd_result is None:
            return
        is_map, result = self._mcd_result
        if not is_map:
            rxes, mcd = result
            path, _ = QFileDialog.getSaveFileName(
                self, "Save RXES-MCD", "rxes_mcd.txt", "Text files (*.txt)"
            )
            if not path:
                return
            try:
                np.savetxt(
                    path,
                    np.column_stack([rxes[:, 0], rxes[:, 1], mcd[:, 1]]),
                )
            except Exception as exc:
                QMessageBox.critical(self, "Save error", str(exc))
                return
        else:
            rxes, mcd, positioner_values = result
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Save RXES-MCD map",
                "rxes_mcd.npz",
                "NumPy archive (*.npz)",
            )
            if not path:
                return
            try:
                np.savez(
                    path,
                    pixel=rxes[0, :, 0],
                    positioner=positioner_values,
                    rxes=rxes[:, :, 1],
                    mcd=mcd[:, :, 1],
                )
            except Exception as exc:
                QMessageBox.critical(self, "Save error", str(exc))
                return
        self.status_bar.showMessage(f"Saved → {path}")


def main():
    """Console entry point for the ``process-images-gui`` application."""
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
