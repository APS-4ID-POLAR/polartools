"""Widgets and helpers shared by polartools.xmcd_gui and polartools.xanes_gui."""

from PyQt6.QtWidgets import (
    QFileDialog,
    QWidget,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
)


class GuiCommonMixin:
    """Source-widget builders and small parsing/line-sync helpers.

    Mixed into the XMCD and XANES GUI ``MainWindow`` classes, which both
    declare the same ``le_*``/``cb_*`` widget attribute names these methods
    build and read.
    """

    # ── Source parameter sub-widgets ────────────────────────────────────────

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

    # ─── Load phase ─────────────────────────────────────────────────────────

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

    # ─── Line ↔ entry synchronization ───────────────────────────────────────

    def _set_line_silent(self, line, pos):
        self._block_line_update = True
        line.setPos(pos)
        self._block_line_update = False

    # ─── Normalize phase ────────────────────────────────────────────────────

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
