"""Post-scan dichroism flyscan viewer.

Loads two scans (e.g. circular-right and circular-left polarizations),
displays their per-image scalar maps side-by-side, and computes a
``(CL - CR)`` difference map interpolated onto a common ``(x, y)`` grid.
Either polarization may be disabled, in which case only that map is drawn
and the difference panel is left blank -- so the viewer also covers the
single-scan case.

Eiger and Vortex each get a tab, with their own ROI and their own copy of
the three maps. The reduction itself -- position stream to per-frame x, y
and I0, detector file to per-frame ROI sum -- is
:mod:`polartools.process_flyscan`.

Settings live in a JSON file rather than in ``QSettings``: the folder,
ROIs, scan numbers and colour limits worth keeping are the ones for the
analysis directory being worked in, and a second experiment should not
inherit them. The file is read at startup -- ``./flyscan_gui_dichro.config.json`` unless
``--config`` names another -- but nothing is written until the Autosave
button in the status bar is pressed, so launching in a directory leaves
no trace in it. The destination can be changed while running; neither it
nor the Autosave state carries over to the next run, by design, so a
session never writes somewhere forgotten. Run with::

    flyscan-dichro-gui [--folder NdFeB] [--cr 490] [--cl 491] [--config FILE]
"""

import argparse
import json
import os
import re
import sys

import h5py
import hdf5plugin  # noqa: F401  -- registers HDF5 codecs (LZ4, Bitshuffle, ...)
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter, median_filter

# Before matplotlib, which picks its Qt binding from whatever is already
# imported -- see the longer note in image_viewer.py.
from PyQt6 import QtCore, QtWidgets

from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from .process_flyscan import (
    process_images,
    process_position_stream,
    process_vortex,
)
from .dichro_maps import interior_columns

#: How the difference colour scale is chosen. "Interior" exists because a
#: thick sample's absorption edges carry a far larger (and spurious)
#: difference than its interior, so a whole-map scale hides the domains; it
#: measures the middle third of the scan in x, which the alignment scan puts
#: inside the sample.
SCALE_MODES = ("Whole map", "Interior", "Symlog")
CMAPS = ("RdBu_r", "seismic", "bwr")

DEFAULT_EIGER_ROI = ((769, 832), (400, 463))  # ((x0, x1), (y0, y1))
DEFAULT_VORTEX_ROI = (800, 900)  # (channel_start, channel_stop)
PROCESS_BATCH = 100

CONFIG_FILENAME = "flyscan_gui_dichro.config.json"
FOLDER_HISTORY_MAX = 20

# Polarization keys; the difference is computed as POL_KEYS[1] - POL_KEYS[0].
POL_KEYS = ("cr", "cl")
POL_LABEL = {"cr": "CR", "cl": "CL"}

_FNAME_PATTERNS = (
    "{}/{}/scan_{:06d}.h5",
    "{}/{}/scan_{:d}.h5",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_h5(folder, stream, scan_number):
    for fmt in _FNAME_PATTERNS:
        path = fmt.format(folder, stream, scan_number)
        if os.path.exists(path):
            return path, fmt
    return None, None


def _parse_scan_list(text):
    """Parse a scan spec into a list of scan numbers.

    Accepts comma/whitespace-separated integers and inclusive ranges
    ``a-b`` (optionally with a step ``a-b:s``). Examples::

        "490"              -> [490]
        "490, 492, 494"    -> [490, 492, 494]
        "754-765"          -> [754, 755, ..., 765]
        "754-764:2"        -> [754, 756, 758, 760, 762, 764]
    """
    scans = []
    for token in re.split(r"[,\s]+", text.strip()):
        if not token:
            continue
        m = re.fullmatch(r"(\d+)-(\d+)(?::(\d+))?", token)
        if m:
            a = int(m.group(1))
            b = int(m.group(2))
            step = int(m.group(3)) if m.group(3) else 1
            if step <= 0:
                raise ValueError(f"step must be > 0 in '{token}'")
            scans.extend(range(a, b + 1, step))
        else:
            scans.append(int(token))
    return scans


def _offset_scan_spec(text, delta):
    """Shift every scan number in a spec by ``delta``, keeping its syntax.

    The field accepts more than a bare number, so stepping it cannot just add
    to an int -- the ranges and separators have to survive::

        "490"            -> "491"
        "490, 492, 494"  -> "491, 493, 495"
        "754-764:2"      -> "755-765:2"

    The step of a stepped range is a stride, not a scan number, so it is left
    alone. Returns ``None`` if the text is not a spec this can step, or if the
    shift would take any scan below zero -- callers leave the field untouched
    rather than mangle what the user typed.
    """
    if not text.strip():
        return None
    out = []
    for part in re.split(r"([,\s]+)", text.strip()):
        if not part or re.fullmatch(r"[,\s]+", part):
            out.append(part)
            continue
        m = re.fullmatch(r"(\d+)-(\d+)(?::(\d+))?", part)
        if m:
            a = int(m.group(1)) + delta
            b = int(m.group(2)) + delta
            if a < 0 or b < 0:
                return None
            out.append(f"{a}-{b}" + (f":{m.group(3)}" if m.group(3) else ""))
        elif re.fullmatch(r"\d+", part):
            v = int(part) + delta
            if v < 0:
                return None
            out.append(str(v))
        else:
            return None
    return "".join(out)


def _sum_scan_arrays(arrays):
    """Truncate a list of 1-D arrays to the shortest length and sum them."""
    n = min(a.size for a in arrays)
    return np.sum(
        np.stack([np.asarray(a[:n], dtype=np.float64) for a in arrays]), axis=0
    )


def _mean_scan_arrays(arrays):
    """Truncate a list of 1-D arrays to the shortest length and average them."""
    n = min(a.size for a in arrays)
    return np.mean(
        np.stack([np.asarray(a[:n], dtype=np.float64) for a in arrays]), axis=0
    )


def _read_eiger_frame(folder, scan_number, index=0):
    path, _ = _resolve_h5(folder, "eiger", scan_number)
    if path is None:
        return None
    with h5py.File(path, "r") as f:
        return f["entry/data/data"][index].astype(np.float32)


def _read_vortex_spectrum(folder, scan_number, index=0):
    path, _ = _resolve_h5(folder, "vortex", scan_number)
    if path is None:
        return None
    with h5py.File(path, "r") as f:
        frame = f["entry/data/data"][index]
    return frame.astype(np.float64).sum(axis=0)


def _get_frame_count(folder, stream, scan_number):
    path, _ = _resolve_h5(folder, stream, scan_number)
    if path is None:
        return 0
    try:
        with h5py.File(path, "r") as f:
            return int(f["entry/data/data"].shape[0])
    except Exception:
        return 0


def _default_config_path():
    return os.path.join(os.getcwd(), CONFIG_FILENAME)


def _load_config(path):
    try:
        with open(path, "r") as f:
            data = json.load(f)
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _save_config(cfg, path):
    """Write ``cfg`` to ``path``; return ``None``, or why it failed.

    Reported rather than swallowed, unlike the load side: the destination
    is one the user picked, so a directory they cannot write to has to say
    so instead of dropping every later change in silence.
    """
    try:
        with open(path, "w") as f:
            json.dump(cfg, f, indent=2)
    except Exception as e:
        return repr(e)
    return None


def _align_and_normalize(x, y, z, i0=None):
    """Trim ``x``/``y``/``z`` to a common length (dropping the first position
    sample) and optionally divide by ``i0`` element-wise. Used for the raw
    per-scan maps where no interpolation is needed.
    """
    xi = np.asarray(x[1:])
    yi = np.asarray(y[1:])
    zi = np.asarray(z, dtype=np.float64).copy()
    n = min(xi.size, yi.size, zi.size)
    xi, yi, zi = xi[:n], yi[:n], zi[:n]
    if i0 is not None and i0.size:
        i0i = np.asarray(i0[1 : 1 + n] if i0.size >= n + 1 else i0[:n])
        if i0i.size == n:
            with np.errstate(divide="ignore", invalid="ignore"):
                zi /= np.where(i0i != 0, i0i, np.nan)
    return xi, yi, zi


def _is_i0_usable(i0):
    if i0 is None or i0.size == 0:
        return False
    return bool(np.any(np.isfinite(i0)) and np.nanmax(np.abs(i0)) > 0)


def _trim_for_grid(xs, ys, sums, i0s, n_skip=3, normalize=True):
    """Trim per-image arrays for interpolation.

    Drops the first position sample (skipped trigger), the next
    ``n_skip`` aligned points (edge artifacts), and truncates to a common
    length in case image and position counts differ by more than one.
    Returns ``(xi, yi, zi)`` where ``zi`` is optionally normalized by I0,
    or ``(None, None, None)`` if there are too few samples left.
    """
    n = min(xs.size - 1, ys.size - 1, sums.size, i0s.size - 1) - n_skip
    if n <= 0:
        return None, None, None
    s = 1 + n_skip
    xi = np.asarray(xs[s : s + n])
    yi = np.asarray(ys[s : s + n])
    zi = np.asarray(sums[n_skip : n_skip + n], dtype=np.float64).copy()
    if normalize:
        i0i = np.asarray(i0s[s : s + n], dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            zi = zi / np.where(i0i != 0, i0i, np.nan)
    return xi, yi, zi


def _grid_pair(cr, cl, n_skip=3, grid_n=200, normalize=True):
    """Interpolate two scans' ``(xs, ys, sums, i0s)`` tuples onto a common grid.

    Returns ``(XG, YG, zcr, zcl, diff)``, or ``(None,)*5`` if either scan
    is empty or the (x, y) overlap is degenerate.
    """
    a = _trim_for_grid(*cr, n_skip=n_skip, normalize=normalize)
    b = _trim_for_grid(*cl, n_skip=n_skip, normalize=normalize)
    if a[0] is None or b[0] is None:
        return None, None, None, None, None
    xi_cr, yi_cr, zi_cr = a
    xi_cl, yi_cl, zi_cl = b
    xlo = max(xi_cr.min(), xi_cl.min())
    xhi = min(xi_cr.max(), xi_cl.max())
    ylo = max(yi_cr.min(), yi_cl.min())
    yhi = min(yi_cr.max(), yi_cl.max())
    if xhi <= xlo or yhi <= ylo:
        return None, None, None, None, None
    xg = np.linspace(xlo, xhi, grid_n)
    yg = np.linspace(ylo, yhi, grid_n)
    XG, YG = np.meshgrid(xg, yg)
    zcr_g = griddata((xi_cr, yi_cr), zi_cr, (XG, YG), method="linear")
    zcl_g = griddata((xi_cl, yi_cl), zi_cl, (XG, YG), method="linear")
    return XG, YG, zcr_g, zcl_g, zcl_g - zcr_g


# ---------------------------------------------------------------------------
# Background worker
# ---------------------------------------------------------------------------


class ProcessWorker(QtCore.QThread):
    finished_ok = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, func, args, kwargs):
        super().__init__()
        self._func = func
        self._args = args
        self._kwargs = kwargs

    def run(self):
        try:
            result = self._func(*self._args, **self._kwargs)
        except Exception as e:
            self.failed.emit(repr(e))
            return
        self.finished_ok.emit(result)


# ---------------------------------------------------------------------------
# Canvas widgets
# ---------------------------------------------------------------------------


class _Canvas(FigureCanvas):
    def __init__(self):
        fig = Figure(constrained_layout=True)
        super().__init__(fig)
        self.ax = fig.add_subplot(111)


class _MapCanvas(FigureCanvas):
    """Single-axis map canvas with a persistent title and one colorbar."""

    def __init__(self, title=""):
        fig = Figure(constrained_layout=True)
        super().__init__(fig)
        self.ax = fig.add_subplot(111)
        self._title = title
        self.ax.set_title(title)
        self.ax.set_aspect("equal")
        self._cbar = None

    def clear_axis(self, message=None):
        self.ax.clear()
        self.ax.set_title(self._title)
        self.ax.set_aspect("equal")
        if self._cbar is not None:
            try:
                self._cbar.remove()
            except Exception:
                # Matplotlib raises assorted types when the colorbar's axes
                # are already gone, which is the normal case on a redraw.
                # The point is to end up without one either way.
                pass
            self._cbar = None
        if message:
            self.ax.text(
                0.5,
                0.5,
                message,
                transform=self.ax.transAxes,
                ha="center",
                va="center",
                color="#888",
            )

    def set_colorbar(self, mappable):
        if self._cbar is not None:
            try:
                self._cbar.remove()
            except Exception:
                # As in clear_axis: dropping the old colorbar is best-effort,
                # and it is replaced on the next line regardless.
                pass
        self._cbar = self.figure.colorbar(mappable, ax=self.ax)


def _make_map_panel(canvas, parent, extra_controls=None):
    toolbar = NavigationToolbar(canvas, parent)
    box = QtWidgets.QWidget()
    v = QtWidgets.QVBoxLayout(box)
    v.setContentsMargins(0, 0, 0, 0)
    v.addWidget(toolbar)
    if extra_controls is not None:
        v.addWidget(extra_controls)
    v.addWidget(canvas, 1)
    return box


def _make_preview_panel(canvas, controls, parent, bottom=None):
    toolbar = NavigationToolbar(canvas, parent)
    box = QtWidgets.QWidget()
    v = QtWidgets.QVBoxLayout(box)
    v.setContentsMargins(0, 0, 0, 0)
    v.addWidget(toolbar)
    v.addWidget(controls)
    v.addWidget(canvas, 1)
    if bottom is not None:
        v.addWidget(bottom)
    return box


def _hline():
    """A thin sunken horizontal separator for use in the left-bar layout."""
    line = QtWidgets.QFrame()
    line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
    line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
    return line


def _make_frame_slider():
    """Slider + label for stepping through frames of a scan."""
    widget = QtWidgets.QWidget()
    h = QtWidgets.QHBoxLayout(widget)
    h.setContentsMargins(0, 0, 0, 0)
    label = QtWidgets.QLabel("Frame: 0 / 0")
    label.setMinimumWidth(90)
    slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
    slider.setMinimum(0)
    slider.setMaximum(0)
    slider.setValue(0)
    slider.setEnabled(False)
    h.addWidget(label)
    h.addWidget(slider, 1)
    return widget, slider, label


# ---------------------------------------------------------------------------
# Difference-map controls (shared by Eiger and Vortex tabs)
# ---------------------------------------------------------------------------


class _DiffControls(QtWidgets.QWidget):
    """Strip controlling interpolation and difference-map colour limits.

    Emits :attr:`changed` whenever any value changes; the host tab connects
    that to its redraw routine.
    """

    changed = QtCore.pyqtSignal()
    auto_requested = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.n_skip = QtWidgets.QSpinBox()
        self.n_skip.setRange(0, 10000)
        self.n_skip.setValue(3)
        self.grid_n = QtWidgets.QSpinBox()
        self.grid_n.setRange(20, 2000)
        self.grid_n.setValue(200)
        self.vmin = QtWidgets.QDoubleSpinBox()
        self.vmin.setDecimals(6)
        self.vmin.setRange(-1e12, 1e12)
        self.vmin.setValue(-1e-3)
        self.vmax = QtWidgets.QDoubleSpinBox()
        self.vmax.setDecimals(6)
        self.vmax.setRange(-1e12, 1e12)
        self.vmax.setValue(1e-3)
        # Cap the double-spinbox width — a ±1e12 range with 6 decimals auto-sizes
        # to ~170px each and makes the diff controls row force the top-right
        # grid cell wide. Realistic dichroism vmin/vmax fit in far less.
        self.vmin.setMaximumWidth(85)
        self.vmax.setMaximumWidth(85)
        self.auto_btn = QtWidgets.QPushButton("Auto")
        self.sym_chk = QtWidgets.QCheckBox("sym")
        self.sym_chk.setChecked(True)
        self.sym_chk.setToolTip(
            "Auto picks limits symmetric about zero.\n"
            "Uncheck to let it set vmin and vmax independently, which is what "
            "you want\nwhen the difference sits on an offset instead of being "
            "centred on zero."
        )

        # A thick sample absorbs so strongly that the difference at its edges
        # dwarfs the dichroism inside it, and a scale set from the whole map
        # renders the interior as flat white. These pick what the scale is set
        # from, and a colour map that changes fast enough near zero to show it.
        self.scale_combo = QtWidgets.QComboBox()
        self.scale_combo.addItems(SCALE_MODES)
        self.scale_combo.setToolTip(
            "Whole map: Auto uses percentiles of everything (the default).\n"
            "Interior: Auto uses only the middle third of the scan in x, so\n"
            "  weak domains fill the colour range and the absorption edges\n"
            "  saturate. The sample is centred by the alignment scan, so that\n"
            "  band is interior by construction.\n"
            "Symlog: linear near zero, logarithmic beyond -- shows interior\n"
            "  and edges at once. vmin/vmax do not apply."
        )
        self.cmap_combo = QtWidgets.QComboBox()
        self.cmap_combo.addItems(CMAPS)
        self.cmap_combo.setToolTip(
            "seismic changes colour about 3.6x faster near zero than RdBu_r\n"
            "(CIELAB dE per unit value), so weak signal is easier to see."
        )

        h = QtWidgets.QHBoxLayout(self)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QtWidgets.QLabel("n_skip:"))
        h.addWidget(self.n_skip)
        h.addSpacing(8)
        h.addWidget(QtWidgets.QLabel("grid:"))
        h.addWidget(self.grid_n)
        h.addSpacing(8)
        h.addWidget(QtWidgets.QLabel("vmin:"))
        h.addWidget(self.vmin)
        h.addWidget(QtWidgets.QLabel("vmax:"))
        h.addWidget(self.vmax)
        h.addWidget(self.auto_btn)
        h.addWidget(self.sym_chk)
        h.addSpacing(8)
        h.addWidget(QtWidgets.QLabel("scale:"))
        h.addWidget(self.scale_combo)
        h.addWidget(QtWidgets.QLabel("cmap:"))
        h.addWidget(self.cmap_combo)
        h.addStretch(1)

        for w in (self.n_skip, self.grid_n, self.vmin, self.vmax):
            w.valueChanged.connect(lambda *_: self.changed.emit())
        self.auto_btn.clicked.connect(self.auto_requested)
        # Toggling re-runs Auto, so the effect shows immediately rather than
        # waiting for the next time the button is pressed.
        self.sym_chk.toggled.connect(lambda *_: self.auto_requested.emit())
        self.cmap_combo.currentIndexChanged.connect(
            lambda *_: self.changed.emit()
        )
        self.scale_combo.currentIndexChanged.connect(self._scale_changed)
        self._sync_enabled()

    def _scale_changed(self):
        self._sync_enabled()
        # Symlog derives its own limits; the others need Auto re-run so the
        # change is visible straight away rather than at the next button press.
        if self.scale_mode() == "Symlog":
            self.changed.emit()
        else:
            self.auto_requested.emit()

    def _sync_enabled(self):
        """Grey out the manual limits when Symlog is driving the scale."""
        manual = self.scale_mode() != "Symlog"
        for w in (self.vmin, self.vmax, self.auto_btn, self.sym_chk):
            w.setEnabled(manual)

    def scale_mode(self):
        """Which data the colour scale is taken from."""
        return self.scale_combo.currentText()

    def cmap(self):
        return self.cmap_combo.currentText()

    def symmetric(self):
        """Whether Auto should force the limits symmetric about zero."""
        return self.sym_chk.isChecked()

    def set_vlim(self, vmin, vmax):
        self.vmin.blockSignals(True)
        self.vmax.blockSignals(True)
        try:
            self.vmin.setValue(float(vmin))
            self.vmax.setValue(float(vmax))
        finally:
            self.vmin.blockSignals(False)
            self.vmax.blockSignals(False)
        self.changed.emit()


class _DiffOptions(QtWidgets.QGroupBox):
    """Left-bar selector for interpolation / diff formula / smoothing."""

    changed = QtCore.pyqtSignal()

    INTERP_ITEMS = (
        ("Linear", "linear"),
        ("Cubic", "cubic"),
        ("Nearest", "nearest"),
    )
    DIFF_ITEMS = (
        ("CL − CR", "sub"),
        ("Normalized (2·Δ/Σ)", "norm_mean"),
        ("Asymmetry (Δ/Σ)", "asymmetry"),
    )
    SMOOTH_ITEMS = (
        ("None", "none"),
        ("Gaussian", "gaussian"),
        ("Median", "median"),
    )

    def __init__(self, parent=None):
        super().__init__("Diff options", parent)
        self.setStyleSheet("QGroupBox::title { font-weight: bold; }")

        self.interp_combo = QtWidgets.QComboBox()
        for label, _ in self.INTERP_ITEMS:
            self.interp_combo.addItem(label)

        self.diff_combo = QtWidgets.QComboBox()
        for label, _ in self.DIFF_ITEMS:
            self.diff_combo.addItem(label)

        self.smooth_combo = QtWidgets.QComboBox()
        for label, _ in self.SMOOTH_ITEMS:
            self.smooth_combo.addItem(label)

        self.param_spin = QtWidgets.QDoubleSpinBox()
        self.param_spin.setRange(0.1, 50.0)
        self.param_spin.setDecimals(1)
        self.param_spin.setSingleStep(0.5)
        self.param_spin.setValue(1.5)
        self.param_spin.setEnabled(False)

        form = QtWidgets.QFormLayout(self)
        form.setContentsMargins(6, 6, 6, 6)
        form.setVerticalSpacing(4)
        self._param_label = QtWidgets.QLabel("σ / kernel")
        form.addRow("Interp", self.interp_combo)
        form.addRow("Diff", self.diff_combo)
        form.addRow("Smooth", self.smooth_combo)
        form.addRow(self._param_label, self.param_spin)

        self.interp_combo.currentIndexChanged.connect(self._emit)
        self.diff_combo.currentIndexChanged.connect(self._emit)
        self.smooth_combo.currentIndexChanged.connect(self._on_smooth_changed)
        self.param_spin.valueChanged.connect(self._emit)

    def _emit(self, *_):
        self.changed.emit()

    def _on_smooth_changed(self, *_):
        mode = self.smooth_mode()
        self.param_spin.setEnabled(mode != "none")
        if mode == "median":
            self._param_label.setText("kernel")
        elif mode == "gaussian":
            self._param_label.setText("σ (px)")
        else:
            self._param_label.setText("σ / kernel")
        self.changed.emit()

    @staticmethod
    def _combo_key(items, combo):
        idx = combo.currentIndex()
        if 0 <= idx < len(items):
            return items[idx][1]
        return items[0][1]

    @staticmethod
    def _set_combo_key(items, combo, key):
        for i, (_, k) in enumerate(items):
            if k == key:
                combo.blockSignals(True)
                try:
                    combo.setCurrentIndex(i)
                finally:
                    combo.blockSignals(False)
                return

    def interp_method(self):
        return self._combo_key(self.INTERP_ITEMS, self.interp_combo)

    def diff_mode(self):
        return self._combo_key(self.DIFF_ITEMS, self.diff_combo)

    def smooth_mode(self):
        return self._combo_key(self.SMOOTH_ITEMS, self.smooth_combo)

    def smooth_param(self):
        return float(self.param_spin.value())

    def state(self):
        return {
            "interp": self.interp_method(),
            "diff_mode": self.diff_mode(),
            "smooth_mode": self.smooth_mode(),
            "smooth_param": self.smooth_param(),
        }

    def apply_state(self, state):
        """Restore from a saved dict; each key is optional."""
        if not isinstance(state, dict):
            return
        blocked = self.blockSignals(True)
        try:
            if "interp" in state:
                self._set_combo_key(
                    self.INTERP_ITEMS, self.interp_combo, str(state["interp"])
                )
            if "diff_mode" in state:
                self._set_combo_key(
                    self.DIFF_ITEMS, self.diff_combo, str(state["diff_mode"])
                )
            if "smooth_mode" in state:
                self._set_combo_key(
                    self.SMOOTH_ITEMS,
                    self.smooth_combo,
                    str(state["smooth_mode"]),
                )
            if "smooth_param" in state:
                try:
                    self.param_spin.setValue(float(state["smooth_param"]))
                except (TypeError, ValueError):
                    # Restored from a hand-editable JSON file, so the value
                    # may be junk. Keep the widget's default rather than
                    # refusing to open the window.
                    pass
            self._on_smooth_changed()
        finally:
            self.blockSignals(blocked)


class _PairVlim(QtWidgets.QGroupBox):
    """Shared vmin / vmax + Auto button for the CR and CL panels.

    When ``manual`` is unchecked (default) matplotlib auto-scales each panel
    just as before. When checked, the CR and CL panels use the spinbox limits
    so the two are directly comparable on the same colour scale.
    """

    changed = QtCore.pyqtSignal()
    auto_requested = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__("CR / CL v-limits", parent)
        self.setStyleSheet("QGroupBox::title { font-weight: bold; }")

        self.manual_chk = QtWidgets.QCheckBox("Manual")
        self.vmin = QtWidgets.QDoubleSpinBox()
        self.vmin.setDecimals(6)
        self.vmin.setRange(-1e12, 1e12)
        self.vmin.setValue(0.0)
        self.vmin.setMaximumWidth(90)
        self.vmax = QtWidgets.QDoubleSpinBox()
        self.vmax.setDecimals(6)
        self.vmax.setRange(-1e12, 1e12)
        self.vmax.setValue(1.0)
        self.vmax.setMaximumWidth(90)
        self.auto_btn = QtWidgets.QPushButton("Auto")

        self.vmin.setEnabled(False)
        self.vmax.setEnabled(False)

        form = QtWidgets.QFormLayout(self)
        form.setContentsMargins(6, 6, 6, 6)
        form.setVerticalSpacing(4)
        form.addRow(self.manual_chk)
        form.addRow("vmin", self.vmin)
        form.addRow("vmax", self.vmax)
        form.addRow(self.auto_btn)

        self.manual_chk.toggled.connect(self._on_manual_toggled)
        self.vmin.valueChanged.connect(self._emit_if_manual)
        self.vmax.valueChanged.connect(self._emit_if_manual)
        self.auto_btn.clicked.connect(self.auto_requested)

    def _on_manual_toggled(self, checked):
        self.vmin.setEnabled(checked)
        self.vmax.setEnabled(checked)
        self.changed.emit()

    def _emit_if_manual(self, *_):
        if self.manual_chk.isChecked():
            self.changed.emit()

    def is_manual(self):
        return self.manual_chk.isChecked()

    def limits(self):
        """Return ``(vmin, vmax)`` if manual is on, else ``(None, None)``."""
        if self.manual_chk.isChecked():
            return self.vmin.value(), self.vmax.value()
        return None, None

    def set_limits(self, vmin, vmax, activate=True):
        """Push computed values into the spinboxes.

        When ``activate`` is true (default) the manual checkbox is also
        turned on so the redraw picks up the new limits.
        """
        self.vmin.blockSignals(True)
        self.vmax.blockSignals(True)
        try:
            self.vmin.setValue(float(vmin))
            self.vmax.setValue(float(vmax))
        finally:
            self.vmin.blockSignals(False)
            self.vmax.blockSignals(False)
        if activate and not self.manual_chk.isChecked():
            self.manual_chk.setChecked(True)  # emits changed via toggle
        else:
            self.changed.emit()

    def state(self):
        return {
            "manual": self.manual_chk.isChecked(),
            "vmin": self.vmin.value(),
            "vmax": self.vmax.value(),
        }

    def apply_state(self, state):
        if not isinstance(state, dict):
            return
        blocked = self.blockSignals(True)
        try:
            if "vmin" in state:
                try:
                    self.vmin.setValue(float(state["vmin"]))
                except (TypeError, ValueError):
                    # Junk in the settings file; keep the current limit.
                    pass
            if "vmax" in state:
                try:
                    self.vmax.setValue(float(state["vmax"]))
                except (TypeError, ValueError):
                    # Junk in the settings file; keep the current limit.
                    pass
            if "manual" in state:
                checked = bool(state["manual"])
                self.manual_chk.setChecked(checked)
                self.vmin.setEnabled(checked)
                self.vmax.setEnabled(checked)
        finally:
            self.blockSignals(blocked)


# ---------------------------------------------------------------------------
# Map rendering
# ---------------------------------------------------------------------------


def _draw_pair(
    map_canvases,
    data,
    n_skip,
    normalize=True,
    scatter=True,
    vmin=None,
    vmax=None,
):
    """Render CR and CL scalar maps into their respective canvases.

    ``map_canvases`` is ``{'cr': _MapCanvas, 'cl': _MapCanvas}``. Missing
    ``data`` keys leave their panel blank. When ``vmin``/``vmax`` are given
    they override matplotlib's auto-scaling (same value for CR and CL so the
    two panels stay comparable). Returns the ``trimmed`` dict so the caller
    can reuse it for the difference panel.
    """
    trimmed = {}
    for key in POL_KEYS:
        entry = data.get(key)
        canvas = map_canvases[key]
        if entry is None:
            canvas.clear_axis(message=f"(no {POL_LABEL[key]} data)")
            continue
        xi, yi, zi = _trim_for_grid(*entry, n_skip=n_skip, normalize=normalize)
        if xi is None:
            canvas.clear_axis(message=f"(no {POL_LABEL[key]} after trim)")
            continue
        # tricontourf raises on non-finite z (e.g. NaN from a zero-I0 frame
        # when normalizing), which otherwise silently falls back to scatter
        # for whichever panel happens to hit it -- drop those points up
        # front so both panels honor the selected style the same way.
        finite = np.isfinite(xi) & np.isfinite(yi) & np.isfinite(zi)
        if not finite.all():
            xi, yi, zi = xi[finite], yi[finite], zi[finite]
        if xi.size == 0:
            canvas.clear_axis(message=f"(no finite {POL_LABEL[key]} samples)")
            continue
        trimmed[key] = (xi, yi, zi)
        canvas.clear_axis()
        if scatter:
            sc = canvas.ax.scatter(xi, yi, c=zi, s=8, vmin=vmin, vmax=vmax)
        else:
            # tricontourf with `levels=int` picks contour boundaries from the
            # data range, ignoring vmin/vmax for level placement. Build the
            # levels explicitly so the colorbar spans exactly [vmin, vmax]
            # when the user set manual limits; `extend='both'` keeps values
            # outside that range visible in the extremal colours.
            if vmin is not None and vmax is not None and vmax > vmin:
                contour_kwargs = dict(
                    levels=np.linspace(vmin, vmax, 51),
                    extend="both",
                )
            else:
                contour_kwargs = dict(levels=50)
            try:
                sc = canvas.ax.tricontourf(xi, yi, zi, **contour_kwargs)
            except Exception:
                sc = canvas.ax.scatter(
                    xi,
                    yi,
                    c=zi,
                    s=8,
                    vmin=vmin,
                    vmax=vmax,
                )
        canvas.ax.set_xlabel("x")
        canvas.ax.set_ylabel("y")
        canvas.set_colorbar(sc)

    if trimmed:
        all_x = np.concatenate([t[0] for t in trimmed.values()])
        all_y = np.concatenate([t[1] for t in trimmed.values()])
        for canvas in map_canvases.values():
            canvas.ax.set_xlim(all_x.min(), all_x.max())
            canvas.ax.set_ylim(all_y.min(), all_y.max())

    for canvas in map_canvases.values():
        canvas.draw_idle()
    return trimmed


def _combine_diff(zcr_g, zcl_g, diff_mode):
    """Combine gridded CR/CL according to ``diff_mode``."""
    if diff_mode == "norm_mean":
        with np.errstate(divide="ignore", invalid="ignore"):
            summed = zcl_g + zcr_g
            return np.where(summed != 0, 2.0 * (zcl_g - zcr_g) / summed, np.nan)
    if diff_mode == "asymmetry":
        with np.errstate(divide="ignore", invalid="ignore"):
            summed = zcl_g + zcr_g
            return np.where(summed != 0, (zcl_g - zcr_g) / summed, np.nan)
    return zcl_g - zcr_g


def _smooth_diff(diff, smooth_mode, smooth_param):
    """Apply the chosen filter to ``diff``, preserving NaN masking."""
    if smooth_mode == "gaussian":
        sigma = max(0.1, float(smooth_param))
        mask = (~np.isnan(diff)).astype(np.float64)
        filled = np.where(np.isnan(diff), 0.0, diff)
        num = gaussian_filter(filled, sigma=sigma, mode="nearest")
        denom = gaussian_filter(mask, sigma=sigma, mode="nearest")
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.where(denom > 1e-6, num / denom, np.nan)
        return out
    if smooth_mode == "median":
        # Force odd, >= 3
        k = max(3, int(round(smooth_param)))
        if k % 2 == 0:
            k += 1
        return median_filter(diff, size=k, mode="nearest")
    return diff


def _draw_diff(
    diff_canvas,
    trimmed,
    grid_n,
    vmin,
    vmax,
    interp_method="linear",
    diff_mode="sub",
    smooth_mode="none",
    smooth_param=1.5,
    cmap="RdBu_r",
    scale_mode="Whole map",
):
    """Render the CL - CR difference map from the trimmed CR/CL arrays.

    Returns ``(diff, mean, x)`` -- the gridded difference, the gridded mean
    intensity, and the grid's x axis. The mean is needed to find the absorbing
    interior when the colour scale is taken from it; all three are ``None``
    when there is nothing to draw.
    """
    diff = None
    mean = None
    xg = None
    if "cr" in trimmed and "cl" in trimmed:
        xi_cr, yi_cr, zi_cr = trimmed["cr"]
        xi_cl, yi_cl, zi_cl = trimmed["cl"]
        xlo = max(xi_cr.min(), xi_cl.min())
        xhi = min(xi_cr.max(), xi_cl.max())
        ylo = max(yi_cr.min(), yi_cl.min())
        yhi = min(yi_cr.max(), yi_cl.max())
        if xhi > xlo and yhi > ylo:
            xg = np.linspace(xlo, xhi, grid_n)
            yg = np.linspace(ylo, yhi, grid_n)
            XG, YG = np.meshgrid(xg, yg)
            xg = XG[0]
            zcr_g = griddata(
                (xi_cr, yi_cr), zi_cr, (XG, YG), method=interp_method
            )
            zcl_g = griddata(
                (xi_cl, yi_cl), zi_cl, (XG, YG), method=interp_method
            )
            diff = _combine_diff(zcr_g, zcl_g, diff_mode)
            mean = 0.5 * (zcr_g + zcl_g)
            if smooth_mode != "none":
                diff = _smooth_diff(diff, smooth_mode, smooth_param)
            diff_canvas.clear_axis()
            limits = _scale_kwargs(diff, scale_mode, vmin, vmax)
            pcm = diff_canvas.ax.pcolormesh(
                XG, YG, diff, shading="auto", cmap=cmap, **limits
            )
            diff_canvas.ax.set_xlabel("x")
            diff_canvas.ax.set_ylabel("y")
            diff_canvas.set_colorbar(pcm)
        else:
            diff_canvas.clear_axis(message="(no overlap)")
    else:
        diff_canvas.clear_axis(message="(need CR and CL)")

    if trimmed:
        all_x = np.concatenate([t[0] for t in trimmed.values()])
        all_y = np.concatenate([t[1] for t in trimmed.values()])
        diff_canvas.ax.set_xlim(all_x.min(), all_x.max())
        diff_canvas.ax.set_ylim(all_y.min(), all_y.max())

    diff_canvas.draw_idle()
    return diff, mean, xg


def _scale_kwargs(diff, scale_mode, vmin, vmax):
    """pcolormesh colour arguments for the selected scale mode.

    "Whole map" and "Interior" both end up as plain vmin/vmax -- they differ
    only in what Auto measured to produce them -- so only Symlog needs a norm
    of its own, derived from the data rather than the spin boxes.
    """
    if scale_mode != "Symlog":
        return dict(vmin=vmin, vmax=vmax)
    finite = diff[np.isfinite(diff)]
    if finite.size == 0:
        return dict(vmin=vmin, vmax=vmax)
    top = float(np.max(np.abs(finite)))
    if top <= 0:
        return dict(vmin=vmin, vmax=vmax)
    thresh = float(np.percentile(np.abs(finite), 99)) / 10.0
    if not np.isfinite(thresh) or thresh <= 0:
        thresh = top / 100.0
    return dict(
        norm=SymLogNorm(linthresh=thresh, linscale=0.6, vmin=-top, vmax=top)
    )


def _draw_maps(
    map_canvases,
    diff_canvas,
    data,
    n_skip,
    grid_n,
    vmin,
    vmax,
    normalize=True,
    scatter=True,
    interp_method="linear",
    diff_mode="sub",
    smooth_mode="none",
    smooth_param=1.5,
    pair_vmin=None,
    pair_vmax=None,
    cmap="RdBu_r",
    scale_mode="Whole map",
):
    trimmed = _draw_pair(
        map_canvases,
        data,
        n_skip=n_skip,
        normalize=normalize,
        scatter=scatter,
        vmin=pair_vmin,
        vmax=pair_vmax,
    )
    diff, mean, xg = _draw_diff(
        diff_canvas,
        trimmed,
        grid_n=grid_n,
        vmin=vmin,
        vmax=vmax,
        interp_method=interp_method,
        diff_mode=diff_mode,
        smooth_mode=smooth_mode,
        smooth_param=smooth_param,
        cmap=cmap,
        scale_mode=scale_mode,
    )
    return trimmed, diff, mean, xg


# ---------------------------------------------------------------------------
# Eiger tab
# ---------------------------------------------------------------------------


class EigerTab(QtWidgets.QWidget):
    process_requested = QtCore.pyqtSignal()
    redraw_requested = QtCore.pyqtSignal()
    frame_requested = QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._preview = None
        self._rect_patch = None
        self._reset_view_next = True
        self._last_diff = None
        self._last_mean = None
        self._last_grid_x = None
        self._last_trimmed = {}
        self._vlim_mode = None  # formula the colour limits were computed for
        self._current_frame = 0

        self.x_cen = QtWidgets.QSpinBox()
        self.x_cen.setMaximum(1_000_000)
        self.y_cen = QtWidgets.QSpinBox()
        self.y_cen.setMaximum(1_000_000)
        self.width = QtWidgets.QSpinBox()
        self.width.setMaximum(1_000_000)
        self.height = QtWidgets.QSpinBox()
        self.height.setMaximum(1_000_000)
        (dx0, dx1), (dy0, dy1) = DEFAULT_EIGER_ROI
        self.x_cen.setValue((dx0 + dx1) // 2)
        self.y_cen.setValue((dy0 + dy1) // 2)
        self.width.setValue(dx1 - dx0)
        self.height.setValue(dy1 - dy0)
        for s in (self.x_cen, self.y_cen, self.width, self.height):
            s.valueChanged.connect(self._update_overlay)

        self.process_btn = QtWidgets.QPushButton("Process")
        self.process_btn.clicked.connect(self.process_requested)

        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.addRow("x_cen", self.x_cen)
        form.addRow("y_cen", self.y_cen)
        form.addRow("width", self.width)
        form.addRow("height", self.height)
        form_widget = QtWidgets.QWidget()
        form_widget.setLayout(form)

        self.diff_options = _DiffOptions()
        self.diff_options.changed.connect(self.redraw_requested)

        self.pair_vlim = _PairVlim()
        self.pair_vlim.changed.connect(self.redraw_requested)
        self.pair_vlim.auto_requested.connect(self._auto_pair_vlim)

        left = QtWidgets.QWidget()
        left.setMaximumWidth(220)
        left_v = QtWidgets.QVBoxLayout(left)
        left_v.setContentsMargins(0, 0, 0, 0)
        left_v.addWidget(form_widget)
        left_v.addWidget(_hline())
        left_v.addWidget(self.process_btn)
        left_v.addWidget(_hline())
        left_v.addWidget(self.diff_options)
        left_v.addWidget(_hline())
        left_v.addWidget(self.pair_vlim)
        left_v.addStretch(1)

        self.scale_combo = QtWidgets.QComboBox()
        self.scale_combo.addItems(["Linear", "Log"])
        self.vmax_spin = QtWidgets.QDoubleSpinBox()
        self.vmax_spin.setRange(5.0, 100.0)
        self.vmax_spin.setDecimals(1)
        self.vmax_spin.setSingleStep(0.5)
        self.vmax_spin.setValue(99.5)
        self.vmax_spin.setSuffix(" %")
        self.scale_combo.currentIndexChanged.connect(self._render_preview)
        self.vmax_spin.valueChanged.connect(self._render_preview)

        ctrl = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(ctrl)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QtWidgets.QLabel("Scale:"))
        h.addWidget(self.scale_combo)
        h.addSpacing(12)
        h.addWidget(QtWidgets.QLabel("vmax:"))
        h.addWidget(self.vmax_spin)
        h.addStretch(1)

        self.preview_canvas = _Canvas()
        self.map_canvases = {
            "cr": _MapCanvas("CR"),
            "cl": _MapCanvas("CL"),
        }
        self.diff_canvas = _MapCanvas("CL - CR")
        self.diff_controls = _DiffControls()
        self.diff_controls.changed.connect(self.redraw_requested)
        self.diff_controls.auto_requested.connect(self._auto_vlim)

        slider_widget, self.frame_slider, self.frame_label = (
            _make_frame_slider()
        )
        self.frame_slider.valueChanged.connect(self._on_slider_changed)

        # Six columns, so the two rows can divide the width differently:
        #   TL Frame 0 preview (2) | TR CL - CR, with diff controls (4)
        #   BL CR (3)              | BR CL (3)
        grid = QtWidgets.QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(4)
        grid.addWidget(
            _make_preview_panel(
                self.preview_canvas,
                ctrl,
                self,
                bottom=slider_widget,
            ),
            0,
            0,
            1,
            2,
        )
        grid.addWidget(
            _make_map_panel(
                self.diff_canvas,
                self,
                self.diff_controls,
            ),
            0,
            2,
            1,
            4,
        )
        grid.addWidget(
            _make_map_panel(self.map_canvases["cr"], self), 1, 0, 1, 3
        )
        grid.addWidget(
            _make_map_panel(self.map_canvases["cl"], self), 1, 3, 1, 3
        )
        for i in range(2):
            grid.setRowStretch(i, 1)
        for i in range(6):
            grid.setColumnStretch(i, 1)
        grid_widget = QtWidgets.QWidget()
        grid_widget.setLayout(grid)

        outer = QtWidgets.QHBoxLayout(self)
        outer.addWidget(left)
        outer.addWidget(grid_widget, 1)

    def get_roi(self):
        xc, yc = self.x_cen.value(), self.y_cen.value()
        w, h = self.width.value(), self.height.value()
        return (
            (xc - w // 2, xc + (w - w // 2)),
            (yc - h // 2, yc + (h - h // 2)),
        )

    def set_preview(self, frame, frame_count=0):
        self._preview = frame
        self._current_frame = 0
        self._reset_view_next = True
        if frame is not None:
            ny, nx = frame.shape
            for spin, upper in (
                (self.x_cen, nx),
                (self.width, nx),
                (self.y_cen, ny),
                (self.height, ny),
            ):
                spin.setMaximum(upper)
        self._set_frame_range(frame_count)
        self._render_preview()

    def _set_frame_range(self, frame_count):
        max_index = max(0, frame_count - 1)
        self.frame_slider.blockSignals(True)
        try:
            self.frame_slider.setMaximum(max_index)
            self.frame_slider.setValue(0)
            self.frame_slider.setEnabled(frame_count > 1)
        finally:
            self.frame_slider.blockSignals(False)
        self.frame_label.setText(f"Frame: 0 / {max_index}")

    def _on_slider_changed(self, value):
        max_index = self.frame_slider.maximum()
        self.frame_label.setText(f"Frame: {value} / {max_index}")
        self.frame_requested.emit(value)

    def update_preview_frame(self, frame, index):
        """Refresh the preview image without resetting the current zoom."""
        self._preview = frame
        self._current_frame = index
        self._render_preview()

    def _render_preview(self):
        ax = self.preview_canvas.ax
        preserve = (
            (not self._reset_view_next)
            and ax.has_data()
            and self._preview is not None
        )
        if preserve:
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.clear()
        self._rect_patch = None
        frame = self._preview
        if frame is None:
            ax.set_title("(no Eiger data loaded)")
            self._reset_view_next = False
            self.preview_canvas.draw_idle()
            return
        pct = self.vmax_spin.value()
        vmax = float(np.percentile(frame, pct))
        if vmax <= 0:
            vmax = float(frame.max()) or 1.0
        if self.scale_combo.currentText() == "Log":
            positive = frame[frame > 0]
            vmin = float(positive.min()) if positive.size else 1.0
            if vmax <= vmin:
                vmax = vmin * 10.0
            ax.imshow(
                frame,
                origin="lower",
                cmap="viridis",
                norm=LogNorm(vmin=vmin, vmax=vmax),
            )
        else:
            ax.imshow(frame, origin="lower", vmin=0, vmax=vmax, cmap="viridis")
        ax.set_title(f"Eiger frame {self._current_frame}")
        self._rect_patch = Rectangle(
            (0, 0), 1, 1, edgecolor="red", facecolor="none", lw=1.5
        )
        ax.add_patch(self._rect_patch)
        if preserve:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        self._reset_view_next = False
        self._update_overlay()

    def _update_overlay(self):
        if self._preview is None or self._rect_patch is None:
            return
        (x0, x1), (y0, y1) = self.get_roi()
        self._rect_patch.set_xy((x0, y0))
        self._rect_patch.set_width(x1 - x0)
        self._rect_patch.set_height(y1 - y0)
        self.preview_canvas.draw_idle()

    def clear_maps(self):
        for canvas in self.map_canvases.values():
            canvas.clear_axis()
            canvas.draw_idle()
        self.diff_canvas.clear_axis()
        self.diff_canvas.draw_idle()
        self._last_diff = None
        self._last_mean = None
        self._last_grid_x = None
        self._last_trimmed = {}
        self._vlim_mode = None  # formula the colour limits were computed for

    def draw_maps(self, data, normalize=True, scatter=True):
        pair_vmin, pair_vmax = self.pair_vlim.limits()
        trimmed, diff, mean, grid_x = _draw_maps(
            self.map_canvases,
            self.diff_canvas,
            data,
            n_skip=self.diff_controls.n_skip.value(),
            grid_n=self.diff_controls.grid_n.value(),
            vmin=self.diff_controls.vmin.value(),
            vmax=self.diff_controls.vmax.value(),
            normalize=normalize,
            scatter=scatter,
            interp_method=self.diff_options.interp_method(),
            diff_mode=self.diff_options.diff_mode(),
            smooth_mode=self.diff_options.smooth_mode(),
            smooth_param=self.diff_options.smooth_param(),
            pair_vmin=pair_vmin,
            pair_vmax=pair_vmax,
            cmap=self.diff_controls.cmap(),
            scale_mode=self.diff_controls.scale_mode(),
        )
        self._last_trimmed = trimmed
        self._last_diff = diff
        self._last_mean = mean
        self._last_grid_x = grid_x

        # The colour limits belong to whichever formula produced them. An
        # asymmetry map runs to ~1e-2 where a raw CL-CR runs to ~1e1, so
        # carrying the old limits across a formula change renders the new map
        # as one flat colour -- which looks like the switch did nothing.
        # Rescale on the first draw and whenever the formula changes; set_vlim
        # emits `changed`, so the host redraws with the new limits.
        mode = self.diff_options.diff_mode()
        if diff is not None and mode != self._vlim_mode:
            self._vlim_mode = mode
            self._auto_vlim()

    def _interior_diff(self):
        """The difference over the central third of the scan in x.

        None when there is no grid to slice, or when the band leaves too
        little to measure -- callers then fall back to the whole map rather
        than scaling to a handful of pixels.
        """
        grid_x = getattr(self, "_last_grid_x", None)
        if grid_x is None or self._last_diff is None:
            return None
        mask = interior_columns(grid_x)
        if mask.sum() < 5:
            return None
        inner = self._last_diff[:, mask]
        return inner if np.isfinite(inner).sum() >= 10 else None

    def _auto_vlim(self):
        """Scale the difference colour bar to the data.

        Symmetric by default, which keeps zero in the middle of the diverging
        colour map. With `sym` unchecked the two ends are set independently
        from the 1st and 99th percentiles, for a difference that sits on an
        offset -- note that zero is then no longer the centre colour.
        """
        if self._last_diff is None:
            return
        source = self._last_diff
        if self.diff_controls.scale_mode() == "Interior":
            # Falls back to the whole map when no interior can be identified,
            # which is the right answer for a sample that does not block the
            # beam enough to have one.
            inner = self._interior_diff()
            if inner is not None:
                source = inner
        finite = source[np.isfinite(source)]
        if finite.size == 0:
            return
        if self.diff_controls.symmetric():
            vmax = float(np.nanpercentile(np.abs(finite), 98))
            if vmax == 0:
                vmax = float(np.nanmax(np.abs(finite))) or 1.0
            self.diff_controls.set_vlim(-vmax, vmax)
            return
        vmin = float(np.nanpercentile(finite, 1))
        vmax = float(np.nanpercentile(finite, 99))
        if vmax <= vmin:  # flat or near-flat map
            vmin = float(np.nanmin(finite))
            vmax = float(np.nanmax(finite))
        if vmax <= vmin:
            vmax = vmin + 1.0
        self.diff_controls.set_vlim(vmin, vmax)

    def _auto_pair_vlim(self):
        zs = [
            t[2]
            for t in self._last_trimmed.values()
            if t is not None and t[2] is not None
        ]
        if not zs:
            return
        finite = np.concatenate(zs)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return
        vmin = float(np.nanpercentile(finite, 2))
        vmax = float(np.nanpercentile(finite, 98))
        if vmax <= vmin:
            vmax = vmin + 1.0
        self.pair_vlim.set_limits(vmin, vmax)


# ---------------------------------------------------------------------------
# Vortex tab
# ---------------------------------------------------------------------------


class VortexTab(QtWidgets.QWidget):
    process_requested = QtCore.pyqtSignal()
    redraw_requested = QtCore.pyqtSignal()
    frame_requested = QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._preview = None
        self._span_patch = None
        self._reset_view_next = True
        self._last_diff = None
        self._last_mean = None
        self._last_grid_x = None
        self._last_trimmed = {}
        self._vlim_mode = None  # formula the colour limits were computed for
        self._current_frame = 0

        self.e_start = QtWidgets.QSpinBox()
        self.e_start.setMaximum(1_000_000)
        self.e_stop = QtWidgets.QSpinBox()
        self.e_stop.setMaximum(1_000_000)
        self.e_start.setValue(DEFAULT_VORTEX_ROI[0])
        self.e_stop.setValue(DEFAULT_VORTEX_ROI[1])
        for s in (self.e_start, self.e_stop):
            s.valueChanged.connect(self._update_overlay)

        self.process_btn = QtWidgets.QPushButton("Process")
        self.process_btn.clicked.connect(self.process_requested)

        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.addRow("e_start (channel)", self.e_start)
        form.addRow("e_stop (channel)", self.e_stop)
        form_widget = QtWidgets.QWidget()
        form_widget.setLayout(form)

        self.diff_options = _DiffOptions()
        self.diff_options.changed.connect(self.redraw_requested)

        self.pair_vlim = _PairVlim()
        self.pair_vlim.changed.connect(self.redraw_requested)
        self.pair_vlim.auto_requested.connect(self._auto_pair_vlim)

        left = QtWidgets.QWidget()
        left.setMaximumWidth(220)
        left_v = QtWidgets.QVBoxLayout(left)
        left_v.setContentsMargins(0, 0, 0, 0)
        left_v.addWidget(form_widget)
        left_v.addWidget(_hline())
        left_v.addWidget(self.process_btn)
        left_v.addWidget(_hline())
        left_v.addWidget(self.diff_options)
        left_v.addWidget(_hline())
        left_v.addWidget(self.pair_vlim)
        left_v.addStretch(1)

        self.scale_combo = QtWidgets.QComboBox()
        self.scale_combo.addItems(["Linear", "Log"])
        self.scale_combo.currentIndexChanged.connect(self._render_preview)

        ctrl = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(ctrl)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(QtWidgets.QLabel("Y scale:"))
        h.addWidget(self.scale_combo)
        h.addStretch(1)

        self.preview_canvas = _Canvas()
        self.map_canvases = {
            "cr": _MapCanvas("CR"),
            "cl": _MapCanvas("CL"),
        }
        self.diff_canvas = _MapCanvas("CL - CR")
        self.diff_controls = _DiffControls()
        self.diff_controls.changed.connect(self.redraw_requested)
        self.diff_controls.auto_requested.connect(self._auto_vlim)

        slider_widget, self.frame_slider, self.frame_label = (
            _make_frame_slider()
        )
        self.frame_slider.valueChanged.connect(self._on_slider_changed)

        # Six columns, so the two rows can divide the width differently:
        #   TL Frame 0 preview (2) | TR CL - CR, with diff controls (4)
        #   BL CR (3)              | BR CL (3)
        grid = QtWidgets.QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(4)
        grid.addWidget(
            _make_preview_panel(
                self.preview_canvas,
                ctrl,
                self,
                bottom=slider_widget,
            ),
            0,
            0,
            1,
            2,
        )
        grid.addWidget(
            _make_map_panel(
                self.diff_canvas,
                self,
                self.diff_controls,
            ),
            0,
            2,
            1,
            4,
        )
        grid.addWidget(
            _make_map_panel(self.map_canvases["cr"], self), 1, 0, 1, 3
        )
        grid.addWidget(
            _make_map_panel(self.map_canvases["cl"], self), 1, 3, 1, 3
        )
        for i in range(2):
            grid.setRowStretch(i, 1)
        for i in range(6):
            grid.setColumnStretch(i, 1)
        grid_widget = QtWidgets.QWidget()
        grid_widget.setLayout(grid)

        outer = QtWidgets.QHBoxLayout(self)
        outer.addWidget(left)
        outer.addWidget(grid_widget, 1)

    def get_roi(self):
        return (self.e_start.value(), self.e_stop.value())

    def set_preview(self, spectrum, frame_count=0):
        self._preview = spectrum
        self._current_frame = 0
        self._reset_view_next = True
        if spectrum is not None:
            nch = spectrum.size
            self.e_start.setMaximum(nch)
            self.e_stop.setMaximum(nch)
        self._set_frame_range(frame_count)
        self._render_preview()

    def _set_frame_range(self, frame_count):
        max_index = max(0, frame_count - 1)
        self.frame_slider.blockSignals(True)
        try:
            self.frame_slider.setMaximum(max_index)
            self.frame_slider.setValue(0)
            self.frame_slider.setEnabled(frame_count > 1)
        finally:
            self.frame_slider.blockSignals(False)
        self.frame_label.setText(f"Frame: 0 / {max_index}")

    def _on_slider_changed(self, value):
        max_index = self.frame_slider.maximum()
        self.frame_label.setText(f"Frame: {value} / {max_index}")
        self.frame_requested.emit(value)

    def update_preview_frame(self, spectrum, index):
        """Refresh the preview spectrum without resetting the current zoom."""
        self._preview = spectrum
        self._current_frame = index
        self._render_preview()

    def _render_preview(self):
        ax = self.preview_canvas.ax
        preserve = (
            (not self._reset_view_next)
            and ax.has_data()
            and self._preview is not None
        )
        if preserve:
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.clear()
        self._span_patch = None
        spectrum = self._preview
        if spectrum is None:
            ax.set_title("(no Vortex data loaded)")
            self._reset_view_next = False
            self.preview_canvas.draw_idle()
            return
        ax.plot(spectrum, lw=0.8)
        ax.set_xlabel("channel")
        ax.set_ylabel("counts (sum over elements)")
        ax.set_title(f"Vortex spectrum, frame {self._current_frame}")
        ax.set_yscale(
            "log" if self.scale_combo.currentText() == "Log" else "linear"
        )
        if preserve:
            ax.set_xlim(xlim)
            if ax.get_yscale() != "log" or ylim[0] > 0:
                ax.set_ylim(ylim)
        self._reset_view_next = False
        self._update_overlay()

    def _update_overlay(self):
        if self._preview is None:
            return
        ax = self.preview_canvas.ax
        if self._span_patch is not None:
            try:
                self._span_patch.remove()
            except Exception:
                # The patch may already have gone with a previous axes
                # clear. A fresh axvspan is drawn just below either way.
                pass
        e0, e1 = self.get_roi()
        self._span_patch = ax.axvspan(e0, e1, alpha=0.25, color="red")
        self.preview_canvas.draw_idle()

    def clear_maps(self):
        for canvas in self.map_canvases.values():
            canvas.clear_axis()
            canvas.draw_idle()
        self.diff_canvas.clear_axis()
        self.diff_canvas.draw_idle()
        self._last_diff = None
        self._last_mean = None
        self._last_grid_x = None
        self._last_trimmed = {}
        self._vlim_mode = None  # formula the colour limits were computed for

    def draw_maps(self, data, normalize=True, scatter=True):
        pair_vmin, pair_vmax = self.pair_vlim.limits()
        trimmed, diff, mean, grid_x = _draw_maps(
            self.map_canvases,
            self.diff_canvas,
            data,
            n_skip=self.diff_controls.n_skip.value(),
            grid_n=self.diff_controls.grid_n.value(),
            vmin=self.diff_controls.vmin.value(),
            vmax=self.diff_controls.vmax.value(),
            normalize=normalize,
            scatter=scatter,
            interp_method=self.diff_options.interp_method(),
            diff_mode=self.diff_options.diff_mode(),
            smooth_mode=self.diff_options.smooth_mode(),
            smooth_param=self.diff_options.smooth_param(),
            pair_vmin=pair_vmin,
            pair_vmax=pair_vmax,
            cmap=self.diff_controls.cmap(),
            scale_mode=self.diff_controls.scale_mode(),
        )
        self._last_trimmed = trimmed
        self._last_diff = diff
        self._last_mean = mean
        self._last_grid_x = grid_x

        # The colour limits belong to whichever formula produced them. An
        # asymmetry map runs to ~1e-2 where a raw CL-CR runs to ~1e1, so
        # carrying the old limits across a formula change renders the new map
        # as one flat colour -- which looks like the switch did nothing.
        # Rescale on the first draw and whenever the formula changes; set_vlim
        # emits `changed`, so the host redraws with the new limits.
        mode = self.diff_options.diff_mode()
        if diff is not None and mode != self._vlim_mode:
            self._vlim_mode = mode
            self._auto_vlim()

    def _interior_diff(self):
        """The difference over the central third of the scan in x.

        None when there is no grid to slice, or when the band leaves too
        little to measure -- callers then fall back to the whole map rather
        than scaling to a handful of pixels.
        """
        grid_x = getattr(self, "_last_grid_x", None)
        if grid_x is None or self._last_diff is None:
            return None
        mask = interior_columns(grid_x)
        if mask.sum() < 5:
            return None
        inner = self._last_diff[:, mask]
        return inner if np.isfinite(inner).sum() >= 10 else None

    def _auto_vlim(self):
        """Scale the difference colour bar to the data.

        Symmetric by default, which keeps zero in the middle of the diverging
        colour map. With `sym` unchecked the two ends are set independently
        from the 1st and 99th percentiles, for a difference that sits on an
        offset -- note that zero is then no longer the centre colour.
        """
        if self._last_diff is None:
            return
        source = self._last_diff
        if self.diff_controls.scale_mode() == "Interior":
            # Falls back to the whole map when no interior can be identified,
            # which is the right answer for a sample that does not block the
            # beam enough to have one.
            inner = self._interior_diff()
            if inner is not None:
                source = inner
        finite = source[np.isfinite(source)]
        if finite.size == 0:
            return
        if self.diff_controls.symmetric():
            vmax = float(np.nanpercentile(np.abs(finite), 98))
            if vmax == 0:
                vmax = float(np.nanmax(np.abs(finite))) or 1.0
            self.diff_controls.set_vlim(-vmax, vmax)
            return
        vmin = float(np.nanpercentile(finite, 1))
        vmax = float(np.nanpercentile(finite, 99))
        if vmax <= vmin:  # flat or near-flat map
            vmin = float(np.nanmin(finite))
            vmax = float(np.nanmax(finite))
        if vmax <= vmin:
            vmax = vmin + 1.0
        self.diff_controls.set_vlim(vmin, vmax)

    def _auto_pair_vlim(self):
        zs = [
            t[2]
            for t in self._last_trimmed.values()
            if t is not None and t[2] is not None
        ]
        if not zs:
            return
        finite = np.concatenate(zs)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return
        vmin = float(np.nanpercentile(finite, 2))
        vmax = float(np.nanpercentile(finite, 98))
        if vmax <= vmin:
            vmax = vmin + 1.0
        self.pair_vlim.set_limits(vmin, vmax)


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------


class MainWindow(QtWidgets.QMainWindow):
    def __init__(
        self,
        initial_folder=None,
        cli_cr=None,
        cli_cl=None,
        config_path=None,
    ):
        super().__init__()
        self.setWindowTitle("Flyscan dichroism viewer")
        # Roughly square content area so the 2x2 grid cells are square-ish
        # (~220px left form + 2 x ~530px cells wide; ~50px top bar +
        # 2 x ~530px cells + status bar tall).
        self.resize(1300, 1150)

        self._folder = ""
        # Per-polarization caches keyed by 'cr' and 'cl'. Each polarization
        # may aggregate several scan numbers (summed for better statistics).
        self._scans = {k: [] for k in POL_KEYS}  # list[int] per key
        self._xs = {
            k: None for k in POL_KEYS
        }  # per-image mean x (summed positions averaged)
        self._ys = {k: None for k in POL_KEYS}  # per-image mean y
        self._i0s = {
            k: None for k in POL_KEYS
        }  # per-image I0 summed across scans
        self._eiger_z = {
            k: None for k in POL_KEYS
        }  # per-image ROI sum summed across scans
        self._vortex_z = {k: None for k in POL_KEYS}
        # Per-scan HDF5 filename format so mixed 6-digit / unpadded scans
        # can be aggregated under the same polarization.
        self._fmt_eiger = {k: [] for k in POL_KEYS}
        self._fmt_vortex = {k: [] for k in POL_KEYS}
        self._workers = {}  # key -> ProcessWorker (keep refs alive)
        self._pending = set()  # in-flight worker keys for current operation
        self._eiger_preview_scan = None
        self._vortex_preview_scan = None
        self._folder_history = []
        self._loading_state = False  # guards auto-save while restoring config
        self._config_path = config_path or _default_config_path()

        self.folder_combo = QtWidgets.QComboBox()
        self.folder_combo.setMinimumWidth(320)
        self.folder_combo.setMaxVisibleItems(20)
        self.folder_combo.setSizeAdjustPolicy(
            QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToContents
        )
        self.folder_combo.currentTextChanged.connect(self._on_folder_selected)
        browse_btn = QtWidgets.QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_folder)

        # CR / CL scan inputs with enable checkboxes. The scan field accepts
        # one or more scan numbers (comma/space separated, plus ranges
        # ``a-b`` and stepped ranges ``a-b:s``); the loader sums them for
        # better statistics.
        self.enable_chk = {}
        self.scan_input = {}
        self.scan_steppers = {}
        # Guards the CR -> CL auto-fill against re-entering itself.
        self._syncing_scans = False
        for k in POL_KEYS:
            chk = QtWidgets.QCheckBox(POL_LABEL[k])
            chk.setChecked(True)
            le = QtWidgets.QLineEdit()
            le.setText("490" if k == "cr" else "491")
            le.setPlaceholderText("e.g. 490 or 754-764:2")
            le.setMinimumWidth(160)

            # A spinbox can't hold this field's syntax, so the steppers are a
            # pair of arrow buttons that add or subtract one from every scan
            # number in whatever is typed. Auto-repeat lets a held button run
            # through a series.
            up = QtWidgets.QToolButton()
            up.setArrowType(QtCore.Qt.ArrowType.UpArrow)
            down = QtWidgets.QToolButton()
            down.setArrowType(QtCore.Qt.ArrowType.DownArrow)
            for btn, step in ((up, 1), (down, -1)):
                btn.setFixedSize(18, 12)
                btn.setAutoRepeat(True)
                btn.setAutoRepeatDelay(400)
                btn.setAutoRepeatInterval(120)
                btn.setToolTip(
                    f"{'Increase' if step > 0 else 'Decrease'}"
                    f" every scan number by 1"
                )
                btn.clicked.connect(
                    lambda _checked=False, key=k, d=step: self._step_scans(
                        key, d
                    )
                )

            self.enable_chk[k] = chk
            self.scan_input[k] = le
            self.scan_steppers[k] = (up, down)
            chk.toggled.connect(le.setEnabled)
            chk.toggled.connect(up.setEnabled)
            chk.toggled.connect(down.setEnabled)

        # Typing in CR carries CL along with it, but only while both are on --
        # unticking CL is how you break the link and drive it by hand.
        self.scan_input["cr"].textEdited.connect(self._sync_cl_from_cr)
        self.enable_chk["cl"].toggled.connect(
            lambda on: self._sync_cl_from_cr() if on else None
        )

        load_btn = QtWidgets.QPushButton("Load")
        load_btn.clicked.connect(self._on_load)

        self.norm_chk = QtWidgets.QCheckBox("Normalize by I0")
        self.norm_chk.setEnabled(False)
        self.norm_chk.stateChanged.connect(self._redraw_current_tab)

        self.style_combo = QtWidgets.QComboBox()
        self.style_combo.addItems(["Scatter", "Tricontourf"])
        self.style_combo.currentIndexChanged.connect(self._redraw_current_tab)

        # CR/CL scan inputs grouped under a shared format hint so the
        # accepted syntaxes are visible without hovering the placeholder.
        scan_group = QtWidgets.QWidget()
        sg_v = QtWidgets.QVBoxLayout(scan_group)
        sg_v.setContentsMargins(0, 0, 0, 0)
        sg_v.setSpacing(2)
        hint = QtWidgets.QLabel(
            "Scans: 490 · 490,492,494 · 754-765 · 754-764:2 (start-end:step)"
        )
        hint.setStyleSheet("QLabel { color: #666; font-size: 10px; }")
        sg_v.addWidget(hint)
        sg_h = QtWidgets.QHBoxLayout()
        sg_h.setContentsMargins(0, 0, 0, 0)
        for k in POL_KEYS:
            sg_h.addWidget(self.enable_chk[k])
            sg_h.addWidget(self.scan_input[k])
            up, down = self.scan_steppers[k]
            arrows = QtWidgets.QVBoxLayout()
            arrows.setContentsMargins(0, 0, 0, 0)
            arrows.setSpacing(0)
            arrows.addWidget(up)
            arrows.addWidget(down)
            sg_h.addLayout(arrows)
            sg_h.addSpacing(6)
        sg_v.addLayout(sg_h)

        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("Folder:"))
        top.addWidget(self.folder_combo, 1)
        top.addWidget(browse_btn)
        top.addSpacing(20)
        top.addWidget(scan_group)
        top.addWidget(load_btn)
        top.addSpacing(20)
        top.addWidget(self.norm_chk)
        top.addWidget(QtWidgets.QLabel("Style:"))
        top.addWidget(self.style_combo)

        top_widget = QtWidgets.QWidget()
        top_widget.setLayout(top)

        self.tabs = QtWidgets.QTabWidget()
        self.eiger_tab = EigerTab()
        self.vortex_tab = VortexTab()
        self.tabs.addTab(self.eiger_tab, "Eiger")
        self.tabs.addTab(self.vortex_tab, "Vortex")
        self.eiger_tab.process_requested.connect(self._process_eiger)
        self.vortex_tab.process_requested.connect(self._process_vortex)
        self.eiger_tab.redraw_requested.connect(self._redraw_current_tab)
        self.vortex_tab.redraw_requested.connect(self._redraw_current_tab)
        self.eiger_tab.frame_requested.connect(self._on_eiger_frame_requested)
        self.vortex_tab.frame_requested.connect(self._on_vortex_frame_requested)
        self.tabs.currentChanged.connect(self._redraw_current_tab)

        central = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(central)
        v.addWidget(top_widget)
        v.addWidget(self.tabs, 1)
        self.setCentralWidget(central)

        self.status = self.statusBar()
        self.status.showMessage("Ready")

        # The settings controls go in the status bar, not the top row: the
        # window is plot-heavy, the top row is already full, and "where am I
        # saving, and am I saving" is the question a status bar exists to
        # answer. Permanent widgets sit to the right of showMessage(), so
        # they survive every status update.
        self.config_label = QtWidgets.QLabel()
        self.config_label.setStyleSheet("QLabel { color: #666; }")
        config_btn = QtWidgets.QPushButton("Settings file...")
        config_btn.setToolTip(
            "Choose where the settings are written.\n"
            "A destination, not a file to open: what is on screen is kept."
        )
        config_btn.clicked.connect(self._choose_config_path)
        self.autosave_btn = QtWidgets.QPushButton("Autosave")
        self.autosave_btn.setCheckable(True)
        self.autosave_btn.setToolTip(
            "Write the settings file on every change, starting now.\n"
            "Off at startup and not remembered between runs, so nothing is\n"
            "written to the working directory unless you ask for it."
        )
        self.autosave_btn.toggled.connect(self._on_autosave_toggled)
        self.status.addPermanentWidget(self.config_label)
        self.status.addPermanentWidget(config_btn)
        self.status.addPermanentWidget(self.autosave_btn)
        self._update_config_label()

        # Restore persisted state before wiring auto-save so re-applied values
        # don't immediately trigger a re-save with the same content.
        self._apply_loaded_config(
            _load_config(self._config_path),
            cli_folder=initial_folder,
            cli_cr=cli_cr,
            cli_cl=cli_cl,
        )

        for spin in (
            self.eiger_tab.x_cen,
            self.eiger_tab.y_cen,
            self.eiger_tab.width,
            self.eiger_tab.height,
            self.vortex_tab.e_start,
            self.vortex_tab.e_stop,
        ):
            spin.valueChanged.connect(self._save_state)
        for k in POL_KEYS:
            self.scan_input[k].editingFinished.connect(self._save_state)
        self.eiger_tab.scale_combo.currentIndexChanged.connect(self._save_state)
        self.vortex_tab.scale_combo.currentIndexChanged.connect(
            self._save_state
        )
        self.eiger_tab.diff_options.changed.connect(self._save_state)
        self.vortex_tab.diff_options.changed.connect(self._save_state)
        self.eiger_tab.pair_vlim.changed.connect(self._save_state)
        self.vortex_tab.pair_vlim.changed.connect(self._save_state)
        self.style_combo.currentIndexChanged.connect(self._save_state)

    # ------------------------------------------------------------ persistence
    def _apply_loaded_config(
        self, cfg, cli_folder=None, cli_cr=None, cli_cl=None
    ):
        self._loading_state = True
        try:
            history = cfg.get("folder_history") or []
            self._folder_history = [p for p in history if isinstance(p, str)]

            folder = cli_folder or cfg.get("folder") or os.getcwd()
            self._folder = folder
            if folder and folder not in self._folder_history:
                self._folder_history.insert(0, folder)
                del self._folder_history[FOLDER_HISTORY_MAX:]
            self._populate_folder_combo()

            eiger_roi = cfg.get("eiger_roi") or {}
            if not any(
                k in eiger_roi for k in ("x_cen", "y_cen", "width", "height")
            ) and all(
                k in eiger_roi for k in ("x_min", "x_max", "y_min", "y_max")
            ):
                try:
                    x0, x1 = int(eiger_roi["x_min"]), int(eiger_roi["x_max"])
                    y0, y1 = int(eiger_roi["y_min"]), int(eiger_roi["y_max"])
                    eiger_roi = {
                        "x_cen": (x0 + x1) // 2,
                        "y_cen": (y0 + y1) // 2,
                        "width": x1 - x0,
                        "height": y1 - y0,
                    }
                except (TypeError, ValueError):
                    eiger_roi = {}
            for name, spin in (
                ("x_cen", self.eiger_tab.x_cen),
                ("y_cen", self.eiger_tab.y_cen),
                ("width", self.eiger_tab.width),
                ("height", self.eiger_tab.height),
            ):
                if name in eiger_roi:
                    try:
                        spin.setValue(int(eiger_roi[name]))
                    except (TypeError, ValueError):
                        # Config file value is not an int; leave this ROI
                        # field at its default and take the rest.
                        pass

            vortex_roi = cfg.get("vortex_roi") or {}
            for name, spin in (
                ("e_start", self.vortex_tab.e_start),
                ("e_stop", self.vortex_tab.e_stop),
            ):
                if name in vortex_roi:
                    try:
                        spin.setValue(int(vortex_roi[name]))
                    except (TypeError, ValueError):
                        # As above: one bad ROI field does not sink the
                        # rest of the restored config.
                        pass

            cfg_scans = cfg.get("scans") or {}
            cli_scan = {"cr": cli_cr, "cl": cli_cl}
            for k in POL_KEYS:
                val = cli_scan.get(k) or cfg_scans.get(k)
                if val is not None:
                    self.scan_input[k].setText(str(val))

            scales = cfg.get("scales") or {}
            for name, combo in (
                ("eiger", self.eiger_tab.scale_combo),
                ("vortex", self.vortex_tab.scale_combo),
            ):
                val = scales.get(name)
                if val:
                    idx = combo.findText(str(val))
                    if idx >= 0:
                        combo.setCurrentIndex(idx)

            diff_opts = cfg.get("diff_options") or {}
            for name, tab in (
                ("eiger", self.eiger_tab),
                ("vortex", self.vortex_tab),
            ):
                state = diff_opts.get(name)
                if isinstance(state, dict):
                    tab.diff_options.apply_state(state)

            pair_vlim = cfg.get("pair_vlim") or {}
            for name, tab in (
                ("eiger", self.eiger_tab),
                ("vortex", self.vortex_tab),
            ):
                state = pair_vlim.get(name)
                if isinstance(state, dict):
                    tab.pair_vlim.apply_state(state)

            style = cfg.get("style")
            if style:
                idx = self.style_combo.findText(str(style))
                if idx >= 0:
                    self.style_combo.setCurrentIndex(idx)
        finally:
            self._loading_state = False

    def _populate_folder_combo(self):
        self.folder_combo.blockSignals(True)
        try:
            self.folder_combo.clear()
            for p in self._folder_history:
                self.folder_combo.addItem(p)
            idx = (
                self.folder_combo.findText(self._folder) if self._folder else -1
            )
            if idx < 0 and self._folder:
                self.folder_combo.insertItem(0, self._folder)
                idx = 0
            if idx >= 0:
                self.folder_combo.setCurrentIndex(idx)
        finally:
            self.folder_combo.blockSignals(False)

    def _add_folder_to_history(self, path):
        if path in self._folder_history:
            self._folder_history.remove(path)
        self._folder_history.insert(0, path)
        del self._folder_history[FOLDER_HISTORY_MAX:]
        self._populate_folder_combo()

    def _update_config_label(self):
        """Show the destination, elided from the left.

        From the left because these paths share their head and differ in
        their tail -- .../NdFeB vs .../test -- so the start is the part
        that can go.
        """
        self.config_label.setText(
            self.config_label.fontMetrics().elidedText(
                self._config_path, QtCore.Qt.TextElideMode.ElideLeft, 320
            )
        )
        self.config_label.setToolTip(self._config_path)

    def _choose_config_path(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Settings file",
            self._config_path,
            "JSON (*.json);;All files (*)",
        )
        if not path:
            return
        self._config_path = path
        self._update_config_label()
        # Only a destination is changed. Reading the file back would throw
        # away the settings on screen, which are the ones being worked with.
        if self.autosave_btn.isChecked():
            self._write_config()
        else:
            self.status.showMessage(f"Settings will be written to {path}")

    def _on_autosave_toggled(self, on):
        if not on:
            self.status.showMessage("Autosave off.")
            return
        # Write at once, so the file matches the screen rather than waiting
        # for whatever spinbox happens to be touched next.
        if self._write_config():
            self.status.showMessage(f"Autosaving to {self._config_path}")

    def _save_state(self):
        if self._loading_state or not self.autosave_btn.isChecked():
            return
        self._write_config()

    def _write_config(self):
        """Write the settings now. Returns whether it worked.

        A failed write switches autosave off instead of failing again on
        every later change: the path is wrong, only the user can fix it, so
        say so once and stop.
        """
        cfg = {
            "folder": self._folder,
            "folder_history": self._folder_history,
            "eiger_roi": {
                "x_cen": self.eiger_tab.x_cen.value(),
                "y_cen": self.eiger_tab.y_cen.value(),
                "width": self.eiger_tab.width.value(),
                "height": self.eiger_tab.height.value(),
            },
            "vortex_roi": {
                "e_start": self.vortex_tab.e_start.value(),
                "e_stop": self.vortex_tab.e_stop.value(),
            },
            "scans": {k: self.scan_input[k].text() for k in POL_KEYS},
            "scales": {
                "eiger": self.eiger_tab.scale_combo.currentText(),
                "vortex": self.vortex_tab.scale_combo.currentText(),
            },
            "diff_options": {
                "eiger": self.eiger_tab.diff_options.state(),
                "vortex": self.vortex_tab.diff_options.state(),
            },
            "pair_vlim": {
                "eiger": self.eiger_tab.pair_vlim.state(),
                "vortex": self.vortex_tab.pair_vlim.state(),
            },
            "style": self.style_combo.currentText(),
        }
        err = _save_config(cfg, self._config_path)
        if err is None:
            return True
        self.status.showMessage(f"Could not write {self._config_path}: {err}")
        self.autosave_btn.blockSignals(True)
        self.autosave_btn.setChecked(False)
        self.autosave_btn.blockSignals(False)
        return False

    def closeEvent(self, event):
        self._save_state()
        super().closeEvent(event)

    # ------------------------------------------------------------------ slots
    def _on_folder_selected(self, path):
        if not path or path == self._folder:
            return
        self._folder = path
        self._add_folder_to_history(path)
        self._save_state()

    def _browse_folder(self):
        start = self._folder or os.getcwd()
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Choose sample folder",
            start,
        )
        if not path:
            return
        self._folder = path
        self._add_folder_to_history(path)
        self._save_state()

    def _enabled_keys(self):
        return [k for k in POL_KEYS if self.enable_chk[k].isChecked()]

    def _reset_key(self, k):
        self._scans[k] = []
        self._xs[k] = self._ys[k] = self._i0s[k] = None
        self._eiger_z[k] = self._vortex_z[k] = None
        self._fmt_eiger[k] = []
        self._fmt_vortex[k] = []

    def _step_scans(self, key, delta):
        """Add ``delta`` to every scan number in one polarization's field."""
        line = self.scan_input[key]
        stepped = _offset_scan_spec(line.text(), delta)
        if stepped is None:
            self.status.showMessage(
                f"Cannot step {POL_LABEL[key]} scans: "
                f"'{line.text()}' is not a scan spec, or the step goes below 0."
            )
            return
        line.setText(stepped)
        if key == "cr":
            self._sync_cl_from_cr()

    def _sync_cl_from_cr(self, *_args):
        """Mirror CR into CL at +1, while both polarizations are enabled.

        Pairs are taken CR then CL on consecutive scans, so CL is almost always
        CR + 1; this keeps the two fields in step whether CR was typed or
        stepped with the arrows. Untick CL to break the link and set it by hand.
        """
        if self._syncing_scans:
            return
        if not (
            self.enable_chk["cr"].isChecked()
            and self.enable_chk["cl"].isChecked()
        ):
            return
        mirrored = _offset_scan_spec(self.scan_input["cr"].text(), 1)
        if mirrored is None:  # half-typed or unparseable -- leave CL alone
            return
        self._syncing_scans = True
        try:
            self.scan_input["cl"].setText(mirrored)
        finally:
            self._syncing_scans = False

    def _on_load(self):
        if not self._folder:
            self.status.showMessage("Pick a sample folder first.")
            return
        enabled = self._enabled_keys()
        if not enabled:
            self.status.showMessage(
                "Enable at least one polarization (CR or CL)."
            )
            return

        any_eiger = False
        any_vortex = False
        usable_i0 = False

        for k in POL_KEYS:
            if k not in enabled:
                self._reset_key(k)

        for k in enabled:
            try:
                requested = _parse_scan_list(self.scan_input[k].text())
            except ValueError as e:
                self.status.showMessage(f"{POL_LABEL[k]}: bad scan list: {e}")
                self._reset_key(k)
                continue
            if not requested:
                self.status.showMessage(f"{POL_LABEL[k]}: no scans specified")
                self._reset_key(k)
                continue

            # Load per-scan position streams; skip any scan whose
            # pos_stream is missing or fails to parse.
            loaded_scans = []
            i0s_list, xs_list, ys_list = [], [], []
            for s in requested:
                pp, pf = _resolve_h5(self._folder, "pos_stream", s)
                if pp is None:
                    self.status.showMessage(
                        f"{POL_LABEL[k]}: pos_stream missing for scan {s} (skipped)"
                    )
                    continue
                try:
                    i0s, xs, ys, _dx, _dy = process_position_stream(
                        s,
                        self._folder,
                        fname_format=pf,
                    )
                except Exception as e:
                    self.status.showMessage(
                        f"{POL_LABEL[k]}: pos_stream {s} failed: {e!r} (skipped)"
                    )
                    continue
                loaded_scans.append(s)
                i0s_list.append(i0s)
                xs_list.append(xs)
                ys_list.append(ys)

            if not loaded_scans:
                self._reset_key(k)
                continue

            self._scans[k] = loaded_scans
            self._i0s[k] = _sum_scan_arrays(i0s_list)
            self._xs[k] = _mean_scan_arrays(xs_list)
            self._ys[k] = _mean_scan_arrays(ys_list)
            self._eiger_z[k] = None
            self._vortex_z[k] = None
            usable_i0 = usable_i0 or _is_i0_usable(self._i0s[k])

            # Per-scan detector file formats. A polarization is considered
            # to have Eiger/Vortex data only if every scan does — otherwise
            # a summed z-array would silently mix subsets of scans.
            efmts, vfmts = [], []
            for s in loaded_scans:
                _, ef = _resolve_h5(self._folder, "eiger", s)
                _, vf = _resolve_h5(self._folder, "vortex", s)
                efmts.append(ef)
                vfmts.append(vf)
            self._fmt_eiger[k] = efmts
            self._fmt_vortex[k] = vfmts
            if all(f is not None for f in efmts):
                any_eiger = True
            if all(f is not None for f in vfmts):
                any_vortex = True

        # Preview comes from the first loaded scan of the first enabled key.
        eiger_preview = None
        eiger_preview_scan = None
        eiger_frame_count = 0
        vortex_preview = None
        vortex_preview_scan = None
        vortex_frame_count = 0
        for k in enabled:
            if not self._scans[k]:
                continue
            first_scan = self._scans[k][0]
            if eiger_preview is None:
                try:
                    f = _read_eiger_frame(self._folder, first_scan)
                    if f is not None:
                        eiger_preview = f
                        eiger_preview_scan = first_scan
                        eiger_frame_count = _get_frame_count(
                            self._folder,
                            "eiger",
                            first_scan,
                        )
                except Exception:
                    # The preview is a convenience for picking an ROI, and
                    # the scan may be mid-write, absent, or compressed with
                    # a codec this machine lacks. Any of those leaves
                    # eiger_preview None, which the panel draws as empty.
                    pass
            if vortex_preview is None:
                try:
                    s = _read_vortex_spectrum(self._folder, first_scan)
                    if s is not None:
                        vortex_preview = s
                        vortex_preview_scan = first_scan
                        vortex_frame_count = _get_frame_count(
                            self._folder,
                            "vortex",
                            first_scan,
                        )
                except Exception:
                    # As for the eiger preview above: no spectrum to show
                    # is a blank panel, not a failure to open the window.
                    pass

        self._eiger_preview_scan = eiger_preview_scan
        self._vortex_preview_scan = vortex_preview_scan

        # Keep both tabs enabled so loading data for only one detector
        # doesn't yank the user off the tab they were viewing.
        self.eiger_tab.set_preview(eiger_preview, frame_count=eiger_frame_count)
        self.eiger_tab.clear_maps()

        self.vortex_tab.set_preview(
            vortex_preview, frame_count=vortex_frame_count
        )
        self.vortex_tab.clear_maps()

        self.norm_chk.setEnabled(usable_i0)
        self.norm_chk.setChecked(usable_i0)

        loaded = [
            f"{POL_LABEL[k]}={len(self._scans[k])}"
            for k in enabled
            if self._scans[k]
        ]
        self.status.showMessage(
            f"Loaded (scans/pol): {', '.join(loaded) if loaded else 'none'} | "
            f"Eiger {'OK' if any_eiger else 'missing'} | "
            f"Vortex {'OK' if any_vortex else 'missing'} | "
            f"I0 {'usable' if usable_i0 else 'absent'}"
        )

        self._save_state()

    # ------------------------------------------------------------ processing
    def _start_worker(self, key, func, per_scan_kwargs, on_done):
        """Run ``func(scan, folder, **kwargs)`` for each scan of ``key`` off
        the GUI thread and sum the returned 1-D arrays element-wise.
        """
        scans = list(self._scans[key])
        folder = self._folder

        def summing():
            results = [
                func(s, folder, **kw) for s, kw in zip(scans, per_scan_kwargs)
            ]
            return _sum_scan_arrays(results)

        worker = ProcessWorker(summing, (), {})
        self._workers[key] = worker
        worker.finished_ok.connect(lambda res, k=key: on_done(k, res))
        worker.failed.connect(lambda msg, k=key: self._on_worker_failed(k, msg))
        worker.start()

    def _on_worker_failed(self, key, msg):
        self.status.showMessage(f"{POL_LABEL[key]}: processing failed: {msg}")
        self._pending.discard(key)
        self.eiger_tab.process_btn.setEnabled(True)
        self.vortex_tab.process_btn.setEnabled(True)

    def _process_eiger(self):
        keys = [k for k in self._enabled_keys() if self._scans[k]]
        if not keys:
            self.status.showMessage("Load at least one scan first.")
            return
        roi = self.eiger_tab.get_roi()
        self.eiger_tab.process_btn.setEnabled(False)
        skipped = [
            k for k in keys if any(f is None for f in self._fmt_eiger[k])
        ]
        runnable = [k for k in keys if k not in skipped]
        if skipped:
            self.status.showMessage(
                f"{', '.join(POL_LABEL[k] for k in skipped)}: Eiger missing "
                f"for some scans; skipping"
            )
        if not runnable:
            self.eiger_tab.process_btn.setEnabled(True)
            return
        self.status.showMessage(
            f"Processing Eiger ROI {roi} for {', '.join(POL_LABEL[k] for k in runnable)}..."
        )
        self._pending = set(runnable)
        for k in runnable:
            per_scan = [
                dict(roi=roi, batch=PROCESS_BATCH, fname_format=fmt)
                for fmt in self._fmt_eiger[k]
            ]
            self._start_worker(k, process_images, per_scan, self._on_eiger_done)

    def _on_eiger_done(self, key, sums):
        self._eiger_z[key] = sums
        self._pending.discard(key)
        if not self._pending:
            self.eiger_tab.process_btn.setEnabled(True)
            self._redraw_current_tab()
            self.status.showMessage("Eiger processed.")

    def _process_vortex(self):
        keys = [k for k in self._enabled_keys() if self._scans[k]]
        if not keys:
            self.status.showMessage("Load at least one scan first.")
            return
        roi = self.vortex_tab.get_roi()
        self.vortex_tab.process_btn.setEnabled(False)
        skipped = [
            k for k in keys if any(f is None for f in self._fmt_vortex[k])
        ]
        runnable = [k for k in keys if k not in skipped]
        if skipped:
            self.status.showMessage(
                f"{', '.join(POL_LABEL[k] for k in skipped)}: Vortex missing "
                f"for some scans; skipping"
            )
        if not runnable:
            self.vortex_tab.process_btn.setEnabled(True)
            return
        self.status.showMessage(
            f"Processing Vortex window {roi} for {', '.join(POL_LABEL[k] for k in runnable)}..."
        )
        self._pending = set(runnable)
        for k in runnable:
            per_scan = [
                dict(roi=roi, batch=PROCESS_BATCH, fname_format=fmt)
                for fmt in self._fmt_vortex[k]
            ]
            self._start_worker(
                k, process_vortex, per_scan, self._on_vortex_done
            )

    def _on_vortex_done(self, key, sums):
        self._vortex_z[key] = sums
        self._pending.discard(key)
        if not self._pending:
            self.vortex_tab.process_btn.setEnabled(True)
            self._redraw_current_tab()
            self.status.showMessage("Vortex processed.")

    def _collect_data(self, z_dict):
        """Return ``{key: (xs, ys, z, i0)}`` for scans whose z is ready and
        whose enable checkbox is on. Missing/unprocessed keys are excluded
        so the map panels draw empty for them.
        """
        out = {}
        for k in POL_KEYS:
            if not self.enable_chk[k].isChecked():
                continue
            if self._xs[k] is None or z_dict[k] is None:
                continue
            out[k] = (self._xs[k], self._ys[k], z_dict[k], self._i0s[k])
        return out

    def _redraw_current_tab(self):
        normalize = self.norm_chk.isChecked()
        scatter = self.style_combo.currentText() == "Scatter"
        idx = self.tabs.currentIndex()
        if idx == 0:
            data = self._collect_data(self._eiger_z)
            self.eiger_tab.draw_maps(data, normalize=normalize, scatter=scatter)
        elif idx == 1:
            data = self._collect_data(self._vortex_z)
            self.vortex_tab.draw_maps(
                data, normalize=normalize, scatter=scatter
            )

    def _on_eiger_frame_requested(self, index):
        if self._eiger_preview_scan is None:
            return
        try:
            frame = _read_eiger_frame(
                self._folder,
                self._eiger_preview_scan,
                index=index,
            )
        except Exception as e:
            self.status.showMessage(
                f"Failed to load Eiger frame {index}: {e!r}"
            )
            return
        if frame is not None:
            self.eiger_tab.update_preview_frame(frame, index)

    def _on_vortex_frame_requested(self, index):
        if self._vortex_preview_scan is None:
            return
        try:
            spectrum = _read_vortex_spectrum(
                self._folder,
                self._vortex_preview_scan,
                index=index,
            )
        except Exception as e:
            self.status.showMessage(
                f"Failed to load Vortex spectrum {index}: {e!r}"
            )
            return
        if spectrum is not None:
            self.vortex_tab.update_preview_frame(spectrum, index)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Flyscan dichroism viewer")
    parser.add_argument(
        "--folder",
        default=None,
        help="Sample folder (containing pos_stream/, eiger/, vortex/)",
    )
    parser.add_argument(
        "--cr",
        default=None,
        help="Initial CR scan(s): e.g. 490, '490,492,494', or '754-764:2'",
    )
    parser.add_argument(
        "--cl", default=None, help="Initial CL scan(s): same syntax as --cr"
    )
    parser.add_argument(
        "--config",
        default=None,
        help=(
            "Settings file to read at startup and autosave to "
            f"(default: ./{CONFIG_FILENAME})"
        ),
    )
    args = parser.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    win = MainWindow(
        initial_folder=args.folder,
        cli_cr=args.cr,
        cli_cl=args.cl,
        config_path=args.config,
    )
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
