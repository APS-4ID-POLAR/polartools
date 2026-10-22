"""A 2D image-stack viewer for the detector files that accompany a master file.

``MyNXWriter`` puts the scan's scalar channels in ``*_master.hdf`` but leaves the
area-detector frames out of it -- ``.../streams/primary/lambda250k`` carries only
``EPOCH``/``time`` and no ``value``.  The pixels are written separately, by one of
two writers, and the master records which one and where.

The pvapy streaming consumer writes a sibling file, named in
``entry/instrument/bluesky/metadata/pva_streaming``::

    file_name: pva_tise2_2_00006.h5
    file_path: ~/6idb-bits/sarah_26_3_2

An areaDetector plugin instead writes one file per detector into a subdirectory
of the master's own, listed in
``entry/instrument/bluesky/metadata/detectors_file_relative_path``::

    eiger: eiger/scan_000310.h5

-- or the literal ``{}`` when the scan wrote no frames at all. Both are resolved
against the master's own directory: ``pva_streaming``'s ``file_path`` and the
companion ``detectors_file_full_path`` name the filesystem the scan was taken
on, which a copy of the data is not. :func:`image_path_for_master` falls back to
deriving ``pva_<stem without "_master">.h5`` when neither is recorded, and opens
the first detector where a master names several.

The image file's layout. The pvapy writer produces all of this::

    entry/data/data                (N, nx, ny)  the stack; ``pixel_order='xy'``,
                                   i.e. ``frame[x, y]`` -- the same convention as
                                   pyqtgraph's default col-major ImageItem
    entry/data/frame_timestamp     (N,)         epoch seconds
    entry/data/rois/<name>         (N, w, h)    pre-baked ROIs; the rectangle is in
                                   the ``x``/``y``/``w``/``h`` attributes
    entry/data/metadata/ca/<pv>    (N,)         one value per frame (eta, delta, cap, ...)

An areaDetector file carries the stack and nothing else from that list -- no
stored ROIs, no frame_timestamp -- but it does bring its own per-frame
attributes, which are read in the same way as the CA metadata above::

    entry/instrument/NDAttributes/NDArrayEpicsTSSec   (N,)  the frame clock,
    entry/instrument/NDAttributes/NDArrayEpicsTSnSec  (N,)  counted from 1990
    entry/instrument/NDAttributes/<name>              (N,)  BeamCenterX, ...

Its frames are stored ``[frame, x, y]`` too, despite declaring no
``pixel_order``: on an EIGER2 X 1M the 38-pixel module gap falls across the
*second* frame axis, splitting it 512/38/512, so that axis is the detector's
height and the first is its width. Nothing needs transposing to display.

The frame count need not match the master scan's point count (the detector
usually starts a little early), so this window is deliberately self-contained:
every per-frame quantity it plots comes from the image file itself.

``entry/data/data`` is compressed with a filter ``hdf5plugin`` registers (blosc
from the pvapy writer, bitshuffle from areaDetector), so that package must be
importable. It is imported here if present; if it is not and the pixel read
fails, :class:`ImageData` raises a :class:`RuntimeError` naming it.

Where the matching ``*_master.hdf`` can be found -- beside the image file, or one
directory up when the frames sit in a per-detector subdirectory -- its primary
stream is offered too, under a ``scan/`` prefix, each frame taking the scan point
nearest it in time, since the detector free-runs and the two files rarely have
the same number of points.

Both ROI plots can **normalize** by any of the same quantities their axes offer
(a CA channel such as a monitor, or another ROI metric): the value axes are
divided by it point by point, after any background subtraction.

Any ROI can be marked as the **background**: its mean intensity per pixel is then
subtracted from every other ROI's sum (as ``background x pixel count``) and from
their mean/min/max, per frame and per file.

An **HKL 3D** button maps the frames into reciprocal space, taking the geometry
from the matching master file (see ``hkl_viewer.py``); it needs ``pyvista`` and
``xrayutilities``, and says so plainly when they are absent.

Several files can be opened at once: they are shown side by side in one window,
sharing a single set of ROIs, one colour scale and one frame control, so scans can
be compared directly. The frame control spans the longest file; a file that has run
out of frames keeps showing its last one.

Run standalone on one or more image files with::

    python image_viewer.py scan_000310.h5 [scan_000311.h5 ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

try:  # Registers the blosc/lz4/bitshuffle filters that entry/data/data needs.
    import hdf5plugin  # noqa: F401
    HAVE_HDF5PLUGIN = True
except ImportError:  # pragma: no cover - depends on the environment
    HAVE_HDF5PLUGIN = False

import h5py

# PyQt6 ahead of pyqtgraph deliberately, and not in alphabetical order.
# pyqtgraph binds to whichever Qt wrapper is already in sys.modules and
# otherwise tries PyQt5 first, so in an environment holding both bindings an
# alphabetised import leaves pyqtgraph on PyQt5 while this module is on PyQt6
# -- two QObject hierarchies in one process, surfacing as opaque TypeErrors
# from pg.RectROI rather than as an import error.
from PyQt6 import QtCore
from PyQt6 import QtGui
from PyQt6 import QtWidgets

import pyqtgraph as pg

IMAGE_DATA_PATH = "entry/data/data"
FRAME_TIMESTAMP_PATH = "entry/data/frame_timestamp"
ROI_GROUP_PATH = "entry/data/rois"
# areaDetector's own per-frame attributes are read just like the CA metadata:
# same shape, same one-value-per-frame meaning, and nothing in them collides.
CA_GROUP_PATHS = (
    "entry/data/metadata/ca",
    "entry/data/metadata/motor_positions",
    "entry/instrument/NDAttributes",
)
EPICS_TS_SEC_PATH = "entry/instrument/NDAttributes/NDArrayEpicsTSSec"
EPICS_TS_NSEC_PATH = "entry/instrument/NDAttributes/NDArrayEpicsTSnSec"
# EPICS time counts from 1990-01-01, UNIX time from 1970-01-01.
EPICS_EPOCH_OFFSET = 631152000.0
PVA_STREAMING_PATH = "entry/instrument/bluesky/metadata/pva_streaming"
DETECTORS_FILE_PATH = "entry/instrument/bluesky/metadata/detectors_file_relative_path"
# Both writers record their file in a "name: value" block. The pvapy consumer
# labels its line "file_name"; an areaDetector writer labels each line with the
# detector's own name, so there the first line simply wins.
#
# Only the relative form of the areaDetector path is read, never the
# detectors_file_full_path beside it. An absolute path does not merely fail
# after the data is copied elsewhere -- while the original mount is still
# visible it succeeds, quietly opening the original file instead of the copy
# that was clicked. A relative path can reach anywhere the full one can,
# including out of the tree via "..", so it gives up nothing.
RECORDED_IMAGE_SOURCES = (
    (PVA_STREAMING_PATH, "file_name"),
    (DETECTORS_FILE_PATH, None),
)
MASTER_PRIMARY_PATH = "entry/instrument/bluesky/streams/primary"
# Master-file channels are offered alongside the image file's own, under this
# prefix so their provenance is plain and they cannot collide.
SCAN_PREFIX = "scan/"

# Refuse rather than thrash: the stack is read into memory in one go (a 441 MB
# 101x1028x1062 eiger stack decodes in well under a second).
MAX_STACK_BYTES = 8 * 1024 ** 3
# Several files are held at once when they are compared side by side.
MAX_TOTAL_BYTES = 8 * 1024 ** 3

HDF5PLUGIN_HINT = (
    "The image stack is compressed and could not be decoded.\n\n"
    "This needs the 'hdf5plugin' package importable in the Python running this "
    "viewer:\n\n    pip install hdf5plugin\n\n"
    "(or conda install -c conda-forge hdf5plugin)"
)


# -- Locating the image file ------------------------------------------------


def _decode(value) -> str:
    """Return *value* -- an h5py scalar that may be bytes or a 0-d array -- as text."""
    value = np.asarray(value).reshape(-1)
    if value.size == 0:
        return ""
    item = value[0]
    if isinstance(item, bytes):
        return item.decode(errors="replace")
    return str(item)


def _block_value(text: str, key: str | None) -> str | None:
    """The value *key* carries in a ``name: value`` block, or the first line's.

    Hand-parsed, because neither viewer depends on PyYAML for the handful of
    metadata blocks it reads. An empty value counts as no value: a scan that
    wrote no frames records the literal ``{}``, which has no colon at all, but
    a truncated ``eiger:`` would otherwise yield an empty relative path, and
    ``master.parent / ""`` is the master's directory rather than a file.
    """
    for line in text.splitlines():
        name, sep, value = line.partition(":")
        if not sep:
            continue
        if key is not None and name.strip() != key:
            continue
        found = value.strip().strip("'\"")
        if found:
            return found
    return None


def _recorded_image_paths(master: Path) -> list[str]:
    """The image-file paths *master* records, in the order worth trying them.

    One open covers both writers' conventions. The browser re-resolves every
    row of the listing whenever the selection changes, and these files sit on
    a beamline mount, so a second open per master is worth avoiding.
    """
    found: list[str] = []
    try:
        with h5py.File(master, "r") as f:
            for path, key in RECORDED_IMAGE_SOURCES:
                item = f.get(path)
                if item is None:
                    continue
                value = _block_value(_decode(item[()]), key)
                if value:
                    found.append(value)
    except OSError:
        return []
    return found


def recorded_detector_file(master: Path, detector: str | None = None) -> str | None:
    """The file *master* records for *detector*, exactly as written.

    The one line of the ``detectors_file_relative_path`` block that *detector*
    labels -- ``"vortex"`` picks ``vortex: vortex/scan_000055.h5`` -- or the
    first line's value when *detector* is ``None``.

    Public, and by name rather than by position, for ``spectrum_viewer``: that
    module wants one detector's line out of a master that may list several,
    where :func:`_recorded_image_paths` wants every source in the order worth
    opening and takes the first line when a block is unlabelled. Reaching
    across for :func:`_block_value` and :func:`_decode` instead would make
    their underscores a lie, and hand-parsing the block a second time in
    another file would duplicate exactly the ``{}`` and truncated-``eiger:``
    reasoning those two were written down once to settle.

    What comes back is a *claim*, not a path that exists: see
    :func:`spectrum_path_for_master` for how far it can be trusted.
    """
    try:
        with h5py.File(master, "r") as f:
            item = f.get(DETECTORS_FILE_PATH)
            if item is None:
                return None
            return _block_value(_decode(item[()]), detector)
    except OSError:
        return None


def image_path_for_master(master: Path) -> Path | None:
    """The image file belonging to *master*, or ``None`` if none was found.

    Prefers what the writer recorded -- a pvapy sibling, or an areaDetector
    file in a per-detector subdirectory -- and falls back to the ``<master stem
    minus "_master">`` -> ``pva_<stem>.h5`` convention. Everything is resolved
    against the master's own directory, so that a copy of the data opens its
    own frames rather than the originals it was copied from. Where a master
    names several detectors, the first is the one opened.
    """
    master = Path(master)
    candidates: list[Path] = []

    for recorded in _recorded_image_paths(master):
        recorded = Path(recorded)
        # Not just tidiness: pathlib's "/" throws the left side away entirely
        # when the right is absolute, so an absolute value here would silently
        # escape the master's directory and defeat the whole policy above.
        if recorded.is_absolute():
            recorded = Path(recorded.name)
        candidates.append(master.parent / recorded)

    stem = master.stem
    if stem.endswith("_master"):
        stem = stem[: -len("_master")]
    for prefix in ("pva_", "pva"):
        for suffix in (".h5", ".hdf5", ".hdf"):
            candidates.append(master.parent / f"{prefix}{stem}{suffix}")

    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def master_path_for_image(image: Path) -> Path | None:
    """The master file this image file belongs to, beside it or one level up.

    An areaDetector writer puts the frames in a per-detector subdirectory of
    the master's own, so ``<dir>/eiger/scan_000310.h5`` belongs to
    ``<dir>/scan_000310_master.hdf``. The image's own directory is searched
    first, so a pvapy sibling still pairs with the master beside it.

    Looking a level up can in principle pair a stray file with a master that
    merely shares its stem; the time-overlap check in
    :func:`read_scan_channels` is what rejects those.
    """
    image = Path(image)
    stem = image.stem
    for prefix in ("pva_", "pva"):
        if stem.startswith(prefix):
            stem = stem[len(prefix):]
            break
    for directory in (image.parent, image.parent.parent):
        for suffix in (".hdf", ".hdf5", ".h5"):
            candidate = directory / f"{stem}_master{suffix}"
            if candidate.is_file():
                return candidate
    return None


def _scan_time_base(group) -> np.ndarray | None:
    """A primary field's EPOCH that is a real, advancing clock.

    Not every field has one: an ``EpicsSignal`` that was never re-read keeps the
    EPICS zero (631152000.0, i.e. 1990) for every point, which places nothing in
    time. One good clock is used for the whole stream, so those fields are
    carried along by it rather than dropped.
    """
    for name in group.keys():
        item = group[name]
        if not isinstance(item, h5py.Group) or "EPOCH" not in item:
            continue
        epoch = np.asarray(item["EPOCH"][()], dtype=float).ravel()
        if (
            epoch.size > 1
            and epoch[0] > 1e9                      # later than 2001, so a real time
            and np.all(np.diff(epoch) >= 0)
            and epoch[-1] > epoch[0]
        ):
            return epoch
    return None


def nearest_index(frame_epoch: np.ndarray, scan_epoch: np.ndarray) -> np.ndarray:
    """For each frame, the scan point closest to it in time.

    Shared with ``hkl_viewer``, which matches the master's diffractometer angles
    onto the frames the same way.
    """
    right = np.clip(np.searchsorted(scan_epoch, frame_epoch), 1, scan_epoch.size - 1)
    left = right - 1
    take_left = frame_epoch - scan_epoch[left] <= scan_epoch[right] - frame_epoch
    return np.where(take_left, left, right)


def _evenly_spaced_index(n_frames: int, n_scan_points: int) -> np.ndarray:
    """Map each frame onto a scan point by position, not by time.

    Used only when the image file has no per-frame timestamp to match against
    (see ``read_scan_channels``), so frame 0 takes scan point 0, the last frame
    takes the last scan point, and the rest are spread evenly in between.
    """
    if n_frames <= 1:
        return np.zeros(max(n_frames, 0), dtype=int)
    positions = np.arange(n_frames) * (n_scan_points - 1) / (n_frames - 1)
    return np.round(positions).astype(int)


def _scan_point_count(group) -> int | None:
    """How many points the master's primary stream has.

    Used in place of ``_scan_time_base`` when there is no frame timestamp to
    read a count off of a real clock: the most common length among the
    stream's own value arrays, since a few fields (e.g. read-once presets) can
    be shorter or longer than the scan itself.
    """
    counts: dict[int, int] = {}
    for name in group.keys():
        item = group[name]
        if not isinstance(item, h5py.Group) or "value" not in item:
            continue
        try:
            size = int(np.asarray(item["value"]).size)
        except (TypeError, ValueError):
            continue
        if size > 1:
            counts[size] = counts.get(size, 0) + 1
    if not counts:
        return None
    return max(counts.items(), key=lambda kv: kv[1])[0]


def read_scan_channels(
    image: Path, frame_epoch, n_frames: int | None = None
) -> dict[str, np.ndarray]:
    """The master file's primary fields, resampled onto the image's frames.

    The detector free-runs, so the two files rarely have the same number of
    points -- 103 frames against 101 scan points is typical, the camera starting
    a moment early. Each frame therefore takes the value of the scan point
    nearest it in time; frames outside the scan hold the first or last point.

    When the image file carries no per-frame timestamp at all, there is no
    clock to match against, so frames are instead spread evenly across the
    master's scan points by position (see ``_evenly_spaced_index``) -- a
    rougher alignment, but the best available without a time base.
    """
    have_epoch = frame_epoch is not None and len(frame_epoch) > 0
    if not have_epoch and not n_frames:
        return {}
    master = master_path_for_image(image)
    if master is None:
        return {}
    try:
        with h5py.File(master, "r") as f:
            group = f.get(MASTER_PRIMARY_PATH)
            if group is None:
                return {}

            if have_epoch:
                frame_epoch = np.asarray(frame_epoch, dtype=float)
                scan_epoch = _scan_time_base(group)
                if scan_epoch is None:
                    return {}
                index = nearest_index(frame_epoch, scan_epoch)

                # Refuse to pair files that do not describe the same stretch
                # of time -- a mismatched master would otherwise be resampled
                # into plausible-looking nonsense.
                spacing = float(np.median(np.diff(scan_epoch))) or 1.0
                if float(np.median(np.abs(frame_epoch - scan_epoch[index]))) > 5 * spacing:
                    return {}
                n_scan_points = scan_epoch.size
            else:
                n_scan_points = _scan_point_count(group)
                if not n_scan_points:
                    return {}
                index = _evenly_spaced_index(n_frames, n_scan_points)

            channels: dict[str, np.ndarray] = {}
            for name in group.keys():
                item = group[name]
                if not isinstance(item, h5py.Group) or "value" not in item:
                    continue
                try:
                    values = np.asarray(item["value"][()], dtype=float).ravel()
                except (TypeError, ValueError):
                    continue
                if values.size == n_scan_points:
                    channels[SCAN_PREFIX + name] = values[index]
            return channels
    except OSError:
        return {}


def is_image_file(path: Path) -> bool:
    """Whether *path* is itself an image file (i.e. it holds an image stack)."""
    try:
        with h5py.File(path, "r") as f:
            return isinstance(f.get(IMAGE_DATA_PATH), h5py.Dataset)
    except OSError:
        return False


# Resolution opens the file to read its metadata, and the browser re-resolves
# the whole selection on every selection change, so remember the answer per
# (path, mtime, size) -- cheap locally, but worth it over a beamline mount.
_RESOLVE_CACHE: dict[tuple, "Path | None"] = {}


def resolve_image_file(path: Path) -> Path | None:
    """The image file to open for *path*, which may be a master or an image file.

    Selecting an image file directly in the file list is the obvious way to ask
    for it, so that opens itself rather than hunting for a partner.
    """
    path = Path(path)
    try:
        stat = path.stat()
        key = (str(path), stat.st_mtime, stat.st_size)
    except OSError:
        key = None
    if key is not None and key in _RESOLVE_CACHE:
        return _RESOLVE_CACHE[key]

    resolved = image_path_for_master(path)
    if resolved is None and is_image_file(path):
        resolved = path
    if key is not None:
        _RESOLVE_CACHE[key] = resolved
    return resolved


# -- Reading the image file -------------------------------------------------


class ImageData:
    """The whole contents of one detector image file, held in memory.

    ``stack`` is ``(n_frames, nx, ny)`` float32 -- indexed ``[frame, x, y]``,
    the pvapy writer's declared ``pixel_order='xy'`` convention and, as the
    module docstring works out from where the module gap falls, areaDetector's
    layout here as well. That is also what pyqtgraph's default (col-major)
    ``ImageItem`` expects, so frames need no transposing to display the right
    way round.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        with h5py.File(self.path, "r") as f:
            dset = f.get(IMAGE_DATA_PATH)
            if dset is None or not isinstance(dset, h5py.Dataset):
                raise RuntimeError(
                    f"{self.path.name} has no image stack at '{IMAGE_DATA_PATH}'."
                )
            if dset.ndim not in (2, 3):
                raise RuntimeError(
                    f"'{IMAGE_DATA_PATH}' is {dset.ndim}-dimensional; "
                    "only 2D images and 3D stacks are supported."
                )
            if dset.nbytes > MAX_STACK_BYTES:
                raise RuntimeError(
                    f"{self.path.name} holds {dset.nbytes / 1024 ** 3:.1f} GB of image "
                    f"data, more than this viewer loads at once "
                    f"({MAX_STACK_BYTES / 1024 ** 3:.0f} GB)."
                )
            try:
                stack = np.asarray(dset[()], dtype=np.float32)
            except OSError as exc:
                if not HAVE_HDF5PLUGIN:
                    raise RuntimeError(HDF5PLUGIN_HINT) from exc
                raise RuntimeError(f"Could not read the image stack: {exc}") from exc

            if stack.ndim == 2:
                stack = stack[np.newaxis, ...]
            self.stack = stack
            n_frames = stack.shape[0]

            self.epoch = self._read_epoch(f, n_frames)
            self.timestamps = (
                None if self.epoch is None else self.epoch - self.epoch[0]
            )
            self.channels = self._read_channels(f, n_frames)
            self.stored_rois = self._read_stored_rois(f)

        # Offered next to the image file's own channels, on the same frames.
        self.channels.update(read_scan_channels(self.path, self.epoch, n_frames))

    @property
    def n_frames(self) -> int:
        return int(self.stack.shape[0])

    @property
    def frame_shape(self) -> tuple[int, int]:
        return int(self.stack.shape[1]), int(self.stack.shape[2])

    @staticmethod
    def _read_epoch(f: h5py.File, n_frames: int) -> np.ndarray | None:
        """Absolute frame times, needed to line the frames up with the scan.

        The pvapy consumer writes them outright. An areaDetector file does not,
        but it stamps every frame with the EPICS time its driver saw, split
        across two attributes, so that is read and shifted onto the UNIX epoch
        instead.
        """
        item = f.get(FRAME_TIMESTAMP_PATH)
        if item is not None:
            stamps = np.asarray(item[()], dtype=float).ravel()
            if stamps.size == n_frames and stamps.size:
                return stamps

        sec = f.get(EPICS_TS_SEC_PATH)
        nsec = f.get(EPICS_TS_NSEC_PATH)
        if sec is None or nsec is None:
            return None
        # Both are uint32 in the file; the seconds alone overflow a float32 and
        # the nanoseconds have to be scaled, so widen before doing any of it.
        seconds = np.asarray(sec[()], dtype=float).ravel()
        nanoseconds = np.asarray(nsec[()], dtype=float).ravel()
        if not n_frames or seconds.size != n_frames or nanoseconds.size != n_frames:
            return None
        epoch = seconds + nanoseconds * 1e-9 + EPICS_EPOCH_OFFSET
        # An attribute the driver never filled in leaves zeroes behind, which
        # place every frame in 1990 -- not a clock at all, so say there is none
        # and let the caller fall back to spreading the frames by position.
        if epoch[0] < 1e9:
            return None
        return epoch

    @staticmethod
    def _read_channels(f: h5py.File, n_frames: int) -> dict[str, np.ndarray]:
        """Per-frame scalar channels (the CA metadata), keyed by name."""
        channels: dict[str, np.ndarray] = {}
        for group_path in CA_GROUP_PATHS:
            group = f.get(group_path)
            if not isinstance(group, h5py.Group):
                continue
            for name in group.keys():
                if name in channels:
                    continue
                item = group.get(name)
                if not isinstance(item, h5py.Dataset):
                    continue
                try:
                    values = np.asarray(item[()], dtype=float).ravel()
                except (OSError, TypeError, ValueError):
                    continue
                if values.size == n_frames:
                    channels[name] = values
        return channels

    @staticmethod
    def _read_stored_rois(f: h5py.File) -> list[dict]:
        """The rectangles saved in ``entry/data/rois`` (x/y/w/h attributes)."""
        group = f.get(ROI_GROUP_PATH)
        if not isinstance(group, h5py.Group):
            return []
        rois = []
        for name in sorted(group.keys()):
            item = group.get(name)
            if not isinstance(item, h5py.Dataset):
                continue  # e.g. the 'info' subgroup
            attrs = item.attrs
            if not all(key in attrs for key in ("x", "y", "w", "h")):
                continue
            rois.append(
                {
                    "name": str(name),
                    "x": float(attrs["x"]),
                    "y": float(attrs["y"]),
                    "w": max(1.0, float(attrs["w"])),
                    "h": max(1.0, float(attrs["h"])),
                    "angle": float(attrs.get("angle", 0.0)),
                }
            )
        return rois


# -- Metrics ----------------------------------------------------------------

METRIC_KEYS = ["frame", "time", "sum", "mean", "min", "max", "comx", "comy"]
SINGLE_FRAME_KEYS = ["proj_x", "proj_y"]
NO_NORMALIZATION = "(none)"

DISPLAY_NAMES = {
    "frame": "Frame",
    "time": "Time (s)",
    "sum": "Sum",
    "mean": "Mean",
    "min": "Min",
    "max": "Max",
    "comx": "CoM X",
    "comy": "CoM Y",
    "proj_x": "Proj X (sum over Y)",
    "proj_y": "Proj Y (sum over X)",
}

AXIS_LABELS = {
    "frame": "Frame index",
    "time": "Time since first frame [s]",
    "sum": "ROI sum",
    "mean": "ROI mean",
    "min": "ROI min",
    "max": "ROI max",
    "comx": "ROI CoM X [px from ROI corner]",
    "comy": "ROI CoM Y [px from ROI corner]",
    "proj_x": "Sum over Y",
    "proj_y": "Sum over X",
    "pixel": "Pixel index within ROI",
}


def add_metric_item(combo: QtWidgets.QComboBox, key: str):
    combo.addItem(DISPLAY_NAMES.get(key, key), key)


def combo_key(combo: QtWidgets.QComboBox, default: str = "frame") -> str:
    data = combo.currentData()
    return str(data) if data is not None else (combo.currentText() or default)


def set_combo_key(combo: QtWidgets.QComboBox, key: str):
    index = combo.findData(key)
    if index >= 0:
        combo.setCurrentIndex(index)


def axis_label(key: str) -> str:
    return AXIS_LABELS.get(key, key)


# -- ROI geometry and per-frame reduction -----------------------------------


def roi_is_rotated(roi) -> bool:
    try:
        return abs(float(roi.angle())) % 360.0 > 1e-6
    except Exception:
        return False


def roi_bbox(roi, shape: tuple[int, int]) -> tuple[int, int, int, int]:
    """The ROI's axis-aligned pixel bounds ``(x0, x1, y0, y1)``, clipped to *shape*."""
    pos, size = roi.pos(), roi.size()
    nx, ny = shape
    x0 = int(np.floor(pos.x()))
    y0 = int(np.floor(pos.y()))
    x1 = int(np.ceil(pos.x() + size.x()))
    y1 = int(np.ceil(pos.y() + size.y()))
    x0 = max(0, min(x0, nx))
    x1 = max(x0, min(x1, nx))
    y0 = max(0, min(y0, ny))
    y1 = max(y0, min(y1, ny))
    return x0, x1, y0, y1


def roi_sub_frame(frame: np.ndarray, roi, image_item=None) -> np.ndarray | None:
    """The pixels of one frame inside *roi*.

    An unrotated ROI is a plain slice, which is both faster and exactly what
    :func:`roi_series` reduces over -- going through ``getArrayRegion`` here
    instead would round fractional edges differently and make the stats line
    disagree with the plotted curve. Only a rotated ROI needs the resampling.
    """
    sub = None
    if roi_is_rotated(roi) and image_item is not None:
        try:
            sub = roi.getArrayRegion(frame, image_item)
            if sub is not None and getattr(sub, "ndim", 0) > 2:
                sub = np.squeeze(sub)
        except Exception:
            sub = None
    if sub is None or getattr(sub, "size", 0) == 0:
        x0, x1, y0, y1 = roi_bbox(roi, (frame.shape[0], frame.shape[1]))
        sub = frame[x0:x1, y0:y1] if x1 > x0 and y1 > y0 else None
    if sub is None or sub.size == 0:
        return None
    return np.asarray(sub, dtype=np.float32)


def roi_pixel_count(roi, shape: tuple[int, int], image_item=None, frame=None) -> int:
    """How many pixels *roi* covers -- the N in a ``background x N`` correction."""
    if roi_is_rotated(roi) and frame is not None:
        sub = roi_sub_frame(frame, roi, image_item)
        return int(sub.size) if sub is not None else 0
    x0, x1, y0, y1 = roi_bbox(roi, shape)
    return max(0, (x1 - x0) * (y1 - y0))


def _empty_series(n_frames: int) -> dict[str, np.ndarray]:
    zeros = np.zeros(n_frames, dtype=float)
    return {key: zeros.copy() for key in ("sum", "mean", "min", "max", "comx", "comy")}


def _reduce_sub(sub: np.ndarray) -> tuple[float, float, float, float, float, float]:
    """``(sum, mean, min, max, comx, comy)`` of one ROI sub-image indexed ``[x, y]``."""
    total = float(np.sum(sub))
    weight = total if total != 0.0 else 1.0
    comx = float((np.arange(sub.shape[0]) @ sub.sum(axis=1)) / weight)
    comy = float((sub.sum(axis=0) @ np.arange(sub.shape[1])) / weight)
    return total, float(np.mean(sub)), float(np.min(sub)), float(np.max(sub)), comx, comy


def roi_series(stack: np.ndarray, roi, image_item=None) -> dict[str, np.ndarray]:
    """Per-frame ROI metrics over the whole stack.

    An unrotated ROI is a plain slice of the stack, so the reductions run over
    all frames at once; a rotated one falls back to ``getArrayRegion`` per frame.
    """
    n_frames = int(stack.shape[0])

    if roi_is_rotated(roi):
        series = _empty_series(n_frames)
        for i in range(n_frames):
            sub = roi_sub_frame(stack[i], roi, image_item)
            if sub is None:
                continue
            values = _reduce_sub(sub)
            for key, value in zip(("sum", "mean", "min", "max", "comx", "comy"), values):
                series[key][i] = value
        return series

    x0, x1, y0, y1 = roi_bbox(roi, (stack.shape[1], stack.shape[2]))
    if x1 <= x0 or y1 <= y0:
        return _empty_series(n_frames)

    sub = stack[:, x0:x1, y0:y1]
    sums = sub.sum(axis=(1, 2), dtype=np.float64)
    weights = np.where(sums != 0.0, sums, 1.0)
    return {
        "sum": sums,
        "mean": sub.mean(axis=(1, 2), dtype=np.float64),
        "min": sub.min(axis=(1, 2)).astype(float),
        "max": sub.max(axis=(1, 2)).astype(float),
        "comx": (sub.sum(axis=2, dtype=np.float64) @ np.arange(x1 - x0)) / weights,
        "comy": (sub.sum(axis=1, dtype=np.float64) @ np.arange(y1 - y0)) / weights,
    }


# -- ROIs, shared across every open image -----------------------------------

ROI_COLORS = [
    (255, 0, 0),
    (0, 255, 0),
    (80, 140, 255),
    (255, 255, 0),
    (255, 0, 255),
    (0, 255, 255),
]

# Curve colours for the per-file traces in a shared ROI plot.
FILE_COLORS = [
    (255, 220, 60),
    (90, 200, 255),
    (255, 120, 120),
    (140, 230, 140),
    (220, 150, 255),
    (255, 175, 80),
]


class ContextRectROI(pg.RectROI):
    """One panel's drawn rectangle for a :class:`SharedRoi`.

    Locking freezes position and size: ``translatable`` stops the body being
    dragged, and :meth:`checkPointMove` -- the hook pyqtgraph's handles consult
    before every resize or rotate -- refuses the move.
    """

    def __init__(self, shared, panel, pos, size, pen):
        super().__init__(pos, size, pen=pen)
        self.shared = shared
        self.panel = panel
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.LeftButton | QtCore.Qt.MouseButton.RightButton)
        self.addRotateHandle([1, 0], [0.5, 0.5])

    def checkPointMove(self, handle, pos, modifiers):
        return not self.shared.locked

    def mouseClickEvent(self, ev):
        if ev.button() == QtCore.Qt.MouseButton.RightButton:
            self.shared.viewer.show_roi_menu(self.shared, QtGui.QCursor.pos())
            ev.accept()
            return
        super().mouseClickEvent(ev)


class SharedRoi(QtCore.QObject):
    """One logical ROI, drawn once per open image and kept in step.

    Every panel shows the same detector, so a rectangle means the same pixels in
    all of them. Dragging any one copy mirrors the geometry onto the others; the
    ROI plots listen to this object rather than to the individual rectangles.
    """

    geometryChanged = QtCore.pyqtSignal(object)
    geometryFinished = QtCore.pyqtSignal(object)

    def __init__(self, viewer, name: str, color, pos, size, angle: float = 0.0):
        super().__init__(viewer)
        self.viewer = viewer
        self.name = name
        self.color = color
        self.locked = False
        self.visible = True
        self.items: list[ContextRectROI] = []
        self._pos = (float(pos[0]), float(pos[1]))
        self._size = (float(size[0]), float(size[1]))
        self._angle = float(angle)
        self._mirroring = False

    # -- geometry, read the way roi_bbox()/roi_series() expect ------------

    def pos(self):
        return pg.Point(*self._pos)

    def size(self):
        return pg.Point(*self._size)

    def angle(self) -> float:
        return self._angle

    # -- items -------------------------------------------------------------

    def add_item(self, panel) -> ContextRectROI:
        item = ContextRectROI(
            self, panel, self._pos, self._size, pg.mkPen(self.color, width=2)
        )
        item.setAngle(self._angle, finish=False)
        item.sigRegionChanged.connect(self._on_item_changed)
        item.sigRegionChangeFinished.connect(self._on_item_finished)
        panel.image_view.addItem(item)
        self.items.append(item)
        self._apply_state(item)
        return item

    def item_for(self, panel) -> ContextRectROI | None:
        for item in self.items:
            if item.panel is panel:
                return item
        return None

    def remove_items(self):
        for item in self.items:
            item.panel.image_view.removeItem(item)
        self.items.clear()

    # -- keeping the copies in step ---------------------------------------

    def _write(self, item):
        item.blockSignals(True)
        try:
            item.setPos(pg.Point(*self._pos), finish=False)
            item.setSize(pg.Point(*self._size), finish=False)
            item.setAngle(self._angle, finish=False)
        finally:
            item.blockSignals(False)

    def _on_item_changed(self, item):
        if self._mirroring:
            return
        self._mirroring = True
        try:
            self._pos = (float(item.pos().x()), float(item.pos().y()))
            self._size = (float(item.size().x()), float(item.size().y()))
            self._angle = float(item.angle())
            for other in self.items:
                if other is not item:
                    self._write(other)
        finally:
            self._mirroring = False
        self.geometryChanged.emit(self)

    def _on_item_finished(self, _item):
        if self._mirroring:
            return
        self.geometryFinished.emit(self)

    def set_geometry(self, x: float, y: float, w: float, h: float, angle: float):
        self._pos = (float(x), float(y))
        self._size = (float(w), float(h))
        self._angle = float(angle)
        self._mirroring = True
        try:
            for item in self.items:
                self._write(item)
        finally:
            self._mirroring = False
        self.geometryChanged.emit(self)
        self.geometryFinished.emit(self)

    def swap_axes(self):
        """Mirror the rectangle across x=y, to follow a Transpose of the images."""
        if roi_is_rotated(self):
            return  # a rotated rectangle does not map onto a simple x/y swap
        (x, y), (w, h) = self._pos, self._size
        self.set_geometry(y, x, h, w, self._angle)

    # -- lock / visibility --------------------------------------------------

    def _apply_state(self, item):
        item.translatable = not self.locked
        item.resizable = not self.locked
        item.rotatable = not self.locked
        item.setAcceptedMouseButtons(
            QtCore.Qt.MouseButton.RightButton
            if self.locked
            else QtCore.Qt.MouseButton.LeftButton | QtCore.Qt.MouseButton.RightButton
        )
        item.setVisible(self.visible)
        for handle in item.getHandles():
            handle.setVisible(self.visible and not self.locked)

    def set_locked(self, locked: bool):
        self.locked = bool(locked)
        for item in self.items:
            self._apply_state(item)

    def set_visible(self, visible: bool):
        self.visible = bool(visible)
        for item in self.items:
            self._apply_state(item)


# -- One image, one panel ----------------------------------------------------


class ImagePanel(QtWidgets.QWidget):
    """The image area for a single file: title, ImageView, crosshair."""

    def __init__(self, viewer, data: ImageData, show_histogram: bool):
        super().__init__()
        self.viewer = viewer
        self.data = data
        self.frame_index = 0
        self._transposed = False
        self._transposed_stack: np.ndarray | None = None
        self._last_cursor: tuple[int, int] | None = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(2)

        self.title_label = QtWidgets.QLabel(data.path.name)
        self.title_label.setToolTip(str(data.path))
        layout.addWidget(self.title_label)

        self.plot_item = pg.PlotItem()
        self.image_view = pg.ImageView(view=self.plot_item)
        self.image_view.ui.roiBtn.hide()
        self.image_view.ui.menuBtn.hide()
        if not show_histogram:
            # One shared colour scale, so only the first panel shows the sidebar.
            # Keep its space reserved though: if this panel's viewport were wider
            # than the first's, the aspect lock would give the two panels
            # different zooms and the images would not line up.
            histogram = self.image_view.ui.histogram
            policy = histogram.sizePolicy()
            policy.setRetainSizeWhenHidden(True)
            histogram.setSizePolicy(policy)
            histogram.hide()
        self.plot_item.setLabel("bottom", "X [pixels]")
        self.plot_item.setLabel("left", "Y [pixels]")
        self.plot_item.getViewBox().setAspectLocked(True)
        layout.addWidget(self.image_view, 1)

        self._setup_hover()

    # -- data ---------------------------------------------------------------

    def stack(self) -> np.ndarray:
        if not self._transposed:
            return self.data.stack
        if self._transposed_stack is None:
            self._transposed_stack = np.ascontiguousarray(
                self.data.stack.transpose(0, 2, 1)
            )
        return self._transposed_stack

    def set_transposed(self, transposed: bool):
        self._transposed = bool(transposed)

    def current_frame(self) -> np.ndarray:
        return self.stack()[self.frame_index]

    def image_item(self):
        return self.image_view.getImageItem()

    def set_frame(self, shared_index: int) -> bool:
        """Point this panel at *shared_index*, clamped to its own last frame.

        Returns whether the panel ran out of frames (i.e. the index was clamped).
        """
        last = self.data.n_frames - 1
        self.frame_index = int(np.clip(shared_index, 0, last))
        return shared_index > last

    def update_title(self, clamped: bool):
        text = f"{self.data.path.name}   [frame {self.frame_index}"
        if clamped:
            text += " — last"
        text += "]"
        self.title_label.setText(text)

    # -- crosshair ----------------------------------------------------------

    def _setup_hover(self):
        pen = pg.mkPen((255, 255, 0, 150), width=1)
        self._hline = pg.InfiniteLine(angle=0, movable=False, pen=pen)
        self._vline = pg.InfiniteLine(angle=90, movable=False, pen=pen)
        self._hover_text = pg.TextItem("", color=(255, 255, 255))
        self._hover_text.setAnchor((0, 1))
        for item in (self._hline, self._vline, self._hover_text):
            item.setZValue(1000)
            self.plot_item.addItem(item, ignoreBounds=True)
        self._mouse_proxy = pg.SignalProxy(
            self.plot_item.scene().sigMouseMoved, rateLimit=60, slot=self._on_mouse_moved
        )

    def _on_mouse_moved(self, event):
        pos = event[0]
        view_box = self.plot_item.getViewBox()
        if not self.plot_item.sceneBoundingRect().contains(pos):
            return
        point = view_box.mapSceneToView(pos)
        x, y = int(np.floor(point.x())), int(np.floor(point.y()))
        frame = self.current_frame()
        if not (0 <= x < frame.shape[0] and 0 <= y < frame.shape[1]):
            self._last_cursor = None
            self._hover_text.setText("")
            return
        self._last_cursor = (x, y)
        self._vline.setPos(x + 0.5)
        self._hline.setPos(y + 0.5)
        self._hover_text.setPos(x, y)
        self.refresh_hover_text()

    def refresh_hover_text(self):
        if self._last_cursor is None:
            return
        x, y = self._last_cursor
        frame = self.current_frame()
        if not (0 <= x < frame.shape[0] and 0 <= y < frame.shape[1]):
            return
        value = float(frame[x, y])
        self._hover_text.setText(f"({x}, {y})  I={value:.6g}")
        self.viewer.show_readout(self, x, y, value)


# -- ROI plot docks ---------------------------------------------------------


class _RoiDockBase(QtWidgets.QDockWidget):
    """Shared series bookkeeping for the two ROI plot docks."""

    def __init__(self, viewer, roi: SharedRoi, title: str):
        super().__init__(title, viewer)
        self.setAllowedAreas(QtCore.Qt.DockWidgetArea.AllDockWidgetAreas)
        self.viewer = viewer
        self.roi = roi

        # Dragging an ROI emits geometryChanged continuously; recomputing every
        # stack on each one would make the drag stutter, so coalesce.
        self._recompute_timer = QtCore.QTimer(self)
        self._recompute_timer.setSingleShot(True)
        self._recompute_timer.setInterval(150)
        self._recompute_timer.timeout.connect(self.recompute)

    def connect_signals(self):
        self.roi.geometryChanged.connect(self._schedule_recompute)
        self.roi.geometryFinished.connect(self._schedule_recompute)
        self.viewer.backgroundChanged.connect(self._schedule_recompute)
        self.viewer.frameChanged.connect(self.on_frame_changed)
        self.viewer.dataChanged.connect(self.recompute)

    def _schedule_recompute(self, *_args):
        self._recompute_timer.start()

    def background_per_pixel(self, panel) -> np.ndarray | None:
        """The background ROI's mean intensity per pixel, per frame, for *panel*.

        ``None`` when no background is set, or when this dock *is* the background
        ROI -- a background does not correct itself.
        """
        background = self.viewer.background_roi()
        if background is None or background is self.roi:
            return None
        return roi_series(
            panel.stack(), background, background.item_for(panel)
        )["mean"]

    def panel_series(self, panel) -> dict[str, np.ndarray]:
        """Every metric for this ROI over one panel's stack.

        With a background ROI set, its per-pixel mean ``b`` is taken off the
        pixel-value metrics: the sum loses ``b x (pixel count)`` and mean/min/max
        each lose ``b``. CoM is deliberately left on the raw pixels -- weighting a
        centroid by background-subtracted values goes unstable wherever the
        corrected sum passes through zero.
        """
        series = roi_series(panel.stack(), self.roi, self.roi.item_for(panel))

        background = self.background_per_pixel(panel)
        if background is not None:
            n_pixels = roi_pixel_count(
                self.roi,
                self.viewer.frame_shape(),
                self.roi.item_for(panel),
                panel.current_frame(),
            )
            series["sum"] = series["sum"] - background * n_pixels
            for key in ("mean", "min", "max"):
                series[key] = series[key] - background

        series["frame"] = np.arange(panel.data.n_frames, dtype=float)
        if panel.data.timestamps is not None:
            series["time"] = panel.data.timestamps
        series.update(panel.data.channels)
        return series

    def normalization_key(self) -> str:
        return combo_key(self.norm_select, NO_NORMALIZATION)

    def normalize(self, values: np.ndarray, series: dict) -> np.ndarray:
        """Divide *values* by the selected normalization channel, point by point.

        Points where the divisor is zero become NaN rather than infinity, so they
        simply leave a gap in the curve.
        """
        key = self.normalization_key()
        if key == NO_NORMALIZATION or key not in series:
            return values
        divisor = np.asarray(series[key], dtype=float)
        length = min(len(values), len(divisor))
        values, divisor = values[:length], divisor[:length]
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(divisor != 0.0, values / divisor, np.nan)

    def normalize_scalar(self, values: np.ndarray, series: dict, index: int) -> np.ndarray:
        """Divide *values* by the normalization channel at one frame.

        Used for the single-frame projections, where the curve runs over pixels
        but the divisor is a per-frame quantity.
        """
        key = self.normalization_key()
        if key == NO_NORMALIZATION or key not in series:
            return values
        divisor = np.asarray(series[key], dtype=float)
        if index >= len(divisor) or divisor[index] == 0.0:
            return np.full_like(values, np.nan)
        return values / divisor[index]

    def normalization_note(self) -> str:
        key = self.normalization_key()
        if key == NO_NORMALIZATION:
            return ""
        return f" / {DISPLAY_NAMES.get(key, key)}"

    def refresh_norm_combo(self, keys: list[str]):
        current = combo_key(self.norm_select, NO_NORMALIZATION)
        self.norm_select.blockSignals(True)
        self.norm_select.clear()
        self.norm_select.addItem(NO_NORMALIZATION, NO_NORMALIZATION)
        for key in keys:
            add_metric_item(self.norm_select, key)
        set_combo_key(self.norm_select, current)
        if self.norm_select.currentIndex() < 0:
            self.norm_select.setCurrentIndex(0)
        self.norm_select.blockSignals(False)

    def background_note(self) -> str:
        """Suffix for axis labels when a background is being subtracted."""
        background = self.viewer.background_roi()
        if background is None:
            return ""
        if background is self.roi:
            return "  [this is the background ROI]"
        return f"  [− {background.name} background]"

    # -- frame scrubber ----------------------------------------------------

    def build_frame_slider(self) -> QtWidgets.QWidget:
        """A frame scrubber for the dock, shown in every mode."""
        widget = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout(widget)
        row.setContentsMargins(0, 0, 0, 0)
        self.frame_label = QtWidgets.QLabel("frame 0")
        self.frame_slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.frame_slider.setMinimum(0)
        self.frame_slider.setMaximum(max(self.viewer.max_frames() - 1, 0))
        self.frame_slider.valueChanged.connect(
            lambda value: self.viewer.set_frame(int(value))
        )
        row.addWidget(self.frame_label)
        row.addWidget(self.frame_slider, 1)
        return widget

    def sync_frame_slider(self, index: int):
        self.frame_slider.blockSignals(True)
        self.frame_slider.setMaximum(max(self.viewer.max_frames() - 1, 0))
        self.frame_slider.setValue(index)
        self.frame_slider.blockSignals(False)
        self.frame_label.setText(f"frame {index}")

    def refresh_combo(self, combo: QtWidgets.QComboBox, keys: list[str]):
        """Rebuild *combo* from *keys*, keeping the current selection if it survives."""
        current = combo_key(combo)
        combo.blockSignals(True)
        combo.clear()
        for key in keys:
            add_metric_item(combo, key)
        set_combo_key(combo, current)
        if combo.currentIndex() < 0:
            combo.setCurrentIndex(0)
        combo.blockSignals(False)

    # Subclass hooks.
    def recompute(self):
        raise NotImplementedError

    def refresh_options(self):
        raise NotImplementedError

    def update_plot(self):
        raise NotImplementedError

    def on_frame_changed(self, index: int):
        raise NotImplementedError

    def retitle(self):
        raise NotImplementedError

    def closeEvent(self, event):
        self.viewer.forget_roi_dock(self)
        super().closeEvent(event)


class RoiPlotDock(_RoiDockBase):
    """The "ROI plot": one curve per open file, so scans can be compared.

    *Offset* stacks the curves -- curve *i* is drawn at ``y + i * offset`` -- which
    separates traces that would otherwise sit on top of each other. 0 overlays them.
    """

    def __init__(self, viewer, roi: SharedRoi):
        super().__init__(viewer, roi, f"ROI plot - {roi.name}")

        container = QtWidgets.QWidget(self)
        layout = QtWidgets.QVBoxLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        self.stats_label = QtWidgets.QLabel("ROI stats: -")
        self.stats_label.setWordWrap(True)
        layout.addWidget(self.stats_label)

        mode_row = QtWidgets.QHBoxLayout()
        mode_row.addWidget(QtWidgets.QLabel("Mode:"))
        self.radio_stack = QtWidgets.QRadioButton("Stack")
        self.radio_single = QtWidgets.QRadioButton("Single frame")
        self.radio_stack.setChecked(True)
        self._mode_group = QtWidgets.QButtonGroup(self)
        self._mode_group.addButton(self.radio_stack, 0)
        self._mode_group.addButton(self.radio_single, 1)
        mode_row.addWidget(self.radio_stack)
        mode_row.addWidget(self.radio_single)
        mode_row.addStretch(1)
        layout.addLayout(mode_row)

        axis_row = QtWidgets.QHBoxLayout()
        self.x_label = QtWidgets.QLabel("X:")
        self.x_select = QtWidgets.QComboBox()
        self.y_select = QtWidgets.QComboBox()
        axis_row.addWidget(self.x_label)
        axis_row.addWidget(self.x_select, 1)
        axis_row.addWidget(QtWidgets.QLabel("Y:"))
        axis_row.addWidget(self.y_select, 1)
        axis_row.addWidget(QtWidgets.QLabel("Norm:"))
        self.norm_select = QtWidgets.QComboBox()
        self.norm_select.setToolTip(
            "Divide Y by this, point by point (after any background subtraction).\n"
            "Same choices as the X and Y axes."
        )
        axis_row.addWidget(self.norm_select, 1)

        axis_row.addWidget(QtWidgets.QLabel("Offset:"))
        self.offset_spin = QtWidgets.QDoubleSpinBox()
        self.offset_spin.setRange(-1e12, 1e12)
        self.offset_spin.setDecimals(3)
        self.offset_spin.setFixedWidth(120)
        self.offset_spin.setKeyboardTracking(False)
        self.offset_spin.setToolTip(
            "Separate the per-file curves: curve i is drawn at y + i × offset.\n"
            "0 overlays them again."
        )
        axis_row.addWidget(self.offset_spin)
        self.auto_offset_btn = QtWidgets.QPushButton("Auto")
        self.auto_offset_btn.setToolTip("Pick an offset from the curves' own span.")
        axis_row.addWidget(self.auto_offset_btn)
        layout.addLayout(axis_row)

        self.plot_item = pg.PlotItem()
        self.plot_item.showGrid(x=True, y=True, alpha=0.3)
        self.plot_widget = pg.PlotWidget(plotItem=self.plot_item)
        layout.addWidget(self.plot_widget, 1)

        self.frame_line = pg.InfiniteLine(angle=90, movable=True, pen="c")
        self.plot_item.addItem(self.frame_line)

        # Kept visible in Single frame mode too -- that is the mode where
        # scrubbing frames matters most, since the plot shows one frame.
        layout.addWidget(self.build_frame_slider())

        self.setWidget(container)

        self.series_by_panel: list[dict[str, np.ndarray]] = []
        self.projections: list[tuple[np.ndarray, np.ndarray]] = []

        self.recompute()
        set_combo_key(self.x_select, "frame")
        set_combo_key(self.y_select, "sum")

        self.connect_signals()
        self.x_select.currentIndexChanged.connect(lambda _: self.update_plot())
        self.y_select.currentIndexChanged.connect(lambda _: self.update_plot())
        self.norm_select.currentIndexChanged.connect(lambda _: self.update_plot())
        self.offset_spin.valueChanged.connect(lambda _: self.update_plot())
        self.auto_offset_btn.clicked.connect(self._auto_offset)
        self.radio_stack.toggled.connect(lambda _: self._on_mode_changed())
        self.frame_line.sigPositionChanged.connect(self._on_line_moved)
        self._on_mode_changed()

    # -- data ---------------------------------------------------------------

    def recompute(self):
        self.series_by_panel = [self.panel_series(p) for p in self.viewer.panels]
        self.refresh_options()
        self.update_plot()

    def refresh_options(self):
        keys = self.viewer.metric_keys()
        self.refresh_combo(self.x_select, keys)
        self.refresh_combo(
            self.y_select, SINGLE_FRAME_KEYS if self._single_frame() else keys
        )
        self.refresh_norm_combo(keys)

    def _single_frame(self) -> bool:
        return self.radio_single.isChecked()

    def _on_mode_changed(self):
        single = self._single_frame()
        self.x_label.setVisible(not single)
        self.x_select.setVisible(not single)
        self.frame_line.setVisible(not single)
        self.refresh_combo(
            self.y_select,
            SINGLE_FRAME_KEYS if single else self.viewer.metric_keys(),
        )
        self.refresh_norm_combo(self.viewer.metric_keys())
        set_combo_key(self.y_select, "proj_x" if single else "sum")
        self.update_plot()

    def _compute_projections(self):
        self.projections = []
        for panel in self.viewer.panels:
            sub = roi_sub_frame(
                panel.current_frame(), self.roi, self.roi.item_for(panel)
            )
            if sub is None:
                self.projections.append((np.zeros(1), np.zeros(1)))
            else:
                # Frames are indexed [x, y]: summing over y leaves a profile along X.
                self.projections.append(
                    (sub.sum(axis=1).astype(float), sub.sum(axis=0).astype(float))
                )

    # -- curves --------------------------------------------------------------

    def _curves(self) -> list[tuple[str, np.ndarray, np.ndarray]]:
        """``(label, x, y)`` for each file, before any offset is applied."""
        curves = []
        if self._single_frame():
            self._compute_projections()
            y_key = combo_key(self.y_select, "proj_x")
            for panel, series, (proj_x, proj_y) in zip(
                self.viewer.panels, self.series_by_panel, self.projections
            ):
                y = proj_x if y_key == "proj_x" else proj_y
                # The projection runs over pixels, so the divisor is this
                # frame's single value rather than a whole series.
                y = self.normalize_scalar(y, series, panel.frame_index)
                curves.append((panel.data.path.name, np.arange(len(y), dtype=float), y))
            return curves

        x_key = combo_key(self.x_select, "frame")
        y_key = combo_key(self.y_select, "sum")
        for panel, series in zip(self.viewer.panels, self.series_by_panel):
            if x_key not in series or y_key not in series:
                continue  # this file does not record that channel
            x = np.asarray(series[x_key], dtype=float)
            y = self.normalize(np.asarray(series[y_key], dtype=float), series)
            length = min(len(x), len(y))
            curves.append((panel.data.path.name, x[:length], y[:length]))
        return curves

    def _auto_offset(self):
        spans = [float(np.ptp(y)) for _label, _x, y in self._curves() if y.size]
        step = max(spans) * 1.1 if spans else 0.0
        if step > 0:
            # Round to 2 significant figures so the axis reads cleanly.
            magnitude = 10 ** (np.floor(np.log10(step)) - 1)
            step = float(np.ceil(step / magnitude) * magnitude)
        self.offset_spin.setValue(step)

    def update_plot(self):
        self.update_stats_label()
        self.plot_item.clear()
        # addLegend() installs on the PlotItem, which clear() empties, so the
        # legend is rebuilt alongside the curves.
        legend = self.plot_item.addLegend(offset=(-10, 10)) if len(self.viewer.panels) > 1 else None

        curves = self._curves()
        offset = self.offset_spin.value()
        for index, (label, x, y) in enumerate(curves):
            color = FILE_COLORS[index % len(FILE_COLORS)]
            shifted = y + index * offset
            item = self.plot_item.plot(
                x, shifted, pen=pg.mkPen(color, width=1), symbol="o", symbolSize=4,
                symbolBrush=color, symbolPen=None,
            )
            if legend is not None:
                suffix = f"  (+{index * offset:g})" if offset and index else ""
                legend.addItem(item, f"{label}{suffix}")

        note = self.normalization_note() + self.background_note()
        if self._single_frame():
            y_key = combo_key(self.y_select, "proj_x")
            self.plot_item.setLabel("bottom", axis_label("pixel"))
            self.plot_item.setLabel("left", axis_label(y_key) + note)
        else:
            self.plot_item.setLabel("bottom", axis_label(combo_key(self.x_select, "frame")))
            self.plot_item.setLabel(
                "left", axis_label(combo_key(self.y_select, "sum")) + note
            )
            self.plot_item.addItem(self.frame_line)
            self._move_line_to_frame(self.viewer.frame_index())

    def update_stats_label(self):
        parts = []
        for panel in self.viewer.panels:
            sub = roi_sub_frame(
                panel.current_frame(), self.roi, self.roi.item_for(panel)
            )
            total = "-" if sub is None else f"{np.sum(sub):.6g}"
            parts.append(f"{panel.data.path.name}: {total}")
        x0, x1, y0, y1 = roi_bbox(self.roi, self.viewer.frame_shape())
        self.stats_label.setText(
            f"{self.roi.name} [x={x0} y={y0} {x1 - x0}x{y1 - y0}]  "
            f"sum at this frame — " + ",  ".join(parts)
        )

    def _x_reference(self) -> np.ndarray:
        """The first file's X values, which the frame marker is positioned on."""
        if not self.series_by_panel:
            return np.zeros(1)
        series = self.series_by_panel[0]
        return np.asarray(
            series.get(combo_key(self.x_select, "frame"), series["frame"]), dtype=float
        )

    def _move_line_to_frame(self, index: int):
        x_data = self._x_reference()
        if len(x_data) == 0:
            return
        index = int(np.clip(index, 0, len(x_data) - 1))
        self.frame_line.blockSignals(True)
        self.frame_line.setPos(float(x_data[index]))
        self.frame_line.blockSignals(False)

    # -- frame sync ----------------------------------------------------------

    def on_frame_changed(self, index: int):
        self.sync_frame_slider(index)
        if self._single_frame():
            self.update_plot()
        else:
            self._move_line_to_frame(index)
            self.update_stats_label()

    def _on_line_moved(self):
        x_data = self._x_reference()
        if len(x_data) == 0:
            return
        # Snap to the nearest sample so the line always lands on a real frame.
        self.viewer.set_frame(
            int(np.argmin(np.abs(x_data - float(self.frame_line.value()))))
        )

    def retitle(self):
        self.setWindowTitle(f"ROI plot - {self.roi.name}")


class Roi2DPlotDock(_RoiDockBase):
    """The "ROI 2D plot" for one file: a scatter point per frame, Z as colour."""

    def __init__(self, viewer, roi: SharedRoi, panel: ImagePanel):
        super().__init__(viewer, roi, f"ROI 2D plot - {roi.name} - {panel.data.path.name}")
        self.panel = panel

        container = QtWidgets.QWidget(self)
        layout = QtWidgets.QVBoxLayout(container)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        self.stats_label = QtWidgets.QLabel("ROI stats: -")
        layout.addWidget(self.stats_label)

        axis_row = QtWidgets.QHBoxLayout()
        self.x_select = QtWidgets.QComboBox()
        self.y_select = QtWidgets.QComboBox()
        self.z_select = QtWidgets.QComboBox()
        self.norm_select = QtWidgets.QComboBox()
        self.norm_select.setToolTip(
            "Divide Y and Z by this, point by point (after any background\n"
            "subtraction). Same choices as the X, Y and Z axes."
        )
        for label, combo in (
            ("X:", self.x_select),
            ("Y:", self.y_select),
            ("Z (colour):", self.z_select),
            ("Norm:", self.norm_select),
        ):
            axis_row.addWidget(QtWidgets.QLabel(label))
            axis_row.addWidget(combo, 1)
        layout.addLayout(axis_row)

        cmap_row = QtWidgets.QHBoxLayout()
        cmap_row.addWidget(QtWidgets.QLabel("Colormap:"))
        self.cmap_select = QtWidgets.QComboBox()
        for name in available_colormaps():
            self.cmap_select.addItem(name, name)
        set_combo_key(self.cmap_select, "viridis")
        cmap_row.addWidget(self.cmap_select)
        cmap_row.addSpacing(12)
        self.z_min_spin = QtWidgets.QDoubleSpinBox()
        self.z_max_spin = QtWidgets.QDoubleSpinBox()
        for label, spin in (("Z min:", self.z_min_spin), ("Z max:", self.z_max_spin)):
            spin.setRange(-1e12, 1e12)
            spin.setDecimals(6)
            spin.setFixedWidth(110)
            spin.setEnabled(False)  # enabled only when Auto is off
            cmap_row.addWidget(QtWidgets.QLabel(label))
            cmap_row.addWidget(spin)
        self.z_auto_check = QtWidgets.QCheckBox("Auto")
        self.z_auto_check.setChecked(True)
        cmap_row.addWidget(self.z_auto_check)
        cmap_row.addStretch(1)
        layout.addLayout(cmap_row)

        self.z_range_label = QtWidgets.QLabel("Z: -")
        layout.addWidget(self.z_range_label)

        self.plot_widget = pg.PlotWidget()
        self.plot_item = self.plot_widget.getPlotItem()
        self.plot_item.showGrid(x=True, y=True, alpha=0.3)
        layout.addWidget(self.plot_widget, 1)

        self.scatter = pg.ScatterPlotItem(size=8, pen=pg.mkPen(None))
        self.highlight = pg.ScatterPlotItem(
            size=14, pen=pg.mkPen("w", width=2), brush=pg.mkBrush(None)
        )
        self.plot_item.addItem(self.scatter)
        self.plot_item.addItem(self.highlight)

        layout.addWidget(self.build_frame_slider())

        self.setWidget(container)

        self.series: dict[str, np.ndarray] = {}
        self.recompute()
        set_combo_key(self.x_select, "frame")
        set_combo_key(self.y_select, "sum")
        set_combo_key(self.z_select, "max")

        self.connect_signals()
        for combo in (self.x_select, self.y_select, self.z_select,
                      self.norm_select, self.cmap_select):
            combo.currentIndexChanged.connect(lambda _: self.update_plot())
        self.z_min_spin.valueChanged.connect(lambda _: self.update_plot())
        self.z_max_spin.valueChanged.connect(lambda _: self.update_plot())
        self.z_auto_check.toggled.connect(self._on_auto_toggled)
        self.update_plot()

    def recompute(self):
        self.series = self.panel_series(self.panel)
        self.refresh_options()
        self.update_plot()

    def refresh_options(self):
        keys = [key for key in self.viewer.metric_keys() if key in self.series]
        for combo in (self.x_select, self.y_select, self.z_select):
            self.refresh_combo(combo, keys)
        self.refresh_norm_combo(keys)

    def series_for(self, key: str) -> np.ndarray:
        values = self.series.get(key)
        if values is None:
            values = self.series.get("frame", np.zeros(1))
        return np.asarray(values, dtype=float)

    def update_stats_label(self):
        sub = roi_sub_frame(
            self.panel.current_frame(), self.roi, self.roi.item_for(self.panel)
        )
        if sub is None:
            self.stats_label.setText("ROI stats: (empty)")
            return
        x0, x1, y0, y1 = roi_bbox(self.roi, self.viewer.frame_shape())
        self.stats_label.setText(
            f"{self.panel.data.path.name} — {self.roi.name} "
            f"[x={x0} y={y0} {x1 - x0}x{y1 - y0}]  "
            f"sum={np.sum(sub):.6g}  min={np.min(sub):.6g}  "
            f"max={np.max(sub):.6g}  mean={np.mean(sub):.6g}"
        )

    def _on_auto_toggled(self, checked: bool):
        self.z_min_spin.setEnabled(not checked)
        self.z_max_spin.setEnabled(not checked)
        self.update_plot()

    def _z_range(self, z_data: np.ndarray) -> tuple[float, float]:
        if self.z_auto_check.isChecked():
            finite = z_data[np.isfinite(z_data)] if z_data.size else z_data
            z_min = float(np.min(finite)) if finite.size else 0.0
            z_max = float(np.max(finite)) if finite.size else 1.0
            for spin, value in ((self.z_min_spin, z_min), (self.z_max_spin, z_max)):
                spin.blockSignals(True)
                spin.setValue(value)
                spin.blockSignals(False)
            return z_min, z_max
        return float(self.z_min_spin.value()), float(self.z_max_spin.value())

    def update_plot(self):
        self.update_stats_label()
        x_key = combo_key(self.x_select, "frame")
        y_key = combo_key(self.y_select, "sum")
        z_key = combo_key(self.z_select, "max")
        # X stays raw -- it is normally a position; Y and Z are the value axes.
        x_data = self.series_for(x_key)
        y_data = self.normalize(self.series_for(y_key), self.series)
        z_data = self.normalize(self.series_for(z_key), self.series)
        length = min(len(x_data), len(y_data), len(z_data))
        if length == 0:
            self.scatter.setData([], [])
            self.highlight.setData([], [])
            return
        x_data, y_data, z_data = x_data[:length], y_data[:length], z_data[:length]

        z_min, z_max = self._z_range(z_data)
        brushes = value_brushes(z_data, z_min, z_max, combo_key(self.cmap_select, "viridis"))
        self.scatter.setData(x=x_data, y=y_data, brush=brushes)

        index = int(np.clip(self.panel.frame_index, 0, length - 1))
        self.highlight.setData(x=[float(x_data[index])], y=[float(y_data[index])])

        mode = "auto" if self.z_auto_check.isChecked() else "manual"
        note = self.normalization_note() + self.background_note()
        self.z_range_label.setText(
            f"Z ({DISPLAY_NAMES.get(z_key, z_key)}{self.normalization_note()}): "
            f"{z_min:.6g} ... {z_max:.6g}  [{mode}]"
        )
        self.plot_item.setLabel("bottom", axis_label(x_key))
        self.plot_item.setLabel("left", axis_label(y_key) + note)

    def on_frame_changed(self, index: int):
        self.sync_frame_slider(index)
        self.update_plot()

    def retitle(self):
        self.setWindowTitle(
            f"ROI 2D plot - {self.roi.name} - {self.panel.data.path.name}"
        )


# -- Colormaps --------------------------------------------------------------

_COLORMAP_CANDIDATES = (
    "viridis", "plasma", "inferno", "magma", "turbo", "cividis",
    "gray", "hot", "jet", "CET-L9",
)
_LUT_CACHE: dict[str, np.ndarray] = {}


def get_colormap(name: str):
    """A pyqtgraph ColorMap by name, or ``None``.

    pyqtgraph only ships some of these itself -- ``gray``, ``hot`` and ``jet``
    come from matplotlib -- and ``get()`` raises rather than returning ``None``
    for a name it cannot find, so try both sources.
    """
    for source in (None, "matplotlib"):
        try:
            cmap = pg.colormap.get(name) if source is None else pg.colormap.get(name, source=source)
        except Exception:
            continue
        if cmap is not None:
            return cmap
    return None


def available_colormaps() -> list[str]:
    """The candidate colormaps this pyqtgraph build can actually supply."""
    names = [name for name in _COLORMAP_CANDIDATES if get_colormap(name) is not None]
    return names or ["viridis"]


def _lookup_table(name: str) -> np.ndarray | None:
    if name not in _LUT_CACHE:
        cmap = get_colormap(name)
        if cmap is None:
            return None
        try:
            _LUT_CACHE[name] = cmap.getLookupTable(nPts=256, alpha=False)
        except Exception:
            return None
    return _LUT_CACHE[name]


def value_brushes(values: np.ndarray, vmin: float, vmax: float, cmap_name: str) -> list:
    """One brush per value, coloured by *cmap_name* over the ``[vmin, vmax]`` range."""
    span = vmax - vmin
    values = np.nan_to_num(np.asarray(values, dtype=float), nan=vmin, posinf=vmax, neginf=vmin)
    fraction = np.clip((values - vmin) / (span if span > 0 else 1.0), 0.0, 1.0)
    lut = _lookup_table(cmap_name)
    if lut is None:
        level = (fraction * 255).astype(int)
        return [pg.mkBrush(int(v), int(v), int(v), 200) for v in level]
    index = (fraction * (len(lut) - 1)).astype(int)
    return [pg.mkBrush(int(r), int(g), int(b), 200) for r, g, b in lut[index]]


# -- The viewer window ------------------------------------------------------


class ImageViewerWindow(QtWidgets.QMainWindow):
    """One or more image stacks side by side, sharing ROIs, levels and frame."""

    frameChanged = QtCore.pyqtSignal(int)
    dataChanged = QtCore.pyqtSignal()
    # Raised when the background ROI is chosen, cleared, or dragged/resized --
    # every other ROI's numbers depend on it.
    backgroundChanged = QtCore.pyqtSignal()

    ROI_DOCK_WIDTH = 240

    def __init__(self, paths, parent=None):
        super().__init__(parent)
        if isinstance(paths, (str, Path)):
            paths = [paths]
        paths = [Path(p) for p in paths]
        if not paths:
            raise RuntimeError("No image file to open.")

        self.datasets = [ImageData(path) for path in paths]
        self._check_total_size()
        self._check_frame_shapes()

        names = ", ".join(d.path.name for d in self.datasets)
        self.setWindowTitle(f"Image viewer - {names}")
        self.resize(1250 if len(self.datasets) == 1 else 1700, 900)

        self._frame_index = 0
        self._transposed = False
        self._rois: list[SharedRoi] = []
        self._roi_docks: list[_RoiDockBase] = []
        self._hkl_windows: list = []
        self._roi_counter = 0
        self._background_roi: SharedRoi | None = None
        self._syncing_roi_table = False
        self._syncing_levels = False
        self._syncing_roi_geometry = False
        self._syncing_views = False
        self._ranged_on_show = False

        self._build_ui()
        self._restore_stored_rois()
        self._show_frame(0, auto_range=True, force_auto_levels=True)

    def _check_total_size(self):
        total = sum(d.stack.nbytes for d in self.datasets)
        if total > MAX_TOTAL_BYTES:
            listing = "\n".join(
                f"  {d.path.name}: {d.stack.nbytes / 1024 ** 3:.2f} GB" for d in self.datasets
            )
            raise RuntimeError(
                f"These {len(self.datasets)} files need "
                f"{total / 1024 ** 3:.1f} GB of memory together, more than this "
                f"viewer loads at once ({MAX_TOTAL_BYTES / 1024 ** 3:.0f} GB):\n\n"
                f"{listing}\n\nSelect fewer files."
            )

    def _check_frame_shapes(self):
        shapes = {d.frame_shape for d in self.datasets}
        if len(shapes) > 1:
            listing = "\n".join(
                f"  {d.path.name}: {d.frame_shape[0]} x {d.frame_shape[1]}"
                for d in self.datasets
            )
            raise RuntimeError(
                "These files have different frame sizes, so they cannot share "
                f"ROIs or be compared side by side:\n\n{listing}"
            )

    # -- construction -----------------------------------------------------

    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QVBoxLayout(central)

        root.addLayout(self._build_control_row())
        root.addLayout(self._build_levels_row())

        self.info_label = QtWidgets.QLabel("")
        root.addWidget(self.info_label)

        # One panel per file, side by side and resizable.
        self.panel_splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.panels: list[ImagePanel] = []
        for index, data in enumerate(self.datasets):
            panel = ImagePanel(self, data, show_histogram=(index == 0))
            self.panels.append(panel)
            self.panel_splitter.addWidget(panel)
        # Equal widths, so the aspect lock gives every panel the same zoom.
        self.panel_splitter.setSizes([1] * len(self.panels))
        self._link_panel_views()
        root.addWidget(self.panel_splitter, 1)

        root.addLayout(self._build_frame_row())

        self._build_roi_dock()
        self._build_view_menu()
        self._wire_level_signals()
        # The panels start greyscale until a colormap is pushed onto them, so
        # apply whatever the combo is showing rather than leaving them bare.
        self._apply_colormap(combo_key(self.cmap_combo, "viridis"))

        self.status = self.statusBar()
        self.status.showMessage(
            " | ".join(str(d.path) for d in self.datasets)
            if len(self.datasets) > 1
            else str(self.datasets[0].path)
        )

    def _link_panel_views(self):
        """Keep the panels showing the same region, without pyqtgraph's view link.

        ``setXLink``/``setYLink`` fight the aspect lock here: panel 0's viewport is
        narrower than the others (it carries the histogram), so each panel's
        aspect correction changes the range, the link pushes that to its
        neighbours, and they correct again -- the range runs away to millions of
        pixels. Copying the range across under a re-entrancy guard applies each
        correction exactly once instead.
        """
        for panel in self.panels:
            panel.plot_item.getViewBox().sigRangeChanged.connect(
                self._on_panel_range_changed
            )

    def _on_panel_range_changed(self, view_box, ranges):
        if self._syncing_views or len(self.panels) < 2:
            return
        self._syncing_views = True
        try:
            for panel in self.panels:
                other = panel.plot_item.getViewBox()
                if other is not view_box:
                    other.setRange(
                        xRange=ranges[0], yRange=ranges[1], padding=0, update=True
                    )
        finally:
            self._syncing_views = False

    def _build_control_row(self) -> QtWidgets.QHBoxLayout:
        row = QtWidgets.QHBoxLayout()

        self.add_roi_btn = QtWidgets.QPushButton("Add ROI")
        self.add_roi_btn.clicked.connect(lambda: self.add_roi())
        row.addWidget(self.add_roi_btn)

        self.hkl_btn = QtWidgets.QPushButton("HKL 3D…")
        self.hkl_btn.setToolTip(
            "Map this file's frames into reciprocal space, using the geometry\n"
            "recorded in the matching master file. Needs pyvista and xrayutilities."
        )
        self.hkl_btn.clicked.connect(self.open_hkl_window)
        row.addWidget(self.hkl_btn)

        row.addWidget(QtWidgets.QLabel("Colormap:"))
        self.cmap_combo = QtWidgets.QComboBox()
        for name in available_colormaps():
            self.cmap_combo.addItem(name, name)
        set_combo_key(self.cmap_combo, "viridis")
        self.cmap_combo.currentIndexChanged.connect(
            lambda _: self._apply_colormap(combo_key(self.cmap_combo, "viridis"))
        )
        row.addWidget(self.cmap_combo)

        self.log_check = QtWidgets.QCheckBox("Log scale")
        self.log_check.toggled.connect(self._on_log_toggled)
        row.addWidget(self.log_check)

        self.aspect_check = QtWidgets.QCheckBox("Aspect 1:1")
        self.aspect_check.setChecked(True)
        self.aspect_check.toggled.connect(self._on_aspect_toggled)
        row.addWidget(self.aspect_check)

        self.transpose_check = QtWidgets.QCheckBox("Transpose")
        self.transpose_check.toggled.connect(self._on_transpose_toggled)
        row.addWidget(self.transpose_check)

        reset_btn = QtWidgets.QPushButton("Reset view")
        reset_btn.clicked.connect(self._reset_view)
        row.addWidget(reset_btn)

        # Closing the ROI dock takes its buttons with it, so there has to be a
        # way back that does not live inside the dock.
        self.panels_btn = QtWidgets.QToolButton()
        self.panels_btn.setText("Panels ▾")
        self.panels_btn.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        row.addWidget(self.panels_btn)

        row.addStretch(1)
        return row

    def _build_levels_row(self) -> QtWidgets.QHBoxLayout:
        """Colour-scale controls: auto, or typed min/max in data units.

        One scale for every panel, so equal brightness means equal counts.
        """
        row = QtWidgets.QHBoxLayout()

        self.auto_levels_check = QtWidgets.QCheckBox("Auto levels")
        self.auto_levels_check.setChecked(True)
        row.addWidget(self.auto_levels_check)

        row.addWidget(QtWidgets.QLabel("Color min:"))
        self.vmin_spin = QtWidgets.QDoubleSpinBox()
        row.addWidget(self.vmin_spin)
        row.addWidget(QtWidgets.QLabel("max:"))
        self.vmax_spin = QtWidgets.QDoubleSpinBox()
        row.addWidget(self.vmax_spin)
        for spin in (self.vmin_spin, self.vmax_spin):
            spin.setRange(-1e12, 1e12)
            spin.setDecimals(3)
            spin.setSingleStep(1.0)
            spin.setFixedWidth(130)
            spin.setKeyboardTracking(False)  # apply on Enter/focus-out, not per digit
            spin.setToolTip(
                "Colour-scale limit in data units (raw counts, even with Log scale on),\n"
                "applied to every image. Typing a value switches off Auto levels."
            )

        full_btn = QtWidgets.QPushButton("Full range")
        full_btn.setToolTip("Set min/max to the combined data range of this frame.")
        full_btn.clicked.connect(self._levels_to_frame_range)
        row.addWidget(full_btn)

        self.levels_hint = QtWidgets.QLabel("")
        row.addWidget(self.levels_hint)
        row.addStretch(1)
        return row

    def _build_frame_row(self) -> QtWidgets.QHBoxLayout:
        row = QtWidgets.QHBoxLayout()
        multi = self.max_frames() > 1

        row.addWidget(QtWidgets.QLabel("Frame:"))
        self.prev_btn = QtWidgets.QPushButton("◀")
        self.prev_btn.setFixedWidth(36)
        self.prev_btn.clicked.connect(lambda: self.set_frame(self._frame_index - 1))
        row.addWidget(self.prev_btn)

        self.frame_spin = QtWidgets.QSpinBox()
        self.frame_spin.setRange(0, max(self.max_frames() - 1, 0))
        self.frame_spin.setEnabled(multi)
        self.frame_spin.valueChanged.connect(self.set_frame)
        row.addWidget(self.frame_spin)

        row.addWidget(QtWidgets.QLabel(f"of {self.max_frames() - 1}"))

        self.next_btn = QtWidgets.QPushButton("▶")
        self.next_btn.setFixedWidth(36)
        self.next_btn.clicked.connect(lambda: self.set_frame(self._frame_index + 1))
        row.addWidget(self.next_btn)

        row.addSpacing(16)
        self.play_btn = QtWidgets.QPushButton("Play")
        self.play_btn.clicked.connect(self.start_playback)
        self.pause_btn = QtWidgets.QPushButton("Pause")
        self.pause_btn.clicked.connect(self.pause_playback)
        self.pause_btn.setEnabled(False)
        self.play_btn.setEnabled(multi)
        row.addWidget(self.play_btn)
        row.addWidget(self.pause_btn)

        row.addWidget(QtWidgets.QLabel("FPS:"))
        self.fps_spin = QtWidgets.QSpinBox()
        self.fps_spin.setRange(1, 60)
        self.fps_spin.setValue(5)
        self.fps_spin.setEnabled(multi)
        self.fps_spin.valueChanged.connect(self._on_fps_changed)
        row.addWidget(self.fps_spin)

        self.replay_check = QtWidgets.QCheckBox("Auto replay")
        self.replay_check.setChecked(True)
        self.replay_check.setEnabled(multi)
        row.addWidget(self.replay_check)

        row.addStretch(1)

        self.play_timer = QtCore.QTimer(self)
        self.play_timer.timeout.connect(self._advance_frame)
        return row

    def _build_roi_dock(self):
        dock = QtWidgets.QDockWidget("ROIs", self)
        dock.setAllowedAreas(QtCore.Qt.DockWidgetArea.AllDockWidgetAreas)
        container = QtWidgets.QWidget(dock)
        container.setMinimumWidth(1)
        layout = QtWidgets.QVBoxLayout(container)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        self.roi_table = QtWidgets.QTableWidget(0, 4)
        self.roi_table.setHorizontalHeaderLabels(["ROI", "Vis", "Lock", "Bkg"])
        self.roi_table.setMinimumWidth(1)
        self.roi_table.setMinimumHeight(90)
        self.roi_table.verticalHeader().setVisible(False)
        self.roi_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.roi_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.roi_table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.roi_table.horizontalHeader().setStretchLastSection(False)
        self.roi_table.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        for column in (1, 2, 3):
            self.roi_table.horizontalHeader().setSectionResizeMode(
                column, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
            )
        self.roi_table.horizontalHeaderItem(3).setToolTip(
            "Use this ROI as the background: its mean intensity per pixel is\n"
            "subtracted from every other ROI's sum, mean, min and max."
        )
        self.roi_table.itemChanged.connect(self._on_roi_table_item_changed)
        self.roi_table.itemDoubleClicked.connect(self._on_roi_table_double_clicked)
        self.roi_table.itemSelectionChanged.connect(self._sync_roi_geometry_panel)
        layout.addWidget(self.roi_table, 1)

        layout.addWidget(self._build_roi_geometry_box())

        # Two columns of short buttons, so the dock stays narrow.
        buttons = [
            ("Add", lambda: self.add_roi()),
            ("Rename", lambda: self.rename_roi()),
            ("Delete", lambda: self.delete_roi()),
            ("Lock all", lambda: self.set_all_locked(True)),
            ("Unlock all", lambda: self.set_all_locked(False)),
            ("Plot", lambda: self.open_roi_plot()),
            ("2D plot", lambda: self.open_roi_2d_plot()),
        ]
        button_grid = QtWidgets.QGridLayout()
        button_grid.setContentsMargins(0, 0, 0, 0)
        button_grid.setSpacing(3)
        for index, (label, slot) in enumerate(buttons):
            button = QtWidgets.QPushButton(label)
            button.clicked.connect(slot)
            button.setMinimumWidth(1)
            button_grid.addWidget(button, index // 2, index % 2)
        layout.addLayout(button_grid)

        dock.setWidget(container)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
        self.roi_dock = dock
        self.resizeDocks([dock], [self.ROI_DOCK_WIDTH], QtCore.Qt.Orientation.Horizontal)

    def _build_roi_geometry_box(self) -> QtWidgets.QGroupBox:
        """Numeric position/size/angle entry for the selected ROI."""
        box = QtWidgets.QGroupBox("Geometry")
        grid = QtWidgets.QGridLayout(box)
        grid.setContentsMargins(6, 6, 6, 6)
        grid.setHorizontalSpacing(4)
        grid.setVerticalSpacing(3)

        specs = [
            ("x", "X:", 0, 0, (-99999.0, 99999.0), ""),
            ("y", "Y:", 0, 2, (-99999.0, 99999.0), ""),
            ("w", "W:", 1, 0, (0.01, 99999.0), ""),
            ("h", "H:", 1, 2, (0.01, 99999.0), ""),
            ("angle", "Angle:", 2, 0, (-360.0, 360.0), " °"),
        ]
        self.roi_geom_spins: dict[str, QtWidgets.QDoubleSpinBox] = {}
        for key, label, row, column, (low, high), suffix in specs:
            spin = QtWidgets.QDoubleSpinBox()
            spin.setRange(low, high)
            spin.setDecimals(2)
            spin.setSingleStep(1.0)
            spin.setSuffix(suffix)
            spin.setKeyboardTracking(False)
            spin.setSizePolicy(
                QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Fixed
            )
            spin.setMaximumWidth(84)
            spin.valueChanged.connect(self._on_roi_geometry_edited)
            grid.addWidget(QtWidgets.QLabel(label), row, column)
            grid.addWidget(spin, row, column + 1)
            self.roi_geom_spins[key] = spin

        self.roi_geom_hint = QtWidgets.QLabel("")
        self.roi_geom_hint.setWordWrap(True)
        self.roi_geom_hint.setMinimumWidth(1)
        box.setToolTip("Position and size in pixels; angle in degrees.")
        grid.addWidget(self.roi_geom_hint, 3, 0, 1, 4)
        return box

    def _build_view_menu(self):
        """A View menu (and Panels button) listing every dock, so a closed one
        can always be brought back."""
        self._view_menu = QtWidgets.QMenu("&View", self)
        self._view_menu.aboutToShow.connect(self._populate_view_menu)
        self.menuBar().addMenu(self._view_menu)
        self.panels_btn.setMenu(self._view_menu)
        self._populate_view_menu()

    def _live_docks(self) -> list[QtWidgets.QDockWidget]:
        """Docks that still exist.

        Not ``findChildren`` -- a plot dock closed a moment ago is deleted
        later (WA_DeleteOnClose), and until then it would still be listed and
        could be resurrected by "Show all panels".
        """
        docks = [getattr(self, "roi_dock", None)]
        docks.extend(self._roi_docks)
        return [dock for dock in docks if dock is not None]

    def _populate_view_menu(self):
        """Rebuild the menu each time -- ROI plot docks come and go."""
        self._view_menu.clear()
        docks = self._live_docks()
        for dock in docks:
            # toggleViewAction() is checkable and restores a closed dock.
            self._view_menu.addAction(dock.toggleViewAction())
        if docks:
            self._view_menu.addSeparator()
        self._view_menu.addAction("Show all panels", self.show_all_docks)

    def show_all_docks(self):
        for dock in self._live_docks():
            dock.setVisible(True)
        if getattr(self, "roi_dock", None) is not None:
            self.resizeDocks([self.roi_dock], [self.ROI_DOCK_WIDTH], QtCore.Qt.Orientation.Horizontal)

    # -- data access ------------------------------------------------------

    def max_frames(self) -> int:
        return max(d.n_frames for d in self.datasets)

    def frame_shape(self) -> tuple[int, int]:
        shape = self.panels[0].stack().shape
        return int(shape[1]), int(shape[2])

    def frame_index(self) -> int:
        return self._frame_index

    def metric_keys(self) -> list[str]:
        """Base metrics plus every channel any open file records."""
        has_time = any(d.timestamps is not None for d in self.datasets)
        keys = [key for key in METRIC_KEYS if key != "time" or has_time]
        channels: set[str] = set()
        for data in self.datasets:
            channels.update(data.channels)
        return keys + sorted(channels)

    # Kept for single-file callers and tests.
    def stack(self) -> np.ndarray:
        return self.panels[0].stack()

    def current_frame(self) -> np.ndarray:
        return self.panels[0].current_frame()

    def image_item(self):
        return self.panels[0].image_item()

    # -- display ----------------------------------------------------------

    def _display_transform(self, frame: np.ndarray) -> np.ndarray:
        if self.log_check.isChecked():
            return np.log1p(np.maximum(frame, 0.0))
        return frame

    def _histogram_levels(self):
        try:
            return self.panels[0].image_view.getHistogramWidget().item.getLevels()
        except Exception:
            return None

    # -- colour-scale levels ----------------------------------------------
    #
    # The image items and the histogram work in the *display* domain (log1p'd
    # when Log scale is on); the spin boxes stay in the *data* domain, so a
    # number typed there always means raw counts. One set of levels drives
    # every panel.

    def _level_to_display(self, value: float) -> float:
        return float(np.log1p(max(value, 0.0))) if self.log_check.isChecked() else float(value)

    def _level_to_data(self, value: float) -> float:
        return float(np.expm1(value)) if self.log_check.isChecked() else float(value)

    def _wire_level_signals(self):
        self.auto_levels_check.toggled.connect(self._on_auto_levels_toggled)
        self.vmin_spin.valueChanged.connect(self._on_level_spin_changed)
        self.vmax_spin.valueChanged.connect(self._on_level_spin_changed)
        try:
            self.panels[0].image_view.getHistogramWidget().item.sigLevelsChanged.connect(
                self._on_histogram_changed
            )
        except Exception:
            pass  # older pyqtgraph: the spin boxes still drive the panels

    def _on_histogram_changed(self, _item=None):
        if self._syncing_levels:
            return
        levels = self._histogram_levels()
        if levels is not None:
            self._apply_levels_to_panels(levels, skip_first=True)
        self._sync_level_spins()

    def _apply_levels_to_panels(self, levels, skip_first: bool = False):
        """Put one pair of display-domain levels on every panel."""
        low, high = float(levels[0]), float(levels[1])
        self._syncing_levels = True
        try:
            for index, panel in enumerate(self.panels):
                if skip_first and index == 0:
                    continue  # it is the one that just changed
                panel.image_view.setLevels(min=low, max=high)
        finally:
            self._syncing_levels = False

    def _on_auto_levels_toggled(self, checked: bool):
        if checked:
            self._show_frame(self._frame_index, force_auto_levels=True)

    def _on_level_spin_changed(self, _value=None):
        if self._syncing_levels:
            return
        if self.auto_levels_check.isChecked():
            # Typing a limit is a request for manual control.
            self.auto_levels_check.blockSignals(True)
            self.auto_levels_check.setChecked(False)
            self.auto_levels_check.blockSignals(False)
        self._apply_manual_levels()

    def _apply_manual_levels(self):
        low = self._level_to_display(self.vmin_spin.value())
        high = self._level_to_display(self.vmax_spin.value())
        if high <= low:
            high = low + 1e-9  # pyqtgraph needs a non-empty range
        self._apply_levels_to_panels((low, high))

    def _sync_level_spins(self):
        """Mirror the shared levels into the spin boxes."""
        levels = self._histogram_levels()
        if levels is None:
            return
        self._syncing_levels = True
        try:
            self.vmin_spin.setValue(self._level_to_data(levels[0]))
            self.vmax_spin.setValue(self._level_to_data(levels[1]))
        finally:
            self._syncing_levels = False

    def _frame_data_range(self) -> tuple[float, float]:
        """The combined min/max of every panel's current frame, in data units."""
        frames = [panel.current_frame() for panel in self.panels]
        return float(min(f.min() for f in frames)), float(max(f.max() for f in frames))

    def _levels_to_frame_range(self):
        low, high = self._frame_data_range()
        self.auto_levels_check.blockSignals(True)
        self.auto_levels_check.setChecked(False)
        self.auto_levels_check.blockSignals(False)
        self._syncing_levels = True
        try:
            self.vmin_spin.setValue(low)
            self.vmax_spin.setValue(high)
        finally:
            self._syncing_levels = False
        self._apply_manual_levels()

    def _update_levels_hint(self):
        low, high = self._frame_data_range()
        label = "frame data" if len(self.panels) == 1 else "all frames"
        self.levels_hint.setText(f"({label}: {low:.6g} … {high:.6g})")

    # -- frames -------------------------------------------------------------

    def _show_frame(self, index: int, auto_range: bool = False, force_auto_levels: bool = False):
        index = int(np.clip(index, 0, self.max_frames() - 1))
        self._frame_index = index

        clamped = [panel.set_frame(index) for panel in self.panels]
        displays = [self._display_transform(panel.current_frame()) for panel in self.panels]

        if force_auto_levels or self.auto_levels_check.isChecked():
            low = float(min(d.min() for d in displays))
            high = float(max(d.max() for d in displays))
            if high <= low:
                high = low + 1e-9
            levels = (low, high)
            auto_histogram = True
        else:
            # Hold the levels explicitly, otherwise stepping frames or playing
            # back would rescale the contrast out from under the user.
            levels = self._histogram_levels()
            auto_histogram = False
            if levels is None:
                levels = (float(min(d.min() for d in displays)),
                          float(max(d.max() for d in displays)))

        self._syncing_levels = True
        try:
            for panel, display in zip(self.panels, displays):
                panel.image_view.setImage(
                    display,
                    autoLevels=False,
                    autoRange=auto_range,
                    autoHistogramRange=auto_histogram,
                    levels=levels,
                )
        finally:
            self._syncing_levels = False

        for panel, was_clamped in zip(self.panels, clamped):
            panel.update_title(was_clamped)

        self.frame_spin.blockSignals(True)
        self.frame_spin.setValue(index)
        self.frame_spin.blockSignals(False)
        self.prev_btn.setEnabled(index > 0)
        self.next_btn.setEnabled(index < self.max_frames() - 1)

        self._update_info_label()
        self._update_levels_hint()
        self._sync_level_spins()
        for panel in self.panels:
            panel.refresh_hover_text()
        self.frameChanged.emit(index)

    def _refresh_image(self):
        self._show_frame(self._frame_index)

    def _update_info_label(self):
        nx, ny = self.frame_shape()
        parts = [f"{nx} x {ny} px", f"frame {self._frame_index} of {self.max_frames() - 1}"]
        if len(self.panels) == 1 and self.datasets[0].timestamps is not None:
            parts.append(f"t = {self.datasets[0].timestamps[self.panels[0].frame_index]:.3f} s")
        elif len(self.panels) > 1:
            parts.append(f"{len(self.panels)} files, ROIs shared")
        self.info_label.setText("   |   ".join(parts))

    def _apply_colormap(self, name: str):
        cmap = get_colormap(name)
        if cmap is None:
            self.statusBar().showMessage(f"Colormap {name!r} is not available.")
            return
        for panel in self.panels:
            panel.image_view.setColorMap(cmap)

    def _on_log_toggled(self, _checked: bool):
        # The display domain changes wholesale, so the old levels are meaningless.
        self._show_frame(self._frame_index, force_auto_levels=True)

    def _on_aspect_toggled(self, checked: bool):
        for panel in self.panels:
            panel.plot_item.getViewBox().setAspectLocked(checked)

    def _on_transpose_toggled(self, checked: bool):
        for panel in self.panels:
            panel.set_transposed(checked)
        for roi in self._rois:
            roi.swap_axes()
        self._show_frame(self._frame_index, auto_range=True, force_auto_levels=True)
        self._sync_roi_geometry_panel()
        self.dataChanged.emit()

    def _reset_view(self):
        self.panels[0].plot_item.getViewBox().autoRange()

    def showEvent(self, event):
        super().showEvent(event)
        if not self._ranged_on_show:
            self._ranged_on_show = True
            # Auto-ranging during __init__ happens before Qt has laid the panels
            # out, so it measures a near-zero viewport and, with the aspect lock,
            # produces a wildly wrong range. Do it once the geometry is real.
            QtCore.QTimer.singleShot(0, self._reset_view)

    def show_readout(self, panel: ImagePanel, x: int, y: int, value: float):
        prefix = "" if len(self.panels) == 1 else f"{panel.data.path.name}:  "
        self.statusBar().showMessage(
            f"{prefix}x={x}  y={y}  value={value:.6g}   frame {panel.frame_index}"
        )

    # -- frame navigation -------------------------------------------------

    def set_frame(self, index: int):
        index = int(np.clip(index, 0, self.max_frames() - 1))
        if index == self._frame_index:
            return
        self._show_frame(index)

    def start_playback(self):
        if self.max_frames() < 2:
            return
        self.play_timer.start(max(1, int(1000 / max(1, self.fps_spin.value()))))
        self.play_btn.setEnabled(False)
        self.pause_btn.setEnabled(True)

    def pause_playback(self):
        self.play_timer.stop()
        self.play_btn.setEnabled(self.max_frames() > 1)
        self.pause_btn.setEnabled(False)

    def _on_fps_changed(self, fps: int):
        if self.play_timer.isActive():
            self.play_timer.start(max(1, int(1000 / max(1, fps))))

    def _advance_frame(self):
        nxt = self._frame_index + 1
        if nxt >= self.max_frames():
            if not self.replay_check.isChecked():
                self.pause_playback()
                return
            nxt = 0
        self._show_frame(nxt)

    # -- ROIs ---------------------------------------------------------------

    def _restore_stored_rois(self):
        """Recreate the rectangles the writer saved in ``entry/data/rois``.

        These describe the detector's own ROIs rather than anything the user
        drew here, so they start locked -- nudging one by accident would
        silently stop it matching the ROI the beamline recorded. With several
        files open the first file's set is used; they are shared by all panels.
        """
        for spec in self.datasets[0].stored_rois:
            self.add_roi(
                pos=(spec["x"], spec["y"]),
                size=(spec["w"], spec["h"]),
                angle=spec["angle"],
                name=spec["name"],
                select=False,
                locked=True,
            )
        if self._rois:
            self.roi_table.selectRow(0)

    def add_roi(
        self, pos=None, size=None, angle=0.0, name=None, select=True, locked=False
    ) -> SharedRoi:
        nx, ny = self.frame_shape()
        if pos is None or size is None:
            side_x, side_y = max(8, nx // 8), max(8, ny // 8)
            pos = (nx // 2 - side_x // 2, ny // 2 - side_y // 2)
            size = (side_x, side_y)
        self._roi_counter += 1
        if name is None:
            name = f"ROI{self._roi_counter}"

        color = ROI_COLORS[(self._roi_counter - 1) % len(ROI_COLORS)]
        roi = SharedRoi(self, name, color, pos, size, angle)
        for panel in self.panels:
            roi.add_item(panel)
        roi.geometryChanged.connect(self._on_roi_geometry_changed)
        self._rois.append(roi)
        self._append_roi_row(roi)
        roi.set_locked(locked)
        self._sync_roi_row(roi)
        if select:
            self.roi_table.selectRow(len(self._rois) - 1)
        return roi

    # -- ROI table ----------------------------------------------------------

    def _append_roi_row(self, roi: SharedRoi):
        self._syncing_roi_table = True
        row = self.roi_table.rowCount()
        self.roi_table.insertRow(row)

        name_item = QtWidgets.QTableWidgetItem(roi.name)
        name_item.setForeground(QtGui.QBrush(QtGui.QColor(*roi.color)))
        self.roi_table.setItem(row, 0, name_item)

        states = ((1, roi.visible), (2, roi.locked), (3, roi is self._background_roi))
        for column, checked in states:
            item = QtWidgets.QTableWidgetItem()
            item.setFlags(
                (item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable) & ~QtCore.Qt.ItemFlag.ItemIsSelectable
            )
            item.setCheckState(QtCore.Qt.CheckState.Checked if checked else QtCore.Qt.CheckState.Unchecked)
            item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            self.roi_table.setItem(row, column, item)
        self._syncing_roi_table = False

    def _sync_roi_row(self, roi: SharedRoi):
        """Push *roi*'s name/visibility/lock state into its table row."""
        if roi not in self._rois:
            return
        row = self._rois.index(roi)
        is_background = roi is self._background_roi
        self._syncing_roi_table = True
        self.roi_table.item(row, 0).setText(
            f"{roi.name}  (bkg)" if is_background else roi.name
        )
        states = ((1, roi.visible), (2, roi.locked), (3, is_background))
        for column, checked in states:
            item = self.roi_table.item(row, column)
            if item is not None:
                item.setCheckState(QtCore.Qt.CheckState.Checked if checked else QtCore.Qt.CheckState.Unchecked)
        self._syncing_roi_table = False

    def _on_roi_table_item_changed(self, item: QtWidgets.QTableWidgetItem):
        if self._syncing_roi_table:
            return
        row, column = item.row(), item.column()
        if not (0 <= row < len(self._rois)) or column not in (1, 2, 3):
            return
        checked = item.checkState() == QtCore.Qt.CheckState.Checked
        roi = self._rois[row]
        if column == 1:
            self.set_roi_visible(roi, checked)
        elif column == 2:
            self.set_roi_locked(roi, checked)
        else:
            self.set_background_roi(roi if checked else None)

    def _on_roi_table_double_clicked(self, item: QtWidgets.QTableWidgetItem):
        if item.column() == 0:
            self.open_roi_plot(self._rois[item.row()])

    # -- ROI lock / visibility --------------------------------------------

    def set_roi_locked(self, roi: SharedRoi, locked: bool):
        """Freeze or release the ROI's position and size in every panel.

        A locked ROI keeps accepting right-clicks so its menu (and the unlock
        entry in it) stays reachable.
        """
        roi.set_locked(locked)
        self._sync_roi_row(roi)
        if roi is self.selected_roi():
            self._sync_roi_geometry_panel()

    def set_roi_visible(self, roi: SharedRoi, visible: bool):
        """Show or hide the ROI. Its plots keep updating -- they follow geometry."""
        roi.set_visible(visible)
        self._sync_roi_row(roi)

    def set_all_locked(self, locked: bool):
        for roi in self._rois:
            self.set_roi_locked(roi, locked)

    # -- background ROI -----------------------------------------------------

    def background_roi(self) -> SharedRoi | None:
        return self._background_roi

    def set_background_roi(self, roi: SharedRoi | None):
        """Designate (or clear) the ROI whose per-pixel mean is subtracted.

        Only one ROI can be the background, so picking a new one releases the
        previous. Every open plot recomputes.
        """
        if roi is self._background_roi:
            return
        previous = self._background_roi
        if previous is not None:
            previous.geometryChanged.disconnect(self._on_background_geometry_changed)
        self._background_roi = roi
        if roi is not None:
            # Moving or resizing the background changes every other ROI's
            # correction, so the plots have to follow it live.
            roi.geometryChanged.connect(self._on_background_geometry_changed)
        for other in (previous, roi):
            if other is not None:
                self._sync_roi_row(other)
        if roi is None:
            self.statusBar().showMessage("Background subtraction off.")
        else:
            self.statusBar().showMessage(
                f"{roi.name} is the background: its mean per pixel is subtracted "
                "from the other ROIs."
            )
        self.backgroundChanged.emit()

    def _on_background_geometry_changed(self, _roi):
        self.backgroundChanged.emit()

    def selected_roi(self) -> SharedRoi | None:
        row = self.roi_table.currentRow()
        if 0 <= row < len(self._rois):
            return self._rois[row]
        return None

    # -- numeric ROI geometry ---------------------------------------------

    def _sync_roi_geometry_panel(self):
        """Show the selected ROI's geometry; disable the fields while it is locked."""
        if self._syncing_roi_geometry:
            return
        roi = self.selected_roi()
        editable = roi is not None and not roi.locked
        for spin in self.roi_geom_spins.values():
            spin.setEnabled(editable)

        if roi is None:
            self.roi_geom_hint.setText("No ROI selected.")
            return
        self.roi_geom_hint.setText(
            f"{roi.name} is locked — untick Lock to edit." if roi.locked else roi.name
        )

        self._syncing_roi_geometry = True
        try:
            pos, size = roi.pos(), roi.size()
            self.roi_geom_spins["x"].setValue(float(pos.x()))
            self.roi_geom_spins["y"].setValue(float(pos.y()))
            self.roi_geom_spins["w"].setValue(float(size.x()))
            self.roi_geom_spins["h"].setValue(float(size.y()))
            self.roi_geom_spins["angle"].setValue(float(roi.angle()))
        finally:
            self._syncing_roi_geometry = False

    def _on_roi_geometry_edited(self, _value=None):
        if self._syncing_roi_geometry:
            return
        roi = self.selected_roi()
        if roi is None or roi.locked:
            return
        self._syncing_roi_geometry = True
        try:
            roi.set_geometry(
                self.roi_geom_spins["x"].value(),
                self.roi_geom_spins["y"].value(),
                self.roi_geom_spins["w"].value(),
                self.roi_geom_spins["h"].value(),
                self.roi_geom_spins["angle"].value(),
            )
        finally:
            self._syncing_roi_geometry = False

    def _on_roi_geometry_changed(self, roi: SharedRoi):
        """Keep the geometry fields in step when an ROI is dragged."""
        if self._syncing_roi_geometry:
            return
        if roi is self.selected_roi():
            self._sync_roi_geometry_panel()

    # -- ROI lifecycle ------------------------------------------------------

    def _require_roi(self, roi=None) -> SharedRoi | None:
        roi = roi or self.selected_roi()
        if roi is None:
            self.statusBar().showMessage("No ROI selected -- use 'Add ROI' first.")
        return roi

    def rename_roi(self, roi=None):
        roi = self._require_roi(roi)
        if roi is None:
            return
        name, ok = QtWidgets.QInputDialog.getText(
            self, "Rename ROI", "Name:", text=roi.name
        )
        if not ok or not name.strip():
            return
        roi.name = name.strip()
        self._sync_roi_row(roi)
        for dock in self._roi_docks:
            if dock.roi is roi:
                dock.retitle()

    def delete_roi(self, roi=None):
        roi = self._require_roi(roi)
        if roi is None:
            return
        for dock in [d for d in self._roi_docks if d.roi is roi]:
            dock.close()
        if roi is self._background_roi:
            roi.geometryChanged.disconnect(self._on_background_geometry_changed)
            self._background_roi = None
            self.backgroundChanged.emit()
        index = self._rois.index(roi)
        roi.remove_items()
        self._rois.pop(index)
        self.roi_table.removeRow(index)
        self._sync_roi_geometry_panel()

    # -- ROI docks ----------------------------------------------------------

    def _existing_roi_dock(self, roi, dock_type, panel=None):
        for dock in self._roi_docks:
            if dock.roi is not roi or not isinstance(dock, dock_type):
                continue
            if panel is None or getattr(dock, "panel", None) is panel:
                return dock
        return None

    def open_roi_plot(self, roi=None):
        """One plot for this ROI holding every file's curve, for comparison."""
        roi = self._require_roi(roi)
        if roi is None:
            return
        existing = self._existing_roi_dock(roi, RoiPlotDock)
        if existing is not None:
            existing.setVisible(True)
            existing.raise_()
            return
        self._add_roi_dock(RoiPlotDock(self, roi), QtCore.Qt.DockWidgetArea.BottomDockWidgetArea)

    def open_roi_2d_plot(self, roi=None):
        """One 2D plot per open file -- the axes are per-file, so they cannot share."""
        roi = self._require_roi(roi)
        if roi is None:
            return
        for panel in self.panels:
            existing = self._existing_roi_dock(roi, Roi2DPlotDock, panel)
            if existing is not None:
                existing.setVisible(True)
                existing.raise_()
                continue
            self._add_roi_dock(
                Roi2DPlotDock(self, roi, panel), QtCore.Qt.DockWidgetArea.BottomDockWidgetArea
            )

    def _add_roi_dock(self, dock: _RoiDockBase, area):
        dock.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._roi_docks.append(dock)
        self.addDockWidget(area, dock)
        dock.show()

    def forget_roi_dock(self, dock: _RoiDockBase):
        if dock in self._roi_docks:
            self._roi_docks.remove(dock)

    def show_roi_menu(self, roi: SharedRoi, global_pos):
        menu = QtWidgets.QMenu()
        actions = [
            ("Open ROI plot", lambda: self.open_roi_plot(roi)),
            ("Open ROI 2D plot", lambda: self.open_roi_2d_plot(roi)),
            (None, None),
            (
                "Unlock ROI" if roi.locked else "Lock ROI",
                lambda: self.set_roi_locked(roi, not roi.locked),
            ),
            (
                "Hide ROI" if roi.visible else "Show ROI",
                lambda: self.set_roi_visible(roi, not roi.visible),
            ),
            (
                "Clear background ROI"
                if roi is self._background_roi
                else "Use as background ROI",
                lambda: self.set_background_roi(
                    None if roi is self._background_roi else roi
                ),
            ),
            (None, None),
            ("Rename ROI", lambda: self.rename_roi(roi)),
            ("Delete ROI", lambda: self.delete_roi(roi)),
        ]
        for label, slot in actions:
            if label is None:
                menu.addSeparator()
            else:
                menu.addAction(label, slot)
        if roi in self._rois:
            self.roi_table.selectRow(self._rois.index(roi))
        menu.exec(global_pos)

    # -- HKL 3D -------------------------------------------------------------

    def open_hkl_window(self):
        """Open the reciprocal-space map for the panel currently in front.

        hkl_viewer is imported here, not at module scope, so that this viewer
        keeps working on an installation without pyvista or xrayutilities.
        """
        try:
            from hkl_viewer import HklWindow
        except ImportError as exc:
            QtWidgets.QMessageBox.warning(
                self,
                "HKL 3D unavailable",
                "The 3D reciprocal-space map needs packages that are not "
                f"installed here:\n\n    {exc}\n\n"
                "Install them with:\n\n    pip install pyvista pyvistaqt xrayutilities",
            )
            return

        panel = self.panels[0]
        try:
            window = HklWindow(panel.data.path, parent=self)
        except (RuntimeError, OSError) as exc:
            QtWidgets.QMessageBox.critical(self, "Could not build the HKL map", str(exc))
            return
        window.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        window.destroyed.connect(lambda _=None, w=window: self._forget_hkl_window(w))
        # Qt destroys a window that nothing holds a reference to.
        self._hkl_windows.append(window)
        window.show()

    def _forget_hkl_window(self, window):
        self._hkl_windows = [w for w in self._hkl_windows if w is not window]

    # -- teardown ---------------------------------------------------------

    def closeEvent(self, event):
        self.play_timer.stop()
        super().closeEvent(event)


def main():
    if len(sys.argv) < 2:
        print(f"usage: {Path(sys.argv[0]).name} <image_file.h5> [more.h5 ...]")
        return 2
    app = QtWidgets.QApplication(sys.argv)
    try:
        window = ImageViewerWindow([Path(arg) for arg in sys.argv[1:]])
    except (RuntimeError, OSError) as exc:
        QtWidgets.QMessageBox.critical(None, "Could not open image file", str(exc))
        return 1
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
