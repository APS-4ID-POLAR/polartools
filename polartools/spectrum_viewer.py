"""A fluorescence-spectrum viewer for the Vortex files beside a master file.

The softglue ("sg") edition. Identical to :mod:`spectrum_viewer` apart from
what the "Combine" box does and what the monitor list offers -- see **The sgz
layout** below. It is a copy and not a subclass because almost every change is
inside a long method rather than at a seam one could override, and because the
original has to keep running untouched while this one is tried on live data.

The multi-element silicon-drift detector writes its MCA spectra the way the
area detectors write their frames: not into ``*_master.hdf`` at all, but into a
file of its own in a subdirectory of the master's, listed in
``entry/instrument/bluesky/metadata/detectors_file_relative_path``::

    vortex: vortex/scan_000055.h5

The cube inside is ``entry/data/data``, shaped ``(frames, elements, bins)``::

    entry/data/data   (800, 7, 4096)  uint32, gzip, chunks (1, 7, 4096)

Two things about that shape are not constants and must be read from the file.
The element count is one: ``/home/beams/POLAR/fluo_data/vortex/scan_000122.h5``,
from the 2025-1 run, is ``(240, 4, 4096)`` and ``uint16`` with no compression.
The dead-time bookkeeping is the other, and it is worse than a shape change --
the names *and the channel numbering* differ between the two vintages::

    2026-3   entry/instrument/NDAttributes/CHAN1DTFactor .. CHAN7DTFactor   1-based
    2025-1   entry/instrument/NDAttributes/ICR_0/OCR_0   .. ICR_3/OCR_3     0-based

so both are probed for, in that order, and a file carrying neither simply does
not offer the correction. ``corrected = raw * DTFactor`` is the XIA convention
and is confirmed rather than assumed here: ``(1 - 1/DTFactor) * 100`` reproduces
the ``CHAN{n}DTPercent`` the file records beside it to within rounding (1.0368
against 3.5448%). Where only ICR and OCR are recorded the same factor is
``ICR / OCR``.

Why the correction is on by default, against this package's general habit of
showing what was recorded. On scan 55 element 3 runs 8 to 33 percent dead while
element 6 runs 0.2 percent. Draw seven raw curves in a window whose whole
purpose is comparing elements and the plot says "element 3 is weak", which is a
wrong conclusion about the hardware stated confidently. The checkbox turns it
off; the default is the one that is not misleading.

**Finding the file.** The recorded path is tried *second*, which is the reverse
of :func:`image_viewer.image_path_for_master` and needs saying why. In
``2026-3/Jeremy_26-3/data/DefaultSample`` there are 97 masters and 66 vortex
files, and the recorded path resolves for **14** of them. From scan 76 the
IOC's file template drifted out of step with the writer::

    master says   vortex/scan76.hdf5
    disk holds    vortex/scan_000076.h5

for 52 masters, with ``detectors_file_full_path`` carrying the same wrong name.
The ``<master dir>/vortex/<stem minus "_master">.h5`` convention resolves all
66. So the convention leads and the recorded line is the fallback, for the
writers that do not follow it. Both are resolved against the master's *own*
directory for the reason :func:`image_viewer.image_path_for_master` gives: an
absolute path does not merely fail after the data is copied, it succeeds,
opening the original instead of the copy that was clicked.

**The cube is never loaded whole.** 8000 frames is 917 MB as uint32, and
several scans overlay in one window. Worse, re-reading it is not something a
frame selector can do per keystroke: it is gzip-compressed, so inflating it
costs ~1.3 s however it is sliced (64, 256 or 1024 frames at a time, or in one
bulk read -- all within 10% of each other, because the inflate is the cost and
not the I/O).

So one streamed pass at open builds a **cumulative sum per block of frames**:
with ``B = ceil(n_frames / PYRAMID_BLOCKS)``, ``_cum[i]`` is the sum of
everything before block ``i``. Any frame range is then ``_cum[j] - _cum[i]``
plus at most two short ragged ends. Measured on the 8000-frame worst case:
building both pyramids 2.2 s against 1.4 s for a plain total -- near enough
free -- 29 MB held per pyramid instead of 917 MB, and an arbitrary range
answered in **81 us** instead of 1.3 s. That is what makes dragging the frame
spin boxes possible at all.

The accumulator is ``float64``, deliberately and not idly: numpy promotes a
``uint32`` sum to ``uint64``, and ``int64 += uint64`` resolves to ``float64``,
which ``+=`` then refuses under the ``same_kind`` casting rule. An ``int64``
accumulator raises ``UFuncTypeError`` on the first block. float64 takes any
source dtype without thinking, is exact for integers to 2**53 against a maximum
bin count here of 2.4e6, and is what dead-time correction and the plot both
want anyway.

**The energy axis is nominal.** Nothing in either file records a calibration.
The defaults here (10 eV/bin, 0 eV offset) come from reading known lines off
scan 55 -- Fe Ka at bin 635, the elastic peak at bin 1117, consistent over a
5 keV span -- so they are close, but they are an inference and the axis says
so. Type over them when the real gain is known.

**Frames are not scan points, but the detector records which is which.**
Division is the wrong instrument, and scan 125 shows why: it declares
``num_points`` 100, was aborted after 9, and still wrote 800 frames. Divide by
the declaration and the note reads a tidy "8 frames per point", which is wrong
and does not look wrong; divide by the 9 rows the primary stream really holds
and the ratio is 88.9. So the count is read from the stream, never from the
metadata. The per-frame timestamps cannot rescue the mapping either --
``NDArrayEpicsTSSec`` and ``nSec`` are all zero in these files, so
:mod:`image_viewer`'s ``nearest_index``/``_scan_time_base`` machinery has
nothing to work from.

What does rescue it is ``ArrayCounter``, which counts frames *within* a point
and resets to 1 at every one, so ``cumsum(counter == 1) - 1`` is an exact
frame-to-point map however ragged the scan -- see
:meth:`SpectrumData.point_index`. On scan 125 it finds **ten** complete
80-frame blocks against the 9 rows the scan recorded: the detector finished a
point the scan was aborted before writing down. So anything pairing frames
with stream rows takes the shorter of the two and says so. The 2025-1 vintage
records no ``ArrayCounter`` and is 1:1 anyway, where even division stands in.

**The second pane answers the question the first one cannot:** how the
intensity in one emission line moves as the scan proceeds. It is off until
"Plot ROI" is pressed, and then takes a third of the width -- off because it
costs a pass over the cube that a window opened to look at a spectrum should
not pay for unasked.

The region of interest is **a centre in keV and a half-width in eV**, shown as
the grey band on the spectrum and draggable there. It opens on the tallest bin
of everything loaded, one default half-width either side; "Peak" re-centres it
on the tallest bin of what is drawn inside the bins shown, which is how you
put it on a line after zooming to one. Be clear about what the opening guess
finds: on these files the tallest peak is the **elastic** line at the incident
energy, 11.20 keV, not a fluorescence line -- Fe Ka on scan 55 is 71 percent
of it -- and the elements do not even agree, element 1 peaking at bin 635
where the sum peaks at 1120.

Held as an energy and not as bins, which is the whole reason it is in keV:
correct the nominal gain below and the ROI stays on the energy it names, so a
band that walks away from a peak is telling you the calibration you typed is
wrong. The price is that a gain edit now moves the derived bins, so it
re-reads the cube where before it only relabelled an axis.

Beside the spectrum, then, **one curve per element** -- never a sum or a mean
across them, since which element sees a line is the comparison this pane
exists to make -- each the mean of the frames on show at each scan point, and
then the mean of the scans on show. Repeated scans are the thing worth
averaging; overlaying them instead would put seven curves times three scans on
one axis. A point the frame range only partly covers is a mean over fewer
frames, which is correct and is counted in the status line, and a point it
misses entirely is NaN rather than zero. The x axis defaults to the field the
scan's own hints nominate (``energy`` here, ``Samp_x`` on the ``lup`` scans),
which is also why it is read from the master rather than assumed: the motor
scans record no ``energy`` column at all. Moving the band costs one streamed
pass over the cube (~1 s on 8000 frames), cached on the bin pair, so every
element tick, dead-time toggle, frame range and axis change afterwards
replots from memory.

The monitor the trend is divided by is chosen by **liveness, not by name**.
Every 2026-3 scan here records ``bI0`` as identically zero while ``gI0`` and
``hI0`` run, so a picker that went down a fixed list of names would open on a
divisor of zero and show a blank plot. A candidate is offered only where every
loaded scan records at least one positive reading of it -- see
:meth:`SpectrumViewerWindow._usable_monitors`.

**The sgz layout.** The softglue scans write more than one vortex frame per
scan point -- 88 of them on the ``qxscan`` runs this was written for -- and
those frames **alternate**, because the softglue toggles between two states
frame by frame. Measured on scan 597, point 0, totalled over all elements::

    14830 16614 15031 16543 15007 16831 14920 16457 14953 16555 14677 16762
    odd  (files 1,3,5..) 14827.9 +- 192.7     even (2,4,6..) 16528.8 +- 203.0
    asymmetry -5.42%; of 400 random 44/44 splits of these same frames the
    widest apart was 548 and the median 112

1701 counts apart with 200 of scatter inside each half, against a shuffle
that never managed a third of it: the split is the measurement, not noise.
Every element shows it, from -2.6% on element 1 to -8.5% on element 5, and
the whole 62-point scan comes to -4.83% (-4.88% dead-time corrected; every
number here is raw counts). So where a scan has more
than one frame per point the **Combine** box offers their difference as well
as their sum -- ``odd - even``, numbering the frames within their point from
1, which is the convention that fixes the sign and is repeated on every label
that carries it. ``ArrayCounter`` gives that numbering exactly, running 1..88
and restarting at each point; where a file lacks it, position within the point
stands in. "Odd only" and "Even only" are there because a difference is only
believable once both halves have been looked at.

This is also what ``vortex_total_el1`` and ``vortex_total_el7`` are doing in
the monitor list, where :mod:`spectrum_viewer` offers nothing of the kind.
Both are computed from the cube rather than read from the master, which is
what lets them differ between the odd frames and the even ones -- a recorded
channel has one reading per scan point and cannot. Dividing the **sum** by
one of them is still not normalisation: it is one region of this spectrum
over another, and it flattens the very trend being looked for. Dividing the
**difference** is the asymmetry, which is the quantity an XMCD measurement is
for. So they are offered, last, where a click can reach them, and never
chosen for you: ``gI0`` and ``hI0`` still open these scans.

Run standalone on one or more spectrum files, or let ``scan_viewer_sgz``'s
"Load Spectrum..." button open them for the scans it lists.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from math import ceil
from pathlib import Path

import h5py
import numpy as np

# Before matplotlib, which picks its Qt binding from whatever is already
# imported -- see the longer note in image_viewer.py.
from PyQt6 import QtCore
from PyQt6 import QtWidgets

import matplotlib
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

from .image_viewer import recorded_detector_file

SPECTRUM_DATA_PATH = "entry/data/data"
NDATTR_PATH = "entry/instrument/NDAttributes"
# Frames within the current scan point, resetting to 1 at each one. The only
# exact frame-to-point map these files carry: see the module docstring.
ARRAY_COUNTER_ATTR = "ArrayCounter"
METADATA_PATH = "entry/instrument/bluesky/metadata"
# Points taken, counted from the primary stream rather than read from
# "entry/instrument/bluesky/metadata/num_points", which is what the plan asked
# for and not what it got: see "Frames are not scan points" above.
PRIMARY_STREAM_PATH = "entry/instrument/bluesky/streams/primary"
# The detector subdirectories to look in, which is also the label whose line is
# read out of the master's detectors block. One entry today; a list because the
# next multi-element detector will want its own and nothing else need change.
SPECTRUM_SUBDIRS = ("vortex",)
SPECTRUM_SUFFIXES = (".h5", ".hdf5", ".hdf")

# Dead-time factors, newest convention first. Each entry is (template, base):
# the template is formatted with the channel number, counted from *base*, and
# a pair means the factor is the first divided by the second (ICR/OCR).
DT_SOURCES = (
    (("CHAN{n}DTFactor",), 1),
    (("ICR_{n}", "OCR_{n}"), 0),
)

# How many blocks the cumulative pyramid is cut into. 128 keeps it near 29 MB
# on the 8000-frame worst case (128 x 7 x 4096 x 8 bytes), and leaves the
# ragged ends at most 2 x 63 frames -- about 20 ms to read. The common small
# scans here (14 to 800 frames) get B = 1 to 7, so most ranges are exact with
# no ragged read at all.
PYRAMID_BLOCKS = 128
# Frames read at once while building, when B is larger than this. Bounds the
# transient: 512 x 7 x 4096 uint32 is 59 MB.
READ_CHUNK = 512
# Several files are held at once when they are overlaid, and each costs its
# pyramid rather than its cube -- so this is generous.
MAX_TOTAL_FRAMES = 200_000

# Inferred from known lines on scan 55, not read from the file. See the module
# docstring: close enough to find a peak by, not a calibration.
DEFAULT_EV_PER_BIN = 10.0
DEFAULT_EV_OFFSET = 0.0
# A gain of zero makes the bin -> keV map constant and its inverse undefined,
# which matplotlib's scale machinery turns into a garbage top axis or a raised
# exception on the next draw. The spin box simply cannot reach it.
MIN_EV_PER_BIN = 0.01

# Redraw this long after the last keystroke. Long enough that holding a spin
# box's arrow down does not queue a read per repeat, short enough to feel live.
REDRAW_DELAY_MS = 150

# The share of the total counts the opening bin range is sized to hold, and the
# fraction of that range added at each end so the edge peaks are not flush
# against the frame. See infer_bin_range for why it is counts and not peak
# height that decides.
BIN_MASS_QUANTILE = 0.99
BIN_RANGE_MARGIN = 0.05

# Headroom above the tallest visible bin, and below the smallest on a log axis:
# a factor, since a log axis has to be padded multiplicatively to be padded at
# all. Deliberately small -- this is to keep the peak off the frame, not to
# leave the decorative empty band a default autoscale would.
Y_MARGIN = 1.05
Y_LOG_MARGIN = 1.6

# What the ROI trend is divided by, best first, out of the fields the master
# actually records. Two vintages with no name in common: 2026-3 writes
# bI0/gI0/hI0, 2025-1 writes 4idgI0, so scan_viewer's REDUCE_MONITORS -- which
# knows only the first set -- would leave the 2025-1 files with nothing.
#
# All beam monitors but the last. Being on this list is not enough to be
# offered -- see _usable_monitors. bI0 stays at the head of it because it is
# the upstream monitor where one is running, but every 2026-3 scan here
# records it as identically zero, so a picker that went by name alone would
# open on a divisor of zero and show a blank trend. gI0 and hI0 are the live
# ones on those files.
#
# The monitors read straight out of the master. The order is the priority:
# whatever is first and live is what the window opens on, so gI0/hI0 keep
# these scans. The detector's own recorded vortex_total_roi1 is deliberately
# *not* here -- it is a fixed hardware window whose bounds these files record
# nowhere, so it is a divisor nobody can describe. What the monitor list
# offers instead is the two computed monitors below.
MONITOR_CANDIDATES = ("bI0", "gI0", "hI0", "4idgI0", "Time1", "Time")

# The monitors this viewer computes rather than reads: name -> the element
# each is taken from, counting from 1 as the detector does. Each is that
# element of the cube summed over a window chosen in its own monitor-ROI
# dialog.
#
# Not in MONITOR_CANDIDATES, and gated differently: nothing records them, so
# what decides whether one is offered is the cube having that element with
# counts in it, not the master having a field of that name.
#
# Named "el" and not "roi". The obvious names would have been
# vortex_total_roi1 and vortex_total_roi7, and the first of those is what
# this monitor was called until 2026-10-18 -- but in the IOC
# vortex_total_roiN means "ROI N summed over all seven elements", not
# "element N". Element 1 got away with the borrowed name because only stats
# plugin 1 reads, so the two coincide (element 1 is 99.93% of the recorded
# vortex_total_roi1 on scan 597). Element 7 would not: configure ROI 7 and
# the master grows a channel of that name meaning something else entirely.
#
# Why these two and not all seven: element 1 carries the signal, so it is the
# monitor that matches what the recorded channel always was, and element 7 is
# the one that carries almost none of it, so dividing by it is not dividing
# the line by itself -- see the warning in MonitorRoiWindow about what each of
# them costs. Both open on the same band (DEFAULT_MONITOR_BINS); what
# distinguishes them is the element, not the window. Adding another is a line
# here, but six more would bury gI0 in the list.
COMPUTED_MONITORS = {"vortex_total_el1": 1, "vortex_total_el7": 7}

# The spectrum window's title. The browser's own "4ID HDF5 File Viewer" is
# left alone; this is the window every softglue feature lives in.
WINDOW_TITLE = "4ID HDF5 Vortex softglue viewer"

# Stands in for an element number in the trend's per-curve dict when the
# elements are summed. Negative so that it cannot collide with one, and so
# that sorted() puts the single summed curve first.
SUM_KEY = -1

# The trend's X axis when no scan axis can be read, and the label it carries.
TREND_POINT_AXIS = "scan point"
TREND_FRAME_AXIS = "frame"

# What the Combine box does with the frames of a scan point, as
# (menu label, part, what the Y axis is then counting) -- "part" being the
# word range_sum and the trend both take. One table rather than three, so
# that a fifth contrast cannot be added to the menu and forgotten on an axis.
# Sum first, so the box opens on what every other viewer here shows.
COMBINE_MODES = (
    ("Sum", "all", "counts"),
    ("Odd - even", "difference", "odd - even counts"),
    ("Odd only", "odd", "odd-numbered counts"),
    ("Even only", "even", "even-numbered counts"),
)
# Modes that subtract one half of the frames from the other, and so can go
# negative. Collected rather than compared against "difference" in five
# places, since a second contrast would otherwise have to find all five.
SIGNED_PARTS = frozenset({"difference"})

# Half-width of the ROI, in eV either side of its centre. 200 eV is about one
# FWHM of the peaks in these files -- the elastic line runs 310 to 500 eV wide
# and the 2025-1 fluorescence line 150 -- so the opening window holds the peak
# and not much else: 39 to 53 percent of the counts in the displayed range,
# from 1.8 percent of its width.
DEFAULT_ROI_WIDTH_EV = 200.0

# The window both computed monitors open on, and the range their dialog
# plots, in MCA bins. Clamped to the file's bin count on use; converted to keV
# through the gain and offset in force, so a recalibration moves the ends.
#
# 2.000-18.000 keV at the default 10 eV/bin: above the zero-energy wall that
# the first few bins carry (bin 2 to 4 on these files, 5.5 to 8.7 percent of
# the tallest bin, and not signal), below the empty tail past bin ~1700.
# Measured on scan 597, it holds 96.5 percent of element 1's counts and 95.7
# percent of element 7's, so either monitor is near enough a total.
#
# A fixed band and not a peak the viewer finds. Opening on the lowest peak of
# each element put element 1's window on the 7.41 keV line being measured, on
# every 2026-3 scan -- dividing the measurement by itself and reading exactly
# 0.000% asymmetry where the raw figure is -2.048%. A band wide enough to be
# a total cannot do that to any one line.
DEFAULT_MONITOR_BINS = (200, 1800)

# The ROI the last window was left on, as bins, for the next one to open on.
# Module level because the scan browser opens one window per scan
# (scan_viewer.py:6585) and they share nothing else; the process is the
# session, which is the scope wanted -- nothing is written to disk, so a
# fresh start is a restart away.
#
# Bins and not keV, which is the opposite of how a window holds its own ROI
# (see SpectrumViewerWindow.__init__) and is deliberate in both places. Within
# a window the ROI is an energy, so recalibrating walks the band off a peak
# that moved -- the symptom worth seeing. Across windows what is wanted is the
# same *line*, and the next scan may be a different detector gain, so bins are
# the thing that carries. The two agree whenever the gain and offset are left
# alone, which is every normal session.
#
# Only ROIs a user *chose* land here: _set_roi_kev is the single funnel for the
# Peak button, the spin boxes and a drag of the band, while the opening default
# is assigned directly in _build_roi_row and does not pass through it.
_REMEMBERED_ROI_BINS: tuple[int, int] | None = None

# How far two scans' X values may differ and still be averaged together, as a
# fraction of one step. Relative and not absolute for the reason
# scan_viewer.same_grid gives: no one figure suits an energy in keV and a motor
# position in mm at once. Half a step apart is two different grids.
SAME_X_STEPS = 0.5


# -- Locating the spectrum file ---------------------------------------------


def spectrum_path_for_master(master: Path) -> Path | None:
    """The spectrum file belonging to *master*, or ``None`` if none was found.

    The ``<master dir>/vortex/<stem minus "_master">.h5`` convention is tried
    first and the path the master records only second -- the reverse of
    :func:`image_viewer.image_path_for_master`, because here the recorded value
    is wrong far more often than it is right (52 of 66 in the 2026-3 run; see
    the module docstring). Everything is resolved against the master's own
    directory either way.

    Two things keep the cost of that second chance off the common case, because
    reading the recorded value means opening the master with h5py and this
    function is asked of every row of the listing. The subdirectory is stat-ed
    first: a run that took no fluorescence at all has no ``vortex/``, so 1254
    TiSe2 masters answer "no spectrum" on one stat each instead of on an open
    each -- 5.9 s of a listing refresh that was being spent to learn nothing.
    And the convention candidates are tested as they are built rather than
    collected first, so the 67 scans that do have a spectrum never open the
    master either. Only a master with a ``vortex/`` beside it and no file of
    the expected name in it is worth asking what it thinks it wrote.

    The stat gate is also the column's definition rather than merely an
    optimisation: "a spectrum in subdirectory ``vortex``" is the question being
    answered, so no ``vortex/``, no tick, whatever the master records.
    """
    master = Path(master)

    stem = master.stem
    if stem.endswith("_master"):
        stem = stem[: -len("_master")]

    present = [d for d in SPECTRUM_SUBDIRS if (master.parent / d).is_dir()]

    for subdir in present:
        for suffix in SPECTRUM_SUFFIXES:
            candidate = master.parent / subdir / f"{stem}{suffix}"
            if candidate.is_file():
                return candidate

    for subdir in present:
        recorded = recorded_detector_file(master, subdir)
        if not recorded:
            continue
        recorded = Path(recorded)
        # pathlib's "/" throws the left side away when the right is absolute,
        # so an absolute value here would escape the master's directory and
        # defeat the policy above -- the same guard, for the same reason, as
        # image_viewer.image_path_for_master.
        if recorded.is_absolute():
            recorded = Path(recorded.name)
        candidate = master.parent / recorded
        if candidate.is_file():
            return candidate
    return None


def is_spectrum_file(path: Path) -> bool:
    """Whether *path* is itself a spectrum file (a 3-D MCA cube).

    This one opens the file, so it is asked only of a file selected directly --
    never of every row of a listing. :func:`resolve_spectrum_file` says why.
    """
    try:
        with h5py.File(path, "r") as f:
            dset = f.get(SPECTRUM_DATA_PATH)
            return isinstance(dset, h5py.Dataset) and dset.ndim == 3
    except OSError:
        return False


# Keyed by (path, mtime, size), like image_viewer's: the browser re-resolves
# the whole selection whenever it changes, and these sit on a beamline mount.
_RESOLVE_CACHE: dict[tuple, "Path | None"] = {}


def resolve_spectrum_file(path: Path) -> Path | None:
    """The spectrum file to open for *path*, master or spectrum file.

    Deliberately a ``stat`` and not an open. The browser calls this once per
    row per refresh, and over the real mount opening each candidate to check
    its dataset costs 6.2 ms against 7 us for :meth:`Path.is_file` -- 0.6 s
    added to a listing of 97 scans, and far worse cold over NFS. So the
    column's tick means "the file is there", and whether it holds a cube this
    viewer can read is settled by :class:`SpectrumData`, which raises. That is
    the same division :func:`image_viewer.image_path_for_master` makes, and for
    the same reason.

    A file selected directly is the exception: nothing else would resolve it,
    and there is one of it rather than a listing's worth. A master is never
    that exception, and the test below says so explicitly, because ``.hdf`` is
    one of :data:`SPECTRUM_SUFFIXES` and the masters are ``*_master.hdf``.
    Without the ``_master`` test every master that did not resolve would be
    opened to be asked whether it is itself a cube -- 4.4 ms a row across 1254
    TiSe2 masters, every one of which answers no. It would also be the one way
    this function could put a tick on a row that has no spectrum, since an
    image master with a 3-D ``entry/data/data`` would answer yes.
    """
    path = Path(path)
    try:
        stat = path.stat()
        key = (str(path), stat.st_mtime, stat.st_size)
    except OSError:
        key = None
    if key is not None and key in _RESOLVE_CACHE:
        return _RESOLVE_CACHE[key]

    resolved = spectrum_path_for_master(path)
    if (
        resolved is None
        and not path.stem.endswith("_master")
        and path.suffix in SPECTRUM_SUFFIXES
        and is_spectrum_file(path)
    ):
        resolved = path
    if key is not None:
        _RESOLVE_CACHE[key] = resolved
    return resolved


def master_path_for_spectrum(spectrum: Path) -> Path | None:
    """The master whose scan wrote *spectrum*, one directory up.

    ``<dir>/vortex/scan_000055.h5`` belongs to ``<dir>/scan_000055_master.hdf``.
    Only wanted for the scan point count, so a miss costs nothing but the
    "frames per point" note.
    """
    spectrum = Path(spectrum)
    candidate = spectrum.parent.parent / f"{spectrum.stem}_master.hdf"
    return candidate if candidate.is_file() else None


# Keyed by (path, mtime, size), like _RESOLVE_CACHE: a listing re-asks about
# every row whenever the directory is refreshed, and a growing file during a
# live scan falls out of the cache on its own.
_MULTI_FRAME_CACHE: dict[tuple, bool] = {}


def is_multi_frame_scan(spectrum: Path) -> bool:
    """Whether a cube holds more than one frame per scan point.

    The sgz softglue layout, and the one condition under which this edition's
    extra features mean anything: with one frame per point there are no two
    halves for the Combine box to tell apart.

    One open, and it reads the ``ArrayCounter`` and nothing else -- a vector
    of small ints beside a cube of hundreds of megabytes. That is what makes
    it affordable once per row of a directory listing, where constructing a
    :class:`SpectrumData` is not: that builds the whole summary pyramid, which
    is a streamed pass over every frame.

    The answer is :attr:`SpectrumData.multi_frame`'s, computed the same way --
    the median point's frame count, so that an aborted scan's short trailing
    point does not speak for the scan. ``False`` for the 2025-1 vintage, which
    records no counter and is one frame per point throughout, and ``False``
    for anything unreadable: this decides a colour, and a colour is not worth
    an exception.
    """
    spectrum = Path(spectrum)
    try:
        stat = spectrum.stat()
        key = (str(spectrum), stat.st_mtime, stat.st_size)
    except OSError:
        return False
    if key in _MULTI_FRAME_CACHE:
        return _MULTI_FRAME_CACHE[key]

    answer = False
    try:
        with h5py.File(spectrum, "r") as f:
            group = f.get(NDATTR_PATH)
            node = None if group is None else group.get(ARRAY_COUNTER_ATTR)
            if node is not None:
                counter = np.asarray(node[()]).ravel()
                starts = counter == 1
                # The same guard point_index uses: a file whose first frame is
                # not the start of a point cannot be cut into points at all.
                if starts.size and starts.any() and starts[0]:
                    index = np.cumsum(starts) - 1
                    answer = float(np.median(np.bincount(index))) >= 2.0
    except (OSError, AttributeError, ValueError, IndexError):
        answer = False

    _MULTI_FRAME_CACHE[key] = answer
    return answer


# -- What the master says the scan was --------------------------------------
#
# Read here rather than imported, which is the part worth explaining:
# ``scan_viewer`` already does all of this and better, in read_scan_info and
# read_stream, but it imports *this* module to open its "Load Spectrum..."
# windows, so importing it back would be a cycle. This module already opens
# the master for the point count and already owns PRIMARY_STREAM_PATH, so the
# little that the trend pane needs lands here instead: the scanned axis, the
# field names, and the values of a named field. The hint-resolution rule below
# is deliberately the same one as scan_viewer's _resolve_hint_field, and
# should be changed in step with it.


def _yaml_block_items(text: str, key: str | None) -> list[str]:
    """Items of a one-level YAML list, as bluesky's metadata writes them.

    No PyYAML, for the same reason :mod:`image_viewer` hand-parses the
    detectors block: these are two fixed shapes written by one writer, not
    arbitrary documents. ``motors`` is a tagged flat list::

        !!python/tuple
        - energy

    and the ``dimensions`` entry of ``hints`` nests one level further, each
    dimension pairing its axes with the stream they are in::

        dimensions:
        - !!python/tuple
          - - energy
          - primary

    With *key* the block under ``key:`` is read and the doubly-dashed axis
    lines are taken; without it the whole text is read as a flat list.
    """
    lines = str(text).splitlines()
    if key is not None:
        try:
            start = next(
                i for i, line in enumerate(lines) if line.strip() == f"{key}:"
            )
        except StopIteration:
            return []
        # The block ends at the next unindented key -- "monitor: ..." here.
        rest = lines[start + 1 :]
        end = next(
            (
                i
                for i, line in enumerate(rest)
                if line[:1] not in (" ", "-", "") and ":" in line
            ),
            len(rest),
        )
        lines = rest[:end]
        # "- - energy" is the axis; "- primary" beside it is the stream.
        return [
            line.strip()[4:].strip()
            for line in lines
            if line.strip().startswith("- - ")
        ]
    return [
        line.strip()[2:].strip()
        for line in lines
        if line.strip().startswith("- ") and not line.strip().startswith("- - ")
    ]


def _meta_text(master, key: str) -> str:
    """One metadata value as text, or ``""``."""
    try:
        with h5py.File(master, "r") as f:
            meta = f.get(METADATA_PATH)
            node = None if meta is None else meta.get(key)
            if node is None:
                return ""
            value = node[()]
    except (OSError, AttributeError, ValueError, IndexError):
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def scan_axis_field(master: Path | None, fields: list[str]) -> str | None:
    """The field this scan stepped, out of those *fields* actually holds.

    The scan's own hints first and its ``motors`` metadata second -- the order
    ``scan_viewer.scan_x_field`` uses -- so an energy scan opens on ``energy``
    and a ``lup`` on ``Samp_x`` without either being written down here. That
    matters more than it looks: the sample-motor scans record no ``energy``
    field at all, so an axis hardcoded to energy would be dead on them.

    ``None`` when the scan names no axis these fields carry, which is also
    what a missing master gives.
    """
    if master is None:
        return None
    hints = _meta_text(master, "hints")
    names = _yaml_block_items(hints, "dimensions") + _yaml_block_items(
        _meta_text(master, "motors"), None
    )
    for name in names:
        # A motor records its readback under its own name and its demand under
        # a suffixed one; which exist depends on the motor.
        for candidate in (name, f"{name}_user_setpoint", f"{name}_setpoint"):
            if candidate in fields:
                return candidate
    return None


def primary_fields(master: Path | None) -> list[str]:
    """Names in the master's primary stream that carry a 1-D ``value``."""
    if master is None:
        return []
    try:
        with h5py.File(master, "r") as f:
            grp = f.get(PRIMARY_STREAM_PATH)
            if grp is None:
                return []
            found = []
            for name in sorted(grp.keys()):
                node = grp[name]
                value = node.get("value") if hasattr(node, "get") else None
                if value is not None and value.ndim == 1:
                    found.append(name)
            return found
    except (OSError, AttributeError, ValueError, IndexError):
        return []


def read_primary_field(master: Path | None, name: str) -> np.ndarray | None:
    """One primary-stream field as float64, or ``None`` if it is not there."""
    if master is None or not name:
        return None
    try:
        with h5py.File(master, "r") as f:
            grp = f.get(PRIMARY_STREAM_PATH)
            node = None if grp is None else grp.get(name)
            value = None if node is None else node.get("value")
            if value is None:
                return None
            return np.asarray(value[()], dtype=float).ravel()
    except (OSError, AttributeError, ValueError, IndexError):
        return None


def infer_bin_range(total: np.ndarray) -> tuple[int, int]:
    """The bin range worth opening on, from where the counts actually are.

    *total* is one spectrum summed over every frame and element. The range
    returned holds :data:`BIN_MASS_QUANTILE` of its counts, centred by taking
    the same tail off each end, then padded by :data:`BIN_RANGE_MARGIN`.

    Counts and not peak height, because the two other rules that suggest
    themselves were measured against real files and neither survived. A
    threshold at a fixed fraction of the maximum depends on how sharp the
    tallest peak happens to be: at 1% of peak it gives bins 0-1173 on scan 55
    but 0-2240 on scan 76, whose pileup peak at twice the elastic energy
    clears the bar, and a useless 292-397 on the 4-element file, which is one
    narrow line. A noise-floor rule -- the last bin standing clear of the
    median of the dead high end -- is worse still, returning 0-3400 for every
    2026 file, because those carry a real low-level continuum all the way out
    and so the floor is never reached. Count mass is indifferent to both: 97
    to 98 percent of every 2026 file's counts sit below the elastic peak and
    99.9 percent below the pileup peak, and the 4-element file concentrates
    at 297-662. One rule, five files, five sensible answers.

    The last bin is dropped before anything is measured. It is not a spectrum
    bin but the accumulator everything above the ADC range piles into, holding
    17454 counts on scan 55 while bins 4090-4094 hold none, so leaving it in
    anchors every inferred range to the right-hand edge of the axis.
    """
    body = np.asarray(total, dtype=float)
    if body.size < 2:
        return 0, max(0, body.size - 1)
    body = body[:-1]
    last = body.size - 1

    mass = body.sum()
    if mass <= 0:
        return 0, last
    cumulative = np.cumsum(body) / mass
    tail = (1.0 - BIN_MASS_QUANTILE) / 2.0
    low = int(np.searchsorted(cumulative, tail))
    high = int(np.searchsorted(cumulative, 1.0 - tail))

    pad = max(1, int(round((high - low) * BIN_RANGE_MARGIN)))
    low = max(0, low - pad)
    high = min(last, high + pad)
    # A single spike holds all the mass at one bin, and an empty-but-for-noise
    # file can put both ends in the same place: fall back to everything rather
    # than open on a degenerate axis.
    return (low, high) if high > low else (0, last)


def stack_mean(pieces: list[np.ndarray]) -> np.ndarray:
    """Mean of several 1-D arrays of unequal length, padding with NaN.

    The same function, and the same reasoning, as ``scan_viewer.stack_mean``.
    It is copied rather than imported because ``scan_viewer`` imports *this*
    module, so the dependency cannot run the other way.

    Row by row, never by interpolation: repeated scans step the same points,
    and one stopped part way is then simply short. A short scan leaves NaN
    where it did not reach, which is not the same as zero -- the mean over the
    tail is over the scans that got there, and a point no scan reached comes
    back NaN and is left off the plot.
    """
    if not pieces:
        return np.empty(0)
    rows = max(p.size for p in pieces)
    stack = np.full((len(pieces), rows), np.nan)
    for row, piece in enumerate(pieces):
        stack[row, : piece.size] = piece
    counted = np.count_nonzero(np.isfinite(stack), axis=0)
    # nansum of an all-NaN column is 0, so the where() is what keeps an
    # unreached point out of the plot instead of drawing it on the axis.
    return np.where(
        counted > 0, np.nansum(stack, axis=0) / np.maximum(counted, 1), np.nan
    )


def scan_number(path: Path) -> str:
    """The scan number in a file name: ``000597`` from ``scan_000597.h5``.

    The same function as ``scan_viewer.scan_number``, copied for the reason
    :func:`stack_mean` gives. The last run of digits wins, so a prefix's own
    digits are ignored.
    """
    stem = path.stem
    if stem.endswith("_master"):
        stem = stem[: -len("_master")]
    digits = re.findall(r"\d+", stem)
    return digits[-1] if digits else stem


def file_prefix(path: Path) -> str:
    """The part of a file name in front of the scan number: ``scan``."""
    stem = path.stem
    if stem.endswith("_master"):
        stem = stem[: -len("_master")]
    digits = list(re.finditer(r"\d+", stem))
    if not digits:
        return stem
    return stem[: digits[-1].start()].rstrip("_-. ") or stem


def scan_tags(paths: list[Path]) -> dict[str, str]:
    """A short name per path, keyed by ``str(path)``: enough to tell them apart.

    The scan number alone normally does it. Where two series in one directory
    collide on it, every name gains its series.
    """
    numbers = [scan_number(path) for path in paths]
    if len(set(numbers)) != len(numbers):
        return {
            str(path): f"{file_prefix(path)} {number}"
            for path, number in zip(paths, numbers)
        }
    return {str(path): number for path, number in zip(paths, numbers)}


def _safe_name(text) -> str:
    """Make *text* safe to drop into a suggested filename."""
    return re.sub(r"[^A-Za-z0-9.-]", "_", str(text))


def write_aligned_columns(
    path: Path, names: list[str], columns: list[np.ndarray], header: list[str]
) -> None:
    """Write *columns* side by side, space-aligned, under a commented header.

    ``scan_viewer.write_columns`` written for whitespace instead of CSV,
    because the beamline's own ``.dat`` files are SPEC-shaped -- a ``#``
    header, an ``#L`` line of column names, space-separated numbers -- and a
    file dropped beside them should read the same way. Every name here is one
    word, which is what lets the names be whitespace-separated at all.

    Each column is padded to the wider of its name and its widest value, so
    the table lines up for an eye as well as for :func:`numpy.loadtxt`.

    Columns of different lengths are squared off with blanks at the bottom. A
    value that is NaN is written ``nan`` instead: a curve normalized by a
    monitor that read zero has real holes in it, and "nothing was measured
    here" must not look like "this curve stopped earlier".
    """
    rows = max((column.size for column in columns), default=0)
    cells = [
        [
            f"{column[row]:.8g}" if row < column.size else ""
            for row in range(rows)
        ]
        for column in columns
    ]
    widths = [
        max(len(name), max((len(cell) for cell in column), default=0))
        for name, column in zip(names, cells)
    ]
    with open(path, "w") as stream:
        for line in header:
            stream.write(f"# {line}\n")
        stream.write(
            "#L "
            + "  ".join(
                name.ljust(w) for name, w in zip(names, widths)
            ).rstrip()
            + "\n"
        )
        for row in range(rows):
            stream.write(
                "   "
                + "  ".join(
                    column[row].rjust(w) for column, w in zip(cells, widths)
                ).rstrip()
                + "\n"
            )


def same_grid(first: np.ndarray, second: np.ndarray) -> bool:
    """Whether two X columns step the same points, over the length they share.

    :data:`SAME_X_STEPS` of a step is the tolerance. Shaped after
    ``scan_viewer.same_grid``; the point of it here is only to decide whether
    averaging the two is honest, so a disagreement is reported rather than
    acted on.
    """
    n = min(first.size, second.size)
    if n < 1:
        return False
    a, b = np.asarray(first[:n], float), np.asarray(second[:n], float)
    both = np.isfinite(a) & np.isfinite(b)
    if not both.any():
        return False
    step = (
        float(np.median(np.abs(np.diff(a[np.isfinite(a)])))) if n > 1 else 0.0
    )
    if step <= 0:
        # One point, or a column that does not move: fall back to exact.
        return bool(np.allclose(a[both], b[both]))
    return float(np.max(np.abs(a[both] - b[both]))) <= SAME_X_STEPS * step


def apply_bold_defaults(app) -> None:
    """Bold text everywhere this edition draws: Qt widgets and plots alike.

    Both halves are needed and neither reaches the other --
    ``QApplication.setFont`` knows nothing about matplotlib, and matplotlib's
    rcParams know nothing about Qt. Called at application start rather than at
    import, so that importing this module for a test leaves the process's
    fonts alone.

    ``axes.labelweight`` and ``axes.titleweight`` are set as well as
    ``font.weight``: matplotlib keeps separate weights for those two and they
    do not follow the general one, so without them the axis labels come back
    normal under bold ticks.
    """
    if app is not None:
        font = app.font()
        font.setBold(True)
        app.setFont(font)
    # matplotlib.rcParams and not pyplot's: this module builds Figures
    # directly, and importing pyplot here would bring its global figure
    # registry into a Qt application that manages its own canvases.
    matplotlib.rcParams.update(
        {
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
        }
    )


# Room for the frame and the two arrows, which belong to the style and are
# not in any string's width. Measured against this application's style at
# both weights; a few pixels over is whitespace, a few under eats a digit.
SPIN_CHROME_PX = 34


def fit_spin_width(spin, widest: str) -> None:
    """Size a spin box to the longest text it can ever show.

    Measured in the font the box is actually using, not set to a constant:
    this edition runs bold (:func:`apply_bold_defaults`), and a width that
    fits "from 40.950 keV" in the normal face loses the suffix in the bold
    one. Call it after :meth:`setRange`, :meth:`setPrefix` and
    :meth:`setSuffix` -- *widest* is built from those.
    """
    spin.setFixedWidth(
        spin.fontMetrics().horizontalAdvance(widest) + SPIN_CHROME_PX
    )


# -- Reading the spectrum file ----------------------------------------------


class SpectrumData:
    """One MCA cube, summarised so that any frame range is cheap to total.

    Never holds the cube. What it holds is a cumulative sum per block of
    frames, raw and dead-time-corrected, each ``(n_blocks + 1, n_elements,
    n_bins)`` float64 -- see the module docstring for the measurements that
    shape this.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        self.master = master_path_for_spectrum(self.path)

        with h5py.File(self.path, "r") as f:
            dset = f.get(SPECTRUM_DATA_PATH)
            if dset is None or not isinstance(dset, h5py.Dataset):
                raise RuntimeError(
                    f"{self.path.name} has no spectrum cube at "
                    f"'{SPECTRUM_DATA_PATH}'."
                )
            if dset.ndim != 3:
                raise RuntimeError(
                    f"'{SPECTRUM_DATA_PATH}' in {self.path.name} is "
                    f"{dset.ndim}-dimensional; a spectrum cube is "
                    "(frames, elements, bins)."
                )
            self.n_frames, self.n_elements, self.n_bins = (
                int(n) for n in dset.shape
            )
            if self.n_frames > MAX_TOTAL_FRAMES:
                raise RuntimeError(
                    f"{self.path.name} holds {self.n_frames} frames, more than "
                    f"this viewer summarises at once ({MAX_TOTAL_FRAMES})."
                )
            self._factors = self._read_dt_factors(f)
            self._counter = self._read_array_counter(f)
            self._build_pyramid(dset)

        self.n_points = self._read_num_points()
        # One ROI trend per slot, each keyed by the bins it was taken over. A
        # single pass answers every element selection and both dead-time
        # variants, so the cache turns element ticking into no file reads at
        # all; see roi_trend.
        self._trends: dict[str, tuple[tuple[int, int], tuple]] = {}
        # The odd-numbered frames' own pyramid, built on first use and kept.
        # A second streamed pass (0.98 s on scan 597's 5456 frames) that a
        # window left on "Sum" never pays; see _odd_pyramid.
        self._cum_odd_raw: np.ndarray | None = None
        self._cum_odd_dt: np.ndarray | None = None

    # -- building -----------------------------------------------------------

    def _read_dt_factors(self, f) -> np.ndarray | None:
        """Per-frame, per-element dead-time factors, or ``None`` if unrecorded.

        Both vintages' names are tried in turn (see the module docstring). A
        source counts only if every element has it -- a half-populated set
        would correct some elements and not others, which is worse than
        correcting none, since the point of the window is comparing them.
        """
        group = f.get(NDATTR_PATH)
        if group is None:
            return None
        for templates, base in DT_SOURCES:
            names = [
                [tpl.format(n=i + base) for tpl in templates]
                for i in range(self.n_elements)
            ]
            if not all(name in group for pair in names for name in pair):
                continue
            columns = []
            for pair in names:
                values = [
                    np.asarray(group[name][()], dtype=np.float64)
                    for name in pair
                ]
                if len(values) == 1:
                    factor = values[0]
                else:
                    # ICR / OCR. A frame with no counts out would divide by
                    # zero; it also has nothing to correct, so it stays at 1.
                    icr, ocr = values
                    factor = np.divide(
                        icr, ocr, out=np.ones_like(icr), where=ocr > 0
                    )
                columns.append(factor)
            factors = np.stack(columns, axis=1)
            if factors.shape[0] != self.n_frames:
                continue
            # A factor below 1 would mean more counts out than in. Seen as
            # exact zeros on frames the detector never armed for, which would
            # otherwise erase them.
            return np.where(factors >= 1.0, factors, 1.0)
        return None

    def _read_array_counter(self, f) -> np.ndarray | None:
        """Per-frame position within its scan point, or ``None``.

        The detector's own count of frames since the point began. It restarts
        at 1 at every point, which is what makes it an exact frame-to-point
        map however ragged the scan -- the only one these files carry, since
        the per-frame timestamps are all zero. The 2025-1 vintage does not
        record it.
        """
        group = f.get(NDATTR_PATH)
        node = None if group is None else group.get(ARRAY_COUNTER_ATTR)
        if node is None:
            return None
        counter = np.asarray(node[()]).ravel()
        if counter.size != self.n_frames:
            return None
        return counter

    def _build_pyramid(self, dset):
        """One streamed pass: cumulative sums per block, raw and corrected."""
        self._block = max(1, ceil(self.n_frames / PYRAMID_BLOCKS))
        n_blocks = ceil(self.n_frames / self._block)
        shape = (n_blocks + 1, self.n_elements, self.n_bins)
        self._cum_raw = np.zeros(shape, dtype=np.float64)
        self._cum_dt = (
            np.zeros(shape, dtype=np.float64)
            if self._factors is not None
            else None
        )

        for b in range(n_blocks):
            lo = b * self._block
            hi = min(lo + self._block, self.n_frames)
            raw, corrected = self._sum_slice(dset, lo, hi)
            self._cum_raw[b + 1] = self._cum_raw[b] + raw
            if self._cum_dt is not None:
                self._cum_dt[b + 1] = self._cum_dt[b] + corrected

    def _sum_slice(
        self, dset, lo: int, hi: int, mask: np.ndarray | None = None
    ):
        """Total frames ``[lo, hi)``, raw and corrected, reading in chunks.

        *mask*, where given, is a per-frame boolean over the whole scan; only
        the frames it selects are counted. It is sliced to the chunk and
        applied after the read, because the frames wanted are every other one
        and h5py would turn a fancy index into that many separate reads.
        """
        raw = np.zeros((self.n_elements, self.n_bins), dtype=np.float64)
        corrected = (
            np.zeros((self.n_elements, self.n_bins), dtype=np.float64)
            if self._factors is not None
            else None
        )
        for start in range(lo, hi, READ_CHUNK):
            stop = min(start + READ_CHUNK, hi)
            block = dset[start:stop]
            factors = (
                None if self._factors is None else self._factors[start:stop]
            )
            if mask is not None:
                keep = mask[start:stop]
                if not keep.any():
                    continue
                block = block[keep]
                if factors is not None:
                    factors = factors[keep]
            raw += block.sum(axis=0)
            if corrected is not None:
                # einsum, not block * factors[:, :, None]: the latter would
                # build a float64 copy of the chunk (59 MB) just to sum it.
                corrected += np.einsum(
                    "ijk,ij->jk", block, factors, optimize=True
                )
        return raw, corrected

    def _read_num_points(self) -> int | None:
        """Scan points actually taken, from the length of the primary stream.

        Every field of that stream carries one row per point, so the first in
        sorted order answers for all of them and the choice is only made sorted
        so that two runs agree. The declared ``num_points`` is deliberately not
        consulted; an aborted scan declares the points it meant to take.
        """
        if self.master is None:
            return None
        try:
            with h5py.File(self.master, "r") as f:
                grp = f.get(PRIMARY_STREAM_PATH)
                if grp is None:
                    return None
                for name in sorted(grp.keys()):
                    value = (
                        grp[name].get("value")
                        if hasattr(grp[name], "get")
                        else None
                    )
                    if value is not None and value.shape:
                        return int(value.shape[0])
                return None
        except (OSError, AttributeError, ValueError, IndexError):
            return None

    # -- asking ------------------------------------------------------------

    @property
    def have_dead_time(self) -> bool:
        return self._cum_dt is not None

    @property
    def frames_per_point(self) -> int | None:
        """Frames per scan point, where the two divide evenly.

        ``None`` when they do not -- an aborted scan writes a ratio that means
        nothing (scan 125: 800 frames, 100 points declared, 9 taken).
        """
        if not self.n_points or self.n_points <= 0:
            return None
        if self.n_frames % self.n_points:
            return None
        return self.n_frames // self.n_points

    @property
    def point_index(self) -> np.ndarray | None:
        """Which scan point each frame belongs to, as ``(n_frames,)`` ints.

        ``ArrayCounter`` first, because it is the only exact answer: it resets
        to 1 at every point, so ``cumsum(counter == 1) - 1`` numbers the points
        however ragged the scan. On scan 125 -- declared 100 points, aborted
        after 9, 800 frames written -- it finds ten complete 80-frame blocks,
        where dividing by either point count gives nonsense.

        Even division second, for the 2025-1 vintage, which records no counter
        and is 1:1 anyway. ``None`` when neither works, and then the caller
        plots against frame number rather than inventing a mapping.
        """
        counter = self._counter
        if counter is not None:
            starts = counter == 1
            # A file whose first frame is not the start of a point would need
            # a partial leading block this cannot size; division is no better
            # informed, but at least it does not silently mislabel it.
            if starts.any() and starts[0]:
                return np.cumsum(starts) - 1
        per_point = self.frames_per_point
        if per_point:
            return np.arange(self.n_frames) // per_point
        return None

    @property
    def multi_frame(self) -> bool:
        """Whether this scan took more than one vortex frame per point.

        The sgz layout, and the one condition under which the odd/even
        contrasts mean anything: with one frame per point there are no two
        halves to tell apart.

        Judged on the *median* point rather than the smallest. An aborted scan
        leaves a short trailing point -- scan 125 writes ten 80-frame blocks
        against the nine rows its stream recorded -- and that one ragged end
        says nothing about how the scan was run.
        """
        index = self.point_index
        if index is None or index.size == 0:
            return False
        return float(np.median(np.bincount(index))) >= 2.0

    @property
    def odd_frames(self) -> np.ndarray | None:
        """Mask of the odd-numbered frames, counting from 1 within their point.

        Frame 1 of a point is odd, frame 2 even, and so on -- the user's
        numbering, and the only thing that fixes the sign of a difference, so
        it is repeated on every label that carries one.

        ``ArrayCounter`` answers exactly, being a 1-based position that
        restarts at each point (1..88 on the sgz scans). Where a file lacks it
        the position within an evenly divided point stands in, and where
        neither is available there is no numbering to speak of and this is
        ``None`` -- which only happens on scans that are not ``multi_frame``
        anyway.
        """
        counter = self._counter
        if counter is not None:
            starts = counter == 1
            if starts.any() and starts[0]:
                return (counter % 2) == 1
        per_point = self.frames_per_point
        if per_point:
            # Position 0 is frame 1, so an even position is an odd frame.
            return (np.arange(self.n_frames) % per_point) % 2 == 0
        return None

    def roi_trend(
        self, lo: int, hi: int, slot: str = "roi"
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """Counts in bins ``[lo, hi]`` per frame, as ``(n_frames, n_elements)``.

        Raw and dead-time corrected, both from the one pass -- the correction
        is a per-frame, per-element factor, constant across bins, so summing
        the bins first and scaling after is the same number as scaling each bin
        and summing, for a fraction of the arithmetic.

        The pyramid cannot serve this: it is cumulative *per block*, so it has
        no frame interiors left to report. So this reads the cube, in
        ``READ_CHUNK`` blocks like :meth:`_sum_slice` -- 0.08 s for 800 frames,
        1.0 s for 8000 -- and the result is cached on ``(lo, hi)``. One array
        answers every element selection and both dead-time variants, which is
        what keeps ticking an element off a redraw rather than a re-read.

        *slot* is which cache the answer is kept in. Several callers want
        different ranges at the same time -- the ROI band, and one band per
        computed monitor -- and a single slot would make every redraw read
        the cube again for each of them. A slot each, rather than a dict of
        the last few ranges: a drag makes a new key per step, so an eviction
        policy would throw away the other bands' arrays on every one of
        them. Nothing is evicted, so the cost is one such array per slot in
        use -- 0.6 MB per slot on scan 597, and only for a monitor that has
        actually been divided by.
        """
        lo = max(0, min(int(lo), self.n_bins - 1))
        hi = max(lo, min(int(hi), self.n_bins - 1))
        held = self._trends.get(slot)
        if held is not None and held[0] == (lo, hi):
            return held[1]

        raw = np.zeros((self.n_frames, self.n_elements), dtype=np.float64)
        corrected = np.zeros_like(raw) if self._factors is not None else None
        with h5py.File(self.path, "r") as f:
            dset = f[SPECTRUM_DATA_PATH]
            for start in range(0, self.n_frames, READ_CHUNK):
                stop = min(start + READ_CHUNK, self.n_frames)
                block = dset[start:stop, :, lo : hi + 1].sum(axis=2)
                raw[start:stop] = block
                if corrected is not None:
                    corrected[start:stop] = block * self._factors[start:stop]

        self._trends[slot] = ((lo, hi), (raw, corrected))
        return self._trends[slot][1]

    def _odd_pyramid(self):
        """Build, once, the same cumulative pyramid over odd frames alone.

        A second streamed pass -- 0.98 s over scan 597's 5456 frames, against
        1.27 s for the pass at open -- paid the first time a contrast is asked
        for and never by a window that stays on "Sum".

        Only the odd half is stored: the even half is the whole minus this, so
        the file is read once more rather than twice and one array pair is
        kept rather than two.

        The mask has to be carried in rather than taken as block parity.
        ``PYRAMID_BLOCKS = 128`` makes the block 43 frames long on a
        5456-frame scan, which is odd, so from the second block on the
        positions within a block and the frames within a point disagree.
        """
        if self._cum_odd_raw is not None:
            return
        mask = self.odd_frames
        if mask is None:
            raise ValueError(
                f"{self.path.name} records no frame numbering within its scan "
                "points, so its frames cannot be split into odd and even."
            )
        n_blocks = self._cum_raw.shape[0] - 1
        self._cum_odd_raw = np.zeros_like(self._cum_raw)
        self._cum_odd_dt = (
            None if self._cum_dt is None else np.zeros_like(self._cum_dt)
        )
        with h5py.File(self.path, "r") as f:
            dset = f[SPECTRUM_DATA_PATH]
            for b in range(n_blocks):
                lo = b * self._block
                hi = min(lo + self._block, self.n_frames)
                raw, corrected = self._sum_slice(dset, lo, hi, mask)
                self._cum_odd_raw[b + 1] = self._cum_odd_raw[b] + raw
                if self._cum_odd_dt is not None:
                    self._cum_odd_dt[b + 1] = self._cum_odd_dt[b] + corrected

    def range_sum(
        self, first: int, last: int, dead_time: bool, part: str = "all"
    ) -> np.ndarray:
        """Total frames *first* to *last* inclusive, as ``(n_elements, n_bins)``.

        Whole blocks come out of the pyramid by subtraction; the partial block
        at each end is read from the file. Clamped rather than refused, because
        one spin-box range is shared by scans of different length.

        *part* picks which frames of each scan point are counted, numbering
        them from 1: ``"all"``, ``"odd"``, ``"even"``, or ``"difference"`` for
        odd minus even. Everything but ``"all"`` needs :attr:`odd_frames`, and
        the first such call builds :meth:`_odd_pyramid`; after that they are
        all subtractions and cost what ``"all"`` costs.
        """
        lo = max(0, int(first))
        hi = min(self.n_frames, int(last) + 1)
        if hi <= lo:
            return np.zeros((self.n_elements, self.n_bins), dtype=np.float64)
        if part not in ("all", "odd", "even", "difference"):
            raise ValueError(f"range_sum: no such part {part!r}")

        if part != "all":
            self._odd_pyramid()
        if part == "even":
            return self.range_sum(first, last, dead_time) - self.range_sum(
                first, last, dead_time, "odd"
            )
        if part == "difference":
            odd = self.range_sum(first, last, dead_time, "odd")
            return odd - (self.range_sum(first, last, dead_time) - odd)

        odd_only = part == "odd"
        use_dt = dead_time and self._cum_dt is not None
        if odd_only:
            cum = self._cum_odd_dt if use_dt else self._cum_odd_raw
        else:
            cum = self._cum_dt if use_dt else self._cum_raw
        block = self._block
        # The first block boundary at or after lo, and the last at or before hi.
        i = ceil(lo / block)
        j = hi // block

        if i >= j:
            # No whole block in between, so there is nothing to subtract.
            return self._read_range(lo, hi, use_dt, odd_only)

        total = cum[j] - cum[i]
        if lo < i * block:
            total = total + self._read_range(lo, i * block, use_dt, odd_only)
        if j * block < hi:
            total = total + self._read_range(j * block, hi, use_dt, odd_only)
        return total

    def _read_range(
        self, lo: int, hi: int, dead_time: bool, odd_only: bool = False
    ) -> np.ndarray:
        """Sum frames ``[lo, hi)`` straight from the file."""
        mask = self.odd_frames if odd_only else None
        with h5py.File(self.path, "r") as f:
            raw, corrected = self._sum_slice(
                f[SPECTRUM_DATA_PATH], lo, hi, mask
            )
        return corrected if dead_time and corrected is not None else raw


# -- The windows ------------------------------------------------------------


class MonitorRoiWindow(QtWidgets.QDialog):
    """Picks the energy window one computed monitor counts in.

    Opened when a monitor from :data:`COMPUTED_MONITORS` is chosen, because
    such a monitor is not a channel read from the master: it is one element of
    this cube, summed over a window, and a window has to be chosen somewhere.
    The master's own ``vortex_total_roi1`` is a fixed hardware region whose
    bounds these files record nowhere, so there is nothing to inherit and the
    choice has to be made here.

    **One instance per computed monitor**, each holding its own band, rather
    than one that retitles itself when the monitor changes. A single instance
    would have to re-seed its two boxes on every switch, and those writes go
    through :meth:`_set_range` while the viewer still holds the *old*
    monitor's range -- which is the value the guard there compares against, so
    the reseed would write one monitor's band into the other's. Two instances
    have no such path.

    Non-modal, and it drives the parent's redraw on every change: the way to
    use it is to drag the band and watch the trend next door move. The range
    outlives the window -- closing it keeps what was picked -- because it
    belongs to the viewer, not to this dialog.

    Shows its own element alone, over the frames the main window is showing,
    and on a log axis by default: picking a region out of a fluorescence
    spectrum means seeing the three decades under the line, which is exactly
    what a linear axis hides.
    """

    def __init__(self, parent, monitor: str):
        super().__init__(parent)
        self.viewer = parent
        self.monitor = monitor
        self.element = COMPUTED_MONITORS[monitor]
        self.setWindowTitle(f"Monitor ROI - {monitor} (element {self.element})")
        # Not setModal(False) and left at that: a dialog with only a close
        # button still takes the Qt default of deleting nothing on close, so
        # the instance -- and the range it is showing -- survives being shut.
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowType.Window)
        self.span = None
        self._syncing = False
        self._build_ui()
        self.refresh()

    # -- construction -------------------------------------------------------

    def _build_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Monitor range"))

        lo_kev, hi_kev = self._range
        self.lo_spin = QtWidgets.QDoubleSpinBox()
        self.lo_spin.setPrefix("from ")
        self.hi_spin = QtWidgets.QDoubleSpinBox()
        self.hi_spin.setPrefix("to ")
        for spin, value in ((self.lo_spin, lo_kev), (self.hi_spin, hi_kev)):
            spin.setDecimals(3)
            spin.setSingleStep(0.01)
            spin.setSuffix(" keV")
            spin.setKeyboardTracking(False)
            spin.setRange(*self._kev_limits())
            fit_spin_width(
                spin, f"{spin.prefix()}{spin.maximum():.3f}{spin.suffix()}"
            )
            # Set before connecting: these are the opening values, not an edit.
            spin.setValue(value)
            spin.valueChanged.connect(self._on_spin)
            row.addWidget(spin)

        # Two ends and not a centre with a half-width, which is how the signal
        # ROI next door is written: Full is one press away and reaches the
        # whole spectrum, and "20.475 keV +- 20475 eV" is not a region
        # anybody reads.
        lo_bin, hi_bin = self.viewer._default_monitor_bins()
        tip = (
            (
                f"The window element {self.element}'s counts are summed over to\n"
                f"make '{self.monitor}'. Drag the band or its edges, or type the\n"
                f"two ends. It opens on bins {lo_bin}-{hi_bin}, the same band on\n"
                "either element -- wide enough to be a near-total (96% of this\n"
                "element's counts) rather than a line divided by itself.\n"
                "Every region of this detector alternates with the frames, so\n"
                "normalising shrinks or inverts the odd - even contrast wherever\n"
                "this band is put. There is no quiet corner to move it to; what\n"
                "the two elements offer is a choice between two distortions,\n"
                "measured on scan 597 at the 7.41 keV line, which reads -2.05%\n"
                "raw:\n"
            )
            + (
                "  element 1 is where the signal is, so its own alternation is\n"
                "  small (-2.4% over the opening band, -2.4% whole spectrum).\n"
                "  The band does contain the line, but the line is a small part\n"
                "  of it, so the -2.05% reads as +0.37% here and +0.33% against\n"
                "  the whole spectrum.\n"
                if self.element == 1
                else "  element 7 carries almost none of the line, so dividing by\n"
                "  it is not dividing the signal by itself -- but element 7\n"
                "  alternates far harder than element 1 (-7.4% over the opening\n"
                "  band, -7.1% whole spectrum). The 7.41 keV line's -2.05% reads\n"
                "  as +6.62% here and +6.19% against the whole spectrum: the\n"
                "  sign flips. This is a second measurement, not a reference.\n"
            )
            + (
                "That is the measurement and not a fault. Pick the region you\n"
                "believe carries beam and not signal, read the rest off the\n"
                "trend, and say in the figure which monitor it was."
            )
        )
        self.lo_spin.setToolTip(tip)
        self.hi_spin.setToolTip(tip)

        full_btn = QtWidgets.QPushButton("Full")
        full_btn.setToolTip(
            "The whole spectrum -- the nearest honest thing to a 'total'\n"
            "when the recorded channel's own window is not in the file, and\n"
            "the fallback when no narrower band looks defensible. A small\n"
            f"step from what the window opens on: bins {lo_bin}-{hi_bin}\n"
            "already hold 96% of this element, so this adds the noise below\n"
            "and the empty tail above rather than a different measurement."
        )
        full_btn.clicked.connect(self._set_full)
        row.addWidget(full_btn)

        self.log_box = QtWidgets.QCheckBox("Log Y")
        self.log_box.setChecked(True)
        self.log_box.setToolTip(
            "On by default here, unlike the spectrum pane: this window is\n"
            "for choosing a region out of the background, and on a linear\n"
            "axis the background is a flat line along the bottom."
        )
        self.log_box.toggled.connect(self.refresh)
        row.addWidget(self.log_box)
        row.addStretch(1)
        layout.addLayout(row)

        self.figure = Figure(figsize=(7, 3.6), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot(111)
        # Built in refresh() and not here: every redraw clears the axes, which
        # takes the secondary axis with it. The main pane rebuilds its own the
        # same way, for the same reason.
        self.energy_axis = None
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.resize(760, 420)

    @property
    def _range(self) -> tuple[float, float]:
        """This monitor's band, which lives on the viewer and not here.

        Read through rather than cached, so that the one on the viewer stays
        the single copy: two dialogs are open at once and only the viewer
        knows which monitor the trend is actually dividing by.
        """
        return self.viewer.monitor_ranges_kev[self.monitor]

    def _kev_limits(self) -> tuple[float, float]:
        """The keV the first and last bin resolve to, as the boxes' range.

        Read from the gain and offset as they are now, so a recalibration
        moves the ends of the box with everything else.
        """
        ends = self.viewer._bin_to_kev([0, self.viewer.max_bins - 1])
        return float(ends[0]), float(ends[-1])

    # -- keeping the three views on one range -------------------------------

    def _on_spin(self, *_):
        self._set_range(self.lo_spin.value(), self.hi_spin.value())

    def _set_full(self):
        self._set_range(*self._kev_limits())

    def _on_span_select(self, low: float, high: float):
        """A drag of the band, which is drawn in bins."""
        lo, hi = int(round(low)), int(round(high))
        if hi <= lo:
            # A click, not a drag. Redraw to put the band back rather than
            # collapse the monitor onto one bin -- the same case, and the same
            # answer, as the signal ROI's own span handler.
            self.refresh()
            return
        self._set_range(
            float(self.viewer._bin_to_kev(lo)),
            float(self.viewer._bin_to_kev(hi)),
        )

    def _set_range(self, lo_kev: float, hi_kev: float):
        """The one place the monitor range changes, boxes and band included.

        Guarded like :meth:`SpectrumViewerWindow._set_roi_kev`, and for the
        same reason: the two boxes are written one after the other, so without
        it the first write would come back here pairing a new low end with the
        old high one and redraw a range nobody asked for on the way past.
        """
        lo_kev, hi_kev = float(lo_kev), float(hi_kev)
        if lo_kev > hi_kev:
            lo_kev, hi_kev = hi_kev, lo_kev
        if self._syncing or (lo_kev, hi_kev) == self._range:
            return

        self._syncing = True
        try:
            self.lo_spin.setValue(lo_kev)
            self.hi_spin.setValue(hi_kev)
            # Read back what the boxes took: they round to three decimals and
            # a drag arrives with far more precision than they show, so this
            # is what stops the trend's label naming a range nowhere on screen.
            # Written under this dialog's own key, which is the whole of what
            # keeps the two computed monitors' bands apart.
            self.viewer.monitor_ranges_kev[self.monitor] = (
                self.lo_spin.value(),
                self.hi_spin.value(),
            )
        finally:
            self._syncing = False
        self.refresh()
        self.viewer._schedule_redraw()

    def _make_span(self):
        """Rebuild the draggable band, as the spectrum pane rebuilds its own.

        The artist does not survive ``ax.clear()`` and the old selector's
        canvas callbacks do, so both halves are needed; see
        :meth:`SpectrumViewerWindow._make_roi_span`.
        """
        if self.span is not None:
            self.span.disconnect_events()
        self.span = SpanSelector(
            self.ax,
            self._on_span_select,
            "horizontal",
            interactive=True,
            drag_from_anywhere=True,
            useblit=False,
            props={"facecolor": "0.35", "alpha": 0.20},
            handle_props={"color": "0.25"},
        )
        self.span.extents = self.viewer.monitor_bins(self.monitor)

    # -- drawing ------------------------------------------------------------

    def refresh(self, *_):
        """Redraw this monitor's element over the frames the main window shows.

        Cheap enough to do on every parent redraw: the pyramid answers a frame
        range without touching the cube, so this is arithmetic on arrays
        already in memory.
        """
        viewer = self.viewer
        first, last = viewer.first_spin.value(), viewer.last_spin.value()
        dead_time = viewer.deadtime_box.isChecked()
        log = self.log_box.isChecked()

        self.ax.clear()
        drawn = 0
        for data in viewer.data:
            if (
                data.n_elements < self.element
                or min(last, data.n_frames - 1) < first
            ):
                continue
            # part="all" whatever the Combine box says: the band is being put
            # on a spectrum, and the odd and even halves of one are the same
            # spectrum at half the counts. Which frames the monitor's *value*
            # is taken over is decided in _monitor_scale, not here.
            totals = data.range_sum(first, last, dead_time)
            y = totals[self.element - 1]
            self.ax.plot(np.arange(y.size), y, lw=0.8, label=data.path.name)
            drawn += 1

        self.ax.set_xlabel("MCA bin")
        self.ax.set_ylabel(f"element {self.element} counts")
        self.ax.set_yscale("log" if log else "linear")
        lo, hi = viewer.monitor_bins(self.monitor)
        # The default band, widened to hold whatever the band is now. Pinned
        # at the default alone, pressing Full would set the band to the whole
        # spectrum and leave both of its edges off screen, with nothing on the
        # plot to say what had happened. Widening keeps every band visible and
        # still opens on the 1600 bins that carry the counts instead of the
        # 4096 the file is long. refresh() runs on every parent redraw, so
        # this decides which range is restored, not whether the toolbar's own
        # zoom can be used in between.
        lo_default, hi_default = viewer._default_monitor_bins()
        self.ax.set_xlim(min(lo_default, lo), max(hi_default, hi))
        if drawn > 1:
            self.ax.legend(fontsize="small")
        lo_kev, hi_kev = self._range
        # The monitor's name in the title and not just the bins: two of these
        # windows can be open at once, and a screenshot of one has to say
        # which element it chose a band on.
        self.ax.set_title(
            f"{self.monitor}: bins {lo}-{hi}, {lo_kev:.3f}-{hi_kev:.3f} keV"
        )

        # Rebuilt after the clear above, through the parent's own map, so this
        # window and the spectrum pane cannot disagree about an energy.
        self.energy_axis = self.ax.secondary_xaxis(
            "top", functions=viewer._energy_map()
        )
        self.energy_axis.set_xlabel(
            "nominal energy [keV] - inferred, not recorded"
        )
        self._make_span()
        self.canvas.draw_idle()

    def closeEvent(self, event):
        """Let the window shut without taking the chosen range with it."""
        self.hide()
        event.ignore()


class SpectrumViewerWindow(QtWidgets.QMainWindow):
    """Plots one or more MCA spectra, overlaid, against bin number.

    Several scans share one frame range and one set of element ticks, because
    the reason to open them together is to compare them and a comparison needs
    the same reduction on both sides. A scan shorter than the range gets what
    it has, and the status bar says which were clamped.
    """

    def __init__(self, paths, parent=None):
        super().__init__(parent)
        paths = [
            Path(p)
            for p in ([paths] if isinstance(paths, (str, Path)) else paths)
        ]
        if not paths:
            raise RuntimeError("No spectrum file to open.")

        self.data = [SpectrumData(p) for p in paths]
        total = sum(d.n_frames for d in self.data)
        if total > MAX_TOTAL_FRAMES:
            raise RuntimeError(
                f"The selected files hold {total} frames between them, more "
                f"than this viewer summarises at once ({MAX_TOTAL_FRAMES})."
            )

        self.max_elements = max(d.n_elements for d in self.data)
        self.max_frames = max(d.n_frames for d in self.data)
        self.max_bins = max(d.n_bins for d in self.data)
        self.any_dead_time = any(d.have_dead_time for d in self.data)

        # Primary-stream columns, by (scan index, field). One open per field
        # per file, not one per redraw -- ticking an element must not go near
        # the disk. Filled first by the monitor survey just below.
        self._field_cache: dict[tuple[int, str], np.ndarray | None] = {}
        self._discover_trend_fields()

        names = ", ".join(d.path.name for d in self.data)
        self.setWindowTitle(f"{WINDOW_TITLE} - {names}")
        # Two plots side by side, so wider than the single pane wanted.
        self.resize(1400, 680)

        # One timer for every control: a redraw is a redraw, whichever box
        # started it, and restarting a single-shot timer is how the keystrokes
        # of a held-down arrow collapse into one.
        self._redraw_timer = QtCore.QTimer(self)
        self._redraw_timer.setSingleShot(True)
        self._redraw_timer.setInterval(REDRAW_DELAY_MS)
        self._redraw_timer.timeout.connect(self._redraw)

        # Measured once, from every frame and every element, and kept. The bin
        # range must not follow the selection -- a window that re-scaled itself
        # each time an element was unticked would make two curves impossible to
        # compare, which is the one thing the element toggles are for.
        combined = self._combined_total()
        self.default_bin_range = infer_bin_range(combined)

        # The ROI the trend pane integrates, held as the energy it names and
        # not as the bins it currently covers. Which way round that goes is a
        # decision and not an accident: the bins are re-derived whenever the
        # gain or offset below change, so correcting a nominal calibration
        # leaves the ROI on the line it was set to and walks the shaded band
        # away from a peak that is no longer there -- which is the symptom
        # worth seeing. Holding bins instead would keep the band on the peak
        # and silently change which energy it stood for.
        #
        # That is the choice *within* a window. What crosses to the next one
        # is bins -- see _REMEMBERED_ROI_BINS, which explains why the two
        # answers differ and why neither is a mistake.
        #
        # Failing a remembered one, it opens one half-width either side of
        # the tallest bin on screen. The maximum is over every element
        # summed, so on these files it is
        # the elastic peak at the incident energy rather than a fluorescence
        # line -- element 1 of scan 55 peaks at bin 635 while the sum peaks at
        # 1120. The "Peak" button is how it gets moved to a line.
        #
        # The centre itself is set in _build_roi_row, not here: either way it
        # ends as a keV, and converting a bin to one goes through the gain
        # and offset boxes, which do not exist until _build_ui runs.
        self.roi_center_kev = 0.0
        self.roi_width_ev = DEFAULT_ROI_WIDTH_EV
        self.roi_span = None
        self._syncing_roi = False

        # The window each computed monitor counts in, by monitor name, kept
        # here rather than in the dialogs that edit them: closing one of those
        # must not take its range with it, and _monitor_scale has to be able
        # to read it whether the window was ever opened or not. In keV, like
        # the ROI, so a change of gain moves them with everything else. Both
        # open on DEFAULT_MONITOR_BINS -- see _build_roi_row, which fills
        # these once the gain and offset boxes they convert through exist.
        #
        # One entry per monitor although both open on the same band: they are
        # bands on different elements, and the point of the dialog is to move
        # one of them. Switching monitors and switching back finds each where
        # it was left.
        self.monitor_ranges_kev: dict[str, tuple[float, float]] = {}
        # Whether the status bar has yet said that this monitor is still on
        # the band it opened with. One note per monitor per window, and only
        # once it matters -- see the end of _redraw_trend.
        self._monitor_default_note: dict[str, bool] = {}
        self.monitor_windows: dict[str, MonitorRoiWindow] = {}
        # What the spectrum pane last drew, summed: the "Peak" button searches
        # this, so that it finds the peak of what is on screen rather than of
        # a selection nobody is looking at. Filled by every _redraw.
        self._shown_total = combined

        self._build_ui()
        self._redraw()

        # The scans hint the computed monitor, so the combo usually opens on
        # it and the trend is already being divided by a window nobody chose.
        # Showing that window at once is the honest opening state: the full
        # spectrum is a placeholder, not a decision, and this is the only
        # thing that says so.
        if self._computed_monitor():
            self._open_monitor_window()

    def _discover_trend_fields(self):
        """What the masters offer the trend: an X axis, and monitors.

        The *intersection* across loaded scans, not the union, because one
        combo drives every curve: a field only one of them records would plot
        that one and silently drop the others. A scan whose master is missing
        contributes nothing, which empties the intersection and leaves the
        trend on its frame axis -- the honest answer, since without the master
        there is nothing to call a scan point.
        """
        per_scan = [set(primary_fields(d.master)) for d in self.data]
        common = set.intersection(*per_scan) if per_scan else set()
        self.axis_fields = sorted(common)
        self.monitor_fields = self._usable_monitors(common)
        # The first scan's hint speaks for the window: several scans opened
        # together are opened to be compared, which means they are the same
        # kind of scan.
        self.hinted_axis = scan_axis_field(
            self.data[0].master, self.axis_fields
        )

    def _usable_monitors(self, common: set) -> list[str]:
        """The monitors worth offering: recorded ones, then computed ones.

        Being recorded is not enough. Every 2026-3 scan here writes ``bI0`` as
        identically zero -- the channel is in the file, the beamline was not
        using it -- and since a zero divisor is read as "no beam" and blanked,
        picking it by name would open the trend on nothing at all. So each
        candidate is read once and kept only where every loaded scan has some
        positive reading for it, which leaves ``gI0``/``hI0`` on 2026-3 and
        ``4idgI0`` on 2025-1. The read is cached, so it is also the read the
        first draw would have done.

        The computed monitors are gated on a different thing, because the
        master records nothing to gate them on: what has to exist is the
        *element*, with counts in it, in every loaded scan. That is what keeps
        ``vortex_total_el7`` off the list for the 2025-1 files, whose detector
        has four elements. Same shape of rule as above -- every scan, and
        positive -- so a window of mixed detectors offers only what all of
        them can answer.

        Computed monitors come last, which is what keeps ``gI0``/``hI0`` the
        opening choice: the combo opens on the first entry.
        """
        usable = []
        for name in MONITOR_CANDIDATES:
            if name not in common:
                continue
            values = [
                self._field_values(i, name) for i in range(len(self.data))
            ]
            if all(v is not None and np.any(v > 0) for v in values):
                usable.append(name)
        for name, element in COMPUTED_MONITORS.items():
            if all(d.n_elements >= element for d in self.data) and np.any(
                self._combined_total(element) > 0
            ):
                usable.append(name)
        return usable

    def _combined_total(self, element: int | None = None) -> np.ndarray:
        """Every loaded scan summed over all frames, on one axis.

        Over all elements, or over the one 1-based *element* if given --
        which is how the monitor's opening window is chosen, since the
        monitor reads one element and the sum of four is a different
        spectrum with different lines in it.

        Dead-time corrected where the file allows it, so the opening range
        matches the curves first drawn. Cheap either way: the pyramid answers
        a full range without touching the cube.
        """
        total = np.zeros(self.max_bins, dtype=np.float64)
        for data in self.data:
            summed = data.range_sum(0, data.n_frames - 1, data.have_dead_time)
            if element is None:
                total[: data.n_bins] += summed.sum(axis=0)
            elif data.n_elements >= element:
                # A detector with fewer elements contributes nothing rather
                # than raising: _computed_scale already has to tolerate one,
                # and a window of mixed detectors should still open.
                total[: data.n_bins] += summed[element - 1]
        return total

    # -- construction -------------------------------------------------------

    def _build_ui(self):
        panel = QtWidgets.QWidget()
        self.setCentralWidget(panel)
        layout = QtWidgets.QVBoxLayout(panel)

        layout.addLayout(self._build_frame_row())
        layout.addLayout(self._build_element_row())
        layout.addLayout(self._build_range_row())
        # Built before the ROI row although it is shown below it: the ROI is
        # held in keV, and its opening centre is a bin converted through these
        # two boxes, so they have to exist first.
        axis_row = self._build_axis_row()
        layout.addLayout(self._build_roi_row())
        layout.addLayout(axis_row)

        self.figure = Figure(figsize=(7, 4.5), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot(111)
        # Tracks the bottom axis through zoom and pan on its own; set_functions
        # is what a change of gain goes through, rather than rebuilding it.
        self.energy_axis = self.ax.secondary_xaxis(
            "top", functions=self._energy_map()
        )
        self.energy_axis.set_xlabel(
            "nominal energy [keV] - inferred, not recorded"
        )
        self.toolbar = NavigationToolbar2QT(self.canvas, panel)

        # A second figure, not a second axes on the first: a splitter resizes
        # widgets, and self.figure already carries one axes plus the secondary
        # energy axis that _redraw rebuilds. The same reasoning, and the same
        # shape, as scan_viewer's trend pane.
        self.trend_figure = Figure(figsize=(5, 4.5), constrained_layout=True)
        self.trend_canvas = FigureCanvasQTAgg(self.trend_figure)
        self.trend_ax = self.trend_figure.add_subplot(111)
        self.trend_toolbar = NavigationToolbar2QT(self.trend_canvas, panel)

        # 2:1 stretch so that the share survives a window resize; the sizes
        # themselves are set once, the first time the trend pane is shown
        # (see _redraw_trend), so that a drag of the handle is not undone by
        # the next redraw.
        self.plot_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        panes = []
        for toolbar, canvas in (
            (self.toolbar, self.canvas),
            (self.trend_toolbar, self.trend_canvas),
        ):
            box = QtWidgets.QVBoxLayout()
            box.setContentsMargins(0, 0, 0, 0)
            box.addWidget(toolbar)
            box.addWidget(canvas, 1)
            pane = QtWidgets.QWidget()
            pane.setLayout(box)
            self.plot_split.addWidget(pane)
            panes.append(pane)
        self.trend_pane = panes[1]
        # Hidden until "Plot ROI" is pressed. Not merely cosmetic: the trend
        # costs a streamed pass over the cube, about a second on the 8000-frame
        # scan, and a window opened to look at a spectrum should not pay it.
        self.trend_pane.setVisible(False)
        self.plot_split.setStretchFactor(0, 2)
        self.plot_split.setStretchFactor(1, 1)
        layout.addWidget(self.plot_split, 1)

        self.status = self.statusBar()

    def _build_frame_row(self):
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Frames"))

        self.first_spin = QtWidgets.QSpinBox()
        self.last_spin = QtWidgets.QSpinBox()
        for spin in (self.first_spin, self.last_spin):
            spin.setRange(0, self.max_frames - 1)
            spin.setKeyboardTracking(False)
            spin.valueChanged.connect(self._schedule_redraw)
        self.last_spin.setValue(self.max_frames - 1)
        self.first_spin.valueChanged.connect(self._sync_frame_slider)
        row.addWidget(self.first_spin)
        row.addWidget(QtWidgets.QLabel("to"))
        row.addWidget(self.last_spin)

        row.addWidget(QtWidgets.QLabel(self._frames_note()))

        # Which frames of each scan point count, in the row that already says
        # how many there are per point -- the note above is the reason this
        # box is enabled, so they belong side by side.
        row.addWidget(QtWidgets.QLabel("Combine"))
        self.combine_combo = QtWidgets.QComboBox()
        for label, part, _ in COMBINE_MODES:
            self.combine_combo.addItem(label, part)
        usable = self._combine_available()
        self.combine_combo.setEnabled(usable)
        self.combine_combo.setToolTip(
            "What to do with the several frames of one scan point.\n"
            "Frames are numbered from 1 within their point, so frame 1 is odd:\n"
            '"Odd - even" is (frames 1, 3, 5 ...) minus (frames 2, 4, 6 ...).\n'
            "The softglue scans alternate frame by frame -- scan 597 point 0\n"
            "runs 5.4% apart, against 0.4% for a typical random split of\n"
            "the same frames -- and that difference is the measurement.\n"
            "Odd only and Even only show the two halves the difference is of;\n"
            "a difference is worth believing once both have been looked at."
            if usable
            else "Needs more than one vortex frame per scan point, in every scan\n"
            "loaded here. These have one, so there are no halves to tell apart."
        )
        self.combine_combo.currentIndexChanged.connect(self._schedule_redraw)
        row.addWidget(self.combine_combo)

        all_btn = QtWidgets.QPushButton("All frames")
        all_btn.clicked.connect(self._select_all_frames)
        row.addWidget(all_btn)

        # Drags the chosen window of frames through the scan without changing
        # its width, so what it does is set by the boxes: leave them a frame
        # apart and it steps frame by frame, set them twenty apart and it walks
        # a scan point at a time. One control for both readings of "slide
        # through frames", and nothing to switch between them. At full width
        # there is nowhere to slide to, which is why the boxes come first.
        self.frame_slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.frame_slider.setRange(0, max(0, self.max_frames - 1))
        self.frame_slider.setToolTip(
            "Slide the selected frames through the scan, keeping their count.\n"
            "Narrow the range above first -- at full width it cannot move."
        )
        self.frame_slider.valueChanged.connect(self._on_frame_slider)
        row.addWidget(self.frame_slider, 1)
        return row

    def _frames_note(self) -> str:
        """ "800 frames, 20 per scan point", where the two divide evenly."""
        if len(self.data) != 1:
            return f"of {self.max_frames}"
        per_point = self.data[0].frames_per_point
        if per_point is None:
            return f"of {self.max_frames}"
        return f"of {self.max_frames}, {per_point} per scan point"

    def _combine_available(self) -> bool:
        """Whether every loaded scan can be split into odd and even frames.

        Every one, not any: a single Combine box drives all the curves, so a
        window holding one sgz scan and one ordinary scan would be asking for
        a difference the second cannot answer. Disabled is the honest state
        there, and the tooltip says which condition failed.
        """
        return bool(self.data) and all(
            d.multi_frame and d.odd_frames is not None for d in self.data
        )

    def _combine_part(self) -> str:
        """The frames the Combine box is asking for: see :data:`COMBINE_MODES`.

        ``"all"`` whenever the box is disabled, so that every caller can ask
        for a part without first asking whether parts exist here.
        """
        if not self.combine_combo.isEnabled():
            return "all"
        return self.combine_combo.currentData() or "all"

    @staticmethod
    def _count_label(part: str) -> str:
        """What a Y axis is counting under *part*, from :data:`COMBINE_MODES`."""
        for _, name, label in COMBINE_MODES:
            if name == part:
                return label
        return "counts"

    @staticmethod
    def _mode_label(part: str) -> str:
        """The Combine box's own wording for *part*, for titles and messages."""
        for label, name, _ in COMBINE_MODES:
            if name == part:
                return label
        return part

    def _build_element_row(self):
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Elements"))

        self.element_boxes = []
        for i in range(self.max_elements):
            box = QtWidgets.QCheckBox(str(i + 1))
            box.setChecked(True)
            box.toggled.connect(self._schedule_redraw)
            self.element_boxes.append(box)
            row.addWidget(box)

        for label, slot in (("All", True), ("None", False)):
            btn = QtWidgets.QPushButton(label)
            btn.clicked.connect(
                lambda _=False, on=slot: self._set_all_elements(on)
            )
            row.addWidget(btn)

        self.separate_box = QtWidgets.QCheckBox("Separate")
        self.separate_box.setToolTip(
            "Draw one curve per element instead of their sum.\n"
            "On by default, where there is a single scan to do it for:\n"
            "the elements are what a single scan is worth looking at, and\n"
            "their sum hides an element that is dead or badly out of step.\n"
            "Several scans times seven elements is a legend nobody can read,\n"
            "so with more than one loaded this is off and cannot be set."
        )
        # Checked and enabled together: checking a box that is ignored would
        # say the curves are separate when _redraw has already overruled it.
        self.separate_box.setEnabled(len(self.data) == 1)
        self.separate_box.setChecked(len(self.data) == 1)
        self.separate_box.toggled.connect(self._schedule_redraw)
        row.addWidget(self.separate_box)

        self.deadtime_box = QtWidgets.QCheckBox("Dead-time corrected")
        self.deadtime_box.setChecked(self.any_dead_time)
        self.deadtime_box.setEnabled(self.any_dead_time)
        self.deadtime_box.setToolTip(
            "Scale each frame by its DTFactor (or ICR/OCR).\n"
            "On by default: the elements run at very different dead times\n"
            "-- 8-33% against 0.2% on scan 55 -- so uncorrected curves make\n"
            "a busy element look like a weak one."
            if self.any_dead_time
            else "This file records no dead-time factors."
        )
        self.deadtime_box.toggled.connect(self._schedule_redraw)
        row.addWidget(self.deadtime_box)

        self.log_box = QtWidgets.QCheckBox("Log Y")
        self.log_box.setChecked(True)
        self.log_box.setToolTip(
            "On by default: the elastic peak is some 1000x the fluorescence\n"
            "lines, which on a linear scale are a flat line on the axis."
        )
        self.log_box.toggled.connect(self._schedule_redraw)
        row.addWidget(self.log_box)

        row.addStretch(1)
        return row

    def _build_range_row(self):
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Bins shown"))

        low, high = self.default_bin_range
        self.xmin_spin = QtWidgets.QSpinBox()
        self.xmax_spin = QtWidgets.QSpinBox()
        for spin, value in ((self.xmin_spin, low), (self.xmax_spin, high)):
            spin.setRange(0, self.max_bins - 1)
            spin.setKeyboardTracking(False)
            # Set before connecting: these are the opening values, not an edit.
            spin.setValue(value)
            spin.valueChanged.connect(self._schedule_redraw)
        self.xmin_spin.setToolTip(
            "The bin range drawn. It is yours once set -- nothing but these\n"
            "boxes, the buttons beside them and the toolbar's zoom changes it,\n"
            "so elements and frames can be switched without the axis moving."
        )
        self.xmax_spin.setToolTip(self.xmin_spin.toolTip())
        row.addWidget(self.xmin_spin)
        row.addWidget(QtWidgets.QLabel("to"))
        row.addWidget(self.xmax_spin)

        auto_btn = QtWidgets.QPushButton("Auto")
        auto_btn.setToolTip(
            f"Back to the opening range, {low}-{high}: the bins holding "
            f"{BIN_MASS_QUANTILE:.0%} of\nthe counts of every frame and element, "
            "plus a margin."
        )
        auto_btn.clicked.connect(
            lambda: self._set_bin_range(*self.default_bin_range)
        )
        row.addWidget(auto_btn)

        full_btn = QtWidgets.QPushButton("Full")
        full_btn.setToolTip(
            f"All {self.max_bins} bins, including the last, which is not a\n"
            "spectrum bin but where everything over the ADC range piles up."
        )
        full_btn.clicked.connect(
            lambda: self._set_bin_range(0, self.max_bins - 1)
        )
        row.addWidget(full_btn)

        row.addStretch(1)
        return row

    def _build_roi_row(self):
        """The ROI, in keV, and what the pane beside the spectrum plots it against."""
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("ROI"))

        # The opening ROI: wherever the last window of this session was left,
        # else the tallest bin of everything loaded, inside the range the
        # spectrum opens on. Done here rather than in __init__ because either
        # way it ends as a keV, and converting a bin to one goes through the
        # gain and offset boxes, which _build_ui creates just before calling
        # this.
        #
        # The remembered bins are taken only if they fit this spectrum. A
        # detector with fewer channels than the one before would otherwise
        # have its ROI clamped to the last bin and open on a sliver at the top
        # of the range -- a worse answer than finding the peak again.
        remembered = _REMEMBERED_ROI_BINS
        if (
            remembered is not None
            and 0 <= remembered[0] < remembered[1] <= self.max_bins - 1
        ):
            lo, hi = remembered
            # The inverse of _roi_bins, as _on_span_select already writes it.
            # Round-trips exactly at the default 10 eV/bin and zero offset:
            # the centre lands on a multiple of 5 eV and the half-width on a
            # whole eV, both of which the boxes' 3 and 0 decimals hold.
            self.roi_center_kev = float(self._bin_to_kev((lo + hi) / 2.0))
            self.roi_width_ev = (
                (float(self._bin_to_kev(hi)) - float(self._bin_to_kev(lo)))
                * 1000.0
                / 2.0
            )
        else:
            self.roi_center_kev = self._peak_kev(
                self._shown_total, *self.default_bin_range
            )

        # Both computed monitors open on DEFAULT_MONITOR_BINS, the same band
        # on two different elements. Here for the same reason and in the same
        # place as the ROI above: they are keV ranges too.
        #
        # A fixed band rather than a peak found per element, which is what
        # this did before. The reasoning that drove the peak-finding still
        # holds, and having a wide default does not soften it: every region
        # of this detector alternates with the frames (measured on scan 597:
        # element 1 -2.4%, element 7 -7.4% over this band), so there is no
        # quiet corner on either element and no window here is a reference.
        # What the fixed band buys is that neither monitor opens *on the line
        # being measured*. The per-element lowest peak did exactly that on
        # element 1 on every 2026-3 file, dividing the measurement by itself
        # for a 0.000% asymmetry; 200-1800 holds 96% of either element, so no
        # single line can dominate it. The default is still a starting point
        # and not an answer: it is named in the y-axis label and the export
        # header (_monitor_label), and _redraw_trend says in the status bar,
        # once per monitor, that nobody chose it.
        #
        # Every computed monitor gets a window, not just the offered ones:
        # _monitor_scale can be asked for any of them by a caller that is not
        # the combo, and a missing key would be a KeyError rather than a
        # fallback.
        lo_bin, hi_bin = self._default_monitor_bins()
        for name in COMPUTED_MONITORS:
            self.monitor_ranges_kev[name] = (
                float(self._bin_to_kev(lo_bin)),
                float(self._bin_to_kev(hi_bin)),
            )
            self._monitor_default_note[name] = True

        self.roi_center_spin = QtWidgets.QDoubleSpinBox()
        self.roi_center_spin.setDecimals(3)
        self.roi_center_spin.setSingleStep(0.01)
        self.roi_center_spin.setSuffix(" keV")

        self.roi_width_spin = QtWidgets.QDoubleSpinBox()
        self.roi_width_spin.setDecimals(0)
        self.roi_width_spin.setSingleStep(10.0)
        self.roi_width_spin.setPrefix("± ")
        self.roi_width_spin.setSuffix(" eV")

        # Sized to the widest reading rather than to the range, whose ±1000
        # keV and 100000 eV are guards against a recalibration (see
        # _set_roi_spin_limits) and not values anybody types: fitting those
        # would give two boxes three times the width they need.
        fit_spin_width(self.roi_center_spin, "-40.950 keV")
        fit_spin_width(self.roi_width_spin, "± 20475 eV")

        self._set_roi_spin_limits()
        for spin, value in (
            (self.roi_center_spin, self.roi_center_kev),
            (self.roi_width_spin, self.roi_width_ev),
        ):
            spin.setKeyboardTracking(False)
            # Set before connecting: these are the opening values, not an edit.
            spin.setValue(value)
            spin.valueChanged.connect(self._on_roi_spin)
        # What the boxes actually took, rounded to their decimals, so the
        # model and the display start out agreeing -- as _set_roi_kev keeps
        # them from here on.
        self.roi_center_kev = self.roi_center_spin.value()
        self.roi_width_ev = self.roi_width_spin.value()
        self.roi_center_spin.setToolTip(
            "The region the pane on the right integrates, shaded on the\n"
            "spectrum. Drag the band or its edges to move it; the boxes and\n"
            "the band are the same region written two ways.\n"
            "Held as an energy, not as bins: correct the gain below and the\n"
            "ROI stays on the energy it names, so a band that walks off a\n"
            "peak is telling you the nominal calibration was wrong.\n"
            "Opens on the tallest bin of everything loaded -- which on these\n"
            "files is the elastic peak, not a fluorescence line."
        )
        self.roi_width_spin.setToolTip(self.roi_center_spin.toolTip())
        row.addWidget(self.roi_center_spin)
        row.addWidget(self.roi_width_spin)

        peak_btn = QtWidgets.QPushButton("Peak")
        peak_btn.setToolTip(
            "Centre the ROI on the tallest bin of the curves now drawn,\n"
            "inside the bins now shown. Zoom onto a line and press this to\n"
            "put the ROI on it. The width is left alone."
        )
        peak_btn.clicked.connect(self._center_on_peak)
        row.addWidget(peak_btn)

        self.plot_roi_btn = QtWidgets.QPushButton("Plot ROI")
        self.plot_roi_btn.setCheckable(True)
        self.plot_roi_btn.setToolTip(
            "Show the pane that plots the ROI across the scan.\n"
            "Off at first because it costs a pass over the cube -- about a\n"
            "second on the 8000-frame scan -- which a window opened to look\n"
            "at a spectrum should not pay for unasked."
        )
        self.plot_roi_btn.toggled.connect(self._schedule_redraw)
        row.addWidget(self.plot_roi_btn)

        self.export_btn = QtWidgets.QPushButton("Export...")
        self.export_btn.setToolTip(
            "Write the plotted ROI trend to a text file: one row per scan\n"
            "point, and all four Combine modes in it -- Sum, Odd - even,\n"
            "Odd only, Even only -- whichever one is on screen. Normalised\n"
            "and raw side by side, with the monitor counts that divided\n"
            "them, so any of it can be recomputed or undone.\n"
            "Goes to an 'export' folder beside the scan data, named for the\n"
            "scan: data_000597.txt."
        )
        # Enabled with the trend and not before: there is nothing to write
        # until the reduce has been asked for, and the button saying so is
        # clearer than a message afterwards. Same rule as 'Separate'.
        self.export_btn.setEnabled(self.plot_roi_btn.isChecked())
        self.plot_roi_btn.toggled.connect(self.export_btn.setEnabled)
        self.export_btn.clicked.connect(self._export_trend)
        row.addWidget(self.export_btn)

        row.addSpacing(16)
        row.addWidget(QtWidgets.QLabel("against"))
        self.x_combo = QtWidgets.QComboBox()
        # Point and frame first: they need no master and always work. The
        # recorded fields follow, and one of them is normally what opens --
        # see the default chosen below.
        self.x_combo.addItem(TREND_POINT_AXIS, ("point", ""))
        self.x_combo.addItem(TREND_FRAME_AXIS, ("frame", ""))
        # The index is kept as the items go in rather than looked up after
        # with findData, which does not match a tuple stored as user data.
        default_index = 0
        for name in self.axis_fields:
            if name == self.hinted_axis:
                default_index = self.x_combo.count()
            self.x_combo.addItem(name, ("field", name))
        self.x_combo.setCurrentIndex(default_index)
        self.x_combo.setToolTip(
            "What the ROI value is plotted against. Opens on the field the\n"
            "scan's own hints name -- energy for an energy scan, the motor\n"
            "for a lup -- so nothing here is hardcoded to one kind of scan.\n"
            "'frame' is every frame as written; the rest are one point each,\n"
            "averaging the frames the detector took there."
        )
        # Allowed to shrink, and to take the row's spare width when there is
        # any. Left alone it asks for its longest field name as a *minimum* --
        # 302 px on a scan with long names -- which is three buttons' worth of
        # window that no other row needs, and it was setting the whole
        # window's minimum width. It still opens at full width; it can now be
        # narrowed, and the names it elides are in the list it drops down.
        self.x_combo.setSizeAdjustPolicy(
            QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.x_combo.setMinimumContentsLength(12)
        self.x_combo.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed,
        )
        self.x_combo.currentIndexChanged.connect(self._schedule_redraw)
        row.addWidget(self.x_combo, 1)

        self.normalise_box = QtWidgets.QCheckBox("/ monitor")
        self.monitor_combo = QtWidgets.QComboBox()
        self.monitor_combo.addItems(self.monitor_fields)
        have = bool(self.monitor_fields)
        self.normalise_box.setChecked(have)
        self.normalise_box.setEnabled(have)
        self.monitor_combo.setEnabled(have)
        self.normalise_box.setToolTip(
            "Divide the ROI value by a beam monitor, so the trend shows the\n"
            "sample and not the ring.\n"
            "On by default, which is a departure from the other viewers --\n"
            "they open on raw recorded numbers. One click turns it off.\n"
            "The vortex_total_el* monitors are not read from the master:\n"
            "each is that element of this cube, summed over a window you\n"
            "pick in its Monitor ROI window, and totalled over the same\n"
            "frames each trend point totals -- that parity's frames under\n"
            "Odd only, Even only and Odd - even. That per-frame divisor is\n"
            "the whole point of them: a recorded monitor has one reading\n"
            "per scan point and cannot tell the odd frames from the even.\n"
            "element 1 carries the signal, so on these scans its window\n"
            "opens on the very line being measured; element 7 never does,\n"
            "but alternates harder and flips the sign of the asymmetry.\n"
            "Neither is a quiet reference -- read the Monitor ROI tooltip.\n"
            "A monitor the file records as all zeros is not offered, which\n"
            "bI0 is on these scans, and an element the detector does not\n"
            "have is not offered either, which element 7 is on 2025-1."
            if have
            else "No live beam monitor among the fields these files record."
        )
        self.monitor_combo.setToolTip(self.normalise_box.toolTip())

        self.monitor_roi_btn = QtWidgets.QPushButton("Monitor ROI...")
        self.monitor_roi_btn.setToolTip(
            "Pick the energy window the selected monitor is summed over.\n"
            "Enabled only while one of the computed monitors is selected,\n"
            "since those are the ones this viewer computes rather than\n"
            "reads. Each keeps its own window: switching between them and\n"
            "back finds each band where it was left."
        )
        self.monitor_roi_btn.clicked.connect(self._open_monitor_window)

        self.normalise_box.toggled.connect(self._schedule_redraw)
        # Before the redraw connection, so that picking the computed monitor
        # has its window open by the time the trend is drawn -- the window is
        # what the user then drags, and a redraw behind it would be the only
        # feedback that anything happened.
        self.monitor_combo.currentIndexChanged.connect(self._on_monitor_changed)
        self.monitor_combo.currentIndexChanged.connect(self._schedule_redraw)
        row.addWidget(self.normalise_box)
        row.addWidget(self.monitor_combo)
        row.addWidget(self.monitor_roi_btn)
        # Set from the opening selection rather than assumed off: the combo
        # may well open on the computed monitor, since the scans hint it.
        self.monitor_roi_btn.setEnabled(
            have and self.monitor_combo.currentText() in COMPUTED_MONITORS
        )

        row.addStretch(1)
        return row

    def _build_axis_row(self):
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Energy axis:"))

        self.gain_spin = QtWidgets.QDoubleSpinBox()
        self.gain_spin.setRange(MIN_EV_PER_BIN, 1000.0)
        self.gain_spin.setDecimals(3)
        self.gain_spin.setValue(DEFAULT_EV_PER_BIN)
        self.gain_spin.setSuffix(" eV/bin")
        self.gain_spin.setKeyboardTracking(False)
        self.gain_spin.valueChanged.connect(self._schedule_redraw)
        row.addWidget(self.gain_spin)

        self.offset_spin = QtWidgets.QDoubleSpinBox()
        self.offset_spin.setRange(-100000.0, 100000.0)
        self.offset_spin.setDecimals(1)
        self.offset_spin.setValue(DEFAULT_EV_OFFSET)
        self.offset_spin.setPrefix("offset ")
        self.offset_spin.setSuffix(" eV")
        self.offset_spin.setKeyboardTracking(False)
        self.offset_spin.valueChanged.connect(self._schedule_redraw)
        row.addWidget(self.offset_spin)

        note = QtWidgets.QLabel(
            "nominal - no calibration is recorded in the file"
        )
        note.setToolTip(
            "The defaults were read off known lines on scan 55 (Fe Ka at bin\n"
            "635, elastic at 1117). Type over them when the gain is known."
        )
        row.addWidget(note)
        row.addStretch(1)
        return row

    # -- controls -----------------------------------------------------------

    def _select_all_frames(self):
        self._set_frames(0, self.max_frames - 1)

    def _set_frames(self, first: int, last: int):
        """Write both frame boxes and the slider at once, then redraw once.

        Signals are blocked rather than left to fire because each box would
        otherwise schedule its own redraw and drag the slider halfway through
        the pair, with the second box not yet written.
        """
        for spin, value in ((self.first_spin, first), (self.last_spin, last)):
            spin.blockSignals(True)
            spin.setValue(value)
            spin.blockSignals(False)
        self._sync_frame_slider()
        self._schedule_redraw()

    def _sync_frame_slider(self, *_):
        """Park the slider on the first selected frame, without a redraw."""
        self.frame_slider.blockSignals(True)
        self.frame_slider.setValue(self.first_spin.value())
        self.frame_slider.blockSignals(False)

    def _on_frame_slider(self, value: int):
        """Move the selected frames to start at *value*, keeping their count.

        Clamped at the far end rather than allowed to shrink: a window that
        narrowed as it reached the last frame would quietly change what is
        being summed just as the interesting part of the scan arrived.
        """
        width = self.last_spin.value() - self.first_spin.value()
        first = max(0, min(value, self.max_frames - 1 - width))
        self._set_frames(first, min(first + width, self.max_frames - 1))

    def _set_bin_range(self, low: int, high: int):
        for spin, value in ((self.xmin_spin, low), (self.xmax_spin, high)):
            spin.blockSignals(True)
            spin.setValue(value)
            spin.blockSignals(False)
        self._schedule_redraw()

    def _on_xlim_changed(self, ax):
        """Follow a toolbar zoom or pan back into the bin boxes.

        Without this the boxes and the axis would disagree the moment anyone
        used the magnifier, and the next redraw -- any element toggled -- would
        throw the zoom away and snap back to whatever the boxes still said.
        Writing them here instead makes the zoom the new range, and the boxes
        go on being the one place the range is kept. Signals are blocked: this
        is a redraw reporting where it landed, not a request for another.
        """
        low, high = sorted(ax.get_xlim())
        low = max(0, min(self.max_bins - 1, int(round(low))))
        high = max(0, min(self.max_bins - 1, int(round(high))))
        if high <= low:
            return
        for spin, value in ((self.xmin_spin, low), (self.xmax_spin, high)):
            spin.blockSignals(True)
            spin.setValue(value)
            spin.blockSignals(False)

    # -- the ROI, which has two views ---------------------------------------

    def _set_roi_spin_limits(self):
        """Ranges for the two ROI boxes -- deliberately not the spectrum's span.

        Tying the centre to the bins would mean that dropping the gain to
        1 eV/bin, which puts the last bin at 4 keV, silently rewrote an 11.2
        keV ROI as 4.095. The ROI names an energy whether or not this detector
        can currently reach it; it is :meth:`_roi_bins` that clamps, and it
        clamps the bins.
        """
        self.roi_center_spin.setRange(-1000.0, 1000.0)
        # One eV, not zero: a zero half-width is a region of no width, and
        # while _roi_bins would widen it to a single bin, a box reading "± 0"
        # over a one-bin window states something that is not true.
        self.roi_width_spin.setRange(1.0, 100000.0)

    def _default_monitor_bins(self) -> tuple[int, int]:
        """:data:`DEFAULT_MONITOR_BINS`, clamped to this file's bin count.

        Both computed monitors open on it, their dialog plots it, and
        :meth:`_monitor_window_is_default` compares against it -- one place,
        so a detector with fewer bins than the constant assumes cannot be
        clamped in two of the three and not the third.
        """
        last = max(1, self.max_bins - 1)
        lo, hi = DEFAULT_MONITOR_BINS
        lo = max(0, min(int(lo), last - 1))
        hi = max(lo + 1, min(int(hi), last))
        return lo, hi

    def monitor_bins(self, monitor: str) -> tuple[int, int]:
        """The inclusive bin window one computed monitor is summed over.

        Derived from :attr:`monitor_ranges_kev` on every use, exactly as
        :meth:`_roi_bins` is derived from the ROI boxes: the range is an
        energy, so a recalibration moves the monitor's bins without anyone
        reopening the dialog.
        """
        lo_kev, hi_kev = self.monitor_ranges_kev[monitor]
        lo = int(round(float(self._kev_to_bin(lo_kev))))
        hi = int(round(float(self._kev_to_bin(hi_kev))))
        lo = max(0, min(lo, self.max_bins - 1))
        hi = max(0, min(hi, self.max_bins - 1))
        if hi <= lo:
            # The range fell off an end, which a large enough gain does. Keep
            # a bin rather than divide the trend by an empty sum.
            hi = min(lo + 1, self.max_bins - 1)
            lo = max(0, hi - 1)
        return lo, hi

    def _computed_monitor(self) -> bool:
        """Whether the monitor in force is one this viewer computes.

        Normalisation has to be on as well: with the box unticked there is no
        divisor at all, and the window that chooses its range has nothing to
        change.
        """
        return (
            self.normalise_box.isChecked()
            and self.monitor_combo.isEnabled()
            and self.monitor_combo.currentText() in COMPUTED_MONITORS
        )

    def _on_monitor_changed(self, *_):
        """Show the selected monitor's ROI window; hide the other one.

        Hidden and not destroyed: the ranges live on the viewer, so nothing is
        lost either way, but keeping the instances means each monitor's chosen
        window is still there after a detour through the other one.
        """
        name = self.monitor_combo.currentText()
        wanted = name in COMPUTED_MONITORS
        self.monitor_roi_btn.setEnabled(
            wanted and self.monitor_combo.isEnabled()
        )
        # Hide every other computed monitor's window first, so that switching
        # from one to the other cannot leave two up: only one of them is
        # dividing the trend, and two open bands would both look live.
        for other, window in self.monitor_windows.items():
            if other != name:
                window.hide()
        if wanted:
            self._open_monitor_window()

    def _open_monitor_window(self, *_):
        """Show the selected monitor's ROI window, building it the first time.

        One instance per computed monitor, raised rather than duplicated if it
        is already up -- pressing the button twice, or reselecting the
        monitor, must not leave two windows editing one range.
        """
        name = self.monitor_combo.currentText()
        if name not in COMPUTED_MONITORS:
            return
        window = self.monitor_windows.get(name)
        if window is None:
            window = MonitorRoiWindow(self, name)
            self.monitor_windows[name] = window
        window.refresh()
        window.show()
        window.raise_()
        window.activateWindow()

    def _roi_bins(self) -> tuple[int, int]:
        """The inclusive bin window the centre and half-width currently name.

        Derived on every use and never stored, which is what makes the ROI an
        energy: change the gain and this returns different bins for the same
        two boxes.

        Rounded to nearest at each end rather than floor/ceil as
        :mod:`image_viewer` does for a pixel box, because this window is
        specified in energy and centred on a value: half a bin of slop at each
        end is what "nearest" means, and it keeps the window symmetric. The
        window is inclusive of both ends, so it holds one bin more than the
        width alone suggests -- 41 bins for ± 200 eV at 10 eV/bin.
        """
        half = self.roi_width_ev / 1000.0
        low = int(round(float(self._kev_to_bin(self.roi_center_kev - half))))
        high = int(round(float(self._kev_to_bin(self.roi_center_kev + half))))
        low = max(0, min(low, self.max_bins - 1))
        high = max(0, min(high, self.max_bins - 1))
        if high <= low:
            # The whole window fell off one end of the spectrum, which a gain
            # small enough does. Keep a bin rather than hand an empty slice to
            # the reader.
            if high < self.max_bins - 1:
                high = low + 1
            else:
                low = high - 1
        return low, high

    def _peak_kev(self, total: np.ndarray, low: int, high: int) -> float:
        """The energy of the tallest bin in ``[low, high]``.

        Not smoothed, and not interpolated between bins. Neither sibling
        viewer smooths before a peak search either, and there is nothing to
        smooth away here: this runs on a spectrum summed over every frame and
        element, whose peak carries 213180 counts on scan 55. A 10 eV bin is
        also far finer than the 310-500 eV peaks it is finding.
        """
        body = np.asarray(total, dtype=float)
        if body.size == 0:
            return 0.0
        low = max(0, min(int(low), body.size - 1))
        high = max(low, min(int(high), body.size - 1))
        peak = low + int(np.argmax(body[low : high + 1]))
        return float(self._bin_to_kev(peak))

    def _center_on_peak(self):
        """The "Peak" button: re-centre on the tallest bin now on screen."""
        self._set_roi_kev(
            self._peak_kev(
                self._shown_total,
                self.xmin_spin.value(),
                self.xmax_spin.value(),
            ),
            self.roi_width_ev,
        )

    def _set_roi_kev(
        self, center: float, width: float, update_span: bool = True
    ):
        """Move the ROI, and bring both views with it.

        The one place ``roi_center_kev``/``roi_width_ev`` change, and so the
        one place that records a *chosen* ROI for the next window.

        The re-entrancy guard is not about the span, which does not re-fire
        ``onselect`` when its extents are set programmatically: it is about
        the two spin boxes, which are written one after the other, so without
        it the first write would come straight back here pairing the new
        centre with the old width -- a region nobody asked for, recomputed
        from the file on its way past.
        """
        center = float(center)
        width = float(width)
        if self._syncing_roi:
            return
        if (center, width) == (self.roi_center_kev, self.roi_width_ev):
            # Nothing to move, but still a choice: pressing Peak when the band
            # is already on the peak, or typing back the number in the box,
            # says "this one" as plainly as a drag does, and the next window
            # should open there. So this returns early from the work and not
            # from the remembering -- the only writes that are not choices are
            # the re-entrant ones caught just above.
            self._remember_roi()
            return

        self._syncing_roi = True
        try:
            self.roi_center_spin.setValue(center)
            self.roi_width_spin.setValue(width)
            # Read back rather than keep what came in: the boxes round to
            # their decimals, and a drag arrives with far more precision than
            # they show. Taking theirs is what stops the title and the derived
            # bins from being computed off a number nowhere on screen.
            self.roi_center_kev = self.roi_center_spin.value()
            self.roi_width_ev = self.roi_width_spin.value()
            if update_span and self.roi_span is not None:
                self.roi_span.extents = self._roi_bins()
        finally:
            self._syncing_roi = False
        # After the read-back above, so what is carried is the ROI the boxes
        # settled on and not the sub-bin number a drag arrived with: the next
        # window must open where this one looks, not half a bin away.
        self._remember_roi()
        self._schedule_redraw()

    def _remember_roi(self):
        """Carry the current ROI to the next window of this session.

        Bins, because that is what survives a change of gain between two
        scans; see :data:`_REMEMBERED_ROI_BINS` for why that is the opposite
        of how a window holds its own ROI, and why both are right.
        """
        global _REMEMBERED_ROI_BINS
        _REMEMBERED_ROI_BINS = self._roi_bins()

    def _on_roi_spin(self, *_):
        self._set_roi_kev(
            self.roi_center_spin.value(), self.roi_width_spin.value()
        )

    def _on_span_select(self, low: float, high: float):
        """A drag of the shaded band on the spectrum, which is drawn in bins."""
        lo, hi = int(round(low)), int(round(high))
        if hi <= lo:
            # A click, not a drag: SpanSelector reports a zero-width span and
            # hides the band on its way out. Redraw to put the band back
            # rather than collapse the ROI onto a single bin.
            self._schedule_redraw()
            return
        # Back out through the same map the band was placed with, so a drag
        # and a typed value mean the same thing. Not back into the span: it is
        # already showing the drag, to sub-bin precision, and the redraw this
        # schedules rebuilds it on the rounded bins anyway.
        center = float(self._bin_to_kev((lo + hi) / 2.0))
        span_kev = float(self._bin_to_kev(hi)) - float(self._bin_to_kev(lo))
        self._set_roi_kev(center, span_kev * 1000.0 / 2.0, update_span=False)

    def _make_roi_span(self):
        """Rebuild the draggable band. Called at the end of every redraw.

        Both halves are necessary. The artist does not survive ``ax.clear()``,
        so there is nothing to keep; and the old selector's canvas callbacks
        do survive it, so without the disconnect the handler count climbed
        from 16 to 56 over six redraws.

        Nothing is needed to keep the toolbar and the band apart:
        ``SpanSelector.ignore`` consults ``canvas.widgetlock``, which zoom and
        pan hold while they are armed.
        """
        if self.roi_span is not None:
            self.roi_span.disconnect_events()
        self.roi_span = SpanSelector(
            self.ax,
            self._on_span_select,
            "horizontal",
            interactive=True,
            drag_from_anywhere=True,
            # Blitting caches the axes background, which a redraw underneath
            # invalidates; the band is cheap enough to draw with the rest.
            useblit=False,
            # Grey, not a colour: every colour here is one of the seven
            # elements' curves, and a band in tab:orange reads as element 2
            # having gone very wide.
            props={"facecolor": "0.35", "alpha": 0.20},
            handle_props={"color": "0.25"},
        )
        self.roi_span.extents = self._roi_bins()

    def _set_all_elements(self, on: bool):
        for box in self.element_boxes:
            box.blockSignals(True)
            box.setChecked(on)
            box.blockSignals(False)
        self._schedule_redraw()

    def _schedule_redraw(self, *_):
        self._redraw_timer.start()

    def _bin_to_kev(self, bins):
        """Bin number -> keV, through the gain and offset boxes as they are now.

        Reads the boxes on every call rather than closing over their values,
        so an edit of either takes effect without anything being rebuilt.
        Accepts a scalar or an array: matplotlib hands the secondary axis a
        whole tick vector, while the ROI asks about one bin at a time.
        """
        gain, offset = self.gain_spin.value(), self.offset_spin.value()
        return (offset + np.asarray(bins, dtype=float) * gain) / 1000.0

    def _kev_to_bin(self, kev):
        """keV -> bin number, fractional. The inverse of :meth:`_bin_to_kev`.

        The gain cannot be zero -- the spin box starts at
        :data:`MIN_EV_PER_BIN` -- so this is always defined.
        """
        gain, offset = self.gain_spin.value(), self.offset_spin.value()
        return (np.asarray(kev, dtype=float) * 1000.0 - offset) / gain

    def _energy_map(self):
        """The bin <-> keV pair the secondary axis is built from.

        The same two methods the ROI converts through, handed to matplotlib as
        a pair. One map and not two: the shaded band, the top axis and the
        bins the trend integrates over are then incapable of disagreeing.
        """
        return self._bin_to_kev, self._kev_to_bin

    # -- drawing ------------------------------------------------------------

    def _selected_elements(self) -> list[int]:
        return [
            i for i, box in enumerate(self.element_boxes) if box.isChecked()
        ]

    def _redraw(self):
        first = self.first_spin.value()
        last = self.last_spin.value()
        if last < first:
            self.status.showMessage("The last frame is before the first.")
            return

        xmin = self.xmin_spin.value()
        xmax = self.xmax_spin.value()
        if xmax <= xmin:
            self.status.showMessage(
                "The last bin shown is not after the first."
            )
            return

        # An empty selection is drawn, not short-circuited. Returning early
        # after a bare clear() left the axes without the bin range, the log
        # scale or the energy axis, so unticking the last element collapsed X
        # to 0-1 and flipped Y to linear under a Log Y box that was still
        # ticked -- the plot appeared to have changed in ways nothing asked
        # for. Falling through draws no curves and keeps the frame, and the
        # count of what was actually drawn settles the status line below. That
        # also covers a case the early return got wrong: elements ticked that
        # none of these files have (element 7 against the 4-element vintage)
        # drew nothing while the status still reported a frame range.
        chosen = self._selected_elements()

        dead_time = self.deadtime_box.isChecked()
        separate = self.separate_box.isChecked() and len(self.data) == 1
        part = self._combine_part()
        signed = part in SIGNED_PARTS
        bins = np.arange(self.max_bins)

        self.ax.clear()
        clamped: list[str] = []
        uncorrected: list[str] = []
        skipped: list[str] = []
        # Scans whose chosen half holds no frame of this range -- "Even only"
        # over frame 0 alone, say. Drawn as the row of zeros it honestly is,
        # but named, because a flat line at zero otherwise reads as a
        # detector that saw nothing rather than as a window with nothing in
        # it. Only the parity modes can produce it.
        empty_half: list[str] = []
        # Only the part of each curve inside the bin range, for the Y scaling
        # below. Collected here because this is the one place the curves exist
        # as arrays; reading them back off the Line2D objects afterwards would
        # work but says less about where they came from.
        shown: list[np.ndarray] = []
        # The same curves summed, at full length, which is what the "Peak"
        # button searches. Full length and not the shown slice because the
        # button hands it bin numbers, and those index this directly.
        #
        # Always the *sum*, whatever the Combine box says: Peak is for putting
        # the ROI on a line, and the tallest point of a difference is wherever
        # the two halves disagree most, which is not where the line is. The
        # pyramid answers a full sum in 0.009 s, so carrying both costs
        # nothing worth saving.
        total = np.zeros(self.max_bins)
        drawn = 0

        def draw(y, label, y_sum):
            nonlocal drawn
            self.ax.plot(bins[: y.size], y, lw=0.8, label=label)
            shown.append(y[xmin : xmax + 1])
            total[: y_sum.size] += y_sum
            drawn += 1

        for data in self.data:
            usable = [i for i in chosen if i < data.n_elements]
            if not usable:
                # It has none of the ticked elements -- a 4-element file in a
                # window whose boxes go to 7. Named in the status line, or the
                # count of scans and the count of curves disagree with nothing
                # on screen to say why.
                skipped.append(data.path.name)
                continue
            if last >= data.n_frames:
                clamped.append(data.path.name)
            if dead_time and not data.have_dead_time:
                uncorrected.append(data.path.name)
            if part != "all":
                half = data.odd_frames[first : min(last, data.n_frames - 1) + 1]
                # A difference needs both halves; the two single-half modes
                # need their own. min() of the two covers all three.
                present = (
                    min(int(half.sum()), int((~half).sum()))
                    if part == "difference"
                    else int(half.sum() if part == "odd" else (~half).sum())
                )
                if present == 0:
                    empty_half.append(data.path.name)

            totals = data.range_sum(first, last, dead_time, part)
            sums = (
                totals
                if part == "all"
                else data.range_sum(first, last, dead_time)
            )
            if separate:
                for i in usable:
                    draw(totals[i], f"element {i + 1}", sums[i])
            else:
                draw(
                    totals[usable].sum(axis=0),
                    data.path.name if len(self.data) > 1 else "sum",
                    sums[usable].sum(axis=0),
                )

        # Only when something was drawn: with every element unticked the
        # "Peak" button would otherwise search a row of zeros and answer with
        # the first bin of the view. Keeping the last real spectrum answers
        # with the peak that was on screen a moment ago, which is the one the
        # user is still looking at the band on.
        if drawn:
            self._shown_total = total

        self.ax.set_xlabel("MCA bin")
        self.ax.set_ylabel(self._count_label(part))
        # The boxes, never an autoscale. Letting matplotlib choose would move
        # the axis whenever an element was unticked, and autoscaling to the
        # last non-zero bin would always find bin 4095, where everything above
        # the ADC range piles up.
        self.ax.set_xlim(xmin, xmax)

        log = self.log_box.isChecked()
        if log:
            # An empty bin is 0, which a log scale cannot place. Clipping is
            # what hides it rather than dropping the whole curve.
            self.ax.set_yscale("log", nonpositive="clip")
        else:
            self.ax.set_yscale("linear")

        # Y from the bins on screen, not from the whole spectrum. matplotlib's
        # autoscale looks at every point of every line, so once the bin range
        # is narrowed the axis would still be stretched by whatever sits
        # outside the view -- the elastic peak, or bin 4095's overflow pile --
        # and the small peaks the range was narrowed to show would be flat
        # along the bottom. Nothing is dropped by this: the curves still run
        # the full width, the frame is just put around the chosen part.
        inside = np.concatenate(shown) if shown else np.empty(0)
        # Bins a log axis cannot place because the difference put them below
        # zero. Counted before they are dropped and reported in the status
        # line: a clipped curve looks like a curve that simply stops, and on a
        # difference that is half the measurement going quietly missing.
        hidden = int((inside < 0).sum()) if log and signed else 0
        if log:
            inside = inside[inside > 0]
        if inside.size:
            top = float(inside.max())
            if log:
                self.ax.set_ylim(
                    float(inside.min()) / Y_LOG_MARGIN, top * Y_LOG_MARGIN
                )
            elif signed:
                # Zero has to stay on the axis -- it is the line the two
                # halves are being compared across -- but the floor cannot
                # stay *at* zero or the negative lobe would be cut off. So
                # the frame spans zero and whichever ends are real, padded
                # the same fraction at both.
                low = min(0.0, float(inside.min()))
                high = max(0.0, top)
                pad = (high - low) * (Y_MARGIN - 1.0) or 1.0
                self.ax.set_ylim(low - pad, high + pad)
            elif top > 0:
                # Not when the curve is flat at zero, which a parity mode can
                # produce from a range holding none of its half: zero to zero
                # is a singular axis and matplotlib says so on the terminal.
                self.ax.set_ylim(0.0, top * Y_MARGIN)
            if signed and not log:
                self.ax.axhline(0.0, color="0.5", lw=0.6, zorder=0.5)

        if drawn:
            self.ax.legend(fontsize="small")
        self.ax.grid(True, which="major", alpha=0.3)

        # The secondary axis was built before this clear() and went with it.
        self.energy_axis = self.ax.secondary_xaxis(
            "top", functions=self._energy_map()
        )
        self.energy_axis.set_xlabel(
            "nominal energy [keV] - inferred, not recorded"
        )

        # clear() takes the axes' callbacks with it -- checked on matplotlib
        # 3.9.1 rather than assumed -- so the zoom watcher is reconnected here
        # every redraw instead of once at build time. After the set_xlim above,
        # or this redraw's own limits would come back in as if someone had just
        # zoomed to them.
        self.ax.callbacks.connect("xlim_changed", self._on_xlim_changed)

        # After the clear and after set_xlim: the band is an artist on these
        # axes and went with the one, and it is placed in data coordinates
        # that the other has just settled.
        self._make_roi_span()

        self.canvas.draw_idle()
        # The same flag the spectrum pane drew with, passed in rather than
        # recomputed, so the two panes cannot come to differ about it.
        notes = self._redraw_trend(chosen, dead_time, part, separate)
        if empty_half:
            notes.insert(
                0,
                f'no frames for "{self._mode_label(part)}" in {first}-{last}: '
                + ", ".join(empty_half),
            )
        if hidden:
            notes.insert(0, f"log Y hides {hidden} negative bins")
        self.status.showMessage(
            self._status_text(first, last, clamped, uncorrected, skipped, notes)
            if drawn
            else "Nothing to draw: no element selected that these files have."
        )

        # Last, and only while it is on screen: the frame range, the dead-time
        # tick and the gain all move what that window draws, and it reads them
        # off this one. Cheap -- the pyramid answers a frame range without
        # touching the cube -- but pointless for a window nobody is looking at.
        for window in self.monitor_windows.values():
            if window.isVisible():
                window.refresh()

    # -- the trend pane -----------------------------------------------------

    def _field_values(self, index: int, name: str) -> np.ndarray | None:
        """One master field, read once per scan and kept.

        Every tick, every drag and every redraw asks for the same two columns;
        without this they would each be a file open.
        """
        key = (index, name)
        if key not in self._field_cache:
            self._field_cache[key] = read_primary_field(
                self.data[index].master, name
            )
        return self._field_cache[key]

    def _monitor_label(self, monitor: str) -> str:
        """How the Y axis names the divisor.

        A computed monitor carries its element and its window, because a
        figure saved from the toolbar is the only record of which region was
        chosen -- the range is not written to the file, and with two of these
        on offer the name alone does not say which band was in force.
        """
        if monitor not in COMPUTED_MONITORS:
            return monitor
        lo, hi = self.monitor_ranges_kev[monitor]
        element = COMPUTED_MONITORS[monitor]
        return f"{monitor} (el {element}, {lo:.3f}-{hi:.3f} keV)"

    def _recorded_scale(self, index, data, monitor):
        """A monitor read from the master, spread from its points over frames.

        One value per scan point, so every frame of a point gets the same
        number. ``(scale, note)``, or ``(None, why)`` when the pairing cannot
        be made.
        """
        mon = self._field_values(index, monitor)
        if mon is None:
            return (
                None,
                f"{data.path.name}: no {monitor} recorded, not normalized",
            )
        index_of = data.point_index
        if index_of is None:
            return (
                None,
                f"{data.path.name}: frames cannot be paired with {monitor}",
            )

        scale = np.full(data.n_frames, np.nan)
        inside = index_of < mon.size
        scale[inside] = mon[index_of[inside]]
        note = (
            None
            if inside.all()
            else f"{data.path.name}: {int((~inside).sum())} frames past the "
            f"last recorded {monitor}"
        )
        return scale, note

    def _computed_scale(self, data, monitor, part, dead_time, window=None):
        """A computed monitor, reduced over the frames the trend reduces.

        A monitor in :data:`COMPUTED_MONITORS` is not a recorded channel: it
        is that element of this cube over the window its own monitor-ROI
        dialog chose. That makes it a per-*frame* number, and the granularity
        it is reduced to before dividing is the whole of whether the Combine
        box still says anything.

        Reduced to one value per scan point **per parity**, over exactly the
        frames on show, and then spread back over those frames. *window* is
        the shown frame range as ``(first, stop)`` inclusive, defaulting to
        all of them. It has to be known here: a point only half inside the
        range would otherwise have its signal totalled over the shown frames
        and its monitor over all of them, which is a 0.1-0.2% error on the
        partial point of scan 597 at frames 10-5455.

        Reducing this way is what lets :meth:`_trend_points` total the half's
        spectra, total this monitor over those same frames, and divide the
        two totals. That is now the arithmetic the pane performs, rather than
        a different one that happened to agree with it.

        What the rule does to a difference is worth knowing before reading
        one. Either computed monitor is an element of the cube the signal
        came out of, so it alternates itself, and dividing each half by its
        own parity carries that alternation into the result. Measured on
        scan 597 at the 7.41 keV line, which reads -2.048% raw: element 1
        alternates -2.414% over the band both monitors open on and -2.379%
        over the whole spectrum, so dividing by it gives +0.366% and
        +0.333% -- the asymmetry pulled most of the way out, since element 1
        is where the signal is. Element 7 holds almost none of the line, but
        it alternates about three times harder (-7.396% over that band,
        -7.134% whole) and so overshoots the other way, to +6.621% and
        +6.189%: the sign flips.
        Dividing both halves by the point's whole monitor would leave it at
        -2.053%, but that is not what was asked for: "sum ROI from odd
        numbered file, and the same for the even numbered". Neither element
        is a quiet reference; read the raw columns the export writes beside
        the normalized ones to see what the division did.

        Spreading back over frames rather than returning one value per point
        is what lets :meth:`_trend_points` take this exactly as it takes the
        values: it totals the divisor over a point's frames and divides the
        totalled signal by it, once. So Sum divides each point by all of its
        frames, and Odd only, Even only and Odd - even divide each half by
        its own half, without a second reduction table to keep in step with
        the first.

        **A mean here, and not a sum** -- even though the pane now plots
        sums. What this returns is a per-*frame* divisor, and the reduction
        totals it over the point's ``n`` frames to reach ``n*m``, which is
        the denominator wanted. Summing here as well would make it
        ``n**2 * m``. The mean is not a leftover of the days when the pane
        averaged, and turning it into a sum is not the tidy-up it looks like.
        """
        element = COMPUTED_MONITORS[monitor]
        if data.n_elements < element:
            return None, (
                f"{data.path.name}: no element {element} to take "
                f"{monitor} from, not normalized"
            )
        lo, hi = self.monitor_bins(monitor)
        # A cache slot per monitor, so that neither range evicts the other's
        # nor the ROI's -- each eviction is a full pass over the cube on every
        # redraw, and with two computed monitors a single shared "monitor"
        # slot would thrash between their two bands.
        raw, corrected = data.roi_trend(lo, hi, slot=f"monitor:{monitor}")
        use = corrected if dead_time and corrected is not None else raw
        per_frame = use[:, element - 1]

        index_of = data.point_index
        if index_of is None:
            # No frame-to-point map, so the frame is the point and there is
            # nothing coarser to average over. _trend_points plots by frame in
            # this case too, so the two still agree.
            return per_frame, None

        points = int(index_of[-1]) + 1
        odd = data.odd_frames if part != "all" else None
        shown = np.zeros(data.n_frames, dtype=bool)
        first, stop = (0, data.n_frames - 1) if window is None else window
        shown[first : stop + 1] = True
        halves = [shown] if odd is None else [odd & shown, ~odd & shown]
        scale = np.full(data.n_frames, np.nan)
        for keep in halves:
            idx = index_of[keep]
            n = np.bincount(idx, minlength=points).astype(float)
            summed = np.bincount(idx, weights=per_frame[keep], minlength=points)
            # NaN where this parity has no frame at a point, matching what
            # _trend_points does with the values themselves: nothing was
            # measured there, so there is no monitor reading either.
            mean = np.where(n > 0, summed / np.maximum(n, 1.0), np.nan)
            scale[keep] = mean[idx]
        return scale, None

    def _monitor_scale(
        self, index, data, monitor, part="all", dead_time=False, window=None
    ):
        """The per-frame divisor for one scan. ``(scale, note)``.

        It does not divide anything. *scale* is one number per frame of the
        whole scan, for the caller to hand to :meth:`_trend_points`, which
        totals it over a point's frames and divides the totalled signal by
        the total once. Numerator and denominator reduced over the same
        frames and divided at the end -- not divided frame by frame and then
        added, which would multiply every point by its own frame count.

        ``(None, why)`` when nothing can be divided: no such channel, no
        frame-to-point map, or a monitor that is zero throughout. The caller
        then plots raw and says so, and -- this is why both panes come
        through here -- the exported ``mon_`` column drops out at exactly the
        same moment the normalized one does.

        *window*, the shown frame range as ``(first, stop)`` inclusive, is
        passed on to :meth:`_computed_scale` so that a point only half inside
        the range divides by the monitor of the frames that are shown. A
        recorded monitor has one value per point however many frames are on
        show, so it ignores the window.

        Where the divisor comes from is the difference between the two kinds
        of monitor -- see :meth:`_recorded_scale` and :meth:`_computed_scale`.
        """
        if monitor in COMPUTED_MONITORS:
            scale, note = self._computed_scale(
                data, monitor, part, dead_time, window
            )
        else:
            scale, note = self._recorded_scale(index, data, monitor)
        if scale is None:
            return None, note

        # A zero or negative monitor is a point that took no beam. NaN leaves
        # a gap in the curve; dividing would throw it off the top of the axis
        # and rescale everything else into a flat line.
        scale = np.where(scale > 0, scale, np.nan)
        if not np.any(np.isfinite(scale)):
            # Nothing would be left to plot. The picker only offers monitors
            # with beam in them, so this is a scan whose monitor died partway
            # through a window opened on several -- or, for the computed one,
            # a monitor window sitting on a part of the spectrum with no
            # counts in it at all.
            return (
                None,
                f"{data.path.name}: {monitor} is zero throughout, not normalized",
            )
        return scale, note

    @staticmethod
    def _frame_pairs(values, index_of, odd, first):
        """Each odd frame against the even frame beside it, within its point.

        The per-frame reading of "odd - even", and the only one the frame axis
        can carry: a point's mean is not a frame. Pairing is by position
        within the half, so frame 1 goes with frame 2, 3 with 4, and a point
        whose halves came out unequal drops its unpartnered tail rather than
        pairing across a point boundary.

        Returns ``(x, y)`` with *x* in the scan's own frame numbers.
        """
        rows = np.arange(values.shape[0])
        xs, ys = [], []
        # At most a few dozen points, so a loop costs less than the index
        # arithmetic that would avoid it, and says plainly what it pairs.
        for p in np.unique(index_of):
            here = index_of == p
            o = rows[here & odd]
            e = rows[here & ~odd]
            k = min(o.size, e.size)
            if k:
                xs.append(o[:k] + first)
                ys.append(values[o[:k]] - values[e[:k]])
        if not xs:
            return np.empty(0), np.empty((0, values.shape[1]))
        return np.concatenate(xs).astype(float), np.concatenate(ys)

    def _trend_points(
        self, index, data, values, first, kind, field, part="all", scale=None
    ):
        """``(x, y, xlabel, notes)`` for one scan: the ROI against the axis.

        *values* is per frame and already cut to the frames on show; *first*
        says which frame of the scan its first row is, so the frames can be
        put back against the points they belong to. *y* comes back per scan
        point, **totalled** over the shown frames the detector took there --
        one frame on every motor scan and on the whole 2025-1 vintage, 20 or
        80 on the 2026-3 energy scans, and fewer than that at either end of a
        restricted frame range.

        *scale* is the per-frame divisor from :meth:`_monitor_scale`, cut to
        the same frames as *values*, or ``None`` to total without dividing.
        It is totalled over each point's frames alongside the values and
        divided in once at the end, so what comes back is a ratio of two
        totals over the same frames. Passing values that have already been
        divided, and a *scale* as well, would divide twice.

        *part* is the Combine box: with anything but ``"all"`` the frames of a
        point are split by :attr:`SpectrumData.odd_frames` and each half
        totalled on its own, the difference being of the two sums -- which is
        what the spectrum pane next door has always shown, so the two panes
        now agree. The cost is that a point which ended up with 44 odd frames
        and 43 even ones reports the extra frame as signal in the *raw*
        difference; a note below says so, and the normalized difference
        divides each half by its own monitor and is unaffected.

        The result is as long as the whole scan, not as long as the part of it
        that has frames: a point outside the window comes back NaN. That is
        what keeps several scans averageable point by point, and what keeps
        the y values paired with the x column the master recorded.
        """
        n_frames = values.shape[0]
        index_all = data.point_index
        full_mask = data.odd_frames if part != "all" else None
        odd = None if full_mask is None else full_mask[first : first + n_frames]

        if scale is not None and (kind == "frame" or index_all is None):
            # One row per frame on these axes, so there is no point's worth of
            # frames to total over and the divisor is simply applied. Clearing
            # *scale* afterwards is what stops the reduction below -- which
            # these branches do not reach today, but would if either grew a
            # path to it -- from dividing a second time.
            values = values / scale[:, None]
            scale = None

        if kind == "frame":
            # Numbered from *first*, not from zero: these are the scan's own
            # frame numbers, the ones the two frame boxes are set to.
            frames = np.arange(first, first + n_frames, dtype=float)
            if odd is None:
                return frames, values, TREND_FRAME_AXIS, []
            if part == "difference":
                x, y = self._frame_pairs(
                    values, index_all[first : first + n_frames], odd, first
                )
                return x, y, TREND_FRAME_AXIS, []
            keep = odd if part == "odd" else ~odd
            return frames[keep], values[keep], TREND_FRAME_AXIS, []

        if index_all is None:
            return (
                np.arange(first, first + n_frames, dtype=float),
                values,
                TREND_FRAME_AXIS,
                [
                    f"{data.path.name}: no frame-to-point map, so plotted by frame"
                ],
            )

        notes = []
        points = int(index_all[-1]) + 1
        index_of = index_all[first : first + n_frames]

        def totals_over(keep):
            """Total per scan point over the frames *keep* selects.

            Divided by the total of *scale* over those same frames when there
            is one. NaN where no frame of a point is in the window, or where
            the monitor was not known over the frames that are.
            """
            idx = index_of if keep is None else index_of[keep]
            rows = values if keep is None else values[keep]
            n = np.bincount(idx, minlength=points).astype(float)
            # bincount per element rather than np.add.at, which is an order
            # slower; there are at most seven of them.
            summed = np.stack(
                [
                    np.bincount(idx, weights=rows[:, e], minlength=points)
                    for e in range(values.shape[1])
                ],
                axis=1,
            )
            # NaN and not zero where the window holds no frame of a point:
            # nothing was measured there, and a zero would be read as a
            # measurement of no counts -- and would drag the mean over scans
            # down with it. A half that came out empty is the same case.
            if scale is None:
                return np.where(n[:, None] > 0, summed, np.nan), n
            den = np.bincount(
                idx,
                weights=scale if keep is None else scale[keep],
                minlength=points,
            )
            # A NaN denominator compares False here, which is what is wanted:
            # the monitor was not known over those frames, so the ratio is not
            # known either. The substitution has to sit *inside* the division
            # -- moved outside it, a NaN denominator would hand back the
            # undivided total, which plots as a spike of raw counts.
            good = (n > 0) & (den > 0)
            return (
                np.where(
                    good[:, None],
                    summed / np.where(good, den, 1.0)[:, None],
                    np.nan,
                ),
                n,
            )

        lopsided = 0
        if odd is None:
            grouped, counts = totals_over(None)
        elif part == "difference":
            odd_sum, odd_n = totals_over(odd)
            even_sum, even_n = totals_over(~odd)
            # NaN propagates on its own where either half is empty, which is
            # the behaviour wanted: a point measured on one side only is not a
            # difference, and plotting it as the half it has would read as a
            # huge asymmetry.
            grouped = odd_sum - even_sum
            counts = odd_n + even_n
            lopsided = int(
                np.count_nonzero((odd_n > 0) & (even_n > 0) & (odd_n != even_n))
            )
        else:
            grouped, counts = totals_over(odd if part == "odd" else ~odd)

        whole = np.bincount(index_all, minlength=points).astype(float)
        if odd is not None and part in ("odd", "even"):
            # Half the frames are the whole of this half, or every point would
            # be reported as partly averaged.
            whole = np.bincount(
                index_all[full_mask if part == "odd" else ~full_mask],
                minlength=points,
            ).astype(float)
        partial = int(np.count_nonzero((counts > 0) & (counts < whole)))
        if partial:
            # Not an error: the total of the frames that are in the window is
            # exactly what was asked for. But on totals this is a step and no
            # longer a wobble -- a point built from one frame sits eighty
            # times lower than its neighbours built from eighty, rather than
            # merely nine times noisier. The normalized curve divides by the
            # monitor over those same frames and is level regardless.
            notes.append(
                f"{data.path.name}: {partial} point(s) total fewer frames "
                f"than the rest, so their raw values sit low in proportion"
            )
        if lopsided and scale is None:
            # On totals, a point with one more odd frame than even carries
            # that whole frame's counts into the difference as signal -- of
            # order a thousand counts against a real asymmetry of tens. The
            # partial test above cannot see it: 44 and 43 frames is 87 of 87.
            # Raw only, hence the scale check: the normalized difference
            # divides each half by its own half's monitor, and the frame
            # count cancels.
            notes.append(
                f"{data.path.name}: {lopsided} point(s) hold unequal numbers "
                f"of odd and even frames, so their raw difference carries the "
                f"extra frame as signal; the normalised difference does not"
            )

        if kind == "point":
            return (
                np.arange(points, dtype=float),
                grouped,
                TREND_POINT_AXIS,
                notes,
            )

        x = self._field_values(index, field)
        if x is None:
            return (
                np.arange(points, dtype=float),
                grouped,
                TREND_POINT_AXIS,
                notes + [f"{data.path.name}: no {field} recorded"],
            )
        # The detector and the scan need not have written the same number of
        # points: scan 125 was aborted after 9 rows with ten complete 80-frame
        # blocks already on disk. Pair what pairs and say what was left over,
        # rather than let numpy broadcast two different scans together.
        n = min(points, x.size)
        if points != x.size:
            notes.append(
                f"{data.path.name}: {points} detector points against "
                f"{x.size} recorded, pairing the first {n}"
            )
        return x[:n], grouped[:n], field, notes

    def _trend_settings(self) -> dict:
        """The widget state one trend pass reads, gathered once.

        A snapshot and not a set of live reads, so that a pass which outlives
        the moment it started -- the export, whose file dialog runs a nested
        event loop -- cannot have the numbers change under it and end up
        disagreeing with its own header.
        """
        normalize = (
            self.normalise_box.isChecked() and self.normalise_box.isEnabled()
        )
        kind, field = self.x_combo.currentData()
        low, high = self._roi_bins()
        return {
            "first": self.first_spin.value(),
            "last": self.last_spin.value(),
            "kind": kind,
            "field": field,
            "normalize": normalize,
            "monitor": self.monitor_combo.currentText() if normalize else "",
            "low": low,
            "high": high,
        }

    def _trend_curves(self, chosen, dead_time, part, separate, settings):
        """The curves the trend pane would draw, reduced but not drawn.

        ``(mean_x, curves, xlabel, notes, used)``. *curves* maps an element
        index -- or :data:`SUM_KEY` for the one summed curve -- to its y
        values, already averaged over the scans on show. *used* is how many
        scans went into that average.

        Split out of :meth:`_redraw_trend` so the export writes the numbers
        the pane plots instead of deriving them a second time. The per-parity
        rule this reaches through :meth:`_monitor_scale` and
        :meth:`_trend_points` is the one piece of arithmetic here that is
        both easy to get wrong and expensive to get wrong, so it exists once.

        Calling it repeatedly is cheap: :meth:`SpectrumData.roi_trend` caches
        on the bin range per slot, so the eight passes the export makes --
        four Combine modes, raw and normalized -- read the cube once between
        them and are masks and :func:`numpy.bincount` thereafter.
        """
        first, last = settings["first"], settings["last"]
        kind, field = settings["kind"], settings["field"]
        normalize, monitor = settings["normalize"], settings["monitor"]
        low, high = settings["low"], settings["high"]

        notes: list[str] = []
        xlabel = TREND_POINT_AXIS
        # Element number -> that element's curve from each scan that has it,
        # to be averaged below. A dict and not an array because the scans need
        # not have the same number of elements, and element 5 of a 4-element
        # file is not a zero to be averaged in: it is a scan with nothing to
        # say about element 5.
        per_element: dict[int, list[np.ndarray]] = {}
        xs: list[np.ndarray] = []
        used = 0
        # Of the scans that actually contributed, not of everything loaded:
        # a scan with none of its frames in the window says nothing about the
        # scale the others are averaged on.
        per_point: set = set()

        for index, data in enumerate(self.data):
            usable = [i for i in chosen if i < data.n_elements]
            if not usable:
                continue
            stop = min(last, data.n_frames - 1)
            if stop < first:
                notes.append(f"{data.path.name}: no frames in {first}-{last}")
                continue
            raw, corrected = data.roi_trend(low, high)
            values = corrected if dead_time and corrected is not None else raw
            scale = None
            if normalize:
                # Built over the whole scan, not the cut: _monitor_scale pairs
                # each row with a point through the scan's full-length
                # point_index. The shown range goes with it, so that a point
                # only half inside the window divides by the monitor of the
                # frames it was actually totalled over. Both it and the values
                # are cut to the window below, in step.
                scale, note = self._monitor_scale(
                    index, data, monitor, part, dead_time, (first, stop)
                )
                if note:
                    notes.append(note)
            x, y, xlabel, scan_notes = self._trend_points(
                index,
                data,
                values[first : stop + 1],
                first,
                kind,
                field,
                part,
                None if scale is None else scale[first : stop + 1],
            )
            notes.extend(scan_notes)
            xs.append(x)
            if separate:
                for i in usable:
                    per_element.setdefault(i, []).append(y[:, i])
            else:
                # Summed here, per scan, and averaged across scans below --
                # the same order the spectrum pane sums in, and the order
                # that keeps a scan missing an element out of the others'
                # sum. NaN propagates through the sum on purpose.
                per_element.setdefault(SUM_KEY, []).append(
                    y[:, usable].sum(axis=1)
                )
            used += 1
            per_point.add(data.frames_per_point)

        if xs:
            step_off = [
                self.data[i].path.name
                for i, x in enumerate(xs[1:], start=1)
                if not same_grid(xs[0], x)
            ]
            if step_off:
                # Averaged anyway -- the points still pair up by index, which
                # is what the average is over -- but an x that is the mean of
                # two different grids is not a position any scan visited, and
                # saying nothing would hide that.
                notes.append(
                    "x values differ between scans, averaged by point: "
                    + ", ".join(step_off)
                )
        if not normalize and len(per_point - {None}) > 1:
            # Totals, so a scan taking 80 frames a point reads eighty times a
            # scan taking one, and the mean below lands wherever the longer
            # one is rather than between them. Harmless while the pane
            # averaged; worth saying now. Normalized curves are ratios and
            # carry no such scale, which is why this is raised only for raw.
            notes.append(
                "scans hold different numbers of frames per point "
                f"({', '.join(str(n) for n in sorted(per_point - {None}))}), "
                "so their raw totals are on different scales and the mean "
                "across scans leans to the longer -- normalise to compare them"
            )
        curves = {i: stack_mean(per_element[i]) for i in sorted(per_element)}
        return stack_mean(xs), curves, xlabel, notes, used

    def _monitor_trend(self, dead_time, part, settings):
        """The monitor itself, totalled per scan point over *part*'s frames.

        The divisor is pushed through :meth:`_trend_points` exactly as the
        data is -- same function, same *part*, same window, same mask, and
        the same :func:`numpy.bincount` over the same array -- so the column
        written beside the normalized values is the denominator that actually
        divided them. Dividing the raw column by this one reproduces the
        normalized column by construction rather than by agreement.

        Totalled, and so no longer the monitor's own reading at a point: on a
        recorded monitor of ``m`` per point it is ``n*m`` over the point's
        ``n`` shown frames. That is the number the raw total was divided by,
        which is what this column is for; the per-point reading is in the
        master file.

        One value per scan point, averaged over the scans on show. ``None``
        when nothing is being normalized, or when no scan can offer the
        monitor.

        *part* is ``"all"``, ``"odd"`` or ``"even"`` only. There is no
        divisor for ``"difference"``: each half is divided by its own
        parity's monitor, so the honest answer for that mode is the odd and
        even columns side by side.
        """
        if not settings["normalize"]:
            return None
        first, last = settings["first"], settings["last"]
        monitor = settings["monitor"]
        per_scan = []
        for index, data in enumerate(self.data):
            stop = min(last, data.n_frames - 1)
            if stop < first:
                continue
            # The same call _trend_curves makes, so the zero guard and the
            # "zero throughout" refusal are the same ones -- a point that took
            # no beam is a hole in this column too rather than the number that
            # was never divided by, and a monitor that never had any drops
            # this column instead of leaving raw numbers under a _norm name.
            scale, _note = self._monitor_scale(
                index, data, monitor, part, dead_time, (first, stop)
            )
            if scale is None:
                continue
            # No *scale* of its own: this column is the divisor, so it is
            # totalled over each point's frames exactly as the values are.
            # That is what makes raw / mon reproduce the normalized column.
            _x, y, _label, _notes = self._trend_points(
                index,
                data,
                scale[first : stop + 1, None],
                first,
                settings["kind"],
                settings["field"],
                part,
            )
            per_scan.append(y[:, 0])
        return stack_mean(per_scan) if per_scan else None

    def _redraw_trend(
        self, chosen, dead_time, part="all", separate=True
    ) -> list[str]:
        """Draw the ROI against the scan axis. Returns notes for the status bar.

        With *separate* -- the "Separate" tick, as the spectrum pane read it
        -- one curve per ticked element, each the mean over the frames on show
        and then over the scans on show. Kept as the default because the
        elements are different detector channels looking at the same sample
        from different angles, and which of them sees a line is a comparison
        a sum throws away. Repeated scans are the thing worth averaging, and
        they are averaged here rather than overlaid -- seven elements times
        three scans would be twenty-one curves, and the scans are already
        overlaid in the spectrum pane next door.

        Without it, one curve labelled ``"sum"``: the ticked elements added
        together, matching the single summed curve the spectrum pane draws in
        the same state. Summed across elements first and averaged across
        scans after, which is the order the pane next door uses, and which
        stops a scan that is missing an element from quietly contributing a
        sum of fewer channels. A NaN in any ticked element makes the sum NaN,
        because the sum of an unknown is unknown.

        The reducing is :meth:`_trend_curves`, which the export calls too;
        what is left here is the pane's own business -- whether it is on
        screen, how wide it opens, and the drawing. *part*, the Combine box,
        goes straight through.

        No recompute/replot split is needed here: :meth:`SpectrumData.roi_trend`
        is cached on the ROI bins, so a change of element, dead time, monitor
        or axis finds the array already there and never opens the file. Only
        moving the ROI -- or the gain or offset under it, which moves the bins
        the keV resolves to -- costs a pass.
        """
        plotting = self.plot_roi_btn.isChecked()
        was_visible = self.trend_pane.isVisible()
        self.trend_pane.setVisible(plotting)
        if not plotting:
            # Before anything is read. This early return is the whole point of
            # the button: roi_trend streams the cube, about a second on the
            # 8000-frame scan, and a window opened to look at a spectrum
            # should not pay for a plot nobody asked for.
            return []
        if not was_visible:
            # Only on the edge, so that a handle the user has since dragged is
            # not shoved back to 2:1 by the next redraw. setMinimumWidth(0)
            # first, or Qt's implicit minimum for the toolbar can overrule the
            # split and the pane comes up wider than a third.
            width = self.plot_split.width()
            if width > 0:
                self.trend_pane.setMinimumWidth(0)
                self.plot_split.setSizes(
                    [2 * width // 3, width - 2 * width // 3]
                )

        settings = self._trend_settings()
        normalize, monitor = settings["normalize"], settings["monitor"]
        low, high = settings["low"], settings["high"]
        mean_x, curves, xlabel, notes, used = self._trend_curves(
            chosen, dead_time, part, separate, settings
        )

        self.trend_ax.clear()
        drawn = 0
        for i, mean_y in curves.items():
            # The element may be missing from the longest scan, leaving its
            # mean shorter than the shared x. Pair what pairs.
            n = min(mean_x.size, mean_y.size)
            if n == 0:
                continue
            # Markers say where the points actually are, which on a scan of
            # forty is the useful part; on eight thousand frames they are a
            # solid bar and slow to draw.
            self.trend_ax.plot(
                mean_x[:n],
                mean_y[:n],
                lw=0.9,
                marker="." if n <= 400 else None,
                ms=3,
                label="sum" if i == SUM_KEY else f"element {i + 1}",
            )
            drawn += 1

        self.trend_ax.set_xlabel(xlabel)
        # "sum", because each point is a total over its frames: the number
        # grows with the frames on show, and two scans of different dwell are
        # not comparable raw -- which is what the monitor is for, and why the
        # divisor is named beside it. Under a contrast it says which frames
        # it is a sum of, since the difference of two sums is the one number
        # on this pane a reader could otherwise take for a count.
        label = f"sum ROI {self._count_label(part)}"
        self.trend_ax.set_ylabel(
            f"{label} / {self._monitor_label(monitor)}" if normalize else label
        )
        # Both readings of the ROI, so that neither can be doubted: the energy
        # it was set to, and the bins that energy currently resolves to.
        title = (
            f"ROI {self.roi_center_kev:.3f} ± "
            f"{self.roi_width_ev / 1000.0:.3f} keV = bins {low}-{high}"
        )
        if part != "all":
            # On the title and not only on the axis, because this pane is
            # saved as a figure from the toolbar and a difference that lost
            # its caption would be read as a count.
            title += ", " + self._mode_label(part).lower()
        if used > 1:
            title += f", mean of {used} scans"
        self.trend_ax.set_title(title, fontsize="small")
        self.trend_ax.grid(True, which="major", alpha=0.3)
        if drawn:
            self.trend_ax.legend(fontsize="small")
        self.trend_canvas.draw_idle()
        # Said once per monitor, the first time a window actually divides by a
        # monitor band nobody picked. The default is wide enough to be a
        # near-total, so it no longer divides one line by itself -- but it
        # does contain the ROI, and every other line of the element with it,
        # which is worth one sentence before a number is believed. The y-label
        # names the window every time; this names the *choice*, which the
        # label cannot. Not said at construction, where it would be noise to
        # anyone who never normalises, and not said again once it has been
        # read. Per monitor and not once overall, because switching to the
        # other one lands on a second unchosen band.
        if (
            normalize
            and monitor in COMPUTED_MONITORS
            and self._monitor_default_note.get(monitor)
            and self._monitor_window_is_default(monitor)
        ):
            self._monitor_default_note[monitor] = False
            lo_bin, hi_bin = self._default_monitor_bins()
            lo_kev, hi_kev = self.monitor_ranges_kev[monitor]
            notes.append(
                f"{monitor} window defaulted to bins {lo_bin}-{hi_bin} "
                f"({lo_kev:.3f}-{hi_kev:.3f} keV) of element "
                f"{COMPUTED_MONITORS[monitor]} -- "
                "check it is not the line you are measuring"
            )
        return notes

    def _monitor_window_is_default(self, monitor: str) -> bool:
        """Whether one monitor's window is still :data:`DEFAULT_MONITOR_BINS`.

        Compared in bins and not in keV, which is both simpler -- no slop,
        because :meth:`monitor_bins` is the same rounding the note quotes --
        and the only version that cannot print a false bin range. A user who
        moved the window and happened to land back on it is told nothing:
        they have seen the dialog, which is the point of the note.

        A recalibration therefore silences the note, because the band is held
        in keV and a new gain walks it off bins 200-1800 -- the same way it
        walks the signal ROI off its line. The window is still one nobody
        chose, so this loses a true warning; accepted because the y-label and
        the export header name the window on every single redraw, and someone
        who has just retyped the gain is already looking at the energy axis.
        """
        return self.monitor_bins(monitor) == self._default_monitor_bins()

    # -- the export ---------------------------------------------------------

    def _export_directory(self) -> Path:
        """Where the save dialog opens: ``export`` beside the scan data.

        The sample directory, so the exports sit next to ``vortex`` rather
        than inside it -- the spectrum files are the instrument's, the export
        is the user's. Derived from the data and not remembered in QSettings:
        ``scan_viewer`` tried a remembered export directory and dropped it
        again, because a run's exports belong with the run.
        """
        master = self.data[0].master
        base = (
            master.parent
            if master is not None
            else self.data[0].path.parent.parent
        )
        return base / "export"

    def _export_folder(self) -> Path:
        """:meth:`_export_directory`, made -- or its parent if it cannot be.

        A directory that is not ours to write in is said once and then
        stepped out of, so the dialog opens on a path that is there and the
        user can steer it anywhere from. The proposal roots here are
        read-only, so this is a live case and not a defensive one.
        """
        folder = self._export_directory()
        try:
            folder.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            self.status.showMessage(
                f"Could not make {folder}  -  {exc.strerror or exc}"
            )
            return folder.parent
        return folder

    def _export_stem(self) -> str:
        """``data_000597``, or ``data_000597_000598`` for two scans."""
        paths = [data.path for data in self.data]
        tags = scan_tags(paths)
        parts = [tags[str(path)] for path in paths]
        # A dozen scans loaded at once would otherwise name every one of them.
        if len(parts) > 3:
            parts = [parts[0], "to", parts[-1]]
        return "data_" + "_".join(_safe_name(part) for part in parts)

    def _export_settings(self):
        """:meth:`_trend_settings`, with the frame axis swapped for a point one.

        ``(settings, note)``. The four Combine modes do not share rows on the
        frame axis -- Sum has every frame, Odd and Even half each, Odd - even
        one row per pair -- so a table with one x column can only be honest
        per scan point. The pane is left alone; the file says what it did.
        """
        settings = self._trend_settings()
        if settings["kind"] != "frame":
            return settings, None
        if self.hinted_axis in self.axis_fields:
            settings["kind"], settings["field"] = "field", self.hinted_axis
            against = self.hinted_axis
        else:
            settings["kind"], settings["field"] = "point", ""
            against = TREND_POINT_AXIS
        return settings, (
            f"the pane was showing the frame axis; exported against "
            f"{against}, since the four Combine modes do not share rows "
            f"per frame"
        )

    def _normalizable(self, dead_time, settings) -> bool:
        """Whether any scan on show can actually be divided by the monitor.

        The picker offers monitor channels that *exist*, not channels with
        beam in them over the window chosen, so "/ monitor" can be ticked and
        every scan still come back undivided -- a monitor ROI parked on an
        empty part of the spectrum is the easy way to do it. The export asks
        this before it names a column ``_norm``, so that a column with that
        name is always a divided one and never a copy of the raw beside it.

        Asked with ``part="all"``: a parity is a subset of the frames, so if
        no frame in the window has a positive monitor, no half has one
        either. Cheap to ask -- :meth:`SpectrumData.roi_trend` caches the
        monitor window, and the export is about to read it anyway.
        """
        if not settings["normalize"]:
            return False
        first, last = settings["first"], settings["last"]
        for index, data in enumerate(self.data):
            stop = min(last, data.n_frames - 1)
            if stop < first:
                continue
            scale, _note = self._monitor_scale(
                index,
                data,
                settings["monitor"],
                "all",
                dead_time,
                (first, stop),
            )
            if scale is not None:
                return True
        return False

    def _export_header(self, settings, separate, dead_time, axis_note, xlabel):
        """The provenance lines above the columns.

        Enough to say what the file is without the window that wrote it:
        which scans, which frames, which ROI in both units, and -- because
        these are the two things the files themselves do not record -- the
        energy map the bins were resolved through and the monitor window the
        divisor was taken over.
        """
        paths = [data.path for data in self.data]
        frames = [data.n_frames for data in self.data]
        per_point = {data.frames_per_point for data in self.data}
        low, high = settings["low"], settings["high"]
        elements = [i + 1 for i in self._selected_elements()]
        lines = [
            "scan_viewer_sgz export -- the ROI trend behind Plot ROI",
            f"written = {datetime.now().astimezone().isoformat(timespec='seconds')}",
            f"directory = {paths[0].parent}",
            f"files = {', '.join(path.name for path in paths)}",
            f"frames = {settings['first']}-{settings['last']} of "
            f"{', '.join(str(n) for n in frames)}"
            + (
                f", {per_point.pop()} per scan point"
                if len(per_point) == 1 and None not in per_point
                else ""
            ),
            f"roi = {self.roi_center_kev:.3f} +- "
            f"{self.roi_width_ev / 1000.0:.3f} keV = bins {low}-{high}",
            f"energy axis = {self.gain_spin.value():.3f} eV/bin, offset "
            f"{self.offset_spin.value():.1f} eV (nominal, not recorded)",
            "elements = "
            + (
                ", ".join(str(e) for e in elements)
                + (" separate" if separate else " summed")
                if elements
                else "(none)"
            ),
            f"dead time = {'corrected' if dead_time else 'uncorrected'}",
        ]
        if self._normalizable(dead_time, settings):
            lines.append(
                f"monitor = {self._monitor_label(settings['monitor'])}"
            )
            lines.append(
                "values = raw and normalised; both are totals over each "
                "point's frames, not means, and so is every mon_ column"
            )
            # The one thing a reader of an odd - even column has to know, and
            # the file is where they will be reading it. It is a different
            # thing for each kind of monitor, so it is not one sentence: the
            # computed one is taken from the same detector and alternates
            # with the signal, a recorded one is a single number per point
            # and cannot.
            lines.append(
                "note = oddeven_norm is exactly odd_norm - even_norm, and has "
                "no single divisor of its own, so there is no mon_oddeven"
            )
            if settings["monitor"] in COMPUTED_MONITORS:
                # What this costs is different for each element, so the note
                # is too. Measured at the 7.41 keV line of scan 597, which
                # reads -2.05% raw. Both monitors open on the same near-total
                # band, so the element is what separates them: element 1
                # alternates -2.41% over it and reads +0.37%, element 7
                # alternates -7.40% and reads +6.62% -- the sign flips. The
                # band does contain the ROI on either element, but it is 96%
                # of the element, so element 1 no longer divides the line by
                # itself the way its old one-peak window did.
                element = COMPUTED_MONITORS[settings["monitor"]]
                lines.append(
                    "note = each parity is divided by its own parity's "
                    f"monitor; this one is element {element} of the same "
                    "detector, which alternates with the frames itself, so "
                    "the odd/even asymmetry below is the ROI's against that "
                    "element's and not against a steady beam"
                )
                lines.append(
                    "note = both computed monitors open on the same band, so "
                    "the element is the choice: element 1 carries the signal "
                    "and pulls the asymmetry towards its own, element 7 "
                    "holds almost none of the line but alternates about "
                    "three times harder and reverses the sign. Neither is a "
                    "quiet reference -- there is no region of this detector "
                    "that does not alternate, and the opening band contains "
                    "the ROI on either element"
                )
            else:
                # These were equal to each other while the pane reduced by
                # the mean: one reading per point, spread over the frames and
                # averaged straight back out. On totals each is that reading
                # times the frames its own half contributed, so they part
                # company at exactly the ragged points -- and that is the
                # division doing its job, not an inconsistency.
                lines.append(
                    "note = one recorded reading per scan point, so mon_odd "
                    "and mon_even are that reading times the frames each half "
                    "totalled: equal wherever the halves are, and normalising "
                    "then leaves the odd/even asymmetry alone; where they "
                    "differ, normalising divides the extra frame back out"
                )
        elif settings["normalize"]:
            # Asked for and not delivered. Saying "(none)" here would be true
            # of the numbers and a lie about the session: the box was ticked,
            # a monitor was named, and it came back empty over the frames
            # shown. The file has to say which one failed, or the next reader
            # re-picks the same window and gets the same silence.
            lines.append(
                f"monitor = {self._monitor_label(settings['monitor'])} "
                "-- zero over every frame shown, so nothing was divided"
            )
            lines.append(
                "values = raw only; totals over each point's frames, not means"
            )
        else:
            lines.append("monitor = (none)")
            lines.append(
                "values = raw only; totals over each point's frames, not means"
            )
        lines.append(f"x axis = {xlabel}")
        if axis_note:
            lines.append(f"note = {axis_note}")
        return lines

    def _export_columns(self, chosen, dead_time, separate, settings):
        """``(names, columns, notes, xlabel)``: the whole table.

        All four Combine modes, raw and normalized, whatever the box is set
        to -- one file per scan that answers every later question, rather
        than four exports to remember to take. It costs one pass over the
        cube between them: :meth:`SpectrumData.roi_trend` caches on the bin
        range, so the seven passes after the first are masks and
        :func:`numpy.bincount`.
        """
        raw_settings = dict(settings, normalize=False, monitor="")
        # A ``_norm`` column is a promise that something was divided. When no
        # scan on show can be, the normalised pass would write the raw totals
        # a second time under that name, beside the all-NaN mon_ columns that
        # _monitor_trend correctly declines to emit -- the same numbers
        # labelled two different things. Drop the pass instead;
        # :meth:`_export_header` asks the same question and says so.
        normalize = self._normalizable(dead_time, settings)
        names: list[str] = ["point"]
        columns: list[np.ndarray] = []
        notes: list[str] = []
        seen: set[str] = set()
        xlabel = TREND_POINT_AXIS
        mean_x = np.empty(0)

        for label, part, _unit in COMBINE_MODES:
            tag = _safe_name(
                label.replace(" - ", "").replace(" only", "")
            ).lower()
            for suffix, these in (("raw", raw_settings), ("norm", settings)):
                if suffix == "norm" and not normalize:
                    continue
                x, curves, xlabel, mode_notes, _used = self._trend_curves(
                    chosen, dead_time, part, separate, these
                )
                if x.size > mean_x.size:
                    mean_x = x
                for note in mode_notes:
                    # The same note comes back from all eight passes.
                    if note not in seen:
                        seen.add(note)
                        notes.append(note)
                for i, y in curves.items():
                    stem = "sum" if i == SUM_KEY else f"el{i + 1}"
                    names.append(f"{tag}_{stem}_{suffix}")
                    columns.append(y)

        for part in ("all", "odd", "even"):
            column = self._monitor_trend(dead_time, part, settings)
            if column is not None:
                names.append(
                    {"all": "mon_sum", "odd": "mon_odd", "even": "mon_even"}[
                        part
                    ]
                )
                columns.append(column)

        points = np.arange(
            max((c.size for c in columns), default=mean_x.size), dtype=float
        )
        columns.insert(0, points)
        if settings["kind"] == "field" and mean_x.size:
            names.insert(1, _safe_name(settings["field"]))
            columns.insert(1, mean_x)
        return names, columns, notes, xlabel

    def _export_trend(self):
        """Write the trend behind Plot ROI to an ASCII table."""
        if not self.plot_roi_btn.isChecked():
            self.status.showMessage(
                "Nothing is plotted -- press Plot ROI first."
            )
            return
        chosen = self._selected_elements()
        if not chosen:
            self.status.showMessage(
                "No element is ticked, so there is nothing to export."
            )
            return
        # Snapshotted before the dialog, which runs a nested event loop: a
        # redraw behind it would otherwise leave the numbers disagreeing
        # with the header written above them.
        settings, axis_note = self._export_settings()
        dead_time = self.deadtime_box.isChecked()
        separate = self.separate_box.isChecked() and len(self.data) == 1

        names, columns, notes, xlabel = self._export_columns(
            chosen, dead_time, separate, settings
        )
        if not columns or max(c.size for c in columns) == 0:
            self.status.showMessage(
                "The trend is empty, so there is nothing to export."
            )
            return

        default = self._export_folder() / (self._export_stem() + ".txt")
        path, _filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export the ROI trend",
            str(default),
            "Text (*.txt);;All files (*)",
        )
        if not path:  # cancelled
            return
        target = Path(path)
        if not target.suffix:
            target = target.with_suffix(".txt")

        header = self._export_header(
            settings, separate, dead_time, axis_note, xlabel
        )
        try:
            write_aligned_columns(target, names, columns, header + notes)
        except OSError as exc:
            self.status.showMessage(
                f"Could not write {target}  -  {exc.strerror or exc}"
            )
            return
        rows = max(c.size for c in columns)
        self.status.showMessage(
            f"Wrote {target}  -  {len(columns)} columns, {rows} rows"
        )

    def _status_text(
        self, first, last, clamped, uncorrected, skipped, notes
    ) -> str:
        parts = [f"Frames {first}-{last}"]
        if len(self.data) > 1:
            parts.append(f"{len(self.data)} scans")
        if clamped:
            parts.append(
                "clamped to their own last frame: " + ", ".join(clamped)
            )
        if uncorrected:
            parts.append("no dead-time factors in: " + ", ".join(uncorrected))
        if skipped:
            parts.append("no selected element in: " + ", ".join(skipped))
        parts.extend(notes)
        return ".  ".join(parts)

    def closeEvent(self, event):
        # A debounced redraw that fires after WA_DeleteOnClose has taken the
        # canvas raises "wrapped C/C++ object has been deleted". The span's
        # canvas callbacks are the same hazard by another route.
        self._redraw_timer.stop()
        if self.roi_span is not None:
            self.roi_span.disconnect_events()
            self.roi_span = None
        super().closeEvent(event)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__.strip().splitlines()[0])
        print("\nusage: spectrum_viewer_sgz.py <file> [...]")
        return 2

    resolved: list[Path] = []
    for name in argv:
        path = resolve_spectrum_file(Path(name))
        if path is None:
            print(f"No spectrum file found for {name}")
            return 1
        if path not in resolved:
            resolved.append(path)

    app = QtWidgets.QApplication(sys.argv[:1])
    # Before the window, so that every widget in it is built against the bold
    # font rather than relaid out afterwards.
    apply_bold_defaults(app)
    window = SpectrumViewerWindow(resolved)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
