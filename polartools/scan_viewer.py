#!/usr/bin/env python
"""The scan browser, with the softglue spectrum window behind its Spc button.

Everything in the browser itself -- the file tree, the scan table, the plots,
the reduce tab -- is :mod:`polartools._scan_viewer_base`, unchanged and
imported. The only difference is which spectrum viewer its "Load Spectrum"
button opens: :mod:`polartools.spectrum_viewer`, the softglue edition, which
adds the Combine box for scans that took more than one vortex frame per point
and offers ``vortex_total_el1`` and ``vortex_total_el7`` as monitors.

The browser itself gains five things, all through the same subclass: bold
text throughout; a blue scan number on any scan that wrote more than one
vortex frame per point -- the sgz softglue layout, and the scans this
edition's Combine box has something to say about; a Find row in the Baseline
window that jumps to a field by part of its name; a directory listing that
appears at once instead of after every file has been read, which on a
2714-scan directory is the difference between 0.05 s and 82 s; and an Auto
refresh button that adds scans to the listing as they land.

Listing a scan whose file is still being written is the one place this
edition deliberately shows something the plain browser does not: that file
is there, so it gets a row with its scan number and blank cells, rather than
being dropped for not opening yet.

This is a rebind and not a fork on purpose. ``_scan_viewer_base`` is 466 KB,
it is what the beamline runs every day, and it looks ``SpectrumViewerWindow``
up as a module global at the moment the button is pressed -- so swapping that
one name is the whole of the difference, and a copy of the file would only
mean every later fix to the browser had to be made twice. Run it the same
way::

    scan-viewer [directory]

Both viewers can be open at once; they share no state.
"""

import contextlib
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from PyQt6 import QtCore
from PyQt6 import QtGui
from PyQt6 import QtWidgets

from . import _scan_viewer_base as _base
from . import spectrum_viewer

# Blue, and nothing else about the row: the scan number is the cell the eye
# goes to down a listing, and a colour says "look here" without taking a
# column away from the eight already in the table. Not the "Spc" tick, which
# answers a different question -- that one says a spectrum exists, this says
# the spectrum has halves.
SOFTGLUE_COLOR = QtGui.QColor("blue")

# How long one chunk of the fill pass may hold the event loop, in
# milliseconds. Time and not a row count, because the cost per row is all
# first-touch I/O on the beamline mount and varies by two orders of magnitude
# between a cold file and a cached one: measured 29 ms per row cold against
# 0.6 ms warm. A deadline is responsive in the first case and barely yields in
# the second.
#
# Small, because the deadline is tested *after* a row and so can only ever
# add rows beyond the first. A row that costs 20 ms overruns any budget below
# 20 ms by the same amount whatever the number is, while a budget large
# enough to fit two such rows simply stalls for both: at 30 ms the watchdog
# measured a 45 ms median stall on a cube-heavy directory, because every
# slice took exactly two 22 ms rows. At 5 ms it takes one, and a warm
# directory still batches eight.
FILL_SLICE_MS = 5

# Reading the first 64 KB of a scan file in a thread, purely to throw the
# bytes away, before the fill pass opens it with h5py.
#
# It looks pointless and it is worth 4.6x. The cost of a row is not the data
# but the open: on GPFS a cold h5py open is ~29 ms of round trips, of which
# 4.55 ms is make_fid alone, and the file is 3.7 MB of which HDF5 reads a few
# kilobytes. Threading h5py itself gains nothing -- it holds a global lock,
# and 8 threads measured slightly *slower* (23.1 -> 24.3 ms/file). A plain
# read() releases the GIL, so sixteen of them overlap their waits and HDF5
# then finds its metadata already in the page cache. Measured cold on 261
# untouched files, three disjoint thirds of one directory:
#
#     no prefetch                    29.39 ms/row
#     64 KB x 16 threads              6.42 ms/row
#     256 KB x 16 threads             6.59 ms/row
#
# 256 KB buying nothing over 64 KB is the confirmation that this is latency
# being overlapped and not bytes being moved, which is why the size below is
# small and the thread count is what matters.
PREFETCH_THREADS = 16
PREFETCH_BYTES = 64 * 1024
# How many reads to keep in flight ahead of the fill. Bounded rather than
# submitting every row at once so that cancelling is just "stop submitting",
# and so that scrolling to the far end of a long listing does not leave two
# thousand queued reads of rows nobody is looking at any more.
PREFETCH_AHEAD = 64
# When to conclude the prefetch is not earning its place, and stop: once this
# many rows are filled, if they are averaging less than this long each.
#
# The test is the speed of the *fill*, not of the reads. A Refresh of a
# directory already listed is the case that matters: read_scan_info and the
# two resolvers all memoise on (path, mtime, size), so every row is a dict
# lookup and the whole pass takes 0.1 s -- but nothing has re-read the files,
# so the page cache may well have lost them and a prefetch read still
# measures a cold 23 ms. Timing the reads therefore concludes "cold, keep
# going" and queues hundreds of reads of bytes no one will look at. Timing
# the rows sees 0.05 ms a row and stops, which is the right answer.
#
# Measured row costs are three orders of magnitude apart either side of this
# line: 22 ms cold on a cube-heavy directory, 5 ms cold on a plain one,
# 0.05 ms on a Refresh.
PREFETCH_PROBE = 32
PREFETCH_IDLE_S = 0.001

# How often the Auto refresh button looks for new scans, in milliseconds.
#
# Set by how soon someone wants to see a scan appear, not by what the looking
# costs: one :func:`_scandir_hdf_files` is 0.09-0.16 ms on the 82-file live
# directory and 1.8 ms on the 2714-file one, so even the big directory spends
# four hundredths of a percent of a core on this.
#
# What the number cannot promise is set by the mount, and the two here differ.
# ``/gdata`` is GPFS, whose metadata is coherent across nodes, so a file
# written by the acquisition host really is visible within the interval.
# ``/home/beams/POLAR/data`` is NFSv4 mounted with no ``actimeo``, where the
# kernel serves a cached directory listing for ``acdirmin``..``acdirmax`` --
# 30 to 60 s by default -- so detection there takes ~30 s whatever this says.
# That is also why this polls instead of using QFileSystemWatcher: inotify is
# told about writes by the local kernel and never hears about another node's
# at all, so the watcher would sit silent while the poll did all the work.
REFRESH_INTERVAL_MS = 5000


def _scandir_hdf_files(directory: Path) -> list[Path]:
    """The HDF5 files directly inside *directory*, sorted by name.

    Replaces :func:`_base.list_hdf_files`, which asks ``Path.is_file()``
    of every entry and so pays a stat syscall each. ``os.scandir`` gets the
    same answer out of the directory block it had to read anyway: on the
    2714-file Kb_comm directory, 0.832 s against 0.011 s, with the two lists
    comparing equal.

    ``os.path.splitext`` rather than ``Path.suffix`` only to avoid building
    2714 Path objects to ask one question of each. The two agree everywhere,
    including on a name that is all suffix -- ".hdf" has no suffix to either.

    Errors are left to propagate, as the base version leaves them: a
    directory that has gone away should say so rather than list as empty.
    """
    with os.scandir(directory) as entries:
        return sorted(
            Path(entry.path)
            for entry in entries
            if os.path.splitext(entry.name)[1].lower() in _base.HDF_SUFFIXES
            and entry.is_file()
        )


def _warm(path: Path):
    """Pull the head of *path* into the page cache; discard what comes back.

    Runs in a prefetch thread, and deliberately touches nothing else -- no
    h5py, no Qt, no shared state. It cannot change what any row ends up
    saying, only how fast the fill pass gets to say it, which is what makes
    threading safe here at all. Failures are silent because there is nothing
    to report: the fill pass opens the same file a moment later and will
    raise, or not, on its own account.

    Only the master, although the fill opens the vortex cube too and the cube
    is the dearer of the two. Warming it was tried twice and lost twice:

      * in this pool, 1.6x slower -- sixteen threads reading two files each
        saturate the mount and the fill's own opens queue behind them. On a
        702-scan directory, 15.9 s to 26 s and the 99th-percentile stall from
        28 ms to 290 ms.
      * in a separate 4-thread pool eight rows ahead, so as not to take
        threads from the masters: still 1.1x slower, 38.4 ms a row against
        42.7 ms, measured on disjoint cold halves of one 704-scan directory.

    Which is worth knowing because it is the opposite of what the same trick
    does for masters (4.6x faster). The likely reason is that 64 KB of the
    head is the wrong 64 KB of a cube: the masters are 3.7 MB and HDF5 finds
    its metadata in the part that was warmed, while a cube is far larger and
    what ``is_multi_frame_scan`` reads need not be near the front.
    """
    try:
        with open(path, "rb") as handle:
            handle.read(PREFETCH_BYTES)
    except OSError:
        pass


def _cheap_scan_info(_path: Path) -> dict:
    """Stand in for ``read_scan_info`` during the cheap first pass.

    The listing comprehension in the base ``_reload_file_list`` asks
    ``read_scan_info(p)["is_scan"]`` of every file, which is one HDF5 open per
    file before a single row is drawn. This says yes without looking.

    Optimistic, and checked: across five data cycles no directory held a
    top-level HDF file that was not a scan. One that did would be listed and
    then hidden a moment later by the fill pass, which asks the real question.

    Only ``is_scan`` is in the dict on purpose. Anything else the base class
    grows a use for here will raise KeyError in testing rather than quietly
    read as a default.
    """
    return {"is_scan": True}


@contextlib.contextmanager
def _module_global(module, name: str, value):
    """Swap one of *module*'s globals for the body of the block, then put it back."""
    original = getattr(module, name)
    setattr(module, name, value)
    try:
        yield
    finally:
        setattr(module, name, original)


class SgHdfViewer(_base.HdfViewer):
    """The browser, in bold, with softglue scan numbers in blue.

    A subclass and not an edit to ``_scan_viewer_base.py``: that file is what the
    beamline runs every day and stays byte-for-byte identical. Both additions
    are things the base class has no hook for, so each overrides the smallest
    method that can carry it.
    """

    def __init__(self, *args, **kwargs):
        # Set before super(), which calls _reload_file_list before this
        # constructor's own body would have run.
        self._fill_timer = None
        self._fill_paths: list[Path | None] = []
        self._fill_done = bytearray()
        self._fill_row = 0
        self._filled = 0
        # Rows whose file could not be read as a scan. They stay on screen
        # showing their scan number -- see _fill_row_now -- and the auto
        # refresh re-reads them, so this is a worklist as well as a count.
        self._bare: set[int] = set()
        self._loading_cheap = False
        self._prefetch_pool: ThreadPoolExecutor | None = None
        self._prefetch_sent = bytearray()
        self._prefetch_futures: list = []
        self._fill_seconds = 0.0
        # The auto refresh, likewise set before super() because
        # _reload_file_list writes _listed_names. The QTimer itself cannot be
        # built here -- QTimer(self) needs a QObject, and self is not one
        # until super().__init__ has run -- so it is made on first use, the
        # same shape as _fill_timer above.
        self._refresh_timer = None
        self._listed_names: set[str] = set()
        self._suppress_replot = False
        # Also before super(), so every widget the browser builds is laid out
        # against the bold font rather than resized under it afterwards.
        # _base.main() builds the QApplication and calls exec_() itself,
        # with no seam between the two to set a font in -- this constructor is
        # that seam, and it keeps main() as the one place the startup
        # sequence lives.
        spectrum_viewer.apply_bold_defaults(QtWidgets.QApplication.instance())
        super().__init__(*args, **kwargs)
        self._add_auto_refresh_button()

    # -- listing in two passes ---------------------------------------------

    def _reload_file_list(self):
        """List the directory at once, then read the files in the background.

        The base class reads every file before it draws anything: three HDF5
        opens a row, 71% of it waiting on ``make_fid``, all of it before the
        window exists. On the 2714-scan Kb_comm directory that is 82 s during
        which there is nothing on screen to look at and no event loop to
        answer a click.

        So the work is split in two. The cheap pass is the base method with
        the expensive questions stubbed out -- it draws every row with its
        scan number, which comes free from the filename, and leaves the other
        seven cells blank: 0.042 s for the same 2714 rows. The fill pass then
        answers those cells a slice at a time under a live event loop.

        Blank rather than a placeholder dash: a dash is what the Img and Spc
        columns say when a file genuinely has no image, so using it for "not
        read yet" would be the table stating something false for twelve
        seconds. Nothing is the one thing that is not a wrong answer.
        """
        self._stop_fill()
        self._loading_cheap = True
        try:
            with _module_global(_base, "read_scan_info", _cheap_scan_info):
                super()._reload_file_list()
        finally:
            self._loading_cheap = False
        if self._directory is None:
            return
        self._start_fill()
        if self._suppress_replot:
            # An auto refresh reload, where this would cost the user their
            # zoom: _update_plot clears every axis it draws on, so replotting
            # redraws the same curves at the default limits. The selection is
            # provably the one already on screen -- the probe reloads only
            # when a name has been *added*, and _select_rows restores by name
            # -- so what is drawn is still right, and these widgets are stale
            # only against changes inside the selected file. A master cannot
            # have any: it is written once, by NXWriter's thread after the
            # run stops, and _setup_file_io refuses to overwrite one.
            #
            # Do not "simplify" this to the argument that an unchanged
            # selection has nothing new to say. That is false -- the hkl
            # label, the stream combo and the baseline window all rekey on
            # (path, mtime, size) and would happily return new content -- and
            # believing it is how the flag gets deleted and the zoom starts
            # vanishing every five seconds.
            return
        # Suppressed during the cheap pass, so run it once now that the real
        # read_scan_info is back -- otherwise a selection restored across a
        # Refresh would plot from the stub's answer and never be corrected.
        self._on_selection_changed()

    def _scan_row(self, path: Path) -> list[str]:
        """The eight cells of one row; during the cheap pass, only the first.

        An instance method overriding a ``@staticmethod``, which works because
        the base calls it as ``self._scan_row(path)``.
        """
        if not self._loading_cheap:
            return super()._scan_row(path)
        return [_base.scan_number(path)] + [""] * (len(_base.FILE_COLUMNS) - 1)

    def _on_selection_changed(self):
        """Nothing during the cheap pass; :meth:`_reload_file_list` calls it after.

        The stub ``read_scan_info`` is in force while the cheap pass runs, and
        this is the one thing the base method does afterwards that would
        believe it.
        """
        if self._loading_cheap:
            return
        super()._on_selection_changed()

    # -- filling the rows in the background --------------------------------

    def _start_fill(self):
        """Begin the pass that reads the files and fills in what they say.

        Restarting from scratch on every reload is deliberate: the base class
        builds fresh items each time, so there is nothing of a previous pass
        to resume, and a directory change must not keep filling rows that have
        since become other scans.
        """
        table = self.file_table
        self._fill_paths = []
        for row in range(table.rowCount()):
            item = table.item(row, 0)
            data = (
                item.data(QtCore.Qt.ItemDataRole.UserRole)
                if item is not None
                else None
            )
            self._fill_paths.append(Path(data) if data else None)
        self._fill_done = bytearray(len(self._fill_paths))
        self._prefetch_sent = bytearray(len(self._fill_paths))
        self._fill_row = 0
        self._filled = 0
        self._bare = set()
        self._fill_seconds = 0.0
        # Every file that now has a row, readable or not, which is what the
        # auto refresh diffs against to decide something has arrived. Set
        # before the empty-listing return below, so that an empty directory
        # leaves an empty set rather than the previous directory's names --
        # watching a fresh proposal directory fill up is the likeliest use of
        # the refresh, and it would otherwise diff against the wrong thing.
        self._listed_names = {p.name for p in self._fill_paths if p is not None}
        if not self._fill_paths:
            return
        self._prefetch_pool = ThreadPoolExecutor(
            PREFETCH_THREADS, thread_name_prefix="sgz-prefetch"
        )
        if self._fill_timer is None:
            self._fill_timer = QtCore.QTimer(self)
            self._fill_timer.setInterval(0)
            self._fill_timer.timeout.connect(self._fill_next_slice)
        self._fill_timer.start()

    def _stop_fill(self):
        """Halt the fill and drop the prefetch pool.

        ``cancel_futures`` matters on the way out: the pool's threads are not
        daemons, so a queue of reads left in it would hold the interpreter
        open after the window has closed.
        """
        if self._fill_timer is not None:
            self._fill_timer.stop()
        self._drop_prefetch_pool()

    def _drop_prefetch_pool(self):
        """Shut the prefetch pool down, if there is one, and forget its work."""
        pool, self._prefetch_pool = self._prefetch_pool, None
        if pool is not None:
            pool.shutdown(wait=False, cancel_futures=True)
        self._prefetch_futures = []

    def closeEvent(self, event):
        """Stop the background work before the base class remembers and closes."""
        self._stop_fill()
        if self._refresh_timer is not None:
            self._refresh_timer.stop()
        super().closeEvent(event)

    def _visible_rows(self) -> tuple[int, int]:
        """The span of rows on screen, as inclusive indices."""
        table = self.file_table
        first = table.rowAt(0)
        if first < 0:
            return 0, -1
        last = table.rowAt(table.viewport().height() - 1)
        if last < 0:
            last = table.rowCount() - 1
        return first, last

    def _next_row(self) -> int | None:
        """The next row to read: what is on screen first, then straight down.

        Re-asked every slice rather than fixed at the start, so that scrolling
        to the bottom of a long listing -- where the newest scans are, and so
        where someone opening this directory is usually headed -- redirects
        the reading to the rows actually being looked at instead of making
        them wait out the two thousand above.
        """
        first, last = self._visible_rows()
        for row in range(
            max(first, 0), min(last, len(self._fill_paths) - 1) + 1
        ):
            if not self._fill_done[row]:
                return row
        while (
            self._fill_row < len(self._fill_paths)
            and self._fill_done[self._fill_row]
        ):
            self._fill_row += 1
        if self._fill_row < len(self._fill_paths):
            return self._fill_row
        return None

    def _prefetch_from(self, start: int):
        """Keep reads in flight for the rows about to be read.

        Scans forward from *start* for rows not yet sent, and gives up after a
        bounded look so that a nearly-finished listing does not walk its whole
        length every slice looking for work.
        """
        pool = self._prefetch_pool
        if pool is None:
            return
        if (
            self._filled >= PREFETCH_PROBE
            and self._fill_seconds / self._filled < PREFETCH_IDLE_S
        ):
            # The rows are coming back from memo caches. Nothing to overlap,
            # so drop the pool for the rest of this directory rather than
            # read a thousand files that no one is waiting on.
            self._drop_prefetch_pool()
            return
        self._prefetch_futures = [
            f for f in self._prefetch_futures if not f.done()
        ]
        room = PREFETCH_AHEAD - len(self._prefetch_futures)
        for row in range(
            start, min(start + 4 * PREFETCH_AHEAD, len(self._fill_paths))
        ):
            if room <= 0:
                break
            if self._prefetch_sent[row] or self._fill_done[row]:
                continue
            path = self._fill_paths[row]
            if path is None:
                continue
            self._prefetch_sent[row] = 1
            self._prefetch_futures.append(pool.submit(_warm, path))
            room -= 1

    def _fill_next_slice(self):
        """Read rows until the deadline, then give the event loop its turn."""
        deadline = QtCore.QElapsedTimer()
        deadline.start()
        table = self.file_table
        row = self._next_row()
        if row is None:
            self._finish_fill()
            return
        self._prefetch_from(row)
        started = time.perf_counter()
        while row is not None:
            self._fill_row_now(table, row)
            self._fill_done[row] = 1
            self._filled += 1
            if deadline.hasExpired(FILL_SLICE_MS):
                break
            row = self._next_row()
        # Only the reading, not the prefetch bookkeeping above it: this is
        # the number the prefetch gate is asking about.
        self._fill_seconds += time.perf_counter() - started
        self.status.showMessage(
            f"{self._scan_count()} scan(s) in {self._directory}. "
            f"Reading {self._filled}/{len(self._fill_paths)}…"
        )

    def _fill_row_now(self, table, row: int):
        """Read one file and put what it says into its row.

        One row is one slice's worth of work and cannot be broken up further,
        because nothing yields inside an HDF5 open. On a directory of vortex
        cubes that open is the whole cost and the floor on how long a slice
        holds the event loop: cold, with the master prefetched, 3.0 ms median
        for ``read_scan_info`` against 37.8 ms for ``is_multi_frame_scan``.

        Testing the deadline between those two steps instead -- so that a
        slice could end between the master and the cube -- was tried, and
        moves the median stall only from 41.0 ms to 37.9 ms, since prefetch
        has already made the master half small. Not worth carrying a
        half-filled row for.
        """
        path = self._fill_paths[row]
        if path is None:
            return
        if not _base.read_scan_info(path)["is_scan"]:
            # Almost always a master still being written, and left on screen
            # on purpose: its scan number from the cheap pass, the other
            # seven cells blank. A file that exists and is not listed is the
            # more surprising of the two, and at the bottom of a filled table
            # a bare row reads as "that one has just finished".
            #
            # The window is real. The master does not exist during the scan
            # at all -- _setup_file_io raises FileExistsError if it does --
            # and is written afterwards by NXWriter's thread, which waits for
            # the area detector's own files before it finishes. Until it
            # does, h5py cannot open what is there: measured on prefixes of a
            # finished master, every one from 1% to 99% reads as is_scan
            # False and only the whole file reads as a scan.
            #
            # Which is also why the row can be filled once and left: there is
            # no half-read state to correct. The Range column would be the
            # travel *so far* if a growing file could be opened -- it is
            # values[0] to values[-1] of the primary stream, not the planned
            # range -- but one never can, so the cells this row eventually
            # gets are the finished scan's.
            self._bare.add(row)
            return
        self._bare.discard(row)
        cells = self._scan_row(path)
        # From 1: column 0 is the scan number, which the cheap pass already
        # wrote, and rewriting it would undo the blue below on a refresh.
        for column in range(1, len(cells)):
            item = table.item(row, column)
            if item is not None:
                item.setText(cells[column])
        self._color_row_now(table, row, path)

    def _color_row_now(self, table, row: int, path: Path):
        """Blue one row's scan number if its cube holds halves.

        Folded into the fill rather than a pass of its own, as it used to be:
        it needs the same ``resolve_spectrum_file`` the "Spc" cell just asked
        for, and one pass over the listing is one set of opens instead of two.
        Both that and :func:`spectrum_viewer.is_multi_frame_scan` cache on
        (path, mtime, size), so the second ask costs a dict lookup and a
        refresh of an unchanged directory reopens nothing at all.
        """
        cube = spectrum_viewer.resolve_spectrum_file(path)
        if cube is None or not spectrum_viewer.is_multi_frame_scan(cube):
            return
        item = table.item(row, 0)
        if item is None:
            return
        item.setForeground(QtGui.QBrush(SOFTGLUE_COLOR))
        item.setToolTip(
            "More than one vortex frame per scan point: a potential sgz\n"
            "softglue scan. Load its spectrum to split the frames into\n"
            "odd and even with the Combine box."
        )

    def _scan_count(self) -> int:
        """Rows that turned out to be scans -- the listing less the bare ones.

        Counts the rows not yet read as scans, so during the fill it starts
        at the whole listing and settles downwards. That matches the
        "Reading n/m" it sits beside, which is also about the pass so far.
        """
        return len(self._fill_paths) - len(self._bare)

    def _finish_fill(self):
        """Drop the machinery and leave the status bar as the base class has it."""
        self._stop_fill()
        self.status.showMessage(
            f"{self._scan_count()} scan(s) in {self._directory}."
        )

    # -- watching the directory for new scans ------------------------------

    def _add_auto_refresh_button(self):
        """Put an Auto refresh toggle in the toolbar, after Refresh.

        A checkable QPushButton because that is the house idiom for a toggle
        here -- ``grid_btn`` in the base class is the same thing -- and there
        is no QCheckBox anywhere in ``_scan_viewer_base.py`` to match instead.

        The row is found by asking which layout holds a button we know rather
        than by position, so that a toolbar edit upstream cannot quietly drop
        this somewhere odd; and if the toolbar is ever rebuilt in a way that
        breaks the search, it says so here instead of raising StopIteration
        out of a constructor three frames up.
        """
        row = next(
            (
                layout
                for layout in self.findChildren(QtWidgets.QHBoxLayout)
                if layout.indexOf(self.export_btn) != -1
            ),
            None,
        )
        if row is None:
            raise RuntimeError(
                "_scan_viewer_base._build_ui no longer puts export_btn in a "
                "QHBoxLayout; the Auto refresh button has nowhere to go."
            )

        self.auto_btn = QtWidgets.QPushButton("Auto refresh")
        self.auto_btn.setCheckable(True)
        self.auto_btn.setToolTip(
            f"Look for new scans every {REFRESH_INTERVAL_MS // 1000} s and add\n"
            "them to the list, scrolling down to follow them if the list is\n"
            "already at the bottom. Nothing is selected or plotted.\n"
            "A scan whose file is still being written is listed with its\n"
            "number alone; the rest of its row fills in when the file is\n"
            "complete."
        )
        self.auto_btn.toggled.connect(self._toggle_auto_refresh)
        row.insertWidget(row.indexOf(self.image_btn), self.auto_btn)

    def _toggle_auto_refresh(self, on: bool):
        """Start or stop the watch, and look once immediately on starting.

        Looking straight away matters because the button is usually pressed
        *after* something has landed -- the list looks out of date, so you
        reach for it -- and waiting out a full interval first would make it
        seem not to have worked.
        """
        if not on:
            if self._refresh_timer is not None:
                self._refresh_timer.stop()
            return
        if self._refresh_timer is None:
            self._refresh_timer = QtCore.QTimer(self)
            self._refresh_timer.setInterval(REFRESH_INTERVAL_MS)
            self._refresh_timer.timeout.connect(self._probe_directory)
        self._refresh_timer.start()
        self._probe_directory()

    def _probe_directory(self):
        """Look for scans that have arrived, and for bare rows that can now be read.

        Everything here runs unattended while the user is in the hutch, and
        an exception out of a Qt timer slot does not reach a traceback: PyQt5
        calls ``abort()``, and nothing in these five modules installs an
        excepthook. So each way out is a guard rather than something to be
        caught afterwards, and the body that is left runs inside one
        ``except OSError``.
        """
        if self._loading_cheap:
            # read_scan_info is the stub inside that window, and it holds
            # nothing but "is_scan".
            return
        if self._directory is None:
            # os.scandir(None) lists the *current working directory* instead
            # of raising, which would diff against our names, find them all
            # new, and reload every five seconds for the life of the window.
            return
        if self._fill_timer is not None and self._fill_timer.isActive():
            # A fill is still running. Reloading now would throw it away and
            # start again -- a cold 2714-scan directory takes ~12 s, so at
            # this interval it would never finish. The "is not None" is load
            # bearing: an empty directory leaves _fill_timer unbuilt, and
            # waiting on an empty directory is exactly what this is for.
            return
        if QtWidgets.QApplication.activeModalWidget() is not None:
            # Not correctness -- a dialog's nested event loop runs timers
            # quite happily -- but a 12 s fill underneath an open dialog is
            # not what anyone wants. Misses native and portal file dialogs,
            # which Qt does not own.
            return

        try:
            names = {path.name for path in _scandir_hdf_files(self._directory)}
            if names - self._listed_names:
                self._reload_for_new_scans()
                return
            self._refill_bare_rows()
        except OSError as exc:
            # _scandir_hdf_files lets a vanished directory raise, which is
            # right from a button press and fatal from a timer. Stop rather
            # than fail silently every five seconds, and untick the button so
            # that the state on screen is the truth.
            self.auto_btn.setChecked(False)
            self.status.showMessage(f"Auto refresh stopped: {exc}")

    def _reload_for_new_scans(self):
        """Rebuild the listing for an arrived scan, following it if we were at the end.

        Additions only. The user asked to be shown new scans, and looking at
        nothing else is also what makes the suppressed replot safe by
        construction rather than by a judgement that has to stay right: no
        listed name has gone away, ``_select_rows`` restores by name, so the
        selection after the rebuild is necessarily the one before it.
        A deleted file therefore leaves a stale row until the next manual
        Refresh, which is a far smaller surprise than a plot disappearing on
        its own.
        """
        bar = self.file_table.verticalScrollBar()
        # Before the rebuild: straight after one, maximum() is still the old
        # value -- measured 44 where 45 was right -- so asking afterwards
        # reads the wrong number. ">=" rather than "==" because a table too
        # short to scroll has value and maximum both 0, and following the
        # newest scan is the right default for a directory just opened.
        at_bottom = bar.value() >= bar.maximum()

        self._suppress_replot = True
        try:
            self._reload_file_list()
        finally:
            self._suppress_replot = False

        if at_bottom:
            # Synchronous: this executes the pending item layout, so the fill
            # pass's first slice already sees the bottom rows as the visible
            # ones and reads the new scan before the rest of the listing.
            self.file_table.scrollToBottom()

    def _refill_bare_rows(self):
        """Re-read the rows whose file could not be read, in place.

        No reload, so no rebuild, no flicker, and nothing done to the
        selection or the scroll position -- the row simply stops being blank.
        Normally there are none, or the one scan that has just finished.

        Cheap because ``read_scan_info`` caches its *failures* too, on
        (path, mtime, size): a file that is not changing is a dict lookup
        behind one stat, so a permanently unreadable stray costs nothing
        after the first look. Only a file that is growing re-keys, which is
        the master being written, and that is one real open per interval
        until it is done.

        A bare row that is selected keeps its stale Img/Spc buttons until it
        is clicked again. Refreshing them means _on_selection_changed, which
        replots every selected scan and so would take the zoom off the others
        in a multi-row selection -- a worse trade than a button that is right
        one click later.
        """
        table = self.file_table
        was = len(self._bare)
        for row in sorted(self._bare):
            self._fill_row_now(table, row)
        if len(self._bare) != was:
            # One of them turned out to be a scan after all, so the count in
            # the status bar is now one short of what the table shows.
            self.status.showMessage(
                f"{self._scan_count()} scan(s) in {self._directory}."
            )

    # -- finding a field in the baseline table -----------------------------

    def _build_baseline_panel(self):
        """The base panel, with a Find row added under the file combo.

        An insert into the built panel rather than a rewritten
        ``_build_baseline_panel``: the base method is 20 lines that would
        then exist twice, and every column or tooltip added to the baseline
        table upstream would have to be added here too.

        Index 1 puts the row between the file combo and the table, which is
        where the eye goes after choosing a file and before reading down the
        rows -- and keeps it out of the way of the table's own stretch.
        """
        panel = super()._build_baseline_panel()
        panel.layout().insertLayout(1, self._build_baseline_find_row(panel))
        return panel

    def _build_baseline_find_row(self, panel):
        """The Find row: a box, two step buttons and a match count."""
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Find field:"))

        self.baseline_find = QtWidgets.QLineEdit()
        self.baseline_find.setPlaceholderText(
            'part of a field name, e.g. "tth" or "temp"'
        )
        self.baseline_find.setClearButtonEnabled(True)
        # On every keystroke, not on Enter alone. A baseline runs to a few
        # hundred fields, which is well past what anyone pages through by
        # eye, so the box is the way to the row -- and one that answers only
        # once the whole name is typed reads as broken while you type it.
        self.baseline_find.textChanged.connect(self._baseline_find_here)
        self.baseline_find.returnPressed.connect(self._baseline_find_next)
        row.addWidget(self.baseline_find, 1)

        for text, tip, slot in (
            ("▲", "Previous match (Shift+Enter)", self._baseline_find_prev),
            ("▼", "Next match (Enter)", self._baseline_find_next),
        ):
            button = QtWidgets.QToolButton()
            button.setText(text)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            row.addWidget(button)

        self.baseline_find_label = QtWidgets.QLabel("")
        # Wide enough for "no match" and for a three-digit "nnn / nnn", so
        # that stepping through the matches does not shuffle the buttons
        # sideways under the pointer.
        self.baseline_find_label.setMinimumWidth(90)
        row.addWidget(self.baseline_find_label)

        # returnPressed cannot tell the two apart, so the backwards step
        # needs a shortcut of its own.
        back = QtGui.QShortcut(
            QtGui.QKeySequence("Shift+Return"), self.baseline_find
        )
        back.activated.connect(self._baseline_find_prev)

        find = QtGui.QShortcut(QtGui.QKeySequence.StandardKey.Find, panel)
        find.setContext(QtCore.Qt.ShortcutContext.WidgetWithChildrenShortcut)
        find.activated.connect(self._focus_baseline_find)
        return row

    def _focus_baseline_find(self):
        """Put the cursor in the Find box and select what is already there."""
        self.baseline_find.setFocus(QtCore.Qt.FocusReason.ShortcutFocusReason)
        self.baseline_find.selectAll()

    def _baseline_matches(self):
        """Rows whose Field cell contains the query, in the order displayed.

        Walked afresh on every search rather than kept: the table sorts on a
        header click, which reorders the rows behind our back and emits
        nothing this subclass is subscribed to, so a remembered list would
        quietly start pointing at the wrong fields. A few hundred rows of
        ``in`` costs well under a millisecond, which is nothing against the
        keystroke that asked for it.
        """
        text = self.baseline_find.text().strip().lower()
        if not text:
            return []
        table = self.baseline_table
        return [
            row
            for row in range(table.rowCount())
            if (item := table.item(row, 0)) is not None
            and text in item.text().lower()
        ]

    def _baseline_find_here(self):
        """Show the match for a changed query, keeping the row if it still fits."""
        self._baseline_step(None)

    def _baseline_find_next(self):
        """Jump to the next match below the current row, wrapping at the end."""
        self._baseline_step(1)

    def _baseline_find_prev(self):
        """Jump to the previous match above the current row, wrapping at the top."""
        self._baseline_step(-1)

    def _baseline_step(self, step):
        """Move to a matching row: *step* None re-answers a changed query."""
        rows = self._baseline_matches()
        if not rows:
            # An empty box is not a failed search -- it is the state the
            # window opens in, and saying "no match" there would be a
            # complaint about not having been asked anything yet.
            self.baseline_find_label.setText(
                "no match" if self.baseline_find.text().strip() else ""
            )
            return

        current = self.baseline_table.currentRow()
        if step is None:
            # Narrowing a query that still describes the row in front of you
            # should not throw the view back to the top: that row is what you
            # were typing towards.
            target = current if current in rows else rows[0]
        elif step > 0:
            target = next((row for row in rows if row > current), rows[0])
        else:
            target = next(
                (row for row in reversed(rows) if row < current), rows[-1]
            )

        self.baseline_table.setCurrentCell(target, 0)
        self.baseline_table.selectRow(target)
        self.baseline_table.scrollToItem(
            self.baseline_table.item(target, 0),
            QtWidgets.QAbstractItemView.ScrollHint.PositionAtCenter,
        )
        self.baseline_find_label.setText(
            f"{rows.index(target) + 1} / {len(rows)}"
        )

    def _update_baseline(self):
        """Refill the table as the base class does, then re-answer the query.

        Without this the count would still describe the file just switched
        away from, and the selected row would be whatever sat at that index
        in a different scan's baseline.
        """
        super()._update_baseline()
        # Guarded because the base constructor connects this to the file
        # combo from inside _build_baseline_panel, which is to say before the
        # Find row that super() call returns into existence.
        if getattr(self, "baseline_find", None) is not None:
            self._baseline_find_here()


# Rebound before main() builds anything, and by name rather than by editing
# the base's import line, so that _scan_viewer_base.py stays byte-for-byte
# the file upstream maintains. HdfViewer has to go this way because main()
# looks it up as a module global at the moment it builds the window.
#
# The first two are no-ops here and kept only so this block still says what
# the browser needs: packaging dropped the plain spectrum viewer, so the base
# already imports both names from the softglue module these assign. Upstream
# they did real work. resolve_spectrum_file goes with the window either way --
# it is how the browser decides a scan has a spectrum at all, and the two have
# to agree about which file that is.
_base.SpectrumViewerWindow = spectrum_viewer.SpectrumViewerWindow
_base.resolve_spectrum_file = spectrum_viewer.resolve_spectrum_file
_base.HdfViewer = SgHdfViewer
# Same output, 75x less time, and the one caller is the listing this edition
# is trying to make fast. Rebound rather than subclassed because it is a
# module-level function the base _reload_file_list looks up as a global.
_base.list_hdf_files = _scandir_hdf_files


def main():
    """Console entry point for the ``scan-viewer`` application."""
    # The rebinds above have already run, at import; the base's main() builds
    # the QApplication, reads sys.argv[1] as the initial directory and does
    # not return.
    return _base.main()


if __name__ == "__main__":
    sys.exit(main())
