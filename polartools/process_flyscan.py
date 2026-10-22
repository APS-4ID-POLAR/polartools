from dataclasses import dataclass
import os
import re

from h5py import File
import hdf5plugin  # noqa: F401  -- registers the HDF5 codecs these files use
import numpy as np
import matplotlib.pyplot as plt

FNAME_FORMAT = "{}/{}/scan_{:06d}.h5"

# The eigerTrig divider is DivByN with N = 5 x N_ckIM and shares its reset
# with the ckIM divider, so it emits its first pulse only after 5 ckIM ticks:
# ckIM 1-5 carry eigerTrig 0 and belong to no exposure, and detector frame k
# is exposed during the ckIM ticks carrying eigerTrig value k + FIRST_TRIGGER.
# Both counters are recorded in the position stream (columns 1 and 2), and
# this holds absolutely -- the counters are zeroed together at each scan.
FIRST_TRIGGER = 1

# ckIM ticks per eigerTrig pulse, i.e. position samples per exposure. Only a
# fallback: the value actually used is measured per scan by `full_group_size`.
CKIM_PER_TRIGGER = 5


@dataclass
class PositionGroups:
    """Per-trigger reduction of a position stream.

    Every array has one entry per trigger group, in file order. ``trigs``
    holds each group's trigger value, which is what associates a group with
    a detector frame (see :func:`groups_for_frames`) — the association is
    made by the hardware trigger counter, not by array position, so a
    missing group drops one frame instead of shifting all the later ones.
    """

    i0s: np.ndarray
    xs: np.ndarray
    ys: np.ndarray
    dxs: np.ndarray
    dys: np.ndarray
    trigs: np.ndarray  # eigerTrig value of each group
    counts: np.ndarray  # rows actually recorded in each group
    coverage: np.ndarray  # ckIM ticks each group spans
    partial: np.ndarray  # bool: group spans less than a full exposure


@dataclass
class PartialGroup:
    """One flagged trigger group, and how it differs from a full exposure."""

    group: int  # index into the PositionGroups / TickResult arrays
    trig: int  # its eigerTrig value
    frame: int  # detector frame it holds; -1 if none maps to it
    coverage: int  # ckIM ticks it spans — how much of the exposure
    expected: int  # ckIM ticks in a full exposure
    count: int  # rows actually recorded (< coverage if any dropped)
    x: float
    y: float
    span: float  # distance the motor covered during the group
    full_span: float  # typical span of a full group, for comparison
    last: bool = False  # the scan's final group, where the detector may
    # itself have been disarmed mid-exposure

    @property
    def dropped(self):
        """Ticks inside this group that the MCS never recorded."""
        return self.coverage - self.count


def full_group_size(coverage):
    """ckIM ticks in a full exposure: the modal value, ties resolved upwards.

    Group 0 spans whatever lies between the two dividers' resets — normally 5
    ticks, but 1, 3 and 37 all occur — so it is not allowed to define what
    "full" means.
    """
    coverage = np.asarray(coverage)
    if coverage.size == 0:
        return 0
    body = coverage[1:] if coverage.size > 1 else coverage
    vals, freq = np.unique(body, return_counts=True)
    return int(vals[freq == freq.max()].max())


def group_coverage(ckim, starts, n):
    """ckIM ticks spanned by each trigger group.

    Measured from the ckIM counter, not from the number of rows, so that a
    sample the MCS failed to record does not look like a short exposure: a
    group whose middle row is missing still spans its full five ticks. Only a
    group that really covers less of the exposure — because the stream
    started or stopped inside it, or lost its leading tick — comes out short.
    """
    ends = np.concatenate((starts[1:], [n]))
    return (np.asarray(ckim)[ends - 1] - np.asarray(ckim)[starts] + 1).astype(
        np.int64
    )


def flag_partial_groups(coverage):
    """Mark groups whose position samples don't span their whole exposure.

    This is a statement about the position record, never about the detector:
    the Eiger is triggered independently and integrates for its full frame
    time regardless. What these groups lack is position data covering all of
    that time, so their mean x/y and I0 describe only part of the travel.

    They are the pre-trigger group, a final group cut off when the position
    stream stopped mid-exposure, and any group that lost its leading ckIM
    tick. A group that merely dropped an interior sample still spans its
    whole exposure and is not flagged — see :func:`group_coverage`.

    Called live as well as post-scan, so with only a couple of closed groups
    the verdict can change as more data arrives.
    """
    coverage = np.asarray(coverage)
    if coverage.size == 0:
        return np.zeros(0, dtype=bool)
    return coverage < full_group_size(coverage)


def describe_partial_groups(groups, n_frames=0):
    """One :class:`PartialGroup` per flagged group, for display.

    ``groups`` is anything carrying the per-group arrays ``xs``, ``ys``,
    ``dxs``, ``dys``, ``trigs``, ``counts`` and ``partial`` — a
    :class:`PositionGroups` or a live ``TickResult``. ``n_frames`` is how
    many detector frames exist, used to say which frame each flagged group
    holds (``-1`` when the group is not plotted at all, as for the
    pre-trigger group).
    """
    partial = np.asarray(groups.partial, dtype=bool)
    if partial.size == 0 or not partial.any():
        return []
    counts = np.asarray(groups.counts)
    coverage = np.asarray(groups.coverage)
    expected = full_group_size(coverage)
    spans = np.hypot(np.asarray(groups.dxs), np.asarray(groups.dys))
    full = spans[~partial]
    full_span = float(np.median(full)) if full.size else float("nan")

    frames, gidx = groups_for_frames(groups.trigs, n_frames)
    frame_of = dict(zip(gidx.tolist(), frames.tolist()))

    return [
        PartialGroup(
            group=int(g),
            trig=int(groups.trigs[g]),
            frame=int(frame_of.get(int(g), -1)),
            coverage=int(coverage[g]),
            expected=expected,
            count=int(counts[g]),
            x=float(groups.xs[g]),
            y=float(groups.ys[g]),
            span=float(spans[g]),
            full_span=full_span,
            last=bool(g == partial.size - 1),
        )
        for g in np.flatnonzero(partial)
    ]


def groups_for_frames(trigs, n_frames):
    """Match detector frames to position groups by trigger value.

    Detector frame ``k`` is exposed during the ckIM ticks carrying eigerTrig
    value ``k + FIRST_TRIGGER`` (see that constant). Returns
    ``(frames, groups)``: the frame indices that have a matching trigger
    group, and the group index each maps to. Frames whose trigger group is
    absent — because the position stream stopped first, or started late, or
    because a group is missing from the middle — are simply left out, instead
    of silently shifting the alignment of everything after them.
    """
    trigs = np.asarray(trigs)
    frames = np.arange(int(n_frames))
    if trigs.size == 0 or frames.size == 0:
        return frames[:0], np.zeros(0, dtype=np.intp)
    # Count from the first group's trigger value rather than from a literal 1,
    # so this still works when the MCS counter was not zeroed at the start of
    # the scan (scan 12 in the 26-3 test folder starts at 834).
    want = FIRST_TRIGGER + frames
    # Don't assume the column is sorted: it is non-decreasing in a normal
    # scan, but a counter that wraps mid-scan is not. Sorting first also
    # makes duplicates resolve to their earliest group.
    order = np.argsort(trigs, kind="stable")
    ordered = trigs[order]
    pos = np.searchsorted(ordered, want)
    ok = pos < ordered.size
    ok[ok] = ordered[pos[ok]] == want[ok]
    return frames[ok], order[pos[ok]]


def group_for_frame(trigs, frame):
    """Index of the position group holding detector frame ``frame``, or None."""
    frames, groups = groups_for_frames(trigs, int(frame) + 1)
    if frames.size == 0 or frames[-1] != frame:
        return None
    return int(groups[-1])


def reduce_position_stream(scan_number, folder, fname_format=FNAME_FORMAT):
    """Aggregate the position stream into per-image x, y, and I0 values.

    The position stream is sampled at a higher rate than the detector
    triggers, so multiple position samples share the same trigger index.
    This function groups samples by trigger and reduces each group to a
    single value per image.

    Parameters
    ----------
    scan_number : int
        Scan number used to build the HDF5 file name.
    folder : str
        Root folder containing the ``pos_stream`` subfolder.
    fname_format : str, optional
        Format string with three placeholders: ``folder``, stream name,
        and scan number. Defaults to :data:`FNAME_FORMAT`.

    Returns
    -------
    groups : PositionGroups
        Per-trigger-group arrays: ``i0s`` (range of column 0 within the
        group), ``xs``/``ys`` (mean motor position, columns 3 and 4),
        ``dxs``/``dys`` (position range within the group — a motor-jitter
        diagnostic), plus ``trigs``, ``counts`` and ``partial``, which are
        what let a caller match groups to detector frames and spot groups
        that only saw part of an exposure.

    Notes
    -----
    Assumes the trigger column (column 2) is monotonically non-decreasing,
    which is true for normal flyscan acquisition. The reduction uses
    ``np.add.reduceat`` and friends so the cost is O(n_samples) regardless
    of the number of images.
    """
    fname = fname_format.format(folder, "pos_stream", scan_number)
    pos_raw = File(fname)["entry/data/data"][()]

    # The hardware sometimes records trailing samples after the scan
    # finishes. Truncate at the first place where the sample counter
    # (column 1) stops advancing. If the counter advances monotonically
    # all the way to the end, there are no trailing samples to drop.
    plateau = np.where((pos_raw[1:, 1] - pos_raw[:-1, 1]) == 0)[0]
    last_point = plateau[0] + 1 if plateau.size else pos_raw.shape[0]
    pos = pos_raw[:last_point, :]

    # Column 2 is the detector trigger index. Samples with the same
    # trigger index belong to the same image exposure.
    trig = pos[:, 2]

    # `starts` holds the first sample index of each trigger group; together
    # with `counts` it defines the contiguous slices used by `reduceat`.
    starts = np.concatenate(([0], np.where(np.diff(trig) != 0)[0] + 1))
    counts = np.diff(np.concatenate((starts, [len(trig)])))

    # Per-image mean position and per-image position range.
    xs = np.add.reduceat(pos[:, 3], starts) / counts
    ys = np.add.reduceat(pos[:, 4], starts) / counts
    dxs = np.maximum.reduceat(pos[:, 3], starts) - np.minimum.reduceat(
        pos[:, 3], starts
    )
    dys = np.maximum.reduceat(pos[:, 4], starts) - np.minimum.reduceat(
        pos[:, 4], starts
    )

    # I0 is recorded as a running counter; the per-image value is the
    # difference between its max and min within the trigger window.
    i0s = np.maximum.reduceat(pos[:, 0], starts) - np.minimum.reduceat(
        pos[:, 0], starts
    )

    coverage = group_coverage(pos[:, 1], starts, len(trig))
    return PositionGroups(
        i0s=i0s,
        xs=xs,
        ys=ys,
        dxs=dxs,
        dys=dys,
        trigs=trig[starts].astype(np.int64),
        counts=counts,
        coverage=coverage,
        partial=flag_partial_groups(coverage),
    )


def process_position_stream(scan_number, folder, fname_format=FNAME_FORMAT):
    """``(i0s, xs, ys, dxs, dys)`` from :func:`reduce_position_stream`.

    Kept for callers that predate :class:`PositionGroups`; new code should
    use :func:`reduce_position_stream` so it also gets the trigger values
    and the partial-group flags.
    """
    g = reduce_position_stream(scan_number, folder, fname_format=fname_format)
    return g.i0s, g.xs, g.ys, g.dxs, g.dys


@dataclass
class ScanGeometry:
    """Beam and detector geometry recorded for a scan.

    Beam centre and detector distance come from the Eiger file's per-frame
    NDAttributes (read at frame 0); wavelength and energy come from the
    bluesky metadata in the scan's master file, which is the only place they
    are written. Any field that could not be read is None.
    """

    beam_center_x: float = None
    beam_center_y: float = None
    beam_center_pv: float = None  # the 4idgSoftX:Eiger:Center PV
    detector_distance: float = None
    wavelength: float = None
    wavelength_units: str = ""
    energy: float = None
    energy_units: str = ""


def _master_path(folder, scan_number):
    for fmt in ("{}/scan_{:06d}_master.hdf", "{}/scan_{:d}_master.hdf"):
        path = fmt.format(folder, scan_number)
        if os.path.exists(path):
            return path
    return None


def _beam_from_master(path):
    """Pull the beam block out of the master file's diffractometer metadata."""
    with File(path, "r") as f:
        key = "entry/instrument/bluesky/metadata/diffractometers"
        if key not in f:
            return {}
        blob = f[key][()]
    if isinstance(blob, bytes):
        blob = blob.decode("utf-8", errors="replace")
    try:
        import yaml

        doc = yaml.safe_load(blob) or {}
        for entry in doc.values():
            beam = (entry or {}).get("beam")
            if beam:
                return beam
        return {}
    except Exception:
        # No PyYAML, or the blob isn't what we expect: pick the scalars out
        # directly. `wavelength_PV` and `wavelength_units` don't match, since
        # the key must be followed immediately by the colon.
        out = {}
        for field in ("wavelength", "energy"):
            m = re.search(rf"^\s*{field}:\s*([-+0-9.eE]+)\s*$", blob, re.M)
            if m:
                out[field] = float(m.group(1))
            m = re.search(rf"^\s*{field}_units:\s*(\S+)\s*$", blob, re.M)
            if m:
                out[field + "_units"] = m.group(1)
        return out


def read_scan_geometry(folder, scan_number, fname_format=FNAME_FORMAT):
    """Beam centre, wavelength and detector distance for one scan.

    Reads whatever it can and leaves the rest None — the two files are
    written by different parts of the acquisition and either may be absent.
    """
    geo = ScanGeometry()

    path = fname_format.format(folder, "eiger", scan_number)
    if not os.path.exists(path):
        path = "{}/{}/scan_{:d}.h5".format(folder, "eiger", scan_number)
    if os.path.exists(path):
        try:
            with File(path, "r") as f:
                attrs = f["entry/instrument/NDAttributes"]
                for name, field in (
                    ("BeamCenterX", "beam_center_x"),
                    ("BeamCenterY", "beam_center_y"),
                    ("BeamCenterPV", "beam_center_pv"),
                    ("DetectorDistance", "detector_distance"),
                ):
                    if name in attrs and attrs[name].shape[0]:
                        setattr(geo, field, float(attrs[name][0]))
        except (OSError, KeyError):
            pass

    master = _master_path(folder, scan_number)
    if master is not None:
        try:
            beam = _beam_from_master(master)
        except (OSError, KeyError):
            beam = {}
        for field in ("wavelength", "energy"):
            if beam.get(field) is not None:
                setattr(geo, field, float(beam[field]))
            units = beam.get(field + "_units")
            if units:
                setattr(geo, field + "_units", str(units))
    return geo


def process_images_v0(scan_number, folder, fname_format=FNAME_FORMAT):
    """Sum a fixed ROI across all Eiger images in a scan (eager version).

    Loads the entire image stack into memory before slicing. Kept as a
    reference implementation; use :func:`process_images` for large scans
    to avoid running out of memory.

    Parameters
    ----------
    scan_number : int
        Scan number used to build the HDF5 file name.
    folder : str
        Root folder containing the ``eiger`` subfolder.
    fname_format : str, optional
        Format string for the HDF5 path. Defaults to :data:`FNAME_FORMAT`.

    Returns
    -------
    sums : ndarray, shape (n_images,)
        ROI intensity sum for each image.
    """
    fname = fname_format.format(folder, "eiger", scan_number)
    images = File(fname)["entry/data/data"][()]
    rois = ((678, 888), (333, 480))
    sums = images[:, rois[1][0] : rois[1][1], rois[0][0] : rois[0][1]].sum(
        axis=(1, 2)
    )
    return sums


def process_images(
    scan_number,
    folder,
    roi=((678, 888), (333, 480)),
    batch=100,
    fname_format=FNAME_FORMAT,
):
    """Sum a fixed ROI across all Eiger images, reading in batches.

    Reads only the ROI window of ``batch`` frames at a time, so peak
    memory is bounded by ``batch * roi_height * roi_width`` regardless of
    the total number of images in the scan.

    Parameters
    ----------
    scan_number : int
        Scan number used to build the HDF5 file name.
    folder : str
        Root folder containing the ``eiger`` subfolder.
    roi : tuple of tuple, optional
        ROI as ``((col_start, col_stop), (row_start, row_stop))``. The
        column range is the x extent and the row range is the y extent.
    batch : int, optional
        Number of frames to read per HDF5 access. Larger values reduce
        per-call overhead at the cost of higher peak memory.
    fname_format : str, optional
        Format string for the HDF5 path. Defaults to :data:`FNAME_FORMAT`.

    Returns
    -------
    sums : ndarray, shape (n_images,)
        ROI intensity sum for each image, stored as float64.
    """
    fname = fname_format.format(folder, "eiger", scan_number)
    # Unpack the ROI; note `roi[0]` is the x (column) range and `roi[1]`
    # is the y (row) range, matching NumPy axis order (rows, cols).
    r0, r1 = roi[1]
    c0, c1 = roi[0]
    with File(fname, "r") as f:
        dset = f["entry/data/data"]
        n = dset.shape[0]
        sums = np.empty(n, dtype=np.float64)
        # Stream the stack in `batch`-sized chunks. h5py reads only the
        # ROI bytes from disk, and the sum keeps memory flat across loops.
        for i in range(0, n, batch):
            j = min(i + batch, n)
            sums[i:j] = dset[i:j, r0:r1, c0:c1].sum(axis=(1, 2))
        return sums


def process_vortex(
    scan_number, folder, roi=(800, 900), batch=100, fname_format=FNAME_FORMAT
):
    """Sum a Vortex fluorescence ROI across all spectra in a scan.

    Reads the Vortex spectrum stack in batches and integrates a single
    energy window (channel range) per spectrum. Memory use is bounded by
    ``batch * n_channels_per_spectrum``.

    Parameters
    ----------
    scan_number : int
        Scan number used to build the HDF5 file name.
    folder : str
        Root folder containing the ``vortex`` subfolder.
    roi : tuple of int, optional
        Energy channel range ``(start, stop)`` to integrate over.
    batch : int, optional
        Number of spectra to read per HDF5 access.
    fname_format : str, optional
        Format string for the HDF5 path. Defaults to :data:`FNAME_FORMAT`.

    Returns
    -------
    sums : ndarray, shape (n_spectra,)
        Integrated counts in the energy window for each spectrum.
    """
    fname = fname_format.format(folder, "vortex", scan_number)
    e0, e1 = roi
    with File(fname, "r") as f:
        dset = f["entry/data/data"]
        n = dset.shape[0]
        sums = np.empty(n, dtype=np.float64)
        # Sum over both the detector-element axis (axis 1) and the
        # selected energy-channel slice (axis 2) for each spectrum.
        for i in range(0, n, batch):
            j = min(i + batch, n)
            sums[i:j] = dset[i:j, :, e0:e1].sum(axis=(1, 2))
    return sums


def plot_data(x, y, z, i0=None, scatter=False, **kwargs):
    """Plot a scalar field ``z(x, y)`` from a flyscan as scatter or contour.

    Drops the first sample and, if needed, the last ``z`` value so the
    three arrays align. Optionally normalizes by an incident-intensity
    monitor ``i0``.

    Parameters
    ----------
    x, y : array_like
        Per-image x and y positions, typically the means returned by
        :func:`process_position_stream`.
    z : array_like
        Per-image scalar values (e.g. ROI sum from :func:`process_images`).
        May be one shorter than ``x``/``y`` when image triggers and
        position samples differ by one; in that case it is used as-is.
    i0 : array_like or None, optional
        Per-image incident intensity. When provided, ``z`` is divided by
        ``i0`` element-wise before plotting.
    scatter : bool, optional
        If True, use ``plt.scatter`` with one marker per point. If False
        (default), use ``plt.tricontourf`` for a filled contour plot.
    **kwargs
        Forwarded to the underlying matplotlib call.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    cbar : matplotlib.colorbar.Colorbar
    """
    # Drop the first position sample; the detector skips the first trigger.
    xi = x[1:]
    yi = y[1:]

    # `z` may already match `yi` in length; otherwise trim its trailing element.
    zi = z[:-1] if z.size != yi.size else z

    if i0 is not None:
        # Align i0 to z the same way before normalizing.
        i0i = i0[:-1] if i0.size != zi.size else i0
        zi /= i0i

    fig, ax = plt.subplots()

    if scatter:
        plt.scatter(xi, yi, c=zi, **kwargs)
    else:
        plt.tricontourf(xi, yi, zi, **kwargs)

    cbar = plt.colorbar()

    return fig, ax, cbar
