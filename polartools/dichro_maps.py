"""Shared helpers for dichroism difference maps from flyscan pairs.

A dichro measurement is a pair of flyscans taken at the two phase-plate
settings; the map of interest is the difference between them. These functions
find the pairs in a scan range, put the two members on a common grid, and plot
a series of the resulting maps.

Factored out of ``data/experiment_plots.ipynb`` so several notebooks can share
one copy. The reduction itself lives in :mod:`process_flyscan`.
"""

import warnings
from pathlib import Path

import h5py
import hdf5plugin  # noqa: F401  -- registers the HDF5 codecs these files use
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.signal import savgol_filter

from .process_flyscan import (
    groups_for_frames,
    process_images,
    reduce_position_stream,
)

# Phase-plate PZT. The two members of a pair differ in this and nothing else.
SWITCH_KEY = "pr2_pzt_localdc"


def baseline_value(scan, key, folder):
    """One signal from a scan's baseline stream, or None if unavailable.

    The baseline is the trustworthy source for per-scan device values. The
    ``diffractometers`` metadata blob also carries a mono energy, but it can be
    stale -- it read 6.200 keV for a whole energy series -- so labels must not
    come from there.
    """
    path = Path(folder) / f"scan_{scan:06d}_master.hdf"
    if not path.exists():
        return None
    try:
        with h5py.File(path, "r") as f:
            node = f"entry/instrument/bluesky/streams/baseline/{key}/value"
            if node not in f:
                return None
            return float(np.atleast_1d(f[node][()])[0])
    except OSError:
        # Still being written -- HDF5 cannot take a read lock. Treat as absent;
        # a re-run picks it up.
        return None


def eiger_frame_count(scan, folder):
    """Frames in a scan's Eiger file, without reading any of them."""
    with h5py.File(Path(folder) / "eiger" / f"scan_{scan:06d}.h5", "r") as f:
        return f["entry/data/data"].shape[0]


def find_pairs(scans, label_key, folder, require_positions=True, atol=0.0):
    """Consecutive scans sharing a label value but differing in phase plate.

    Returns ``[(label, scan_a, scan_b), ...]`` with ``scan_a`` the lower number.
    Pairs without a position stream are dropped unless asked for.

    ``atol`` is how far the label may differ between the two members and still
    count as one setting. Motor readbacks jitter -- ``huber_hp_z`` gives 0.0500
    for one member and 0.0501 for the other -- so those series need a
    tolerance, well below the spacing between settings. The default of 0 means
    exact equality, right for setpoint-valued series.
    """
    seen = []
    for scan in scans:
        label = baseline_value(scan, label_key, folder)
        switch = baseline_value(scan, SWITCH_KEY, folder)
        if label is not None and switch is not None:
            seen.append((scan, label, switch))
    pairs, i = [], 0
    while i < len(seen) - 1:
        (sa, la, wa), (sb, lb, wb) = seen[i], seen[i + 1]
        if sb == sa + 1 and abs(la - lb) <= atol and wa != wb:
            if not require_positions or all(
                (Path(folder) / "pos_stream" / f"scan_{s:06d}.h5").exists()
                for s in (sa, sb)
            ):
                pairs.append(((la + lb) / 2, sa, sb))
            i += 2
        else:
            i += 1
    return pairs


def normalised_map(scan, roi, folder, batch=200):
    """``(x, y, z/I0)`` for one flyscan.

    Detector frame k is exposed during the position samples carrying eigerTrig
    k+1, and ``groups_for_frames`` matches them by that trigger value, so a
    frame whose trigger group is missing drops out instead of shifting every
    later frame along by one.
    """
    z = process_images(scan, folder, roi=roi, batch=batch)
    groups = reduce_position_stream(scan, folder)
    frames, gi = groups_for_frames(groups.trigs, z.size)
    i0 = groups.i0s[gi]
    with np.errstate(divide="ignore", invalid="ignore"):
        return (
            groups.xs[gi],
            groups.ys[gi],
            z[frames] / np.where(i0 != 0, i0, np.nan),
        )


def absorbance_map(
    scan, roi, folder, batch=200, reference="open", open_pct=99.0
):
    """``(x, y, mu*t)`` for one flyscan.

    ``reference="open"`` (default) gives **true absorbance**,
    ``mu*t = ln(I_open / I)``: zero where the beam is clear, rising to the
    sample's own absorbance. ``I_open`` is the ``open_pct`` percentile of
    ``I/I0`` over the scan, which is the unobstructed level provided the scan
    actually runs off the sample at one end -- if it does not, the reference
    lands on the thinnest part of the sample and every value is offset low.

    ``reference="i0"`` uses ``ln(I0 / I)`` literally. ``I`` is a detector sum
    and ``I0`` a monitor counter, so that carries a constant ln(gain) offset
    (about -11.9 in the open beam here).

    The offset cancels in a CL - CR *difference* either way, but **not** in
    their sum, so anything normalising by ``mu_L t + mu_R t`` needs the open
    reference to mean anything.

    Absorbance, not transmission, is what is linear in path length, so this is
    the quantity a tomographic reconstruction wants.
    """
    x, y, ratio = normalised_map(scan, roi, folder, batch)
    with np.errstate(divide="ignore", invalid="ignore"):
        positive = np.where(ratio > 0, ratio, np.nan)
        if reference == "open":
            open_level = np.nanpercentile(positive, open_pct)
            if not np.isfinite(open_level) or open_level <= 0:
                open_level = 1.0
            return x, y, -np.log(positive / open_level)
        if reference == "i0":
            return x, y, -np.log(positive)
    raise ValueError(f"reference must be 'open' or 'i0', got {reference!r}")


def dichroic_ratio(diff, mean, floor_frac=0.1):
    """``(mu_L t - mu_R t) / (mu_L t + mu_R t)`` from a difference and mean.

    The denominator is twice the mean absorbance, which goes to zero wherever
    the beam is clear -- the ratio is 0/0 in the open beam and would amplify
    noise without bound there. Points whose sum falls below ``floor_frac`` of
    its 99th percentile are masked out instead of returned as huge numbers.

    Being a ratio of absorbances it divides out the sample thickness, so it
    measures the magnetisation itself rather than magnetisation x path length.
    """
    total = 2.0 * np.asarray(mean, dtype=float)
    scale = np.nanpercentile(np.abs(total), 99)
    if not np.isfinite(scale) or scale <= 0:
        return np.full_like(total, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(
            np.abs(total) > floor_frac * scale,
            np.asarray(diff, dtype=float) / total,
            np.nan,
        )


def pair_difference(
    scan_a, scan_b, roi, folder, grid_n=200, quantity="transmission"
):
    """Regrid both scans of a pair and return ``(XG, YG, difference, mean)``.

    The two scans don't visit the same points -- the stage free-runs -- so each
    is interpolated onto a common grid before subtracting. The grid spans only
    the overlap of the two scans' travel, so neither is extrapolated.

    ``quantity`` selects what is differenced: ``"transmission"`` gives
    ``I/I0``, ``"absorbance"`` gives ``mu*t = ln(I0/I)``. The transform is
    applied per point *before* interpolating, because a log of an interpolated
    value is not the interpolation of the logs.
    """
    reader = absorbance_map if quantity == "absorbance" else normalised_map
    xa, ya, za = reader(scan_a, roi, folder)
    xb, yb, zb = reader(scan_b, roi, folder)
    xg = np.linspace(max(xa.min(), xb.min()), min(xa.max(), xb.max()), grid_n)
    yg = np.linspace(max(ya.min(), yb.min()), min(ya.max(), yb.max()), grid_n)
    XG, YG = np.meshgrid(xg, yg)
    ga = griddata((xa, ya), za, (XG, YG), method="linear")
    gb = griddata((xb, yb), zb, (XG, YG), method="linear")
    return XG, YG, gb - ga, (ga + gb) / 2


def average_pairs(pairs, roi, folder, grid_n=200):
    """Average several dichro pairs measured at the same position.

    ``pairs`` is ``[(scan_a, scan_b), ...]``. Every pair is put on one common
    grid -- the intersection of all their travel -- before differencing, since
    each is otherwise regridded to its own overlap and they share no sample
    points. Each pair's own common-mode offset is removed first: those differ
    from pair to pair and would bias the mean.

    Returns ``(XG, YG, mean_difference, mean_transmission, stack)``, where
    ``stack`` holds the individual offset-removed differences so the caller can
    see how far the repeats agree.
    """
    loaded = []
    for scan_a, scan_b in pairs:
        xa, ya, za = normalised_map(scan_a, roi, folder)
        xb, yb, zb = normalised_map(scan_b, roi, folder)
        loaded.append((xa, ya, za, xb, yb, zb))

    x0 = max(max(m[0].min(), m[3].min()) for m in loaded)
    x1 = min(min(m[0].max(), m[3].max()) for m in loaded)
    y0 = max(max(m[1].min(), m[4].min()) for m in loaded)
    y1 = min(min(m[1].max(), m[4].max()) for m in loaded)
    XG, YG = np.meshgrid(
        np.linspace(x0, x1, grid_n), np.linspace(y0, y1, grid_n)
    )

    diffs, means = [], []
    for xa, ya, za, xb, yb, zb in loaded:
        ga = griddata((xa, ya), za, (XG, YG), method="linear")
        gb = griddata((xb, yb), zb, (XG, YG), method="linear")
        d = gb - ga
        diffs.append(d - np.nanmedian(d))
        means.append((ga + gb) / 2)
    stack = np.stack(diffs)
    with warnings.catch_warnings():  # rim cells lie outside every map
        warnings.simplefilter("ignore", RuntimeWarning)
        return (
            XG,
            YG,
            np.nanmean(stack, axis=0),
            np.nanmean(np.stack(means), axis=0),
            stack,
        )


def averaged_series(groups, roi, folder, grid_n=200, verbose=True):
    """Build plottable entries from ``[(title, [(a, b), ...]), ...]``.

    One entry per group, holding the mean of that group's repeats. The offsets
    are already removed inside :func:`average_pairs`, so ``offset`` is zero.
    """
    series = []
    for i, (title, pairs) in enumerate(groups):
        XG, YG, diff, mean, stack = average_pairs(pairs, roi, folder, grid_n)
        entry = dict(
            label=float(i),
            a=pairs[0][0],
            b=pairs[-1][1],
            XG=XG,
            YG=YG,
            diff=diff,
            mean=mean,
            offset=0.0,
            stack=stack,
            title=f"{title}  ({len(pairs)}\u00d7)",
            amplitude=float(np.nanstd(diff)),
        )
        series.append(entry)
        if verbose:
            msg = f"  {title:<24} {len(pairs)} pair(s)"
            if len(pairs) > 1:
                # Spread between repeats, halved, estimates one map's noise.
                # Measured robustly: griddata throws a handful of wild values
                # onto the rim of the interpolated area, and a plain std is
                # dominated by them -- for #778-781 four rim pixels inflate it
                # from 8 to 104.
                d = stack[0] - stack[1]
                d = d[np.isfinite(d)]
                single = (
                    1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)
                )
                msg += (
                    f"   single-map noise {single:6.2f}"
                    f"   expected for the mean {single / np.sqrt(len(pairs)):6.2f}"
                )
            print(msg)
    return series


def build_series(pairs, roi, folder, grid_n=200, verbose=True):
    """Difference map for every pair, ordered by the series parameter.

    Sorted by label rather than scan number: settings are not always visited in
    order, and a trend drawn in acquisition order would double back on itself.
    A pair that fails to load -- a scan still being written -- is skipped rather
    than losing the ones that did.
    """
    series = []
    for entry in pairs:
        # A pair may carry an explicit panel title, for groups whose members
        # differ in more than one setting and so have no single label.
        label, a, b = entry[:3]
        title = entry[3] if len(entry) > 3 else None
        try:
            XG, YG, diff, mean = pair_difference(a, b, roi, folder, grid_n)
        except (OSError, KeyError, ValueError) as exc:
            print(f"  {label:>8.3f}  #{a}/{b}:  skipped ({type(exc).__name__})")
            continue
        offset = float(np.nanmedian(diff))
        series.append(
            dict(
                label=label,
                a=a,
                b=b,
                XG=XG,
                YG=YG,
                diff=diff,
                mean=mean,
                offset=offset,
                title=title,
                amplitude=float(np.nanstd(diff - offset)),
            )
        )
        if verbose:
            print(
                f"  {label:>8.3f}  #{a}/{b}:  offset {offset:>10.2f}   "
                f"structure std {series[-1]['amplitude']:>8.2f}"
            )
    series.sort(key=lambda s: s["label"])
    return series


def plot_difference_series(
    series, title, label_fmt="{:.3f}", ncol=3, clip=99.0, share=True
):
    """Difference maps on one shared symmetric scale, common offsets removed.

    The two phase-plate settings don't transmit equally, so every difference
    sits on a large constant offset that would swamp the structure. It is
    instrumental, so it is subtracted per panel.

    ``share`` links the panel axes, which is what makes a series comparable
    point to point. Turn it off when the members sit at different places on
    the sample: shared limits would then span the whole excursion and squash
    every map into a sliver.
    """
    nrow = int(np.ceil(len(series) / ncol))
    deviation = np.concatenate(
        [(s["diff"] - s["offset"]).ravel() for s in series]
    )
    v = np.nanpercentile(np.abs(deviation), clip)
    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(3.6 * ncol, 3.5 * nrow),
        constrained_layout=True,
        sharex=share,
        sharey=share,
    )
    axes = np.atleast_1d(axes).ravel()
    for ax, s in zip(axes, series):
        pcm = ax.pcolormesh(
            s["XG"],
            s["YG"],
            s["diff"] - s["offset"],
            shading="auto",
            cmap="RdBu_r",
            vmin=-v,
            vmax=v,
        )
        ax.set_aspect("equal")
        heading = s.get("title") or label_fmt.format(s["label"])
        ax.set_title(f"{heading}\n#{s['a']}/{s['b']}", fontsize=8)
    for ax in axes[len(series) :]:
        ax.set_visible(False)
    fig.colorbar(
        pcm,
        ax=axes.tolist(),
        shrink=0.5,
        label="difference b − a (I0-normalised, offset removed)",
    )
    fig.suptitle(title, fontsize=13)
    return fig


def plot_stitched(
    series,
    title,
    label_fmt="{:.4f}",
    clip=99.0,
    width=5.0,
    max_height=34.0,
    cmap="seismic",
):
    """Lay every map of a series at its own coordinates on a single axes.

    For a set of maps taken at stepped sample positions: rather than one panel
    each, they are drawn where they actually sit, so the strip reads as one
    picture of the sample. The tiles are positioned by the interferometer
    coordinates, which do follow the coarse stage -- a 10 um move in
    ``huber_hp_y`` shows up as about +10000 here.

    Each tile still has its own common-mode offset removed. Without that the
    instrumental offset, which differs from pair to pair, would put a visible
    step at every tile boundary.

    The default colour map is ``seismic`` rather than the ``RdBu_r`` used for
    panel figures. Most of a dichroic difference sits well inside the clipped
    range -- the median deviation is under a third of the 99th percentile --
    and ``RdBu_r`` spends that whole span near white, which washes the domains
    out. ``seismic`` changes lightness about four times faster at zero, so the
    structure carries. It saturates harder at the ends, which also makes
    per-scan artefact bands obvious.
    """
    deviation = np.concatenate(
        [(s["diff"] - s["offset"]).ravel() for s in series]
    )
    v = np.nanpercentile(np.abs(deviation), clip)

    x0 = min(s["XG"].min() for s in series)
    x1 = max(s["XG"].max() for s in series)
    y0 = min(s["YG"].min() for s in series)
    y1 = max(s["YG"].max() for s in series)
    height = min(max_height, width * (y1 - y0) / (x1 - x0))

    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    for s in series:
        pcm = ax.pcolormesh(
            s["XG"],
            s["YG"],
            s["diff"] - s["offset"],
            shading="auto",
            cmap=cmap,
            vmin=-v,
            vmax=v,
        )
        heading = s.get("title") or label_fmt.format(s["label"])
        ax.text(
            x1,
            float(np.mean(s["YG"])),
            f" {heading}\n #{s['a']}/{s['b']}",
            va="center",
            ha="left",
            fontsize=7,
            color="0.3",
        )
    ax.set_aspect("equal")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel("x (interferometer units, ~1 unit = 1 nm)")
    ax.set_ylabel("y (interferometer units, ~1 unit = 1 nm)")
    fig.colorbar(
        pcm,
        ax=ax,
        shrink=0.25,
        pad=0.16,
        label="difference b \u2212 a (I0-normalised, offset removed)",
    )
    ax.set_title(title, fontsize=11)
    return fig


def write_gif(
    series,
    path,
    title,
    label_fmt="{:.1f}",
    clip=99.0,
    fps=2.5,
    cmap="RdBu_r",
    width=5.4,
    pingpong=True,
    dpi=100,
    scale="frame",
):
    """Animate a series of difference maps as a GIF, one frame per member.

    All frames share one colour scale and one set of axis limits, so what moves
    between frames is the data and not the rendering. ``pingpong`` plays the
    sweep forwards then backwards, which avoids the jump a looping monotonic
    sweep otherwise makes when it wraps from the last setting back to the first.

    The writer is driven directly rather than through ``FuncAnimation``: that
    calls the draw function an extra time for setup, and the correspondence
    between requested frames and grabbed images is then not one to one -- the
    first angle of a sweep goes missing. Here one explicit grab follows each
    draw.

    Returns the path written.
    """
    from matplotlib.animation import PillowWriter

    # One scale for the whole animation, but chosen so no single frame blows
    # through it. Pooling every frame's pixels and taking one percentile lets
    # the many weak frames set the limit, and the few strong ones then clip
    # badly -- in the 807-842 sweep that cost 9% of the pixels at xr 27.5 while
    # every other frame lost under 1%. Taking the largest per-frame percentile
    # instead caps the worst frame at the same (100 - clip)%.
    per_frame = [
        np.nanpercentile(np.abs(s["diff"] - s["offset"]), clip) for s in series
    ]
    if scale == "frame":
        v = float(np.max(per_frame))
    elif scale == "pooled":
        v = float(
            np.nanpercentile(
                np.abs(
                    np.concatenate(
                        [(s["diff"] - s["offset"]).ravel() for s in series]
                    )
                ),
                clip,
            )
        )
    else:
        raise ValueError(f"scale must be 'frame' or 'pooled', got {scale!r}")
    x0 = min(s["XG"].min() for s in series)
    x1 = max(s["XG"].max() for s in series)
    y0 = min(s["YG"].min() for s in series)
    y1 = max(s["YG"].max() for s in series)

    order = list(range(len(series)))
    if pingpong and len(series) > 2:
        order += order[-2:0:-1]

    fig, ax = plt.subplots(
        figsize=(width, width * (y1 - y0) / (x1 - x0) + 1.2),
        constrained_layout=True,
    )
    ax.set_aspect("equal")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel("x (interferometer units, ~1 unit = 1 nm)")
    ax.set_ylabel("y (interferometer units, ~1 unit = 1 nm)")
    # The series name goes on the figure so it can use the full width; only the
    # short per-frame label sits on the axes, clear of the colour bar.
    fig.suptitle(title, fontsize=9)

    first = series[0]
    mesh = ax.pcolormesh(
        first["XG"],
        first["YG"],
        first["diff"] - first["offset"],
        shading="auto",
        cmap=cmap,
        vmin=-v,
        vmax=v,
    )
    fig.colorbar(
        mesh,
        ax=ax,
        shrink=0.8,
        label="difference b \u2212 a (I0-normalised, offset removed)",
    )

    writer = PillowWriter(fps=fps)
    with writer.saving(fig, str(path), dpi=dpi):
        for i in order:
            s = series[i]
            for coll in list(ax.collections):
                # A fresh mesh per frame: pairs are regridded separately, so the
                # meshes are not guaranteed to share a shape.
                coll.remove()
            ax.pcolormesh(
                s["XG"],
                s["YG"],
                s["diff"] - s["offset"],
                shading="auto",
                cmap=cmap,
                vmin=-v,
                vmax=v,
            )
            label = s.get("title") or label_fmt.format(s["label"])
            ax.set_title(f"{label}   #{s['a']}/{s['b']}", fontsize=11)
            writer.grab_frame()
    plt.close(fig)
    return path


def line_cut(entry, y0, robust=8.0):
    """One row of a map at the y nearest ``y0``.

    Returns ``(y_actual, x, i_a, i_b, mean, diff, kept)``. The two polarisations
    are recovered from the stored mean and difference (``diff = i_b - i_a``),
    so no re-read is needed.

    ``griddata`` leaves a few wild values on the rim of the interpolated area --
    a single such column sits at the left edge of the #778-781 map, where the
    difference reads -12709 against a row median of -35. They are dropped here
    rather than allowed to dominate the derivative; ``kept`` is the mask.
    """
    yv = entry["YG"][:, 0]
    j = int(np.argmin(np.abs(yv - y0)))
    x = entry["XG"][j].astype(float)
    diff = entry["diff"][j].astype(float) - entry["offset"]
    mean = entry["mean"][j].astype(float)

    kept = np.isfinite(x) & np.isfinite(diff) & np.isfinite(mean)
    for row in (diff, mean):
        good = row[kept]
        if good.size:
            med = np.median(good)
            mad = 1.4826 * np.median(np.abs(good - med))
            if mad > 0:
                kept &= np.abs(row - med) <= robust * mad
    x, diff, mean = x[kept], diff[kept], mean[kept]
    return float(yv[j]), x, mean - diff / 2, mean + diff / 2, mean, diff, kept


def _derivative(x, y, window=11, order=2):
    """d/dx, lightly smoothed. Falls back to a plain gradient if too short."""
    raw = np.gradient(y, x)
    if y.size < window + 2:
        return raw, raw
    step = float(np.median(np.diff(x)))
    return raw, savgol_filter(y, window, order, deriv=1, delta=step)


def peak_fwhm(x, y, k):
    """Full width at half maximum of the peak of ``y`` at index ``k``.

    Walks out to each side until the curve falls through half the peak height
    and interpolates the crossing. Returns ``(fwhm, x_lo, x_hi)``, with NaNs if
    the peak does not come back down inside the cut.
    """
    half = y[k] / 2.0
    sign = np.sign(y[k]) or 1.0
    i = k
    while i > 0 and sign * y[i] > sign * half:
        i -= 1
    j = k
    while j < y.size - 1 and sign * y[j] > sign * half:
        j += 1
    if sign * y[i] > sign * half or sign * y[j] > sign * half:
        return float("nan"), float("nan"), float("nan")
    lo = np.interp(
        half, [y[i], y[i + 1]][:: int(sign)], [x[i], x[i + 1]][:: int(sign)]
    )
    hi = np.interp(
        half, [y[j], y[j - 1]][:: int(sign)], [x[j], x[j - 1]][:: int(sign)]
    )
    return float(hi - lo), float(lo), float(hi)


def plot_line_cut(
    entry, y0, title=None, window=11, figsize=(11, 3.3), clip=99.0
):
    """Map, the dichroic difference along one row, and its x-derivative.

    The derivative peaks where the difference crosses between domains, so its
    extrema locate the walls and their width sets an upper bound on the probe
    size convolved with the true wall width.
    """
    y_at, x, _i_a, _i_b, _mean, diff, kept = line_cut(entry, y0)
    fig, ax = plt.subplots(
        1, 3, figsize=figsize, gridspec_kw=dict(width_ratios=[1, 1.35, 1.35])
    )

    v = np.nanpercentile(np.abs(entry["diff"] - entry["offset"]), clip)
    ax[0].pcolormesh(
        entry["XG"],
        entry["YG"],
        entry["diff"] - entry["offset"],
        cmap="RdBu_r",
        vmin=-v,
        vmax=v,
        shading="auto",
    )
    ax[0].axhline(y_at, color="k", lw=1.2, ls="--")
    ax[0].set_title(f"{entry['title']}\ncut at y = {y_at:.0f}", fontsize=8)
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")
    ax[0].set_aspect("equal")

    ax[1].axhline(0, color="0.7", lw=0.8)
    ax[1].plot(x, diff, lw=1.2, color="tab:red")
    ax[1].set_ylabel("difference  (I$_{CL}$ - I$_{CR}$)")
    ax[1].set_title("dichroic difference", fontsize=8)

    raw, smooth = _derivative(x, diff, window)
    ax[2].axhline(0, color="0.7", lw=0.8)
    ax[2].plot(x, raw, lw=0.6, color="tab:red", alpha=0.25)
    ax[2].plot(x, smooth, lw=1.4, color="tab:red")
    k = int(np.abs(smooth).argmax())
    fwhm, lo, hi = peak_fwhm(x, smooth, k)
    ax[2].plot(x[k], smooth[k], "ko", ms=4)
    if np.isfinite(fwhm):
        ax[2].hlines(smooth[k] / 2, lo, hi, color="k", lw=1.1)
        ax[2].vlines([lo, hi], 0, smooth[k] / 2, color="k", lw=0.7, ls=":")
        ax[2].annotate(
            f"FWHM {fwhm:.0f}",
            ((lo + hi) / 2, smooth[k] / 2),
            textcoords="offset points",
            xytext=(0, -12),
            ha="center",
            fontsize=7,
        )
    ax[2].annotate(
        f"{smooth[k]:+.3f}\nat x = {x[k]:.0f}",
        (x[k], smooth[k]),
        textcoords="offset points",
        xytext=(6, -4),
        fontsize=7,
    )
    ax[2].set_ylabel("d(difference)/dx")
    ax[2].set_title(
        f"derivative (Savitzky-Golay, {window} pts;\nraw in pale)", fontsize=8
    )

    for a in ax[1:]:
        a.set_xlabel("x")
        a.grid(alpha=0.25)
    dropped = int((~kept).sum())
    if title:
        fig.suptitle(
            title + (f"   ({dropped} rim point(s) dropped)" if dropped else ""),
            fontsize=10,
        )
    fig.tight_layout()
    return fig, dict(
        y=y_at,
        x=x,
        diff=diff,
        ddiff=smooth,
        dropped=dropped,
        peak=float(smooth[k]),
        peak_x=float(x[k]),
        fwhm=fwhm,
    )


def edge_position(x, y, window=11):
    """Sub-pixel x of the steepest point of a step in ``y(x)``.

    The derivative peak is located on the smoothed curve, then refined by a
    parabola through the peak sample and its two neighbours, so the result is
    not quantised to the 50 nm grid step.

    Returns ``(x_edge, slope)``; ``slope`` carries the sign of the step.
    """
    _, g = _derivative(x, y, window)
    a = np.abs(g)
    k = int(np.argmax(a))
    if 0 < k < a.size - 1:
        denom = a[k - 1] - 2 * a[k] + a[k + 1]
        if denom != 0:
            shift = 0.5 * (a[k - 1] - a[k + 1]) / denom
            if abs(shift) <= 1:
                return float(np.interp(k + shift, np.arange(x.size), x)), float(
                    g[k]
                )
    return float(x[k]), float(g[k])


def _subpixel(v, k):
    """Refine the index of a peak in ``v`` at ``k`` by a three-point parabola."""
    if 0 < k < v.size - 1:
        den = v[k - 1] - 2 * v[k] + v[k + 1]
        if den != 0:
            return k + 0.5 * (v[k - 1] - v[k + 1]) / den
    return float(k)


def _resample(series, step=50.0, field="diff"):
    """Every frame's map on one common grid (the intersection of their travel)."""
    lo_x = max(e["XG"].min() for e in series)
    hi_x = min(e["XG"].max() for e in series)
    lo_y = max(e["YG"].min() for e in series)
    hi_y = min(e["YG"].max() for e in series)
    xg = np.arange(lo_x, hi_x, step)
    yg = np.arange(lo_y, hi_y, step)
    X, Y = np.meshgrid(xg, yg)
    out = []
    for e in series:
        z = e[field] - (e["offset"] if field == "diff" else 0.0)
        f = RegularGridInterpolator(
            (e["YG"][:, 0], e["XG"][0]),
            z,
            bounds_error=False,
            fill_value=np.nan,
        )
        out.append(f(np.stack([Y, X], -1)))
    return xg, yg, np.asarray(out)


def _xcorr2(ref, img):
    """(dx, dy, peak) in grid cells, from a normalised 2D cross-correlation."""
    a = np.nan_to_num(ref - np.nanmean(ref))
    b = np.nan_to_num(img - np.nanmean(img))
    a /= np.linalg.norm(a) or 1.0
    b /= np.linalg.norm(b) or 1.0
    c = np.fft.fftshift(
        np.real(np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))))
    )
    ky, kx = np.unravel_index(np.argmax(c), c.shape)
    cy, cx = np.array(c.shape) // 2
    # F(a).conj(F(b)) peaks at minus the displacement of b relative to a --
    # verified against synthetic shifts -- so negate to return the displacement
    # itself, which is what the caller subtracts off.
    return (
        cx - _subpixel(c[ky], kx),
        cy - _subpixel(c[:, kx], ky),
        float(c[ky, kx]),
    )


def align_series(
    series, method="map", y0=None, step=50.0, reference=None, max_shift=4000.0
):
    """Line frames up in x, undoing the xeryon stage runout.

    Rotating the xeryon translates the sample sideways, so the frames of an
    animation jump about. Each frame is registered against a reference frame and
    its x coordinates offset by the measured displacement. **Only the coordinates
    move** -- the data are never resampled, so no interpolation blur is added.

    ``method="map"`` cross-correlates the whole difference map, ``method="cut"``
    only the row nearest ``y0``. The map is the safer default: it uses every
    pixel, and it also measures the vertical component, which is the assumption
    the correction rests on. For scans 807-878 that assumption holds well --
    dy rms 34-36 nm against dx rms 326-373 nm.

    The reference defaults to the middle frame, which correlates better with both
    ends of a sweep than either end does with the other. ``max_shift`` leaves a
    frame alone if its measured displacement is implausibly large.

    Returns ``(aligned_series, report)``.
    """
    if reference is None:
        reference = len(series) // 2
    if method == "map":
        _xg, _yg, C = _resample(series, step)
        meas = [_xcorr2(C[reference], img) for img in C]
        dxs = np.array([m[0] * step for m in meas])
        dys = np.array([m[1] * step for m in meas])
        quality = [m[2] for m in meas]
    elif method == "cut":
        if y0 is None:
            raise ValueError("method='cut' needs y0")
        cuts = [line_cut(e, y0) for e in series]
        lo = max(c[1].min() for c in cuts)
        hi = min(c[1].max() for c in cuts)
        xg = np.arange(lo, hi, step)
        P = np.array([np.interp(xg, c[1], c[5]) for c in cuts])
        ref = P[reference] - P[reference].mean()
        ref /= np.linalg.norm(ref) or 1.0
        n = int(max_shift / step)
        dxs, quality = [], []
        for p in P:
            b = p - p.mean()
            b /= np.linalg.norm(b) or 1.0
            lags = np.arange(-n, n + 1)
            c = np.array(
                [
                    np.dot(
                        ref[max(0, l) : len(ref) + min(0, l)],
                        b[max(0, -l) : len(b) + min(0, -l)],
                    )
                    for l in lags
                ]
            )
            k = int(np.argmax(c))
            dxs.append((n - _subpixel(c, k)) * step)  # same sign convention
            quality.append(float(c[k]))
        dxs = np.array(dxs)
        dys = np.full(dxs.shape, np.nan)
    else:
        raise ValueError(f"method must be 'map' or 'cut', got {method!r}")

    dxs = np.where(np.abs(dxs) <= max_shift, dxs, 0.0)
    out, report = [], []
    for e, dx, dy, q in zip(series, dxs, dys, quality):
        new = dict(e)
        new["XG"] = e["XG"] - dx
        new["x_shift"] = float(dx)
        out.append(new)
        report.append(
            dict(
                label=e["label"],
                a=e["a"],
                b=e["b"],
                dx=float(dx),
                dy=float(dy),
                quality=float(q),
            )
        )
    return out, report


# ---------------------------------------------------------------------------
# Stitching adjacent tiles along x
# ---------------------------------------------------------------------------


def overlap_report(tile_a, tile_b, field="diff"):
    """Compare two x-adjacent tiles where their travel overlaps.

    The overlap is the one place the same sample is measured twice, so it is
    the only direct check that the two tiles agree -- on level, and on whether
    they are placed correctly relative to each other.

    Returns a dict with the overlap span, the median ``b - a`` offset, and the
    correlation between the two copies. A high correlation with a non-zero
    offset means the tiles are aligned but sit on different baselines; a low
    correlation means the x placement itself is suspect.
    """
    lo = max(tile_a["XG"].min(), tile_b["XG"].min())
    hi = min(tile_a["XG"].max(), tile_b["XG"].max())
    if hi <= lo:
        return dict(
            overlaps=False, span=0.0, n=0, offset=None, corr=None, lo=lo, hi=hi
        )

    ylo = max(tile_a["YG"].min(), tile_b["YG"].min())
    yhi = min(tile_a["YG"].max(), tile_b["YG"].max())
    xs = np.linspace(lo, hi, 60)
    ys = np.linspace(ylo, yhi, 120)
    X, Y = np.meshgrid(xs, ys)

    sampled = []
    for tile in (tile_a, tile_b):
        f = RegularGridInterpolator(
            (tile["YG"][:, 0], tile["XG"][0]),
            tile[field],
            bounds_error=False,
            fill_value=np.nan,
        )
        sampled.append(f(np.stack([Y, X], -1)))

    good = np.isfinite(sampled[0]) & np.isfinite(sampled[1])
    if good.sum() < 10:
        return dict(
            overlaps=True,
            span=hi - lo,
            n=int(good.sum()),
            offset=None,
            corr=None,
            lo=lo,
            hi=hi,
        )
    va, vb = sampled[0][good], sampled[1][good]
    return dict(
        overlaps=True,
        span=hi - lo,
        n=int(good.sum()),
        offset=float(np.median(vb - va)),
        corr=float(np.corrcoef(va, vb)[0, 1]),
        lo=lo,
        hi=hi,
    )


def stitch_x(
    tiles, fields=("diff", "mean"), step=None, seam=None, level_match=False
):
    """Join x-adjacent tiles into one map, using each point only once.

    Tiles are placed by their own interferometer coordinates -- nothing is
    fitted or slid -- and the duplicated strip where two tiles overlap is
    resolved by cutting at its midpoint, so each tile contributes the half
    nearer its own centre and no sample point is counted twice.

    ``seam`` overrides the cut position (a list of x values, one per gap).
    ``level_match`` shifts each tile after the first by the median difference
    measured in its overlap with the previous one; off by default, because a
    step at the seam is usually worth seeing rather than hiding.

    Returns ``(XG, YG, {field: array}, info)``.
    """
    tiles = sorted(tiles, key=lambda t: float(t["XG"].mean()))
    if step is None:
        step = float(np.median([np.diff(t["XG"][0]).mean() for t in tiles]))

    shifts = [0.0] * len(tiles)
    if level_match:
        for i in range(1, len(tiles)):
            rep = overlap_report(tiles[i - 1], tiles[i], fields[0])
            if rep["offset"] is not None:
                shifts[i] = shifts[i - 1] - rep["offset"]

    # Cut each overlap down the middle unless told otherwise.
    seams = (
        list(seam)
        if seam is not None
        else [
            0.5 * (tiles[i]["XG"].max() + tiles[i + 1]["XG"].min())
            for i in range(len(tiles) - 1)
        ]
    )

    x0 = min(t["XG"].min() for t in tiles)
    x1 = max(t["XG"].max() for t in tiles)
    ylo = max(t["YG"].min() for t in tiles)
    yhi = min(t["YG"].max() for t in tiles)
    xs = np.arange(x0, x1 + 0.5 * step, step)
    ys = np.linspace(ylo, yhi, max(t["YG"].shape[0] for t in tiles))
    XG, YG = np.meshgrid(xs, ys)

    # Each output column belongs to exactly one tile.
    edges = [-np.inf, *seams, np.inf]
    out = {f: np.full(XG.shape, np.nan) for f in fields}
    owned = []
    for i, tile in enumerate(tiles):
        col = (xs >= edges[i]) & (xs < edges[i + 1])
        owned.append(int(col.sum()))
        if not col.any():
            continue
        sub = np.stack([YG[:, col], XG[:, col]], -1)
        for f in fields:
            g = RegularGridInterpolator(
                (tile["YG"][:, 0], tile["XG"][0]),
                tile[f],
                bounds_error=False,
                fill_value=np.nan,
            )
            out[f][:, col] = g(sub) + (shifts[i] if f == fields[0] else 0.0)

    info = dict(
        step=step,
        seams=seams,
        columns=owned,
        shifts=shifts,
        x_range=(x0, x1),
        y_range=(ylo, yhi),
        overlaps=[
            overlap_report(tiles[i], tiles[i + 1], fields[0])
            for i in range(len(tiles) - 1)
        ],
    )
    return XG, YG, out, info


def interior_columns(x, lo_frac=1.0 / 3.0, hi_frac=2.0 / 3.0):
    """Columns in the central band of a scan -- by default its middle third.

    Used to set colour limits from the sample's interior rather than from its
    absorption edges, where a sub-pixel CR/CL position mismatch across a large
    intensity step makes a difference many times any real dichroism (8.5x on
    the #946-949 tiles). Scaling to the whole map therefore renders the
    interior flat, hiding the domains.

    The middle third is interior by construction: the sample is centred in the
    scan, which is what the short absorption alignment scan before each angle
    establishes. Taking a fixed fraction is deliberate -- it clips the same way
    on every scan, where locating the edges from the intensity shifts with the
    contrast and can creep back towards them on a weakly absorbing sample.

    Returns a boolean mask over ``x``.
    """
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return np.zeros(0, dtype=bool)
    lo, hi = float(np.min(x)), float(np.max(x))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.ones(x.size, dtype=bool)
    span = hi - lo
    mask = (x >= lo + lo_frac * span) & (x <= lo + hi_frac * span)
    return mask if mask.any() else np.ones(x.size, dtype=bool)


# ---------------------------------------------------------------------------
# Correcting the beam shift between the two polarisations
# ---------------------------------------------------------------------------
#
# Moving the phase plate moves the beam a little, and the beam also wanders
# between one scan line and the next. Where the sample is steeply absorbing --
# at a slanted edge, say -- a shift of a few tens of nm changes the measured
# transmission far more than the dichroism does, so subtracting CL from CR
# without aligning them first leaves a large false signal along every edge and
# a stripe on every line.
#
# Measured on #1374/1375: the *global* shift is only ~18 nm and correcting it
# alone changes nothing, but correcting each line separately drops the edge
# artefact from 0.048 to 0.004 and nearly halves the line-to-line striping.


def shift_rows(image, x, dx):
    """Shift each row of ``image`` along ``x`` by ``dx`` (scalar or per row).

    Positive ``dx`` moves the row's content towards +x. Points shifted in from
    beyond the measured range become NaN rather than being extrapolated.
    """
    image = np.asarray(image, dtype=float)
    x = np.asarray(x, dtype=float)
    dx = (
        np.full(image.shape[0], float(dx))
        if np.isscalar(dx)
        else np.asarray(dx, dtype=float)
    )
    out = np.full_like(image, np.nan)
    for j in range(image.shape[0]):
        step = dx[j] if np.isfinite(dx[j]) else 0.0
        out[j] = np.interp(x - step, x, image[j], left=np.nan, right=np.nan)
    return out


def line_shifts(ref, moving, x, max_shift=1500.0, min_points=40):
    """Per-row displacement of ``moving`` relative to ``ref``, in x.

    One normalised cross-correlation per row, refined by a parabola through
    the peak. Normalising matters: a plain dot product shrinks as the overlap
    shrinks with lag, so it always peaks at zero and reports no shift at all.

    Rows with too little overlap give NaN, which :func:`shift_rows` treats as
    "leave this row alone".
    """
    ref = np.asarray(ref, dtype=float)
    moving = np.asarray(moving, dtype=float)
    x = np.asarray(x, dtype=float)
    step = float(np.median(np.diff(x)))
    n = max(1, int(max_shift / step))
    lags = np.arange(-n, n + 1)

    out = np.full(ref.shape[0], np.nan)
    for j in range(ref.shape[0]):
        a, b = ref[j], moving[j]
        if (np.isfinite(a) & np.isfinite(b)).sum() < min_points:
            continue
        c = np.full(lags.size, np.nan)
        for i, lag in enumerate(lags):
            u = a[max(0, lag) : len(a) + min(0, lag)]
            v = b[max(0, -lag) : len(b) + min(0, -lag)]
            m = np.isfinite(u) & np.isfinite(v)
            if m.sum() < min_points:
                continue
            uu = u[m] - u[m].mean()
            vv = v[m] - v[m].mean()
            denom = uu.std() * vv.std()
            if denom > 0:
                c[i] = float((uu * vv).mean() / denom)
        if not np.isfinite(c).any():
            continue
        k = int(np.nanargmax(c))
        if (
            0 < k < c.size - 1
            and np.isfinite(c[k - 1])
            and np.isfinite(c[k + 1])
        ):
            den = c[k - 1] - 2 * c[k] + c[k + 1]
            if den != 0:
                k = k + 0.5 * (c[k - 1] - c[k + 1]) / den
        out[j] = (n - k) * step
    return out


def realign_pair(
    diff, mean, x, mode="line", manual=0.0, max_shift=1500.0, shifts=None
):
    """Re-derive a pair's difference with the two polarisations aligned in x.

    Takes the difference and mean a pair was already reduced to -- the two
    polarisations are recovered exactly as ``mean -/+ diff/2`` -- so this needs
    no re-reading of the scans and can be re-run interactively.

    ``mode`` is ``"none"``, ``"global"`` (one shift for the whole map, the
    median of the per-line measurements) or ``"line"``. ``manual`` adds a fixed
    offset in nm on top, for tuning by eye.

    Returns ``(diff, mean, shifts)``.
    """
    diff = np.asarray(diff, dtype=float)
    mean = np.asarray(mean, dtype=float)
    cr = mean - 0.5 * diff
    cl = mean + 0.5 * diff

    if shifts is None and mode != "none":
        shifts = line_shifts(cr, cl, x, max_shift=max_shift)
    if mode == "none":
        shifts = np.zeros(diff.shape[0])
    else:
        shifts = np.asarray(shifts, dtype=float)
        if mode == "global":
            median = np.nanmedian(shifts)
            shifts = np.full(diff.shape[0], 0.0 if np.isnan(median) else median)
        elif mode != "line":
            raise ValueError(f"mode must be none/global/line, got {mode!r}")

    total = np.where(np.isfinite(shifts), shifts, 0.0) + float(manual)
    if not np.any(total):
        return diff, mean, shifts
    cl = shift_rows(cl, x, -total)
    return cl - cr, 0.5 * (cl + cr), shifts


def tile_pairs(pairs, roi, folder, grid_n=200, open_pct=99.0, open_level=None):
    """Reduce several CR/CL pairs against **one** open-beam reference.

    Referencing each scan to its own high percentile treats the datasets
    differently, in two ways that both corrupt a stitch:

    * Within a pair, the difference picks up ``ln(ref_CL / ref_CR)`` -- a
      constant that has nothing to do with dichroism and differs from tile to
      tile. On #1959-1962 that alone put the two tiles at -0.0053 and +0.0007.
    * Across tiles, a scan that never runs off the sample references itself to
      its own thinnest point and reports every absorbance too low by a
      constant. #1959 peaked at 2648 where #1961 plateaued at 3651, so tile A
      came out 0.349 low throughout -- about 12%.

    Transmission cannot exceed the unobstructed beam, so the largest
    percentile across all the scans is the best available estimate of it, and
    using that one value for every scan makes the tiles directly comparable:
    the reference cancels exactly in every difference, and absorbances share a
    zero. Pass ``open_level`` to override it.

    Each scan is read once, however many pairs mention it. Returns
    ``(tiles, open_level)``; each tile carries ``XG``, ``YG``, ``diff``,
    ``mean``, ``a``, ``b``.
    """
    raw = {}
    for a, b in pairs:
        for scan in (a, b):
            if scan not in raw:
                raw[scan] = normalised_map(scan, roi, folder)

    if open_level is None:
        levels = [np.nanpercentile(r, open_pct) for _x, _y, r in raw.values()]
        levels = [v for v in levels if np.isfinite(v) and v > 0]
        if not levels:
            raise ValueError("no usable open-beam level in any scan")
        open_level = float(max(levels))

    tiles = []
    for a, b in pairs:
        xa, ya, ra = raw[a]
        xb, yb, rb = raw[b]
        with np.errstate(divide="ignore", invalid="ignore"):
            # Same convention as absorbance_map: a non-positive ratio has no
            # logarithm and becomes NaN rather than being clipped to a floor.
            za = -np.log(np.where(ra > 0, ra, np.nan) / open_level)
            zb = -np.log(np.where(rb > 0, rb, np.nan) / open_level)
        xg = np.linspace(
            max(xa.min(), xb.min()), min(xa.max(), xb.max()), grid_n
        )
        yg = np.linspace(
            max(ya.min(), yb.min()), min(ya.max(), yb.max()), grid_n
        )
        XG, YG = np.meshgrid(xg, yg)
        ga = griddata((xa, ya), za, (XG, YG), method="linear")
        gb = griddata((xb, yb), zb, (XG, YG), method="linear")
        tiles.append(
            dict(XG=XG, YG=YG, diff=gb - ga, mean=0.5 * (ga + gb), a=a, b=b)
        )
    return tiles, float(open_level)


def common_y_range(first_scans, folder, group=4):
    """The y interferometer range shared by several angle groups.

    The stage drifts in y over a long tomographic series -- about 1.56 µm
    between #1716 (25°) and #2479 (203.8°) -- so the slice measured at one
    angle is not the slice measured at another. Only the overlap of every
    angle can be reconstructed; outside it, different angles are looking at
    different parts of the sample.

    Reads position streams only, never detector frames, so it is cheap. Pass
    the *first* scan of each group; ``group`` scans are examined from each.
    Returns ``(lo, hi)``, or ``None`` if no scan could be read.
    """
    lo, hi = -np.inf, np.inf
    seen = False
    for first in first_scans:
        for scan in range(int(first), int(first) + group):
            try:
                ys = reduce_position_stream(scan, folder).ys
            except Exception:
                continue  # a missing or aborted scan just drops out
            if ys.size == 0:
                continue
            lo = max(lo, float(ys.min()))
            hi = min(hi, float(ys.max()))
            seen = True
    if not seen or not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return None
    return lo, hi
