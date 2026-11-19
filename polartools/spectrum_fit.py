"""Gaussian peak fitting for the spectrum viewer's trend.

The trend in :mod:`spectrum_viewer` totals a box ROI per frame. This module
offers the alternative of fitting several Gaussians and a background over a
window of bins, and taking the *area* of the chosen peaks instead -- for lines
that overlap, or sit on a background a box cannot take out.

**The shape is fitted once, the areas are a projection.** The centres and
widths are fitted (with lmfit) to a summed spectrum -- many frames, both
halves of a contrast -- and then held fixed. With the shape fixed, the
least-squares areas of any spectrum ``y`` over the same window are linear in
``y``: ``areas = pinv(G) @ y``, where ``G`` holds one unit-area Gaussian per
peak and the background's columns. :func:`area_weights` returns the row of
``pinv(G)`` for the peaks wanted, so that one frame's fitted area is a dot
product with its spectrum.

That linearity is the whole design. The area of a sum of frames is the sum
of the frames' areas, exactly as for an ROI sum, so the viewer can compute a
per-frame area and hand it to the same machinery that totals ROI counts over
points, halves and monitor-normalised ratios. A free fit at every point
would follow a line that moves, and would lose that: the halves and the
monitors could no longer share the reduction, and it is three orders of
magnitude slower. Measured on a 2026-3 dichro scan's Pt line, the fixed-shape
areas agree with free per-frame fits to 0.09 % (median), and cost 0.07 ms
for 312 frames against 140 ms. Shape-from-the-sum is also the physics a
contrast wants: polarization changes intensities, not the detector's line
shape. It is the wrong tool for the elastic and Compton lines of an energy
scan, which move with the incident energy.

**The background costs noise.** Its columns are estimated from the wings of
the window, which hold few counts. On that same Pt line, an isolated one,
the point-to-point noise of the XMCD/XAS ratio is 9.0e-3 from an ROI sum,
9.9e-3 from a fitted area with no background, and 1.26e-2 with a constant or
linear one. Fitting pays where peaks overlap or the background slopes; on a
clean line a box ROI is as good.

The least-squares projection is unweighted on purpose. Poisson weights
would be better per spectrum, but they depend on ``y``, and the projection
then stops being linear -- the property everything above rests on. The
shape fit itself is Poisson-weighted, since it is done once on a spectrum of
many counts.

No Qt here, so that all of it can be tested on synthetic spectra.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from lmfit.models import ConstantModel, GaussianModel, LinearModel
from scipy.signal import find_peaks, savgol_filter

SQRT_2PI = np.sqrt(2.0 * np.pi)
FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))

BACKGROUNDS = ("none", "constant", "linear")

# Bounds on a fitted width, in bins: below half a bin a Gaussian is a spike
# the bins cannot describe, and the upper bound is set per window.
MIN_SIGMA_BINS = 0.5
# How far a centre may move from where it started: this many starting
# widths, and at least CENTRE_SLACK_BINS. Unbounded, a weak peak can walk
# under a strong one, and two Gaussians on one line split its area at random
# -- seen on a softglue scan, where a third peak started at 7.9 keV ended at
# 7.35, under the 7.41 keV line.
CENTRE_SLACK_SIGMAS = 3.0
CENTRE_SLACK_BINS = 5.0
# How wide a peak may grow, as a multiple of its starting width. Unbounded,
# a Gaussian left with nothing to fit spreads out to soak up a neighbour's
# tail or the background -- seen growing from 12 to 1566 eV FWHM.
MAX_SIGMA_GROWTH = 4.0
# The narrowest a found peak may be, as a fraction of the most prominent
# one's width: every line is broadened by the same detector, so a curvature
# bump far narrower than the main line is noise, not a line.
MIN_RELATIVE_WIDTH = 1.0 / 3.0


@dataclass
class PeakSpec:
    """One peak's starting values, in bins."""

    centre: float
    sigma: float
    vary_centre: bool = True


@dataclass
class ShapeFit:
    """A fitted shape over bins ``lo..hi`` inclusive.

    *centres* and *sigmas* are in bins; *areas* are the fit's own, on the
    spectrum it was made from. ``success`` False means the fit raised or did
    not converge, and *message* says why.
    """

    lo: int
    hi: int
    background: str
    centres: np.ndarray = field(default_factory=lambda: np.empty(0))
    sigmas: np.ndarray = field(default_factory=lambda: np.empty(0))
    areas: np.ndarray = field(default_factory=lambda: np.empty(0))
    redchi: float = float("nan")
    success: bool = False
    message: str = ""

    @property
    def n_peaks(self) -> int:
        return int(self.centres.size)


def unit_gaussian(x: np.ndarray, centre: float, sigma: float) -> np.ndarray:
    """A Gaussian of unit area, so a coefficient on it is an area."""
    return np.exp(-0.5 * ((x - centre) / sigma) ** 2) / (sigma * SQRT_2PI)


def background_columns(x: np.ndarray, background: str) -> list[np.ndarray]:
    """The background's columns of the design matrix.

    The slope column is centred and scaled to the window, so that it is not
    nearly parallel to the constant one and ``pinv`` stays well conditioned.
    """
    if background == "none":
        return []
    columns = [np.ones_like(x)]
    if background == "linear":
        half = max((x[-1] - x[0]) / 2.0, 1.0)
        columns.append((x - (x[0] + x[-1]) / 2.0) / half)
    elif background != "constant":
        raise ValueError(f"no such background {background!r}")
    return columns


def design_matrix(
    lo: int, hi: int, centres, sigmas, background: str
) -> np.ndarray:
    """``(hi - lo + 1, n_peaks + n_background)``: peaks first, then background."""
    x = np.arange(lo, hi + 1, dtype=float)
    columns = [unit_gaussian(x, c, s) for c, s in zip(centres, sigmas)]
    columns += background_columns(x, background)
    return np.stack(columns, axis=1)


def fit_shape(
    y: np.ndarray,
    lo: int,
    hi: int,
    peaks: list[PeakSpec],
    background: str = "linear",
    shared_width: bool = False,
) -> ShapeFit:
    """Fit centres, widths and areas of *peaks* to ``y[lo:hi + 1]``.

    *y* is a whole spectrum, indexed by bin. Poisson-weighted, with the
    centres kept inside the window and near where they started
    (:data:`CENTRE_SLACK_SIGMAS`), and, with *shared_width*, every width
    tied to the first -- the detector's resolution changes little across a
    few hundred eV. Never raises: a failure comes back as ``success`` False.
    """
    result = ShapeFit(lo=int(lo), hi=int(hi), background=background)
    if not peaks:
        result.message = "no peaks to fit"
        return result
    if hi - lo + 1 < 3 * len(peaks) + 2:
        result.message = "window too narrow for that many peaks"
        return result

    x = np.arange(lo, hi + 1, dtype=float)
    data = np.asarray(y[lo : hi + 1], dtype=float)
    if not np.all(np.isfinite(data)) or data.sum() <= 0:
        result.message = "no counts in the window"
        return result
    weights = 1.0 / np.sqrt(np.maximum(data, 1.0))

    model = None
    for i in range(len(peaks)):
        part = GaussianModel(prefix=f"p{i}_")
        model = part if model is None else model + part
    if background == "constant":
        model = model + ConstantModel(prefix="bg_")
    elif background == "linear":
        model = model + LinearModel(prefix="bg_")
    elif background != "none":
        raise ValueError(f"no such background {background!r}")

    params = model.make_params()
    floor = float(np.percentile(data, 10)) if background != "none" else 0.0
    max_sigma = max(float(hi - lo), MIN_SIGMA_BINS * 2)
    for i, peak in enumerate(peaks):
        centre = float(np.clip(peak.centre, lo, hi))
        sigma = float(np.clip(peak.sigma, MIN_SIGMA_BINS, max_sigma))
        widest = min(max_sigma, max(MAX_SIGMA_GROWTH * sigma, 2.0))
        height = max(float(data[int(round(centre)) - lo]) - floor, 1.0)
        slack = max(CENTRE_SLACK_SIGMAS * sigma, CENTRE_SLACK_BINS)
        params[f"p{i}_center"].set(
            value=centre,
            min=max(lo, centre - slack),
            max=min(hi, centre + slack),
            vary=peak.vary_centre,
        )
        params[f"p{i}_sigma"].set(value=sigma, min=MIN_SIGMA_BINS, max=widest)
        params[f"p{i}_amplitude"].set(value=height * sigma * SQRT_2PI, min=0.0)
        if shared_width and i > 0:
            params[f"p{i}_sigma"].set(expr="p0_sigma")
    if background == "constant":
        params["bg_c"].set(value=floor)
    elif background == "linear":
        params["bg_intercept"].set(value=floor)
        params["bg_slope"].set(value=0.0)

    try:
        fit = model.fit(data, params, x=x, weights=weights)
    except Exception as exc:  # lmfit raises a variety; all mean "no fit"
        result.message = f"fit failed: {exc}"
        return result

    n = len(peaks)
    result.centres = np.array(
        [fit.params[f"p{i}_center"].value for i in range(n)]
    )
    result.sigmas = np.array(
        [fit.params[f"p{i}_sigma"].value for i in range(n)]
    )
    result.areas = np.array(
        [fit.params[f"p{i}_amplitude"].value for i in range(n)]
    )
    result.redchi = float(fit.redchi)
    # lmfit's own success flag also drops for "could not estimate
    # error-bars", which says nothing about the parameters, and the areas
    # use no error bars. So a fit is good when it ran to finite numbers;
    # lmfit's message is kept as a warning.
    result.success = bool(
        np.all(np.isfinite(result.centres))
        and np.all(np.isfinite(result.sigmas))
        and np.isfinite(result.redchi)
    )
    result.message = "" if fit.success else (fit.message or "")
    if not result.success:
        result.message = result.message or "the fit gave non-finite values"
    return result


def area_weights(shape: ShapeFit, used) -> np.ndarray:
    """Weights over bins ``lo..hi`` whose dot product with a spectrum is its area.

    The summed area of the peaks indexed by *used*, solved by unweighted
    least squares with *shape*'s centres and widths held fixed. Linear in the
    spectrum, so totals over frames commute with it -- see the module
    docstring. Zeros when nothing is used.
    """
    width = shape.hi - shape.lo + 1
    used = [i for i in used if 0 <= i < shape.n_peaks]
    if not used:
        return np.zeros(width)
    projector = np.linalg.pinv(
        design_matrix(
            shape.lo, shape.hi, shape.centres, shape.sigmas, shape.background
        )
    )
    return projector[used].sum(axis=0)


def decompose(y: np.ndarray, shape: ShapeFit) -> tuple[np.ndarray, np.ndarray]:
    """``(coefficients, matrix)``: *y*'s fixed-shape least-squares solve.

    The same projection the trend uses, so the curves drawn from it and the
    areas listed beside them are the numbers the trend sums. Coefficients
    are the peaks' areas first, then the background's.
    """
    matrix = design_matrix(
        shape.lo, shape.hi, shape.centres, shape.sigmas, shape.background
    )
    data = np.asarray(y[shape.lo : shape.hi + 1], dtype=float)
    return np.linalg.pinv(matrix) @ data, matrix


def suggest_peaks(
    y: np.ndarray,
    lo: int,
    hi: int,
    max_peaks: int = 6,
    prominence: float = 0.05,
    smooth: int = 7,
) -> list[PeakSpec]:
    """Starting peaks for ``y[lo:hi + 1]``, the most prominent first.

    Found as peaks of the curvature ``-y''`` rather than of ``y``, because
    the lines worth fitting are the ones that overlap, and a weak line on a
    strong one's tail is a shoulder, not a maximum: ``find_peaks`` on the
    counts misses it outright. Every Gaussian, shoulder or not, is a
    maximum of ``-y''``. The derivative is taken through a Savitzky-Golay
    filter *smooth* bins wide, since a raw second difference of counts is
    mostly noise.

    *prominence* is a fraction of the largest curvature in the window, and
    a candidate narrower than :data:`MIN_RELATIVE_WIDTH` of the most
    prominent peak is dropped as a ripple.
    A Gaussian's curvature changes sign at one sigma either side of its
    centre, so half the run over which ``-y''`` stays positive is the
    starting width -- overlapping lines narrow it, and the fit corrects it.
    """
    data = np.asarray(y[lo : hi + 1], dtype=float)
    smooth = max(5, smooth | 1)
    if data.size < smooth or not np.any(data > 0):
        return []
    curvature = -savgol_filter(data, smooth, 3, deriv=2)
    top = curvature.max()
    if top <= 0:
        return []
    found, props = find_peaks(curvature, prominence=prominence * top)
    if not found.size:
        return []
    order = np.argsort(props["prominences"])[::-1]

    def width(p):
        left, right = p, p
        while left > 0 and curvature[left - 1] > 0:
            left -= 1
        while right < curvature.size - 1 and curvature[right + 1] > 0:
            right += 1
        return max((right - left + 1) / 2.0, MIN_SIGMA_BINS)

    sigmas = {int(p): width(p) for p in found[order]}
    floor = MIN_RELATIVE_WIDTH * sigmas[int(found[order[0]])]
    kept = [int(p) for p in found[order] if sigmas[int(p)] >= floor]
    return [
        PeakSpec(centre=float(lo + p), sigma=float(sigmas[p]))
        for p in sorted(kept[:max_peaks])
    ]
