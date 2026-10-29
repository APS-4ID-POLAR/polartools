"""Browse a directory of Bluesky scans and plot them.

For each ``*_master.hdf`` file written by ``MyNXWriter``
(``callbacks/nexus_data_file_writer.py``), this shows:

- a scan listing on the left -- one row per scan giving its number, plan,
  motor, range, points and count time -- filled in as the directory is listed,
  without having to select anything. Only files carrying Bluesky scan metadata
  are listed; the detector image files, whether beside the masters or in a
  per-detector subdirectory, are not. Every column can be
  dragged to a different width. The leading part of whatever is selected -- the
  name up to the scan number, given once however many scans of that series are
  picked -- is shown above the plot. An "Img" tick marks the
  scans whose detector images were found on disk. Click (or
  ctrl/shift-click for several) to open one or more scans; Refresh picks up
  new/changed files without losing the current selection,
- one of the scan's ``entry/instrument/bluesky/streams``, chosen from the
  Stream box, with a field table under it whose X / Y / Norm ticks drive the
  plot -- X and Norm take a single field each, Y any number. An
  ``EPOCH`` entry (seconds past the UNIX epoch, one per point) joins the real
  fields wherever the file records it. Most scans record only ``primary`` and
  the box has nothing to choose between; a dichro scan also writes the XAS and
  XMCD it worked out to ``dichro_monitor``, and is shown on that stream, since
  that reduced result is what it was run for. Picking a stream by hand keeps it
  across scans that have it. X and Y start on what the scan's own ``hints``
  metadata nominates -- the axis it scanned for X (the ``EPOCH`` column for a
  plan that moved nothing) and the one signal per detector worth looking at
  for Y -- so a scan plots itself the way it was meant to be read as soon as
  it is selected. Norm starts clear: the hints name the monitor the scan
  counted against, but normalizing rescales every number on the plot and in
  the statistics table, so a scan opens showing what it recorded and the tick
  is yours to make. The scan
  that leads the selection decides, and decides again whenever the kind of
  scan changes -- the columns on offer, or which of them the scan asks for.
  Step along a run of lineup scans and the ticks made on the first stand for
  the rest, which is the comparison they were made for; open something else
  and all three are read back from that scan's own hints, so nothing you
  ticked an hour ago, on a scan since forgotten, is still deciding what you
  are looking at now. (``hints``
  describes ``primary``; a reduced stream brings its own defaults instead.
  The stream itself is the exception, and is kept when picked by hand.)
  Curves are overlaid across every selected file, optionally normalized by
  dividing each Y field by the monitor field -- or, with "log(mon/det)"
  ticked, by taking log(Norm / Y) instead, which is the absorption a
  measurement through the sample asks for and is what makes a transmission
  XAS scan readable straight off ``primary`` rather than only through the
  reduced panels below. The two forms are the same button in both places,
  so a curve and a panel of one measurement cannot come out in different
  units; it starts unticked, and a plain ratio is what a Norm tick gives
  until it is set. Pressed with no Norm ticked it ticks the monitor the scan
  hints, which is the one place that hint is read -- the Norm column itself
  still starts clear, so a scan opens showing what it recorded, but the
  button a transmission measurement wants is one press rather than two.
  A colour map normalizes the same way. A field a stream keeps on its
  own scale is drawn against a right-hand Y axis coloured to match -- XMCD,
  being a difference often a percent of the XAS it rides on, would otherwise
  be flattened into the axis. "Right scale" sends whole scans there by hand,
  for the other reason two curves fail to compare: two runs of one channel
  offset from each other by more than either one's own swing sit at two
  heights in a frame stretched to hold both, so their shapes -- which is what
  the comparison is about -- are read nowhere near each other. Given its own
  scale each fills the frame and the two superimpose. Sending every scan
  there would leave nothing to compare against, so that does nothing. A
  magneto-chiral scan (``macros_other_scans.py``) steps a second positioner -- the magnet current
  -- inside the scanned one, so its reduced stream holds one XAS/XMCD pair per
  step, which drawn as a single curve runs back and forth across the plot. It
  is split into one curve per step instead, and with more than one Y field
  each gets a plot of its own, stacked above the next and sharing the X axis;
  the "Average repeats" button gathers the steps that sat at one field (a
  pattern of [1, -1, -1, 1] measures each twice) into a curve each. A
  fixed-polarization scan has no reduced stream at all -- it steps the
  polarization and the magnet as two grid axes of its own, so ``primary``
  holds every setting itself -- and is split the same way, along both axes at
  once, from the raster the scan declares. Declares and the data bears out: a
  magneto-chiral scan's ``primary`` carries sub-measurements its declared
  shape does not describe, so a period taken there on trust would draw two
  polarizations as one curve. Either way those curves are the raw absorption,
  and what the measurement was for is the four contrasts the "XMCD panels"
  button reduces them to: XAS, XMCD, XMuCD and XNCD against energy, four
  panels stacked, every selected scan overlaid in its own colour. The points
  taken at one energy are grouped into the four (polarization, field)
  quadrants, which separates the terms of ``mu(p, B) = A + C*p + D*B + E*p*B``
  -- the grouping, the tolerances and the four combinations are those of the
  beamline's own ``xmucd_reduce.py``, reproduced here rather than imported
  because it lives in a per-cycle data directory. XNCD is the field-even
  combination and should sit at the noise floor, so the status bar reports
  its rms beside XMuCD's: when it does not, the polarization average is not
  cancelling XMCD and the magneto-chiral signal is contaminated. The reduction
  reads ``primary`` whatever the Stream box shows -- a dichro scan opens on
  its reduced stream, which records no retarder column -- and the axis label
  names the stream it used. The Detector and Monitor boxes beside the button
  say which pair to divide, a choice the Y and Norm ticks cannot stand in for
  since a scan's hinted detector is not always one that was live, and
  "log(mon/det)" switches between a transmission measurement and a plain
  ratio -- the same button, and the same state, that the Norm tick above
  uses, so the panels open on a plain det/mon and a transmission run is
  ticked into the log once for both. "Remove spikes" drops the points the beam was gone for, in the
  panels and in the curves alike. The monitor reads a fraction of a percent
  of its usual at such a point and every channel goes with it, so it holds no
  measurement -- yet it plots like any other, and reduced it is worse than
  that, since the reduction averages the raw points sharing an energy and one
  dead point of four becomes the whole energy's answer: on scan 239 one such
  point stretches the XMuCD panel from 0.026 to 0.511, twenty times the real
  figure and the largest feature on it. The test is the monitor against its
  own scan's median, which is not a judgement call -- over the 988 scans here
  that carry one, the worst live point sits at 0.84 of its median and the
  best dead one at 0.02, so any cut inside that fortyfold gap draws the same
  line. Only lone points go: where the beam was down for a stretch, as in
  scan 244's 265 consecutive points, that is the run and not a spike, and the
  scan goes on showing it. What goes leaves a gap rather than a guess,
  nothing being interpolated into a place nothing was measured, and because
  the points are dropped before the reduction rather than after it the panels
  keep their energies -- the affected bin simply averages the points that did
  have beam. Reading the monitor is also what keeps a real one-point feature:
  index 80 of scans 239 and 242 is a thirteen-sigma excursion at one energy,
  but its monitor is healthy and all four of its raw points agree, which
  makes it a monochromator glitch and the filter leaves it alone. A colour
  map is left out, a raster's wild cells being what "Clip map scale" below is
  for. A two-motor grid scan is instead drawn as a colour map over the two
  motors, one per ticked Y field, side by side on shared axes and each with
  its own colour scale. Several selected grid scans over the same two motors
  go on one map: scan 39 of Ryan_26-3/sro was stopped at h = 0.985 and scan
  41 started again there, so apart they are two halves of one picture. A
  cell more than one of them measured -- the 52 of that shared row -- is
  their average, and the heading says how many there were; scans over other
  motors are left out with a note. A dichro mesh
  scan's XAS and XMCD cover the same area but not remotely the same range, so
  one shared scale would wash the XMCD out. Those two are mapped off the
  scan's reduced stream, which holds one point per mesh position: primary has
  both polarizations at every position, so a map of it can only show one of
  them. The reduced stream carries averaged readbacks and no setpoint
  alongside, so where they wander off the raster -- scan 280's by a part in
  two hundred of a step -- they are gathered back onto the one the scan
  declares. (A grid scan declares itself in
  its ``hints``, which is also what keeps a magneto-chiral scan -- two
  positioners as well, but the inner one a setting rather than a coordinate --
  on the split curves above.) Where one or two wild cells -- a detector
  glitch, a corner where the beam ran off the sample -- stretch a scale over a
  range nothing else occupies and leave the rest of the map one flat shade,
  the "Clip map scale" button keeps the most extreme cells at each end out of
  the scale: they saturate, the colour bar grows an arrow at each end it cut,
  and the statistics below still report the true extremes. The box beside sets
  how much goes, starting at half a percent -- a spike or two needs no more
  than that, a map still washed out can take a good deal more -- and each map
  is cut by its own cells, so clipping a spiky one leaves the one beside it
  alone. The amount is asked for rather than worked out because the rules that
  judge outliers for themselves mistake a lineup peak for one. The colour box
  beside that picks the map they are drawn in; left on "Automatic" each panel
  chooses for itself, which means a diverging seismic for a signed channel
  like a mesh scan's XMCD -- where zero is a value the sample has and which
  side of it a cell falls on is the result -- and an ordinary viridis for
  everything else, XAS included, that merely runs from less to more. A
  diverging map is read by where its middle colour falls, so "Center on 0"
  widens such a scale to whichever end reaches further and puts white at
  zero; fitted to the data instead, white lands wherever the range happens to
  centre, drawing signal as though there were none. It applies only to the
  maps actually drawn in a diverging scale, so the ordinary ones are left as
  they are. The
  "2D map"
  button turns the map off to plot the raster as an ordinary X-Y curve, zoomable by dragging (via the
  plot toolbar's rectangle-zoom tool; "Log X" and "Log Y" beside it put
  either axis on a log scale, leaving out what is at or below zero), with a
  peak-statistics table (max,
  min, half-max center and FWHM) below the plot for each curve shown, and
  a "Show dY/dX" toggle that overlays each curve's numerical derivative on
  a further Y axis, stepped outward when a right-hand axis is already there.
  That toggle gives the table a second row per curve, directly under the
  first, with the same statistics taken on the derivative -- which is what
  carries the answer on a knife edge or an absorption edge, where the curve
  rises to a plateau and so has no peak of its own to measure,
- an "Integrate" toggle giving the table two further columns, the
  trapezoidal area under each curve and the same above a straight line
  joining its first and last points. On repeats of one rocking curve the
  area compares them where the peak height does not -- scans 32 to 38 of
  TiSe2 spread 21% in height against 15% in area -- and the gap between
  the two columns is the sloping background, which there is 9 to 15% and
  a different share in every scan. It works on the reduced panels too,
  where the integral of XMCD over the edge is the sum-rule quantity,
- a "Fit" toggle with a shape beside it -- Gauss, Lorentzian, Lorentzian
  squared or pseudo-Voigt, the last a Gauss and a Lorentzian blended at
  one shared width in a proportion fitted along with the rest and
  reported as "Fit eta", 1 the pure Lorentzian and 0 the pure Gauss --
  fitting ``peak + m*x + c`` to every curve on screen, drawing
  the model over the data in the curve's own colour and giving the table
  five further columns: the fitted centre, FWHM, height and area, and an
  R2, with eta a sixth for the one shape that has one.
  Where the "Center"/"FWHM" columns interpolate half-max crossings between
  two points, the fit uses every point at once, which is what the repeats
  ask for: over scans 32 to 38 the crossing width spreads 78% of its
  smallest value and the fitted one 62%. "Fit area" is likewise the
  integrated intensity of the peak alone -- height x FWHM x a constant of
  the shape, exact over the whole line rather than a trapezoid over the
  points -- so it is not the "Net area" of "Integrate" and does not match
  it (2219 against 2791 on scan 32). It is reckoned over the whole line
  on purpose, so that the same peak swept wider reports the same
  intensity; the Lorentzian pays for that with 6 to 11% of its area in
  tails past the ends of the scan, the pseudo-Voigt with whatever share
  of that its fitted eta asks for -- 4 to 7% at a half -- where
  the Gaussian has nothing worth printing and its square under 1%. The
  background is always a sloping line because on those scans it
  genuinely slopes -- 1754 counts at the left edge against 476 at the
  right on scan 32 -- and fitting it
  lifts the mean R2 of the whole curve from 0.957 to 0.980 for a Gauss
  (0.976 for a flat offset, so it is the slope and not merely a pedestal
  that is being bought). The shape barely moves the centre (0.0005 degrees
  across all four on scan 32) and moves the width a great deal (0.039 to
  0.045), so it is a width control more than a position one. Read R2
  first, and note that it is reckoned on the peak alone with the fitted
  line taken off both sides: against the raw curve the line explains most
  of the variance by itself and pure noise scores 0.6 to 0.9. A peak model
  fitted to an absorption edge returns a narrow peak and a centre to five
  figures with no complaint -- scan 171 does, for all four shapes -- and
  only the R2, 0.11 there, says the model describes nothing. A fit further
  gone than that is declined rather than drawn, and declining is common:
  of the 625 curves the two test directories here plot, 263 blank, nearly
  all of them steps, plateaus, monotonic rises and scans too short for
  five parameters. A peak the scan stepped over and caught in a single
  point is declined too (scan 25), since one sample fixes no width.
  While "Show dY/dX" is on the fit moves onto the derivative instead of
  the curve, which is what makes a step measurable: a knife edge or an
  absorption edge holds no peak and is declined, while its derivative is
  the beam profile or the edge itself. Scans 101 to 106 of the tests
  directory fit at R2 0.99 that way and not at all on the curve, and 64
  of the 613 curves here fit only in that direction. Going the other way
  is a trap the R2 catches: the derivative of a peak is a lobe up then a
  lobe down, so ``peak + line`` takes one lobe, displaces the centre and
  flips the height's sign between identical repeats (scans 32 to 38, R2
  0.52 to 0.71). Differentiating amplifies noise too, so a derivative
  ragged against its own step size is declined -- scan 171's is,
- a "Multi" toggle beside "Fit", and a count beside that, fitting several
  peaks at once over one shared sloping line, three parameters apiece.
  The three cases it was built for are all in TiSe2: two peaks barely
  overlapping (scan 109, which one peak will not fit at all and two fit
  at R2 0.92), two overlapping badly (scan 34, whose second peak at
  44.7406 is a shoulder on the flank of the first and is no local
  maximum, 0.982 to 0.992), and the derivative of a peak (scan 31 with
  dY/dX, 0.641 to 0.990) -- which is also the way out of the one-lobe
  trap above, since two lobes are two peaks. Seeded from
  ``scipy.signal.find_peaks`` run over the data and its negative, so a
  dip seeds as readily as a rise, and any shortfall taken from the
  largest residual, which is what finds scan 34's shoulder. Each peak
  gets a row of its own, numbered from the left, and the curve's own row
  blanks: a sum of peaks has no one centre to report. The count has to be
  the user's, because no finder supplies it -- one candidate on scan 34
  where there are two, six on scan 109 where there are two -- and R2
  cannot supply it either, the models being nested so that it only ever
  rises (scan 32, 0.978 on one peak and 0.998 on two). Nor can the gates:
  forced to two, a single peak is divided rather than refused, scan 32
  coming back as a broad base of FWHM 0.0914 under a sharp core of
  0.0347 where the peak is really 0.0447. Off by default, and meant to be
  switched on for a curve that needs it rather than left on,
- a "Parts" toggle beside the count, drawing what the fit is made of as
  well as the fit: the sloping background in dash-dot, and with "Multi"
  on each component in dots, every one of them on the background rather
  than from zero, so a component lies over the stretch of curve it is
  claiming. Drawn from zero they would be nowhere near the data wherever
  the background is not near it either -- scan 109's is some 800 counts.
  This is the answer to the one thing R2 and the gates cannot tell you,
  that a component is not really there: a divided single peak shows as
  two dotted lines over the same bump (scan 32 at two peaks, the broad
  base and the sharp core), where a real pair shows one line over each
  (scan 109). A single peak adds only the background, its one component
  over that background being the fitted line already drawn. Nothing is
  exported or tabulated that was not already: the pieces are the fit the
  five columns and the header name, split up,
- a "Plot" toggle beside "Fit", opening a second plot to the right of the
  first that trends one fitted quantity -- height, FWHM, area, centre, or
  a pseudo-Voigt's eta -- selected in the combo
  next to it, against the scan number or against a baseline reading the
  selected scans share, chosen in the combo after it. It needs "Fit", and
  greys out without it, because fitted quantities are all it
  draws -- never the crossing width or the trapezoid area, which would
  quietly mix two different measurements into one line. The baseline axis
  is what turns a row of repeats into the curve they were measured to
  draw: scans 32 to 39 step huber_euler_ana_eta from 0 to 90 degrees, and
  the fitted area against eta is that sweep rather than eight points in
  scan order. Only fields that moved are offered, most-moved first, judged
  against how far each drifted inside a scan -- otherwise the list is led
  by beam diagnostics that wander further than anything deliberate and say
  nothing. One series per Y field, taken from
  whatever the fit ran on, so with dY/dX shown it trends the derivative's
  fit and the axis label says so, and one per peak as well while "Multi"
  is on -- ordered by centre, which is what lets a line follow the same
  physical peak through a run of repeats. The ten knife-edge repeats at
  scans 101 to 127 are what it is for: their beam size against scan
  number, 0.0198 to 0.0217 over the run. Changing the peak count partway
  through such a run moves every value by more than that drift (scan 101
  reads 0.0210 at one peak and 0.0234 at two), so a series has to be made
  at one setting. Curves with no scan number of their own --
  averages and combinations, which keep only their first member's path --
  are left out and the status bar says how many,
- an "Export…" at the end of that row writing the trend to a file:
  one row per scan, the chosen axis in the first column and then every
  fitted quantity rather than the one being plotted, since a peak's
  numbers are read together and exporting under each in turn would make
  four files to join back up. One set of columns per series, named as the
  trend names its lines, so "Multi" gives each peak its own. A declined
  fit writes ``nan``, which is the gap in the plot beside it. With
  "Integrate" on the curve's own two trapezoids come as well, as
  ``curve_area``/``curve_net_area`` -- one pair per Y field and not per
  peak, since they measure a curve rather than a peak of it. The header
  is the one the "Export…" at the top writes, that one holding the curves
  themselves,
- an "Invert" button turning whole scans upside down, ``y -> -y``, picked
  from a menu of the selected scans exactly as "Right scale" is. A channel
  read at the opposite lock-in phase comes back as a dip where its partner
  is a peak, which neither compares by eye nor measures: the statistics
  table finds a "Center" and an "FWHM" on a peak and leaves them blank on a
  dip. Negated, the two lie on top of each other and both are measured. The
  negation also throws the scan clear of the one it is to be compared with
  -- 212 and 213 both read LockAC at 2.6e6, so inverting one spreads the
  shared scale over 5.2e6 and leaves each 4.2e4 swing at 0.8% of it -- so
  while nothing else holds the right-hand scale the inverted scan takes it,
- an "Add" button drawing the sum of two selected scans, or their
  difference either way round, as a third curve of its own beside the two it
  came from -- a repeat summed for counting statistics, a pair differenced
  to see what moved between them. The second scan is read onto the first
  one's X by interpolation rather than added point by point: scans 171 and
  173 are 109 energies each and differ by up to 4e-05 keV, the readback
  sampled at two different moments, so the points that share an index are
  not readings of the same energy. Where the two do not overlap the sum is
  blank and the line simply stops. A scan ticked in "Invert" is inverted
  first, so inverting one and summing is the same curve as differencing. A
  difference is far smaller than either scan it came from -- 212 - 213 in
  LockAC spans -25852..5567 where both parents sit at 2.6e6, a flat line on
  their axis -- so it too takes the right-hand scale while that scale is
  free, and "Right scale" lists it either way to move it back,
- a "Baseline…" button opening a separate window with the
  ``entry/instrument/bluesky/streams/baseline`` stream for one selected file at
  a time, as a before/after table (the baseline stream records each field once
  before the scan and once after),
- a "Load Image" button that opens the detector frames belonging to the
  selected file(s) in a separate 2D image window (see ``image_viewer.py``).
  Selecting several files puts their images side by side in one window,
  sharing ROIs, so scans can be compared. The frames are not in the master
  file at all -- depending on which writer took them they live either in a
  sibling ``pva*.h5`` or in a per-detector subdirectory such as
  ``eiger/scan_000310.h5``, and the master records which,
- a "Load Spectrum" button that opens the multi-element fluorescence spectra
  belonging to the selected file(s) in a separate window (see
  ``spectrum_viewer.py``). Those live in a ``vortex`` subdirectory beside the
  master, as a ``(frames, elements, bins)`` MCA cube, and the ``Spc`` column
  of the listing says which scans have one. Selecting several overlays them on
  one pair of axes, which is the comparison worth making for a curve where
  side by side is the one worth making for an image,
- an "Export…" button writing out whatever is plotted, which kind of file
  decided by the name given. A ``.csv``, ``.dat`` or ``.txt`` name writes the
  numbers: an X and a Y column per curve, curves drawn against one sampling
  sharing its X column, a ``dy[..]`` column beside each Y while the
  derivative is shown, and a commented header naming the scans, the stream,
  the ticks and -- reducing -- the channels, the tolerances and how well each
  scan reduced, so the file can be read without this window. Columns of
  different lengths are squared off with blanks; a reading that came back
  NaN is written ``nan``, which is not the same thing. Any name matplotlib
  can render (``.png``, ``.pdf``, ``.svg`` and more) saves the picture
  instead, as it stands on screen, zoom and all. A colour map's numbers are
  left out -- it saves only as a picture. Files are offered in an ``export``
  folder under the directory being browsed, named ``data_000171_000173`` or
  ``figure_000171_000173`` after the scans in the plot. The convention
  followed, down to the header, is ``xmucd_live.py``'s own export.

Each ``primary``/``baseline`` field is a subgroup holding ``EPOCH``, ``time``
and (for a scalar field) ``value`` datasets -- confirmed by inspecting a real
master file -- rather than a plain top-level dataset the way
``entry/data`` stores its (smaller) set of fields. ``lambda250k`` itself has
no ``value`` (its image lives elsewhere, linked separately) and is skipped.

``image_viewer.py`` and ``spectrum_viewer.py`` must sit beside this file -- the
image and spectrum windows are imported from them.

Run with the ``6idb-bits`` environment's Python, with a display available::

    DISPLAY=localhost:10.0 \
        /home/beams18/USER6IDB/.conda/envs/6idb-bits/bin/python scan_viewer.py \
        [initial_directory]
"""

from __future__ import annotations

import csv
import re
import sys
import warnings
from datetime import datetime
from functools import partial
from pathlib import Path

import h5py
import numpy as np

# Before matplotlib, which picks its Qt binding from whatever is already
# imported -- see the longer note in image_viewer.py.
from PyQt6 import QtCore
from PyQt6 import QtWidgets

from matplotlib import rcParams
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter
from matplotlib.transforms import Bbox
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from .image_viewer import ImageViewerWindow
from .image_viewer import resolve_image_file

# The softglue spectrum viewer, which is the only one packaged here. Upstream
# imported the plain spectrum_viewer and scan_viewer.py rebound both of these
# names to the softglue fork before building a window; pointing the import
# straight at it leaves that rebinding a no-op and drops 2185 lines of module
# that nothing could reach.
from .spectrum_viewer import SpectrumViewerWindow
from .spectrum_viewer import resolve_spectrum_file

# numpy 2 renamed trapz to trapezoid and removed the old spelling. Nothing
# here pins numpy, so both are still in the field: 1.26 has only trapz and
# 2.x only trapezoid. Bound once, by name, rather than at each call site.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

STREAMS_PATH = "entry/instrument/bluesky/streams"
PRIMARY_STREAM = "primary"
PRIMARY_PATH = f"{STREAMS_PATH}/{PRIMARY_STREAM}"
METADATA_PATH = "entry/instrument/bluesky/metadata"
# Offered alongside the real fields: each primary field carries its own EPOCH
# (seconds past the UNIX epoch, one per scan point), so absolute time can be an
# axis like anything else.
EPOCH_FIELD = "EPOCH"
BASELINE_PATH = f"{STREAMS_PATH}/baseline"
# Two points, before and after, rather than a curve -- it has its own window.
NON_PLOT_STREAMS = {"baseline"}

# The reduced stream a dichro scan writes beside primary: the XAS and XMCD
# worked out from the polarization states measured at each energy, so it holds
# a quarter of primary's points and none of its raw channels. Its field names
# are fixed by the writer, and the scan's own ``hints`` describe primary rather
# than this stream, so what to plot of it is spelled out here instead.
DICHRO_STREAM = "dichro_monitor"
STREAM_DEFAULTS = {
    DICHRO_STREAM: {
        "x": ("dichro_positioner1",),
        "y": ("dichro_xas", "dichro_xmcd"),
        # XAS and XMCD are worked out from monitor-normalized counts already.
        "norm": (),
        # XMCD is a difference of two absorption measurements, often a percent
        # or less of the XAS it rides on; sharing one scale would flatten it
        # into the axis, so it is drawn against its own on the right.
        "right": ("dichro_xmcd",),
        # The same field again, for the same reason seen from the colour map's
        # side: XMCD is signed, so zero is a value the sample has rather than
        # the bottom of the range, and which side of it a cell falls on is the
        # result. Its map is drawn in a diverging scale centred there. XAS is
        # not on this list -- it is an absorption that merely happens to come
        # out negative, with no meaning at zero, and a diverging scale would
        # put white at an arbitrary one of its values.
        #
        # Named rather than detected. Of the 22 dichro scans in the two test
        # directories only 7 have XMCD straddling zero: 8 are wholly positive
        # and 2 wholly negative. All 22 want the same treatment, so a "does it
        # cross zero?" test would get it wrong more often than right.
        "diverging": ("dichro_xmcd",),
        # A magneto-chiral scan (``macros_other_scans.py``) steps a second
        # positioner -- the magnet current -- inside the scanned one, and this
        # stream holds one XAS/XMCD pair per step. Drawn as a single curve it
        # runs back and forth across the plot; split on this field it is one
        # curve per step, which is what the scan measured.
        "split": "dichro_positioner2",
        # On a mesh scan the same two positioner fields are the raster's own
        # axes, outermost first -- whatever the motors behind them are called.
        "grid_axes": ("dichro_positioner1", "dichro_positioner2"),
    }
}
# How near two steps' readbacks must sit to count as the same setting when
# "Average repeats" gathers them. Kepco readbacks wander by a few mA, so the
# steps sitting at one field are gathered by value rather than by equality.
# The same 0.03 A as process_data.py and xmucd_reduce.py use on the same
# readback. A scan stepping two positioners must agree on both of them.
SPLIT_ATOL = 0.03
# The four-contrast reduction, transcribed from xmucd_reduce.py rather than
# imported: that module lives in a per-cycle data directory and pulls in
# polartools and pyplot behind it, so a hard-wired path would stop working at
# the next cycle. The arithmetic is checked against it -- see reduce_dichro.
#
# Binning tolerances for the three grid axes. The magnet keeps SPLIT_ATOL,
# which is already the 0.03 A those references use on this same readback.
ENERGY_ATOL = 0.0001
PR_ATOL = 0.1
# The panels, in the order they are drawn. "XMuCD" is the magneto-chiral term;
# "XNCD" is field-even and so should sit at the noise floor -- where it does
# not, the polarization average is not cancelling XMCD and XMuCD is
# contaminated, which is what makes it worth a panel of its own.
CONTRASTS = ("xas", "xmcd", "xmucd", "xncd")
CONTRAST_LABELS = {
    "xas": "XAS",
    "xmcd": "XMCD",
    "xmucd": "XMuCD",
    "xncd": "XNCD",
}
CONTRAST_NOTES = {
    "xas": "absorption",
    "xmcd": "odd in helicity and field",
    "xmucd": "magneto-chiral, field-odd",
    "xncd": "field-even -- expect noise",
}
# The two ways the same four contrasts are measured here. A dichro scan steps
# the retarder against the magnet and gives four groups per energy; a lock-in
# scan has the retarder oscillating, so the helicity difference arrives
# ready-made on the AC channel and only the field is stepped -- two groups.
# "lockin_field", the third mode xmucd_reduce.py knows about, is deliberately
# absent: no scan here is one, and the two lock-in modes differ only in which
# of XMuCD and XNCD comes off the DC pair, with nothing in the file to say
# which. Guessing it would swap two panels in silence.
MODE_DICHRO = "dichro"
MODE_LOCKIN_POL = "lockin_pol"
MODE_LABELS = {
    MODE_DICHRO: "dichro, 4-quadrant",
    MODE_LOCKIN_POL: "lock-in, polarization modulated",
}
# Which raw channel each panel came off, per mode: "DC" for a panel in
# absorption units, "AC" for one in the lock-in's own. The AC difference is
# used exactly as the scaler wrote it -- neither normalized by the monitor nor
# put through the transmission log -- so the two are not comparable in
# magnitude, which is why the panels say which they are.
CONTRAST_SOURCES = {
    MODE_DICHRO: {key: "DC" for key in CONTRASTS},
    MODE_LOCKIN_POL: {
        "xas": "DC",
        "xmucd": "DC",
        "xncd": "AC",
        "xmcd": "AC",
    },
}
# The hints value that says a scan locked in -- qxscan_current_pr(lockin=True)
# writes it. Read only in the positive direction; see reduce_channels.
LOCKIN_HINT = "lockin"
# The three columns the reduction groups by. Named rather than discovered:
# they are the beamline's own devices, and the same names process_data.py
# works from.
ENERGY_FIELD = "energy"
PR_FIELD = "pr2_pzt_localdc"
MAGNET_FIELD = "bmag_kepco_current"
REDUCE_FIELDS = (ENERGY_FIELD, PR_FIELD, MAGNET_FIELD)
# The lock-in scaler's own channels, named by the scaler and so without a
# device prefix. The difference the lock-in measured is LockAC - LockACoff;
# LockDC is the ordinary absorption alongside it.
DC_FIELD = "LockDC"
AC_FIELD = "LockAC"
ACOFF_FIELD = "LockACoff"
LOCKIN_FIELDS = (DC_FIELD, AC_FIELD, ACOFF_FIELD)
TEMPERATURE_FIELD = "temp_336_4idg_loop1"
# How far the cryostat may drift with the scans still one measurement. The
# allowance cannot go below the drift it has to survive, which is 0.91 K
# between 235 and 236 -- a single measurement written to two files.
TEMPERATURE_ATOL = 2.0
# What a lock-in scan is grouped by when its field was stepped between files
# rather than inside one -- see :func:`lockin_groups`, which explains the
# tolerances. Deliberately not the field itself.
GROUP_FIELDS = {
    "bmag_sx": 0.01,
    "bmag_sy": 0.01,
    "bmag_srot": 0.1,
    TEMPERATURE_FIELD: TEMPERATURE_ATOL,
}
# Where the sample sat, for telling repeats of one measurement from a series
# that stepped something between them -- see :func:`scan_groups`. One set per
# endstation, the two in this data sharing no sample stage, and a scan has to
# carry a whole set before it can be grouped at all: a stage nobody listed
# here would otherwise leave every scan standing on it looking alike. The
# tolerances are the measured gaps. A series steps bmag_sy by 0.2 mm (174,
# 175, 176) or huber_euler_chi by 90 deg (741 to 783); repeats hold every one
# of these to four decimals and move only the temperature, by up to 0.9 K.
SAMPLE_STAGES = (
    {"bmag_sx": 0.01, "bmag_sy": 0.01, "bmag_srot": 0.1},
    {
        "huber_euler_x": 0.01,
        "huber_euler_y": 0.01,
        "huber_euler_z": 0.01,
        "huber_euler_chi": 0.05,
        "huber_euler_phi": 0.05,
    },
)
# How far two scans' X may drift apart and still be the same points measured
# twice, as a fraction of one step. Past half a step, point i of one scan is
# nearer point i + 1 of the other and averaging them row by row would smear
# the curve rather than clean it. Every repeat set in this data sits between
# 0.03 and 0.07 of a step; scans 4 and 6 sit at exactly 1.0, sharing a plan,
# a motor, a range and 50 points but stepping them in a different order.
SAME_X_STEPS = 0.5
# What the Detector and Monitor boxes fall back on, best first, when the scan's
# own hints name nothing that works. They are needed more often than not: a
# hint can be wrong in two ways that look identical on an empty plot, and both
# happen here. Scan 96 hints the monitor "Time", which no stream records --
# Time1 and Time2 are what it writes -- and the scans before 33 hint bI0, which
# read a flat zero for that whole run. Seeding from the hints alone leaves 13
# of the 14 reducible scans in the test directory opening on empty panels.
#
# The 4id<x>I0/I1 and eiger/vortex entries are the beamline's modern naming
# (2026 cycles onwards): the ion chambers at the b/g/h stations and the two
# counting detectors the mesh scans use for fluorescence and scattering. Scans
# 753 and 754 of NdFeB_26-1 record no detector hint at all, so the lists
# below are their only guide; without the modern names reduce_defaults falls
# through to channels[0] and picks 4idbI0 twice -- a dead channel over itself.
REDUCE_DETECTORS = (
    "eiger_stats1_total",
    "vortex_total_roi1",
    "vortex_total_roi2",
    "4idbI1",
    "4idgI1",
    "4idhI1",
    "bI1",
    "bI0",
    "gI1",
    "APD",
)
REDUCE_MONITORS = (
    "4idgI0",
    "4idbI0",
    "4idhI0",
    "bI0",
    "gI0",
    "hI0",
    "Time1",
)
# Where the "Clip map scale" box starts: on a 1500-cell raster, the eight most
# extreme cells either way. Small enough that a map with nothing wrong with it
# barely changes, which is what a starting point should do; one or two cells
# are what usually wreck a scale. Type over it when more needs to go -- see
# clipped_limits for why the amount is asked for rather than worked out.
COLOR_CLIP_PERCENT = 0.5
# Half the data from each end is the whole of it, and the two percentiles meet
# at the median. Stop just short, where there is still a scale left to draw.
COLOR_CLIP_MAX = 49.0
# Below this share of its own median the monitor is not dim, it is off: the
# beam was gone for that point. Unlike COLOR_CLIP_PERCENT above this is not a
# starting point to type over, because there is nothing to tune. Over the 988
# scans in both data directories that carry a monitor, the worst live point
# sits at 0.8359 of its scan's median and the best dead one at 0.0208 -- a gap
# of forty, with nothing anywhere inside it. Any value from 0.05 to 0.5 draws
# the same line; this one sits at the dead end of the gap so that a real but
# weakly illuminated point is never taken for a dropout.
BEAM_GONE_FRACTION = 0.05
# The colour maps the box offers, sequential first and then diverging. Short
# on purpose: matplotlib ships upwards of a hundred and a list that long is
# harder to pick from than it is useful. These cover the two things a map can
# be asked to say.
COLOR_MAPS = (
    "viridis",
    "plasma",
    "inferno",
    "magma",
    "cividis",
    "turbo",
    "gray",
    "jet",
    "seismic",
    "coolwarm",
    "bwr",
    "RdBu_r",
)
# Those of COLOR_MAPS built around a middle colour rather than running dark to
# bright. The distinction is the whole of what "Center on 0" acts on: a
# sequential map says "more", and has nothing to centre; a diverging one says
# "which side", which is only true where the middle of the scale is zero.
# RdBu rather than RdBu_r would run blue high, which reads backwards.
DIVERGING_COLOR_MAPS = frozenset({"seismic", "coolwarm", "bwr", "RdBu_r"})
# Field names whose signed meaning the stream config cannot convey -- a
# per-cell XMCD worked out from a dichro mesh scan lives in the primary
# stream, not the reduced one, and its field name is the only hint left that
# it is signed. XAS is deliberately absent: see the diverging-list note at
# :data:`STREAM_DEFAULTS`.
DIVERGING_FIELDS = frozenset({"xmcd", "xmucd", "xncd"})
# What the box starts on, and what it means: let each map choose, by the field
# it is drawing. Naming a map in the box instead applies it to every panel.
AUTO_COLOR_MAP = "Automatic"
DEFAULT_COLOR_MAP = "viridis"
DICHRO_COLOR_MAP = "seismic"
# How many curves the peak-statistics table shows before it starts to scroll.
# Deliberately few: the plot above it is what the window is for. Counted in
# curves rather than rows because with dY/dX shown each curve brings a second
# row for its derivative, and the pair should be as visible as the single row
# was -- see HdfViewer._size_stats_table.
STATS_TABLE_CURVES = 3
# What marks a "Multi" component's row as belonging to the fitted row above
# it. On screen it is just an indent; in the export it is also the one thing
# that says whose peak a line describes, since "peak 1 of 2" on its own
# names no curve once several are plotted. No curve label can begin with it
# -- they are built from scan numbers and field names -- so it discriminates
# as well as it indents.
STATS_PEAK_INDENT = "    "
STATS_COLUMNS = [
    "Curve",
    "Max Y",
    "@X (max)",
    "Min Y",
    "@X (min)",
    "Center",
    "FWHM",
]
# Appended to those, and only while "Integrate" is ticked. Kept off by
# default because the table shares its height with the plot, and two more
# columns is a cost every other scan would pay for a number it never asked
# for -- see HdfViewer._set_stats_columns.
STATS_AREA_COLUMNS = ["Area", "Net area"]
# The same two quantities as columns of the fit export, kept here beside
# the table's names for them so the two cannot drift apart. Named away
# from the "area" the fitted quantities already spell, which is the
# model's height x fwhm x shape integral and a different number -- see
# fit_export_columns and HdfViewer._export_fits.
CURVE_AREA_EXPORT = [("curve_area", "area"), ("curve_net_area", "net_area")]
# Appended after those, and only while "Fit" is ticked, for the same reason.
# "Center"/"FWHM" above are half-max crossings of the data, and "Area"/"Net
# area" are trapezoids over the points drawn; these five come from the
# fitted model and will not agree with either -- see peak_fit.
STATS_FIT_COLUMNS = [
    "Fit center",
    "Fit FWHM",
    "Fit height",
    "Fit area",
    "R2",
]
# Slipped in among those, and only while the shape being fitted has a
# mixing to report -- the pseudo-Voigt and nothing else, see PEAK_MIX.
# Conditional for the reason the area pair is: a column costs table
# height on every scan, and for three shapes of the four this one would
# be blank for ever. It sits beside "Fit FWHM" because the two are the
# halves of one description -- the width says how broad the peak is and
# the mixing says what shape it is -- and because R2 stays last, which
# is where a reader goes looking for it.
STATS_MIX_COLUMN = "Fit eta"
# What "Plot" can trend across a selection, mapped to the key the fit
# returns it under. The display names are the table headings above,
# letter for letter, so the combo and the column it came from say the
# same word. "Fit height" is first because it is the default.
#
# The first three are the intensity, the width and the integrated
# intensity originally asked for. The centre was added afterwards, on
# request: against scan number it is a drift plot, but against one of
# the baseline fields the axis combo offers it is the useful one -- the
# peak position as a function of whatever was stepped, which is what
# scans 32 to 39 are stepping huber_euler_ana_eta to measure. The
# mixing is last and is blank for every shape that has none.
TREND_QUANTITIES = {
    "Fit height": "fit_height",
    "Fit FWHM": "fit_fwhm",
    "Fit area": "fit_area",
    "Fit center": "fit_center",
    "Fit eta": "fit_eta",
}
# The trend's X axis when it is not a baseline field, and the default.
TREND_X_SCAN = "Scan number"
# What counts as "deliberately moved" when deciding which baseline fields are
# worth offering as that axis. A field's noise is taken as how far it drifted
# within a scan -- it is read before and after, so every scan carries its own
# estimate -- plus a floor of BASELINE_TREND_NOISE of its own magnitude, for
# the setpoints that never drift at all because they are digital. A field is
# offered when its motion across the selection beats that by
# BASELINE_TREND_SNR.
#
# Both halves earn their keep, measured on scans 32 to 39, where
# huber_euler_ana_eta is stepped 0 to 90 degrees and 136 of the 519 numeric
# fields move at all. Without the floor the list is headed by energy and
# mono_th wobbling in their fifth decimal -- clean, since a setpoint does not
# drift, but four parts in a billion of their own value. Without the ratio it
# is headed by the gsydor beam diagnostics, which swing further than anything
# else here and are noise; they are why a plain "did it change" filter cannot
# do this. Together they leave 19 fields with eta third, and dropped every
# gsydor field on all six selections tried.
BASELINE_TREND_NOISE = 1e-4
BASELINE_TREND_SNR = 10.0
HDF_SUFFIXES = (".hdf", ".h5", ".hdf5")
FILE_COLUMNS = [
    "Scan",
    "Img",
    "Spc",
    "Plan",
    "Motor",
    "Range",
    "Pts",
    "Count [s]",
]
# Starting widths. Every column is user-resizable from here; the motor name is
# the one that can run long, so it starts wide and is elided rather than
# allowed to shove the other columns aside.
FILE_COLUMN_WIDTHS = [58, 38, 38, 72, 130, 118, 40, 66]
# Asked of the list rather than written down, so that inserting a column moves
# them instead of leaving one quietly pointing at its neighbour.
IMAGE_COLUMN = FILE_COLUMNS.index("Img")
SPECTRUM_COLUMN = FILE_COLUMNS.index("Spc")
# The glyph marks a file of some kind, not an image: Spc uses it too.
HAS_FILE, NO_FILE = "✓", "–"
FIELD_COLUMNS = ["Field", "X", "Y", "Norm"]
FIELD_X, FIELD_Y, FIELD_NORM = 1, 2, 3


def list_hdf_files(directory: Path) -> list[Path]:
    """Return the HDF5 files directly inside *directory*, sorted by name."""
    return sorted(
        p
        for p in directory.iterdir()
        if p.is_file() and p.suffix.lower() in HDF_SUFFIXES
    )


def _is_time_like(name: str) -> bool:
    return "time" in name.lower() or name == EPOCH_FIELD


def scan_streams(path: Path) -> list[str]:
    """The streams of *path* that can be plotted as curves, ``primary`` first.

    Everything the writer put under ``streams`` that holds at least one
    value-bearing field, less the ones that are not curves at all
    (:data:`NON_PLOT_STREAMS`). A dichro scan adds :data:`DICHRO_STREAM`.
    """
    try:
        with h5py.File(path, "r") as f:
            group = f.get(STREAMS_PATH)
            if group is None:
                return []
            names = [
                name
                for name in group.keys()
                if name not in NON_PLOT_STREAMS
                and isinstance(group[name], h5py.Group)
                and any(
                    isinstance(field, h5py.Group) and "value" in field
                    for field in group[name].values()
                )
            ]
    except OSError:
        return []
    # primary is the raw scan and the one every file has, so it leads.
    return sorted(names, key=lambda name: (name != PRIMARY_STREAM, name))


def stream_fields(path: Path, stream: str) -> list[str]:
    """Field names under *stream* that carry a ``value`` dataset.

    h5py iterates group keys in name order (no ``track_order`` here), so this
    is also the order ``default_x_field`` treats as "last".
    """
    try:
        with h5py.File(path, "r") as f:
            group = f.get(f"{STREAMS_PATH}/{stream}")
            if group is None:
                return []
            names = [
                name
                for name in group.keys()
                if isinstance(group[name], h5py.Group)
                and "value" in group[name]
            ]
            if _epoch_source(group) is not None:
                names.append(EPOCH_FIELD)
            return names
    except OSError:
        return []


def _epoch_source(group) -> str | None:
    """The field whose EPOCH stands for the scan's clock.

    Every field records its own EPOCH and they disagree by a few seconds --
    each device is stamped when it was read -- so one has to be picked. The
    first in the group's (alphabetical) order is used, consistently, rather
    than silently mixing time bases between axes.
    """
    for name in group.keys():
        item = group[name]
        if isinstance(item, h5py.Group) and "EPOCH" in item:
            return name
    return None


def default_x_field(fields: list[str]) -> str | None:
    """The last field, skipping trailing time-like ones (e.g. ``scaler_time``).

    Only a fallback: :func:`scan_x_field` asks the scan what it moved.
    """
    for name in reversed(fields):
        if not _is_time_like(name):
            return name
    return fields[-1] if fields else None


def _resolve_hint_field(name: str, fields: list[str]) -> str | None:
    """The recorded field a hinted name stands for, or ``None`` for none.

    A motor is hinted under its own name but may only have been written as a
    setpoint, so the readback is preferred and the setpoint accepted after it.
    ``time`` is Bluesky's pseudo-axis for a plan that moved nothing; here that
    is the synthesized :data:`EPOCH_FIELD`.
    """
    if not name:
        return None
    if name == "time":
        return EPOCH_FIELD if EPOCH_FIELD in fields else None
    for candidate in (name, f"{name}_user_setpoint", f"{name}_setpoint"):
        if candidate in fields:
            return candidate
    return None


def stream_default_fields(
    stream: str, key: str, fields: list[str]
) -> list[str]:
    """A stream's own default ``x``/``y``/``norm``/``right`` fields that *fields* holds.

    Empty for a stream with no entry in :data:`STREAM_DEFAULTS`, which is every
    stream whose contents the scan's ``hints`` already describe.
    """
    wanted = STREAM_DEFAULTS.get(stream, {}).get(key, ())
    return [name for name in wanted if name in fields]


def scan_x_field(
    path: Path, fields: list[str], stream: str = PRIMARY_STREAM
) -> str | None:
    """The field to plot against: the axis this scan actually scanned.

    Taken from the scan's own ``hints`` -- ``dimensions`` names the axis of
    each scanned dimension -- and failing that from the ``motors`` metadata,
    rather than guessed from the field ordering. Returns ``None`` when the scan
    names no axis that was recorded. A stream the hints do not describe uses
    its own default instead (see :data:`STREAM_DEFAULTS`).
    """
    if stream in STREAM_DEFAULTS:
        default = stream_default_fields(stream, "x", fields)
        return default[0] if default else None
    info = read_scan_info(path)
    for name in (*info["hints"]["dimensions"], *info["motors"]):
        resolved = _resolve_hint_field(name, fields)
        if resolved is not None:
            return resolved
    return None


def scan_y_fields(
    path: Path, fields: list[str], stream: str = PRIMARY_STREAM
) -> list[str]:
    """The detectors the scan itself nominates to plot, from its ``hints``.

    The hints name the one signal per detector worth looking at --
    ``eiger_stats1_total`` out of the eiger's three dozen columns -- so they,
    rather than nothing, are what the Y ticks start on. Returns an empty list
    when the scan hints at nothing, or at fields it did not record. A stream
    the hints do not describe uses its own defaults -- for a dichro scan, the
    XAS and XMCD it was run to measure.
    """
    if stream in STREAM_DEFAULTS:
        return stream_default_fields(stream, "y", fields)
    chosen: list[str] = []
    for name in read_scan_info(path)["hints"]["detectors"]:
        resolved = _resolve_hint_field(name, fields)
        if resolved is not None and resolved not in chosen:
            chosen.append(resolved)
    return chosen


def default_monitor_field(
    fields: list[str], stream: str = PRIMARY_STREAM
) -> str:
    """What the Norm tick starts on -- nothing, unless the stream asks.

    A scan's ``hints`` nominate a monitor and this deliberately passes it over.
    Normalizing rescales every number on the plot and in the statistics table,
    so a scan opens showing what it recorded, with Norm left a tick away. A
    stream spec may still name one, that being a fact about how the stream was
    computed rather than a guess at intent -- though ``dichro_monitor``, whose
    values are normalized already, names none either, so in practice the
    column now always starts clear.
    """
    if stream in STREAM_DEFAULTS:
        default = stream_default_fields(stream, "norm", fields)
        return default[0] if default else ""
    return ""


def scan_number(path: Path) -> str:
    """The scan number in a file name: ``00008`` from ``tise2_2_00008_master.hdf``.

    The last run of digits wins, so the prefix's own digits (``tise2_2``) are
    ignored and an image file reduces to the same number as its master, whether
    it is named ``pva_tise2_2_00008.h5`` or ``eiger/scan_000310.h5``.
    """
    stem = path.stem
    if stem.endswith("_master"):
        stem = stem[: -len("_master")]
    digits = re.findall(r"\d+", stem)
    return digits[-1] if digits else stem


def file_prefix(path: Path) -> str:
    """The part of a file name in front of the scan number: ``tise2_2``.

    Everything from the scan number onwards is dropped, so the files of one
    series all reduce to the same string.
    """
    stem = path.stem
    if stem.endswith("_master"):
        stem = stem[: -len("_master")]
    digits = list(re.finditer(r"\d+", stem))
    if not digits:
        return stem
    return stem[: digits[-1].start()].rstrip("_-. ") or stem


def scan_tags(paths: list[Path]) -> dict[str, str]:
    """A short name per path, keyed by ``str(path)``: enough to tell them apart.

    The scan number alone normally does it. Where it does not -- two series in
    one directory whose numbers collide, :func:`scan_number` taking the last
    run of digits either way -- every name gains its series. The file listing
    gets away with the bare number because its other six columns disambiguate;
    a bare menu of numbers would not.
    """
    numbers = [scan_number(path) for path in paths]
    if len(set(numbers)) != len(numbers):
        return {
            str(path): f"{file_prefix(path)} {number}"
            for path, number in zip(paths, numbers)
        }
    return {str(path): number for path, number in zip(paths, numbers)}


def _as_text(value) -> str:
    """An h5py scalar -- bytes, numpy scalar or 0-d array -- as plain text."""
    value = np.asarray(value).reshape(-1)
    if value.size == 0:
        return ""
    item = value[0]
    return (
        item.decode(errors="replace") if isinstance(item, bytes) else str(item)
    )


def _meta_text(group, key: str) -> str:
    item = group.get(key)
    return "" if item is None else _as_text(item[()])


def _yaml_items(text: str) -> list[str]:
    """The entries of a simple YAML block sequence (``- name`` per line).

    ``motors`` and ``detectors`` are stored as YAML text (``motors`` even keeps a
    ``!!python/tuple`` tag), so pull the items out directly rather than take a
    PyYAML dependency for two one-line lists.
    """
    return [
        line[2:].strip().strip("'\"")
        for line in text.splitlines()
        if line.startswith("- ")
    ]


def _yaml_block_items(text: str, key: str) -> list[str]:
    """The entries of the block sequence under the top-level ``key:`` of *text*.

    ``hints`` is a small YAML mapping rather than a bare list, so
    :func:`_yaml_items` -- which takes every ``- `` line in the document --
    would mix its sections together. Collect only the lines belonging to one
    key, stopping at the next top-level one.
    """
    items: list[str] = []
    collecting = False
    for line in text.splitlines():
        if not line.strip():
            continue
        if not collecting:
            collecting = line.rstrip() == f"{key}:"
            continue
        if line.startswith("- "):
            items.append(line[2:].strip().strip("'\""))
        elif not line[:1].isspace():
            break  # a new top-level key ends the block
    return items


def _yaml_scalar(text: str, key: str) -> str:
    """The value of a top-level ``key: value`` line in *text*."""
    prefix = f"{key}:"
    for line in text.splitlines():
        if line.startswith(prefix):
            return line[len(prefix) :].strip().strip("'\"")
    return ""


def _hint_dimensions(text: str) -> list[str]:
    """The axes named by the ``dimensions`` hint, outermost first.

    Bluesky hints one ``(fields, stream)`` tuple per scanned dimension, which
    PyYAML writes two or three levels deep with ``!!python/tuple`` tags::

        dimensions:
        - !!python/tuple
          - - bmag_sx
            - bmag_sy
          - primary

    Only the field names are wanted, so the trailing stream name -- the last
    leaf of each tuple -- is dropped. A plan that moved nothing hints the
    pseudo-axis ``time`` (see :func:`_resolve_hint_field`).
    """
    groups: list[list[str]] = []
    collecting = False
    for line in text.splitlines():
        if not line.strip():
            continue
        if not collecting:
            collecting = line.rstrip() == "dimensions:"
            continue
        if line.startswith("- "):
            groups.append([])  # one more (fields, stream) tuple
        elif not line[:1].isspace():
            break  # a new top-level key ends the block
        item = line.strip()
        while item.startswith("- "):
            item = item[2:].strip()
        if item and not item.startswith("!!") and groups:
            groups[-1].append(item.strip("'\""))
    return [name for group in groups for name in group[:-1]]


# Reading this opens the file, and the browser re-reads the whole selection on
# every selection change, so remember it per (path, mtime, size).
# The diffractometer's own reading of where it was pointing, per scan point.
# Older scans do not record it.
HUBER_HKL_FIELDS = ("huber_euler_h", "huber_euler_k", "huber_euler_l")

_SCAN_HKL_CACHE: dict[tuple, tuple | None] = {}
_SCAN_INFO_CACHE: dict[tuple, dict] = {}
_REDUCE_CACHE: dict[tuple, list[str]] = {}
_BASELINE_CACHE: dict[tuple, dict[str, tuple[float, float]]] = {}
# Keyed on the selection, since the ranking in baseline_trend_fields is of
# one selection against itself and every redraw asks for it again. That
# function says what the key is and why.
_TREND_FIELD_CACHE: dict[tuple, list[str]] = {}


def scan_hkl(path: Path):
    """Mean and spread of huber_euler_h/k/l over the scan, or None if not recorded.

    The spread matters as much as the mean: on a rocking curve one of the three
    sweeps across the peak while the others hold still, and a bare average
    would read like a fixed position.
    """
    path = Path(path)
    try:
        stat = path.stat()
        key = (str(path), stat.st_mtime, stat.st_size)
    except OSError:
        key = None
    if key is not None and key in _SCAN_HKL_CACHE:
        return _SCAN_HKL_CACHE[key]

    result = None
    if set(HUBER_HKL_FIELDS) <= set(stream_fields(path, "primary")):
        try:
            data = read_stream(path, list(HUBER_HKL_FIELDS))
        except OSError:
            data = {}
        gathered = []
        for name in HUBER_HKL_FIELDS:
            values = np.asarray(data.get(name, ()), dtype=float).ravel()
            values = values[np.isfinite(values)]
            if values.size == 0:
                gathered = []
                break
            gathered.append((float(values.mean()), float(values.std())))
        result = tuple(gathered) if gathered else None
    if key is not None:
        _SCAN_HKL_CACHE[key] = result
    return result


def format_scan_hkl(values, decimals: int = 4) -> str:
    """``H 1.0135  K -0.4926 ± 0.0386  L 7.4840`` for one scan."""
    smallest = 0.5 * 10.0 ** (-decimals)
    parts = []
    for name, (mean, spread) in zip("HKL", values):
        text = f"{name} {mean:.{decimals}f}"
        if spread >= smallest:
            text += f" ± {spread:.{decimals}f}"
        parts.append(text)
    return "   ".join(parts)


def _dwell(seconds: np.ndarray) -> float | None:
    """The counting time per point in *seconds*, as the one number it mostly is.

    The commonest reading, not the first. A scan whose dwell is not the same
    at every point is rare -- 32 of the 1012 scans in the two test
    directories -- and in every one of those it is a stray point or two
    rather than a change of plan: scan 83 reads 0.2 s at 44 of its 45 points
    and 0.0 at the first, scan 172 reads 5.0 at 108 of 109. The first point
    is the one that strays, being the one the scaler was still settling into,
    so reading the dwell off the front reports 0.0 for a scan that counted
    for 0.2 s throughout -- which is what 18 of those 1012 would have done.

    No rounding before counting: the repeated readings are bit-identical in
    all 231 scans here that record one, so the only thing rounding could
    merge is strays, which are outvoted anyway.
    """
    seconds = seconds[np.isfinite(seconds)]
    if seconds.size == 0:
        return None
    values, counts = np.unique(seconds, return_counts=True)
    return float(values[np.argmax(counts)])


def read_scan_info(path: Path) -> dict:
    """What the scan was: plan, points, motors and the range they covered.

    The motor names come from the ``motors`` metadata, and each one's actual
    travel is read back from its own field in the primary stream -- the
    requested range in ``plan_args`` is not necessarily what the motor did.
    """
    path = Path(path)
    try:
        stat = path.stat()
        key = (str(path), stat.st_mtime, stat.st_size)
    except OSError:
        key = None
    if key is not None and key in _SCAN_INFO_CACHE:
        return _SCAN_INFO_CACHE[key]

    info = {
        "file": path.name,
        "is_scan": False,
        "plan": "",
        "points": "",
        "motors": [],
        "shape": [],
        "detectors": [],
        # What the scan says is worth looking at: the axes it scanned, the one
        # signal per detector to plot, and the monitor it counted against.
        # These are field names (``eiger_stats1_total``), unlike ``detectors``
        # above, which names the devices (``eiger``).
        "hints": {
            "dimensions": [],
            "detectors": [],
            "monitor": "",
            "gridding": "",
            "polarization": "",
        },
        "ranges": [],  # (motor, first, last, n)
        "count_time": None,
        "started": "",
        "plan_args": "",
        "uid": "",
    }
    try:
        with h5py.File(path, "r") as f:
            meta = f.get(METADATA_PATH)
            if meta is not None:
                info["is_scan"] = True
                info["plan"] = _meta_text(meta, "plan_name")
                info["points"] = _meta_text(meta, "num_points")
                info["motors"] = _yaml_items(_meta_text(meta, "motors"))
                # grid_scan records the raster size here; a 1-D plan has none.
                info["shape"] = [
                    int(n)
                    for n in _yaml_items(_meta_text(meta, "shape"))
                    if n.lstrip("-").isdigit()
                ]
                info["detectors"] = _yaml_items(_meta_text(meta, "detectors"))
                hints = _meta_text(meta, "hints")
                info["hints"] = {
                    "dimensions": _hint_dimensions(hints),
                    "detectors": _yaml_block_items(hints, "detectors"),
                    # Read because the file says it, but nothing ticks Norm
                    # from it -- see :func:`default_monitor_field`.
                    "monitor": _yaml_scalar(hints, "monitor"),
                    # ``grid_scan`` and its relatives declare their raster
                    # here; ``list_grid_scan``, which the magneto-chiral
                    # macros use, does not -- see :func:`is_mesh_scan`.
                    "gridding": _yaml_scalar(hints, "gridding"),
                    # ``lockin``, ``dichro`` or ``fixed``. Only the first is
                    # acted on -- see :func:`reduce_mode`.
                    "polarization": _yaml_scalar(hints, "polarization"),
                }
                info["plan_args"] = _meta_text(meta, "plan_args")
                info["uid"] = _meta_text(meta, "run_start_uid")

            primary = f.get(PRIMARY_PATH)
            if primary is not None:
                for motor in info["motors"]:
                    sub = primary.get(motor)
                    if sub is None or "value" not in sub:
                        continue
                    values = np.asarray(sub["value"][()], dtype=float).ravel()
                    if values.size:
                        info["ranges"].append(
                            (
                                motor,
                                float(values[0]),
                                float(values[-1]),
                                int(values.size),
                            )
                        )
                # The dwell per point, off the scaler's own clock channel.
                # That channel is named for the scaler it belongs to --
                # scaler1_time and scaler2_time here -- and never plainly
                # "scaler_time", which is what this used to look for and why
                # the column stood empty on all 232 scans in the test
                # directory and all 779 in lineup. Lowest-numbered wins: the
                # two carry the same dwell, differing only on the stray
                # points _dwell outvotes, in 12 of 231 scans.
                #
                # Not the "Time1" channel beside it, which is the same dwell
                # counted in ticks of the scaler's 10 MHz clock: 5e7 where
                # scaler1_time reads 5. Seconds are what the column asks for.
                clocks = sorted(
                    name
                    for name in primary
                    if name.startswith("scaler") and name.endswith("_time")
                )
                if clocks and "value" in primary[clocks[0]]:
                    info["count_time"] = _dwell(
                        np.asarray(
                            primary[clocks[0]]["value"][()], dtype=float
                        ).ravel()
                    )
                for name in primary:
                    field = primary[name]
                    stamp = (
                        field.get("time")
                        if isinstance(field, h5py.Group)
                        else None
                    )
                    if stamp is not None and "start_time_iso" in stamp.attrs:
                        info["started"] = _as_text(
                            stamp.attrs["start_time_iso"]
                        )[:19]
                        break
    except OSError:
        # A file still being written, or not HDF5 at all, contributes no
        # timestamp. The listing shows the row with blank cells rather than
        # dropping the scan, which is the documented behaviour.
        pass

    if key is not None:
        _SCAN_INFO_CACHE[key] = info
    return info


def _grid_axis_candidates(available: list[str], motor: str) -> list[str]:
    """Fields that could supply this motor's grid axis, best first.

    The readback comes first, because it is where the motor actually was. It is
    only usable when every row repeats its coordinates exactly, which a readback
    that wanders by a few counts will not; the ``*_user_setpoint`` companion
    holds the demanded value and always does. :func:`read_grid` tries them in
    this order and keeps the first that reconstructs a raster.
    """
    names = (motor, f"{motor}_user_setpoint", f"{motor}_setpoint")
    return [name for name in names if name in available]


def inner_axis_fields(path: Path, stream: str) -> list[tuple[str, str, str]]:
    """The scan's inner mesh axes, as ``(motor, readback, setpoint)`` each.

    The axes a scan steps inside its outermost one: for a magneto-chiral scan
    the magnet current, for a fixed-polarization one the polarization as well.
    They are what the stream has to be cut along to get a curve per setting
    rather than one running back and forth, and the scan names them itself.

    Both fields are wanted because they answer different questions. The
    readback says where the positioner actually sat, which is what a curve is
    labelled with; the setpoint says where it was *told* to sit, which is the
    only one of the two that can be trusted to repeat exactly -- see
    :func:`confirms_period`.

    Returns ``[]`` unless every axis yields both, which is what keeps a
    reduced stream on the split its own spec names: ``dichro_monitor`` calls
    these fields ``dichro_positioner1`` and ``2``, so no motor name reaches
    them.
    """
    info = read_scan_info(path)
    shape = info["shape"]
    axes = list(info["hints"]["dimensions"]) or list(info["motors"])
    # An ordinary scan lists its count-time signals among its motors, so the
    # names alone would suggest a raster where there is none. The declared
    # shape is what says how many axes were really stepped, and one that does
    # not account for every name means the names cannot be trusted either.
    if len(shape) < 2 or len(axes) != len(shape):
        return []
    available = stream_fields(path, stream)
    resolved = []
    for motor in axes[1:]:
        names = _grid_axis_candidates(available, motor)
        setpoints = [name for name in names if name != motor]
        if not names or not setpoints:
            return []
        resolved.append((motor, names[0], setpoints[0]))
    return resolved


def confirms_period(columns: list[np.ndarray], steps: int) -> bool:
    """Whether *columns* really do hold one value per step of the raster.

    A scan's declared shape says how many points make up one row of its mesh,
    but it does not always describe every repeat the stream holds: a dichro
    scan measures each mesh point at two polarizations, so its ``primary``
    carries twice the rows the shape accounts for and splitting on the
    declared count would put both polarizations in one curve. Asking the data
    whether the period holds is what tells the two apart, and it costs a scan
    that does step a plain raster nothing.

    Setpoints rather than readbacks: a readback is mid-ramp on the rows a
    supply spent reaching its setting, so it is unequal to itself within a
    step by amounts far past any sensible tolerance -- the magnet in scan 96
    reads 2 A either side of a step that is in fact one setting.
    """
    for column in columns:
        read_something = False
        for step in range(steps):
            sample = column[step::steps]
            sample = sample[np.isfinite(sample)]
            if not sample.size:
                continue
            read_something = True
            if len(np.unique(sample)) != 1:
                return False
        # A column that was never written holds nothing to confirm the period
        # against, which is not the same as confirming it.
        if not read_something:
            return False
    return True


def is_mesh_scan(path: Path) -> bool:
    """Whether the scan declares itself a rectilinear mesh.

    ``grid_scan`` and its relatives write ``gridding: rectilinear`` into their
    hints; ``list_grid_scan``, which the magneto-chiral macros use to cycle a
    magnet current inside the scanned axis, does not. Both plans leave two
    positioners in the reduced dichro stream, so this is what says whether its
    inner axis is a coordinate to raster over or a handful of settings to draw
    a curve each for.
    """
    return read_scan_info(path)["hints"].get("gridding") == "rectilinear"


def has_dichro_stream(path: Path) -> bool:
    """True when the scan carries a dichro_monitor stream with XAS and XMCD.

    Written by the beamline when the retarder was stepped: four points of
    primary (two helicities, two field directions) reduce to one XAS and one
    XMCD per energy, and the writer puts both into this stream. The 1D
    reducer here needs the retarder *and* the magnet column to run; where
    the magnet was held still the four-quadrant reduction cannot (one
    quadrant missing), yet the file already has the result the user came
    for, so this flag lets the viewer open on the panels anyway.
    """
    fields = stream_fields(path, DICHRO_STREAM)
    return "dichro_xas" in fields and "dichro_xmcd" in fields


def scan_kind(path: Path, fields: list[str], stream: str) -> tuple:
    """The stream, its columns, and what the scan's hints make of X, Y and Norm.

    This is what says whether two scans are the same kind of measurement, and
    so whether ticks made on one still describe the other -- see
    :meth:`HdfViewer._refresh_field_lists`. Two scans are the same kind when
    the same columns are on offer and the scans ask for the same three of
    them: the ticks stand until the scan on screen says something different
    from the one before it.

    Drawing the line here rather than at the plan keeps both halves honest.
    Scans 353 and 232 in the lineup directory are both ``rel_grid_scan`` over
    ``huber_euler_z`` recording the same columns, so a plan-and-motors rule
    would call them one kind -- but 232 asks for the Eiger's third ROI and
    353 for its second, and inheriting the third there would quietly draw the
    wrong map. Conversely a scan repeated over a narrower range asks for
    exactly what the one before it did, and re-reading the hints there would
    undo a correction just made by hand.

    The last four elements are the seed itself, so a caller that is about to
    seed from hints can unpack them rather than compute them twice.
    """
    return (
        stream,
        tuple(fields),
        scan_x_field(path, fields, stream) or default_x_field(fields) or "",
        tuple(scan_y_fields(path, fields, stream)),
        default_monitor_field(fields, stream),
    )


def reduce_mode(path: Path) -> str:
    """Which of the two reductions a scan asks for.

    :data:`MODE_LOCKIN_POL` where the scan's ``polarization`` hint reads
    ``lockin``, :data:`MODE_DICHRO` for everything else. The hint is read in
    that one direction only: ``lockin`` is written by the macro that turns the
    lock-in on and means it, while ``fixed`` -- which is what the great
    majority of scans say, every plain ``lup`` among them -- means nothing at
    all. So a positive reading switches modes and a negative one decides
    nothing, which is what :func:`reduce_channels` goes on to do with columns.
    """
    return (
        MODE_LOCKIN_POL
        if read_scan_info(path)["hints"]["polarization"] == LOCKIN_HINT
        else MODE_DICHRO
    )


def reduce_channels(path: Path) -> list[str]:
    """The channels a four-contrast reduction can run on, or ``[]`` for none.

    What is asked of ``primary`` depends on :func:`reduce_mode`. A dichro scan
    must record every one of :data:`REDUCE_FIELDS`, the columns deciding it
    rather than the hint -- see there for why the hint cannot refuse. A
    lock-in scan must record the energy and all of :data:`LOCKIN_FIELDS`, and
    neither of the other two: the retarder is oscillating rather than stepped,
    so it is never read per point, and the field is commonly stepped between
    two files rather than inside one -- scans 235 and 236 are one measurement
    at +30.02 A and -29.96 A. Where it is missing here it is taken from the
    baseline instead, so asking for it would turn away the usual case.

    Neither mode goes on the plan name, every reducible scan here being a
    ``list_grid_scan`` along with fifty that are not.

    Always ``primary``, whatever the Stream box shows. A reduced stream holds
    a quarter of the points with the quadrants already averaged away, so there
    is nothing left there to reduce -- yet scan 98 opens on ``dichro_monitor``
    and is one of the scans this was written for. Gating on the stream on
    screen would have hidden the button on exactly that scan.

    Everything but the synthesized ``EPOCH`` is offered, since which channel
    was the detector is the user's to say and the hints get it wrong often
    enough to matter -- see :func:`reduce_defaults`.

    :meth:`HdfViewer._update_plot` asks on every redraw and the answer opens
    the file, so it is remembered per (path, mtime, size) as
    :func:`read_scan_info` is.
    """
    path = Path(path)
    try:
        stat = path.stat()
        key = (str(path), stat.st_mtime, stat.st_size)
    except OSError:
        key = None
    if key is not None and key in _REDUCE_CACHE:
        return _REDUCE_CACHE[key]

    fields = stream_fields(path, PRIMARY_STREAM)
    present = set(fields)
    # A pol-stepped mesh (dichro map) does not record the energy or the
    # magnet because it steps neither -- the retarder is all that moves per
    # cell -- so the four-column requirement cannot stand. It is reducible
    # all the same, into one XAS and one XMCD map per selection; see
    # :func:`read_dichro_grid`.
    enough = (
        ENERGY_FIELD in present and set(LOCKIN_FIELDS) <= present
        if reduce_mode(path) == MODE_LOCKIN_POL
        else (
            set(REDUCE_FIELDS) <= present
            or (
                PR_FIELD in present
                and is_mesh_scan(path)
                and len(read_scan_info(path)["motors"]) == 2
            )
        )
    )
    channels = (
        [name for name in fields if name != EPOCH_FIELD] if enough else []
    )
    if key is not None:
        _REDUCE_CACHE[key] = channels
    return channels


def reduce_defaults(
    path: Path, channels: list[str], transmission: bool = True
) -> tuple[str, str]:
    """A detector and monitor to open the reduction on: the scan's, then a live pair.

    The hints are asked first -- they name the signal the scan was run to
    watch and the channel it counted against -- but only as a pair, since half
    a hint is not a reduction the scan described. Whatever comes out of that
    is then tried, and :data:`REDUCE_DETECTORS` and :data:`REDUCE_MONITORS`
    walked until one pair gives a finite absorption, the monitor varied first
    as the usual culprit.

    That last step is not belt and braces. A hint naming a channel that is
    recorded but dead opens the scan on four empty panels, which reads exactly
    like a scan with nothing in it, and on the test directory only 1 of the 14
    reducible scans survives its own hints. Probing first, all 14 open on the
    pair the offline reduction was in fact run with.

    A lock-in scan asks the hints for nothing but the monitor. Its detector is
    the lock-in's own :data:`DC_FIELD`, which is the channel the AC difference
    beside it was measured against; the scan still hints the ordinary
    detectors, and reducing on one of those would pair the modulated
    difference with a channel that never saw the modulation.
    """
    if not channels:
        return "", ""
    info = read_scan_info(path)
    hinted_detector = next(
        (name for name in info["hints"]["detectors"] if name in channels), ""
    )
    monitor_hint = info["hints"]["monitor"]
    hinted_monitor = monitor_hint if monitor_hint in channels else ""
    paired = bool(hinted_detector and hinted_monitor)
    lockin = reduce_mode(path) == MODE_LOCKIN_POL

    def order(first: str, fallbacks: tuple[str, ...]) -> list[str]:
        names: list[str] = []
        for name in ((first,) if first else ()) + fallbacks:
            if name in channels and name not in names:
                names.append(name)
        return names

    detectors = order(
        DC_FIELD if lockin else (hinted_detector if paired else ""),
        REDUCE_DETECTORS,
    )
    monitors = order(
        hinted_monitor if lockin or paired else "", REDUCE_MONITORS
    )
    if not detectors or not monitors:
        return channels[0], channels[0]
    try:
        data = read_stream(
            path, sorted(set(detectors + monitors)), PRIMARY_STREAM
        )
    except OSError:
        return detectors[0], monitors[0]
    for detector in detectors:
        for monitor in monitors:
            if (
                detector == monitor
                or detector not in data
                or monitor not in data
            ):
                continue
            if np.isfinite(
                absorption(data[detector], data[monitor], transmission)
            ).any():
                return detector, monitor
    # Nothing alive on any pair. Open on the best guess and let the status bar
    # say so: the boxes are there to be moved, and a dead I0 is a fact about
    # the scan rather than a reason to refuse it.
    return detectors[0], monitors[0]


def split_axes(
    path: Path, stream: str
) -> tuple[list[str], list[str], list[str]]:
    """What this stream splits into curves: axis names, values, setpoints.

    Two shapes of scan arrive here. A magneto-chiral scan is reduced by the
    writer into a stream that names its own split field, and the scan's own
    motor list supplies the name to label it with. A fixed-polarization scan
    has no reduced stream at all -- it steps its positioners in ``primary``
    itself -- so the axes come from the raster the scan declares, which is the
    only reason a setpoint column is wanted as well: see
    :func:`confirms_period`.

    Empty lists where there is nothing to split along, a mesh among them: its
    inner axis sweeps a coordinate, and cutting it into one curve per position
    would give fifty of them.
    """
    if is_mesh_scan(path):
        return [], [], []
    spec = STREAM_DEFAULTS.get(stream, {})
    if spec.get("split"):
        # The stream calls it dichro_positioner2; the scan knows its name.
        motors = read_scan_info(path)["motors"]
        return [motors[-1] if motors else spec["split"]], [spec["split"]], []
    axes = inner_axis_fields(path, stream)
    return (
        [motor for motor, _value, _setpoint in axes],
        [value for _motor, value, _setpoint in axes],
        [setpoint for _motor, _value, setpoint in axes],
    )


# How much wider the step from one mesh position to the next must be than
# the widest gap left inside one of them before the readings are read as a
# raster at all. A factor rather than a plain comparison, because the cuts
# are made at the widest gaps there are and so are always wider than what
# they leave behind: without one, twenty-five numbers with no levels in
# them whatever would gather into a five-by-five raster that was never
# measured.
#
# Where to draw it is measured from both sides, and the two are far apart.
# Numbers with no structure hardly reach two: an even ramp gives exactly
# one, and of 4000 random draws the worst was 2.6 over 25 points and 1.2
# over a hundred or more. Real meshes clear it by a long way -- scan 280 by
# 185, its wander a part in two hundred of the step -- but not limitlessly:
# scan 283, a 41-column raster over the same 0.2 mm, manages 15.6, and a
# finer one would manage less, the wander being the motor's and the step
# the scan's. Five leaves room on both sides.
GRID_LEVEL_MARGIN = 5.0


def _grid_levels(values, count) -> np.ndarray | None:
    """*values* gathered onto *count* levels, or ``None`` if they hold none.

    A reduced stream's positioner is an average of the readbacks taken while
    the point was being measured, so the readings meant to sit at one
    position differ in the fourth decimal and :func:`np.unique` counts
    twenty-three levels where the scan declares five (scan 280). Sorting
    and cutting at the
    *count* - 1 widest gaps gathers them again; each reading comes back as
    its own level's mean, which is as near the position as readbacks can
    say -- a shade off the setpoints a primary-stream map of the same scan
    is drawn against.
    """
    if count < 2 or count > len(values):
        return None
    order = np.argsort(values, kind="stable")
    gaps = np.diff(values[order])
    cuts = np.sort(np.argsort(gaps, kind="stable")[-(count - 1) :])
    within = np.delete(gaps, cuts)
    if len(within) and gaps[cuts].min() <= GRID_LEVEL_MARGIN * within.max():
        return None
    starts = np.zeros(len(values), dtype=bool)
    starts[cuts + 1] = True
    level = np.empty(len(values), dtype=int)
    level[order] = np.cumsum(starts)
    return (np.bincount(level, weights=values) / np.bincount(level))[level]


def _snap_grid_axes(slow, fast, count, shape) -> tuple | None:
    """Both axis arrays gathered onto the raster *shape* declares."""
    n_points = min(len(slow), len(fast), count)
    if len(shape) != 2 or shape[1] < 2:
        return None
    # The declared width says how many columns there are to find. How many
    # rows were reached is what the points themselves say, so a scan stopped
    # part way gathers into the rows it got to rather than the twenty-one it
    # set out to measure -- scan 283, interrupted four rows in -- and its
    # last row keeps whatever holes it really has.
    rows = _grid_levels(slow[:n_points], -(-n_points // shape[1]))
    columns = _grid_levels(fast[:n_points], shape[1])
    return None if rows is None or columns is None else (rows, columns)


def _build_grid(slow, fast, values, shape, slow_name, fast_name) -> dict | None:
    """Lay *values* onto the raster described by the two axis arrays, or ``None``."""
    n_points = min(len(slow), len(fast), len(values))
    if n_points < 4:
        return None
    slow, fast, values = slow[:n_points], fast[:n_points], values[:n_points]

    # return_inverse gives each point's row/column straight off, so a snaked
    # raster (alternate rows reversed) lands in the right place unaided.
    y_axis, rows = np.unique(slow, return_inverse=True)
    x_axis, columns = np.unique(fast, return_inverse=True)
    if len(x_axis) < 2 or len(y_axis) < 2:
        return None

    complete = len(x_axis) * len(y_axis) == n_points
    # A scan stopped part way has fewer rows than were asked for, so compare
    # against the declared raster rather than expecting a full rectangle: the
    # fast axis must be the declared width, the slow axis no taller.
    declared = (
        len(shape) == 2 and len(x_axis) == shape[1] and len(y_axis) <= shape[0]
    )
    if not (complete or declared):
        return None

    grid = np.full((len(y_axis), len(x_axis)), np.nan)
    grid[rows, columns] = values
    return {
        "x_name": fast_name,
        "y_name": slow_name,
        "x": x_axis,
        "y": y_axis,
        "z": grid,
        "points": n_points,
    }


def merge_grids(grids: list[dict]) -> dict:
    """Several rasters over the same two motors laid onto one, or the only one.

    A map stopped part way is finished by a second scan started where it
    left off -- scan 39 of Ryan_26-3/sro reached h = 0.985 of 1.3, and scan
    41 began again there, at half the step -- so the two belong on one pair
    of axes. Each axis takes every level any of the scans has; levels nearer
    than a hundredth of the finest step are one, which absorbs what two
    setpoint lists differ by in the last bits. A cell measured more than
    once is the mean of its finite values, and ``averaged`` counts them (the
    52 cells of the h = 0.985 row both scans measured). Uneven steps are
    left as they are: :func:`axis_edges` sizes each cell from its own
    neighbours.
    """
    if len(grids) == 1:
        return grids[0]
    steps = [np.diff(grid[key]) for grid in grids for key in ("x", "y")]
    positive = np.concatenate([step[step > 0] for step in steps])
    tolerance = 0.01 * positive.min() if positive.size else 0.0

    def union(key):
        levels = np.sort(np.concatenate([grid[key] for grid in grids]))
        return levels[np.r_[True, np.diff(levels) > tolerance]]

    def index(levels, values):
        return np.searchsorted(levels, values + tolerance, side="right") - 1

    x_axis, y_axis = union("x"), union("y")
    total = np.zeros((len(y_axis), len(x_axis)))
    count = np.zeros((len(y_axis), len(x_axis)), dtype=int)
    for grid in grids:
        rows, columns = np.meshgrid(
            index(y_axis, grid["y"]), index(x_axis, grid["x"]), indexing="ij"
        )
        finite = np.isfinite(grid["z"])
        np.add.at(total, (rows[finite], columns[finite]), grid["z"][finite])
        np.add.at(count, (rows[finite], columns[finite]), 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(count > 0, total / count, np.nan)
    return {
        "x_name": grids[0]["x_name"],
        "y_name": grids[0]["y_name"],
        "x": x_axis,
        "y": y_axis,
        "z": z,
        "points": sum(grid["points"] for grid in grids),
        "averaged": int(np.count_nonzero(count > 1)),
        "scans": len(grids),
    }


def read_grid(
    path: Path,
    value_field: str,
    monitor_field: str = "",
    stream: str = PRIMARY_STREAM,
    transmission: bool = False,
) -> dict | None:
    """*value_field* laid out over a two-motor grid scan, or ``None``.

    ``None`` means "not a raster": a two-motor ``inner_product`` scan moves both
    motors together, so its axes have as many levels as there are points and the
    product test in :func:`_build_grid` fails, which is exactly what tells the
    two plans apart without relying on the plan name.

    A reduced stream carries the same raster at one point per mesh position --
    a dichro mesh scan's XAS and XMCD maps -- under its own positional axis
    names, which the stream spec supplies. There it is read only for a scan
    that declares itself a mesh: a magneto-chiral scan leaves two positioners
    in the same stream, and its inner one is a setting, not a coordinate.
    Those axis fields are averaged readbacks with no setpoint beside them, so
    they are gathered back onto the declared raster first -- see
    :func:`_snap_grid_axes`.

    With *monitor_field* the cells are divided by it, or -- *transmission* --
    taken as ``log(monitor/value)`` instead, the two forms :func:`absorption`
    offers. The same choice the XMCD panels make, so a map of a transmission
    measurement reads in the same units its panels would.
    """
    info = read_scan_info(path)
    if len(info["motors"]) != 2:
        return None

    available = stream_fields(path, stream)
    axes = STREAM_DEFAULTS.get(stream, {}).get("grid_axes")
    if axes:
        if not is_mesh_scan(path):
            return None
        # One field per scanned axis, so there is nothing to try in turn; the
        # motors keep their own names for the labels.
        slow_options = [axes[0]] if axes[0] in available else []
        fast_options = [axes[1]] if axes[1] in available else []
        labels = tuple(info["motors"])
    else:
        slow_options = _grid_axis_candidates(available, info["motors"][0])
        fast_options = _grid_axis_candidates(available, info["motors"][1])
        labels = ()
    if not slow_options or not fast_options or value_field not in available:
        return None

    wanted = set(slow_options) | set(fast_options) | {value_field}
    if monitor_field:
        wanted.add(monitor_field)
    try:
        data = read_stream(path, list(wanted), stream)
    except OSError:
        return None
    if value_field not in data:
        return None

    values = data[value_field]
    monitor = data.get(monitor_field) if monitor_field else None
    if monitor is not None:
        length = min(len(values), len(monitor))
        values = absorption(values[:length], monitor[:length], transmission)

    # Readback first, setpoint only if the readback will not form a raster.
    for slow_name in slow_options:
        for fast_name in fast_options:
            if slow_name not in data or fast_name not in data:
                continue
            grid = _build_grid(
                data[slow_name],
                data[fast_name],
                values,
                info["shape"],
                labels[0] if labels else slow_name,
                labels[1] if labels else fast_name,
            )
            if grid is None and axes:
                # A reduced stream has no setpoint column to try next -- it
                # carries the two averaged positioners and nothing else --
                # so readings that will not form a raster as they stand are
                # gathered onto the one the scan declares instead. Second
                # rather than first for the same reason the setpoint is:
                # a reading that already lands on the raster is better
                # evidence of it than one binned there. The dichro mesh
                # scans whose positioners come out exact (156, 168) never
                # reach here; scan 280's, which wander, are mapped by it.
                snapped = _snap_grid_axes(
                    data[slow_name],
                    data[fast_name],
                    len(values),
                    info["shape"],
                )
                if snapped is not None:
                    grid = _build_grid(
                        snapped[0],
                        snapped[1],
                        values,
                        info["shape"],
                        labels[0],
                        labels[1],
                    )
            if grid is not None:
                return grid
    return None


def read_dichro_grid(
    path: Path,
    detector_field: str,
    monitor_field: str,
    stream: str = PRIMARY_STREAM,
    transmission: bool = False,
) -> list[tuple[str, dict]] | None:
    """XAS and XMCD maps from a two-motor mesh scan with polarization pairing.

    A dichro mesh steps the retarder between its two settings at every grid
    cell -- :data:`PR_FIELD` toggles, the magnet stays put -- so each cell
    carries one measurement per helicity. The absorption is the usual
    detector/monitor (or its log in transmission), averaged over the plus
    and minus retarder groups to give XAS, and differenced to give XMCD.
    This is :func:`reduce_dichro`'s split applied per cell instead of per
    energy, with the field axis dropped because a mesh scan has none to
    split on.

    Returns ``[("xas", grid), ("xmcd", grid)]`` in the shape
    :meth:`HdfViewer._plot_grids` already consumes, or ``None`` when the
    scan is not a dichro mesh -- no second motor, no retarder, only one
    retarder setting, or missing detector/monitor columns.
    """
    if not detector_field or not monitor_field:
        return None
    info = read_scan_info(path)
    if len(info["motors"]) != 2 or not is_mesh_scan(path):
        return None

    available = stream_fields(path, stream)
    axes = STREAM_DEFAULTS.get(stream, {}).get("grid_axes")
    if axes:
        slow_options = [axes[0]] if axes[0] in available else []
        fast_options = [axes[1]] if axes[1] in available else []
        labels = tuple(info["motors"])
    else:
        slow_options = _grid_axis_candidates(available, info["motors"][0])
        fast_options = _grid_axis_candidates(available, info["motors"][1])
        labels = ()
    if (
        not slow_options
        or not fast_options
        or PR_FIELD not in available
        or detector_field not in available
        or monitor_field not in available
    ):
        return None

    wanted = (
        set(slow_options)
        | set(fast_options)
        | {detector_field, monitor_field, PR_FIELD}
    )
    try:
        data = read_stream(path, list(wanted), stream)
    except OSError:
        return None
    pr = data.get(PR_FIELD)
    detector = data.get(detector_field)
    monitor = data.get(monitor_field)
    if pr is None or detector is None or monitor is None:
        return None
    size = min(pr.size, detector.size, monitor.size)
    if size < 4:
        return None

    retarders = collect_uniques(pr[:size], PR_ATOL)
    if retarders.size != 2:
        return None
    # Same test as :func:`reduce_dichro`: "plus" is the larger retarder
    # readback, so the two halves of the split are named the same way the
    # 1D reduction names them.
    pr_plus = np.isclose(pr[:size], retarders.max(), atol=PR_ATOL)
    mu = absorption(detector[:size], monitor[:size], transmission)

    # Try the readback first, then its setpoint, the same order
    # :func:`read_grid` tries them in -- a readback that already lands on
    # the raster is better evidence of it than one gathered onto it.
    for slow_name in slow_options:
        for fast_name in fast_options:
            if slow_name not in data or fast_name not in data:
                continue
            slow = data[slow_name][:size]
            fast = data[fast_name][:size]
            result = _dichro_grid_from_cells(
                slow,
                fast,
                mu,
                pr_plus,
                info["shape"],
                labels[0] if labels else slow_name,
                labels[1] if labels else fast_name,
            )
            if result is not None:
                return result
    return None


def _dichro_grid_from_cells(
    slow: np.ndarray,
    fast: np.ndarray,
    mu: np.ndarray,
    pr_plus: np.ndarray,
    shape: tuple,
    slow_name: str,
    fast_name: str,
) -> list[tuple[str, dict]] | None:
    """Group *mu* by cell and retarder half, and lay the result on a raster.

    The grid is read off the same ``np.unique(..., return_inverse=True)``
    pair :func:`_build_grid` uses, so a snaking raster lands in the right
    place without an explicit un-snake.
    """
    # Reusing :func:`_build_grid`'s rectangle test rather than copying it:
    # it already knows how to accept a scan stopped part way through, and
    # the test is the same for either half of the split.
    y_axis, rows = np.unique(slow, return_inverse=True)
    x_axis, columns = np.unique(fast, return_inverse=True)
    if len(x_axis) < 2 or len(y_axis) < 2:
        return None
    # A dichro mesh visits every cell twice, so the point count is twice
    # the cell count rather than equal to it. A scan stopped part way
    # through still has to look like one, so compare against the declared
    # raster the same way :func:`_build_grid` does.
    pairs = len(mu) // 2
    complete = len(x_axis) * len(y_axis) == pairs
    declared = (
        len(shape) == 2 and len(x_axis) == shape[1] and len(y_axis) <= shape[0]
    )
    if not (complete or declared):
        return None

    plus_sum = np.zeros((len(y_axis), len(x_axis)))
    plus_count = np.zeros_like(plus_sum)
    minus_sum = np.zeros_like(plus_sum)
    minus_count = np.zeros_like(plus_sum)
    finite = np.isfinite(mu)
    plus = finite & pr_plus
    minus = finite & ~pr_plus
    # np.add.at is the right primitive here: a cell with two plus-samples
    # (one scan looped twice) sums both, where plus_sum[rows, cols] = mu
    # would keep only the last.
    np.add.at(plus_sum, (rows[plus], columns[plus]), mu[plus])
    np.add.at(plus_count, (rows[plus], columns[plus]), 1)
    np.add.at(minus_sum, (rows[minus], columns[minus]), mu[minus])
    np.add.at(minus_count, (rows[minus], columns[minus]), 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        plus_mean = np.where(plus_count > 0, plus_sum / plus_count, np.nan)
        minus_mean = np.where(minus_count > 0, minus_sum / minus_count, np.nan)
    xas = (plus_mean + minus_mean) / 2.0
    xmcd = plus_mean - minus_mean
    if not (np.isfinite(xas).any() and np.isfinite(xmcd).any()):
        return None

    def wrap(z):
        return {
            "x_name": fast_name,
            "y_name": slow_name,
            "x": x_axis,
            "y": y_axis,
            "z": z,
            "points": int(pairs),
        }

    return [("xas", wrap(xas)), ("xmcd", wrap(xmcd))]


def axis_edges(values: np.ndarray) -> np.ndarray:
    """Cell boundaries around the grid's sample positions.

    Midway between neighbours, with the outer two extended by half their own
    spacing, so uneven steps keep their proportions.
    """
    v = np.asarray(values, dtype=float)
    if v.size == 1:
        return np.array([v[0] - 0.5, v[0] + 0.5])
    middles = (v[1:] + v[:-1]) / 2.0
    return np.concatenate(
        ([2 * v[0] - middles[0]], middles, [2 * v[-1] - middles[-1]])
    )


def clipped_limits(
    z: np.ndarray, percent: float = COLOR_CLIP_PERCENT
) -> tuple[float | None, float | None, str]:
    """Colour limits leaving out the most extreme *percent* at each end.

    A handful of wild cells -- a spot where the beam ran off the sample, a
    detector glitch -- stretch a map's colour scale over a range nothing else
    occupies, leaving the rest of it one flat shade. Cutting the scale back to
    where the bulk of the data lives brings that structure back; the cells
    outside simply saturate.

    How much to cut is asked for rather than worked out, which is the point.
    Rules that decide for themselves what is an outlier (Tukey's fences on the
    interquartile range, say) read the top of a sharp peak as one, because a
    peaked map legitimately has a long tail above its upper quartile. On the
    sample-position rasters this beamline lines up on -- ``huber_euler_z`` by
    ``huber_euler_y`` against a detector total, where the peak *is* the result
    -- those fences clip most of the peak away: measured over the twelve of
    them, 3x IQR left scan 468 with 7% of its height, 232 with 14%, 228 with
    32%. A number someone typed can be wrong, but it cannot be surprising,
    and the same number gives the same answer on every map it is applied to.

    *percent* is a percentile, so what it costs depends on how the values are
    distributed, not on how extreme they are: 0.5 on a 1500-cell raster is the
    eight cells at each end whatever they hold. Nothing stops it from being
    large -- a map still washed out after a fraction of a percent may well
    need several -- but past a few percent it is no longer only spikes going.

    Returns the two limits and which ends they actually cut, in the wording
    :meth:`Figure.colorbar` wants for its ``extend`` arrows, so a clipped
    scale says on the bar that it is clipped. ``(None, None, "neither")``
    means leave the limits alone.
    """
    finite = np.asarray(z, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < 3:
        return None, None, "neither"
    low, high = (
        float(v) for v in np.percentile(finite, [percent, 100.0 - percent])
    )
    if not (np.isfinite(low) and np.isfinite(high)) or low >= high:
        # Flat, or so nearly flat that the percentiles coincide.
        return None, None, "neither"
    return low, high, scale_extend(finite, low, high)


def lone_dropouts(monitor: np.ndarray) -> np.ndarray:
    """Mask the points where the beam was gone for that point alone.

    A dropout is a point whose *monitor* has collapsed to less than
    :data:`BEAM_GONE_FRACTION` of the column's own median. Every channel goes
    with it -- on scan 239's row 367 the monitor reads 5372 against a median
    of 1.85e6, ``LockDC`` 5502 against 713000 -- so the point carries no
    measurement at all, and in the reduction it drags a whole energy bin with
    it: 239's XMuCD spans 0.511 with that row in and 0.0262 with it out.

    Unlike :func:`clipped_limits` above this needs no amount typed into it,
    and the difference is worth stating because the two sit side by side. The
    clip is a judgement -- where a peak stops being signal is not written in
    the data. This is not: a monitor measures the beam arriving, and a beam
    that is arriving at 2% of normal has gone. The two populations are forty
    apart with nothing between them, so the threshold is a reading rather than
    a choice.

    Only *isolated* dropouts are masked -- both neighbours must be live. A run
    of them is not a spike but a dump, and the scan should show it: of the 24
    scans in these directories carrying a dead point, 17 carry nothing but
    single ones and 7 were taken while the beam was down, scan 244 for 265
    consecutive rows. Deleting that would hide a quarter of the run behind a
    curve that looked fine. Length, not depth, is what tells the two apart --
    scan 89's one dead row is 5% of a 20-point scan and still a spike.

    Neighbours are the rows either side *in acquisition order*, so this wants
    one file's column at a time: run two files together and the join would
    stand one scan's last row next to another's first.

    What this deliberately will not catch is a spike the monitor knew nothing
    about. Reduced index 80 of scans 239 and 242 (E = 7.26907) is a 13-sigma
    one-point excursion with the monitor at 0.91 of median and all four raw
    rows agreeing to four decimals -- a monochromator glitch, reproduced in
    every repeat. It is data, and it stays.
    """
    values = np.asarray(monitor, dtype=float)
    gone = np.zeros(values.shape, dtype=bool)
    if values.ndim != 1 or values.size < 3:
        return gone
    finite = values[np.isfinite(values)]
    if finite.size < 3:
        return gone
    median = float(np.median(finite))
    if not np.isfinite(median) or median <= 0.0:
        # No scale to judge against: a monitor that is flat zero, or counting
        # down, says nothing about which of its points had beam.
        return gone
    dead = np.isfinite(values) & (values < BEAM_GONE_FRACTION * median)
    padded = np.concatenate(([False], dead, [False]))
    return dead & ~padded[:-2] & ~padded[2:]


def scale_extend(z: np.ndarray, low: float | None, high: float | None) -> str:
    """Which ends of *z* fall outside ``(low, high)``, as colour-bar arrows.

    In the wording :meth:`Figure.colorbar` wants for ``extend``, so a scale
    that does not reach its own data says so on the bar rather than leaving it
    to be noticed. ``None`` for a limit is matplotlib's "fit the data", which
    by definition cuts nothing off at that end.

    Asked separately from the limits because the limits move after they are
    first worked out: :func:`centered_limits` widens a range, and cells that
    saturated against the clipped one need not saturate against the wider. An
    arrow left over from the narrower scale points at data that is not out
    there.
    """
    finite = np.asarray(z, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return "neither"
    cuts = (
        low is not None and low > float(finite.min()),
        high is not None and high < float(finite.max()),
    )
    return {
        (True, True): "both",
        (True, False): "min",
        (False, True): "max",
        (False, False): "neither",
    }[cuts]


def centered_limits(
    z: np.ndarray, low: float | None, high: float | None
) -> tuple[float | None, float | None]:
    """``(low, high)`` widened until zero sits exactly in the middle.

    A diverging colour map is read by where its middle colour falls, and left
    to the data that lands wherever the range happens to centre. On scan 156's
    XMCD the clipped range is -0.0009 to +0.0036, which puts white at +0.0013:
    a map with signal everywhere, drawn so that most of it reads as none. Zero
    is where white belongs, because zero is what the map is a difference from.

    Widening to whichever end reaches further costs the shorter end some of
    its contrast, and buys a map that can be read without consulting the bar.
    Either limit may be ``None``, meaning the data's own extreme.

    Only ever widens, so no cell that was inside the scale leaves it -- but
    the arrows still have to be recomputed, which is :func:`scale_extend`.
    """
    finite = np.asarray(z, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return low, high
    lo = float(finite.min()) if low is None else low
    hi = float(finite.max()) if high is None else high
    reach = max(abs(lo), abs(hi))
    # Every cell is zero, and there is no scale left to centre on it.
    return (-reach, reach) if reach > 0 else (low, high)


def grid_stats(grid: dict) -> dict[str, float | None]:
    """Max/min of a grid and the (x, y) they sit at."""
    z = grid["z"]
    stats = {"max_z": None, "max_at": None, "min_z": None, "min_at": None}
    if not np.any(np.isfinite(z)):
        return stats
    flat_max = int(np.nanargmax(z))
    flat_min = int(np.nanargmin(z))
    for key, flat in (("max", flat_max), ("min", flat_min)):
        row, column = np.unravel_index(flat, z.shape)
        stats[f"{key}_z"] = float(z[row, column])
        stats[f"{key}_at"] = (float(grid["x"][column]), float(grid["y"][row]))
    return stats


def read_stream(
    path: Path, fields: list[str], stream: str = PRIMARY_STREAM
) -> dict[str, np.ndarray]:
    """Read the ``value`` array of each of *fields* from one stream.

    Fields the stream does not hold are left out of the result rather than
    raising, so a selection spanning files that do not all record the same
    thing -- or all the same streams -- still plots what it can.
    """
    out: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as f:
        group = f.get(f"{STREAMS_PATH}/{stream}")
        if group is None:
            return out
        epoch_source = _epoch_source(group)
        for name in fields:
            if name == EPOCH_FIELD:
                if epoch_source is not None:
                    out[name] = np.asarray(
                        group[epoch_source]["EPOCH"][()], dtype=float
                    )
                continue
            item = group.get(name)
            if item is not None and "value" in item:
                out[name] = np.asarray(item["value"][()], dtype=float)
    return out


def mesh_step_count(path: Path, points: int) -> int:
    """How many of a stream's points belong to one row of the scan's mesh.

    A ``grid_scan`` records its raster size in the metadata, innermost axis
    last, and the reduced stream of a magneto-chiral scan holds exactly that
    mesh -- one point per (outer position, inner step). Returns ``1`` for
    anything that is not a mesh of more than one row, which is everything
    there is nothing to split.

    Every axis inside the outermost counts, not just the last: a scan that
    steps a polarization and cycles a magnet current within it comes round to
    its starting setting only after the product of the two.
    """
    shape = read_scan_info(path)["shape"]
    if len(shape) < 2:
        return 1
    steps = int(np.prod(shape[1:]))
    return steps if 1 < steps < points else 1


def step_groups(
    columns: list[np.ndarray], steps: int, atol: float, average: bool
) -> list[tuple[tuple[float, ...], list[int]]]:
    """The curves a mesh's inner axes fall into: a value each, and the steps at it.

    One group per step by default, so a magnet current cycled over
    ``[1, -1, -1, 1]`` gives four curves and the two measurements of each field
    stay apart. With *average* the steps whose readbacks agree to within *atol*
    are gathered instead, and the same pattern gives one curve per field, ``+I``
    and ``-I``. None of this compares readbacks for equality: a step stands at
    the *median* of what it read, which the rows a supply spent ramping or
    settling to its set point cannot drag off the value it held for the rest.

    One *column* per axis, so a step is a position in all of them at once and
    a value is a tuple. Gathering needs every axis to agree: the polarization
    and the magnet of a fixed-polarization scan repeat on different periods,
    and two steps sharing a current but not a polarization are two
    measurements, not a repeat of one.
    """
    values = []
    for step in range(steps):
        position = []
        for column in columns:
            sample = column[step::steps]
            sample = sample[np.isfinite(sample)]
            position.append(
                float(np.median(sample)) if sample.size else float("nan")
            )
        values.append(tuple(position))
    if not average:
        return [(value, [step]) for step, value in enumerate(values)]
    gathered: list[list[int]] = []
    for step, value in enumerate(values):
        for group in gathered:
            if all(
                abs(own - theirs) <= atol
                for own, theirs in zip(value, values[group[0]])
            ):
                group.append(step)
                break
        else:  # no group near enough, so this step starts one
            gathered.append([step])
    return [
        (
            tuple(
                float(np.median([values[step][axis] for step in group]))
                for axis in range(len(columns))
            ),
            group,
        )
        for group in gathered
    ]


def stack_mean(pieces: list[np.ndarray], rows: int) -> np.ndarray:
    """The mean of *pieces* row by row over *rows*, counting only the finite.

    A measurement stopped part way leaves its last rows short, so a row only
    some of the pieces reached is the mean of those and a row none of them
    reached comes back NaN -- never a short piece quietly read as zero.
    Shared by the repeats within one scan and the repeats across several,
    which are the same reduction asked for twice.
    """
    stack = np.full((len(pieces), rows), np.nan)
    for row, piece in enumerate(pieces):
        taken = np.asarray(piece, dtype=float)
        stack[row, : taken.size] = taken
    counted = np.count_nonzero(np.isfinite(stack), axis=0)
    return np.where(
        counted > 0, np.nansum(stack, axis=0) / np.maximum(counted, 1), np.nan
    )


def group_curve(
    x: np.ndarray, y: np.ndarray, steps: int, members: list[int]
) -> tuple[np.ndarray, np.ndarray]:
    """One group's curve: the points its *members* took, one per mesh row.

    A group covering a single step is that step's slice of the stream. One
    covering several -- the repeats within a current pattern -- is their mean
    row by row, the same reduction ``process_data.process_xmucd`` makes of the
    same measurement. A scan stopped part way leaves its last row short, so
    rows a step never reached are left out of the mean rather than counted.
    """
    if len(members) == 1:
        step = members[0]
        return x[step::steps], y[step::steps]
    rows = -(-max(x.size, y.size) // steps)
    return (
        stack_mean([x[step::steps] for step in members], rows),
        stack_mean([y[step::steps] for step in members], rows),
    )


def collect_uniques(values: np.ndarray, atol: float) -> np.ndarray:
    """The distinct settings in *values*, counting any two within *atol* as one.

    ``process_data.collect_uniques`` verbatim. A readback never repeats
    exactly, so the settings have to be found by distance; the first value of
    each cluster is the one kept, which is why the result follows the sorted
    order rather than the order measured.
    """
    unique = np.unique(values)
    kept = [unique[0]]
    for value in unique[1:]:
        if not np.any(np.isclose(value, kept, atol=atol)):
            kept.append(value)
    return np.array(kept)


def absorption(
    detector: np.ndarray, monitor: np.ndarray, transmission: bool
) -> np.ndarray:
    """``log(monitor/detector)`` or ``detector/monitor``, NaN where not finite.

    A dead I0 gives an infinity rather than an exception, and treating that as
    missing drops it out of a group's mean instead of poisoning the mean with
    it. Every point coming back NaN is how a dead channel is recognised --
    see :meth:`HdfViewer._contrast_message`.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.log(monitor / detector) if transmission else detector / monitor
    return np.where(np.isfinite(mu), mu, np.nan)


def _empty_contrasts(mode: str) -> dict:
    """What a reduction returns when it could make nothing of its columns.

    Every key the full result has, so a caller reading ``n_steps`` to say why
    the panels are empty need not first ask whether there is anything there.
    ``_empty_result`` in xmucd_reduce.py.
    """
    return {
        "energy": np.array([]),
        **{key: np.array([]) for key in CONTRASTS},
        "mode": mode,
        "sources": CONTRAST_SOURCES[mode],
        "n_points": 0,
        "n_finite": 0,
        "n_energies": 0,
        "n_complete": 0,
        "n_steps": 0,
    }


def reduce_dichro(
    energy: np.ndarray,
    pr: np.ndarray,
    current: np.ndarray,
    detector: np.ndarray,
    monitor: np.ndarray,
    *,
    transmission: bool = True,
) -> dict:
    """Split raw per-point columns into XAS, XMCD, XMuCD and XNCD.

    ``xmucd_reduce.reduce_columns`` transcribed -- see the note above
    :data:`CONTRASTS` for why it is copied rather than imported. The model is
    ``mu(p, B) = A + C*p + D*B + E*p*B`` with the helicity *p* and the field
    *B* both +/-1, so grouping one energy's points into the four (p, B)
    quadrants separates the terms: XAS is A, XMCD 2E, XMuCD 2D and XNCD 2C.

    "Plus" is the largest magnet current and the largest retarder readback.
    The split stays binary however many settings the readback bins into -- a
    ramping kepco gives five on scan 96 -- because what is being asked is
    which side of the measurement a point is on, not which setting it sat at.

    Returns ``energy`` and one array per entry of :data:`CONTRASTS`, covering
    only the energies whose four groups all hold a finite point, plus
    ``n_points`` (points consumed), ``n_finite`` (of those, how many gave a
    usable absorption -- zero means the detector or the monitor is dead),
    ``n_energies``, ``n_complete`` and ``n_steps`` (distinct field settings).
    """
    columns = [
        np.asarray(column, dtype=float).ravel()
        for column in (energy, pr, current, detector, monitor)
    ]
    # A scan stopped part way leaves its last point half written, which is a
    # column or two shorter than the rest rather than an error.
    size = min(column.size for column in columns)
    empty = _empty_contrasts(MODE_DICHRO)
    if size == 0:
        return empty
    energy, pr, current, detector, monitor = (
        column[:size] for column in columns
    )

    mu = absorption(detector, monitor, transmission)
    on_grid = np.isfinite(energy) & np.isfinite(pr) & np.isfinite(current)
    if not on_grid.any():
        return empty
    energy, pr, current, mu = (
        column[on_grid] for column in (energy, pr, current, mu)
    )

    energies = collect_uniques(energy, ENERGY_ATOL)
    retarders = collect_uniques(pr, PR_ATOL)
    fields = collect_uniques(current, SPLIT_ATOL)
    field_plus = np.isclose(current, fields.max(), atol=SPLIT_ATOL)
    pr_plus = np.isclose(pr, retarders.max(), atol=PR_ATOL)

    out = {key: [] for key in ("energy",) + CONTRASTS}
    for value in energies:
        here = np.isclose(energy, value, atol=ENERGY_ATOL)
        groups = (
            mu[here & pr_plus & field_plus],
            mu[here & pr_plus & ~field_plus],
            mu[here & ~pr_plus & field_plus],
            mu[here & ~pr_plus & ~field_plus],
        )
        # Part way through a scan only some of the grid has been visited, and
        # an energy missing any one of the four quadrants cannot be reduced at
        # all -- every contrast below needs all four.
        if any(not np.isfinite(group).any() for group in groups):
            continue
        pp, pm, mp, mm = (float(np.nanmean(group)) for group in groups)
        out["energy"].append(float(value))
        out["xas"].append(np.mean([pp, pm, mp, mm]))
        out["xmcd"].append(np.mean([pp, mm]) - np.mean([pm, mp]))
        out["xmucd"].append(np.mean([pp, mp]) - np.mean([pm, mm]))
        out["xncd"].append(np.mean([pp, pm]) - np.mean([mp, mm]))

    result = {key: np.array(values) for key, values in out.items()}
    result["mode"] = MODE_DICHRO
    result["sources"] = CONTRAST_SOURCES[MODE_DICHRO]
    result["n_points"] = int(size)
    result["n_finite"] = int(np.isfinite(mu).sum())
    result["n_energies"] = int(energies.size)
    result["n_complete"] = int(result["energy"].size)
    result["n_steps"] = int(fields.size)
    return result


def reduce_lockin(
    energy: np.ndarray,
    current: np.ndarray,
    detector: np.ndarray,
    ac: np.ndarray,
    monitor: np.ndarray,
    *,
    transmission: bool = True,
) -> dict:
    """Split a lock-in scan's columns into XAS, XMCD, XMuCD and XNCD.

    ``xmucd_reduce.reduce_lockin_columns`` transcribed, in its
    :data:`MODE_LOCKIN_POL` sense -- see the note above :data:`CONTRASTS` for
    why it is copied rather than imported, and above :data:`MODE_DICHRO` for
    why the other lock-in sense is not here.

    The lock-in oscillates the retarder, so the helicity difference is
    measured directly on *ac* and only the field is stepped: each energy has
    two groups rather than :func:`reduce_dichro`'s four. Writing the field as
    s = +/-1::

        dc(s) = A + D*s                 ac(s) = k*(C + E*s)
        xas   = mean([dc(+), dc(-)])    xncd  = ac(+) + ac(-)
        xmucd = dc(+) - dc(-)           xmcd  = ac(+) - ac(-)

    *ac* is used exactly as handed in -- the caller has already subtracted
    :data:`ACOFF_FIELD` -- so it carries the lock-in's own units and gain,
    with neither the monitor nor *transmission* applied. The AC-derived
    panels therefore cannot be read against the DC-derived ones by magnitude;
    the returned ``sources`` says which is which.

    Returns what :func:`reduce_dichro` returns. ``n_steps`` must reach 2
    before any energy completes, so a single file at one field gives nothing
    -- which is the usual case here, the field being stepped between files.
    """
    columns = [
        np.asarray(column, dtype=float).ravel()
        for column in (energy, current, detector, ac, monitor)
    ]
    size = min(column.size for column in columns)
    empty = _empty_contrasts(MODE_LOCKIN_POL)
    if size == 0:
        return empty
    energy, current, detector, ac, monitor = (
        column[:size] for column in columns
    )

    mu = absorption(detector, monitor, transmission)
    # An absent LockAC or LockACoff arrives as an all-NaN column rather than
    # an error, and has to stay that way: it means there is no lock-in data
    # here, which the status bar reports, not that the difference was zero.
    ac = np.where(np.isfinite(ac), ac, np.nan)
    on_grid = np.isfinite(energy) & np.isfinite(current)
    if not on_grid.any():
        return empty
    energy, current, mu, ac = (
        column[on_grid] for column in (energy, current, mu, ac)
    )

    energies = collect_uniques(energy, ENERGY_ATOL)
    fields = collect_uniques(current, SPLIT_ATOL)
    # "Plus" is the largest current, as in the dichro reduction. With one
    # setting every point is plus, no energy has the pair it needs, and the
    # loop below yields nothing.
    field_plus = np.isclose(current, fields.max(), atol=SPLIT_ATOL)

    out = {key: [] for key in ("energy",) + CONTRASTS}
    for value in energies:
        here = np.isclose(energy, value, atol=ENERGY_ATOL)
        plus, minus = here & field_plus, here & ~field_plus
        groups = (mu[plus], mu[minus], ac[plus], ac[minus])
        if any(not np.isfinite(group).any() for group in groups):
            continue
        dc_p, dc_m, ac_p, ac_m = (float(np.nanmean(group)) for group in groups)
        out["energy"].append(float(value))
        out["xas"].append(np.mean([dc_p, dc_m]))
        out["xmcd"].append(ac_p - ac_m)
        out["xmucd"].append(dc_p - dc_m)
        out["xncd"].append(ac_p + ac_m)

    result = {key: np.array(values) for key, values in out.items()}
    result["mode"] = MODE_LOCKIN_POL
    result["sources"] = CONTRAST_SOURCES[MODE_LOCKIN_POL]
    result["n_points"] = int(size)
    # Both channels have to be alive for a lock-in scan to reduce, so count
    # where they both are: a live LockDC beside a dead LockAC otherwise reads
    # as a healthy scan right up until no energy ever completes.
    result["n_finite"] = int((np.isfinite(mu) & np.isfinite(ac)).sum())
    result["n_energies"] = int(energies.size)
    result["n_complete"] = int(result["energy"].size)
    result["n_steps"] = int(fields.size)
    return result


def xncd_partner(result: dict) -> str:
    """The contrast XNCD's leakage is read against: the other one off its channel.

    XNCD is even in the field, so the sample cannot make it and whatever
    shows there is the helicity average failing to cancel. The same leakage
    sits in exactly one other panel -- XMuCD in the dichro reduction, where
    both come off the DC channel, and XMCD in the lock-in one, where both
    come off AC. Comparing across the two would set an absorption against a
    raw lock-in difference, and the ratio would mean nothing.
    """
    sources = result["sources"]
    return next(
        key for key in ("xmucd", "xmcd") if sources[key] == sources["xncd"]
    )


def _scalarize(value: np.ndarray):
    """Return a plain float for a numeric size-1 array, else a short text summary.

    Some baseline fields are waveform-valued (e.g. a UB matrix or a
    pixel-center pair) rather than scalar, so ``value_start``/``value_end``
    are not always convertible with a bare ``float()``. Others are size-1
    but string-valued (e.g. a motor direction like ``b"y-"``), which
    ``float()`` also rejects.
    """
    value = np.asarray(value)
    if value.size == 1:
        item = value.reshape(-1)[0]
        if isinstance(item, bytes):
            return item.decode(errors="replace")
        try:
            return float(item)
        except (TypeError, ValueError):
            return str(item)
    return f"array{value.shape}: {np.array2string(value, precision=4, threshold=6)}"


def read_baseline(path: Path) -> dict[str, tuple[object, object]]:
    """Return ``{field: (before, after)}`` from the baseline stream.

    ``before``/``after`` are ``float`` for a scalar field, or a text summary
    for an array-valued one (see :func:`_scalarize`).
    """
    out: dict[str, tuple[object, object]] = {}
    with h5py.File(path, "r") as f:
        group = f.get(BASELINE_PATH)
        if group is None:
            return out
        for name in group.keys():
            sub = group[name]
            if not isinstance(sub, h5py.Group):
                continue
            if "value_start" in sub and "value_end" in sub:
                before = _scalarize(sub["value_start"][()])
                after = _scalarize(sub["value_end"][()])
            elif (
                "value" in sub
                and sub["value"].shape
                and sub["value"].shape[0] >= 2
            ):
                value = np.asarray(sub["value"][()])
                before, after = _scalarize(value[0]), _scalarize(value[-1])
            else:
                continue
            out[name] = (before, after)
    return out


def _baseline_numbers(path: Path) -> dict[str, tuple[float, float]]:
    """Every numeric baseline field of *path*, as its before/after pair.

    Numeric fields only: a baseline holds motor directions and UB matrices
    alongside the readbacks, and :func:`_scalarize` hands those back as text.
    A field the scan did not record, or recorded as text, is left out rather
    than filled in, so a caller can tell "not measured" from a value.

    Remembered per (path, mtime, size) the way :func:`reduce_channels` is,
    and the whole baseline at once rather than the fields asked for:
    :func:`read_baseline` reads all 322 of them whatever is wanted, three
    callers here want three different handfuls, and every redraw asks again
    for every selected scan. Keying on the field list too would open each
    file three times over for one answer. ``OSError`` gives ``{}``, which
    leaves the scan ungroupable rather than taking the window down.

    Both ends are kept rather than averaged away because the gap between
    them is the only noise estimate one scan carries on its own, which is
    what :func:`baseline_trend_fields` ranks by.
    """
    path = Path(path)
    try:
        stat = path.stat()
        key = (str(path), stat.st_mtime, stat.st_size)
    except OSError:
        key = None
    if key is not None and key in _BASELINE_CACHE:
        return _BASELINE_CACHE[key]
    try:
        baseline = read_baseline(path)
    except OSError:
        baseline = {}
    numbers: dict[str, tuple[float, float]] = {}
    for name, pair in baseline.items():
        try:
            numbers[name] = (float(pair[0]), float(pair[1]))
        except (TypeError, ValueError):
            continue
    if key is not None:
        _BASELINE_CACHE[key] = numbers
    return numbers


def baseline_values(path: Path, fields: tuple[str, ...]) -> dict[str, float]:
    """The mean of each named baseline field's before and after reading.

    The mean of the two because they are one quantity read twice, before and
    after; the sample does not move during a scan and the cryostat drifts
    under a kelvin, so the midpoint is the setting better than either end.

    Numeric fields only, and a field the scan did not record is left out
    rather than filled in -- :func:`_baseline_numbers` does the reading and
    the remembering, and says why.
    """
    numbers = _baseline_numbers(path)
    return {
        name: (numbers[name][0] + numbers[name][1]) / 2.0
        for name in fields
        if name in numbers
    }


def baseline_trend_fields(paths: list[Path]) -> list[str]:
    """Baseline fields worth trending against, most deliberate first.

    The answer to "what was stepped across these scans": every numeric field
    the whole selection recorded, whose motion over the selection clears the
    noise bar BASELINE_TREND_SNR and BASELINE_TREND_NOISE describe, best
    first. Fewer than two scans gives ``[]`` -- one point has nothing to be
    plotted against.

    Ordered by how far each moved relative to its own noise rather than
    alphabetically, so the motor somebody deliberately drove comes out at
    the top: scans 32 to 39 put huber_euler_ana_eta third of 19, behind only
    its own setpoint and the huber_hp alias that reads the same encoder.

    The ranking is of one selection against itself, so it says nothing about
    a wider one. Over all 146 TiSe2 scans at once eta falls below the bar
    altogether, correctly -- somewhere in those 146 it was moving *during* a
    scan, and a field that will not hold still inside one scan cannot be the
    axis that scan sits at.
    """
    paths = list(paths)
    if len(paths) < 2:
        return []
    # Keyed on the same (path, mtime, size) the readings underneath are, and
    # not on the names alone: a scan still being written to is one the
    # listing offers like any other, and on the names alone it would keep
    # whatever ranking its first few points happened to earn. Worth a cache
    # at all because the ranking is 14 ms over a 27 scan selection and every
    # redraw asks for it again, where the readings behind it are already
    # remembered and cost nothing the second time.
    stamps = []
    for name in sorted(str(path) for path in paths):
        try:
            stat = Path(name).stat()
            stamps.append((name, stat.st_mtime, stat.st_size))
        except OSError:
            stamps.append((name, None, None))
    key = tuple(stamps)
    if key in _TREND_FIELD_CACHE:
        return _TREND_FIELD_CACHE[key]
    pairs = [_baseline_numbers(path) for path in paths]
    common = set(pairs[0]).intersection(*pairs[1:]) if pairs else set()
    scored: list[tuple[float, str]] = []
    for name in common:
        readings = [reading[name] for reading in pairs]
        if not all(
            np.isfinite(before) and np.isfinite(after)
            for before, after in readings
        ):
            continue
        values = [(before + after) / 2.0 for before, after in readings]
        across = max(values) - min(values)
        if across <= 0.0:
            continue
        within = max(abs(after - before) for before, after in readings)
        floor = within + BASELINE_TREND_NOISE * max(
            abs(value) for value in values
        )
        if floor <= 0.0 or across / floor > BASELINE_TREND_SNR:
            scored.append((across / floor if floor > 0.0 else np.inf, name))
    # By score, then by name, so that the aliases reading one encoder -- and
    # they are legion here -- come out in a stable order rather than
    # whichever the set happened to yield.
    scored.sort(key=lambda entry: (-entry[0], entry[1]))
    fields = [name for _score, name in scored]
    _TREND_FIELD_CACHE[key] = fields
    return fields


def lockin_groups(paths: list[Path]) -> list[list[Path]]:
    """Gather lock-in scans into the sets that reduce together.

    A lock-in scan needs two field settings and a single-pass one has only
    the one: 235 sits at +30.02 A and 236 at -29.96 A, and they are a single
    measurement written to two files. So the field is stepped *between* files
    and the reduction has to join them, which is only sound where the files
    are the same measurement -- same spot on the sample, same temperature.
    That is the key, read from each scan's baseline through
    :data:`GROUP_FIELDS`. A scan that stepped the field inside itself (237,
    238) is left on its own: it already has both settings, and a partner
    would only dilute it.

    Each quantity is clustered across the selection with
    :func:`collect_uniques` at its own tolerance and every scan assigned to
    the nearest surviving value, so the grouping does not depend on which
    scan was looked at first. The tolerances are the measured gaps: within a
    group the positions agree to four decimals and the temperature drifts by
    up to 0.91 K, while between groups the positions differ by 0.02 mm and
    0.5 deg and the temperature by 19.8 K.

    The field itself is deliberately not in the key. The two halves of one
    pair sit at +30.02 and -29.96 A, so binning on the magnitude would split
    the very pair being looked for.

    The temperature tolerance cannot go below the drift it has to survive, so
    a series stepped finer than 2 K would chain into one group. Nothing here
    can tell that from a real pair; the status bar names each group's
    temperature so it can be seen, and selecting fewer scans splits it.

    Groups come back in the order their first scan appears, which is the
    table's own order.
    """
    pooled = [
        path for path in paths if MAGNET_FIELD not in reduce_channels(path)
    ]
    readings = {
        str(path): baseline_values(path, tuple(GROUP_FIELDS)) for path in pooled
    }
    settings: dict[str, np.ndarray] = {}
    for name, atol in GROUP_FIELDS.items():
        measured = [
            reading[name] for reading in readings.values() if name in reading
        ]
        if measured:
            settings[name] = collect_uniques(np.array(measured), atol)

    groups: dict[object, list[Path]] = {}
    for path in paths:
        reading = readings.get(str(path), {})
        if any(name not in reading for name in GROUP_FIELDS):
            # Either the field is already stepped inside this file, or the
            # baseline does not say where the scan was. Neither can be
            # pooled, and a scan reduced alone reports itself.
            key: object = ("alone", str(path))
        else:
            key = tuple(
                int(np.argmin(np.abs(settings[name] - reading[name])))
                for name in GROUP_FIELDS
            )
        groups.setdefault(key, []).append(path)
    return list(groups.values())


def scan_signature(path: Path) -> tuple | None:
    """What a scan has to share with another before it can be a repeat of it.

    The plan and the motors it stepped, taken from the metadata rather than
    the columns: two scans that moved different things are not repeats of one
    another however alike their numbers look. ``None`` where the scan stepped
    nothing -- a ``count`` has no shape to repeat, and the 32 of them in one
    directory here would otherwise read as a single measurement taken 32
    times, their row counts running from 1 to 600.
    """
    info = read_scan_info(path)
    motors = tuple(name for name, *_ in info["ranges"])
    return (info["plan"], motors) if motors else None


def sample_atol(name: str) -> float:
    """How far *name* may move with the sample still counting as unmoved."""
    for stage in SAMPLE_STAGES:
        if name in stage:
            return stage[name]
    return TEMPERATURE_ATOL


def scan_groups(paths: list[Path]) -> list[list[Path]]:
    """Gather the selection into the sets that are one measurement repeated.

    Three things have to agree. The plan and the motors, from
    :func:`scan_signature`. Where the sample sat and how cold it was, read
    from the baseline through :data:`SAMPLE_STAGES` -- which is what tells
    repeats from a series, 174, 175 and 176 being the same energy scan with
    the sample stepped 0.2 mm between them. And the points themselves, which
    only the curves can answer and :func:`same_grid` settles later.

    Each quantity is clustered across the selection with
    :func:`collect_uniques` at its own tolerance and every scan then assigned
    to the nearest surviving value, so the grouping does not depend on which
    scan happened to be looked at first. Groups come back in the order their
    first scan appears, which is the table's own order.

    A scan whose baseline carries no whole stage is left on its own rather
    than pooled on the temperature alone. That is the conservative way round:
    an endstation nobody has listed here draws its scans separately, which is
    what the viewer does today, instead of averaging a series into one curve
    that looks clean and is not one measurement.
    """
    readings: dict[str, dict[str, float]] = {}
    for path in paths:
        for stage in SAMPLE_STAGES:
            reading = baseline_values(path, tuple(stage))
            if len(reading) == len(stage):
                reading.update(baseline_values(path, (TEMPERATURE_FIELD,)))
                readings[str(path)] = reading
                break
    settings: dict[str, np.ndarray] = {}
    for name in {name for reading in readings.values() for name in reading}:
        measured = [
            reading[name] for reading in readings.values() if name in reading
        ]
        settings[name] = collect_uniques(np.array(measured), sample_atol(name))

    groups: dict[object, list[Path]] = {}
    for path in paths:
        signature = scan_signature(path)
        reading = readings.get(str(path))
        if signature is None or reading is None:
            # Nothing to repeat, or nothing saying where it was measured.
            key: object = ("alone", str(path))
        else:
            key = (
                signature,
                tuple(
                    (
                        name,
                        int(np.argmin(np.abs(settings[name] - reading[name]))),
                    )
                    for name in sorted(reading)
                ),
            )
        groups.setdefault(key, []).append(path)
    return list(groups.values())


def curve_step(x: np.ndarray) -> float:
    """The typical distance between one point of *x* and the next.

    The median rather than the span over the count: a list scan need not
    step evenly, and one long jump back to the start would make every real
    step look small beside it. Zero where nothing moved.
    """
    gaps = np.abs(np.diff(np.asarray(x, dtype=float)))
    gaps = gaps[np.isfinite(gaps) & (gaps > 0)]
    return float(np.median(gaps)) if gaps.size else 0.0


def same_grid(x_a: np.ndarray, x_b: np.ndarray) -> bool:
    """Whether two scans measured the same points, to within half a step.

    Compared row by row over the overlap, so a scan stopped part way is
    still the same grid as far as it got. The tolerance is a fraction of the
    scan's own step rather than an absolute distance, there being no one
    figure that suits an energy in keV and an angle in degrees at once; see
    :data:`SAME_X_STEPS` for where the half comes from.
    """
    rows = min(np.size(x_a), np.size(x_b))
    if rows < 2:
        return False
    first = np.asarray(x_a, dtype=float)[:rows]
    second = np.asarray(x_b, dtype=float)[:rows]
    shared = np.isfinite(first) & np.isfinite(second)
    step = max(curve_step(first), curve_step(second))
    if not shared.any() or not step:
        return False
    gap = float(np.max(np.abs(first[shared] - second[shared])))
    return gap <= SAME_X_STEPS * step


def _unit_tag(tags: dict[str, str], paths: list[Path]) -> str:
    """One name for the scans a reduction ran on: ``236``, ``235+236``, ``228-231``.

    A consecutive run collapses to its ends, which is what the macros write:
    ``current_pattern=[1, -1, -1, 1]`` is four scans in a row and naming all
    four would fill the legend. One scan keeps its bare number, so nothing
    the dichro reduction draws is renamed.
    """
    names = [tags[str(path)] for path in paths]
    if len(names) == 1:
        return names[0]
    numbers = [int(name) for name in names if name.isdigit()]
    if len(numbers) == len(names) and numbers == list(
        range(numbers[0], numbers[0] + len(numbers))
    ):
        return f"{names[0]}-{names[-1]}"
    return "+".join(names)


def contrast_units(paths: list[Path]) -> list[tuple[str, list[Path], str]]:
    """What each reduction runs on: ``(tag, paths, mode)``, in the table's order.

    A dichro scan is its own unit, as it always was. Lock-in scans are
    gathered by :func:`lockin_groups`, because the field they need stepped is
    commonly stepped between files rather than inside one.
    """
    tags = scan_tags(paths)
    units: list[tuple[str, list[Path], str]] = []
    grouped: dict[str, list[Path]] = {}
    for group in lockin_groups(
        [path for path in paths if reduce_mode(path) == MODE_LOCKIN_POL]
    ):
        for path in group:
            grouped[str(path)] = group
    for path in paths:
        group = grouped.get(str(path))
        if group is None:
            units.append((tags[str(path)], [path], MODE_DICHRO))
        elif str(group[0]) == str(path):
            # Emitted once, where its first member comes up, so the units
            # follow the order the scans are listed in.
            units.append((_unit_tag(tags, group), group, MODE_LOCKIN_POL))
    return units


def averaged_units(
    units: list[tuple[str, list[Path], str]], paths: list[Path]
) -> tuple[list[tuple[str, list[Path], str]], list[str]]:
    """Join the units that are one measurement repeated into one unit each.

    Pooled and reduced once rather than reduced apart and meaned after: the
    reducers already take the mean of every point that fell in one energy
    and field bin (``np.nanmean`` of each group, ``xmucd_reduce.py:262``),
    so handing them the repeats' columns together *is* the reference's own
    averaging. The two agree point for point where the repeats are even, and
    where they are not -- a scan stopped part way -- pooling is the better
    answer, weighing the points it actually has instead of letting a short
    curve count for as much as a whole one. It also needs no new arithmetic:
    235 and 236 are already joined this way.

    Two units are the same measurement when :func:`scan_groups` puts every
    one of their scans in one group and they reduce in the same mode. The
    grid test that curves-mode averaging applies is deliberately not made
    here: the reduction bins by energy *value*, so repeats that visited
    their points in a different order still pool correctly.

    Returns the units to reduce and the names of those that were joined.
    """
    where = {
        str(path): index
        for index, group in enumerate(scan_groups(paths))
        for path in group
    }
    tags = scan_tags(paths)
    pooled: dict[object, list[tuple[str, list[Path], str]]] = {}
    for unit in units:
        _tag, members, mode = unit
        # Keyed on the set of groups the unit's scans fall in, not on one of
        # them, so a unit the two groupings disagree about -- lockin_groups
        # joining what scan_groups separates -- pools only with a unit that
        # straddles the same pair rather than being forced into either.
        key = (
            mode,
            frozenset(where.get(str(path), str(path)) for path in members),
        )
        pooled.setdefault(key, []).append(unit)
    joined: list[tuple[str, list[Path], str]] = []
    names: list[str] = []
    for members in pooled.values():
        if len(members) == 1:
            joined.append(members[0])
            continue
        together = [path for _tag, unit, _mode in members for path in unit]
        tag = _unit_tag(tags, together)
        joined.append((tag, together, members[0][2]))
        names.append(tag)
    return joined, names


def read_contrast_columns(
    paths: list[Path],
    mode: str,
    detector: str,
    monitor: str,
    drop_spikes: bool = False,
) -> tuple[dict[str, np.ndarray], list[str], int]:
    """The columns one reduction runs on, joined over *paths*, and what is missing.

    Keyed by the reducers' own argument names rather than by field name, so a
    lock-in unit and a dichro scan come back in the same shape and the caller
    need only pick which reducer to hand them to.

    The lock-in ``step`` is the magnet: :data:`MAGNET_FIELD` where the scan
    stepped it within one file, and otherwise a constant column at the
    baseline reading. That substitution is the whole of why 235 and 236
    reduce together -- see :func:`lockin_groups`.

    Each file's columns are cut to its own shortest before being joined. A
    scan stopped part way leaves its last point a column short, and squaring
    that off only at the end would slide every later file against its own
    energies.

    With *drop_spikes* the rows where the beam was gone are left out, counted
    into the third return value. They go before the reducer rather than after
    it because the reducer averages every raw row sharing an energy: by the
    time it has run, one dead row of four has already become the whole bin's
    answer. Taken out here the bin simply averages the three that remain, so
    all 129 energies of scans 239-242 survive having a point removed.

    Dichro is the one case where an energy can go with the row. Its four
    quadrants are each averaged separately and :func:`reduce_dichro` keeps
    only the energies holding all four, so a row alone in its quadrant takes
    that energy with it -- scan 96 falls to 106 energies of 107 if the row at
    E = 7.24344 goes, and scan 98 to 43 of 44. Most quadrants hold two points
    and lose nothing; the lock-in scans above hold four. The energy is dropped
    rather than drawn from three quadrants because three cannot be reduced at
    all: every contrast is a difference between pairs of them.

    Each file is judged on its own column, which is what :func:`lone_dropouts`
    needs -- its own median, and its own rows either side.

    ``OSError`` propagates, as it does out of :func:`read_stream`.
    """
    lockin = mode == MODE_LOCKIN_POL
    fields = (
        [ENERGY_FIELD, *LOCKIN_FIELDS] if lockin else list(REDUCE_FIELDS)
    ) + [detector, monitor]
    parts: dict[str, list[np.ndarray]] = {}
    missing: list[str] = []
    removed = 0
    for path in paths:
        # Always primary, whatever the Stream box says. These are raw
        # per-point columns and a reduced stream holds a quarter of them with
        # the quadrants already averaged away -- yet scan 98 opens on
        # dichro_monitor, so reading the stream on screen would fail on one
        # of the two scans this was written for.
        data = read_stream(path, [*fields, MAGNET_FIELD], PRIMARY_STREAM)
        absent = [name for name in fields if name not in data]
        if lockin and MAGNET_FIELD not in data:
            setting = baseline_values(path, (MAGNET_FIELD,)).get(MAGNET_FIELD)
            if setting is None:
                absent.append(MAGNET_FIELD)
        if absent:
            missing += [name for name in absent if name not in missing]
            continue
        columns = {
            "energy": data[ENERGY_FIELD],
            "detector": data[detector],
            "monitor": data[monitor],
        }
        if lockin:
            columns["ac"] = data[AC_FIELD] - data[ACOFF_FIELD]
            columns["step"] = (
                data[MAGNET_FIELD]
                if MAGNET_FIELD in data
                else np.full(data[ENERGY_FIELD].size, setting)
            )
        else:
            columns["pr"] = data[PR_FIELD]
            columns["step"] = data[MAGNET_FIELD]
        size = min(column.size for column in columns.values())
        keep: slice | np.ndarray = slice(None)
        if drop_spikes:
            gone = lone_dropouts(columns["monitor"][:size])
            removed += int(gone.sum())
            keep = ~gone
        for name, column in columns.items():
            parts.setdefault(name, []).append(column[:size][keep])
    if missing:
        return {}, missing, 0
    return (
        {name: np.concatenate(columns) for name, columns in parts.items()},
        [],
        removed,
    )


def _half_max_crossing(
    x: np.ndarray, y: np.ndarray, half: float, start: int, stop: int, step: int
) -> float | None:
    """Linearly-interpolated X where *y* first crosses *half*, walking from
    *start* to (but not including) *stop* in steps of *step*."""
    i = start
    while i != stop:
        j = i + step
        yi, yj = y[i], y[j]
        if yi == half:
            return x[i]
        if (yi > half) != (yj > half):
            frac = (half - yi) / (yj - yi)
            return x[i] + frac * (x[j] - x[i])
        i = j
    return None


def curve_stats(x: np.ndarray, y: np.ndarray) -> dict[str, float | None]:
    """Max/min (with position), FWHM center/width and area, for one curve.

    ``center``/``fwhm`` come from linearly-interpolated crossings of the
    half-max level (halfway between the curve's min and max) either side of
    the peak -- the same convention tools like SPEC's peak search use --
    rather than a fit, so they are ``None`` when the peak sits at an end of
    the data (no crossing on that side).

    These are the table's "Center"/"FWHM" columns. ``peak_fit`` below
    answers the same two questions by fitting a model instead, and fills
    the separate "Fit center"/"Fit FWHM" columns. The two disagree on real
    data and are both right: the crossings sit on a curve that still has
    its background under it, while the fit takes a sloping background out
    first, so on scan 32 the crossings give 0.0421 where the fit gives
    0.0447.

    ``area``/``net_area`` are trapezoidal integrals over the whole curve,
    the second one above the straight line joining its first and last
    points. Both are signed, so a contrast that runs negative integrates
    negative.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]

    stats: dict[str, float | None] = {
        "max_y": None,
        "max_x": None,
        "min_y": None,
        "min_x": None,
        "center": None,
        "fwhm": None,
        "area": None,
        "net_area": None,
    }
    if x.size == 0:
        return stats

    imax = int(np.argmax(y))
    imin = int(np.argmin(y))
    stats["max_y"], stats["max_x"] = y[imax], x[imax]
    stats["min_y"], stats["min_x"] = y[imin], x[imin]

    half = (y[imax] + y[imin]) / 2.0
    left = _half_max_crossing(x, y, half, imax, 0, -1)
    right = _half_max_crossing(x, y, half, imax, len(x) - 1, 1)
    if left is not None and right is not None:
        stats["center"] = (left + right) / 2.0
        stats["fwhm"] = abs(right - left)

    # Two points to make a trapezoid from. One gets no area rather than a
    # zero: an integral over a single reading is zero in the arithmetic but
    # says nothing, and the table shows None as a blank cell, which is the
    # honest answer. Sorted ascending first, and that is correctness rather
    # than tidiness -- a scan started at the top of an edge steps downwards
    # (171 and 173 both run 7.30399 to 7.19292 keV) and integrating one as
    # read silently returns the negative: 171's LockAC gives -290364 that
    # way against +290364 sorted. Sorted into locals, never in place: the
    # caller's x is shared between a scan's curves and is what the export
    # groups its columns by.
    if x.size >= 2:
        order = np.argsort(x, kind="stable")
        xs, ys = x[order], y[order]
        stats["area"] = float(_trapezoid(ys, xs))
        # Both ends at one X leaves no line to draw -- and no width to
        # integrate over either, so the raw area is already 0. Blank, not a
        # copy of it, so the column never claims a background was removed.
        if xs[-1] != xs[0]:
            slope = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
            stats["net_area"] = float(
                _trapezoid(ys - (ys[0] + (xs - xs[0]) * slope), xs)
            )
    return stats


# The peak shapes "Fit" can use, each normalised to unit height and -- this
# is the part worth getting right -- parameterised by its FWHM directly
# rather than by a sigma or a half-width. A Gaussian's sigma is the FWHM
# over 2.355 and a Lorentzian's gamma is half of it, so carrying the native
# parameter would mean three conversions going in to seed the fit and three
# coming back out to report it, including the squared Lorentzian's
# 2*gamma*sqrt(sqrt(2)-1), which is the kind of constant that gets typed
# wrong once and never noticed. This way the seed and the answer are both
# just the width. Another shape is one more line here.
def _peak_gauss(x: np.ndarray, centre: float, fwhm: float) -> np.ndarray:
    return np.exp(-4.0 * np.log(2.0) * ((x - centre) / fwhm) ** 2)


def _peak_lorentzian(x: np.ndarray, centre: float, fwhm: float) -> np.ndarray:
    return 1.0 / (1.0 + 4.0 * ((x - centre) / fwhm) ** 2)


def _peak_lorentzian2(x: np.ndarray, centre: float, fwhm: float) -> np.ndarray:
    root = 4.0 * (np.sqrt(2.0) - 1.0) * ((x - centre) / fwhm) ** 2
    return 1.0 / (1.0 + root) ** 2


# Lorentzian and Gaussian at one shared FWHM, which is the whole trick of
# the pseudo-Voigt: both are already unit height at the centre and already
# half height at centre +/- fwhm/2, so any mix of them is too, and the
# blend keeps that FWHM without a correction term. *eta* is the Lorentzian
# fraction -- 1 is the pure Lorentzian, 0 the pure Gauss -- and it is
# fitted rather than assumed, which is what PEAK_MIX below arranges.
def _peak_pseudo_voigt_mixed(
    x: np.ndarray, centre: float, fwhm: float, eta: float
) -> np.ndarray:
    return eta * _peak_lorentzian(x, centre, fwhm) + (1.0 - eta) * _peak_gauss(
        x, centre, fwhm
    )


# The same shape with the mixing held at a half: where the fit starts, and
# what it falls back to when a free eta will not stand. Kept as a function
# of three parameters of its own because that is the signature PEAK_SHAPES
# promises and that every shape without a mixing has.
def _peak_pseudo_voigt(x: np.ndarray, centre: float, fwhm: float) -> np.ndarray:
    return _peak_pseudo_voigt_mixed(x, centre, fwhm, 0.5)


def _pseudo_voigt_area(eta: float) -> float:
    """The mix's own integral, at unit height and unit FWHM.

    Linear in *eta*, because the integral of a sum is the sum of the
    integrals: no third closed form is needed beyond the two above, and a
    fitted mixing reaches the "Fit area" column without a quadrature
    anywhere. At eta a half this is the constant the fixed shape carried
    before the mixing was fitted, to the last bit -- halving is exact.
    """
    return float(
        eta * np.pi / 2.0 + (1.0 - eta) * np.sqrt(np.pi / (4.0 * np.log(2.0)))
    )


# Each shape is paired with the integral of its own function over the whole
# line, at unit height and unit FWHM, so that a new shape cannot be added
# with its area quietly left out. These are closed forms, checked against a
# trapezoid of the function itself: the Gaussian and the square agree to
# nine figures, the Lorentzian only to four, the rest being tail that a
# numeric check truncated at 4000 half-widths has still not collected. That
# same heavy tail is what the "Fit area" column carries -- see peak_fit.
PEAK_SHAPES = {
    "Gauss": (_peak_gauss, float(np.sqrt(np.pi / (4.0 * np.log(2.0))))),
    "Lorentzian": (_peak_lorentzian, float(np.pi / 2.0)),
    "Lorentzian squared": (
        _peak_lorentzian2,
        float(np.pi / (4.0 * np.sqrt(np.sqrt(2.0) - 1.0))),
    ),
    # Half of each of the two above it. Its area comes from the same
    # function the fitted mixing uses, so the two cannot drift apart.
    "Pseudo-Voigt": (_peak_pseudo_voigt, _pseudo_voigt_area(0.5)),
}
# The shapes that carry a mixing parameter of their own, paired with the
# area factor as a function of it. Being listed here is what tells both
# fitters to spend a fourth parameter a component on that shape; a shape
# absent from it is solved on three, exactly as all four were before, and
# only a shape listed here is given the "Fit eta" column.
#
# Both fitters also solve the fixed shape and keep whichever scores
# better, so the extra parameter cannot cost anything, and a free eta is
# refused outright unless it lands clear of its two ends: a mixing resting
# on 0 or 1 is the solver saying "pure Gauss" or "pure Lorentzian" rather
# than measuring a blend, and both of those are in the combo already.
#
# Measured on the scans the tooltips name. It is kept where it buys
# something -- scan 31's derivative at two peaks goes 0.996363 to
# 0.997166 at eta 0.690 and 0.802, scan 34 goes 0.998357 to 0.999409 at
# 0.851 and 0.073 -- and refused where it does not: scan 32 at two peaks
# runs its first component to the 1.000 rail, and scan 109 drives one
# component under the step size. Those two fall back to the fixed half,
# report precisely the fit they reported before and leave the eta cell
# blank, which is the honest reading and is why the fallback exists at
# all rather than letting one of the three cases "Multi" is for decline.
PEAK_MIX = {
    "Pseudo-Voigt": (_peak_pseudo_voigt_mixed, _pseudo_voigt_area),
}
# How far clear of either of its two ends a fitted mixing has to land
# before it is read as a measurement rather than as the solver resting
# against a wall. Small on purpose, because the wall is the only thing
# being detected: a bounded least squares that wants to be outside its
# box stops hard against it, and scan 32 at two peaks stops at one ulp
# below 1 -- a plain "strictly less than 1" lets that straight through
# and prints 1.000 as though it had been measured. A mixing genuinely
# near a pure shape is a different statement and is nowhere near this
# close; it is reported as whatever number it is.
PEAK_MIX_RAIL = 1e-6
# Below this many finite points a five-parameter fit is not measuring
# anything: over 200 six-point curves of pure noise, five parameters scored
# a median R2 of 0.71, one in twenty above 0.99 and the best a flat 1.0000.
# Seven leaves two degrees of freedom, which is thin but honest, and the
# sanity gate below is what actually does the rejecting.
PEAK_FIT_MIN_POINTS = 7
# How far above the data's own range a fitted peak may stand before it is
# read as a solver artefact rather than a signal. Generous on purpose --
# a real peak sits inside the data range, so ten times it is already far
# past anything physical -- but finite, which is the point: see peak_fit.
PEAK_FIT_MAX_HEIGHT = 10.0
# How prominent a bump must be, as a fraction of the data's own Y range,
# before find_peaks offers it as a starting point for one of "Multi"'s
# components. A seed and nothing more: a bump it passes over can still be
# fitted -- scan 34's second peak is not a local maximum at all, and is
# found by the padding in _peak_seeds -- and one it offers can still be
# turned down by the gate afterwards.
PEAK_SEED_PROMINENCE = 0.05
# Bounded least squares walks to its answer in much smaller steps than the
# unbounded fit peak_fit does, and with three parameters a component there
# are a good many more of them to take. Generous rather than measured: the
# failure mode is a RuntimeError and a blank row, and even eight parameters
# over fifty points is a few milliseconds.
PEAK_FIT_MULTI_MAXFEV = 60000
# What the count beside "Multi" offers. It starts at two because one peak
# is what the button being off already means, so a "1" in the list would
# be a second way of saying the same thing. Two and three are measured on
# the scans the tooltips name; four is allowed by the same arithmetic and
# held in check by the minimum-points rule, which asks three more readings
# for every component.
PEAK_COUNTS = ["2", "3", "4"]


def _peak_mix_refine(
    xs: np.ndarray,
    ys: np.ndarray,
    mixed,
    start: tuple[float, float, float, float, float],
    span: float,
    spread: float,
    step: float,
) -> tuple[float, float, float, float, float, float] | None:
    """Re-solve one fitted peak with the shape's mixing set free.

    *mixed* is the shape's four-parameter form and *start* is what
    :func:`peak_fit` already has -- height, centre, FWHM, slope, offset.
    The mixing is seeded at the half the fixed shape assumes, so the
    six-parameter solve starts from the five-parameter answer with only
    the mixing left to find.

    Returns the six refined parameters, or ``None`` for "keep the answer
    you had". That happens when the solver gives up, when the refined
    peak fails any of the four geometric tests ``peak_fit`` has just
    applied -- a freed mixing can trade width against shape and walk out
    of them -- or when the mixing does not land strictly between its own
    two ends, which PEAK_MIX explains.

    Unbounded, like the fit it refines. Bounding it would switch
    ``curve_fit`` from Levenberg-Marquardt to a trust region and move
    every other parameter along with the mixing; the ends are tested
    here afterwards instead, exactly as the geometry already is.
    """
    height, centre, fwhm, line_slope, line_offset = start

    def model(grid, height, centre, fwhm, line_slope, line_offset, eta):
        return (
            height * mixed(grid, centre, fwhm, eta)
            + line_slope * grid
            + line_offset
        )

    try:
        # As peak_fit: catch_warnings only saves and restores the filter
        # state, and the simplefilter is what silences the covariance and
        # overflow warnings scipy prints straight to stderr.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, _covariance = curve_fit(
                model,
                xs,
                ys,
                p0=[height, centre, fwhm, line_slope, line_offset, 0.5],
                maxfev=20000,
            )
    except (RuntimeError, ValueError):
        return None
    if not np.all(np.isfinite(popt)):
        return None
    height, centre, fwhm, line_slope, line_offset, eta = (
        float(value) for value in popt
    )
    if not xs[0] <= centre <= xs[-1]:
        return None
    if not 0.0 < fwhm <= span:
        return None
    if abs(height) > PEAK_FIT_MAX_HEIGHT * spread:
        return None
    if fwhm < step:
        return None
    # Nothing holds this solve inside the two ends, so a mixing that
    # wants to be outside them lands outside and is refused here. The
    # same margin as the bounded path all the same, so that the two
    # report the same shapes as measured and the same ones as pure.
    if not PEAK_MIX_RAIL <= eta <= 1.0 - PEAK_MIX_RAIL:
        return None
    return height, centre, fwhm, line_slope, line_offset, eta


def peak_fit(
    x: np.ndarray,
    y: np.ndarray,
    shape: str,
    fwhm_seed: float | None = None,
) -> dict[str, float | None]:
    """Fit *shape* plus a sloping line to one curve.

    The model is ``height * peak(x; center, fwhm) + slope * x + offset``,
    five parameters. The background is always sloping and never optional:
    on the rocking curves in TiSe2 (32 to 38) it runs 1754 counts down to
    476 across scan 32 alone, and fitting it lifts the mean R2 of the whole
    curve from 0.957 to 0.980 for a Gaussian -- 0.976 for a flat offset, so
    it is the slope that is bought and not merely a pedestal. (Those three
    are measured against the raw curve, which is not how ``fit_r2`` below
    is reckoned; they compare background models, not fits.) ``height`` is
    measured above that line, which is why it comes out below the curve's
    own maximum (46664 against 51863 on scan 32).

    Every value is ``None`` when no trustworthy fit came back, which the
    table shows as a blank row and the export writes as ``(none)``. That
    happens a great deal more often than it might seem it should, and
    deliberately: see the gate below.

    *fwhm_seed* is ``curve_stats``' crossing width when it found one. It is
    passed in rather than recomputed because the caller has just worked it
    out for the same curve.
    """
    stats: dict[str, float | None] = {
        "fit_center": None,
        "fit_fwhm": None,
        "fit_height": None,
        "fit_area": None,
        "fit_eta": None,
        "fit_r2": None,
        "fit_slope": None,
        "fit_offset": None,
    }
    entry = PEAK_SHAPES.get(shape)
    if entry is None:
        return stats
    peak, area_factor = entry

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # Both, as curve_stats masks: an averaged curve carries a NaN pad in X
    # and a combined one carries NaN outside the overlap.
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if x.size < PEAK_FIT_MIN_POINTS:
        return stats
    # Sorted into locals, never in place, for the reason curve_stats gives:
    # the caller's x is shared between a scan's curves and is what the
    # export groups its columns by. Sorting matters here too -- 171 and 173
    # step downwards, and the end points below are read as first and last.
    order = np.argsort(x, kind="stable")
    xs, ys = x[order], y[order]
    span = float(xs[-1] - xs[0])
    # Every reading at one motor position: no span to fit along, and the
    # model grid the caller draws would be 400 copies of one point.
    if span <= 0.0:
        return stats

    # Seeded off the straight line through the end points. The largest
    # residual from it gives both the height and the position, and because
    # that residual keeps its sign a dip needs no special case at all --
    # it simply seeds a negative height and fits one.
    slope = (ys[-1] - ys[0]) / span
    offset = ys[0] - slope * xs[0]
    residual = ys - (slope * xs + offset)
    tallest = int(np.argmax(np.abs(residual)))
    width = fwhm_seed if fwhm_seed and fwhm_seed > 0.0 else span / 10.0

    def model(grid, height, centre, fwhm, line_slope, line_offset):
        return (
            height * peak(grid, centre, fwhm) + line_slope * grid + line_offset
        )

    try:
        # catch_warnings on its own only saves and restores the filter
        # state -- it suppresses nothing, and a flat curve prints
        # "OptimizeWarning: Covariance of the parameters could not be
        # estimated" straight to stderr through it. The simplefilter is
        # what silences that, and the RuntimeWarning from an overflowing
        # np.exp on a bad step. This file neither prints nor logs.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, _covariance = curve_fit(
                model,
                xs,
                ys,
                p0=[residual[tallest], xs[tallest], width, slope, offset],
                maxfev=20000,
            )
    except (RuntimeError, ValueError):
        return stats
    height, centre, fwhm, line_slope, line_offset = (float(v) for v in popt)

    # The gate. An unbounded five-parameter fit returns nonsense without
    # raising and without complaining: over 400 curves of pure noise, 76
    # came back with a negative width and 33 with a centre outside the data
    # entirely, and the table would have printed those to six figures. The
    # four tests are geometric and cost nothing.
    #
    # The last two are not tidiness. A tall narrow spike between two
    # readings passes the first two happily, and the plot cannot survive
    # it: matplotlib's relim() takes the model line as readily as the data,
    # so the Y axis scales to the spike and the real curve flattens onto
    # the axis. Before these two the worst model spanned 6.6e28 times the
    # data's own Y range; after them, under ten.
    if not np.all(np.isfinite(popt)):
        return stats
    if not xs[0] <= centre <= xs[-1]:
        return stats
    if not 0.0 < fwhm <= span:
        return stats
    spread = float(ys.max() - ys.min())
    if spread <= 0.0 or abs(height) > PEAK_FIT_MAX_HEIGHT * spread:
        return stats
    # A peak narrower than the step size was never resolved -- the model is
    # interpolating between two readings, not measuring a width.
    step = float(np.median(np.diff(xs)))
    if fwhm < step:
        return stats

    # The mixing, for a shape that has one. Refined from the five
    # parameters above rather than solved alongside them from scratch, so
    # that a shape without a mixing goes through exactly the code it
    # always did, and so that a pseudo-Voigt whose mixing will not stand
    # still reports the fit it would have reported anyway -- the answer
    # already in hand is the fallback, at no extra solve.
    eta: float | None = None
    mixed = PEAK_MIX.get(shape)
    # One more parameter wants one more reading, keeping the two degrees
    # of freedom PEAK_FIT_MIN_POINTS was chosen to leave.
    if mixed is not None and xs.size >= PEAK_FIT_MIN_POINTS + 1:
        refined = _peak_mix_refine(
            xs,
            ys,
            mixed[0],
            (height, centre, fwhm, line_slope, line_offset),
            span,
            spread,
            step,
        )
        if refined is not None:
            height, centre, fwhm, line_slope, line_offset, eta = refined
            area_factor = mixed[1](eta)
    profile = (
        peak(xs, centre, fwhm)
        if eta is None
        else mixed[0](xs, centre, fwhm, eta)
    )

    # R2 against the background-subtracted data, not the raw curve. With a
    # sloping line in the model the line otherwise explains most of the
    # variance by itself and the number stops discriminating: pure noise
    # scored 0.62 to 0.87 that way, which a reader would take for a real
    # fit. Taken over the peak alone, the same noise sits around 0.29 and
    # the rocking curves stay above 0.96.
    base = line_slope * xs + line_offset
    target = ys - base
    ss_total = float(np.sum((target - target.mean()) ** 2))
    if ss_total > 0.0:
        ss_residual = float(np.sum((target - height * profile) ** 2))
        stats["fit_r2"] = 1.0 - ss_residual / ss_total
    stats["fit_center"] = centre
    stats["fit_fwhm"] = fwhm
    stats["fit_height"] = height
    # Blank for every shape without a mixing, and for one whose mixing was
    # refused: in both the model drawn is the three-parameter shape, and a
    # printed 0.5 would read as a measurement when it is an assumption.
    stats["fit_eta"] = eta
    # The integrated intensity of the peak alone, the sloping line taken
    # off: height x FWHM x the shape's own integral, which is exact rather
    # than a sum over the drawn points. Deliberately the whole line and not
    # just the scanned range, so that the same peak measured over a wider
    # sweep reports the same intensity and the repeats stay comparable --
    # which is the whole reason for fitting them. The price is paid by the
    # Lorentzian, 6 to 11% of whose area on scans 32 to 38 lies in tails
    # beyond where anybody actually measured; the Gaussian extrapolates
    # nothing worth printing and its square under 1%. A dip fits a negative
    # height and so reports a negative area, exactly as the height column
    # does, rather than pretending to an intensity it does not have.
    stats["fit_area"] = height * fwhm * area_factor
    stats["fit_slope"] = line_slope
    stats["fit_offset"] = line_offset
    return stats


def fit_columns(mixed: bool) -> list[str]:
    """The fit headings, the mixing's among them only where it is fitted.

    One function so that the headings and the cells under them cannot be
    built from two different ideas of the order -- ``_set_stats_columns``
    and ``_update_stats_table`` both come here.
    """
    if not mixed:
        return list(STATS_FIT_COLUMNS)
    at = STATS_FIT_COLUMNS.index("Fit FWHM") + 1
    return STATS_FIT_COLUMNS[:at] + [STATS_MIX_COLUMN] + STATS_FIT_COLUMNS[at:]


def _blank_stats() -> dict[str, float | None]:
    """A statistics row with every curve column empty.

    For the per-peak rows "Multi" adds: a component of a fit has no curve
    of its own, so it has no maximum, no crossing width and no trapezoid.
    Built by asking ``curve_stats`` about no data rather than by listing
    its keys here, because ``_update_stats_table`` reads those columns by
    subscript -- a row missing one raises rather than showing a blank --
    and a list of them kept by hand is a list that drifts.
    """
    return curve_stats(np.empty(0), np.empty(0))


def _peak_seeds(
    xs: np.ndarray,
    ys: np.ndarray,
    peaks: int,
    step: float,
    peak,
    slope: float,
    offset: float,
    rises_first: bool,
) -> list[tuple[float, float, float]]:
    """Starting height, centre and width for each of *peaks* components.

    ``find_peaks`` is run over the data and over its negative, so a dip
    seeds as readily as a rise. That is the whole of the derivative case:
    the derivative of a peak is a lobe up then a lobe down, and seeding
    only maxima finds one of the two -- which is exactly the one-lobe trap
    the "Fit" tooltip has warned about on scans 32 to 38.

    Candidates are taken in order of prominence, and *rises_first* says
    whether every maximum is ranked ahead of every dip or the two signs
    are pooled. Neither order is right on its own, which is why the
    caller tries both and keeps the fit that scores better; the sort
    below says what each of them gets wrong.

    The widths come from ``find_peaks``' own ``widths`` rather than a
    fraction of the scan, and that is not a refinement. A tenth of the
    range on scan 31's derivative is about twice the true lobe width and
    wider than the gap between the lobes; seeded that way the solver
    walked off and returned no fit at all.
    """
    spread = float(np.max(ys) - np.min(ys))
    reach = float(xs[-1] - xs[0])
    found: list[tuple[int, float, float, float, float]] = []
    for rank, sign in ((0, 1.0), (1, -1.0)):
        # width=1 is the loosest setting that makes find_peaks report a
        # width at all; prominence is what keeps it off the noise.
        index, props = find_peaks(
            sign * ys, prominence=PEAK_SEED_PROMINENCE * spread, width=1
        )
        for i, rise, samples in zip(
            index, props["prominences"], props["widths"]
        ):
            found.append(
                (
                    rank,
                    float(rise),
                    sign * float(rise),
                    float(xs[i]),
                    max(float(samples) * step, 2.0 * step),
                )
            )
    # Each ordering has one shape it gets wrong, and they are not the
    # same shape, which is why both are tried.
    #
    # Pooling the signs on prominence alone sounds even-handed and is
    # not. The valley between two separate peaks is about as prominent as
    # the smaller of the two, so a single pool decides between them on a
    # margin of about a percent: a clean pair of Gaussians at 0.3 and
    # 0.65 loses on exactly that -- valley 59.05 against peak 58.30 --
    # and the fit then hangs a wide negative component off the left edge
    # and scores R2 0.955 with every gate below satisfied.
    #
    # Ranking every rise first settles that the way round the data means
    # it, and starves the derivative case instead, wherever the positive
    # side has a candidate to spare. Scan 35's derivative offers two
    # rises and two dips, and its second lobe is the most prominent thing
    # in the scan: taken rises first it is passed over for a 1.7e5 bump,
    # both components come back positive and R2 falls to 0.415, below
    # even the one-lobe fit this is all meant to replace.
    found.sort(
        key=lambda candidate: (
            candidate[0] if rises_first else 0,
            -candidate[1],
        )
    )
    seeds = [(h, c, w) for _rank, _rise, h, c, w in found[:peaks]]

    # Fewer candidates than asked for is the interesting case rather than
    # an error, and it is why the count has to be the user's to set: scan
    # 34's second peak hides under the falling flank of the first, so it
    # is not a local maximum and no peak finder will ever report it.
    #
    # The shortfall is seeded where the data departs furthest from what
    # has been seeded so far -- peak_fit's own trick of taking the largest
    # residual, applied once per missing component instead of once. Where
    # that residual keeps its sign a dip seeds a negative height and fits
    # one, exactly as it does there. Spreading the spare seeds evenly
    # across the scan instead was tried and is worse than it looks: on
    # scan 34 it lands an even third of the way along, to the *left* of
    # the one real candidate, and the fit settles for a 645-count
    # satellite out at 44.80 (R2 0.9818) rather than the shoulder at
    # 44.74 that is actually there (0.9926).
    width = (
        float(np.median([w for _k, _r, _h, _c, w in found])) if found else 0.0
    )
    if width <= 0.0:
        width = max(reach / 10.0, 2.0 * step)
    while len(seeds) < peaks:
        settled = slope * xs + offset
        for height, centre, fwhm in seeds:
            settled = settled + height * peak(xs, centre, fwhm)
        residual = ys - settled
        apart = int(np.argmax(np.abs(residual)))
        # A residual of nothing left to explain seeds no position worth
        # having, and two components on one point will not converge.
        if not np.isfinite(residual[apart]) or residual[apart] == 0.0:
            seeds.append((spread, float(xs[xs.size // 2]), width))
        else:
            seeds.append((float(residual[apart]), float(xs[apart]), width))
    return seeds


def _peak_solve(
    xs: np.ndarray,
    ys: np.ndarray,
    peaks: int,
    step: float,
    peak,
    seeds: list[tuple[float, float, float]],
    slope: float,
    offset: float,
    mixed=None,
) -> (
    tuple[
        float | None,
        list[tuple[float, float, float, float | None]],
        float,
        float,
    ]
    | None
):
    """Solve one seeding and gate the answer, or ``None`` if it will not do.

    Returns ``(r2, [(centre, fwhm, height, eta), ...], slope, offset)``
    with the components ordered by centre. ``None`` is this seeding
    declined: the solver gave up, or came back with something the gates
    refuse, and the caller is free to try another seeding or to report no
    fit at all.

    *mixed* is the shape's four-parameter form for a shape that has one,
    and giving it is what spends a fourth parameter a component on the
    mixing; left out, every component is solved on three and its ``eta``
    comes back ``None``. The caller solves both ways and keeps the better
    answer -- see PEAK_MIX and ``peak_fit_multi``.

    Expects *xs* sorted, *step* and the data's span and range positive,
    and warnings already silenced -- ``peak_fit_multi`` sees to all four.
    """
    span = float(xs[-1] - xs[0])
    spread = float(ys.max() - ys.min())
    # The stride through the solver's flat parameter vector, which is
    # what the mixing costs: everything below indexes components by it.
    width = 3 if mixed is None else 4
    guess: list[float] = []
    lower: list[float] = []
    upper: list[float] = []
    for height, centre, fwhm in seeds:
        guess += [height, centre, fwhm]
        # Bounded, where peak_fit only checks afterwards. With several
        # components free at once there is always somewhere for a runaway
        # to hide, so the solver is held to the same geometry the gate
        # below tests: inside the data, and no taller than
        # PEAK_FIT_MAX_HEIGHT times the data's range.
        #
        # The width is the one that deliberately falls short of its gate,
        # at a tenth of a step rather than a step. A bound is a place the
        # solver is allowed to sit, so a bound at the gate is a gate that
        # can never fire -- and it did not: scan 171's derivative, which
        # is too ragged for its own step size and which peak_fit declines
        # outright, came back as two components pinned at exactly one
        # step and R2 0.17. Left room beneath it, the solver goes under,
        # the strict test below catches it, and the two paths agree
        # again. The bound is kept off zero only because a width of
        # nothing divides by nothing in every shape in PEAK_SHAPES.
        lower += [-PEAK_FIT_MAX_HEIGHT * spread, xs[0], step / 10.0]
        upper += [PEAK_FIT_MAX_HEIGHT * spread, xs[-1], span]
        if mixed is not None:
            # Seeded at the half the fixed shape assumes, so the mixed
            # solve starts from the shape the fixed one would have
            # fitted, and bounded by the two pure shapes it blends --
            # outside them one of the halves turns negative and the
            # profile stops being a peak. Unlike the width above, these
            # two are deliberately set *at* the gate rather than beyond
            # it, and the margin that lets the gate still fire is kept
            # on the gate's side instead -- see PEAK_MIX_RAIL. A mixing
            # resting on either wall is the solver asking for a pure
            # Gauss or a pure Lorentzian, which is what gets caught.
            guess.append(0.5)
            lower.append(0.0)
            upper.append(1.0)
    guess += [slope, offset]
    lower += [-np.inf, -np.inf]
    upper += [np.inf, np.inf]
    # curve_fit raises ValueError on a seed outside its bounds, and a
    # padded seed on a short scan can be: reach/10 is under one step
    # wherever there are fewer than ten readings.
    guess = [min(max(g, lo), hi) for g, lo, hi in zip(guess, lower, upper)]

    def model(grid, *free):
        total = free[-2] * grid + free[-1]
        for k in range(peaks):
            block = free[width * k : width * k + width]
            total = total + block[0] * (
                peak(grid, block[1], block[2])
                if mixed is None
                else mixed(grid, block[1], block[2], block[3])
            )
        return total

    try:
        popt, _covariance = curve_fit(
            model,
            xs,
            ys,
            p0=guess,
            bounds=(lower, upper),
            maxfev=PEAK_FIT_MULTI_MAXFEV,
        )
    except (RuntimeError, ValueError, np.linalg.LinAlgError):
        return None
    if not np.all(np.isfinite(popt)):
        return None
    line_slope, line_offset = float(popt[-2]), float(popt[-1])

    found: list[tuple[float, float, float, float | None]] = []
    for k in range(peaks):
        height = float(popt[width * k])
        centre = float(popt[width * k + 1])
        fwhm = float(popt[width * k + 2])
        eta = None if mixed is None else float(popt[width * k + 3])
        # peak_fit's four geometric tests, once per component, and the
        # whole seeding declined when any one of them fails. Declined
        # whole rather than per component because the model is drawn as a
        # sum: one runaway takes the Y axis with it however sane its
        # neighbours are, which is the matplotlib relim() failure the
        # height test was put there to stop. A solver sitting exactly on
        # a bound trips the same tests, since the bounds are that
        # geometry.
        if not xs[0] <= centre <= xs[-1]:
            return None
        if not 0.0 < fwhm <= span:
            return None
        if abs(height) > PEAK_FIT_MAX_HEIGHT * spread:
            return None
        if fwhm < step:
            return None
        # And the mixing clear of its two ends, where the four above are
        # the geometry. A component resting on 0 or 1 declines the whole
        # mixed seeding, which costs nothing: peak_fit_multi solves the
        # fixed shape alongside this one and reports that instead. Scan
        # 32 at two peaks is the measured case.
        if eta is not None and not (
            PEAK_MIX_RAIL <= eta <= 1.0 - PEAK_MIX_RAIL
        ):
            return None
        found.append((centre, fwhm, height, eta))
    # By centre, so "peak 1" is the leftmost -- see peak_fit_multi. By
    # the centre alone and not the whole tuple, so that the mixing, which
    # is None for every shape without one, is never compared.
    found.sort(key=lambda component: component[0])

    # R2 against the background-subtracted data, for the reason peak_fit
    # gives: with a sloping line in the model the line explains most of
    # the variance by itself and the number stops discriminating. One R2
    # for the whole model, since the components were fitted together and
    # none of them has an R2 of its own.
    base = line_slope * xs + line_offset
    target = ys - base
    ss_total = float(np.sum((target - target.mean()) ** 2))
    r2: float | None = None
    if ss_total > 0.0:
        shaped = np.zeros_like(xs)
        for centre, fwhm, height, eta in found:
            shaped = shaped + height * (
                peak(xs, centre, fwhm)
                if eta is None
                else mixed(xs, centre, fwhm, eta)
            )
        r2 = 1.0 - float(np.sum((target - shaped) ** 2)) / ss_total
    return r2, found, line_slope, line_offset


def peak_fit_multi(
    x: np.ndarray,
    y: np.ndarray,
    shape: str,
    peaks: int,
) -> dict[str, float | None]:
    """Fit *peaks* copies of *shape* over one shared sloping line.

    The model is ``sum(height_k * peak(x; centre_k, fwhm_k))`` plus
    ``slope * x + offset``: three parameters a component and two for the
    line, which the components share because one background under several
    peaks is what the data has.

    Beside the keys ``peak_fit`` fills, this adds ``fit_peaks`` -- a list
    of one row per component, **ordered by centre**. That ordering is what
    makes a trend across scans mean anything: "peak 1" is then always the
    leftmost, and its line follows the same physical peak through a set of
    repeats even where the two swap which is the taller.

    The four per-peak values on the returned dict itself stay ``None``, so
    the curve's own row shows five blank cells and the numbers appear on
    the peak rows underneath it. ``fit_slope``/``fit_offset`` are kept,
    because the shared line is what ``peak_model`` needs to draw the sum.

    Two things this cannot do, both measured and both in the tooltips.

    R2 cannot choose *peaks* for you. The models are nested, so it only
    ever rises: on scan 32, a single clean rocking curve, one peak scores
    0.97832 and two score 0.99761.

    And a component too many does not announce itself by failing. Asked
    for two, a single peak is *divided* rather than fitted alongside a
    rejected spare, and the division need not be lopsided: scan 32 comes
    back as a broad base (FWHM 0.0914) under a sharp core (0.0347) at
    much the same place, and neither of those is the 0.0447 that one peak
    reports and that the crossings agree with. The derivative of a knife
    edge divides too -- scan 101 grows a negative component at -6.4542
    and its main width moves 0.02095 to 0.02340, 12%, which is more than
    the drift the whole trend exists to show. Every component is therefore
    given a row of its own and the gate below stays geometric: a broad
    base is as sane a shape as a peak, and no test here can know it was
    not wanted.
    """
    stats: dict[str, float | None] = {
        "fit_center": None,
        "fit_fwhm": None,
        "fit_height": None,
        "fit_area": None,
        "fit_eta": None,
        "fit_r2": None,
        "fit_slope": None,
        "fit_offset": None,
        "fit_peaks": None,
    }
    entry = PEAK_SHAPES.get(shape)
    if entry is None or peaks < 1:
        return stats
    peak, area_factor = entry
    mixed = PEAK_MIX.get(shape)

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # Both, as peak_fit and curve_stats mask: an averaged curve carries a
    # NaN pad in X and a combined one carries NaN outside the overlap.
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    # Three parameters a component against peak_fit's seven points for
    # five, written so that one peak gives PEAK_FIT_MIN_POINTS back exactly.
    if x.size < PEAK_FIT_MIN_POINTS + 3 * (peaks - 1):
        return stats
    # Sorted into locals, never in place, for the reason curve_stats gives:
    # the caller's x is shared between a scan's curves and is what the
    # export groups its columns by.
    order = np.argsort(x, kind="stable")
    xs, ys = x[order], y[order]
    span = float(xs[-1] - xs[0])
    spread = float(ys.max() - ys.min())
    if span <= 0.0 or spread <= 0.0:
        return stats
    # The lower bound on every width below, so it has to be positive or
    # the bounds themselves are malformed. Zero where more than half the
    # readings repeat an X position.
    step = float(np.median(np.diff(xs)))
    if step <= 0.0:
        return stats

    try:
        # As peak_fit: catch_warnings only saves and restores the filter
        # state, and it is the simplefilter that silences the covariance
        # and overflow warnings scipy prints straight to stderr. polyfit
        # is inside it too, for its RankWarning on a flat end. This file
        # neither prints nor logs.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # The background is seeded from a line through both ends of
            # the scan rather than peak_fit's chord between the two end
            # points. A derivative sits at about zero at either end, so
            # that chord is mostly noise, and seeding the offset from it
            # put scan 31's fit in the wrong basin. A tenth of the points
            # each side, never fewer than two.
            edge = max(2, xs.size // 10)
            slope, offset = np.polyfit(
                np.concatenate((xs[:edge], xs[-edge:])),
                np.concatenate((ys[:edge], ys[-edge:])),
                1,
            )
            # Both the orderings _peak_seeds offers, keeping whichever
            # fits better, because each of them has a shape it gets
            # wrong and they are not the same shape -- the comment there
            # names the two scans. Comparing R2 between these two is
            # sound where comparing it between peak counts is not: this
            # is one model with one parameter count, solved twice from
            # different starting points, so the better number is the
            # better fit of the same thing and not a reward for more
            # freedom.
            seedings: list[list[tuple[float, float, float]]] = []
            for rises_first in (False, True):
                seeds = _peak_seeds(
                    xs,
                    ys,
                    peaks,
                    step,
                    peak,
                    float(slope),
                    float(offset),
                    rises_first,
                )
                # The two orderings agree whenever find_peaks offers no
                # more candidates than were asked for, which is most
                # scans; solving the same seeding twice is only slow.
                if seeds not in seedings:
                    seedings.append(seeds)
            # And, for a shape with a mixing, every one of those
            # seedings solved twice more: once with the mixing free and
            # once with it held at a half, keeping whichever scores
            # best. The free solves are tried first so that a tie keeps
            # the measured mixing over the assumed one.
            #
            # Scored against each other rather than simply trusted
            # because the free model contains the fixed one: nesting
            # promises that the better optimum exists, not that this
            # solver walked to it, and a freed mixing is one more
            # direction in which to get lost. The caution is the same
            # one the two seedings already get, and it is free -- where
            # the mixing is turned down the fixed answer is in hand.
            attempts = [None]
            # A fourth parameter a component asks a reading a component
            # more, on the same terms the minimum above is written to.
            if (
                mixed is not None
                and xs.size >= PEAK_FIT_MIN_POINTS + 4 * peaks - 3
            ):
                attempts.insert(0, mixed[0])
            best = None
            best_score = -np.inf
            for mix in attempts:
                for seeds in seedings:
                    got = _peak_solve(
                        xs,
                        ys,
                        peaks,
                        step,
                        peak,
                        seeds,
                        float(slope),
                        float(offset),
                        mix,
                    )
                    if got is None:
                        continue
                    # A fit whose R2 could not be formed -- the
                    # background alone accounts for the data -- is still
                    # a fit, and so beats nothing while losing to
                    # anything scored.
                    score = -np.inf if got[0] is None else got[0]
                    if best is None or score > best_score:
                        best, best_score = got, score
    except (RuntimeError, ValueError, np.linalg.LinAlgError):
        return stats
    if best is None:
        return stats
    r2, found, line_slope, line_offset = best

    components: list[dict[str, float | None]] = []
    for centre, fwhm, height, eta in found:
        component = _blank_stats()
        component.update(
            {
                "fit_center": centre,
                "fit_fwhm": fwhm,
                "fit_height": height,
                # Each component's own integrated intensity, on the same
                # terms peak_fit states: height x width x the shape's
                # integral over the whole line, the shared background
                # off. The integral is the mixing's where one was
                # fitted, so two components of one fit can take two
                # different factors -- they are two different shapes.
                "fit_area": height
                * fwhm
                * (area_factor if eta is None else mixed[1](eta)),
                # Blank unless this component's own mixing was fitted
                # and kept, which peak_fit says more about.
                "fit_eta": eta,
                "fit_r2": r2,
                "fit_slope": line_slope,
                "fit_offset": line_offset,
            }
        )
        components.append(component)
    stats["fit_slope"] = line_slope
    stats["fit_offset"] = line_offset
    stats["fit_peaks"] = components
    return stats


# How many points the drawn model is evaluated at. Enough that a narrow
# peak still reads as a curve rather than a polygon.
PEAK_FIT_PLOT_POINTS = 400


def _peak_shaped(
    peak, mixed, stats: dict[str, float | None], x: np.ndarray
) -> np.ndarray:
    """One fitted row's unit-height profile, its own mixing applied.

    *stats* is a whole fit or one of its components. ``fit_eta`` is
    there and not ``None`` only where a mixing was fitted and kept, so a
    shape that has none -- and a pseudo-Voigt whose mixing was refused
    and which fell back to the fixed half -- goes through *peak* exactly
    as it did before there was a mixing to fit.
    """
    eta = stats.get("fit_eta")
    if mixed is not None and eta is not None:
        return mixed[0](x, stats["fit_center"], stats["fit_fwhm"], eta)
    return peak(x, stats["fit_center"], stats["fit_fwhm"])


def peak_model(
    x: np.ndarray, shape: str, stats: dict[str, float | None]
) -> np.ndarray | None:
    """Evaluate over *x* the fit ``peak_fit`` returned in *stats*.

    ``None`` when there is no fit to draw, so the caller can ask without
    checking first. Uses ``.get`` throughout because a row that was not the
    one fitted -- the curve while dY/dX is shown, the derivative while it is
    not -- and any colour map never carry these keys at all.

    A ``peak_fit_multi`` result is drawn as the one summed model its
    components make, not as a line apiece: they were fitted together over
    a single background, and several lines in one colour over one curve
    would read as several fits rather than one.
    """
    entry = PEAK_SHAPES.get(shape)
    if entry is None:
        return None
    peak, _area_factor = entry
    mixed = PEAK_MIX.get(shape)
    components = stats.get("fit_peaks")
    if components:
        total = stats["fit_slope"] * x + stats["fit_offset"]
        for component in components:
            total = total + component["fit_height"] * _peak_shaped(
                peak, mixed, component, x
            )
        return total
    if stats.get("fit_center") is None:
        return None
    return (
        stats["fit_height"] * _peak_shaped(peak, mixed, stats, x)
        + stats["fit_slope"] * x
        + stats["fit_offset"]
    )


def peak_parts(
    x: np.ndarray, shape: str, stats: dict[str, float | None]
) -> list[np.ndarray]:
    """The pieces the fit in *stats* is made of, for "Parts" to draw.

    The fitted background first, then one line a component -- each drawn
    **on** that background rather than from zero, so a component lies
    where the data it is claiming lies. From zero they would be readable
    against each other and nowhere near the curve: scan 109's background
    is about 800 counts, so every component would float a screen below
    its own peak.

    Empty where there is nothing to add to what ``peak_model`` already
    drew. A single-component fit is its background plus its one peak, so
    the component line would be the summed model over again, in the same
    colour, and only the background is worth a line of its own.
    """
    entry = PEAK_SHAPES.get(shape)
    if entry is None or stats.get("fit_slope") is None:
        return []
    peak, _area_factor = entry
    mixed = PEAK_MIX.get(shape)
    base = stats["fit_slope"] * x + stats["fit_offset"]
    components = stats.get("fit_peaks") or []
    parts = [base]
    if len(components) > 1:
        for component in components:
            parts.append(
                base
                + component["fit_height"]
                * _peak_shaped(peak, mixed, component, x)
            )
    return parts


def combine_curves(
    x_a: np.ndarray,
    y_a: np.ndarray,
    x_b: np.ndarray,
    y_b: np.ndarray,
    subtract: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """``A + B``, or ``A - B``, with B read onto A's X by interpolation.

    Two scans of one thing do not step the same points. Scans 171 and 173
    are 109 energies each and differ by up to 4e-05 keV -- the readback,
    sampled at two different moments -- so adding them point by point would
    be adding readings taken at different energies.

    Only across the overlap. Left to itself :func:`numpy.interp` holds B's
    first and last readings flat right across A, drawing a sum out where
    there is no second scan at all; ``left`` and ``right`` make that stretch
    NaN and the line simply stops. On 171 and 173 that is the first point
    and the last, 171 starting above 173 and ending below it.

    The sort is what makes this correct rather than merely tidy.
    :func:`numpy.interp` documents ``xp`` as increasing and says nothing
    whatever when it is not -- and both of those scans run *downwards* in
    energy, as does any scan stepped from the top of an edge. Handed 173 as
    it comes off disk it returns 0 finite points out of 109, so the sum
    would be a blank curve on exactly the pair this was written for; without
    the NaN guards it is instead wrong by 64452 on values around 4e5, in
    silence. An averaged curve's X can also carry the NaN pad
    :func:`group_curve` leaves, which the same clean removes.
    """
    # The reader does not truncate X against Y, so neither pair may be taken
    # at face value; the shortest wins, as reduce_dichro does.
    n_a, n_b = min(x_a.size, y_a.size), min(x_b.size, y_b.size)
    x_a = np.asarray(x_a[:n_a], dtype=float)
    y_a = np.asarray(y_a[:n_a], dtype=float)
    x_b = np.asarray(x_b[:n_b], dtype=float)
    y_b = np.asarray(y_b[:n_b], dtype=float)
    good = np.isfinite(x_b) & np.isfinite(y_b)
    order = np.argsort(x_b[good], kind="stable")
    x_b, y_b = x_b[good][order], y_b[good][order]
    if x_b.size == 0:
        return x_a, np.full(x_a.shape, np.nan)
    onto = np.interp(x_a, x_b, y_b, left=np.nan, right=np.nan)
    return x_a, (y_a - onto if subtract else y_a + onto)


def combine_name(first: str, second: str, subtract: bool) -> str:
    """How a combination reads, in the menus and on the plot alike."""
    return f"{first} {'-' if subtract else '+'} {second}"


def _rms(values) -> float:
    """Root mean square of the finite values, NaN when there are none."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.sqrt(np.mean(values**2))) if values.size else float("nan")


def _safe_name(text) -> str:
    """Make ``text`` safe to drop into a suggested filename."""
    return re.sub(r"[^A-Za-z0-9.-]", "_", str(text))


def trend_points(
    curves: list[dict], against: str, peaks: int
) -> tuple[dict[str, list[tuple[float, dict | None]]], int]:
    """Where each fitted curve sits on the trend axis, and what it fitted.

    ``{series: [(x, that curve's fit, or None), ...]}``, together with how
    many curves had to be left out for want of an x. One series per Y
    field, or per field and component with *peaks* non-zero, named the way
    :meth:`HdfViewer._draw_trend` labels its lines.

    The whole fit rather than one number out of it, because the plot wants
    one quantity and the export beside it wants every quantity, and the
    two must not be able to disagree about which scan sat where. The
    caller picks what it reads out.

    *against* is :data:`TREND_X_SCAN` for the scan number, or the name of
    a baseline field every plotted scan recorded, which is what turns a
    row of repeats into the curve they were measured to draw -- scans 32
    to 39 against huber_euler_ana_eta being the case it was built for. The
    caller offers only fields the selection shares, but a curve that
    cannot supply one is dropped here anyway, since a point with no X has
    nowhere to go.
    """
    series: dict[str, list[tuple[float, dict | None]]] = {}
    skipped = 0
    for curve in curves:
        # An average or a combination keeps only its first member's path
        # -- both say so where they build it -- so it has no scan number
        # of its own and cannot be placed on this axis. The "combo" key
        # is the honest discriminator; curve["path"] looks usable and is
        # not.
        if "combo" in curve:
            skipped += 1
            continue
        path = Path(curve["path"])
        if against == TREND_X_SCAN:
            try:
                position = float(int(scan_number(path)))
            except ValueError:
                # scan_number falls back to the whole stem when a name
                # carries no digits at all.
                skipped += 1
                continue
        else:
            reading = baseline_values(path, (against,))
            if against not in reading:
                # Dropped rather than gapped, which is the opposite of
                # what a declined fit gets: a fit that failed still knows
                # where on this axis it belongs, so its gap lands
                # somewhere meaningful. A scan missing the X reading has
                # nowhere to put one.
                skipped += 1
                continue
            position = reading[against]
        fit = curve.get("fit")
        if not peaks:
            series.setdefault(curve["y_field"], []).append((position, fit))
            continue
        # Every component gets its point even when there is none to give
        # it: a Multi fit is declined whole, so a scan the gate turned
        # down holds no components at all, and leaving it out would join
        # the scans either side of it with a line. The series are built
        # from the count asked for rather than from what came back, so the
        # gap lands in each of them.
        found = (fit or {}).get("fit_peaks") or []
        for index in range(peaks):
            series.setdefault(
                f"{curve['y_field']} peak {index + 1}", []
            ).append((position, found[index] if index < len(found) else None))
    # By the position alone and never by what is beside it: two scans may
    # share one -- scan 39 returns eta to 0 and lands under scan 32 -- and
    # a fit is not a thing that can be compared with another fit. Stable,
    # so a tie keeps the order the curves were plotted in.
    for points in series.values():
        points.sort(key=lambda point: point[0])
    return series, skipped


def fit_export_columns(mixed: bool) -> list[tuple[str, str]]:
    """``(column name, stats key)`` for the fit export, in combo order.

    Every quantity the trend offers, since the export writes all of them
    at once rather than the one selected -- which is most of why it is
    there: a peak's numbers are read together, and exporting under each in
    turn would make four files that then have to be joined back up.

    Named by their stats keys without the prefix all of them share, so a
    column can be matched by eye to the table heading it came from. The
    mixing only where there is one to report, on the same terms as that
    table's own column -- see :func:`fit_columns`.
    """
    return [
        (key.removeprefix("fit_"), key)
        for key in TREND_QUANTITIES.values()
        if mixed or key != "fit_eta"
    ]


def _fit_column(points: list[tuple[float, dict | None]], key: str):
    """One fitted quantity down a trend's points, every gap a NaN.

    A declined fit and a quantity the shape does not have read the same
    here, as they do in the table: nothing was measured, and ``nan`` is
    what :func:`write_columns` writes for that.
    """
    return np.asarray(
        [
            np.nan
            if stats is None or stats.get(key) is None
            else float(stats[key])
            for _position, stats in points
        ],
        dtype=float,
    )


def export_table(
    curves: list[dict], tags: dict[str, str], overlay: bool, contrast: bool
) -> tuple[list[str], list[np.ndarray]]:
    """The columns an export writes: an X and a Y per curve, side by side.

    Curves that were plotted against the *same* X array get one X column
    between them, so three detectors of one scan, or the four contrasts of
    one reduction, come out as the ``energy,xas,xmcd,xmucd,xncd`` block
    xmucd_live.py writes rather than as four repeats of the same energies.

    The curves have to share the array itself, not merely equal values. The
    two are not the same question here: scan 96 split along its inner axis
    gives eight curves whose energies differ in the fifth decimal -- the
    readback, sampled at eight different moments -- and those are eight
    honest X columns. Comparing values would also weld together any two
    one-point scans, and would never collapse an averaged curve, whose X can
    carry a NaN pad and so is not even equal to itself. Reading the column
    once per scan and handing the same array to every field is what
    :meth:`HdfViewer._read_curves` already does, so identity answers exactly
    the question being asked: were these drawn against one sampling.

    Columns are returned at their own lengths, short and long together; it
    is :func:`write_columns` that squares them off. Nothing zips an X against
    its Y, because the reader does not truncate one to the other.
    """

    def name_of(curve: dict) -> str:
        if contrast:
            # The scan tag always, overlaid or not: a panel is named by its
            # contrast and the four of them are only told from another run's
            # by the scan, which is how the legend reads them too. The
            # curve's own tag rather than its path, which for a lock-in
            # group names only the first of the scans that made it.
            return f"{curve['tag']}: {CONTRAST_LABELS[curve['y_field']]}"
        # A curve made from several scans -- a sum, a difference, a group of
        # repeats averaged -- is named by all of them. Its path is only the
        # first, which is the one thing it must not be called: that name is
        # already taken by the column that scan writes itself.
        name = curve.get("combo")
        if name is None:
            name = (
                f"{tags[curve['path']]}: {curve['y_field']}"
                if overlay
                else (curve["y_field"])
            )
            if curve["tag"]:
                name += f" @ {curve['tag']}"
        return name + (" (right)" if curve["on_right"] else "")

    groups: dict[int, list[dict]] = {}
    for curve in curves:
        groups.setdefault(id(curve["x"]), []).append(curve)

    names: list[str] = []
    columns: list[np.ndarray] = []
    for members in groups.values():
        names.append(f"x[{name_of(members[0])}]")
        columns.append(np.asarray(members[0]["x"], dtype=float))
        for curve in members:
            name = name_of(curve)
            names.append(f"y[{name}]")
            columns.append(np.asarray(curve["y"], dtype=float))
            # Absent rather than None in anything stashed before the drawing
            # loop materializes it -- a map, say -- so ask rather than index.
            derivative = curve.get("dy_dx")
            if derivative is not None:
                names.append(f"dy[{name}]")
                columns.append(np.asarray(derivative, dtype=float))
    return names, columns


def write_columns(
    path: Path, names: list[str], columns: list[np.ndarray], header: list[str]
) -> None:
    """Write ``columns`` side by side as commented CSV.

    ``header`` lines come first, each prefixed with ``#``, so the file says
    which scans and which settings produced it -- the convention
    ``write_reduced`` in xmucd_live.py follows.

    Through the :mod:`csv` module and not a plain join, because a split
    curve's name carries the commas of the axes it was split along
    (``pr2_pzt_localdc +27.33, bmag_kepco_current +1.029 #0``) and would
    otherwise break the header into fields that do not line up with the
    numbers under them.

    Columns of different lengths are squared off with blanks at the bottom.
    A value that is NaN is written ``nan`` instead: a curve normalized by a
    monitor that read zero has real holes in it, and "nothing was measured
    here" must not look like "this curve stopped earlier".
    """
    rows = max((column.size for column in columns), default=0)
    with open(path, "w", newline="") as stream:
        for line in header:
            stream.write(f"# {line}\n")
        writer = csv.writer(stream)
        writer.writerow(names)
        for row in range(rows):
            writer.writerow(
                [
                    f"{column[row]:.8g}" if row < column.size else ""
                    for column in columns
                ]
            )


def _migrate_settings(settings: QtCore.QSettings):
    """Carry the remembered directory over from the old settings file, once.

    The one setting this viewer keeps used to live under an organisation
    named for another beamline's account, ``USER6IDB``, which has nothing to
    do with 4-ID and was only ever a string copied along with the code.
    Renaming it on its own would have quietly lost whichever directory was
    remembered, so the old file is read first and its value adopted.

    Then the old file goes, and its directory with it if nothing else is in
    there. Leaving an emptied one behind -- which is all ``clear()`` does --
    would keep the name in sight and leave something for a later version to
    wonder about. Gone is unambiguous, and there is nothing to lose by then:
    the only key still read has just been copied across.

    ``last_export_directory`` is deliberately not brought over. It dates from
    when the export folder was picked by dialog and remembered; exports now
    go to ``export/`` under the directory being browsed -- see
    :meth:`HdfViewer._export_directory` -- so nothing reads that key, and
    carrying it would only move the litter to a new address.
    """
    old = QtCore.QSettings("USER6IDB", "hdf5-viewer")
    if not old.allKeys():
        return
    if not settings.contains("last_directory"):
        remembered = old.value("last_directory", "")
        if remembered:
            settings.setValue("last_directory", remembered)
    stale = Path(old.fileName())
    try:
        stale.unlink()
        stale.parent.rmdir()
    except OSError:
        # Tidying the migrated settings file is best-effort: another session
        # may still hold it, or the directory may not be empty. Its values
        # were copied across above, so leaving it costs nothing.
        pass
    # Written out now rather than left to the event loop, so that a session
    # that ends badly in its first second does not lose what was just read
    # out of a file this has already deleted.
    settings.sync()


class HdfViewer(QtWidgets.QMainWindow):
    """Directory-of-HDF5-files browser with a primary-stream plot and baseline table."""

    def __init__(self, initial_directory: Path | None = None):
        super().__init__()
        self.setWindowTitle("4ID HDF5 File Viewer")
        self.resize(1400, 900)

        # ~/.config/POLAR/scan-viewer.conf on Linux, holding one key,
        # "last_directory" -- see _remember_directory for when it is written.
        self._settings = QtCore.QSettings("POLAR", "scan-viewer")
        _migrate_settings(self._settings)
        self._directory: Path | None = None
        self._fields: list[str] = []
        self._last_paths: tuple[str, ...] = ()
        self._image_windows: list[ImageViewerWindow] = []
        self._spectrum_windows: list[SpectrumViewerWindow] = []
        self._syncing_fields = False
        self._syncing_channels = False
        self._colorbars: list = []
        self._last_mode = ""
        # The X, Y and Norm ticks, and the kind of scan they were made for.
        # They are kept and re-applied while that kind of scan is what is
        # selected, and read afresh from hints when it is not -- see
        # _refresh_field_lists. Held here rather than read back off the table,
        # because the table is emptied whenever the selection is, and letting
        # go of a scan for a moment is not a change of mind about what to plot.
        self._kept: tuple[str, list[str], str] = ("", [], "")
        self._fields_kind: tuple = ()
        # Which stream the field table currently describes, which scan's hints
        # built it, and whether the choice of stream was the viewer's or the
        # user's.
        self._fields_stream = ""
        self._fields_leader = ""
        self._stream_is_auto = True
        self._manual_stream = ""
        self._last_x_field: str | None = None
        self._last_stream: str | None = None
        self._last_monitor_field: str | None = None
        self._last_log_norm: bool | None = None
        # Which scans are drawn against the right-hand scale, by str(path).
        # A hand choice, so it survives a redraw -- but it is pruned to the
        # selection on every one, so it holds only while the scan it was made
        # on is on screen. Unlike the ticks above, which name a column and so
        # can describe the next scan too, this names a scan: once that scan is
        # deselected there is nothing left for it to mean. Paths and not scan
        # numbers: moving to another directory of the same series would
        # otherwise carry the choice onto a different scan 171.
        self._right_paths: set[str] = set()
        self._last_right_scans: frozenset[str] = frozenset()
        # Which scans are drawn upside down, and which pairs are drawn summed
        # or differenced. Both by str(path) and both pruned to the selection
        # on every redraw, for the reasons given just above: they name scans
        # rather than columns, so they mean nothing once the scan is gone.
        self._invert_paths: set[str] = set()
        # (first, second, subtract) triples. A set rather than one choice at
        # a time: with four scans up, 171+173 and 175+177 are two comparisons
        # of one experiment, and the picker is a ticked menu either way.
        self._combine: set[tuple[str, str, bool]] = set()
        # Which of those go to the right-hand scale. Kept apart from
        # _right_paths because a combination is not a scan, and it needs that
        # scale more than either of its parents: 171 and 173 read 3.3e5..4.5e5
        # while their difference reads -4.9e4..-2.0e4, which on one scale is a
        # flat line below the axis.
        self._right_combos: set[tuple[str, str, bool]] = set()
        self._last_inverted: frozenset[str] = frozenset()
        self._last_combined: tuple = ()
        # Which channels the last reduction was drawn from -- part of the
        # same-view test, so that changing one autoscales the panels afresh.
        self._last_reduce_key: tuple = ()
        self._last_unit_key: tuple = ()
        # What the plot is currently showing, for Export to write out, or None
        # while it shows nothing -- see _stash_plot. Kept rather than read back
        # at export time: the curves are locals of _update_plot and die with
        # it, re-reading would reduce a dichro scan a second time, and a file
        # written since the draw would give a file that disagrees with the
        # screen it claims to be.
        self._plotted: dict | None = None

        self._build_ui()

        start_dir = initial_directory
        if start_dir is None:
            remembered = self._settings.value("last_directory", "")
            if remembered and Path(remembered).is_dir():
                start_dir = Path(remembered)
        if start_dir is not None and start_dir.is_dir():
            self._set_directory(start_dir)

    # -- UI construction ------------------------------------------------

    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QVBoxLayout(central)

        top = QtWidgets.QHBoxLayout()
        open_btn = QtWidgets.QPushButton("Open Directory…")
        open_btn.clicked.connect(self._choose_directory)
        refresh_btn = QtWidgets.QPushButton("Refresh")
        refresh_btn.clicked.connect(self._reload_file_list)
        self.image_btn = QtWidgets.QPushButton("Load Image…")
        self.image_btn.setEnabled(False)
        self.image_btn.clicked.connect(self._open_image_viewers)
        self.spectrum_btn = QtWidgets.QPushButton("Load Spectrum…")
        self.spectrum_btn.setEnabled(False)
        self.spectrum_btn.clicked.connect(self._open_spectrum_viewers)
        self.baseline_btn = QtWidgets.QPushButton("Baseline…")
        self.baseline_btn.setToolTip(
            "Show the baseline (before/after) values for one selected file."
        )
        self.baseline_btn.clicked.connect(self._show_baseline_window)
        self.export_btn = QtWidgets.QPushButton("Export…")
        self.export_btn.setEnabled(False)
        self.export_btn.setToolTip(
            "Write what is plotted to a file, which kind decided by the name\n"
            "given. A .csv, .dat or .txt name writes the numbers -- an X and\n"
            "a Y column per curve, with the scans, channels and settings in a\n"
            "commented header, so the file can be read without this window.\n"
            "A .png, .pdf or .svg name saves the picture instead, as it is on\n"
            "screen, zoom and all. A colour map saves only as a picture.\n"
            "Offered in an export folder beside the scans, named data_ or\n"
            "figure_ and the scan numbers.\n"
            "Only active while something is plotted."
        )
        self.export_btn.clicked.connect(self._export)
        self.dir_label = QtWidgets.QLabel("No directory selected")
        top.addWidget(open_btn)
        top.addWidget(refresh_btn)
        top.addWidget(self.image_btn)
        top.addWidget(self.spectrum_btn)
        top.addWidget(self.baseline_btn)
        top.addWidget(self.export_btn)
        top.addWidget(self.dir_label, 1)
        root.addLayout(top)

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        root.addWidget(splitter, 1)

        # Left: the files, each with what its scan was.
        left = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left)
        left_layout.addWidget(
            QtWidgets.QLabel("HDF5 files (ctrl/shift-click for several):")
        )
        self.file_table = QtWidgets.QTableWidget(0, len(FILE_COLUMNS))
        self.file_table.setHorizontalHeaderLabels(FILE_COLUMNS)
        self.file_table.verticalHeader().setVisible(False)
        self.file_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.file_table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.file_table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
        )
        file_header = self.file_table.horizontalHeader()
        # Interactive throughout: every column starts at a set width -- so a
        # long motor name does not push the others around -- and any of them
        # can be dragged wider or narrower. The last one takes up the slack.
        file_header.setStretchLastSection(True)
        for column, width in enumerate(FILE_COLUMN_WIDTHS):
            file_header.setSectionResizeMode(
                column, QtWidgets.QHeaderView.ResizeMode.Interactive
            )
            self.file_table.setColumnWidth(column, width)
        self.file_table.setWordWrap(False)
        self.file_table.setTextElideMode(QtCore.Qt.TextElideMode.ElideRight)
        row_header = self.file_table.verticalHeader()
        row_header.setMinimumSectionSize(16)
        row_header.setDefaultSectionSize(18)
        row_header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
        self.file_table.itemSelectionChanged.connect(self._on_selection_changed)
        left_layout.addWidget(self.file_table)

        # The field chooser sits under the file list and drives all three axes.
        left_splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        left_splitter.addWidget(left)
        left_splitter.addWidget(self._build_field_panel())
        left_splitter.setSizes([520, 420])
        splitter.addWidget(left_splitter)

        splitter.addWidget(self._build_primary_panel())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        # The scan listing gets about a third of the window.
        total = self.width() or 1400
        splitter.setSizes([total // 3, total - total // 3])

        # Built now but shown only on demand, in its own window.
        self._baseline_panel = self._build_baseline_panel()
        self._baseline_window: QtWidgets.QDialog | None = None

        self.status = self.statusBar()

    def _build_field_panel(self) -> QtWidgets.QWidget:
        """One row per field of the chosen stream, with X / Y / Norm ticks.

        X and Norm hold a single choice each -- ticking one clears the previous
        -- while any number of Y fields can be on at once.
        """
        panel = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)

        stream_row = QtWidgets.QHBoxLayout()
        stream_row.setContentsMargins(0, 0, 0, 0)
        stream_row.addWidget(QtWidgets.QLabel("Stream:"))
        self.stream_combo = QtWidgets.QComboBox()
        self.stream_combo.setToolTip(
            "Which of the scan's streams to plot.\n"
            "A dichro scan writes the XAS and XMCD it worked out to\n"
            "dichro_monitor, beside the raw counts in primary, and is\n"
            "shown on that stream unless you pick another."
        )
        self.stream_combo.currentTextChanged.connect(self._on_stream_changed)
        stream_row.addWidget(self.stream_combo, 1)
        layout.addLayout(stream_row)

        layout.addWidget(
            QtWidgets.QLabel("Fields — X (one), Y (any number), Norm (one):")
        )
        self.field_table = QtWidgets.QTableWidget(0, len(FIELD_COLUMNS))
        self.field_table.setHorizontalHeaderLabels(FIELD_COLUMNS)
        self.field_table.verticalHeader().setVisible(False)
        self.field_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.field_table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.NoSelection
        )
        self.field_table.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        for column in range(1, len(FIELD_COLUMNS)):
            self.field_table.horizontalHeader().setSectionResizeMode(
                column, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
            )
        row_header = self.field_table.verticalHeader()
        row_header.setMinimumSectionSize(16)
        row_header.setDefaultSectionSize(18)
        row_header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
        self.field_table.itemChanged.connect(self._on_field_item_changed)
        layout.addWidget(self.field_table)
        return panel

    def _build_primary_panel(self) -> QtWidgets.QWidget:
        panel = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(panel)

        self.file_label = QtWidgets.QLabel("No scan selected")
        self.file_label.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse
        )
        layout.addWidget(self.file_label)

        self.hkl_label = QtWidgets.QLabel("")
        self.hkl_label.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse
        )
        self.hkl_label.setWordWrap(True)
        self.hkl_label.setStyleSheet("color: gray;")
        self.hkl_label.setToolTip(
            "The diffractometer's own huber_euler_h/k/l, averaged over the scan.\n"
            "A spread is shown for any axis that moved -- on a rocking curve\n"
            "one of the three sweeps while the others hold still."
        )
        layout.addWidget(self.hkl_label)

        axis_row = QtWidgets.QHBoxLayout()
        self.axis_label = QtWidgets.QLabel("")
        axis_row.addWidget(self.axis_label, 1)
        self.grid_btn = QtWidgets.QPushButton("2D map")
        self.grid_btn.setCheckable(True)
        self.grid_btn.setChecked(True)
        self.grid_btn.setToolTip(
            "Draw a two-motor grid scan as a colour map -- one per ticked Y\n"
            "field, side by side, each with its own colour scale.\n"
            "Untick to plot it as an ordinary X-Y curve instead.\n"
            "On a dichro mesh scan, choose the dichro_monitor stream to map\n"
            "the XAS and the XMCD themselves: primary holds both\n"
            "polarizations at every position, so a map of it can only show\n"
            "one of them.\n"
            "Only active when the selected scan is a grid."
        )
        self.grid_btn.toggled.connect(self._update_plot)
        axis_row.addWidget(self.grid_btn)

        self.clip_btn = QtWidgets.QPushButton("Clip map scale")
        self.clip_btn.setCheckable(True)
        self.clip_btn.setToolTip(
            "Leave the most extreme cells of each colour map out of its colour\n"
            "scale, so that a few wild ones -- the corner where the beam ran\n"
            "off the sample, a detector glitch -- stop flattening everything\n"
            "else into one shade. Clipped cells saturate, and the colour bar\n"
            "grows an arrow at each end it cut. How much to leave out is the\n"
            "box alongside. The statistics below still report the true maximum\n"
            "and minimum.\n"
            "Only active while a colour map is shown."
        )
        self.clip_btn.toggled.connect(self._update_plot)
        axis_row.addWidget(self.clip_btn)

        self.clip_spin = QtWidgets.QDoubleSpinBox()
        self.clip_spin.setRange(0.0, COLOR_CLIP_MAX)
        self.clip_spin.setDecimals(2)
        self.clip_spin.setSingleStep(0.5)
        self.clip_spin.setSuffix(" %")
        self.clip_spin.setValue(COLOR_CLIP_PERCENT)
        self.clip_spin.setToolTip(
            "How much of each colour map to leave out at each end of its\n"
            "scale. A spike or two needs a fraction of a percent; a map that\n"
            "is still washed out at that can take a good deal more. The two\n"
            "ends are cut by the same amount, each map by its own cells.\n"
            "Only active while the map scale is being clipped."
        )
        # Re-plotting on every keystroke would redraw the rasters once per
        # digit typed, and "12" would clip at 1% on its way to 12%. Wait for
        # Enter, the arrows, or the focus leaving the box.
        self.clip_spin.setKeyboardTracking(False)
        self.clip_spin.valueChanged.connect(self._update_plot)
        axis_row.addWidget(self.clip_spin)

        axis_row.addWidget(QtWidgets.QLabel("Colour:"))
        self.cmap_combo = QtWidgets.QComboBox()
        self.cmap_combo.addItems([AUTO_COLOR_MAP, *COLOR_MAPS])
        self.cmap_combo.setToolTip(
            "Which colour map the 2D maps are drawn in. On Automatic each map\n"
            "picks for itself: a diverging scale for a signed channel like\n"
            "dichro_xmcd, where which side of zero a cell falls on is the\n"
            "result, and an ordinary one everywhere else. Naming a map here\n"
            "applies it to every map shown instead.\n"
            "Only active while a colour map is shown."
        )
        self.cmap_combo.currentTextChanged.connect(self._update_plot)
        axis_row.addWidget(self.cmap_combo)

        self.center_btn = QtWidgets.QPushButton("Center on 0")
        self.center_btn.setCheckable(True)
        self.center_btn.setChecked(True)
        self.center_btn.setToolTip(
            "Put zero in the middle of the colour scale, by widening it to\n"
            "whichever end reaches further. A diverging map is read by where\n"
            "its middle colour falls, and fitted to the data that lands\n"
            "wherever the range happens to centre rather than at zero --\n"
            "which draws signal as though there were none. Clipped cells\n"
            "still saturate and the bar keeps the arrows it has earned.\n"
            "Only active while a diverging colour map is shown."
        )
        self.center_btn.toggled.connect(self._update_plot)
        axis_row.addWidget(self.center_btn)

        layout.addLayout(axis_row)

        # Second row for the curve-side buttons. The one above holds the
        # label and the colour-map controls, which are what belongs with
        # the "what is on the X and Y axes" readout; everything from
        # "Remove spikes" on acts on the curves themselves, in the order it
        # applies to them. Split in two because one row of all of these
        # ran past the right edge of a 1400-pixel window -- 1674 pixels of
        # rigid widgets, see the note by :attr:`reduce_row` below.
        button_row = QtWidgets.QHBoxLayout()

        # First in the row because it is first in the order the buttons along
        # it apply: spikes go, then the repeats inside a scan are averaged,
        # then whole scans, then the inversions and the sums. A point the
        # beam never lit has to leave before anything averages over it.
        self.spike_btn = QtWidgets.QPushButton("Remove spikes")
        self.spike_btn.setCheckable(True)
        self.spike_btn.setToolTip(
            "Drop the single points where the beam was gone. The monitor\n"
            "reads a fraction of a percent of its usual there and every\n"
            "channel goes with it, so the point holds no measurement -- yet\n"
            "it plots like any other, and in the XMCD panels one such point\n"
            "takes a whole energy with it: scan 239's XMuCD spans 0.511 with\n"
            "it and 0.026 without. Removed points leave a gap rather than a\n"
            "guess, and the panels keep all their energies, the reduction\n"
            "averaging whichever points of that energy did have beam.\n"
            "Only runs of a single point go. Where the beam was down for\n"
            "longer -- scan 244 lost 265 points together -- the scan shows\n"
            "it, which is not a spike but the run.\n"
            "Only active while there is a monitor to judge the beam by."
        )
        self.spike_btn.toggled.connect(self._update_plot)
        button_row.addWidget(self.spike_btn)

        self.average_btn = QtWidgets.QPushButton("Average repeats")
        self.average_btn.setCheckable(True)
        self.average_btn.setToolTip(
            "Gather the steps that sat at the same setting into one curve\n"
            "each, averaging the repeats -- a magnet pattern of\n"
            "[1, -1, -1, 1] becomes two curves rather than four. A scan\n"
            "stepping two positioners gathers only steps that agree on both.\n"
            "Only active when the selected scan steps inner positioners."
        )
        self.average_btn.toggled.connect(self._update_plot)
        button_row.addWidget(self.average_btn)

        # Its own button rather than a second meaning for the one above. That
        # one gathers the steps within a scan and is live only on a scan that
        # steps inner positioners; this gathers whole scans and is live only
        # on several. They compose -- the repeats inside each scan collapse
        # first, and the groups left are then matched across the files.
        self.average_scans_btn = QtWidgets.QPushButton("Average scans")
        self.average_scans_btn.setCheckable(True)
        self.average_scans_btn.setToolTip(
            "Gather the selected scans that are one measurement repeated\n"
            "into a single curve each, drawn in place of the repeats. Two\n"
            "scans are the same measurement when they ran the same plan over\n"
            "the same motors, the sample sat in the same place at the same\n"
            "temperature, and they stepped the same points -- 174, 175 and\n"
            "176 look alike but move the sample 0.2 mm between them, so they\n"
            "stay three curves. A scan stopped part way is averaged over the\n"
            "rows it reached. The new curve is named after every scan in it.\n"
            "With the XMCD panels on it joins the repeats into one reduction\n"
            "instead, over all four panels at once: 239 to 242 are one\n"
            "measurement run four times and draw one set of panels rather\n"
            "than four on top of each other.\n"
            "Only active while two or more scans are plotted."
        )
        self.average_scans_btn.toggled.connect(self._update_plot)
        button_row.addWidget(self.average_scans_btn)

        # A QPushButton with a menu rather than image_viewer's QToolButton:
        # this row is six push buttons, and a tool button's frame and its own
        # arrow would make it ragged. setMenu() supplies the indicator.
        self.right_btn = QtWidgets.QPushButton("Right scale")
        self.right_btn.setToolTip(
            "Draw whole scans against the right-hand Y axis, so that two runs\n"
            "of one channel can be compared when an offset between them is\n"
            "wider than either one's own swing: sharing a scale they sit at\n"
            "two heights in a frame stretched to hold both, and on separate\n"
            "ones each fills the frame and the shapes superimpose. The menu\n"
            'lists the scans now selected, and any curve "Add" has made --\n'
            "a difference is the one that needs this most. Sending every\n"
            "scan leaves nothing to compare against and so does nothing.\n"
            "Only active while two or more scans are plotted as curves."
        )
        self._right_menu = QtWidgets.QMenu(self.right_btn)
        # Rebuilt on every opening, as image_viewer's view menu is: the
        # selection changes underneath it and the menu cannot then go stale.
        self._right_menu.aboutToShow.connect(self._populate_right_menu)
        self.right_btn.setMenu(self._right_menu)
        button_row.addWidget(self.right_btn)

        # Both of these are the same widget as "Right scale" above and for
        # the same reason: a menu of the scans now selected, rebuilt as it
        # opens, ticking whole scans rather than fields.
        self.invert_btn = QtWidgets.QPushButton("Invert")
        self.invert_btn.setToolTip(
            "Turn whole scans upside down -- y becomes -y, so a dip reads as\n"
            "a peak of the same size and the statistics below can find its\n"
            "centre and width, which on a dip they leave blank. The menu\n"
            "lists the scans now selected; the ones left unticked are drawn\n"
            "as they were. Negating throws a scan clear of the one it is\n"
            "meant to be compared with, so it takes the right-hand scale\n"
            'when nothing else is using it -- as a tick in "Right scale",\n'
            "to clear there if the shared scale was wanted after all.\n"
            "Only active while scans are plotted as curves."
        )
        self._invert_menu = QtWidgets.QMenu(self.invert_btn)
        self._invert_menu.aboutToShow.connect(self._populate_invert_menu)
        self.invert_btn.setMenu(self._invert_menu)
        button_row.addWidget(self.invert_btn)

        self.add_btn = QtWidgets.QPushButton("Add")
        self.add_btn.setToolTip(
            "Draw the sum, or the difference, of two scans as a third curve\n"
            "of its own, keeping both originals. The second is read onto the\n"
            "first one's X by interpolation -- two runs rarely step quite the\n"
            "same points -- and where they do not overlap it is left blank.\n"
            "Both orders of the subtraction are offered, so either sign can\n"
            'be had. A scan ticked in "Invert" is inverted first. A\n'
            "difference is far smaller than the scans it came from, so the\n"
            "new curve takes the right-hand scale when nothing else is\n"
            'using it; "Right scale" lists it either way, to move it back.\n'
            "Only active while two or more scans are plotted as curves."
        )
        self._add_menu = QtWidgets.QMenu(self.add_btn)
        self._add_menu.aboutToShow.connect(self._populate_add_menu)
        self.add_btn.setMenu(self._add_menu)
        button_row.addWidget(self.add_btn)

        self.derivative_btn = QtWidgets.QPushButton("Show dY/dX")
        self.derivative_btn.setCheckable(True)
        self.derivative_btn.setToolTip(
            "Add the slope of every curve, as a dashed line on a third axis\n"
            "at the right. The statistics below gain a row per derivative,\n"
            "under the curve it came from: on a knife edge or an absorption\n"
            "edge, where the curve itself has no peak to measure, that row is\n"
            "where the edge sits and how wide it is.\n"
            "While this is on, Fit works on the derivative instead of the\n"
            "curve -- which is what makes those edges fittable at all."
        )
        self.derivative_btn.toggled.connect(self._update_plot)
        button_row.addWidget(self.derivative_btn)

        self.integrate_btn = QtWidgets.QPushButton("Integrate")
        self.integrate_btn.setCheckable(True)
        self.integrate_btn.setToolTip(
            "Add the area under every curve to the statistics below, in two\n"
            'columns: "Area" as the curve stands, "Net area" above the\n'
            "straight line joining its first and last points. On a rocking\n"
            "curve measured several times over, the area is what compares\n"
            "the repeats -- scans 32 to 38 of TiSe2 spread 21% in peak\n"
            "height and 15% in area -- and the difference between the two\n"
            "columns is the sloping background, 9 to 15% there and not the\n"
            "same share in any two scans. Trapezoidal, over every point\n"
            'drawn, so a gap left by "Remove spikes" is spanned by a\n'
            "straight line rather than left out.\n"
            "Only active while the scans are drawn as curves or panels."
        )
        self.integrate_btn.toggled.connect(self._update_plot)
        button_row.addWidget(self.integrate_btn)

        layout.addLayout(button_row)

        # The reduction controls go on a row of their own rather than at the
        # end of the one above. That row is full and rigid: two labels,
        # eleven buttons, a spin box and a combo want 1674 pixels, and none
        # of them wraps, so every minimum equals its own sizeHint and the row
        # gives nothing back under pressure -- it already pushes the window
        # past the 1400 asked for at startup. Another button, two labels, two
        # combos and a second button would be some 400 more. They also read
        # as one unit: what to reduce, from which channels, which way up.
        reduce_row = QtWidgets.QHBoxLayout()
        self.xmcd_btn = QtWidgets.QPushButton("XMCD panels")
        self.xmcd_btn.setCheckable(True)
        self.xmcd_btn.setToolTip(
            "Reduce the scan to XAS, XMCD, XMuCD and XNCD against energy,\n"
            "four stacked panels, every selected scan overlaid in its own\n"
            "colour. XNCD is even in the field, so the sample cannot make it\n"
            "-- whatever shows there is the helicity average failing to\n"
            "cancel, and the same leakage is in the XMuCD above it, which is\n"
            "why both are drawn and the status bar gives the two as a pair.\n"
            "A lock-in scan needs the scan at the opposite field selected\n"
            "with it: the field is stepped between the two files, not inside\n"
            "either, and scans at one sample position and temperature are\n"
            "reduced together as one measurement.\n"
            "Only active when the scan records the retarder and the magnet,\n"
            "or the lock-in channels."
        )
        self.xmcd_btn.toggled.connect(self._update_plot)
        reduce_row.addWidget(self.xmcd_btn)

        # Dropdowns rather than the Y and Norm ticks: the reduction needs
        # exactly one channel over exactly one monitor, while Y is
        # multi-choice and Norm deliberately starts clear.
        reduce_row.addWidget(QtWidgets.QLabel("Detector:"))
        self.detector_combo = QtWidgets.QComboBox()
        self.detector_combo.setToolTip(
            "Which channel the absorption is measured on. It starts on the\n"
            "pair the scan hints, where both halves of that pair are recorded\n"
            "and have beam on them, and otherwise walks to the first pair that\n"
            "gives a finite absorption. The hints need that fallback: scan 96\n"
            "hints a monitor no stream records, and the early scans here hint\n"
            "one that read a flat zero for the whole run.\n"
            "Only active while the XMCD panels are shown."
        )
        self.detector_combo.currentTextChanged.connect(self._on_channel_changed)
        reduce_row.addWidget(self.detector_combo)

        reduce_row.addWidget(QtWidgets.QLabel("Monitor:"))
        self.monitor_combo = QtWidgets.QComboBox()
        self.monitor_combo.setToolTip(
            "Which channel the detector is taken against. Which one carries\n"
            "the beam varies between runs -- bI0 reads a flat zero on the\n"
            "early scans here, where gI0 is the live one -- and a dead monitor\n"
            "leaves four empty panels and a word in the status bar.\n"
            "Only active while the XMCD panels are shown."
        )
        self.monitor_combo.currentTextChanged.connect(self._on_channel_changed)
        reduce_row.addWidget(self.monitor_combo)

        # Starts unticked, which is a plain ratio in both the places it now
        # applies. One button and one state covering the panels and the Norm
        # tick alike: it is the same question in both -- how a detector and a
        # monitor become an absorption -- and two buttons that could disagree
        # about it would draw the panels in units the curves below them were
        # not in. The cost is that the panels no longer open on log(mon/det):
        # a transmission measurement is ticked into it once per session.
        self.transmission_btn = QtWidgets.QPushButton("log(mon/det)")
        self.transmission_btn.setCheckable(True)
        self.transmission_btn.setToolTip(
            "Take the absorption as log(mon/det), which is what a measurement\n"
            "through the sample needs. Untick for a plain det/mon instead, the\n"
            "form a fluorescence channel wants -- the same choice pr_setup\n"
            "makes with its transmission flag.\n"
            "Applies to the XMCD panels, and to the ordinary plot wherever a\n"
            "Norm field is ticked: each Y field is then drawn as\n"
            "log(Norm / Y) rather than Y / Norm, a colour map with it.\n"
            "With no Norm ticked this ticks the monitor the scan hints, so a\n"
            "transmission scan takes one press rather than two.\n"
            "Only active while the panels are shown, or there is a monitor to\n"
            "be had -- a Norm tick, or a hint to make one from."
        )
        self.transmission_btn.toggled.connect(self._on_transmission_toggled)
        reduce_row.addWidget(self.transmission_btn)
        reduce_row.addStretch(1)
        layout.addLayout(reduce_row)

        # A row of its own rather than more width on either of the two above:
        # the axis row already wants 1674 pixels against the 1400 the window
        # asks for, and the reduction row reads as one unit about the
        # reduction that fitting is no part of. The cost is the ~30 pixels of
        # height this takes permanently from the canvas, which is the trade
        # made knowingly -- the controls are two widgets wide and could not
        # go anywhere else without crowding a row that has nothing to give.
        fit_row = QtWidgets.QHBoxLayout()
        fit_row.addWidget(QtWidgets.QLabel("Fit:"))
        self.shape_combo = QtWidgets.QComboBox()
        self.shape_combo.addItems(list(PEAK_SHAPES))
        self.shape_combo.setToolTip(
            "Peak shape fitted to every curve, on a sloping background.\n"
            "Gauss falls away fastest, a Lorentzian carries the most in its"
            " tails, and squared sits between them. Pseudo-Voigt blends a"
            " Gauss and a Lorentzian at one width, for a peak whose core"
            " and tails disagree about which it is, and fits the"
            ' proportion as well: "Fit eta" is how much of it is'
            " Lorentzian, 1 all of it and 0 none.\n"
            "That column stands blank where a free mixing would not hold"
            " -- it starts at a half and is kept only if it settles clear"
            ' of either end, since a solver resting on "pure Gauss" or'
            ' "pure Lorentzian" has not measured a blend, and both of'
            " those are in this combo already. The fixed half-and-half"
            " shape is fitted alongside in any case and whichever scores"
            " better is the one drawn, so the extra parameter can only"
            " help: at two peaks scan 34 goes R2 0.9984 to 0.9994 with it"
            " and scan 109 keeps the 0.9131 it had without, to the"
            " figure.\n"
            "The shape barely moves the fitted centre -- 0.0005 degrees"
            " across all four on scan 32 -- but it does move the width and"
            " height, by some 15%.\n"
            "So it is a width-and-height choice far more than a position one.\n"
            "Only active while curves are plotted; a colour map has no axis"
            " to fit along."
        )
        self.shape_combo.currentTextChanged.connect(self._update_plot)
        fit_row.addWidget(self.shape_combo)
        self.fit_btn = QtWidgets.QPushButton("Fit")
        self.fit_btn.setCheckable(True)
        self.fit_btn.setToolTip(
            "Fit the chosen peak shape plus a sloping line to every curve"
            " shown, and draw the model over it.\n"
            "Adds the fitted centre, width, height and R2 to the table.\n"
            "The fitted width disagrees with the FWHM column, and both are"
            " right: FWHM is a half-max crossing of the curve as it stands,"
            " while the fit takes the sloping background out from under it"
            " first -- 0.0421 against 0.0447 on scan 32.\n"
            "The height is likewise measured above that background, so it"
            " sits below the curve's own maximum.\n"
            '"Fit area" is the integrated intensity of the peak alone,'
            " height x width x a constant of the shape, worked out over"
            " the whole line rather than summed over the points drawn. It"
            ' is therefore not "Net area" and will not match it -- 2219'
            " against 2791 on scan 32 -- and because it does not depend on"
            " how wide the scan was, it is the one to compare repeats"
            " with. On a Lorentzian 6 to 11% of it is tail beyond the ends"
            " of the scan, and on a pseudo-Voigt whatever share of that"
            " its fitted eta asks for -- 4 to 7% at a half.\n"
            "Read R2 first. A peak model fitted to something that is not a"
            " peak comes back looking perfectly reasonable and is not: on"
            " the edge in scan 171 all four shapes report a centre to five"
            " figures and a width of a few thousandths, and only the R2 of"
            " 0.11 says the model describes nothing.\n"
            "A fit too far gone even for that is rejected outright, which"
            " leaves the row blank and says so in the status bar.\n"
            "Expect that often on an alignment scan: a step, a plateau or"
            " a monotonic rise holds no peak to find, and a peak the scan"
            " stepped over and caught in a single point (scan 25) is"
            " declined as well, since one sample fixes no width.\n"
            "With dY/dX shown the fit moves onto the derivative: the curve's"
            " own row blanks and the model is drawn over the dashed line"
            " instead. That is what fits a step -- a knife edge through the"
            " beam, or an absorption edge -- since the step itself holds no"
            " peak and is declined, while its derivative is the beam profile"
            " or the edge. Scans 101 to 106 fit at R2 0.99 that way and not"
            " at all on the curve.\n"
            "The converse is a trap: the derivative of a peak is a lobe up"
            " then a lobe down, so the fit takes one lobe, puts the centre"
            " beside the true one, and flips the height's sign between"
            " repeats that are otherwise identical (scans 32 to 38, R2 0.52"
            " to 0.71). Read R2 and the sign together.\n"
            '"Multi" at two peaks is the way out of it: two lobes are two'
            " peaks, and fitting them as such takes those same scans to R2"
            " 0.96 to 0.99, with the heights coming back opposite in sign"
            " as the two lobes are.\n"
            "Scan 34 is the one exception, at 0.90, and for a reason worth"
            " knowing: its curve is itself two peaks, so its derivative"
            " has four lobes and two components cannot be all of them.\n"
            "Differentiating also amplifies noise, so a derivative too ragged"
            " for its own step size is declined outright -- scan 171 is.\n"
            "Only active while curves are plotted; a colour map has no axis"
            " to fit along."
        )
        self.fit_btn.toggled.connect(self._update_plot)
        fit_row.addWidget(self.fit_btn)
        # Beside the fit it modifies rather than further along the row, so
        # that the shape, the count and the button read as one control.
        self.multi_btn = QtWidgets.QPushButton("Multi")
        self.multi_btn.setCheckable(True)
        self.multi_btn.setToolTip(
            "Fit more than one peak at once, over a single shared sloping"
            " background.\n"
            "Three cases it was built for, all in TiSe2: two peaks that"
            " barely overlap (scan 109, which one peak cannot fit at all"
            " and two fit at R2 0.92), two that overlap badly (scan 34,"
            " whose second peak at 44.7406 is a shoulder on the flank of"
            " the first and shows no maximum of its own -- 0.982 to"
            " 0.992), and the derivative of a peak (scan 31 with dY/dX"
            " shown, 0.641 to 0.990).\n"
            'That last one repairs the trap the "Fit" button warns'
            " about: the derivative of a peak is a lobe up then a lobe"
            " down, and fitting two peaks fits both lobes instead of"
            " locking onto one -- scans 32 to 38 go from R2 0.52-0.71 to"
            " 0.96-0.99. Scan 34 is the exception at 0.90, being itself"
            " two peaks, so its derivative has four lobes and two"
            " components cannot cover them.\n"
            "Each peak gets a row of its own in the table, numbered from"
            " the left, and the curve's own row goes blank: a sum of two"
            " peaks has no one centre or width to put there.\n"
            "Leave it off unless the curve in front of you needs it."
            " Asked for two peaks a single peak is divided rather than"
            " refused, and the division looks perfectly reasonable: scan"
            " 32 comes back as a broad base of FWHM 0.0914 beneath a"
            " sharp core of 0.0347, where the peak is really 0.0447 and"
            " neither number is it.\n"
            "R2 cannot warn you of that and must not be read as if it"
            " could -- more peaks always fit better, so it only ever"
            " rises (scan 32, 0.978 on one peak and 0.998 on two). The"
            " heights and widths of the separate rows are what tell a"
            " real second peak from a divided one.\n"
            "Only active while a fit is being made."
        )
        self.multi_btn.toggled.connect(self._update_plot)
        fit_row.addWidget(self.multi_btn)
        self.peaks_combo = QtWidgets.QComboBox()
        self.peaks_combo.addItems(PEAK_COUNTS)
        self.peaks_combo.setToolTip(
            "How many peaks are fitted at once.\n"
            "There is no reading it off the data: the count is the one"
            " thing a peak finder cannot supply. On scan 34 it finds one"
            " peak where there are two, because the second is a shoulder"
            " and not a maximum, and on scan 109 it offers six candidates"
            " where there are two.\n"
            'Two answers all three of the cases "Multi" was built for,'
            " and three does not improve on them. Scan 34 reaches 0.997"
            " from 0.992 only by cutting its taller peak in two, and scan"
            " 109 is declined outright: the third component comes back as"
            " a dip between the two real peaks, and paying for it narrows"
            " the left one to 0.0197 against a step of 0.0250.\n"
            "Each peak costs three more parameters, so more of them need"
            " a longer scan: 7 points for one, then 3 more apiece.\n"
            "Changing the count mid-series moves every number, so a trend"
            " across scans has to be made at one setting -- on the knife"
            " edge, scan 101's width reads 0.0210 at one peak and 0.0234"
            " at two, which is wider than the drift being looked for.\n"
            'Only active while "Multi" is on.'
        )
        self.peaks_combo.currentTextChanged.connect(self._update_plot)
        fit_row.addWidget(self.peaks_combo)
        self.parts_btn = QtWidgets.QPushButton("Parts")
        self.parts_btn.setCheckable(True)
        self.parts_btn.setToolTip(
            "Draw what the fit is made of, as well as the fit itself.\n"
            'The sloping background goes in dash-dot, and with "Multi"'
            " on each peak goes in dots -- drawn on the background, not"
            " from zero, so a component sits over the part of the curve"
            " it is claiming. Scan 109's background is around 800 counts,"
            " and from zero its peaks would be nowhere near the data.\n"
            "This is how to see a component that is not really there."
            ' "Multi" divides a single peak rather than refusing it, and'
            " the division reads perfectly well as a sum -- it is the two"
            " dotted lines that show one of them lying flat along the"
            " background, or two of them stacked on the same bump. On"
            " scan 32 at two peaks they sit at the same place, a broad"
            " base under a sharp core; on scan 109 they stand over"
            " separate peaks, which is what a real pair looks like.\n"
            "A single peak adds the background alone: its one component"
            " plus that background is the fitted line already drawn.\n"
            "All of it is drawn in the curve's own colour and left out of"
            " the legend, which would otherwise grow four entries a"
            " curve.\n"
            "Only active while a fit is being made."
        )
        self.parts_btn.toggled.connect(self._update_plot)
        fit_row.addWidget(self.parts_btn)
        # On the same row as the fit it reads from, which is also the only
        # row with width to spare: these take it to 996 pixels of the 1400
        # the window asks for, where axis_row alone already wants 1674.
        # 709 of the 996 were here before the X axis combo; the 207 it adds
        # are why that combo is pinned to a width a few lines below rather
        # than left to grow to the longest field name it happens to hold,
        # which would put the row over the 1400 on its own.
        fit_row.addWidget(QtWidgets.QLabel("Plot:"))
        self.trend_combo = QtWidgets.QComboBox()
        self.trend_combo.addItems(list(TREND_QUANTITIES))
        self.trend_combo.setToolTip(
            "Which fitted quantity the trend plots, the combo after it"
            " choosing what against.\n"
            "Height is the peak above the fitted background, FWHM its width,"
            " centre where it sits, and area the integrated intensity of the"
            " peak alone.\n"
            "Area is the one to compare repeats with: it does not depend on"
            " how wide the scan was, where height alone misses a peak that"
            " spread out rather than weakened. Centre is the one to watch"
            " against a baseline axis -- a peak moving as something is"
            " driven is the measurement scans 32 to 39 were taken for.\n"
            "Eta is the pseudo-Voigt's fitted mixing, 1 a pure Lorentzian"
            " and 0 a pure Gauss. It is blank for every other shape, and"
            " for a peak whose mixing would not stand.\n"
            "The export at the end of this row writes all of them at once.\n"
            "Only active while the trend is shown."
        )
        self.trend_combo.currentTextChanged.connect(self._update_plot)
        fit_row.addWidget(self.trend_combo)
        fit_row.addWidget(QtWidgets.QLabel("vs"))
        self.trend_x_combo = QtWidgets.QComboBox()
        self.trend_x_combo.addItems([TREND_X_SCAN])
        # The one combo here whose entries are not written into the source,
        # and baseline names run to 37 characters, so it is held to a width
        # rather than left to size itself to its longest entry: it would
        # swallow the row and change width on every selection. The popup
        # still opens as wide as it needs, so nothing becomes unreadable.
        # ...WithIcon because Qt6 removed the plain AdjustToMinimumContentsLength
        # and this is the only survivor that honours setMinimumContentsLength.
        # It reserves width for an icon the combo does not have, which costs a
        # few pixels and changes nothing else.
        self.trend_x_combo.setSizeAdjustPolicy(
            QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.trend_x_combo.setMinimumContentsLength(16)
        self.trend_x_combo.setToolTip(
            "What the trend plots against: the scan number, or a baseline"
            " reading the whole selection shares.\n"
            "Scans 32 to 39 are the case it was built for, stepping"
            " huber_euler_ana_eta from 0 to 90 degrees. Against the scan"
            " number the fitted area is eight points in a row; against eta"
            " it is the curve those eight were measured to draw.\n"
            "Only fields that actually moved are offered, the most moved"
            ' first, and "moved" is judged against how far each drifted'
            " within a scan. Without that second half the list is headed by"
            " the beam diagnostics, which wander further than anything"
            " deliberate here and mean nothing on this axis.\n"
            "The list is built from the selection, so it changes with it. A"
            " selection holding nothing that moved offers the scan number"
            " alone, and a chosen field the next selection does not share"
            " falls back to it.\n"
            "Two scans sitting at the same value put two points on the same"
            " vertical -- scan 39 returns eta to 0 and lands under scan 32 --"
            " which is the data saying the series closed where it opened.\n"
            "Only active while the trend is shown."
        )
        self.trend_x_combo.currentTextChanged.connect(self._update_plot)
        fit_row.addWidget(self.trend_x_combo)
        self.trend_btn = QtWidgets.QPushButton("Plot")
        self.trend_btn.setCheckable(True)
        self.trend_btn.setToolTip(
            "Trend the chosen fitted quantity across the selected scans, in"
            " a second plot sharing the width with the first. The combo"
            " after the quantity chooses the axis: the scan number, or a"
            " baseline reading the selection shares.\n"
            "One line per Y field, over whatever the fit ran on -- so with"
            " dY/dX shown it trends the derivative's fit, and the axis says"
            " so.\n"
            'It plots the fitted numbers only, never the "Center"/"FWHM"'
            ' crossings or the "Net area" trapezoid: mixing a fitted width'
            " and a measured one into one line would compare two different"
            ' measurements. That is why it needs "Fit".\n'
            "Scans 101 to 127 are the case it was built for: a knife edge"
            " driven through the beam, whose derivative gives the beam size,"
            " trended over ten repeats to show it drifting between 0.0198"
            " and 0.0217.\n"
            "A scan whose fit was declined leaves a gap rather than a line"
            " drawn straight across it.\n"
            'With "Multi" on there is a line per peak as well as per'
            ' field, named "peak 1", "peak 2" and so on from the left.'
            " Numbering them by position rather than by height is what"
            " makes a line follow one peak through a set of repeats, even"
            " across the scan where the two swap which is the taller.\n"
            "Make the whole series at one count. Changing it partway moves"
            " every number, by more than the drift above: the same knife"
            " edge reads 0.0210 at one peak and 0.0234 at two.\n"
            "Averaged and combined curves are left out: they keep only their"
            " first member's file name, so they have no scan number of their"
            " own to sit at.\n"
            'Only active while curves are plotted and "Fit" is on.'
        )
        self.trend_btn.toggled.connect(self._update_plot)
        fit_row.addWidget(self.trend_btn)
        self.fit_export_btn = QtWidgets.QPushButton("Export…")
        self.fit_export_btn.setToolTip(
            "Write the fitted numbers of every selected scan to a file,"
            " one row per scan, against the axis chosen beside it.\n"
            "All of them, not only the quantity being plotted: centre,"
            " FWHM, height and area side by side, since those four are"
            " read together and four files would have to be joined back"
            " up. A pseudo-Voigt adds its fitted eta.\n"
            'With "Multi" on each peak gets its own set of columns,'
            " named as the trend names its lines.\n"
            'With "Integrate" on the curve\'s own trapezoids come too,'
            " as curve_area and curve_net_area -- one pair per Y field,"
            " since they measure the curve and not a peak of it.\n"
            "A scan whose fit was declined writes nan, as a gap in the"
            " plot beside it.\n"
            'The header is the one the "Export…" at the top writes, so'
            " the two files say the same thing about where the numbers"
            " came from. That one writes the curves themselves.\n"
            "Only active while the trend is shown."
        )
        self.fit_export_btn.clicked.connect(self._export_fits)
        fit_row.addWidget(self.fit_export_btn)
        # Both a combo and a button carry Qt's grow flag, so without this the
        # two would spread across the whole window. reduce_row ends the same
        # way; axis_row gets the same effect from its stretching label.
        fit_row.addStretch(1)
        layout.addLayout(fit_row)

        canvas_box = QtWidgets.QVBoxLayout()
        self.figure = Figure(figsize=(6, 4))
        self.canvas = FigureCanvasQTAgg(self.figure)
        # One plot normally, and one per Y field when a scan splits into
        # several curves per field and overlaying them all would be a thicket,
        # or when each field gets its own colour map -- see _set_panels, which
        # builds them on demand, stacked for curves and in a row for maps.
        self._panels: list[tuple] = []
        self._panels_across = False
        # Gives the standard matplotlib "zoom to rectangle" (drag to zoom),
        # pan, and home/back/forward view-history controls.
        self.toolbar = NavigationToolbar2QT(self.canvas, panel)
        self._set_panels(1)

        # The toolbar's built-in coordinate label lives inside the QToolBar, so
        # once the readout grows past the space left over Qt shunts it into the
        # overflow area and it vanishes -- which it did as soon as the labels
        # carried field names. Route the message to a label of our own, on the
        # row below, where it has the full width to itself.
        self.coord_label = QtWidgets.QLabel("")
        self.coord_label.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse
        )
        if hasattr(self.toolbar, "locLabel"):
            self.toolbar.locLabel.hide()
        self.toolbar.set_message = self.coord_label.setText

        # Beside the toolbar rather than in the button rows below: they only
        # change how the axes are drawn, as the toolbar's own tools do, and
        # those rows are already as wide as the window allows. Checkable, so
        # the state survives every redraw -- _apply_log_scales puts it back
        # on the fresh axes each draw makes.
        self.log_x_btn = QtWidgets.QPushButton("Log X")
        self.log_x_btn.setCheckable(True)
        self.log_x_btn.setToolTip(
            "Draw the X axis on a logarithmic scale.\n"
            "Points at or below zero are left out while it is on."
        )
        self.log_x_btn.toggled.connect(self._toggle_log_scales)
        self.log_y_btn = QtWidgets.QPushButton("Log Y")
        self.log_y_btn.setCheckable(True)
        self.log_y_btn.setToolTip(
            "Draw the Y axes -- the right-hand ones too -- on a logarithmic\n"
            "scale. Points at or below zero are left out while it is on."
        )
        self.log_y_btn.toggled.connect(self._toggle_log_scales)

        toolbar_row = QtWidgets.QHBoxLayout()
        toolbar_row.setContentsMargins(0, 0, 0, 0)
        toolbar_row.addWidget(self.toolbar)
        toolbar_row.addWidget(self.log_x_btn)
        toolbar_row.addWidget(self.log_y_btn)
        toolbar_row.addWidget(self.coord_label, 1)
        canvas_box.addLayout(toolbar_row)
        canvas_box.addWidget(self.canvas, 1)
        canvas_widget = QtWidgets.QWidget()
        canvas_widget.setLayout(canvas_box)

        # A figure of its own, deliberately not routed through _set_panels:
        # that helper clears self.figure and owns a single set of axes,
        # twins, colour bars and a default format_coord, so a second caller
        # would tear down the first's. This one is built once and its axis
        # cleared per draw.
        trend_box = QtWidgets.QVBoxLayout()
        self.trend_figure = Figure(figsize=(4, 4))
        self.trend_canvas = FigureCanvasQTAgg(self.trend_figure)
        self.trend_ax = self.trend_figure.add_subplot(111)
        self.trend_toolbar = NavigationToolbar2QT(self.trend_canvas, panel)
        # This one keeps matplotlib's own coordinate label, where the main
        # toolbar's is rerouted to a label below it. The reason there was
        # that readouts carrying field names outgrew the QToolBar and Qt
        # hid them in the overflow; a trend reads "scan number, value" and
        # fits. Its save button is also the only way to write this figure
        # out -- Export's savefig takes the main one alone.
        trend_box.addWidget(self.trend_toolbar)
        trend_box.addWidget(self.trend_canvas, 1)
        self.trend_widget = QtWidgets.QWidget()
        self.trend_widget.setLayout(trend_box)
        self.trend_widget.setVisible(False)

        # Equal stretch so that "sharing the space equally" survives a window
        # resize; the sizes themselves are set when the panel is first shown,
        # once rather than on every draw, so that a drag of the handle is not
        # undone by the next redraw.
        self.plot_split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.plot_split.addWidget(canvas_widget)
        self.plot_split.addWidget(self.trend_widget)
        self.plot_split.setStretchFactor(0, 1)
        self.plot_split.setStretchFactor(1, 1)
        layout.addWidget(self.plot_split, 1)

        layout.addWidget(
            QtWidgets.QLabel("Peak statistics (per plotted curve):")
        )
        self.stats_table = QtWidgets.QTableWidget(0, len(STATS_COLUMNS))
        self.stats_table.setHorizontalHeaderLabels(STATS_COLUMNS)
        self.stats_table.horizontalHeader().setStretchLastSection(True)
        self.stats_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        vertical_header = self.stats_table.verticalHeader()
        vertical_header.setMinimumSectionSize(16)
        # Compact rows, so that the few the table is allowed take as little
        # off the plot as possible.
        vertical_header.setDefaultSectionSize(18)
        # defaultSectionSize only wins when the rows are not auto-sized.
        vertical_header.setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.Fixed
        )
        self._size_stats_table(0)
        layout.addWidget(self.stats_table)
        return panel

    def _build_baseline_panel(self) -> QtWidgets.QWidget:
        panel = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(panel)

        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("Baseline (before/after) for file:"))
        self.baseline_combo = QtWidgets.QComboBox()
        self.baseline_combo.currentIndexChanged.connect(self._update_baseline)
        top.addWidget(self.baseline_combo, 1)
        layout.addLayout(top)

        self.baseline_table = QtWidgets.QTableWidget(0, 4)
        self.baseline_table.setHorizontalHeaderLabels(
            ["Field", "Before scan", "After scan", "Δ (after − before)"]
        )
        self.baseline_table.horizontalHeader().setStretchLastSection(True)
        self.baseline_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.baseline_table.setSortingEnabled(True)
        layout.addWidget(self.baseline_table)
        return panel

    # -- Directory / file list -------------------------------------------

    def _choose_directory(self):
        start = str(self._directory) if self._directory else ""
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Choose a directory of HDF5 files", start
        )
        if directory:
            self._set_directory(Path(directory))
            # Remembered here as well as on close, so that a session killed
            # from the terminal, or one that never ends tidily, still comes
            # back where it was asked to go.
            self._remember_directory()

    def _set_directory(self, directory: Path):
        """Browse *directory*, without that alone deciding what is remembered.

        Nothing is written here. This runs for the directory handed in on the
        command line and for the one a test harness constructs the window
        with, and either writing itself down would replace what the last real
        session chose: start the viewer once from a script pointed at
        somewhere else and the next plain start would open there instead.
        """
        self._directory = directory
        self.dir_label.setText(str(directory))
        self._reload_file_list()

    def _remember_directory(self):
        """Store the directory being browsed as the one to open next time."""
        if self._directory is not None:
            self._settings.setValue("last_directory", str(self._directory))

    def closeEvent(self, event):
        """Remember the directory in view, rather than the last one opened.

        Which is not the same thing once a second window is in play. Writing
        at open time instead meant whichever viewer was *started* last had the
        final say for good: open one on lineup and another on tests, close the
        second, work on in the first for an hour, and the next start still
        came up in tests, because the first window never wrote anything again.
        Closing is the moment that actually says where the work was.
        """
        self._remember_directory()
        super().closeEvent(event)

    def _reload_file_list(self):
        """List the directory, filling in each file's scan details as we go.

        The scan columns are read for every file up front rather than on
        selection, so the listing can be scanned for the run you want. Reads are
        cached per (path, mtime, size), so a Refresh only touches new or changed
        files.
        """
        if self._directory is None:
            self.file_table.setRowCount(0)
            return
        selected_names = {path.name for path in self._selected_paths()}
        self.file_table.blockSignals(True)
        self.file_table.setRowCount(0)
        # Only the scans: the detector image files carry no scan metadata, and
        # their frames are reached through "Load Image…" anyway. The listing is
        # not recursive either, so a per-detector subdirectory stays out of it.
        files = [
            p
            for p in list_hdf_files(self._directory)
            if read_scan_info(p)["is_scan"]
        ]
        self.file_table.setRowCount(len(files))
        restore: list[int] = []
        for row, path in enumerate(files):
            for column, text in enumerate(self._scan_row(path)):
                item = QtWidgets.QTableWidgetItem(text)
                if column == 0:
                    item.setData(QtCore.Qt.ItemDataRole.UserRole, str(path))
                elif column in (IMAGE_COLUMN, SPECTRUM_COLUMN):
                    item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.file_table.setItem(row, column, item)
            if path.name in selected_names:
                restore.append(row)
        self._select_rows(restore)
        self.file_table.blockSignals(False)
        self.status.showMessage(f"{len(files)} scan(s) in {self._directory}.")
        self._on_selection_changed()

    def _select_rows(self, rows: list[int]):
        """Select every one of *rows*.

        Not ``selectRow`` in a loop: in ExtendedSelection mode each call
        *replaces* the selection, so restoring a multi-file selection after a
        Refresh would keep only the last row.
        """
        model = self.file_table.model()
        last_column = self.file_table.columnCount() - 1
        selection = QtCore.QItemSelection()
        for row in rows:
            selection.select(model.index(row, 0), model.index(row, last_column))
        self.file_table.selectionModel().select(
            selection,
            QtCore.QItemSelectionModel.SelectionFlag.Select
            | QtCore.QItemSelectionModel.SelectionFlag.Rows,
        )

    @staticmethod
    def _scan_row(path: Path) -> list[str]:
        """One file's cells: its scan number, then what that scan was."""
        info = read_scan_info(path)
        motors = ", ".join(name for name, *_ in info["ranges"]) or ", ".join(
            info["motors"]
        )
        span = ", ".join(
            f"{start:.6g} → {stop:.6g}"
            for _name, start, stop, _n in info["ranges"]
        )
        return [
            scan_number(path),
            HAS_FILE if resolve_image_file(path) is not None else NO_FILE,
            HAS_FILE if resolve_spectrum_file(path) is not None else NO_FILE,
            info["plan"],
            motors,
            span,
            "x".join(str(n) for n in info["shape"])
            if len(info["shape"]) == 2
            else str(info["points"]),
            "" if info["count_time"] is None else f"{info['count_time']:.6g}",
        ]

    def _selected_paths(self) -> list[Path]:
        paths = []
        for row in sorted(
            {index.row() for index in self.file_table.selectedIndexes()}
        ):
            item = self.file_table.item(row, 0)
            if item is not None:
                paths.append(Path(item.data(QtCore.Qt.ItemDataRole.UserRole)))
        return paths

    # -- Selection handling -----------------------------------------------

    def _on_selection_changed(self):
        paths = self._selected_paths()
        self._update_file_label(paths)
        self._update_hkl_label(paths)
        self._update_image_button(paths)
        self._update_spectrum_button(paths)
        if not paths:
            # Only what describes the table, which is about to be emptied.
            # The ticks and the kind they belong to are left alone: selecting
            # the same kind of scan again puts them straight back.
            self._fields = []
            self._fields_stream = ""
            self._syncing_fields = True
            self.field_table.setRowCount(0)
            self._syncing_fields = False
            self.stream_combo.blockSignals(True)
            self.stream_combo.clear()
            self.stream_combo.blockSignals(False)
            self.baseline_combo.clear()
            self._refresh_reduce_channels([])
            self._update_plot()
            self._update_baseline()
            return
        self._refresh_stream_combo(paths)
        self._refresh_field_lists(paths, self.stream_combo.currentText())
        self._refresh_reduce_channels(paths)
        self._refresh_baseline_combo(paths)
        self._auto_check_xmcd()
        self._update_plot()
        self._update_baseline()

    def _auto_check_xmcd(self):
        """Pre-tick XMCD when the stream on screen is ``dichro_monitor``.

        That stream is where the panels or the maps are the point of the
        scan -- primary carries the raw retarder-stepped points, which
        plotted as they stand are a 1D view of the measurement before the
        split. A user who picked that stream wanted the split, so opening
        on the 1D view and asking them to tick again is a round trip.

        Does nothing when XMCD is already on, when the stream is anything
        else, or when the selection holds no scan the split can run on.
        Signals are blocked round the tick: the caller is about to run
        :meth:`_update_plot` already, and two redraws for one gesture
        would strobe the plot area.
        """
        if (
            self.stream_combo.currentText() != DICHRO_STREAM
            or self.xmcd_btn.isChecked()
        ):
            return
        if not any(
            reduce_channels(path) or has_dichro_stream(path)
            for path in self._selected_paths()
        ):
            return
        self.xmcd_btn.blockSignals(True)
        try:
            self.xmcd_btn.setChecked(True)
        finally:
            self.xmcd_btn.blockSignals(False)

    def _refresh_stream_combo(self, paths: list[Path]):
        """Offer the streams the selection holds and settle on one to plot.

        A dichro scan leads with its reduced ``dichro_monitor`` stream, since
        the XAS and XMCD there are what the scan was run for; anything else
        leads with ``primary``. A stream picked by hand is kept as long as the
        selection still offers it.
        """
        names: list[str] = []
        for path in paths:
            for name in scan_streams(path):
                if name not in names:
                    names.append(name)
        if DICHRO_STREAM in names:
            wanted = DICHRO_STREAM
        else:
            wanted = names[0] if names else ""
        if not self._stream_is_auto and self._manual_stream in names:
            wanted = self._manual_stream
        self.stream_combo.blockSignals(True)
        try:
            self.stream_combo.clear()
            self.stream_combo.addItems(names)
            if wanted:
                self.stream_combo.setCurrentText(wanted)
        finally:
            self.stream_combo.blockSignals(False)
        # Nothing to choose between, so nothing to choose with.
        self.stream_combo.setEnabled(len(names) > 1)

    def _on_stream_changed(self, name: str):
        """Follow a stream picked by hand, and stay on it across scans."""
        self._stream_is_auto = False
        self._manual_stream = name
        paths = self._selected_paths()
        if not paths:
            return
        self._refresh_field_lists(paths, name)
        self._auto_check_xmcd()
        self._update_plot()

    def _refresh_field_lists(self, paths: list[Path], stream: str):
        # Union of the stream's fields across every selected file, in the first
        # file's order (which, absent track_order, is alphabetical). A
        # selection with no stream to speak of -- a file holding none the
        # viewer can plot -- leaves the table empty.
        fields: list[str] = []
        seen: set[str] = set()
        if stream:
            for path in paths:
                for name in stream_fields(path, stream):
                    if name not in seen:
                        fields.append(name)
                        seen.add(name)

        # The scan that leads the selection decides all three columns -- the
        # axis it scanned for X and its hinted detectors for Y, with Norm
        # clear unless the stream asks (see default_monitor_field) -- but it
        # decides them only when the kind of scan changes. Within one kind the
        # ticks stand: a run of lineup scans is stepped through one at a time,
        # and having to re-tick the detector on each would make the comparison
        # they were taken for unworkable.
        #
        # What the hints are there for is the move between kinds. A tick made
        # on a mesh scan governing the energy scan opened after it was the old
        # behaviour here, and the trouble with it was that a field ticked an
        # hour ago, on a scan since forgotten, decided what the scan now on
        # screen showed. scan_kind draws the line at the hints themselves: the
        # ticks stand until the scan on screen asks for something the one
        # before it did not.
        #
        # Bringing a second scan alongside the first is not a move between
        # kinds either, even when it is a different kind of scan: the scan the
        # ticks were made for is still on screen, and the whole point of
        # adding the second one is to see the same channels on both. So the
        # ticks also stand while the scan that seeded them is still selected,
        # whichever end the new one was added at.
        kind = scan_kind(paths[0], fields, stream)
        preferred_x, hinted_y, preferred_norm = kind[2], kind[3], kind[4]
        preferred_y = list(hinted_y)
        keeping = kind == self._fields_kind or self._fields_leader in {
            str(path) for path in paths
        }
        if keeping:
            # What was ticked wins, as long as the selection still records it:
            # the kept ticks outlive the table they were made on, and a second
            # scan can be dropped as well as added, so a tick whose field went
            # with it has to give way to the hint. An empty Norm is a choice
            # and stays empty; an empty Y is only ever what is left when every
            # ticked field has gone.
            kept_x, kept_y, kept_norm = self._kept
            live_y = [name for name in kept_y if name in seen]
            preferred_x = kept_x if kept_x in seen else preferred_x
            preferred_y = live_y if live_y or not kept_y else preferred_y
            preferred_norm = kept_norm if kept_norm in seen else ""
        else:
            # Only a re-seed records the kind, so that the ticks stay attached
            # to the kind of scan they were made for rather than to a mixed
            # selection they merely survived.
            self._fields_kind = kind
            self._fields_leader = str(paths[0])

        if fields == self._fields and stream == self._fields_stream:
            # The same fields as the last scan, so the table itself can stand
            # -- but the ticks on it belonged to that scan, and this need not
            # be one of the same kind. An X the new scan does not name is the
            # exception: clearing it would leave nothing to plot against, so
            # it keeps the one it has.
            if preferred_x and preferred_x != self._x_field():
                self._set_checked(FIELD_X, [preferred_x])
            if preferred_y != self._checked_y_fields():
                self._set_checked(FIELD_Y, preferred_y)
            if preferred_norm != self._monitor_field():
                self._set_checked(FIELD_NORM, [preferred_norm])
            self._remember_fields()
            return
        self._fields = fields
        self._fields_stream = stream

        # A new field set, so the table is rebuilt and the ticks go on from
        # scratch. A preference read from hints already names one of *fields*,
        # having been resolved against it; one carried over from the previous
        # selection need not, since the scan that recorded it may be the one
        # just dropped, so both are filtered here rather than only the second.
        current_x = preferred_x if preferred_x in seen else ""
        current_y = {name for name in preferred_y if name in seen}
        current_norm = preferred_norm if preferred_norm in seen else ""

        self._syncing_fields = True
        try:
            self.field_table.setRowCount(0)
            self.field_table.setRowCount(len(fields))
            for row, name in enumerate(fields):
                label = QtWidgets.QTableWidgetItem(name)
                self.field_table.setItem(row, 0, label)
                states = (
                    (FIELD_X, name == current_x),
                    (FIELD_Y, name in current_y),
                    (FIELD_NORM, name == current_norm),
                )
                for column, checked in states:
                    item = QtWidgets.QTableWidgetItem()
                    item.setFlags(
                        (item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
                        & ~QtCore.Qt.ItemFlag.ItemIsSelectable
                    )
                    item.setCheckState(
                        QtCore.Qt.CheckState.Checked
                        if checked
                        else QtCore.Qt.CheckState.Unchecked
                    )
                    item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                    self.field_table.setItem(row, column, item)
        finally:
            self._syncing_fields = False
        self.field_table.resizeColumnsToContents()
        self.field_table.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        self._remember_fields()

    def _remember_fields(self):
        """Hold on to the ticks, so the next scan of this kind opens on them."""
        self._kept = (
            self._x_field(),
            self._checked_y_fields(),
            self._monitor_field(),
        )

    def _refresh_reduce_channels(self, paths: list[Path]):
        """Offer the reducible scan's channels and settle on a detector and monitor.

        A pair chosen by hand stands for as long as it still works: every
        redraw goes through _update_plot rather than through here, so picking
        a channel and then toggling dY/dX or adding a scan alongside keeps it.
        It gives way only when a newly selected scan does not record it, or
        records it dead -- and that is not a nicety either, since the pair
        that reduces is a fact about the run. The scans here before 33 need
        bI1/gI0 and 96 onwards bI1/bI0, and carrying one over to the other
        leaves four empty panels with nothing on screen to say why.

        A looser rule than the one the X, Y and Norm ticks follow, which lets
        go of a hand choice as soon as the scan selected asks for something
        else (see :func:`scan_kind`). A pair here is held on to past that,
        because what a scan asks for is no guide at all to which channels
        reduce -- see :func:`reduce_defaults`.
        """
        # The first scan that can be reduced, not paths[0]: selecting a plain
        # scan alongside a dichro one should still offer the dichro one's
        # channels, and which of the two leads is only where they happen to
        # sit in the listing.
        leader = next((path for path in paths if reduce_channels(path)), None)
        channels = reduce_channels(leader) if leader is not None else []
        transmission = self.transmission_btn.isChecked()
        detector = self.detector_combo.currentText()
        monitor = self.monitor_combo.currentText()
        if channels and not self._pair_works(
            leader, detector, monitor, transmission
        ):
            detector, monitor = reduce_defaults(leader, channels, transmission)
        self._syncing_channels = True
        try:
            for combo, choice in (
                (self.detector_combo, detector),
                (self.monitor_combo, monitor),
            ):
                combo.clear()
                combo.addItems(channels)
                if choice in channels:
                    combo.setCurrentIndex(channels.index(choice))
        finally:
            self._syncing_channels = False

    @staticmethod
    def _pair_works(
        path, detector: str, monitor: str, transmission: bool
    ) -> bool:
        """Whether this scan records both channels and they give any absorption."""
        if path is None or not detector or not monitor or detector == monitor:
            return False
        try:
            data = read_stream(path, [detector, monitor], PRIMARY_STREAM)
        except OSError:
            return False
        if detector not in data or monitor not in data:
            return False
        return bool(
            np.isfinite(
                absorption(data[detector], data[monitor], transmission)
            ).any()
        )

    def _on_channel_changed(self, _name: str = ""):
        """Redraw for a channel picked by hand, not for one the viewer just set."""
        if self._syncing_channels:
            return
        self._update_plot()

    # -- Field choices -------------------------------------------------------

    def _checked_in_column(self, column: int) -> list[str]:
        names = []
        for row in range(self.field_table.rowCount()):
            tick = self.field_table.item(row, column)
            label = self.field_table.item(row, 0)
            if (
                tick is not None
                and label is not None
                and tick.checkState() == QtCore.Qt.CheckState.Checked
            ):
                names.append(label.text())
        return names

    def _set_checked(self, column: int, names: list[str]):
        """Tick exactly *names* in *column*, without counting as a hand choice."""
        wanted = {name for name in names if name}
        self._syncing_fields = True
        try:
            for row in range(self.field_table.rowCount()):
                tick = self.field_table.item(row, column)
                label = self.field_table.item(row, 0)
                if tick is None or label is None:
                    continue
                tick.setCheckState(
                    QtCore.Qt.CheckState.Checked
                    if label.text() in wanted
                    else QtCore.Qt.CheckState.Unchecked
                )
        finally:
            self._syncing_fields = False

    def _x_field(self) -> str:
        names = self._checked_in_column(FIELD_X)
        return names[0] if names else ""

    def _checked_y_fields(self) -> list[str]:
        return self._checked_in_column(FIELD_Y)

    def _monitor_field(self) -> str:
        """The normalization field, or ``""`` for none."""
        names = self._checked_in_column(FIELD_NORM)
        return names[0] if names else ""

    def _hinted_monitor(self) -> str:
        """The monitor the leading scan nominates, if the table offers it.

        What the Norm tick would start on, were it not that
        :func:`default_monitor_field` passes the hint over on purpose so that
        a scan opens showing what it recorded. That policy stands; this reads
        the hint for the one place a monitor is wanted without a tick -- the
        log button, which has no "mon" to take a log against until something
        names one.

        Asked of the leading scan rather than of all of them, and of the names
        the field table was built from rather than of the file, the same two
        economies :meth:`_spike_monitor` makes: it is asked once per redraw,
        and :func:`stream_fields` is not cached.

        And, like that one, falling back to :data:`REDUCE_MONITORS` where the
        hint names nothing the stream carries. This method first refused to,
        on the grounds that a guess is worth making about which channel judges
        the beam and not about which one every number on the plot gets divided
        by. The scans decided it: 279 lineup scans on TiSe2 and 7 here hint
        ``Time``, which is not a field of theirs -- the stream calls them
        ``Time1`` and ``Time2`` -- so the button went grey on 286 scans for
        want of a name, and the complaint it drew was that lineup scans did
        not work. They are not what is special; the hint is.

        What makes the guess safe here is that it is not a guess the user
        cannot see. :meth:`_on_transmission_toggled` writes it into the Norm
        column as a tick, where a wrong channel is on screen and one click
        from being moved -- where :meth:`_spike_monitor` judges the beam
        silently and has to be right. A monitor that can be corrected by
        looking at it is worth offering; the alternative was offering nothing.
        """
        paths = self._selected_paths()
        if not paths:
            return ""
        hint = read_scan_info(paths[0])["hints"]["monitor"]
        if hint in self._fields:
            return hint
        return next(
            (name for name in REDUCE_MONITORS if name in self._fields), ""
        )

    def _on_transmission_toggled(self, checked: bool):
        """Supply the monitor the log needs, where the press did not have one.

        The button is live wherever a monitor can be had, not only where Norm
        is already ticked, because the scans it exists for open with Norm
        clear: a transmission run hints its monitor and the field table starts
        it unticked, so the button sat grey on exactly the measurement it was
        added for, and a control that is dead when wanted is one that cannot
        be found. Pressing it now ticks the hinted monitor itself.

        The tick is made rather than implied. It goes into the Norm column
        where it can be seen, moved to another channel or taken off again, so
        that the field dividing every number on the plot is never one with
        nothing on screen naming it. A Norm already ticked is left alone --
        a choice outranks a hint, and the press was about the form, not the
        field -- and unticking the button changes the form back without
        disturbing the tick it made, which is what lets the two be compared.
        """
        if checked and not self._monitor_field():
            hinted = self._hinted_monitor()
            if hinted:
                self._set_checked(FIELD_NORM, [hinted])
                # Kept as a hand choice is kept: the press was one, and a tick
                # that went away at the next scan of the same kind would take
                # the log the user asked for away with it.
                self._remember_fields()
        self._update_plot()

    def _on_field_item_changed(self, item: QtWidgets.QTableWidgetItem):
        if self._syncing_fields:
            return
        column = item.column()
        # X and Norm are single-choice: ticking one clears whatever held it.
        if (
            column in (FIELD_X, FIELD_NORM)
            and item.checkState() == QtCore.Qt.CheckState.Checked
        ):
            self._syncing_fields = True
            try:
                for row in range(self.field_table.rowCount()):
                    other = self.field_table.item(row, column)
                    if other is not None and other is not item:
                        other.setCheckState(QtCore.Qt.CheckState.Unchecked)
            finally:
                self._syncing_fields = False
        self._remember_fields()
        self._update_plot()

    def _update_hkl_label(self, paths: list[Path]):
        """Where the diffractometer was, averaged over each selected scan."""
        lines = []
        missing = 0
        for path in paths:
            values = scan_hkl(path)
            if values is None:
                missing += 1
                continue
            text = format_scan_hkl(values)
            lines.append(
                text if len(paths) == 1 else f"{scan_number(path)}:  {text}"
            )
        if not lines:
            # Silent when nothing records it, rather than a standing complaint.
            self.hkl_label.setText(
                ""
                if not paths
                else "huber_euler_h/k/l not recorded in this scan"
            )
            return
        if missing:
            lines.append(f"({missing} more without huber_euler_h/k/l)")
        self.hkl_label.setText("        ".join(lines))

    def _update_file_label(self, paths: list[Path]):
        """Name the selected file(s) above the plot, by their leading part.

        Only what comes before the scan number is shown, and a series shared by
        several selected files is named once rather than repeated. The full
        names are on the tooltip.
        """
        if not paths:
            self.file_label.setText("No scan selected")
            self.file_label.setToolTip("")
            return
        prefixes: list[str] = []
        for path in paths:
            prefix = file_prefix(path)
            if prefix not in prefixes:
                prefixes.append(prefix)
        self.file_label.setText(",  ".join(prefixes))
        self.file_label.setToolTip("\n".join(str(path) for path in paths))

    def _right_fields(self, y_fields: list[str]) -> set[str]:
        """Of *y_fields*, those drawn against the right-hand scale.

        A stream nominates these so that a small signal is not flattened
        against a large one sharing the scale -- XMCD against its XAS. With
        nothing left on the main scale there is nothing to be flattened
        against, so in that case they keep it.

        Failing a nomination, two ticked fields split one to each scale. Two
        detectors is the case the second scale exists for: they are rarely
        the same size -- an I0 against a fluorescence channel is orders apart
        -- and sharing a scale draws the smaller as a flat line along the
        bottom. Each scale fits its own and the two shapes superimpose, which
        is the comparison that ticking both was asking for. Three or more
        stay as they were: there are only two scales, and picking which two
        of the three get one would be arbitrary.
        """
        stream = self.stream_combo.currentText()
        names = set(stream_default_fields(stream, "right", y_fields))
        if not names and len(y_fields) == 2:
            names = {y_fields[1]}
        return set() if len(names) == len(y_fields) else names

    def _right_scans(self, paths: list[Path]) -> set[str]:
        """Of *paths*, those sent to the right-hand scale by hand.

        The ticked scans that are still selected. Whether that leaves anything
        on the main scale cannot be settled here -- a selected scan need not
        record the ticked field -- so :meth:`_update_plot` makes the final
        call once it knows which curves there actually are.
        """
        return self._right_paths & {str(path) for path in paths}

    def _invert_scans(self, paths: list[Path]) -> set[str]:
        """Of *paths*, those ticked to be drawn upside down."""
        return self._invert_paths & {str(path) for path in paths}

    def _right_axis_free(self) -> bool:
        """Whether the right-hand scale is standing empty and worth filling.

        Read off the last drawing rather than worked out afresh: that is the
        plot being looked at when the menu is reached for, and every rule
        that fills the scale -- the stream's own fields, the scans ticked for
        it, and the back-off when those would empty the main one -- has
        already been through :meth:`_update_plot` by then.

        One scan does not count. With nothing to compare against, that
        back-off hands the curve straight back to the main scale, and the
        tick would sit where it cannot be seen -- "Right scale" is disabled
        on a single scan -- waiting to move it the moment a second is picked.
        """
        return (
            len(self._selected_paths()) > 1
            and bool(self._plotted)
            and not any(curve["on_right"] for curve in self._plotted["curves"])
        )

    def _combine_keys(self, paths: list[Path]) -> set[tuple[str, str, bool]]:
        """The combinations whose two scans are both still selected."""
        live = {str(path) for path in paths}
        return {
            key for key in self._combine if key[0] in live and key[1] in live
        }

    def _populate_right_menu(self):
        """Rebuild the menu of selected scans each time it is opened."""
        paths = self._selected_paths()
        tags = scan_tags(paths)
        self._right_menu.clear()
        for path in paths:
            key = str(path)
            action = self._right_menu.addAction(tags[key])
            action.setCheckable(True)
            # Before connecting, or setting it would re-enter _update_plot
            # once per scan already on the right.
            action.setChecked(key in self._right_paths)
            action.toggled.connect(partial(self._on_right_toggled, key))
        # The curves "Add" has made belong on this menu too, and want it more
        # than the scans do: 171 and 173 read 3.3e5..4.5e5 and their
        # difference -4.9e4..-2.0e4, which shares a scale with neither.
        combos = sorted(self._combine_keys(paths))
        if combos:
            self._right_menu.addSeparator()
        for key in combos:
            first, second, subtract = key
            action = self._right_menu.addAction(
                combine_name(tags[first], tags[second], subtract)
            )
            action.setCheckable(True)
            action.setChecked(key in self._right_combos)
            action.toggled.connect(partial(self._on_right_combo_toggled, key))

    def _on_right_toggled(self, key: str, checked: bool):
        if checked:
            self._right_paths.add(key)
        else:
            self._right_paths.discard(key)
        self._update_plot()

    def _on_right_combo_toggled(self, key: tuple, checked: bool):
        if checked:
            self._right_combos.add(key)
        else:
            self._right_combos.discard(key)
        self._update_plot()

    def _populate_invert_menu(self):
        """Rebuild the menu of selected scans each time it is opened."""
        paths = self._selected_paths()
        tags = scan_tags(paths)
        self._invert_menu.clear()
        for path in paths:
            key = str(path)
            action = self._invert_menu.addAction(tags[key])
            action.setCheckable(True)
            # Before connecting, or setting it would re-enter _update_plot
            # once per scan already inverted.
            action.setChecked(key in self._invert_paths)
            action.toggled.connect(partial(self._on_invert_toggled, key))

    def _on_invert_toggled(self, key: str, checked: bool):
        if checked:
            self._invert_paths.add(key)
            # Negating a scan throws it a scale-height away from the one it
            # is meant to be compared with: 212 and 213 both read LockAC at
            # 2.6e6, so inverting one puts the pair at +2.6e6 and -2.6e6 and
            # a shared scale spans 5.2e6, in which each curve's own 4.2e4
            # swing is 0.8% -- two flat lines. Its own scale is where it can
            # be read, so it goes there unasked while that scale is free.
            # Not when it is busy: it is then already carrying a comparison,
            # and taking it over would break that one to make this one. This
            # lands as an ordinary tick in "Right scale" and is cleared there.
            if self._right_axis_free():
                self._right_paths.add(key)
        else:
            # The tick stays. The curve is still drawn and still on that
            # scale, so unlike a combination below there is something left
            # for it to describe, and silently moving a curve the user never
            # asked to move would be the surprise.
            self._invert_paths.discard(key)
        self._update_plot()

    def _populate_add_menu(self):
        """Rebuild the menu of pairs each time it is opened.

        Three entries to a pair -- the sum, and the subtraction each way
        round, which is how the sign is chosen. Grouped by pair rather than
        all the sums and then all the differences: the question being asked
        is about two scans, and this way the three answers to it are
        together. A scan is never offered against itself.
        """
        paths = self._selected_paths()
        tags = scan_tags(paths)
        self._add_menu.clear()
        for index, first in enumerate(paths):
            for second in paths[index + 1 :]:
                if self._add_menu.actions():
                    self._add_menu.addSeparator()
                for key in (
                    (str(first), str(second), False),
                    (str(first), str(second), True),
                    (str(second), str(first), True),
                ):
                    action = self._add_menu.addAction(
                        combine_name(tags[key[0]], tags[key[1]], key[2])
                    )
                    action.setCheckable(True)
                    action.setChecked(key in self._combine)
                    action.toggled.connect(
                        partial(self._on_combine_toggled, key)
                    )

    def _on_combine_toggled(self, key: tuple, checked: bool):
        if checked:
            self._combine.add(key)
            # The same rule as inverting, and the plainer case for it: a
            # difference is nothing like the size of either scan it came
            # from. 212 - 213 in LockAC spans -25852..5567 where both parents
            # sit at 2.6e6, so on their scale it is a flat line on the axis.
            # A sum is the size of its parents and would be happy beside
            # them, but it is the same menu entry and the same gesture, and
            # one rule that is always true beats two that have to be told
            # apart. Only while that scale is free, and only the first thing
            # to want it: the second keeps the scale the first one chose.
            if self._right_axis_free():
                self._right_combos.add(key)
        else:
            self._combine.discard(key)
            # Nothing stale left pointing at a curve that no longer exists.
            # Asked for again it is a new combination, and the rule above
            # decides its scale afresh against whatever the plot holds then
            # -- which is not always the same answer, and should not be.
            self._right_combos.discard(key)
        self._update_plot()

    def _map_colormap(self, field: str, stream: str) -> str:
        """Which colour map *field* is drawn in as a map on *stream*.

        The box names one for every panel, or leaves the choice here. Left
        here, a field the stream declares ``diverging`` -- a signed channel,
        where zero means something -- gets a scale built around a middle
        colour, and everything else the ordinary dark-to-bright one.

        Deliberately not :meth:`_right_fields`, though both ask the same table
        about the same field. That one gives up when every ticked field is on
        the right, there being nothing left to separate them from; here a
        lone XMCD is exactly as signed as one beside its XAS.
        """
        chosen = self.cmap_combo.currentText()
        if chosen and chosen != AUTO_COLOR_MAP:
            return chosen
        if stream_default_fields(stream, "diverging", [field]):
            return DICHRO_COLOR_MAP
        # The XAS/XMCD pair from a dichro mesh scan is drawn out of the
        # primary stream, which has no stream config to say which of them
        # is signed; the field name itself carries that.
        if field in DIVERGING_FIELDS:
            return DICHRO_COLOR_MAP
        return DEFAULT_COLOR_MAP

    def _update_axis_label(
        self,
        right: set[str],
        split_axis: str = "",
        grid: dict | None = None,
        right_scans: list[str] = (),
        inverted: list[str] = (),
        added: list[str] = (),
        averaged: list[str] = (),
        reduction: tuple[str, str, bool] | None = None,
        reduction_mode: str = "",
        dropped: int = 0,
        log_norm: bool = False,
    ):
        """Say what is being plotted against what, as the plot ended up drawn.

        *right* is the Y fields that went on the right-hand scale, which is
        not simply :meth:`_right_fields` -- a stacked plot gives every field
        its own scale and uses neither. *right_scans* are the scans sent there
        by hand, named rather than folded into the field list: a field sent
        right one scan at a time is on both scales at once, and "(right)"
        after its name would be a flat contradiction of the curve still drawn
        on the left. Pass *grid* when the fields are drawn as colour maps,
        whose axes are the scan's two motors rather than the ticked X field.

        *inverted*, *added* and *averaged* name the scans drawn upside down,
        the curves made out of two of them, and the groups of repeats drawn
        as one. All three say something about the numbers rather than where
        they sit, so unlike "(right)" they are their own clause: a plot
        whose peak is really a dip has to say so, and so does one whose
        single curve is four scans meaned together.

        Pass *reduction* as ``(detector, monitor, transmission)`` where the
        four contrasts are drawn: the ticked X, Y and Norm then describe a
        plot that is not on screen, and the label has to say what is.
        *reduction_mode* names which reduction, where every unit on screen
        agrees on one -- the same four panels come out of two quite different
        measurements, and the file written from here keeps this line. There
        *averaged* and *dropped* report their empty answers as well, for want
        of a status bar to put them in.

        *log_norm* says the Norm tick was applied as a log rather than as a
        division, which the curves and the maps have to spell out for the
        same reason the panels do: the two forms are the same two channels
        and nothing on the plot itself distinguishes a ratio from its log.
        """
        if reduction is not None:
            detector, monitor, transmission = reduction
            ratio = (
                f"log({monitor} / {detector})"
                if transmission
                else f"{detector} / {monitor}"
            )
            # Naming the stream matters here and nowhere else: scan 98's
            # Stream box says dichro_monitor while the reduction reads
            # primary, and this is the only place that shows.
            text = (
                f"X: {ENERGY_FIELD}"
                f"    panels: {', '.join(CONTRAST_LABELS[k] for k in CONTRASTS)}"
                f"    from {PRIMARY_STREAM}: {ratio}"
            )
            if reduction_mode:
                text += f"    {MODE_LABELS[reduction_mode]}"
            # Which units were repeats pooled into one reduction. Said here
            # rather than only in the export because the panels otherwise
            # look the same either way -- one curve per panel is what a
            # single scan gives too -- and this line is what the export
            # copies into its ``view =``. Both answers, unlike the curves
            # case: there the status bar is free to carry them, here it is
            # already spoken for by the per-unit reduction report, and a
            # ticked button that found nothing alike would otherwise be
            # indistinguishable from a button that did nothing.
            if averaged:
                text += f"    averaged: {', '.join(averaged)}"
            elif (
                self.average_scans_btn.isEnabled()
                and self.average_scans_btn.isChecked()
            ):
                text += "    averaged: nothing alike"
            # Both answers again, for the same reason and more sharply: a
            # dropout is one raw point in several hundred, so the panels look
            # much the same whether one went or none did, and only the XMuCD
            # changes enough to notice. Silence would leave the user unable to
            # tell a scan with no dropouts from a button that failed to run.
            if self.spike_btn.isEnabled() and self.spike_btn.isChecked():
                found = self._spikes_removed(dropped, {})
                text += f"    spikes removed: {found}"
            self.axis_label.setText(text)
            return
        monitor = self._monitor_field()
        y_fields = self._checked_y_fields()
        if grid is not None:
            # Two lines rather than one: the axes names and the colour
            # fields each run long on their own -- motor fields are
            # ``huber_hp_nanox_user_setpoint`` and the dichro grid colours
            # two fields at once -- so one line of both and the modifiers
            # after reaches past the window on an ordinary monitor.
            text = (
                f"X: {grid['x_name']}    Y: {grid['y_name']}\n"
                f"colour: {', '.join(y_fields)}"
            )
            if monitor:
                # "colour" and not "Y": on a map that letter is already one
                # of the two motors, so the ratio has to name what is being
                # coloured rather than borrow an axis name for it.
                text += (
                    f"    as log({monitor} / colour)"
                    if log_norm
                    else f"    normalized by {monitor}"
                )
            # A clipped scale is not showing the full range, which the plot
            # should say rather than leave to be noticed from the bar.
            if self.clip_btn.isChecked() and self.clip_btn.isEnabled():
                text += f"    colour clipped at {self.clip_spin.value():g}%"
            # Which map, only when it is not the one the panels would have
            # picked anyway -- naming viridis every time would be noise.
            chosen = self.cmap_combo.currentText()
            if chosen and chosen != AUTO_COLOR_MAP:
                text += f"    in {chosen}"
            # A centred scale reaches past the data at one end, so the bar is
            # wider than the map's range. Say so for the same reason as clip.
            if self.center_btn.isChecked() and self.center_btn.isEnabled():
                text += "    centered on 0"
            self.axis_label.setText(text)
            return
        x_field = self._x_field() or "(pick one)"
        shown = [
            f"{name} (right)" if name in right else name for name in y_fields
        ]
        text = f"X: {x_field}    Y: {', '.join(shown) or '(none)'}"
        if monitor:
            text += (
                f"    as log({monitor} / Y)"
                if log_norm
                else f"    normalized by {monitor}"
            )
        if right_scans:
            text += f"    right scale: {', '.join(right_scans)}"
        if inverted:
            text += f"    inverted: {', '.join(inverted)}"
        if added:
            text += f"    added: {', '.join(added)}"
        if averaged:
            text += f"    averaged: {', '.join(averaged)}"
        if split_axis:
            text += f"    split by {split_axis}"
        self.axis_label.setText(text)

    def _refresh_trend_x_combo(self, paths: list[Path]):
        """Re-offer the baseline fields this selection has in common.

        The chosen field survives a selection that still offers it and falls
        back to the scan number otherwise, that being the one axis every
        selection can draw. Signals stay blocked throughout: the caller is
        mid-redraw and is about to draw the answer anyway.
        """
        names = [TREND_X_SCAN] + baseline_trend_fields(paths)
        current = self.trend_x_combo.currentText()
        self.trend_x_combo.blockSignals(True)
        self.trend_x_combo.clear()
        self.trend_x_combo.addItems(names)
        self.trend_x_combo.setCurrentText(
            current if current in names else TREND_X_SCAN
        )
        self.trend_x_combo.blockSignals(False)

    def _refresh_baseline_combo(self, paths: list[Path]):
        names = [p.name for p in paths]
        current = self.baseline_combo.currentText()
        self.baseline_combo.blockSignals(True)
        self.baseline_combo.clear()
        self.baseline_combo.addItems(names)
        self.baseline_combo.blockSignals(False)
        if current in names:
            self.baseline_combo.setCurrentText(current)

    # -- Image windows ------------------------------------------------------

    def _image_pairs(self, paths: list[Path]) -> list[tuple[Path, Path | None]]:
        return [(path, resolve_image_file(path)) for path in paths]

    def _resolved_images(self, paths: list[Path]) -> list[Path]:
        """The image files for *paths*, in selection order and without repeats.

        Two masters can point at the same image file, and a selection may mix a
        master with its own image file; either way it should appear once.
        """
        images: list[Path] = []
        for _master, image in self._image_pairs(paths):
            if image is not None and image not in images:
                images.append(image)
        return images

    def _image_label(self, image: Path) -> str:
        """How to name *image* in the tooltips -- short, but saying which detector.

        A detector file is named after its scan, so its bare name says nothing
        the master's row does not already say; the subdirectory it sits in is
        the detector. Shown relative to the directory being browsed where it
        is below it, which turns ``scan_000310.h5`` into ``eiger/scan_000310.h5``.
        """
        if self._directory is not None:
            try:
                return str(image.relative_to(self._directory))
            except ValueError:
                pass
        return image.name

    def _update_image_button(self, paths: list[Path]):
        """Enable "Load Image…" only for selections that have image files."""
        found = self._resolved_images(paths)
        self.image_btn.setEnabled(bool(found))
        if not paths:
            self.image_btn.setToolTip(
                "Select a master file to load its detector images."
            )
        elif not found:
            self.image_btn.setToolTip(
                "No detector image file found for the selected file(s)."
            )
        elif len(found) == 1:
            self.image_btn.setToolTip(f"Open {self._image_label(found[0])}")
        else:
            self.image_btn.setToolTip(
                f"Open {len(found)} images side by side, sharing ROIs:\n  "
                + "\n  ".join(self._image_label(image) for image in found)
            )

    def _open_image_viewers(self):
        """Open the selected files' images -- several go side by side in one window."""
        paths = self._selected_paths()
        pairs = self._image_pairs(paths)
        missing = [master.name for master, image in pairs if image is None]
        images = self._resolved_images(paths)
        if not images:
            self.status.showMessage(
                f"No image file found for: {', '.join(missing)}"
            )
            return

        try:
            window = ImageViewerWindow(images, parent=self)
        except (RuntimeError, OSError) as exc:
            QtWidgets.QMessageBox.critical(
                self, "Could not open images", str(exc)
            )
            return

        window.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        window.destroyed.connect(
            lambda _=None, w=window: self._forget_image_window(w)
        )
        # Qt destroys a window that nothing holds a reference to.
        self._image_windows.append(window)
        window.show()

        message = f"Opened {len(images)} image(s) in one window."
        if missing:
            message += f"  No image file found for: {', '.join(missing)}"
        self.status.showMessage(message)

    def _forget_image_window(self, window):
        # Identity, not equality: this runs from destroyed(), by which point the
        # underlying C++ object is gone and comparing wrappers can raise.
        self._image_windows = [
            w for w in self._image_windows if w is not window
        ]

    # -- Spectrum windows ---------------------------------------------------

    def _spectrum_pairs(
        self, paths: list[Path]
    ) -> list[tuple[Path, Path | None]]:
        return [(path, resolve_spectrum_file(path)) for path in paths]

    def _resolved_spectra(self, paths: list[Path]) -> list[Path]:
        """The spectrum files for *paths*, in selection order and without repeats."""
        spectra: list[Path] = []
        for _master, spectrum in self._spectrum_pairs(paths):
            if spectrum is not None and spectrum not in spectra:
                spectra.append(spectrum)
        return spectra

    def _spectrum_label(self, spectrum: Path) -> str:
        """How to name *spectrum* in the tooltips -- as ``vortex/scan_000055.h5``."""
        if self._directory is not None:
            try:
                return str(spectrum.relative_to(self._directory))
            except ValueError:
                pass
        return spectrum.name

    def _update_spectrum_button(self, paths: list[Path]):
        """Enable "Load Spectrum…" only for selections that have spectrum files."""
        found = self._resolved_spectra(paths)
        self.spectrum_btn.setEnabled(bool(found))
        if not paths:
            self.spectrum_btn.setToolTip(
                "Select a master file to load its fluorescence spectra."
            )
        elif not found:
            self.spectrum_btn.setToolTip(
                "No spectrum file found for the selected file(s)."
            )
        elif len(found) == 1:
            self.spectrum_btn.setToolTip(
                f"Open {self._spectrum_label(found[0])}"
            )
        else:
            self.spectrum_btn.setToolTip(
                f"Overlay {len(found)} spectra in one window:\n  "
                + "\n  ".join(self._spectrum_label(s) for s in found)
            )

    def _open_spectrum_viewers(self):
        """Open the selected scans' spectra -- several overlay in one window.

        Overlaid rather than side by side, unlike the image windows: two
        spectra on one pair of axes is the comparison, where two images on one
        would be a mess.
        """
        paths = self._selected_paths()
        pairs = self._spectrum_pairs(paths)
        missing = [
            master.name for master, spectrum in pairs if spectrum is None
        ]
        spectra = self._resolved_spectra(paths)
        if not spectra:
            self.status.showMessage(
                f"No spectrum file found for: {', '.join(missing)}"
            )
            return

        # Summarising the cubes takes a couple of seconds apiece on the long
        # scans, and the window cannot be shown until it is done.
        QtWidgets.QApplication.setOverrideCursor(
            QtCore.Qt.CursorShape.WaitCursor
        )
        self.status.showMessage(f"Reading {len(spectra)} spectrum file(s)…")
        QtWidgets.QApplication.processEvents()
        try:
            window = SpectrumViewerWindow(spectra, parent=self)
        except (RuntimeError, OSError) as exc:
            QtWidgets.QMessageBox.critical(
                self, "Could not open spectra", str(exc)
            )
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

        window.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        window.destroyed.connect(
            lambda _=None, w=window: self._forget_spectrum_window(w)
        )
        # Qt destroys a window that nothing holds a reference to.
        self._spectrum_windows.append(window)
        window.show()

        message = f"Opened {len(spectra)} spectrum(s) in one window."
        if missing:
            message += f"  No spectrum file found for: {', '.join(missing)}"
        self.status.showMessage(message)

    def _forget_spectrum_window(self, window):
        # Identity, not equality: see _forget_image_window.
        self._spectrum_windows = [
            w for w in self._spectrum_windows if w is not window
        ]

    # -- Plotting -----------------------------------------------------------

    def _set_panels(self, count: int, across: bool = False) -> bool:
        """Give the figure *count* plots, stacked or side by side.

        Curves stack vertically on a shared X axis, one Y field above the
        next. Colour maps go side by side instead, sharing both axes: what two
        of them have in common is the area they were measured over, and the
        point of drawing them together is to compare the same spot.

        Each plot gets two more axes sharing its X, for data needing its own Y
        scale: one for the Y fields a stream puts on the right (XMCD against
        its XAS), one for the dY/dX overlay. They are independent so both can
        be shown at once, in which case the derivative's spine is pushed clear
        of the other's; otherwise whichever is visible sits at the right edge
        on its own.

        Rebuilding throws the axes away, so it happens only when the count or
        the direction actually changes. Returns whether it did.
        """
        count = max(1, count)
        # A single plot is the same plot either way, so do not throw it away
        # -- and the zoom on it -- just because the direction nominally flips.
        across = across and count > 1
        if len(self._panels) == count and across == self._panels_across:
            return False
        # figure.clear() takes the colour bars' own axes with everything else.
        self._colorbars = []
        self.figure.clear()
        self._panels_across = across
        mains = (
            self.figure.subplots(
                1, count, sharex=True, sharey=True, squeeze=False
            )[0]
            if across
            else self.figure.subplots(count, 1, sharex=True, squeeze=False)[
                :, 0
            ]
        )
        self._panels = [(main, main.twinx(), main.twinx()) for main in mains]
        self.ax, self.ax2, self.ax3 = self._panels[0]
        self._twins = tuple(
            twin for _main, *twins in self._panels for twin in twins
        )
        # Kept so each grid's own readout can be swapped in and out.
        self._default_format_coord = [
            main.format_coord for main, *_ in self._panels
        ]
        # The view history still points at the axes just discarded.
        self.toolbar.update()
        return True

    def _apply_log_scales(self, autoscale: bool = False):
        """Put the Log X / Log Y states on the main figure's axes.

        Every draw clears the axes, and clear() puts a linear scale back, so
        this runs on each draw as well as on a press. Non-positive points are
        masked rather than clipped: clipped, a zero is drawn at the bottom of
        the frame and pulls a vertical stroke down to it from its neighbours.
        *autoscale* refits the view, which a press wants -- limits sized for
        a linear axis are no use on a log one, and the reverse -- and a draw
        does not, since it sets its own limits afterwards.

        An axis with nothing positive on it stays linear: a log locator with
        no positive data raises on the next draw. That covers the empty and
        hidden twins, and a channel such as XMCD that is negative throughout.
        """

        def scale(wanted: bool, minpos: float) -> dict:
            if wanted and np.isfinite(minpos):
                return {"value": "log", "nonpositive": "mask"}
            return {"value": "linear"}

        log_x = self.log_x_btn.isChecked()
        log_y = self.log_y_btn.isChecked()
        for main, *twins in self._panels:
            if autoscale:
                for axes in (main, *twins):
                    # A toolbar zoom turns autoscale off, and autoscale_view()
                    # then does nothing: the new log axis would keep a linear
                    # zoom's limits, below zero included. A map is the
                    # exception -- relim() skips its pcolormesh, so a refit
                    # would size it to nothing; its limits stay as drawn.
                    if not main.collections:
                        axes.set_autoscale_on(True)
                    axes.relim()
            # The twins share main's X, so its scale reaches them too.
            main.set_xscale(**scale(log_x, main.dataLim.minposx))
            for axes in (main, *twins):
                axes.set_yscale(
                    **scale(log_y and axes.get_visible(), axes.dataLim.minposy)
                )
                if autoscale:
                    axes.autoscale_view()
        if autoscale:
            for twin in self._twins:
                twin.set_xlim(self.ax.get_xlim(), emit=False)

    def _toggle_log_scales(self):
        self._apply_log_scales(autoscale=True)
        # The view history holds limits sized for the old scale; Home would
        # go back to them.
        self.toolbar.update()
        self.canvas.draw_idle()

    def _draw_trend(
        self,
        curves: list[dict],
        trending: bool,
        quantity: str,
        show_derivative: bool,
        notes: list[str],
        peaks: int = 0,
        against: str = TREND_X_SCAN,
    ) -> list[tuple[str, list[tuple[float, float | None]]]]:
        """Trend one fitted quantity across the selection, beside the plot.

        Returns what it drew, as ``(y field, [(x, value), ...])`` per
        series, for the export header to report. Draws nothing and returns
        nothing when the trend is off.

        *against* is the X axis and *peaks* how many components "Multi"
        fitted, or 0 when it is off; :func:`trend_points` places the points
        and names the series, and says what both mean.

        Non-zero *peaks* gives one line per component, taken in the order
        ``peak_fit_multi`` returns them -- which is by centre, and that
        ordering is the only reason a line here means anything. By height
        the series would swap members at the scan where two peaks cross,
        and draw a step that is in the labelling and not in the sample.
        """
        was_visible = self.trend_widget.isVisible()
        self.trend_widget.setVisible(trending)
        if not trending:
            return []
        # Only on the way in, so that dragging the handle to favour one plot
        # survives every later redraw. Qt gives a freshly shown child a width
        # of its own choosing otherwise, which is rarely half.
        if not was_visible:
            width = self.plot_split.width()
            if width > 0:
                self.trend_widget.setMinimumWidth(0)
                self.plot_split.setSizes([width // 2, width - width // 2])

        key = TREND_QUANTITIES[quantity]
        # Where every point goes and what it fitted, read out here down to
        # the one quantity asked for. The same function the export button
        # at the end of this row calls, so the picture and the file are the
        # same numbers laid out two ways -- see :func:`trend_points`.
        series, skipped = trend_points(curves, against, peaks)

        self.trend_ax.clear()
        for field in sorted(series):
            points = series[field]
            # A declined fit becomes a gap rather than a straight line drawn
            # across it: nan breaks the line where the number is missing,
            # which is the honest picture of a scan that would not fit.
            self.trend_ax.plot(
                [position for position, _stats in points],
                _fit_column(points, key),
                marker="o",
                markersize=4,
                linewidth=1,
                label=field,
            )
        self.trend_ax.set_xlabel(against)
        # Says which data the numbers came off, since with dY/dX shown these
        # are the derivative's fit and nothing on this axis would otherwise
        # distinguish a beam size from a peak width.
        self.trend_ax.set_ylabel(
            f"{quantity} of dY/dX" if show_derivative else quantity
        )
        # Scan numbers are integers, and the default locator is happy to
        # label a run of three scans at 101.5. A baseline reading is not an
        # integer and must not be rounded to one: eta's 0 to 90 would
        # survive it, but the knife edge's 0.0198 to 0.0217 would come back
        # with no labels at all. clear() above resets the locator, so there
        # is nothing to put back on the other branch.
        if against == TREND_X_SCAN:
            self.trend_ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        self.trend_ax.grid(True, alpha=0.3)
        if len(series) > 1:
            self.trend_ax.legend(fontsize=8)
        self.trend_figure.tight_layout()
        self.trend_canvas.draw_idle()

        # A trend that drew nothing looks exactly like a button that did
        # nothing, the trap the despike and averaging notes are there to
        # avoid, so both ways of coming up empty are said out loud.
        if skipped:
            one = skipped == 1
            missing = (
                "no scan number"
                if against == TREND_X_SCAN
                else f"no {against} reading"
            )
            notes.append(
                f"Trend: {skipped} curve{'' if one else 's'} with"
                f" {missing} of {'its' if one else 'their'} own left out."
            )
        if not series:
            # Unless the line above has already said why, in which case
            # repeating it adds a second sentence and no second fact.
            if not skipped:
                notes.append("Trend: nothing with a scan number to plot.")
        elif not any(
            stats is not None and stats.get(key) is not None
            for points in series.values()
            for _position, stats in points
        ):
            notes.append(f"Trend: no curve has a {quantity} to plot.")
        # The one quantity again, not the fits behind it: this is what the
        # header reports, and it reports what was drawn.
        return [
            (
                field,
                [
                    (position, None if stats is None else stats.get(key))
                    for position, stats in series[field]
                ],
            )
            for field in sorted(series)
        ]

    def _split_steps(
        self, path: Path, stream: str, data: dict
    ) -> tuple[int, list, list[str]]:
        """Whether and how this scan's stream splits into several curves.

        Returns ``(points per mesh row, groups, the axes split on)``, and
        ``(1, [], [])`` for a stream that names no axis to split on or a scan
        that is not a mesh -- which is every ordinary scan.
        """
        names, values, setpoints = split_axes(path, stream)
        columns = [data[field] for field in values if field in data]
        if not columns or len(columns) != len(values):
            return 1, [], []
        if not all(np.any(np.isfinite(column)) for column in columns):
            return 1, [], []
        steps = mesh_step_count(path, columns[0].size)
        if steps < 2:
            return 1, [], []
        # Where the axes came from the scan's declared raster rather than from
        # the stream's own spec, the data has to bear that raster out: a dichro
        # scan's primary holds sub-measurements its shape does not describe.
        demanded = [data[field] for field in setpoints if field in data]
        if setpoints and (
            len(demanded) != len(setpoints)
            or not confirms_period(demanded, steps)
        ):
            return 1, [], []
        groups = step_groups(
            columns, steps, SPLIT_ATOL, self.average_btn.isChecked()
        )
        return steps, groups, names

    @staticmethod
    def _group_tags(groups: list, axes: list[str]) -> list[str]:
        """A short name per split group: each axis, and the value it sat at.

        Readbacks differ from step to step, so the values alone normally tell
        the steps apart. Where they do not -- a simulated positioner that never
        moved, or two steps of a cycle that return to the same setting -- the
        step number is added so that no two curves share a label.
        """
        tags = [
            ", ".join(
                f"{axis} {value:+.4g}" for axis, value in zip(axes, values)
            )
            for values, _members in groups
        ]
        if len(set(tags)) != len(tags):
            tags = [
                f"{tag} #{members[0]}"
                for tag, (_values, members) in zip(tags, groups)
            ]
        return tags

    def _spike_monitor(self, path: Path, monitor_field: str) -> str:
        """Which field says whether the beam was there, for one scan's curves.

        Whatever is normalising the curves, where anything is, so that the
        filter and the division agree about which channel the beam is. Where
        nothing is -- the norm box is left unticked by default -- the scan's
        own hinted monitor, which is the channel it counted against and the
        one the reduction would reach for too.

        The hint is only taken when the stream in fact carries it. Scans 45,
        89, 90 and 91 hint ``Time``, which is not a field of theirs, and taken
        on trust that would light the button on four of the very scans it was
        written for and then quietly do nothing to them. Where the hint is no
        use :data:`REDUCE_MONITORS` is walked instead -- the same list in the
        same order as :func:`reduce_defaults`, and for the same reason.

        Asked of the names the field table was built from rather than of the
        file, so that it costs no read: :func:`read_scan_info` is cached but
        :func:`stream_fields` is not, and this is asked per path per redraw.
        That list is the union over the selection, so a name it offers that
        this one file lacks leaves the scan unjudged -- which is what no
        monitor at all would have done anyway.
        """
        if monitor_field:
            return monitor_field
        hint = read_scan_info(path)["hints"]["monitor"]
        if hint in self._fields:
            return hint
        return next(
            (name for name in REDUCE_MONITORS if name in self._fields), ""
        )

    def _read_curves(
        self,
        paths: list[Path],
        stream: str,
        x_field: str,
        y_fields: list[str],
        monitor_field: str,
        transmission: bool = False,
        drop_spikes: bool = False,
    ) -> tuple[list[dict], str, dict[str, int]]:
        """Every curve the selection plots, read and reduced but not yet drawn.

        One curve per (file, Y field) as a rule, and -- where the stream splits
        along a mesh's inner axis -- one per step or distinct value of it on
        top of that. Colours are handed out per file and split group rather
        than per curve, so one magnet current keeps its colour across the Y
        fields it is plotted in. Returns the curves, the axis split on, if
        any -- which is what decides whether the plot is stacked -- and how
        many points the spike filter took from each scan.

        That last is keyed by path and holds an entry only for the scans
        actually judged, so that a scan whose beam never faltered reads as the
        zero it is rather than as one there was no monitor to ask about.

        *transmission* picks which way the monitor is applied: the plain
        ``y/monitor`` a counting channel wants, or ``log(monitor/y)`` for a
        measurement made through the sample. :func:`absorption` does both,
        the same call the XMCD panels reduce with, so a transmission curve
        here and a panel above it cannot come out in different units.

        With *drop_spikes* the monitor is read whether or not it is being
        divided by, because here it is being asked a different question: not
        how bright the beam was, but whether there was any. A point it says
        had none becomes a gap rather than a value, before the split groups
        are averaged, since a dead point meaned into a repeat has already
        done its damage.
        """
        normalize = bool(monitor_field)
        overlay = len(paths) > 1
        curves: list[dict] = []
        split_axis = ""
        removed: dict[str, int] = {}
        colors: dict[tuple, int] = {}
        for path in paths:
            fields_needed = {x_field, *y_fields}
            if normalize:
                fields_needed.add(monitor_field)
            judge = (
                self._spike_monitor(path, monitor_field) if drop_spikes else ""
            )
            if judge:
                fields_needed.add(judge)
            try:
                # Which columns carry the split is a question only this scan
                # answers -- where the stream names no split field of its own,
                # the axes are the scan's -- and asking reads the file, so it
                # belongs where an unreadable one is already handled.
                _names, values, setpoints = split_axes(path, stream)
                fields_needed.update(values, setpoints)
                data = read_stream(path, list(fields_needed), stream)
            except OSError as exc:
                self.status.showMessage(f"Could not read {path.name}: {exc}")
                continue
            if x_field not in data:
                continue
            x = data[x_field]
            monitor = data.get(monitor_field) if normalize else None
            if normalize and monitor is None:
                self.status.showMessage(
                    f"{path.name}: monitor field {monitor_field!r} not found -- plotting raw values."
                )
            # Judged on the monitor as read, not on the normalised curve: the
            # question is whether the beam was there, and the ratio of two
            # dead channels answers a different one. A scan whose monitor the
            # stream does not carry simply keeps all its points.
            gone = None
            if judge and data.get(judge) is not None:
                gone = lone_dropouts(data[judge])
                removed[str(path)] = int(gone.sum())
                if not gone.any():
                    gone = None
            steps, groups, axes = self._split_steps(path, stream, data)
            split_axis = split_axis or ", ".join(axes)
            tags = self._group_tags(groups, axes)
            for y_field in y_fields:
                if y_field not in data:
                    continue
                y = data[y_field]
                if monitor is not None:
                    y = absorption(y, monitor, transmission)
                if gone is not None:
                    # Ragged columns: a scan stopped part way leaves its last
                    # point written into some fields and not others, so the
                    # monitor and this field need not be the same length. They
                    # start at the same row, so the shorter says how far the
                    # judgement reaches and the rest of the curve stands.
                    reach = min(gone.size, y.size)
                    here = np.zeros(y.size, dtype=bool)
                    here[:reach] = gone[:reach]
                    # np.where rather than an assignment: unnormalised, y is
                    # still the array read_stream handed back, and X is shared
                    # across this scan's curves -- the same aliasing the
                    # inversion below is careful about.
                    y = np.where(here, np.nan, y)
                pieces = [("", x, y)]
                if groups:
                    pieces = [
                        (tag, *group_curve(x, y, steps, members))
                        for tag, (_values, members) in zip(tags, groups)
                    ]
                base = f"{path.name}: {y_field}" if overlay else y_field
                for tag, piece_x, piece_y in pieces:
                    # Stacked, each plot's own Y label names the field, so the
                    # legend there need only say which group a curve is.
                    short = (
                        (f"{path.name}: {tag}" if overlay else tag)
                        if tag
                        else base
                    )
                    key = (path, tag) if tag else (path, y_field)
                    curves.append(
                        {
                            "label": f"{base} @ {tag}" if tag else base,
                            "short": short,
                            "y_field": y_field,
                            # Which scan it came from, for the right-hand
                            # scale: that choice is per scan, not per field.
                            "path": str(path),
                            "x": piece_x,
                            "y": piece_y,
                            "color": f"C{colors.setdefault(key, len(colors)) % 10}",
                            "tag": tag,
                        }
                    )
        return curves, split_axis, removed

    def _averaged_curves(
        self, curves: list[dict], paths: list[Path], tags: dict[str, str]
    ) -> tuple[list[dict], list[str]]:
        """The repeats among *curves* gathered into one curve each.

        A group is scans that :func:`scan_groups` already calls one
        measurement, drawing the same field and the same split group, whose
        X columns agree point for point. Its members leave the plot and
        their mean takes their place -- in the first one's position and its
        colour, so that ticking the button does not reshuffle the plot or
        recolour what is left of it.

        Averaged row by row rather than by interpolation: repeats step the
        same points, to within 0.07 of a step everywhere in this data, and a
        scan stopped part way is then simply short. :func:`stack_mean`
        averages the rows it reached and leaves the rest blank, which is
        what the repeats inside a single scan already do.

        Returns the curves to draw and the names of the groups made, for the
        axis label and the export header to say what was averaged.
        """
        where = {
            str(path): index
            for index, group in enumerate(scan_groups(paths))
            for path in group
        }
        buckets: dict[tuple, list[int]] = {}
        for index, curve in enumerate(curves):
            key = (
                where.get(curve["path"], curve["path"]),
                curve["y_field"],
                curve["tag"],
            )
            buckets.setdefault(key, []).append(index)
        made: dict[int, dict] = {}
        gone: set[int] = set()
        names: list[str] = []
        for members in buckets.values():
            # One measurement and one field, and the scans can still have
            # stepped different points: 4 and 6 share a plan, a motor, a
            # range and 50 points, and visit them in a different order. Only
            # the columns say so, and they are read by now. Gathered
            # first-fit against each cluster's first scan, as step_groups
            # gathers the repeats within one.
            clusters: list[list[int]] = []
            for index in members:
                for cluster in clusters:
                    if same_grid(curves[cluster[0]]["x"], curves[index]["x"]):
                        cluster.append(index)
                        break
                else:  # nothing on these points yet, so this starts a cluster
                    clusters.append([index])
            for cluster in clusters:
                if len(cluster) == 1:
                    continue
                group = [curves[index] for index in cluster]
                first = group[0]
                rows = max(np.size(curve["x"]) for curve in group)
                name = _unit_tag(tags, [Path(c["path"]) for c in group])
                base = f"{name}: {first['y_field']}"
                label = f"{base} @ {first['tag']}" if first["tag"] else base
                made[cluster[0]] = {
                    "label": label,
                    "short": f"{name}: {first['tag']}"
                    if first["tag"]
                    else base,
                    "y_field": first["y_field"],
                    # The first member's, so that anything looking a path up
                    # in tags still finds one. Never used to name it.
                    "path": first["path"],
                    "x": stack_mean([c["x"] for c in group], rows),
                    "y": stack_mean([c["y"] for c in group], rows),
                    "color": first["color"],
                    "tag": first["tag"],
                    "panel": first["panel"],
                    # Already settled for every member, so taken from them
                    # rather than worked out again: a repeat sent to the
                    # right-hand scale takes the curve replacing it there.
                    "on_right": any(c["on_right"] for c in group),
                    "combo": label,
                }
                gone.update(cluster[1:])
                names.append(name)
        return (
            [
                made.get(index, curve)
                for index, curve in enumerate(curves)
                if index not in gone
            ],
            names,
        )

    def _combined_curves(
        self,
        curves: list[dict],
        paths: list[Path],
        tags: dict[str, str],
        y_fields: list[str],
        stacked: bool,
    ) -> list[dict]:
        """The sums and differences asked for, as curves of their own.

        Paired by field and by split group, so two scans with three fields
        ticked give three sums, one per field. ``tag`` is "" on an unsplit
        scan and the ordinary case matches at once; where the scans were
        split, only groups that agree on their setpoints are paired, which
        is the only pairing that means anything.

        Built here rather than in the reader because the inversion has to
        happen first -- ticking 173 inverted and asking for 171 + 173 then
        gives exactly 171 - 173, which it would not the other way round.
        """
        by_key = {}
        for curve in curves:
            by_key[(curve["path"], curve["y_field"], curve["tag"])] = curve
        # Carry on from the colours already drawn rather than restarting at
        # C0, or a sum would come out the same colour as one of its parents.
        used = len({curve["color"] for curve in curves})
        combined: list[dict] = []
        missed: list[str] = []
        for key in sorted(self._combine_keys(paths)):
            first, second, subtract = key
            name = combine_name(tags[first], tags[second], subtract)
            pairs = [
                (by_key[(first, field, tag)], by_key[(second, field, tag)])
                for (path, field, tag) in by_key
                if path == first and (second, field, tag) in by_key
            ]
            if not pairs:
                missed.append(name)
                continue
            for curve_a, curve_b in pairs:
                x, y = combine_curves(
                    curve_a["x"],
                    curve_a["y"],
                    curve_b["x"],
                    curve_b["y"],
                    subtract,
                )
                # Named the way _read_curves names an overlaid curve, with
                # the two short scan tags standing in for the file name: two
                # 26-character names in one legend entry would be half the
                # box, and the short form is what the user clicked anyway.
                tag = curve_a["tag"]
                base = f"{name}: {curve_a['y_field']}"
                combined.append(
                    {
                        "label": f"{base} @ {tag}" if tag else base,
                        "short": f"{name}: {tag}" if tag else base,
                        "y_field": curve_a["y_field"],
                        # The first scan's, so that anything looking a path
                        # up in tags still finds one. Never used to name it.
                        "path": first,
                        "x": x,
                        "y": y,
                        "color": f"C{used % 10}",
                        "tag": tag,
                        # Past the loop that settles these for the curves
                        # that were read, so they are settled here instead.
                        "panel": (
                            y_fields.index(curve_a["y_field"]) if stacked else 0
                        ),
                        # Nothing is inherited: a combination starts on the
                        # main scale and is sent right only by being ticked
                        # for it. Neither rule the parents follow carries any
                        # meaning here. A scale chosen for two LockDC
                        # readings says nothing about their difference, which
                        # on 171 and 173 is -4.9e4 where they are 3.3e5..4.5e5
                        # -- a scale that fits the parents draws it flat
                        # along the bottom. And which scan the user sent
                        # right is a choice about two particular curves; a
                        # sum is neither of them. So it gets its own entry in
                        # that menu, which is where the choice belongs: a
                        # parent on the right with its sum on the left is a
                        # reading the two scales are there to allow.
                        "on_right": key in self._right_combos,
                        # What the export calls it: its path names only the
                        # first scan, whose own column has that name already.
                        "combo": f"{base} @ {tag}" if tag else base,
                    }
                )
                used += 1
        # Said once rather than left to be noticed from a curve that never
        # appeared, and the tick stays where it was put. Two ways to get
        # here: a field ticked on one scan and not recorded by the other,
        # and -- scans 96 and 98 -- two scans that split differently, where
        # pairing anything would be pairing two different field settings.
        if missed:
            self.status.showMessage(
                f"{', '.join(missed)}  -  the two draw no curve in common,"
                " nothing to add."
            )
        return combined

    def _read_contrasts(
        self,
        paths: list[Path],
        detector: str,
        monitor: str,
        transmission: bool,
        average: bool = False,
        drop_spikes: bool = False,
    ) -> tuple[list[dict], list[tuple], list[str], int]:
        """Reduce each unit into one curve per contrast, already assigned to panels.

        A unit is what one reduction runs on -- a dichro scan on its own, or
        a group of lock-in scans whose field was stepped between them, as
        :func:`contrast_units` decides. Everything below is per unit rather
        than per scan for that reason: 235 and 236 are one measurement and
        draw one set of curves. With *average* the repeats among those units
        are joined into one apiece as well (:func:`averaged_units`), which
        is the same reduction over more points rather than a second kind of
        curve -- so nothing below this has to know it happened.

        Curve dicts of exactly the shape :meth:`_read_curves` returns, so the
        drawing below takes them unchanged -- which is the point. The legend,
        the zoom restore, the statistics table and the dY/dX overlay all come
        for free, and there is no second copy of that code to drift out of
        step with this one.

        The colour is per unit rather than per curve, so one measurement keeps
        its colour down all four panels and they can be read as a column.

        A unit the reduction cannot make anything of contributes no curves and
        is reported instead of aborting the lot: the reason for overlaying two
        runs is to compare them, and dropping the pair because one has a dead
        I0 would lose the one that does not. The reports are
        ``(tag, result, missing)``, which is what the status bar says.

        With *drop_spikes* the points where the beam was gone are left out of
        every unit, and how many went altogether comes back as the fourth
        value -- a total rather than a tally per unit, because the panels have
        no status bar to spend on it and the axis label says it once.
        """
        curves: list[dict] = []
        reports: list[tuple] = []
        colors: dict[str, int] = {}
        units = contrast_units(paths)
        averaged: list[str] = []
        dropped = 0
        if average:
            units, averaged = averaged_units(units, paths)
        for tag, unit, mode in units:
            try:
                columns, missing, removed = read_contrast_columns(
                    unit, mode, detector, monitor, drop_spikes
                )
            except OSError as exc:
                reports.append((tag, None, [str(exc)]))
                continue
            if missing:
                reports.append((tag, None, missing))
                continue
            dropped += removed
            if mode == MODE_LOCKIN_POL:
                result = reduce_lockin(
                    columns["energy"],
                    columns["step"],
                    columns["detector"],
                    columns["ac"],
                    columns["monitor"],
                    transmission=transmission,
                )
                # Which conditions the group was joined on, carried alongside
                # rather than inferred downstream: the join is the one thing
                # here that can be wrong while looking entirely right, and
                # the status bar has to be able to name it. A copy, so the
                # reducer's own result stays what it returned.
                result = dict(result, conditions=self._unit_conditions(unit))
            else:
                result = reduce_dichro(
                    columns["energy"],
                    columns["pr"],
                    columns["step"],
                    columns["detector"],
                    columns["monitor"],
                    transmission=transmission,
                )
            reports.append((tag, result, []))
            if not result["n_complete"]:
                continue
            color = f"C{colors.setdefault(tag, len(colors)) % 10}"
            for panel, key in enumerate(CONTRASTS):
                label = f"{tag}: {CONTRAST_LABELS[key]}"
                curves.append(
                    {
                        "label": label,
                        # Each panel's Y label names the contrast, so the one
                        # legend on top need only say which scan.
                        "short": tag,
                        "y_field": key,
                        # The unit's first scan. Nothing keyed on this is
                        # live while reducing -- the right-hand scale, the
                        # inversions and the sums are all switched off there
                        # -- and it keeps the curves the shape the drawing
                        # expects.
                        "path": str(unit[0]),
                        "x": result["energy"],
                        "y": result[key],
                        "color": color,
                        "tag": tag,
                        "panel": panel,
                        "on_right": False,
                    }
                )
        return curves, reports, averaged, dropped

    @staticmethod
    def _unit_conditions(paths: list[Path]) -> str:
        """The temperature a lock-in group was joined at, for the status bar.

        Joining scans measured at two temperatures gives a clean-looking
        curve that is half one and half the other, so the figure the grouping
        went on is worth showing rather than leaving to be trusted.
        """
        readings = [
            value
            for path in paths
            if (
                value := baseline_values(path, (TEMPERATURE_FIELD,)).get(
                    TEMPERATURE_FIELD
                )
            )
            is not None
        ]
        return f"{np.mean(readings):.0f} K" if readings else ""

    @staticmethod
    def _spikes_removed(dropped: int, removed: dict[str, int]) -> str:
        """What the filter took, worded once for label, bar and file alike.

        In one place because the three otherwise drift apart, and visibly: the
        axis label is copied verbatim into the export's ``view =`` line, which
        sits directly above the ``spikes =`` line this also writes.

        *removed* is the per-scan count the curves keep and names only the
        scans that lost something; *dropped* the single total the panels have,
        their files being pooled before the cut. Either way nothing found says
        so in words rather than as a nought, since "0" beside a button reads
        as a button that did not run.
        """
        named = "   ".join(
            f"{tag}: {count} point(s)"
            for tag, count in removed.items()
            if count
        )
        if named:
            return named
        return f"{dropped} point(s)" if dropped else "none found"

    @staticmethod
    def _contrast_message(
        reports: list[tuple], detector: str, monitor: str
    ) -> str:
        """What the reduction made of each scan, for the status bar.

        One segment per unit, because the panels overlay several and a count
        that does not say whose it is answers nothing. The failures are
        spelled out rather than left as a zero: a dead channel and a scan
        still on its first field setting both show as four empty panels, and
        which of the two it was decides whether the fix is the Monitor box or
        selecting the scan at the other field.
        """
        parts = []
        for tag, result, missing in reports:
            if result is None:
                parts.append(f"{tag}: no {', '.join(missing)}")
                continue
            text = (
                f"{tag}: {result['n_complete']}/{result['n_energies']} energies"
                f"  -  {result['n_points']} points"
            )
            # Which conditions a lock-in group was pooled on; empty for a
            # dichro scan, which is one file and pools nothing.
            if result.get("conditions"):
                text += f"  -  {result['conditions']}"
            if result["n_complete"]:
                # XNCD is even in the field, so the sample cannot make it:
                # what shows there is the helicity average failing to cancel
                # XMCD, and the same leakage is in one other panel. The two
                # together, not either alone, say whether that panel is
                # signal -- see :func:`xncd_partner` for which, and why it is
                # not always the same one.
                partner = xncd_partner(result)
                text += (
                    f"  -  XNCD rms {_rms(result['xncd']):.3g}"
                    f" vs {CONTRAST_LABELS[partner]} rms"
                    f" {_rms(result[partner]):.3g}"
                )
            elif result["n_points"] and not result["n_finite"]:
                text += "  -  no finite absorption, try another channel"
            elif result["n_steps"] < 2:
                # A dichro scan steps the field itself and has not got there
                # yet; a lock-in scan usually never will, the other setting
                # being a second file that has to be selected alongside.
                text += (
                    "  -  only one field setting, select the scan at the"
                    " opposite field"
                    if result["mode"] == MODE_LOCKIN_POL
                    else "  -  only one field setting so far"
                )
            parts.append(text)
        return f"{detector}/{monitor}  -  " + "   ".join(parts)

    def _export(self):
        """Write what is on screen to a file the user names."""
        # Snapshot before the dialog, as xmucd_live.py's export does and for
        # the same reason: the dialog runs a nested event loop, so a Refresh
        # behind it can redraw, and the file would otherwise not match the
        # numbers in its own header.
        plotted = self._plotted
        if plotted is None:
            self.status.showMessage("Nothing to export -- nothing is plotted.")
            return

        folder = self._export_folder()
        path, _filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export what is plotted",
            str(folder / self._export_name(plotted)),
            "Text data (*.csv *.dat *.txt);;Picture (*.png *.pdf *.svg);;"
            "All files (*)",
        )
        if not path:  # cancelled
            return
        if not Path(path).suffix:
            path += ".csv"

        # Asked of the canvas rather than matched against a list written here:
        # it also renders eps, jpg, tif and several more, and a name this does
        # not recognise falls through to the text branch -- so a .jpg left out
        # of the list would quietly be written as CSV.
        suffix = Path(path).suffix.lstrip(".").lower()
        picture = suffix in self.canvas.get_supported_filetypes()

        # Refused before the name is settled below, so that a map asked for
        # as text never gets as far as offering to replace anything.
        if not picture and plotted["mode"] == "grid":
            self.status.showMessage(
                "A colour map saves as a picture only -- try a .png name."
            )
            return
        target = self._export_target(Path(path), plotted, picture)
        if target is None:  # the rename would have replaced a file
            return

        if picture:
            try:
                self.figure.savefig(target)
            except OSError as exc:
                self.status.showMessage(f"Export failed: {exc.strerror or exc}")
                return
            # Where the toolbar's own Save button remembers its directory. The
            # two save the same picture, so they share the one memory rather
            # than wander apart.
            rcParams["savefig.directory"] = str(target.parent)
            self.status.showMessage(f"Saved the plot to {target}")
            return

        curves = plotted["curves"]
        if not curves:
            self.status.showMessage("Nothing to export -- no curve is plotted.")
            return

        tags = scan_tags(plotted["paths"])
        names, columns = export_table(
            curves,
            tags,
            len(plotted["paths"]) > 1,
            plotted["mode"] == "xmcd",
        )
        try:
            write_columns(target, names, columns, self._export_header(plotted))
        except OSError as exc:
            self.status.showMessage(f"Export failed: {exc.strerror or exc}")
            return
        rows = max((column.size for column in columns), default=0)
        self.status.showMessage(
            f"Exported {len(curves)} curves, {rows} rows, to {target}"
        )

    def _export_fits(self):
        """Write the fitted values of every scan against the trend's axis.

        The trend's picture as a file, and all of it: the plot draws the
        one quantity its combo names, where a peak is read by its four
        numbers together, and exporting under each in turn would make four
        files that then have to be joined back up. "Integrate" adds the
        curve's own trapezoids on the same terms, since they are read
        beside the fit and the table already shows the two together.
        """
        # Snapshotted before the dialog, for the reason _export gives.
        plotted = self._plotted
        if plotted is None or not plotted["trended"]:
            self.status.showMessage(
                "Nothing to export -- the trend is not shown."
            )
            return
        against = plotted["trend_against"]
        series, skipped = trend_points(
            plotted["curves"],
            against,
            plotted["peak_count"] if plotted["multi"] else 0,
        )
        if not series:
            # The trend said so in the status bar when it drew nothing, but
            # that was several clicks ago, and the alternative is a dialog
            # offering to write a file with one empty column in it.
            self.status.showMessage(
                f"Nothing to export -- no curve has a {against} to sit at."
            )
            return

        wanted = fit_export_columns(plotted["fit_shape"] in PEAK_MIX)
        # Each set of columns carried beside the series it describes,
        # because "Integrate" measures a curve where the fit measures a
        # peak. With "Multi" on those are no longer the same thing: a
        # component holds no trapezoid -- peak_fit_multi blanks every
        # curve_stats key on one, which is the statement the table's own
        # blank cells make -- so the pair is taken from the whole curve
        # and written once per Y field rather than copied into each peak.
        # Off, the curve and the fit are one series and the two sets land
        # side by side on it, the second call never made.
        blocks = [(field, series[field], wanted) for field in sorted(series)]
        if plotted["integrated"]:
            whole = (
                trend_points(plotted["curves"], against, 0)[0]
                if plotted["multi"]
                else series
            )
            blocks += [
                (field, whole[field], CURVE_AREA_EXPORT)
                for field in sorted(whole)
            ]
        # One x between the series that share one, as export_table does for
        # the curves: the usual selection is one Y field across a row of
        # scans, which is a single column of x. Two fields measured on
        # different scans do not share it, and each group then carries its
        # own, named after the first series in it so the two can be told
        # apart.
        groups: dict[tuple[float, ...], list[tuple[str, list, list]]] = {}
        for field, points, spec in blocks:
            groups.setdefault(
                tuple(position for position, _stats in points), []
            ).append((field, points, spec))
        names: list[str] = []
        columns: list[np.ndarray] = []
        for positions, members in groups.items():
            names.append(
                f"x[{against}]"
                if len(groups) == 1
                else f"x[{against} @ {members[0][0]}]"
            )
            columns.append(np.asarray(positions, dtype=float))
            for field, points, spec in members:
                for column, key in spec:
                    names.append(f"{column}[{field}]")
                    columns.append(_fit_column(points, key))

        folder = self._export_folder()
        path, _filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export the fitted values",
            str(folder / f"fits_{self._export_scans(plotted)}.csv"),
            "Text data (*.csv *.dat *.txt);;All files (*)",
        )
        if not path:  # cancelled
            return
        if not Path(path).suffix:
            path += ".csv"
        # Nothing to keep honest the way _export_target does: there is one
        # kind of file here, so its prefix cannot come to disagree with
        # what is in it.
        target = Path(path)

        header = self._export_header(plotted) + [
            # What this file is, beside the provenance it shares with the
            # data export: the same fits the header lists one line per
            # curve, laid out one scan per row, which is the shape a plot
            # of them wants.
            f"fits_against = {against}",
        ]
        if plotted["integrated"]:
            # Two columns both named for an area is a question a reader
            # will ask, so it is answered here rather than left to the
            # area= line above, which explains the trapezoid but does not
            # know the fit is beside it.
            header.append(
                'curve_area = the "Integrate" trapezoids of the row that'
                ' was fitted -- the derivative\'s own while "Show dY/dX"'
                " is on, since that is what the fit describes"
                "   --   the area[..] column beside them is the fitted"
                " peak's instead, and the two do not agree"
            )
        if skipped:
            one = skipped == 1
            header.append(
                f"left_out = {skipped} curve{'' if one else 's'} with no"
                f" {against} of {'its' if one else 'their'} own"
            )
        try:
            write_columns(target, names, columns, header)
        except OSError as exc:
            self.status.showMessage(f"Export failed: {exc.strerror or exc}")
            return
        rows = max((column.size for column in columns), default=0)
        # Series and rows rather than curves and rows: the file is one row
        # per scan, so saying "8 curves" of a selection of eight would read
        # as the 8 columns the other export writes.
        self.status.showMessage(
            f"Exported {len(series)} fitted series, {rows} rows, to {target}"
        )

    def _export_folder(self) -> Path:
        """:meth:`_export_directory`, made -- or its parent if it cannot be.

        A directory that is not ours to write in is said once and then
        stepped out of, so that the dialog opens on a path that is there
        and the user can steer it anywhere from.
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

    def _export_directory(self) -> Path:
        """Where the save dialog opens: ``export`` beside the scans.

        Under the directory being browsed rather than one remembered, so a
        run's exports stay with the run and moving to another directory takes
        the exports with it. Nothing is created here -- the caller does that,
        and copes with a directory it may not write in.
        """
        return (self._directory or Path.cwd()) / "export"

    @staticmethod
    def _export_scans(plotted: dict) -> str:
        """``000171_000173``: the scans every exported name is built on."""
        tags = scan_tags(plotted["paths"])
        parts = [tags[str(path)] for path in plotted["paths"]]
        # A whole directory selected would otherwise name every scan in it.
        if len(parts) > 3:
            parts = [parts[0], "to", parts[-1]]
        return "_".join(_safe_name(part) for part in parts)

    @staticmethod
    def _export_stem(plotted: dict, picture: bool) -> str:
        """``data_000171_000173``, or ``figure_...`` for the picture of it."""
        prefix = "figure_" if picture else "data_"
        return prefix + HdfViewer._export_scans(plotted)

    @staticmethod
    def _export_name(plotted: dict) -> str:
        """The name the dialog opens with."""
        # A map has no numbers to offer, so what it offers is a picture.
        picture = plotted["mode"] == "grid"
        return HdfViewer._export_stem(plotted, picture) + (
            ".png" if picture else ".csv"
        )

    def _export_target(
        self, typed: Path, plotted: dict, picture: bool
    ) -> Path | None:
        """The file to write, with ``data_``/``figure_`` kept honest.

        The prefix says what is inside, but a name has to be offered before
        the suffix is settled -- so a suggestion left standing is put right
        here, while a name typed by hand is used exactly as typed. None when
        the rename would replace a file and the user would rather it did not.
        """
        if typed.stem != self._export_stem(plotted, not picture):
            return typed
        target = typed.with_name(
            self._export_stem(plotted, picture) + typed.suffix
        )
        # The dialog asked about replacing the name it showed, not this one.
        if target.exists():
            reply = QtWidgets.QMessageBox.question(
                self,
                "Replace file?",
                f"{target.name} already exists in {target.parent}.\n\n"
                f"The prefix follows what the file holds, so {typed.name} is "
                f"written as {target.name}. Replace it?",
                QtWidgets.QMessageBox.StandardButton.Yes
                | QtWidgets.QMessageBox.StandardButton.No,
            )
            if reply != QtWidgets.QMessageBox.StandardButton.Yes:
                return None
        return target

    @staticmethod
    def _export_header(plotted: dict) -> list[str]:
        """The provenance lines that go above the columns.

        Enough to say what the file is without the window that wrote it:
        which scans, read how, with what ticked, and -- reducing -- from
        which channels and how well the reduction went. The convention is
        ``_export_header`` in xmucd_live.py, and the reduction block below
        keeps its wording so the two kinds of file read alike.
        """
        paths = plotted["paths"]
        lines = [
            "scan_viewer export",
            f"written = {datetime.now().astimezone().isoformat(timespec='seconds')}",
            f"directory = {paths[0].parent}",
            f"files = {', '.join(path.name for path in paths)}",
            f"stream = {plotted['stream']}",
            f"mode = {plotted['mode']}",
            # The line already written above the plot, so the file and the
            # window cannot drift apart in how they describe the same view.
            f"view = {plotted['view']}",
        ]
        if plotted["mode"] == "xmcd":
            detector, monitor, transmission = plotted["reduction"]
            # A lock-in export is a different measurement reduced a different
            # way, so it names its own model and says which panel came off
            # which channel. Only where every unit agrees: a mixed overlay
            # has no one answer and the per-unit lines below carry it instead.
            modes = {
                result["mode"]
                for _tag, result, _missing in plotted["reports"]
                if result is not None
            }
            lockin = modes == {MODE_LOCKIN_POL}
            lines.append(
                "contrasts = "
                + ", ".join(CONTRAST_LABELS[key] for key in CONTRASTS)
            )
            if lockin:
                sources = CONTRAST_SOURCES[MODE_LOCKIN_POL]
                lines += [
                    "model = dc(s) = A + D*s, ac(s) = k*(C + E*s), s = +/-B"
                    "  --  XAS = mean(dc), XMuCD = dc(+) - dc(-),"
                    " XNCD = ac(+) + ac(-), XMCD = ac(+) - ac(-)",
                    "sources = "
                    + ", ".join(
                        f"{CONTRAST_LABELS[key]} {sources[key]}"
                        for key in CONTRASTS
                    )
                    + "  --  AC is the raw lock-in difference, neither"
                    " normalized nor logged, so not comparable with DC",
                    "scans = "
                    + ", ".join(tag for tag, _r, _m in plotted["reports"]),
                    "source = reduce_lockin_columns of xmucd_reduce.py,"
                    " reproduced here",
                ]
            else:
                lines += [
                    "model = mu(p, B) = A + C*p + D*B + E*p*B"
                    "  --  XAS = A, XNCD = 2C, XMuCD = 2D, XMCD = 2E",
                    "source = reduce_columns of xmucd_reduce.py, reproduced here",
                ]
            lines += [
                f"detector = {detector}",
                f"monitor = {monitor}",
                "signal = "
                + (
                    "log(monitor/detector)"
                    if transmission
                    else "detector/monitor"
                ),
            ]
            if lockin:
                lines.append(f"ac = {AC_FIELD} - {ACOFF_FIELD}")
            lines.append(
                f"energy_atol = {ENERGY_ATOL:g}   magnet_atol = {SPLIT_ATOL:g}"
                if lockin
                else f"energy_atol = {ENERGY_ATOL:g}   pr_atol = {PR_ATOL:g}"
                f"   magnet_atol = {SPLIT_ATOL:g}"
            )
            for tag, result, missing in plotted["reports"]:
                if result is None:
                    lines.append(
                        f"{tag} = not reduced, no {', '.join(missing)}"
                    )
                    continue
                lines.append(
                    f"{tag} = {result['n_complete']} of {result['n_energies']}"
                    f" energies, {result['n_points']} points,"
                    f" {result['n_finite']} finite,"
                    f" {result['n_steps']} field settings"
                    + (
                        f", joined at {result['conditions']}"
                        if result.get("conditions")
                        else ""
                    )
                )
                if result["n_complete"]:
                    # The leakage check the panels are read for; keeping it in
                    # the file means an exported curve can be judged without
                    # the window -- see _contrast_message for what it says,
                    # and xncd_partner for why it is not always XMuCD.
                    partner = xncd_partner(result)
                    lines.append(
                        f"{tag} = xncd_rms {_rms(result['xncd']):.6g}"
                        f" vs {partner}_rms {_rms(result[partner]):.6g}"
                    )
        else:
            lines += [
                f"x_field = {plotted['x_field']}",
                f"y_fields = {', '.join(plotted['y_fields'])}",
                "normalized_by = " + (plotted["monitor_field"] or "(none)"),
            ]
            # Named "signal" and worded as the reduction block above words
            # it, so that a transmission curve and a transmission panel are
            # read out of their two files the same way. Only where there is
            # a monitor: with none there is no form to be in.
            if plotted["monitor_field"]:
                lines.append(
                    "signal = "
                    + ("log(monitor/y)" if plotted["log_norm"] else "y/monitor")
                )
            if plotted["split_axis"]:
                lines.append(f"split_by = {plotted['split_axis']}")
            if plotted["right_scans"]:
                lines.append(
                    f"right_scale = {', '.join(plotted['right_scans'])}"
                )
            if plotted["inverted"]:
                lines.append(
                    f"inverted = {', '.join(plotted['inverted'])}"
                    "   --   y columns are -y as recorded"
                )
            if plotted["added"]:
                lines.append(
                    f"added = {', '.join(plotted['added'])}"
                    "   --   second scan interpolated onto the first one's x,"
                    " blank outside the overlap"
                )
            if plotted["averaged"]:
                lines.append(
                    f"averaged = {', '.join(plotted['averaged'])}"
                    "   --   repeats of one measurement, meaned row by row;"
                    " a scan stopped part way counts only where it reached"
                )
        # Written even when nothing went, and outside the branches above
        # because it applies to both: a reader asking why a curve has a gap
        # needs to know the filter ran, and one asking why it has none needs
        # to know it ran and found nothing. Per scan off the curves, one
        # total off the panels, which pool their files before the cut.
        if plotted["despiked"]:
            found = HdfViewer._spikes_removed(
                plotted["dropped"], plotted["removed"]
            )
            lines.append(
                f"spikes = {found}"
                "   --   removed where the monitor fell below"
                f" {BEAM_GONE_FRACTION:g} of its scan's median with both"
                " neighbours alive; a removed point is absent, never"
                " interpolated"
            )
        if plotted["derivative"]:
            lines.append(
                "derivative = the dy[..] columns, np.gradient(y, x) point by point"
            )
        # Written whenever the button was on, the explaining line even when
        # no curve was long enough to give a number -- for the reason the
        # spike line above is. A blank where a figure should be is a curve
        # of under two finite points, not a column that failed to run.
        if plotted["integrated"]:

            def area_of(value: float | None) -> str:
                return "(none)" if value is None else f"{value:.6g}"

            lines.append(
                "area = trapezoidal, x sorted ascending, over every point"
                " plotted"
                "   --   net is the same above the straight line joining the"
                " curve's first and last points; a gap is spanned, not skipped"
            )
            for label, area, net in plotted["areas"]:
                lines.append(
                    f"area[{label}] = {area_of(area)}   --   net {area_of(net)}"
                )
        # As above: the explaining line goes in whenever the button was on.
        # Here a blank means the fit was turned down as not describing a
        # peak, which is a result and worth recording as "(none)" rather
        # than an absent line. Only rows that were offered to the fit get a
        # line at all -- see where record["fits"] is filled.
        if plotted["fitted"]:

            def fit_of(value: float | None) -> str:
                return "(none)" if value is None else f"{value:.6g}"

            # Two shapes of model share these lines, so the explaining one
            # says which was used. Under "Multi" each line below is one
            # component of a joint fit, not a fit of its own: they were
            # found together over one background, so R2 is of the sum and
            # is repeated rather than apportioned, and the slope and offset
            # are the one line they share.
            lines.append(
                f"fit shape = {plotted['fit_shape']}"
                + (
                    f", {plotted['peak_count']} peaks"
                    if plotted["multi"]
                    else ""
                )
                + (
                    "   --   peaks + m*x + c, three parameters each over"
                    " one shared line, least squares, ordered by centre;"
                    " one line below per peak"
                    if plotted["multi"]
                    else "   --   peak + m*x + c, five parameters,"
                    " least squares"
                )
                + "; the height is above that sloping line, so it sits"
                " below Max Y and the width disagrees with the FWHM column;"
                " the area is height x fwhm x the shape's own integral over"
                " the whole line, so it is not the Net area trapezoid and"
                " for a Lorentzian it counts tails past the scan;"
                + (
                    " R2 is of the peaks summed, the line taken off both"
                    " sides, and is repeated on every line rather than"
                    " apportioned between them"
                    if plotted["multi"]
                    else " R2 is of the peak alone, the line taken off both"
                    " sides"
                )
                + (
                    "; eta is the pseudo-Voigt's fitted Lorentzian"
                    " fraction, 1 pure Lorentzian and 0 pure Gauss, and"
                    " a peak reporting one was fitted on a parameter"
                    " more than the count just given -- a peak without"
                    " one is the fixed half-and-half shape, which is"
                    " what was drawn for it"
                    if plotted["fit_shape"] in PEAK_MIX
                    else ""
                )
            )
            for (
                label,
                centre,
                fwhm,
                eta,
                height,
                area,
                r2,
                slope,
                offset,
            ) in plotted["fits"]:
                lines.append(
                    f"fit[{label}] = {fit_of(centre)}"
                    f"   --   fwhm {fit_of(fwhm)}"
                    # Left out rather than written "(none)", so that a
                    # shape with no mixing exports the line it always did.
                    + (f"   --   eta {fit_of(eta)}" if eta is not None else "")
                    + f"   --   height {fit_of(height)}"
                    f"   --   area {fit_of(area)}"
                    f"   --   r2 {fit_of(r2)}"
                    f"   --   slope {fit_of(slope)}"
                    f"   --   offset {fit_of(offset)}"
                )
        # Which quantity the second plot was trending, and nothing more: the
        # numbers it drew are the fit[..] lines above, scan by scan, so a
        # second copy here could only ever come to disagree with them.
        if plotted["trended"]:
            fields = ", ".join(field for field, _points in plotted["trend"])
            lines.append(
                f"trend = {plotted['trend_quantity']}"
                f"   --   against {plotted['trend_against']}"
                f", over {fields or '(nothing)'};"
                " the values are the fit[..] lines above"
            )
        lines += [
            f"curves = {len(plotted['curves'])}",
            "columns = x[..] then y[..] per curve, curves sharing an X sharing"
            " its column",
            "padding = a blank is past the end of a shorter curve, nan a"
            " reading that came back NaN",
        ]
        return lines

    def _stash_plot(self, record: dict | None):
        """Remember what was drawn, for Export, and let its button follow.

        The curve dicts in ``record`` are the drawing's own and not copies.
        That is safe only because every redraw builds them afresh and nothing
        touches them afterwards -- if anything ever did, the file Export
        writes would stop matching the plot it is named after.
        """
        self._plotted = record
        self.export_btn.setEnabled(record is not None)

    def _update_plot(self):
        # Cleared first so that no path out of here -- an unreadable file, a
        # selection with nothing ticked -- can leave the last plot behind to
        # be exported as though it were this one. The two exits that draw
        # something fill it back in.
        self._stash_plot(None)
        paths = self._selected_paths()
        # Pruned here rather than on selection changes: every path into the
        # viewer ends up redrawing, including the one that empties the listing
        # without a selection change to report.
        live = {str(path) for path in paths}
        self._right_paths &= live
        self._invert_paths &= live
        # A combination names two scans, so it survives only while both of
        # them do; and one sent to the right-hand scale is forgotten with it,
        # or asking for the same pair again would bring it back on a scale
        # that was never chosen for this one.
        self._combine = {
            key for key in self._combine if key[0] in live and key[1] in live
        }
        self._right_combos &= self._combine
        stream = self.stream_combo.currentText() or PRIMARY_STREAM
        x_field = self._x_field()
        y_fields = self._checked_y_fields()
        monitor_field = self._monitor_field()
        normalize = bool(monitor_field)

        # Which scans the four contrasts can be had from. Enabled from this
        # rather than from what is on screen, the way grid_btn is enabled from
        # the pre-filter grids below: a button that goes dead the moment it is
        # unticked could never be ticked back. One reducible scan is enough --
        # the others are left out with a note in the status bar.
        reducible = [path for path in paths if reduce_channels(path)]
        # A fixed-field dichro scan has the retarder stepped but the magnet
        # held still, so the four-quadrant reducer finds only two groups and
        # returns nothing. The beamline writer already wrote the result into
        # ``dichro_monitor``, though -- see :func:`has_dichro_stream` -- and
        # that is what the button opens on where the four-quadrant run cannot.
        on_dichro = stream == DICHRO_STREAM
        dichro_stream_scans = (
            [p for p in paths if has_dichro_stream(p)] if on_dichro else []
        )
        self.xmcd_btn.setEnabled(bool(reducible) or bool(dichro_stream_scans))
        reducing = bool(reducible) and self.xmcd_btn.isChecked()
        # Separate flag because reducing stays False here (there is no
        # four-quadrant split to run), yet the button is on and the plot
        # has to react: the stacked-panels rule below uses it to give XAS
        # and XMCD their own axes rather than overlaying them on one.
        dichro_stream_view = (
            self.xmcd_btn.isChecked()
            and not reducing
            and bool(dichro_stream_scans)
        )
        detector_channel = self.detector_combo.currentText()
        monitor_channel = self.monitor_combo.currentText()
        transmission = self.transmission_btn.isChecked()
        # Which way the Norm tick divides, where there is one. The button is
        # live in both cases, so there is no state here to leak from a mode
        # it is not applying to: with no Norm ticked and the panels off it is
        # dead and this is False whatever it still shows.
        log_norm = normalize and transmission
        # Nothing for them to choose while the panels they feed are not drawn.
        # The log is not one of them -- it governs the Norm tick as well, and
        # is enabled below, once the mode it is being asked about is known.
        for widget in (self.detector_combo, self.monitor_combo):
            widget.setEnabled(reducing)
        # The X axis is the reduction's own, not the ticked field: the
        # contrasts exist per energy and nothing else.
        x_label = ENERGY_FIELD if reducing else x_field

        # A two-motor raster can be shown either as colour maps -- one per
        # ticked Y field, side by side -- or, with the 2D map button off, as
        # ordinary curves. Read them either way so the button can tell the
        # user whether there is a raster to switch to. A dichro mesh scan has
        # one in its reduced stream too: XAS and XMCD over the same area.
        # Not while reducing, which would be file reads for a raster that is
        # not going to be drawn.
        grids: list[tuple[str, dict]] = []
        # A dichro mesh scan hits every cell twice -- once per retarder
        # setting -- so the raw per-cell read above has nothing to pick
        # between the two measurements. :func:`read_dichro_grid` splits them
        # the way :func:`reduce_dichro` splits one energy's measurements and
        # gives one XAS and one XMCD map in return. Tried only while
        # reducing: a scan that is not being asked for XMCD has no reason to
        # pay the file read, and the "2D map" button already governs the
        # fallback to per-cell counts.
        # Every selected scan over the same two motors as the first goes on
        # the one map, which is how a raster stopped part way is finished by
        # a second scan -- scans 39 and 41 of Ryan_26-3/sro -- see
        # :func:`merge_grids`. Matched on the motors rather than the grids'
        # axis names, which say readback for one scan and setpoint for the
        # next. The first scan still has to give a map of its own, as it did
        # when it was the only one drawn.
        motors = read_scan_info(paths[0])["motors"] if paths else None
        same_motors = [
            path for path in paths if read_scan_info(path)["motors"] == motors
        ]

        def merged(read) -> tuple[list[tuple[str, dict]], list[Path]]:
            """Each label's maps over *same_motors* merged, and the scans used."""
            by_label: dict[str, list[dict]] = {}
            used = []
            for path in same_motors:
                found = read(path)
                if not found:
                    if not used:
                        return [], []
                    continue
                used.append(path)
                for label, grid in found:
                    by_label.setdefault(label, []).append(grid)
            return [
                (label, merge_grids(found)) for label, found in by_label.items()
            ], used

        grid_paths: list[Path] = []
        dichro_grids = None
        if reducing and paths and is_mesh_scan(paths[0]):
            # Primary rather than the stream on screen for the same reason
            # :func:`reduce_channels` and :meth:`_read_contrasts` do: the raw
            # per-point retarder values the split needs are only in primary,
            # and a scan showing dichro_monitor still has its primary beside
            # it.
            dichro_grids, grid_paths = merged(
                lambda path: (
                    is_mesh_scan(path)
                    and read_dichro_grid(
                        path,
                        detector_channel,
                        monitor_channel,
                        PRIMARY_STREAM,
                        transmission,
                    )
                )
            )
        if dichro_grids:
            grids = dichro_grids
        elif paths and not reducing:
            grids, grid_paths = merged(
                lambda path: [
                    (field, grid)
                    for field in y_fields
                    if (
                        grid := read_grid(
                            path, field, monitor_field, stream, transmission
                        )
                    )
                    is not None
                ]
            )
        self.grid_btn.setEnabled(bool(grids))
        if not self.grid_btn.isChecked():
            grids = []
        # The dichro grid replaces the 1D reduction and so cancels it: the
        # four-panel contrast plot and the two-map one answer the same
        # question, so running both would read the file twice and draw a
        # plot the second one erases. Carried as its own flag because
        # reducing has to go down to keep the panels out of the dispatch
        # below, and nothing else there says a dichro grid is on screen.
        drawing_dichro_grid = bool(dichro_grids) and bool(grids)
        if drawing_dichro_grid:
            reducing = False
        # A mesh scan whose grid cannot be built -- scan 741 of NdFeB_26-1
        # was aborted after seven cells of a 31x21 raster; scan 753 defaults
        # to a dead detector against itself and gives all-NaN absorption --
        # has no map to draw, but the dichro_monitor stream still carries
        # one XAS and one XMCD value per cell that was measured. Fall back
        # to those as stacked 1D curves rather than running the four-quadrant
        # reducer, which has no field axis to split on here and shows four
        # empty panels under "Nothing reduced from this selection".
        elif (
            reducing
            and paths
            and is_mesh_scan(paths[0])
            and bool(dichro_stream_scans)
        ):
            reducing = False
            dichro_stream_view = True
        # Three ways to draw a selection, and the reduction wins where it is
        # asked for. grids is already empty above when it is, so the map
        # controls below go dead through the conditions they already have.
        mode = "xmcd" if reducing else ("grid" if grids else "curves")
        # Live wherever there are two channels to make an absorption out of:
        # the panels, which always have their pair; a curve or a map with a
        # Norm field ticked; and -- since Norm starts clear even on a scan
        # that nominates a monitor -- one whose hint can supply the tick,
        # which pressing the button then makes (see _on_transmission_toggled).
        # Dead only where no monitor is to be had at all, there being nothing
        # then for the log to be the form of: a raw count with nothing to
        # divide it is a choice that would change nothing on screen.
        #
        # The hint is asked for second, so a ticked Norm settles it without a
        # read: this runs on every redraw.
        self.transmission_btn.setEnabled(
            reducing
            or (
                mode in ("curves", "grid")
                and (normalize or bool(self._hinted_monitor()))
            )
        )
        # Something to judge the beam by. The panels always have one -- the
        # monitor channel is what they reduce against -- and a curve has
        # whichever field is normalising it, or failing that the one its scan
        # hints. Set here rather than with the rest below because the read
        # under it asks whether the button is live. A colour map is left out:
        # a raster's wild cells are what "Clip map scale" is for, and the two
        # would otherwise both claim the same spike.
        self.spike_btn.setEnabled(
            mode == "xmcd"
            or (
                mode == "curves"
                and bool(x_field)
                and bool(y_fields)
                and any(
                    self._spike_monitor(path, monitor_field) for path in paths
                )
            )
        )
        # Nothing to compare one scan against on its own, and a colour map has
        # no second Y scale to send it to.
        self.right_btn.setEnabled(
            mode == "curves"
            and len(paths) > 1
            and bool(x_field)
            and bool(y_fields)
        )
        # One scan on its own is worth inverting -- a dip read as a peak is
        # what lets the statistics table measure its centre and width -- so
        # unlike the right-hand scale this does not ask for two. Not while
        # the contrasts are shown: the choice is per scan and would take XAS,
        # an absorption positive by construction, negative along with XMCD.
        self.invert_btn.setEnabled(
            mode == "curves"
            and bool(paths)
            and bool(x_field)
            and bool(y_fields)
        )
        # Two scans to combine, and curves to combine them from.
        self.add_btn.setEnabled(
            mode == "curves"
            and len(paths) > 1
            and bool(x_field)
            and bool(y_fields)
        )
        # Two scans to find repeats among. Live in the panels as well, where
        # it joins the repeats into one reduction instead of one curve: 239
        # to 242 each step the field inside themselves, so each reduces alone
        # and the panels draw four copies of one measurement over each other
        # until this is ticked. Live even where none of them turn out to be
        # repeats, which is a thing only the files can answer and the status
        # bar reports: a button that greys itself out on the answer gives no
        # way to ask the question.
        self.average_scans_btn.setEnabled(
            (
                mode == "curves"
                and len(paths) > 1
                and bool(x_field)
                and bool(y_fields)
            )
            or (mode == "xmcd" and len(reducible) > 1)
        )
        # A colour map's statistics are a brightest and a dimmest cell;
        # there is no axis to integrate along, so the two columns would be
        # blank in every row. Curves and panels both integrate, which is
        # the point of putting the area in curve_stats rather than beside
        # it: the four contrasts get theirs with no code of their own.
        self.integrate_btn.setEnabled(mode != "grid")
        # Same rule and the same reason: a map has no axis to fit along
        # either. The combo follows the button rather than the mode, so the
        # shape cannot be changed while nothing is being fitted.
        self.fit_btn.setEnabled(mode != "grid")
        self.shape_combo.setEnabled(mode != "grid" and self.fit_btn.isChecked())
        # One cascade further down the same chain: the count modifies the
        # fit, so it follows the button that modifies it, which follows
        # the fit the shape already follows.
        self.multi_btn.setEnabled(mode != "grid" and self.fit_btn.isChecked())
        self.peaks_combo.setEnabled(
            self.multi_btn.isEnabled() and self.multi_btn.isChecked()
        )
        # Follows the fit rather than "Multi": a single fit has a
        # background to show even with no components to separate, and it
        # is the same background the five fit columns were reckoned
        # against.
        self.parts_btn.setEnabled(mode != "grid" and self.fit_btn.isChecked())
        # The trend plots the fitted numbers and only those, so it follows
        # the fit rather than the mode: with nothing fitted there is nothing
        # it is allowed to draw. It keeps its tick while greyed out, so that
        # unticking "Fit" and ticking it again brings the trend back as it
        # was, and its combo follows it the way the shape follows the fit.
        self.trend_btn.setEnabled(mode != "grid" and self.fit_btn.isChecked())
        self.trend_combo.setEnabled(
            self.trend_btn.isEnabled() and self.trend_btn.isChecked()
        )
        # Refilled only while the trend is on: its entries are read out of
        # the selected files, and a trend nobody is looking at has no
        # business opening 27 baselines on every redraw. Cached on the
        # selection, so every redraw after the first one is a lookup.
        if self.trend_btn.isChecked():
            self._refresh_trend_x_combo(paths)
        self.trend_x_combo.setEnabled(
            self.trend_btn.isEnabled() and self.trend_btn.isChecked()
        )
        # With the trend rather than with the fit, although there would be
        # numbers to write either way: what it writes them against is the
        # combo above, and that combo is only filled and only selectable
        # while the trend is on. Offering the export without it would be
        # offering a file whose x axis cannot be chosen.
        self.fit_export_btn.setEnabled(
            self.trend_btn.isEnabled() and self.trend_btn.isChecked()
        )
        self.clip_btn.setEnabled(bool(grids))
        # No amount to set when nothing is being clipped by it.
        self.clip_spin.setEnabled(bool(grids) and self.clip_btn.isChecked())
        self.cmap_combo.setEnabled(bool(grids))
        # Nothing to centre on a scale with no middle colour, so the button
        # goes dead unless at least one map on screen is drawn in one.
        self.center_btn.setEnabled(
            bool(grids)
            and any(
                self._map_colormap(field, stream) in DIVERGING_COLOR_MAPS
                for field, _ in grids
            )
        )

        # Read every curve before touching the figure: how many plots the
        # figure needs depends on whether the scan splits into several curves
        # per field, which is only known once its data has been read.
        curves: list[dict] = []
        split_axis = ""
        reports: list[tuple] = []
        averaged: list[str] = []
        # One count for the panels, where the reduction pools several files
        # into one answer and there is nothing finer to report; one per scan
        # for the curves, where each is drawn in its own right.
        dropped = 0
        removed: dict[str, int] = {}
        tags = scan_tags(paths)
        despiking = self.spike_btn.isEnabled() and self.spike_btn.isChecked()
        if reducing:
            curves, reports, averaged, dropped = self._read_contrasts(
                reducible,
                detector_channel,
                monitor_channel,
                transmission,
                self.average_scans_btn.isChecked(),
                despiking,
            )
        elif mode == "curves" and paths and x_field and y_fields:
            curves, split_axis, removed = self._read_curves(
                paths,
                stream,
                x_field,
                y_fields,
                monitor_field,
                transmission,
                despiking,
            )
        # Which reduction the panels are showing, where they all agree on
        # one. Only the lock-in answer is named: the four panels have been
        # the dichro ones since they were added and saying so now would be
        # noise, while a lock-in panel carries the lock-in's own units and an
        # overlay of both kinds has no single answer to give.
        reduced_modes = {
            result["mode"]
            for _tag, result, _missing in reports
            if result is not None
        }
        reduced_mode = (
            MODE_LOCKIN_POL if reduced_modes == {MODE_LOCKIN_POL} else ""
        )
        # split_axis stays empty while reducing, so this goes dead by itself.
        self.average_btn.setEnabled(bool(split_axis))
        if reducing:
            # Both rules below read the ticked Y fields, and while reducing
            # those ticks describe a plot that is not on screen. The panels
            # are the contrasts in their fixed order, already settled by the
            # reader, and nothing goes to a right-hand scale: every contrast
            # has a panel of its own, which is the only thing that scale is
            # for.
            stacked = True
            right_fields, right_scans = set(), set()
        else:
            # Several curves per field and several fields is a thicket on one
            # plot, so each field gets its own, stacked and sharing the X
            # axis. One field, or no split, keeps the single plot and its
            # right-hand scale -- an ordinary dichro energy scan looks as it
            # always did.
            # dichro_stream_view stacks whatever the tag-count rule would
            # say: on a dichro_monitor stream with XMCD on, the point is to
            # see XAS and XMCD on their own panels, same as the four-quadrant
            # path does for ordinary dichro scans.
            stacked = dichro_stream_view or (
                len({curve["tag"] for curve in curves}) > 1
                and len(y_fields) > 1
            )
            right_fields = set() if stacked else self._right_fields(y_fields)
            # Which scale each curve belongs to, settled before anything is
            # drawn: the axis label below, the Y labels and the legend all
            # have to agree with it, and they are written at three different
            # points. The two rules -- the stream's fields and the user's
            # scans -- are a union, so a dichro scan keeps XMCD on the right
            # whatever is ticked here. Stacking drops the field rule (every
            # field already has a plot of its own) but keeps the scan rule,
            # which separates curves within a plot.
            right_scans = self._right_scans(paths)
            for curve in curves:
                curve["panel"] = (
                    y_fields.index(curve["y_field"]) if stacked else 0
                )
                curve["on_right"] = (
                    curve["y_field"] in right_fields
                    or curve["path"] in right_scans
                )
            # _right_fields guards against emptying the main scale by counting
            # ticked fields; the scan rule cannot, because a selected scan need
            # not record the ticked field and so need not contribute a curve.
            # Decide it on the curves that exist, and back off the scan rule
            # alone -- dropping both would take XMCD off the scale its stream
            # chose. Globally rather than per plot: one scan on the right in
            # one plot of a stack and on the left in the next would be
            # unreadable.
            if curves and all(curve["on_right"] for curve in curves):
                right_scans = set()
                for curve in curves:
                    curve["on_right"] = curve["y_field"] in right_fields
        # Negating before combining, so that ticking 000173 inverted and
        # asking for 000171 + 000173 draws exactly 000171 - 000173. The other
        # order would make those two disagree and neither would be wrong.
        # A fresh array rather than y *= -1: the reader hands the same object
        # to more than one curve where a scan is split, and an in-place
        # negation would flip a neighbour too.
        inverted = self._invert_scans(paths) if mode == "curves" else set()
        for curve in curves:
            if curve["path"] in inverted:
                curve["y"] = -np.asarray(curve["y"], dtype=float)
        # Past the negation, so that a scan ticked upside down is averaged
        # upside down and the two buttons read in the order they are written;
        # past the panel and right-scale loop, so the curve that replaces a
        # group takes both from the members rather than work them out again.
        # Before combining: a sum of two scans that were each averaged is the
        # comparison worth having, and the reverse -- averaging a scan with a
        # sum it is itself a parent of -- is not.
        # Already filled in by the reduction above, where averaging joins
        # whole units rather than curves, so not reset here.
        # Gathered rather than shown as they are decided: two buttons can both
        # have something to say about one plot, and whichever spoke last would
        # otherwise be the only one heard.
        notes: list[str] = []
        if mode == "curves" and self.average_scans_btn.isChecked():
            curves, averaged = self._averaged_curves(curves, paths, tags)
            # Both ways round, because the status bar is where the button
            # gives its answer and a message only on the empty answer would
            # be left standing over the next selection, telling the user
            # there was nothing to average beside the curve that was.
            if averaged:
                notes.append(
                    f"Averaged {', '.join(averaged)}"
                    "  -  repeats of one measurement, meaned row by row."
                )
            else:
                notes.append(
                    "No two of the selected scans are the same measurement"
                    "  -  nothing to average."
                )
        # Named per scan rather than totalled: the count answers "is this
        # curve's gap real", which is a question about one scan and not about
        # the selection. Keyed only by the scans there was a monitor to ask
        # about -- see _read_curves -- so an empty answer here means the
        # filter ran and found nothing, never that it had nothing to run on.
        tagged = {
            tags[key]: count for key, count in removed.items() if key in tags
        }
        if despiking and removed:
            notes.append(
                f"Spikes removed  -  {self._spikes_removed(0, tagged)}"
            )
        # Appended past the panel and right-scale loop above -- which is why
        # _combined_curves settles both itself -- and past stacked, so that
        # asking for a sum cannot be what tips the plot into a stack.
        combined = (
            self._combined_curves(curves, paths, tags, y_fields, stacked)
            if mode == "curves"
            else []
        )
        curves += combined
        # What each curve is called, decided once here for all four places it
        # is named. The right-hand marker goes on every one of them: a peak in
        # the statistics table is found by matching the legend entry it
        # belongs to, and with two scans of one field the scale is the only
        # thing telling those entries apart -- the axis colouring below gives
        # up past a single curve. The legend and the table part company only
        # in a stack, where the legend sits on a plot whose own Y label names
        # the field and the table is one flat list over every plot, so the
        # table keeps the long form rather than repeat a string per plot.
        for curve in curves:
            # Which way up comes before which scale: the first says what the
            # numbers are, the second only where they are drawn. A sum of an
            # inverted scan is not marked -- it is a curve in its own right
            # and its name already says which two it came from.
            mark = (
                " (inverted)"
                if curve["path"] in inverted and "combo" not in curve
                else ""
            )
            mark += " (right)" if curve["on_right"] else ""
            curve["legend"] = (
                curve["short"] if stacked else curve["label"]
            ) + mark
            curve["stat"] = curve["label"] + mark
            curve["d_legend"] = f"d({curve['legend']})/d{x_label}"
            curve["d_stat"] = f"d({curve['stat']})/d{x_label}"
        right_tags = [tags[str(p)] for p in paths if str(p) in right_scans]
        # One clause listing everything on that scale rather than two that
        # have to be read together; a combination sent there is named the
        # same way its menu entry is.
        right_tags += [
            combine_name(tags[key[0]], tags[key[1]], key[2])
            for key in sorted(self._right_combos & self._combine_keys(paths))
        ]
        inverted_tags = [tags[str(p)] for p in paths if str(p) in inverted]
        added_tags = [
            combine_name(tags[key[0]], tags[key[1]], key[2])
            for key in sorted(self._combine_keys(paths))
        ]
        self._update_axis_label(
            right_fields,
            split_axis,
            grids[0][1] if grids else None,
            right_tags,
            inverted=inverted_tags,
            added=added_tags,
            averaged=averaged,
            reduction=(
                (detector_channel, monitor_channel, transmission)
                if reducing
                else None
            ),
            reduction_mode=reduced_mode,
            dropped=dropped,
            log_norm=log_norm,
        )

        # Preserve a manual zoom/pan across redraws triggered by checkbox or
        # dY/dX-toggle changes -- but only while the plot still shows the same
        # data plotted the same way (same files, X field and normalization),
        # otherwise old limits would be meaningless on the new data. Selecting
        # a different file is a different scan, so its range is autoscaled
        # afresh even when the X field happens to be the same one.
        path_key = tuple(str(p) for p in paths)
        # Energy does not move when the detector changes, but the four panels'
        # values can change by orders of magnitude, so a Y zoom made on one
        # pair would hold the XAS panel at a range the next pair never
        # reaches. Empty outside the reduction, where it compares equal to
        # itself and so changes nothing.
        reduce_key = (
            (detector_channel, monitor_channel, transmission)
            if reducing
            else ()
        )
        # Which scans went into each curve. Not folded into reduce_key, which
        # Export reads back as the three channel settings. A lock-in group
        # can change membership without the selection changing -- a scan
        # dropping out of a group leaves the rest drawing different numbers
        # on the same energies -- and a zoom kept across that is a zoom onto
        # data that is no longer there.
        unit_key = tuple(tag for tag, _result, _missing in reports)
        same_view = (
            path_key == self._last_paths
            and stream == self._last_stream
            and x_field == self._last_x_field
            and monitor_field == self._last_monitor_field
            # Beside the monitor field rather than folded into reduce_key,
            # which is the panels' own and empty here. A ratio of a few
            # thousand and the log of it are not the same numbers -- a Y
            # zoom made on one is off the top of the other -- so the form
            # counts as a different view exactly as the monitor does.
            and log_norm == self._last_log_norm
            and mode == self._last_mode
            and reduce_key == self._last_reduce_key
            and unit_key == self._last_unit_key
        )
        prev_xlim = prev_ylim = None
        if same_view and self.ax.lines:
            if not self.ax.get_autoscalex_on():
                prev_xlim = self.ax.get_xlim()
            # Sending a scan to the other scale changes what the left one
            # holds, so a Y zoom made while both shared it was sized to span
            # the offset between them and would hold the left scale at a range
            # it no longer needs -- the one thing the move is meant to undo.
            # The right-hand scale is never saved, so it autoscales regardless,
            # and keeping a stale left one would make the two sides disagree.
            # X is deliberately kept: having zoomed into a feature, moving a
            # scan off the scale is how you then compare the shapes there.
            # Inverting and adding join it for the same reason and on the
            # same side. Negating a curve moves its range to the other side
            # of zero and a sum sits an octave above its parts, both while
            # leaving X exactly where it was: a scan is still the same
            # energies upside down, and a sum is drawn on the first scan's
            # own X. A pair rather than the union of the two sets, which
            # would read a combination moved to the right as unchanged --
            # moving one there empties that much of the left scale, so it
            # belongs in the same test as asking for it in the first place.
            combined_key = (
                frozenset(self._combine),
                frozenset(self._right_combos),
            )
            if not self.ax.get_autoscaley_on() and (
                frozenset(right_scans) == self._last_right_scans
                and frozenset(inverted) == self._last_inverted
                and combined_key == self._last_combined
            ):
                prev_ylim = self.ax.get_ylim()
        self._last_right_scans = frozenset(right_scans)
        self._last_inverted = frozenset(inverted)
        self._last_combined = (
            frozenset(self._combine),
            frozenset(self._right_combos),
        )
        self._last_paths = path_key
        self._last_stream = stream
        self._last_x_field = x_field
        self._last_monitor_field = monitor_field
        self._last_log_norm = log_norm
        self._last_mode = mode
        self._last_reduce_key = reduce_key
        self._last_unit_key = unit_key

        # Everything Export needs to describe this plot, settled here where
        # all of it is known. The curves join it at the exit below, once the
        # drawing has worked out their derivatives -- a map has none to add,
        # its numbers being left out of the data export altogether.
        record = {
            "mode": mode,
            "curves": [],
            "paths": paths,
            "stream": stream,
            "x_label": x_label,
            "x_field": x_field,
            "y_fields": y_fields,
            "monitor_field": monitor_field,
            # Which way that monitor was applied. A flag beside the field
            # rather than left to be read out of the view line: a column of
            # numbers whose header says only which channel divided them
            # cannot be told from its own logarithm.
            "log_norm": log_norm,
            "split_axis": split_axis,
            "right_scans": [tags[path] for path in sorted(right_scans)],
            "inverted": inverted_tags,
            "added": added_tags,
            "averaged": averaged,
            # Both shapes, and the flag beside them: a header saying nothing
            # about spikes has to mean the button was off, not that it was on
            # and found none. Only the file can carry that -- by the time it
            # is read the button is long since not on screen.
            "despiked": despiking,
            "dropped": dropped,
            "removed": tagged,
            "derivative": self.derivative_btn.isChecked(),
            # Same reasoning as the spike flag above: the areas are filled
            # in at the exit below, off the same statistics the table got,
            # and the flag is what tells a header with no area lines that
            # the button was off rather than every curve too short to
            # integrate. Tied to the enable rule and not just the tick,
            # because a button left on and then greyed out by switching to
            # a map keeps its tick, and the map path never integrates.
            "integrated": mode != "grid" and self.integrate_btn.isChecked(),
            "areas": [],
            # And again for the fit, tied to the enable rule for the reason
            # the integral's is. The shape is recorded beside the flag as
            # its own field rather than left to be read out of the view
            # line, so the header says which model produced the numbers in
            # a form something other than a person can pick up.
            "fitted": mode != "grid" and self.fit_btn.isChecked(),
            "fit_shape": self.shape_combo.currentText(),
            # The count is recorded beside the shape for the same reason:
            # a table of two peaks apiece and a header that does not say
            # so is a table nobody can reproduce. Read as an integer here
            # so that the one place a combo string becomes a number is
            # this one.
            "multi": mode != "grid" and self.multi_btn.isChecked(),
            "peak_count": int(self.peaks_combo.currentText()),
            # Drawing only, which is why it joins the record and not the
            # header: the pieces are the fit the header already names,
            # split up, and every number in them is a column or an
            # exported line already.
            "fit_parts": mode != "grid" and self.parts_btn.isChecked(),
            "fits": [],
            # The trend plots the fit's own output, so it is off whenever
            # the fit is -- the enable rule greys it out, and this makes the
            # record agree rather than trusting a tick left behind.
            "trended": (
                mode != "grid"
                and self.fit_btn.isChecked()
                and self.trend_btn.isChecked()
            ),
            "trend_quantity": self.trend_combo.currentText(),
            # Read off the widget with the rest of them so that the axis
            # the header names cannot drift from the axis that was drawn.
            "trend_against": self.trend_x_combo.currentText(),
            "trend": [],
            "reduction": reduce_key or None,
            "reports": reports,
            "view": self.axis_label.text(),
        }

        # A colour bar owns its own axes, which ax.clear() does not touch.
        for bar in self._colorbars:
            bar.remove()
        self._colorbars = []
        # A heading over several maps belongs to the figure, so neither
        # ax.clear() nor the next draw takes it away on its own.
        self.figure.suptitle("")
        # Maps go in a row, stacked curves in a column; anything else is the
        # single plot the viewer normally shows.
        if mode == "grid":
            wanted, across = len(grids), True
        elif reducing:
            wanted, across = len(CONTRASTS), False
        else:
            wanted, across = (len(y_fields) if stacked else 1), False
        if self._set_panels(wanted, across):
            # Brand-new axes: whatever was zoomed to belonged to the old ones.
            prev_xlim = prev_ylim = None
        for (main, *_twins), default in zip(
            self._panels, self._default_format_coord
        ):
            # clear() resets the limits to 0..1 before the scale, so on a log
            # axis it warns about the zero; _apply_log_scales puts it back.
            main.set_xscale("linear")
            main.set_yscale("linear")
            main.clear()
            main.format_coord = default
        # Every plot and twin shares one X axis. matplotlib unions dataLim
        # across every axes in that shared group when autoscaling --
        # regardless of each one's own autoscale-on flag -- and clear() does
        # NOT reset dataLim (it's left stale until new data replaces it), so
        # a twin's last (e.g. derivative-off, file-A) data silently
        # polluted self.ax's autoscaled range for the next file. Forcing it
        # back to an empty bbox is what actually excludes it from that
        # union; the twins only ever get their X range explicitly afterward.
        for twin in self._twins:
            twin.set_yscale("linear")
            twin.clear()
            # clear() puts the Y label back on the left, where it lands on top
            # of self.ax's own; a twin's belongs with its ticks on the right.
            twin.yaxis.set_label_position("right")
            twin.yaxis.set_ticks_position("right")
            twin.dataLim = Bbox.null()
            twin.set_autoscalex_on(False)
        # The X and Y ticks are what the reduction does not use, so while it
        # is running an empty one is no reason to refuse to draw.
        if not paths or (not reducing and (not x_field or not y_fields)):
            self.ax.set_title(
                "Select file(s), an X field and at least one Y field."
            )
            for twin in self._twins:
                twin.set_visible(False)
            self.canvas.draw_idle()
            self.stats_table.setRowCount(0)
            # The headings would otherwise survive from the previous draw,
            # leaving thirteen of them stranded over an empty table.
            self._set_stats_columns(False, False)
            self._size_stats_table(0)
            # Neither of these two paths reaches _draw_trend, so the panel
            # would keep the previous selection's trend beside a plot that
            # no longer has anything to do with it.
            self.trend_widget.setVisible(False)
            return

        if grids:
            self._stash_plot(record)
            self.trend_widget.setVisible(False)
            # Dichro grids carry XAS and XMCD, which are already worked out
            # from the ticked detector/monitor pair, so dividing them again
            # in the colour-bar label would misrepresent what the colours
            # mean. The regular grids still name their ratio.
            self._plot_grids(
                paths,
                grid_paths,
                grids,
                "" if drawing_dichro_grid else monitor_field,
                False if drawing_dichro_grid else log_norm,
                stream,
                prev_xlim,
                prev_ylim,
            )
            return

        def scale_label(names: list[str], scans: list[str] = ()) -> str:
            # Each field wrapped in its own log rather than one log over the
            # joined list: "log(gI0 / a, b)" parses as a ratio against two
            # things at once, which is not what either curve is.
            # log_norm already implies a monitor -- it is only ever set
            # where one is ticked -- so it need not ask again.
            text = (
                ", ".join(f"log({monitor_field} / {name})" for name in names)
                if log_norm
                else ", ".join(names)
            )
            if normalize and text and not log_norm:
                text += f" / {monitor_field}"
            # Both scales can hold the same field, whole scans being what is
            # split; then the field name alone names neither of them.
            if text and scans:
                text += f" ({', '.join(scans)})"
            return text

        def side(
            own: list[dict], on_right: bool
        ) -> tuple[list[str], list[str]]:
            """The fields, and the scans, drawn on one plot's one scale.

            Taken from the curves actually drawn rather than by subtracting
            the right-hand fields from the ticked ones: with a scan on each
            scale that arithmetic puts the same field on both, which says
            nothing about either.
            """
            fields: list[str] = []
            scans: list[str] = []
            for curve in own:
                if curve["on_right"] is not on_right:
                    continue
                if curve["y_field"] not in fields:
                    fields.append(curve["y_field"])
                tag = tags.get(curve["path"], "")
                # Only worth naming when the scans were split by hand;
                # otherwise every scan is on the left and listing them all
                # would push the field name off the plot.
                if right_scans and tag and tag not in scans:
                    scans.append(tag)
            return fields, scans

        show_derivative = self.derivative_btn.isChecked()
        # The record's flag rather than the button's, so that the columns,
        # the drawn models and the header lines cannot disagree about
        # whether the fit ran. Counted as we go: a fit that the gate turns
        # down leaves a blank row, which is right but looks exactly like a
        # button that did nothing, so the tally below says which it was.
        fitting = record["fitted"]
        fit_shape = record["fit_shape"]
        multi = record["multi"]
        peak_count = record["peak_count"]
        fit_parts = record["fit_parts"]
        fitted = eligible = 0
        stats_rows: list[tuple[str, dict[str, float | None]]] = []
        for curve in curves:
            main, right, derivative = self._panels[curve["panel"]]
            axis = right if curve["on_right"] else main
            axis.plot(
                curve["x"],
                curve["y"],
                marker="o",
                markersize=3,
                linewidth=1,
                label=curve["legend"],
                color=curve["color"],
            )
            # Differentiated here rather than where the dashed line is drawn
            # further down, so that the gradient is taken once for all its
            # uses -- the dashed line, the statistics row, and now the fit --
            # and the derivative's statistics can sit directly under those of
            # the curve they came from. np.gradient needs two points to
            # difference, and a curve of one is left without a derivative;
            # it raises IndexError rather than returning empty if asked.
            curve["dy_dx"] = None
            d_stats = None
            if show_derivative and curve["x"].size >= 2:
                with np.errstate(divide="ignore", invalid="ignore"):
                    curve["dy_dx"] = np.gradient(curve["y"], curve["x"])
                # Worth a row of its own: on an absorption edge the derivative
                # is the peak, so its maximum is where the curve climbs most
                # steeply -- the edge position -- and its FWHM the edge's width.
                d_stats = curve_stats(curve["x"], curve["dy_dx"])
            stats = curve_stats(curve["x"], curve["y"])
            curve["fit"] = None
            if fitting:
                # Showing dY/dX moves the fit onto it rather than adding a
                # second one, so the curve's own row blanks and its model
                # line goes. That is the point of the pairing: a knife edge
                # or an absorption edge is a step, which holds no peak and
                # is declined, while its derivative is the beam profile or
                # the edge -- scans 101-106 fit at R2 0.99 on the derivative
                # and not at all on the curve. The converse is the trap: the
                # derivative of a peak is a lobe up then a lobe down, so the
                # fit takes one lobe, reports a displaced centre, and flips
                # the height's sign between identical repeats (scans 32-38,
                # R2 0.52 to 0.71). R2 is what tells the two apart.
                target = curve["dy_dx"] if show_derivative else curve["y"]
                into = d_stats if show_derivative else stats
                if target is not None:
                    eligible += 1
                    # The crossing width seeds the fitted one: it has just
                    # been worked out for this same data, and a seed near
                    # the answer is what keeps a five-parameter fit out of
                    # the wrong basin. Both come from `into`, so the seed
                    # always describes whatever is being fitted.
                    #
                    # "Multi" seeds itself instead, from find_peaks over
                    # both signs, and so takes no width from here: one
                    # crossing width describes one peak, and which of
                    # several it would describe is not a question with an
                    # answer. Nothing else about the pairing changes, so
                    # Multi and dY/dX compose -- which is the whole of the
                    # derivative-of-a-peak case the trap above describes.
                    into.update(
                        peak_fit_multi(
                            curve["x"], target, fit_shape, peak_count
                        )
                        if multi
                        else peak_fit(
                            curve["x"], target, fit_shape, into["fwhm"]
                        )
                    )
                    curve["fit"] = into
                    good = np.isfinite(curve["x"]) & np.isfinite(target)
                    model = None
                    # A Multi fit leaves fit_center blank on the parent --
                    # its numbers are on the component rows -- so the guard
                    # asks for either, and peak_model draws whichever it
                    # was given.
                    if into["fit_center"] is not None or into.get("fit_peaks"):
                        span = curve["x"][good]
                        grid = np.linspace(
                            span.min(), span.max(), PEAK_FIT_PLOT_POINTS
                        )
                        model = peak_model(grid, fit_shape, into)
                    if model is not None:
                        fitted += 1
                        # Onto whichever axis carries the data it describes:
                        # the derivative's own twin when that is what was
                        # fitted, otherwise the parent curve's axis, the
                        # expression repeated from the plot above rather
                        # than hoisted. The rule being kept is that a fit
                        # line never exists on an axis carrying no data
                        # line -- the fallback further down that gives
                        # self.ax the curves' own X range is guarded on
                        # self.ax having no lines, and only that rule keeps
                        # it correct. The twin obeys it too: the dashed
                        # derivative is drawn on it further down, before
                        # the fallback is reached. Put these on self.ax to
                        # share one scale, or on an overlay axis of their
                        # own, and that fallback quietly stops firing.
                        #
                        # Evaluated between the data's own ends, never past
                        # them, because the plots share an X axis. Solid,
                        # against the derivative's dashes, so the two read
                        # apart in one colour. The label is underscored so
                        # matplotlib leaves it out of the legend: seven
                        # curves would otherwise read as fourteen entries,
                        # each colour twice.
                        (derivative if show_derivative else axis).plot(
                            grid,
                            model,
                            linewidth=1,
                            color=curve["color"],
                            label=(
                                f"_fit {curve['d_legend']}"
                                if show_derivative
                                else f"_fit {curve['legend']}"
                            ),
                        )
                        # The same line, taken apart, on the same axis and
                        # in the same colour: background dash-dot and the
                        # components dotted, against the model's solid and
                        # the derivative's dashes, so four kinds of line in
                        # one colour still read apart. Underscored for the
                        # reason the model is, and more so -- at three
                        # peaks this is four more lines a curve, which
                        # would be twenty-eight legend entries for the
                        # seven curves the model alone makes fourteen of.
                        if fit_parts:
                            for index, part in enumerate(
                                peak_parts(grid, fit_shape, into)
                            ):
                                (derivative if show_derivative else axis).plot(
                                    grid,
                                    part,
                                    linewidth=1,
                                    color=curve["color"],
                                    linestyle="-." if index == 0 else ":",
                                    label="_fit part",
                                )
            stats_rows.append((curve["stat"], stats))
            if d_stats is not None:
                stats_rows.append((curve["d_stat"], d_stats))
            # One row a component, indented under whichever row was fitted
            # -- the derivative's when dY/dX is on, since that is what was
            # fitted -- and so appended last. Each carries the numbers the
            # parent's five fit cells left blank, and nothing else: a
            # component of a fit has no curve, so no maximum and no
            # trapezoid. They are ordinary (label, dict) pairs, which is
            # why the table below needs no change to show them.
            components = (curve["fit"] or {}).get("fit_peaks") or []
            for index, component in enumerate(components, start=1):
                stats_rows.append(
                    (
                        f"{STATS_PEAK_INDENT}peak {index} of {len(components)}",
                        component,
                    )
                )

        # Only the bottom plot of a stack carries the X label; sharex hides
        # the tick labels of the ones above it.
        self._panels[-1][0].set_xlabel(x_label)
        for index, (main, right, derivative) in enumerate(self._panels):
            own = [curve for curve in curves if curve["panel"] == index]
            left_fields, left_scans = side(own, False)
            if reducing:
                key = CONTRASTS[index]
                # Which channel the panel came off, where every unit on
                # screen is a lock-in one. Two of the four are then the raw
                # lock-in difference and two an absorption, so reading one
                # against the other by eye would be comparing nothing -- the
                # label says so rather than leaving it to be worked out.
                source = (
                    f" ({CONTRAST_SOURCES[reduced_mode][key]})"
                    if reduced_mode
                    else ""
                )
                main.set_ylabel(CONTRAST_LABELS[key] + source)
                # In the corner rather than under the label: four stacked
                # panels have no width to spare, and what the panel means is
                # wanted beside the curve rather than out in the margin.
                main.text(
                    0.99,
                    0.95,
                    CONTRAST_NOTES[key],
                    transform=main.transAxes,
                    ha="right",
                    va="top",
                    fontsize=7,
                    color="0.45",
                )
                # Headroom, so the autoscaled curve stays clear of that note.
                main.set_ymargin(0.18)
                if key != "xas":
                    # A difference is read against zero, which is off the
                    # bottom of a panel autoscaled to a small signal.
                    main.axhline(0.0, linewidth=0.8, color="0.6", zorder=0)
            else:
                main.set_ylabel(
                    scale_label(
                        [y_fields[index]] if stacked else left_fields,
                        left_scans,
                    )
                )
            main.grid(True, alpha=0.3)
            # clear() gives every plot its tick labels back; in a stack only
            # the bottom one, which carries the shared X axis, should have them.
            main.label_outer()
            # Resolve each plot's autoscaled view now instead of leaving it to
            # the deferred draw-time autoscale: the twins share the X axis and,
            # left to their own autoscale-on-draw, can win the shared-axis
            # update with their own (possibly stale/empty) data limits -- that
            # was the cause of the X range getting stuck on the derivative's
            # old range after toggling dY/dX off and then selecting a new file.
            main.relim()

            right_colors = [
                curve["color"] for curve in own if curve["on_right"]
            ]
            right.set_visible(bool(right_colors))
            if right_colors:
                right_names, right_tags_here = side(own, True)
                right.set_ylabel(scale_label(right_names, right_tags_here))
                # Resolve its range now rather than leaving it to the draw:
                # tight_layout below measures the tick labels, and a scale
                # still showing 0..1 would have the real ones laid out against
                # a margin sized for the placeholders.
                right.autoscale_view()
                # With a single curve over there, colouring the axis to match
                # says which of the two scales it is read against. clear() has
                # already put these back to the defaults for any other case.
                if len(right_colors) == 1:
                    right.yaxis.label.set_color(right_colors[0])
                    right.tick_params(axis="y", colors=right_colors[0])

            derivative.set_visible(show_derivative)
            if show_derivative:
                for curve in own:
                    dy_dx = curve["dy_dx"]
                    if dy_dx is None:
                        continue
                    x = curve["x"]
                    finite = np.isfinite(dy_dx)
                    derivative.plot(
                        x[finite],
                        dy_dx[finite],
                        linestyle="--",
                        linewidth=1,
                        color=curve["color"],
                        label=curve["d_legend"],
                    )
                derivative.set_ylabel(f"d/d{x_label}")
            # Both twins put their spine at the right edge, so the
            # derivative's steps outward to clear the stream's own axis when
            # both are shown -- by enough for that axis's tick labels and its
            # own ylabel.
            derivative.spines["right"].set_position(
                ("outward", 60 if show_derivative and right_colors else 0)
            )

        # Autoscale only once every plot's limits have been recomputed: they
        # share the X axis, so each one's range is the union over all of them.
        for main, *_twins in self._panels:
            main.autoscale_view()
        if not self.ax.lines and not self.ax2.lines and curves:
            # The left-hand field is ticked but no selected file records it,
            # so self.ax -- which owns the shared X axis -- has nothing to
            # scale itself by. Give it the range the other curves cover. The
            # right-hand scale is asked too: it shares the X axis and so can
            # supply the range itself, and curves live there now by choice
            # rather than only when a stream nominates the field.
            spans = [curve["x"][np.isfinite(curve["x"])] for curve in curves]
            spans = [span for span in spans if span.size]
            if spans:
                low = min(span.min() for span in spans)
                high = max(span.max() for span in spans)
                pad = 0.05 * (high - low) if high > low else 0.5
                # auto=None leaves the autoscale flag alone. Turning it off --
                # which is what set_xlim does by default -- would have the next
                # redraw read this fallback back as a user zoom and keep it.
                self.ax.set_xlim(low - pad, high + pad, auto=None)

        # One legend, on the top plot. Stacked, each curve's colour means the
        # same thing on every plot, so that one legend reads them all.
        top, *top_twins = self._panels[0]
        handles, labels = top.get_legend_handles_labels()
        for twin in top_twins:
            if twin.get_visible():
                extra_handles, extra_labels = twin.get_legend_handles_labels()
                handles, labels = handles + extra_handles, labels + extra_labels
        # The twins are made after the main plot and so are drawn over it,
        # legend and all. That went unnoticed while the right-hand scale
        # carried a small XMCD trace; a whole scan over there is a curve
        # across the whole frame and would run straight through the box. The
        # transparent patch is what keeps lifting the main plot from hiding
        # the twins' own curves behind its background.
        top.set_zorder(max(twin.get_zorder() for twin in top_twins) + 1)
        top.patch.set_visible(False)
        if handles:
            top.legend(handles, labels, fontsize=8, loc="best")
        self._apply_log_scales()
        if prev_xlim is not None:
            self.ax.set_xlim(prev_xlim)
        if prev_ylim is not None:
            self.ax.set_ylim(prev_ylim)

        # The twins' X autoscale is off (set above); give them self.ax's final
        # range explicitly. Only their Y axes -- the right-hand field's scale
        # and the derivative's -- stay autoscaled. emit=False keeps this from
        # propagating back over the shared X axis to self.ax, which would
        # clear self.ax's autoscale-X flag on every redraw -- the flag is how
        # the next redraw tells a real user zoom from a plain autoscaled view.
        for twin in self._twins:
            twin.set_xlim(self.ax.get_xlim(), emit=False)

        if reducing and not curves:
            # Four labelled but empty panels read as a scan with nothing in
            # it. Say on the plot that it was the reduction and not the scan
            # that came up empty; the status bar says which of the reasons.
            self.ax.set_title(
                "Nothing reduced from this selection -- see the status bar.",
                fontsize=9,
            )

        # The gate turns down a fit that does not describe a peak, which
        # leaves a blank row -- correct, and indistinguishable from a button
        # that did nothing, the same trap the despike and averaging notes
        # above are there to avoid. So the tally is said out loud whenever
        # some curve went unfitted, and -- whatever the count -- whenever
        # dY/dX moved the fit off the curves. The status bar only changes
        # when there is a note to show, so on the knife-edge scans (101-127,
        # a step on bI1 whose derivative is the beam profile) ticking dY/dX
        # turns ten refusals into ten good fits with nothing left to say,
        # and the refusal from the press before would still be sitting
        # there calling them unfitted.
        missed = eligible - fitted
        # "Multi" joins the two reasons a tally is worth printing even when
        # nothing failed: several rows have just appeared under one curve,
        # and the only thing that says how many were asked for -- as against
        # how many the data happened to hold -- is this line.
        if fitting and (missed or show_derivative or multi):
            # Named for what was actually fitted: with dY/dX shown the fit
            # is on the derivatives, and "curves" would point the reader at
            # the rows that are blank by design rather than by failure.
            what = "derivatives" if show_derivative else "curves"
            # The count goes in beside the shape, because "3 of 4" says
            # nothing about how many components each of the three has and
            # the rows below are about to show several apiece.
            how = f"{fit_shape}, {peak_count} peaks" if multi else fit_shape
            notes.append(
                f"Fitted {fitted} of {eligible} {what} ({how})"
                + (f"  -  {missed} gave no usable peak." if missed else ".")
            )

        record["trend"] = self._draw_trend(
            curves,
            record["trended"],
            record["trend_quantity"],
            show_derivative,
            notes,
            peak_count if multi else 0,
            record["trend_against"],
        )

        self.figure.tight_layout()
        self.canvas.draw_idle()

        if reducing:
            self.status.showMessage(
                self._contrast_message(
                    reports, detector_channel, monitor_channel
                )
            )
        elif notes:
            self.status.showMessage("   ".join(notes))
        # The record's flag rather than the button's own, so that the two
        # columns and the header lines can never disagree about whether the
        # integral ran.
        integrating = record["integrated"]
        self._update_stats_table(
            stats_rows,
            # Everything one curve can put in the table: its own row, the
            # derivative's when that is shown, and a row a component under
            # "Multi". Counted from the settings rather than from the rows
            # so that the height does not jump about as fits are declined.
            (2 if show_derivative else 1) + (peak_count if multi else 0),
            integrating,
            fitting,
            # The shape's own property and not the fit's: the column is
            # offered wherever a mixing could have been measured, and
            # stands blank on the rows where it was not. Tying it to the
            # answers instead would make the table grow and shrink a
            # column as scans are clicked through.
            fit_shape in PEAK_MIX,
        )
        # Taken off the rows the table was just given rather than integrated
        # a second time, so the file cannot disagree with the screen. The
        # derivative rows come along, which is wanted: over a whole curve
        # the integral of dY/dX is the rise from end to end, i.e. the edge
        # jump.
        if integrating:
            record["areas"] = [
                (label, stats["area"], stats["net_area"])
                for label, stats in stats_rows
            ]
        # The same again for the fits, with one difference: only the rows
        # that were offered to peak_fit at all. Exactly one of each pair is,
        # since showing dY/dX moves the fit onto the derivative rather than
        # adding a second one -- so the other row carries none of these keys,
        # and writing it out as a line of "(none)" would say "tried and
        # failed" about data that was never tried. Key presence is the only
        # discriminator the rows carry.
        if fitting:
            # Under "Multi" the line a reader wants is the component's, and
            # the parent carries none of those numbers -- writing it would
            # put "(none)" against a fit that succeeded. So a row holding
            # fit_peaks is stepped over and its name kept, to go in front of
            # the indented rows that follow it. A Multi fit the gate turned
            # down holds no fit_peaks and is written out as the blank line
            # it is, exactly as a declined single fit always has been.
            record["fits"] = []
            parent = ""
            for label, stats in stats_rows:
                child = label.startswith(STATS_PEAK_INDENT)
                if not child:
                    parent = label
                if "fit_center" not in stats or stats.get("fit_peaks"):
                    continue
                record["fits"].append(
                    (
                        f"{parent} {label.strip()}" if child else label,
                        stats.get("fit_center"),
                        stats.get("fit_fwhm"),
                        # Beside the width, as in the table: the two
                        # together are the shape, and _export_header
                        # leaves the line out where there is no mixing.
                        stats.get("fit_eta"),
                        stats.get("fit_height"),
                        stats.get("fit_area"),
                        stats.get("fit_r2"),
                        stats.get("fit_slope"),
                        stats.get("fit_offset"),
                    )
                )
        # Last, so that the curves carry the derivatives the loop above put on
        # them: Export writes a dy column beside every y while they are shown.
        record["curves"] = curves
        self._stash_plot(record)

    def _plot_grids(
        self,
        paths,
        grid_paths,
        grids,
        monitor_field,
        log_norm,
        stream,
        prev_xlim,
        prev_ylim,
    ):
        """Draw each ticked field's raster as a colour map, side by side.

        Every map keeps its own colour scale and bar. A dichro mesh scan's
        XMCD is a fraction of a percent of the XAS beside it, and on one
        shared scale it would be a flat wash of the middle colour.

        *grid_paths* are the selected scans the maps were merged from -- see
        :func:`merge_grids` -- and *paths* all those selected, so the heading
        can say which went in, how many cells more than one of them measured,
        and how many were left out for being over other motors.

        *stream* is wanted because the colour map each panel is drawn in
        depends on what the stream says about its field -- the grid dicts
        carry the data and not its name, so the panel loop's own field name is
        the only identity there is.

        *monitor_field* and *log_norm* name the normalization the cells were
        already given by :func:`read_grid`; nothing is divided here. They are
        passed only so that the title, the colour bar and the pointer readout
        say which of the two forms the numbers under them are in.
        """
        for twin in self._twins:
            twin.set_visible(False)
        first = grids[0][1]
        names = [path.name for path in grid_paths]
        if len(names) > 3:
            names = [names[0], f"{len(names) - 1} more"]
        named = (
            f"{' + '.join(names)}  —  "
            f"{len(first['x'])} x {len(first['y'])} grid"
        )
        notes = []
        if first.get("averaged"):
            notes.append(
                f"[{first['averaged']} cell(s) averaged over "
                f"{first['scans']} scans]"
            )
        left_out = len(paths) - len(grid_paths)
        if left_out:
            # Two rasters over different motors cannot share one pair of
            # axes, nor can a scan that is not a raster join one.
            notes.append(f"({left_out} selected scan(s) not on this map)")
        missing = int(np.count_nonzero(~np.isfinite(first["z"])))
        if missing:
            notes.append(f"[{missing} point(s) not measured]")
        heading = "   ".join([named, *notes])
        # One map labels itself; several need the heading above them all so
        # that each can say which field it is.
        if len(grids) > 1:
            self.figure.suptitle(heading, fontsize=9)
        # A single map is titled on its own axes, which are narrower than the
        # figure by the width of the colour bar, so one long line runs out from
        # over the plot and into the bar's power-of-ten header. The notes go on
        # a second line instead. The suptitle above has the full width.
        axes_heading = "\n".join([named, "   ".join(notes)]) if notes else named

        clip = self.clip_spin.value() if self.clip_btn.isChecked() else None
        meshes = []
        stats_rows = []
        for (main, *_twins), (value_field, grid) in zip(self._panels, grids):
            if not monitor_field:
                label = value_field
            elif log_norm:
                label = f"log({monitor_field} / {value_field})"
            else:
                label = f"{value_field} / {monitor_field}"
            x_edges, y_edges = axis_edges(grid["x"]), axis_edges(grid["y"])
            # Each map's scale is worked out from its own cells, so clipping
            # one spiky map does not touch the well-behaved one beside it.
            low, high, extend = (
                clipped_limits(grid["z"], clip)
                if clip is not None
                else (None, None, "neither")
            )
            cmap = self._map_colormap(value_field, stream)
            if self.center_btn.isChecked() and cmap in DIVERGING_COLOR_MAPS:
                low, high = centered_limits(grid["z"], low, high)
                # The widened scale cuts off less than the clipped one did,
                # so the arrows have to be asked again rather than carried.
                extend = scale_extend(grid["z"], low, high)
            meshes.append(
                (
                    main,
                    label,
                    extend,
                    main.pcolormesh(
                        x_edges,
                        y_edges,
                        grid["z"],
                        shading="flat",
                        vmin=low,
                        vmax=high,
                        cmap=cmap,
                    ),
                )
            )
            main.set_xlabel(grid["x_name"])
            main.set_ylabel(grid["y_name"])
            main.set_title(
                label if len(grids) > 1 else axes_heading, fontsize=9
            )
            main.grid(False)
            # Side by side the maps share both axes, so only the leftmost
            # needs its Y ticks and label.
            main.label_outer()

            # Autoscale cannot be relied on here: pcolormesh adds a
            # *collection*, which relim() does not look at, and the twin axes
            # sharing this X suppresses the draw-time rescale -- between them
            # the X range was left at the empty default while Y came out
            # right. Set both ends from the cell edges instead.
            main.set_xlim(x_edges[0], x_edges[-1])
            main.set_ylim(y_edges[0], y_edges[-1])

            # Report the value under the pointer too. The stock readout gives
            # only x and y (and, because of the twin axes, a second bare
            # coordinate pair), so replace it outright rather than append.
            def format_coord(
                x, y, grid=grid, label=label, x_edges=x_edges, y_edges=y_edges
            ):
                parts = [
                    f"{grid['x_name']}={x:.6g}",
                    f"{grid['y_name']}={y:.6g}",
                ]
                column = int(np.searchsorted(x_edges, x) - 1)
                row = int(np.searchsorted(y_edges, y) - 1)
                if 0 <= column < len(grid["x"]) and 0 <= row < len(grid["y"]):
                    value = grid["z"][row, column]
                    parts.append(
                        f"{label}="
                        + (
                            "no data"
                            if not np.isfinite(value)
                            else f"{value:.6g}"
                        )
                    )
                return "   ".join(parts)

            main.format_coord = format_coord
            stats_rows.append((label, grid_stats(grid)))

        self._apply_log_scales()
        if prev_xlim is not None:
            self.ax.set_xlim(prev_xlim)
        if prev_ylim is not None:
            self.ax.set_ylim(prev_ylim)
        self.ax2.set_xlim(self.ax.get_xlim(), emit=False)

        # Lay the axes out *before* adding the colour bars. tight_layout()
        # afterwards gives each bar's strip back to its plot and the bars end
        # up drawn on top of the maps. The rect leaves room for the heading,
        # which tight_layout() does not account for by itself.
        self.figure.tight_layout(
            rect=(0, 0, 1, 0.94) if len(grids) > 1 else (0, 0, 1, 1)
        )
        for main, label, extend, mesh in meshes:
            bar = self.figure.colorbar(mesh, ax=main, extend=extend)
            # XMCD runs to a few thousandths, and spelled out in full those
            # tick labels are wider than the bar's whole strip -- they ran off
            # the edge of the figure. A shared power above the bar keeps them
            # to a couple of characters. Ordinary ranges are untouched.
            compact = ScalarFormatter(useMathText=True)
            compact.set_powerlimits((-3, 4))
            bar.ax.yaxis.set_major_formatter(compact)
            # A single map is headed by the file, so its bar says which field
            # the colours are. Side by side each map is already titled with
            # its own, and a second copy on every bar only runs into the map
            # to its right.
            if len(grids) == 1:
                bar.set_label(label)
            self._colorbars.append(bar)
        self.canvas.draw_idle()

        def at(point):
            return "" if point is None else f"{point[0]:.6g}, {point[1]:.6g}"

        # A raster has no axis to integrate along, nor one to fit a peak
        # along, so both sets of extra columns are put away here whatever
        # the buttons say -- they are dead on a map in any case.
        self._set_stats_columns(False, False)
        self.stats_table.setRowCount(len(stats_rows))
        for row, (label, stats) in enumerate(stats_rows):
            values = [
                label,
                "" if stats["max_z"] is None else f"{stats['max_z']:.6g}",
                at(stats["max_at"]),
                "" if stats["min_z"] is None else f"{stats['min_z']:.6g}",
                at(stats["min_at"]),
                "",
                "",
            ]
            for column, text in enumerate(values):
                self.stats_table.setItem(
                    row, column, QtWidgets.QTableWidgetItem(text)
                )
        self.stats_table.resizeColumnsToContents()
        self._size_stats_table(len(stats_rows))

    def _size_stats_table(self, rows: int, per_curve: int = 1):
        """Size the statistics table for *rows*, *per_curve* of them to a curve.

        Two limits, whichever is the lower. ``STATS_TABLE_CURVES`` curves'
        worth is as much of the window as the table may ever take -- and with
        dY/dX shown a curve is *per_curve* rows, so the budget grows with it
        rather than halving how many curves fit. Below that the table is only
        as tall as what it is actually holding, so a one-curve scan, or
        turning the derivative back off, hands the space to the plot.
        """
        header = self.stats_table.horizontalHeader().height()
        row_height = self.stats_table.verticalHeader().defaultSectionSize()
        shown = min(rows, STATS_TABLE_CURVES * per_curve)
        # The columns are sized to their contents and the curve labels run
        # long, so past about nine columns a horizontal scrollbar appears --
        # and it takes its height from inside the viewport, which is most of
        # one of the three rows this is trying to show. Counted rather than
        # left to clip the bottom row to a sliver.
        scrollbar = self.stats_table.horizontalScrollBar()
        overflow = scrollbar.sizeHint().height() if scrollbar.isVisible() else 0
        self.stats_table.setMaximumHeight(
            header
            + shown * row_height
            + overflow
            + 2 * self.stats_table.frameWidth()
        )

    def _set_stats_columns(
        self, areas: bool, fits: bool = False, mixed: bool = False
    ):
        """Show the area and fit columns, or put the table back to seven.

        Called by every fill path -- the curves below, the maps above, and
        the bail-out that empties the table -- because they write their
        rows separately, and a column count changed in one and not the
        other puts every cell under the wrong heading rather than failing
        outright.

        Areas first, then fits, and the rows must be built in that same
        order. Both can be on at once, which is where getting it wrong
        would show.
        """
        headings = (
            STATS_COLUMNS
            + (STATS_AREA_COLUMNS if areas else [])
            + (fit_columns(mixed) if fits else [])
        )
        if self.stats_table.columnCount() != len(headings):
            self.stats_table.setColumnCount(len(headings))
        self.stats_table.setHorizontalHeaderLabels(headings)

    def _update_stats_table(
        self,
        rows: list[tuple[str, dict[str, float | None]]],
        per_curve: int = 1,
        areas: bool = False,
        fits: bool = False,
        mixed: bool = False,
    ):
        def fmt(value: float | None) -> str:
            return "" if value is None else f"{value:.6g}"

        self._set_stats_columns(areas, fits, mixed)
        self.stats_table.setRowCount(len(rows))
        for row, (label, stats) in enumerate(rows):
            values = [
                label,
                fmt(stats["max_y"]),
                fmt(stats["max_x"]),
                fmt(stats["min_y"]),
                fmt(stats["min_x"]),
                fmt(stats["center"]),
                fmt(stats["fwhm"]),
            ]
            if areas:
                values += [fmt(stats["area"]), fmt(stats["net_area"])]
            if fits:
                # .get, not subscripting: a dY/dX row never went through
                # peak_fit and genuinely has none of these keys, so [] would
                # raise on the first scan with the derivative shown. A blank
                # cell is the right answer for such a row in any case.
                values += [
                    fmt(stats.get("fit_center")),
                    fmt(stats.get("fit_fwhm")),
                ]
                # In the order fit_columns puts the headings, which is
                # why both of them ask it rather than agreeing by hand.
                if mixed:
                    values.append(fmt(stats.get("fit_eta")))
                values += [
                    fmt(stats.get("fit_height")),
                    fmt(stats.get("fit_area")),
                    fmt(stats.get("fit_r2")),
                ]
            for col, text in enumerate(values):
                self.stats_table.setItem(
                    row, col, QtWidgets.QTableWidgetItem(text)
                )
        self.stats_table.resizeColumnsToContents()
        self._size_stats_table(len(rows), per_curve)

    # -- Baseline -------------------------------------------------------

    def _show_baseline_window(self):
        """Open (or raise) the separate baseline window."""
        if self._baseline_window is None:
            window = QtWidgets.QDialog(self)
            window.setWindowTitle("Baseline (before / after)")
            window.resize(760, 560)
            layout = QtWidgets.QVBoxLayout(window)
            layout.setContentsMargins(6, 6, 6, 6)
            layout.addWidget(self._baseline_panel)
            self._baseline_window = window
        # Show first: _update_baseline() does nothing while the window is
        # hidden, so filling the table before showing it would leave it empty.
        self._baseline_window.show()
        self._baseline_window.raise_()
        self._baseline_window.activateWindow()
        self._refresh_baseline_combo(self._selected_paths())
        self._update_baseline()

    def _update_baseline(self):
        window = getattr(self, "_baseline_window", None)
        if window is None or not window.isVisible():
            return  # nothing on screen to fill in; reading the file can wait
        name = self.baseline_combo.currentText()
        self.baseline_table.setRowCount(0)
        if not name:
            return
        path = next((p for p in self._selected_paths() if p.name == name), None)
        if path is None:
            return
        try:
            baseline = read_baseline(path)
        except OSError as exc:
            self.status.showMessage(f"Could not read baseline of {name}: {exc}")
            return

        self.baseline_table.setSortingEnabled(False)
        self.baseline_table.setRowCount(len(baseline))
        for row, field in enumerate(sorted(baseline)):
            before, after = baseline[field]
            is_scalar = isinstance(before, float) and isinstance(after, float)
            delta = f"{after - before:.6g}" if is_scalar else ""
            values = [
                field,
                f"{before:.6g}" if isinstance(before, float) else str(before),
                f"{after:.6g}" if isinstance(after, float) else str(after),
                delta,
            ]
            for col, text in enumerate(values):
                self.baseline_table.setItem(
                    row, col, QtWidgets.QTableWidgetItem(text)
                )
        self.baseline_table.setSortingEnabled(True)
        self.baseline_table.resizeColumnsToContents()


def main():
    initial_directory = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    app = QtWidgets.QApplication(sys.argv)
    window = HdfViewer(initial_directory)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
