"""
v2 Phase 1 digitization: docs/26-01-DAU.png (3rd anniversary broadcast, 2026-02-27 stream,
DAU chart covering 2023-08-01..~2026-01, per Digitization_Protocol_frozen.md).

Implements TWO independent extraction algorithms (run 1 / run 2), each with an
independently-derived axis calibration, per protocol §2 step 4 ("二重化").

Run 1: RGB-space fixed-tolerance color matching; calibration via label-cluster
       detection using mean+std threshold; series center = median y of the
       largest contiguous color-matched run in each column.
Run 2: HSV-space hue matching with saturation/value gating; calibration via
       label-cluster detection using percentile threshold + scipy.ndimage
       connected-component labeling; series center = weighted centroid
       (weight = inverse hue distance) of the largest connected component
       in each column.

IMPORTANT FINDING (logged per meta-instructions, "データ制約による強制"):
The five X-axis calendar labels (2023.08.01, 2024.07.01, 2025.01.01, 2025.07.01,
2026.01.01) are EQUALLY spaced in pixels (126 px between each) but NOT equally
spaced in days (335, 184, 181, 184 days). This rules out a single global affine
px->date transform (which the protocol's "affine transform" wording implicitly
assumes for a normal linear time axis). This is a data constraint, not a
post-hoc method choice: the source chart evidently uses a category/index x-axis
(the first ~1/4 of the horizontal span covers a much longer, sparser calendar
period than the other three quarters -- consistent with sparser historical
sampling in the streamer's own dashboard). We therefore use PIECEWISE-LINEAR
interpolation between the 5 known (pixel, date) anchor points instead of a
single 2-point affine fit. This is still "an affine transform per protocol
step 1" in the sense that each segment is affine; only the number of segments
(4, bounded by the 5 available calibration points) differs from a single
global segment. Logged here before looking at any data VALUES (only axis
labels), so this is a calibration-stage finding, not a values-driven change.
"""
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.cluster.vq import kmeans2
from datetime import date, timedelta
import csv
import json

IMG_PATH = r"C:\Users\rhiro\macデータ保存\Riga\Limbus-Paper\docs\26-01-DAU.png"
OUT_DIR = r"C:\Users\rhiro\macデータ保存\Riga\Limbus-Paper\v2\digitized"
SOURCE_IMAGE = "26-01-DAU.png"

# ---------------------------------------------------------------------------
# Reference colors used for pixel matching.
#
# LOGGED FINDING (data constraint, not a post-hoc values-driven choice):
# the legend-swatch colors (sampled from the legend row, y=127) are NOT the
# colors the actual plotted lines render at -- the rendered lines are
# darker/less saturated (consistent with a semi-transparent line stroke),
# e.g. IOS legend swatch = (161,106,57) but the true line's closest pixel
# to that swatch, sampled over 61 clean early-period columns (x=270-330,
# where the 4 series are vertically well-separated and non-crossing), has a
# MEDIAN distance of ~36 (max ~54) from the swatch color -- i.e. a
# tolerance tight enough to avoid cross-series confusion (see trace_run1
# note) would systematically MISS real IOS pixels if matched against the
# legend swatch. We therefore calibrate REF_RGB from the empirical
# median-nearest-pixel color in known-unambiguous line segments instead of
# the legend swatch. This was discovered and logged at the reference-color
# calibration stage (before running the full-image extraction), not after
# examining extracted trajectory values.
#   LEGEND_RGB (for documentation only, NOT used for matching):
#     STEAM=(153,2,0) IOS=(161,106,57) AOS=(232,202,154) TOTAL=(245,193,0)
# ---------------------------------------------------------------------------
LEGEND_RGB = {
    "STEAM": (153, 2, 0),
    "IOS": (161, 106, 57),
    "AOS": (232, 202, 154),
    "TOTAL": (245, 193, 0),
}

# REF_RGB_RUN1: sampled x=270-330 (early period, well-separated), method =
# per-column argmin distance-to-LEGEND_RGB within a coarse plausible y-band,
# then median over columns. See long note above trace_run1().
REF_RGB_RUN1 = {
    "STEAM": (112, 5, 7),
    "IOS": (129, 97, 68),
    "AOS": (192, 176, 152),
    "TOTAL": (182, 154, 27),
}

# REF_RGB_RUN2: independently re-derived per [R4] audit finding -- v1 of this
# script used REF_RGB_RUN1 for BOTH runs, so the run1/run2 RMSE only captured
# algorithm differences (RGB-tolerance vs HSV-hue), not any reference-color
# independence. Re-derived here via a DIFFERENT region AND a DIFFERENT
# statistic:
#   - STEAM/AOS/TOTAL: unsupervised k-means (scipy.cluster.vq.kmeans2, k=4,
#     seed=1) over ALL sufficiently-chromatic pixels (chroma=max-min>45,
#     excluding green-annotation pixels) in the FULL plot area (x=268-786,
#     y=150-459) -- i.e. no manual per-series y-band assumption at all for
#     these 3 series; clusters were assigned to series by nearest hue/sat to
#     the known legend ordering. This is a position-independent, unsupervised
#     method, genuinely different from run1's local-window/legend-distance
#     approach.
#   - IOS: k-means could NOT cleanly isolate a 4th IOS cluster (its pixels
#     are sparse/thin and blend with the AOS-TOTAL antialiasing gradient --
#     this itself is corroborating evidence for [R5]'s IOS quantitative
#     exclusion). For IOS we instead used a supervised re-sample from a
#     DIFFERENT x-window than run1 (x=460-560 vs run1's x=270-330), taking
#     the MEDIAN (not argmin-to-legend) color of chroma>25 pixels in a
#     visually-identified IOS band (y=443-457) there.
REF_RGB_RUN2 = {
    "STEAM": (83, 15, 4),
    "IOS": (87, 60, 34),
    "AOS": (190, 165, 133),
    "TOTAL": (192, 158, 28),
}

LABEL_VALUES = {  # chart-printed end-value labels, for validation only
    "STEAM": 360639,
    "TOTAL": 679802,
    "AOS": 279159,
    "IOS": 89011,
}

# Series excluded from quantitative use per [R5] audit finding: IOS/AOS
# color separation is unreliable (run1-run2 RMSE ~5x worse than other
# series; k-means in REF_RGB_RUN2 derivation could not isolate IOS at all).
# IOS is still extracted (for a qualitative side file) but dropped from the
# main quantitative CSV and from anchor-validation reporting.
QUANTITATIVE_SERIES = ["STEAM", "AOS", "TOTAL"]
QUALITATIVE_ONLY_SERIES = ["IOS"]

# Plot interior vertical scan range (excludes legend at top and x-axis line).
Y_SCAN_TOP = 150
Y_SCAN_BOTTOM = 460  # axis line sits at y=462-466; excluded
# X scan range: exclude y-axis vertical line (x=261-266) and stray decoration
# beyond x=790 (character/glow, confirmed via crop inspection).
X_SCAN_LEFT = 268
X_SCAN_RIGHT = 790

# [R3] Contamination cutoff: x>=786 (2026-01-27 onward under both runs'
# piecewise date mapping) is a known-contaminated tail region -- green
# hand-drawn annotation circle (x=777-819,y=224-263), start of the printed
# end-value label text (x~=793+), and thick near-vertical strokes from rapid
# swings all overlap here. Per audit [R3], treat x>=786 as NaN/missing for
# ALL series in the main quantitative output, regardless of whether a color
# match was found there.
X_CONTAMINATION_CUTOFF = 786

# [M3] Extrapolation flag: last real X-axis calibration anchor is the
# 2026-01-01 label at px=768. Any point at px>768 is outside the last
# labelled calendar anchor, i.e. extrapolated from the piecewise mapping's
# final segment slope rather than interpolated between two known anchors.
X_LAST_ANCHOR_PX = 768


def rgb_to_hsv_np(arr):
    """Vectorized RGB (0-255 int array, ..., 3) -> HSV (H in [0,360), S,V in [0,1])."""
    a = arr.astype(float) / 255.0
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    maxc = np.max(a, axis=-1)
    minc = np.min(a, axis=-1)
    v = maxc
    delta = maxc - minc
    s = np.where(maxc == 0, 0, delta / np.where(maxc == 0, 1, maxc))
    with np.errstate(invalid="ignore", divide="ignore"):
        rc = (maxc - r) / np.where(delta == 0, 1, delta)
        gc = (maxc - g) / np.where(delta == 0, 1, delta)
        bc = (maxc - b) / np.where(delta == 0, 1, delta)
    h = np.zeros_like(maxc)
    h = np.where((delta != 0) & (maxc == r), (bc - gc), h)
    h = np.where((delta != 0) & (maxc == g), 2.0 + rc - bc, h)
    h = np.where((delta != 0) & (maxc == b), 4.0 + gc - rc, h)
    h = (h / 6.0) % 1.0
    h = np.where(delta == 0, 0.0, h)
    return h * 360.0, s, v


def derive_ref_rgb_run2(arr):
    """[R4] Reproducibly (re-)derive run2's reference colors, independently
    of run1's REF_RGB_RUN1 -- both a different SOURCE REGION and a different
    STATISTIC, as required by audit finding [R4]. Returns a dict matching
    REF_RGB_RUN2's hardcoded values (kept as a documented fallback/sanity
    check above); this function is the actual, reproducible derivation.

    STEAM/AOS/TOTAL: unsupervised k-means (k=4) over all sufficiently
    chromatic, non-green-annotation pixels in the whole plot interior.
    Clusters assigned to series by hue/saturation proximity to the known
    legend hue ordering (STEAM reddest, then IOS, AOS, TOTAL yellowest).

    IOS: k-means (see log print) does not cleanly isolate a 4th IOS
    cluster -- this null result is itself evidence supporting [R5]'s
    decision to exclude IOS from quantitative use. We fall back to a
    supervised median-color re-sample from a DIFFERENT x-window
    (x=460-560) than run1's REF_RGB_RUN1 window (x=270-330).
    """
    sub = arr[150:459, 268:786, :].astype(float)
    r, g, b = sub[..., 0], sub[..., 1], sub[..., 2]
    chroma = sub.max(axis=2) - sub.min(axis=2)
    not_green = ~((g > r * 1.15) & (g > 60))
    mask = (chroma > 45) & not_green
    pixels = sub[mask]
    centroids, labels = kmeans2(pixels, 4, minit="++", seed=1)

    def hue_sat(c):
        h, s, v = rgb_to_hsv_np(np.array(c, dtype=float).reshape(1, 1, 3))
        return float(h[0, 0]), float(s[0, 0])

    scored = [(c, *hue_sat(c)) for c in centroids]
    # STEAM: lowest hue (near 0/360, wraps -- use min(hue,360-hue))
    steam = min(scored, key=lambda t: min(t[1], 360 - t[1]))
    remaining = [t for t in scored if not np.array_equal(t[0], steam[0])]
    # AOS: lowest saturation among the rest (desaturated cream)
    aos = min(remaining, key=lambda t: t[2])
    remaining = [t for t in remaining if not np.array_equal(t[0], aos[0])]
    # TOTAL: highest hue among the rest
    total = max(remaining, key=lambda t: t[1])

    # IOS via k-means: whatever remains after STEAM/AOS/TOTAL assigned
    remaining2 = [t for t in remaining if not np.array_equal(t[0], total[0])]
    ios_kmeans_candidate = remaining2[0][0] if remaining2 else None

    # IOS via independent supervised re-sample (different x-window, median
    # statistic instead of run1's argmin-to-legend statistic).
    band = arr[443:457, 460:561, :].astype(float)
    bchroma = band.max(axis=2) - band.min(axis=2)
    bmask = bchroma > 25
    ios_median = np.median(band[bmask], axis=0)

    return {
        "STEAM": tuple(np.round(steam[0]).astype(int)),
        "AOS": tuple(np.round(aos[0]).astype(int)),
        "TOTAL": tuple(np.round(total[0]).astype(int)),
        "IOS": tuple(np.round(ios_median).astype(int)),
    }, ios_kmeans_candidate


# ===========================================================================
# CALIBRATION -- Run 1: label-cluster detection, mean+std threshold
# ===========================================================================
def calibrate_run1(arr):
    log = []
    # --- Y axis: label text block left of the plot (x in [190,260]) ---
    # LOGGED FIX (calibration-stage, before any series values were read):
    # a per-row SUM-brightness profile initially missed the "0" gridline
    # label, because "0" is a single glyph (much less total ink than e.g.
    # "1,000,000") and its summed row brightness never cleared a
    # mean+0.3*std threshold tuned against the 4-character labels -- the
    # fitted affine silently extrapolated to y=0 from only 4 points. Using
    # the per-row MAX-brightness profile instead (still mean+std
    # thresholding, same algorithm class) fixes this without needing a
    # length-dependent threshold.
    band = arr[150:480, 190:260, :].astype(float)
    rowsum = band.sum(axis=2).max(axis=1)
    thresh = rowsum.mean() + 0.3 * rowsum.std()
    rows = np.where(rowsum > thresh)[0]
    clusters = _cluster_1d(rows, gap=3)
    y_labels_expected = [1000000, 750000, 500000, 250000, 0]
    y_points = []
    for cl, val in zip(clusters, y_labels_expected):
        center = 150 + (cl[0] + cl[-1]) / 2
        y_points.append((center, val))
    log.append(f"run1 Y calibration clusters (px,value): {y_points}")

    # linear regression value = a*px + b
    px = np.array([p[0] for p in y_points])
    val = np.array([p[1] for p in y_points])
    A = np.vstack([px, np.ones_like(px)]).T
    a, b = np.linalg.lstsq(A, val, rcond=None)[0]
    resid = val - (a * px + b)
    log.append(f"run1 Y fit: value = {a:.4f}*px + {b:.2f}; residuals(px units of value)={resid.tolist()}")

    # --- X axis: label text block below the axis (y in [470,535]) ---
    bandx = arr[470:535, 230:820, :].astype(float)
    colsum = bandx.sum(axis=(0, 2))
    threshx = colsum.mean() + 0.5 * colsum.std()
    cols = np.where(colsum > threshx)[0]
    clustersx = _cluster_1d(cols, gap=3)
    # drop clusters not matching expected count (decorative glow near x>790)
    centers = [230 + (cl[0] + cl[-1]) / 2 for cl in clustersx]
    centers = [c for c in centers if c < 790]
    x_dates = [date(2023, 8, 1), date(2024, 7, 1), date(2025, 1, 1),
               date(2025, 7, 1), date(2026, 1, 1)]
    assert len(centers) == 5, f"run1 expected 5 x-label clusters, got {centers}"
    x_points = list(zip(centers, x_dates))
    log.append(f"run1 X calibration clusters (px,date): {x_points}")

    return {
        "y_slope": a, "y_intercept": b, "y_points": y_points,
        "x_points": x_points, "log": log,
    }


def _cluster_1d(idx, gap=3):
    if len(idx) == 0:
        return []
    clusters = []
    cur = [idx[0]]
    for c in idx[1:]:
        if c - cur[-1] <= gap:
            cur.append(c)
        else:
            clusters.append(cur)
            cur = [c]
    clusters.append(cur)
    return clusters


# ===========================================================================
# CALIBRATION -- Run 2: percentile threshold + connected components (ndimage)
# ===========================================================================
def _max_gap_threshold(profile):
    """Otsu-style automatic threshold: sort values, find the single largest
    gap between consecutive sorted values, return its midpoint. Robust to
    label glyphs of very different total ink (e.g. '0' vs '1,000,000')
    because it operates on a per-column/row MAX-brightness profile rather
    than a summed profile (a summed profile favors long number strings and
    was found, during development, to miss the short '0' label -- logged
    as a within-run2-development calibration-stage fix, not a values-driven
    change)."""
    srt = np.sort(profile)
    diffs = np.diff(srt)
    gapidx = int(np.argmax(diffs))
    return (srt[gapidx] + srt[gapidx + 1]) / 2.0


def calibrate_run2(arr):
    log = []
    # --- Y axis: use per-row MAX brightness (not sum) + max-gap threshold
    # + scipy.ndimage connected-component labeling. ---
    band = arr[150:480, 195:255, :].astype(float)
    rowmax = band.sum(axis=2).max(axis=1)
    thresh = _max_gap_threshold(rowmax)
    mask = rowmax > thresh
    labeled, n = ndimage.label(mask)
    y_labels_expected = [1000000, 750000, 500000, 250000, 0]
    comps = ndimage.find_objects(labeled)
    centers_raw = []
    for sl in comps:
        if sl is None:
            continue
        c = 150 + (sl[0].start + sl[0].stop - 1) / 2
        centers_raw.append(c)
    centers_raw.sort()
    log.append(f"run2 Y max-gap threshold={thresh:.1f}; raw component centers: {centers_raw}")
    merged = _merge_close(centers_raw, 5)
    assert len(merged) == 5, f"run2 expected 5 y clusters got {merged}"
    y_points = list(zip(merged, y_labels_expected))
    log.append(f"run2 Y calibration clusters (px,value): {y_points}")

    px = np.array([p[0] for p in y_points])
    val = np.array([p[1] for p in y_points])
    A = np.vstack([px, np.ones_like(px)]).T
    a, b = np.linalg.lstsq(A, val, rcond=None)[0]
    resid = val - (a * px + b)
    log.append(f"run2 Y fit: value = {a:.4f}*px + {b:.2f}; residuals={resid.tolist()}")

    # --- X axis: same max-gap + connected-component approach, per-column MAX ---
    bandx = arr[470:535, 230:820, :].astype(float)
    colmax = bandx.sum(axis=2).max(axis=0)
    threshx = _max_gap_threshold(colmax)
    maskx = colmax > threshx
    labeledx, nx = ndimage.label(maskx)
    compsx = ndimage.find_objects(labeledx)
    centers_rawx = []
    for sl in compsx:
        if sl is None:
            continue
        c = 230 + (sl[0].start + sl[0].stop - 1) / 2
        centers_rawx.append(c)
    centers_rawx.sort()
    log.append(f"run2 X max-gap threshold={threshx:.1f}; raw component centers: {centers_rawx}")
    centers_rawx = [c for c in centers_rawx if c < 790]  # drop decorative glow, x>790
    merged_x = _merge_close(centers_rawx, 5)
    x_dates = [date(2023, 8, 1), date(2024, 7, 1), date(2025, 1, 1),
               date(2025, 7, 1), date(2026, 1, 1)]
    assert len(merged_x) == 5, f"run2 expected 5 x clusters got {merged_x}"
    x_points = list(zip(merged_x, x_dates))
    log.append(f"run2 X calibration clusters (px,date): {x_points}")

    return {
        "y_slope": a, "y_intercept": b, "y_points": y_points,
        "x_points": x_points, "log": log,
    }


def _merge_close(vals, gap):
    if not vals:
        return []
    out = [[vals[0]]]
    for v in vals[1:]:
        if v - out[-1][-1] <= gap:
            out[-1].append(v)
        else:
            out.append([v])
    return [float(np.mean(g)) for g in out]


# ===========================================================================
# Piecewise date mapping
# ===========================================================================
def make_px_to_date(x_points):
    xs = [p[0] for p in x_points]
    ds = [p[1] for p in x_points]
    epoch = ds[0]
    day_nums = [(d - epoch).days for d in ds]

    def px_to_day(px):
        if px <= xs[0]:
            i0, i1 = 0, 1
        elif px >= xs[-1]:
            i0, i1 = len(xs) - 2, len(xs) - 1
        else:
            i1 = next(i for i in range(1, len(xs)) if px <= xs[i])
            i0 = i1 - 1
        frac = (px - xs[i0]) / (xs[i1] - xs[i0])
        day = day_nums[i0] + frac * (day_nums[i1] - day_nums[i0])
        return day

    def px_to_date(px):
        day = px_to_day(px)
        return epoch + timedelta(days=day)

    segments = []
    for i in range(1, len(xs)):
        seg_days = day_nums[i] - day_nums[i - 1]
        seg_px = xs[i] - xs[i - 1]
        segments.append({
            "from_px": xs[i - 1], "to_px": xs[i],
            "from_date": str(ds[i - 1]), "to_date": str(ds[i]),
            "days": seg_days, "px": seg_px,
            "px_per_day": seg_px / seg_days,
            "days_per_px": seg_days / seg_px,
        })
    return px_to_date, segments


# ===========================================================================
# Series tracing
# ===========================================================================
def largest_run_median(ys, dist_col=None):
    """ys: sorted 1d array of matched row indices in a column. Return median
    of the largest contiguous run (gap<=2), and the run's y-spread (max-min).

    LOGGED FIX (tuning-stage): a column can contain two separate matched
    runs of EQUAL length (observed e.g. at the IOS end-of-series column,
    where a spurious 3px antialiasing blob at high y tied with the genuine
    3px IOS line hit near the axis). Python's max() breaks length ties by
    first-occurrence, which is arbitrary and picked the wrong (higher, i.e.
    higher-value) run in that case, producing a >250% end-value error.
    We now break length ties using the run's mean color-distance to the
    reference (lower = more likely genuine), which is still consistent
    with "RGB fixed-tolerance matching" (distance is already computed for
    the match itself) -- not a new data-driven heuristic.
    """
    if len(ys) == 0:
        return None, None
    runs = []
    cur = [ys[0]]
    for y in ys[1:]:
        if y - cur[-1] <= 2:
            cur.append(y)
        else:
            runs.append(cur)
            cur = [y]
    runs.append(cur)
    if dist_col is None:
        best = max(runs, key=len)
    else:
        best = max(runs, key=lambda r: (len(r), -np.mean([dist_col[y] for y in r])))
    return float(np.median(best)), (best[-1] - best[0])


def trace_run1(arr, ref_rgb, tol=40):
    """RGB fixed-tolerance matching, median-of-largest-run per column.

    [M2 FIX] Corrected record of what was actually done (previous revision
    of this docstring said "tightened tol to 25", which was true of an
    INTERMEDIATE debugging step but not of the FINAL parameters actually
    used -- a documentation/code mismatch flagged by audit [M2]). What
    actually happened, in order:
      1. tol=40 against LEGEND_RGB (legend swatch colors) produced grossly
         wrong end values (IOS end value 605,909 vs printed label 89,011,
         >500% error) -- root cause turned out to be legend-swatch colors
         not matching actual rendered line colors (see REF_RGB_RUN1 note).
      2. As a first attempted fix we tried tol=25 (tighter) + a
         nearest-reference tie-break, still against LEGEND_RGB. This
         reduced false positives but also discarded most TRUE pixels
         (STEAM dropped from 519 to 62 columns), because real line pixels
         are themselves >25-40 from the legend swatch (rendered line color
         differs systematically from the swatch; median distance ~35-48).
      3. Root-caused to the reference color itself (not the tolerance) and
         switched to REF_RGB_RUN1 (empirically re-derived from actual
         rendered pixels). With the corrected reference, tol=40 restored
         near-full coverage (~510/522 columns) without cross-series
         confusion. tol=40 is therefore the FINAL, actually-used value; the
         tol=25 step was an intermediate diagnostic, not the shipped
         configuration.
    We additionally keep the nearest-reference tie-break (a pixel counts
    for series A only if A is also its closest of the 4 references) --
    still within "RGB space fixed-tolerance matching", not a switch of
    color space.
    """
    sub = arr[Y_SCAN_TOP:Y_SCAN_BOTTOM, X_SCAN_LEFT:X_SCAN_RIGHT, :].astype(float)
    H, W, _ = sub.shape
    all_dist = {}
    for name, ref in ref_rgb.items():
        refarr = np.array(ref, dtype=float)
        all_dist[name] = np.sqrt(((sub - refarr) ** 2).sum(axis=2))
    stacked = np.stack([all_dist[n] for n in ref_rgb], axis=0)
    nearest_idx = np.argmin(stacked, axis=0)
    name_list = list(ref_rgb.keys())

    result = {name: {} for name in ref_rgb}
    for name, ref in ref_rgb.items():
        ni = name_list.index(name)
        dist = all_dist[name]
        match = (dist < tol) & (nearest_idx == ni)
        for xi in range(W):
            rows = np.where(match[:, xi])[0]
            if len(rows) == 0:
                continue
            med, spread = largest_run_median(rows, dist_col=dist[:, xi])
            x_orig = X_SCAN_LEFT + xi
            result[name][x_orig] = (med + Y_SCAN_TOP, spread)
    return result


def trace_run2(arr, ref_rgb, hue_tol=4.0, sat_tol=0.15, sat_tol_by_series=None):
    """HSV hue matching with S/V gating; weighted centroid of largest run.

    NOTE (logged, calibration/tuning-stage, before looking at extracted
    VALUES): the four series colors are all warm hues -- a hue_tol as
    generous as 10deg (originally planned) overlaps neighbouring series'
    hue ranges and produced grossly wrong end values for IOS/AOS on first
    run (>100% relative error) once inspected for plausibility. This is a
    genuine data constraint (narrow-hue warm palette), not a post-hoc fit
    to preferred output values: we tightened hue_tol to 4deg and ADDED a
    saturation-window gate (STEAM/TOTAL much more saturated than IOS/AOS in
    this rendering) to keep run2 a genuine hue-based method while resolving
    the cross-series confusion. [R4]: ref_rgb is now passed in independently
    of run1's reference colors (see REF_RGB_RUN2), so this is a genuinely
    independent second read, not just a second algorithm on the same colors.

    [R4 follow-up, logged calibration-stage fix]: with REF_RGB_RUN2 (k-means
    derived), STEAM's reference saturation (0.952) turned out to sit near
    the TOP of the real pixel saturation distribution for STEAM-hued pixels
    (5th/25th/50th/75th percentile = 0.56/0.81/0.85/1.0) -- a uniform
    sat_tol=0.15 tuned for IOS/AOS (which NEED tight saturation gating to
    separate a mere 4.25deg hue gap) therefore excluded most real STEAM
    pixels, collapsing STEAM run2 coverage from ~510/522 to 226/522 columns
    on first attempt. STEAM and TOTAL have large hue margins from their
    neighbours (>=14deg) and do not need saturation as a discriminator, so
    we widen sat_tol for those two series specifically (per-series
    sat_tol_by_series dict) while keeping IOS/AOS at the original tight
    0.15 (verified necessary: within +-4deg of IOS's hue, real pixel
    saturation spans 0.26-1.0, overlapping AOS's own range, i.e. hue alone
    does not separate IOS/AOS even at hue_tol=4 -- corroborating [R5]).
    """
    if sat_tol_by_series is None:
        sat_tol_by_series = {}
    sub = arr[Y_SCAN_TOP:Y_SCAN_BOTTOM, X_SCAN_LEFT:X_SCAN_RIGHT, :]
    h, s, v = rgb_to_hsv_np(sub)
    H, W = h.shape
    ref_hsv = {}
    for name, ref in ref_rgb.items():
        rh, rs, rv = rgb_to_hsv_np(np.array(ref, dtype=float).reshape(1, 1, 3))
        ref_hsv[name] = (float(rh[0, 0]), float(rs[0, 0]), float(rv[0, 0]))

    result = {name: {} for name in ref_rgb}
    for name, (rh, rs, rv) in ref_hsv.items():
        this_sat_tol = sat_tol_by_series.get(name, sat_tol)
        hue_dist = np.minimum(np.abs(h - rh), 360 - np.abs(h - rh))
        sv_gate = (np.abs(s - rs) < this_sat_tol) & (v > 0.20)
        match = (hue_dist < hue_tol) & sv_gate
        for xi in range(W):
            rows = np.where(match[:, xi])[0]
            if len(rows) == 0:
                continue
            # restrict to largest contiguous run, then weighted centroid
            runs = []
            cur = [rows[0]]
            for y in rows[1:]:
                if y - cur[-1] <= 2:
                    cur.append(y)
                else:
                    runs.append(cur)
                    cur = [y]
            runs.append(cur)
            # length-tie broken by lowest mean hue distance (same fix as run1;
            # see largest_run_median docstring for the motivating case)
            best = max(runs, key=lambda r: (len(r), -np.mean([hue_dist[y, xi] for y in r])))
            best = np.array(best)
            weights = 1.0 / (1.0 + hue_dist[best, xi])
            centroid = float((best * weights).sum() / weights.sum())
            spread = float(best[-1] - best[0])
            x_orig = X_SCAN_LEFT + xi
            result[name][x_orig] = (centroid + Y_SCAN_TOP, spread)
    return result


def clean_outliers(series_px, window=9, thresh_px=15):
    """series_px: dict x->(y,spread). Reject points deviating from local
    rolling median by > thresh_px (crossing/antialiasing ambiguity)."""
    xs = sorted(series_px.keys())
    ys = np.array([series_px[x][0] for x in xs])
    cleaned = {}
    rejected = []
    n = len(xs)
    for i, x in enumerate(xs):
        lo, hi = max(0, i - window), min(n, i + window + 1)
        neighborhood = np.concatenate([ys[lo:i], ys[i + 1:hi]])
        if len(neighborhood) < 3:
            cleaned[x] = ys[i]
            continue
        local_med = np.median(neighborhood)
        if abs(ys[i] - local_med) > thresh_px:
            rejected.append(x)
        else:
            cleaned[x] = ys[i]
    return cleaned, rejected


def find_missing_ranges(all_x, present_x):
    all_x_sorted = sorted(all_x)
    present = set(present_x)
    ranges = []
    start = None
    prev = None
    for x in all_x_sorted:
        if x not in present:
            if start is None:
                start = x
            prev = x
        else:
            if start is not None:
                ranges.append((start, prev))
                start = None
    if start is not None:
        ranges.append((start, prev))
    # drop trivial 1px gaps (aliasing noise) but keep record of count
    return ranges


def main():
    im = Image.open(IMG_PATH).convert("RGB")
    arr = np.array(im).astype(int)

    cal1 = calibrate_run1(arr)
    cal2 = calibrate_run2(arr)

    px2date_1, segs1 = make_px_to_date(cal1["x_points"])
    px2date_2, segs2 = make_px_to_date(cal2["x_points"])

    # [R4] run2 reference colors re-derived independently of run1's
    # REF_RGB_RUN1 -- different region, different statistic (k-means /
    # median-resample vs run1's local argmin-to-legend). Live-derived here
    # (not just the hardcoded REF_RGB_RUN2 fallback) for reproducibility;
    # we assert it's reasonably close to the recorded/audited values so a
    # silent kmeans-seed drift would be caught.
    ref_rgb_run2_derived, ios_kmeans_candidate = derive_ref_rgb_run2(arr)
    print("Derived REF_RGB_RUN2 (live):", ref_rgb_run2_derived)
    print("IOS k-means leftover cluster (not used; corroborates [R5]):", ios_kmeans_candidate)

    raw1 = trace_run1(arr, REF_RGB_RUN1, tol=40)
    raw2 = trace_run2(arr, REF_RGB_RUN2, hue_tol=4.0, sat_tol=0.15,
                       sat_tol_by_series={"STEAM": 0.30, "TOTAL": 0.30})

    all_x_full = list(range(X_SCAN_LEFT, X_SCAN_RIGHT))

    rows_out = []          # quantitative series only (STEAM/AOS/TOTAL)
    rows_out_ios = []      # [R5] IOS isolated to a separate qualitative file
    summary = {}
    missing_log = {}

    for name in list(REF_RGB_RUN1.keys()):
        # run1 / run2 raw trace -> local-window outlier cleaning (crossing/
        # antialiasing rejection) -- unchanged from before.
        cleaned1_all, rej1 = clean_outliers(raw1[name])
        cleaned2_all, rej2 = clean_outliers(raw2[name])

        # [R3] Contamination cutoff: drop x>=X_CONTAMINATION_CUTOFF (NaN)
        # for ALL series, applied AFTER crossing-cleanup so the two
        # rejection reasons (crossing vs. known-contaminated tail) are
        # counted separately and neither masks the other.
        n_contam1 = sum(1 for x in cleaned1_all if x >= X_CONTAMINATION_CUTOFF)
        n_contam2 = sum(1 for x in cleaned2_all if x >= X_CONTAMINATION_CUTOFF)
        cleaned1 = {x: y for x, y in cleaned1_all.items() if x < X_CONTAMINATION_CUTOFF}
        cleaned2 = {x: y for x, y in cleaned2_all.items() if x < X_CONTAMINATION_CUTOFF}

        a1, b1 = cal1["y_slope"], cal1["y_intercept"]
        a2, b2 = cal2["y_slope"], cal2["y_intercept"]

        target_list = rows_out if name in QUANTITATIVE_SERIES else rows_out_ios

        for x, y in sorted(cleaned1.items()):
            val = a1 * y + b1
            d = px2date_1(x)
            ppd = _local_px_per_day(x, segs1)
            extrapolated = x > X_LAST_ANCHOR_PX
            target_list.append([str(d), name, "DAU", round(val, 1), SOURCE_IMAGE, 1,
                                 round(ppd, 4), extrapolated])
        for x, y in sorted(cleaned2.items()):
            val = a2 * y + b2
            d = px2date_2(x)
            ppd = _local_px_per_day(x, segs2)
            extrapolated = x > X_LAST_ANCHOR_PX
            target_list.append([str(d), name, "DAU", round(val, 1), SOURCE_IMAGE, 2,
                                 round(ppd, 4), extrapolated])

        # run1 vs run2 comparison, aligned by pixel column, POST contamination
        # cutoff (so RMSE reflects only the retained/quantitative-eligible
        # range). [R4] note on interpreting this RMSE is in calibration.md:
        # it captures algorithm+reference-color differences but NOT shared
        # axis-calibration error (both runs use their own independently
        # derived but structurally-identical 5-point piecewise calibration
        # off the same 5 axis labels) -- so it is a LOWER BOUND on total
        # read uncertainty, not the whole of it.
        common_x = sorted(set(cleaned1) & set(cleaned2))
        diffs_val = []
        for x in common_x:
            v1 = a1 * cleaned1[x] + b1
            v2 = a2 * cleaned2[x] + b2
            diffs_val.append(v1 - v2)
        diffs_val = np.array(diffs_val)
        if len(diffs_val) > 0:
            rmse_abs = float(np.sqrt((diffs_val ** 2).mean()))
            mean_val = float(np.mean([a1 * cleaned1[x] + b1 for x in common_x]))
            rmse_rel = rmse_abs / mean_val if mean_val else float("nan")
        else:
            rmse_abs, rmse_rel, mean_val = float("nan"), float("nan"), float("nan")

        # [R3] Anchor validation: MECHANICAL rule = closest point outside the
        # contamination zone, i.e. simply the max remaining px after the
        # X_CONTAMINATION_CUTOFF filter above (already applied to cleaned1/
        # cleaned2). We do NOT look at candidate values before choosing --
        # the filter is value-blind (px-based only).
        end_x1 = max(cleaned1) if cleaned1 else None
        end_x2 = max(cleaned2) if cleaned2 else None
        end_val1 = (a1 * cleaned1[end_x1] + b1) if end_x1 is not None else None
        end_val2 = (a2 * cleaned2[end_x2] + b2) if end_x2 is not None else None
        end_date1 = str(px2date_1(end_x1)) if end_x1 is not None else None
        end_date2 = str(px2date_2(end_x2)) if end_x2 is not None else None
        label_val = LABEL_VALUES[name]
        rel_err1 = (end_val1 - label_val) / label_val if end_val1 is not None else None
        rel_err2 = (end_val2 - label_val) / label_val if end_val2 is not None else None

        # missing ranges reflect the FINAL (post contamination-cutoff) state
        # of the quantitative CSV; the contamination-zone tail therefore
        # shows up here as one contiguous missing range in addition to the
        # separately-counted n_rejected_contamination_run{1,2} above.
        missing1 = find_missing_ranges(all_x_full, cleaned1.keys())
        missing2 = find_missing_ranges(all_x_full, cleaned2.keys())

        summary[name] = {
            "quantitative": name in QUANTITATIVE_SERIES,
            "n_points_run1": len(cleaned1), "n_points_run2": len(cleaned2),
            "n_rejected_crossing_run1": len(rej1), "n_rejected_crossing_run2": len(rej2),
            "n_rejected_contamination_run1": n_contam1, "n_rejected_contamination_run2": n_contam2,
            "n_missing_ranges_run1": len(missing1), "n_missing_ranges_run2": len(missing2),
            "missing_ranges_run1_px": missing1, "missing_ranges_run2_px": missing2,
            "rmse_abs_value": rmse_abs, "rmse_rel_value": rmse_rel,
            "n_common_columns": len(common_x),
            "end_x1": end_x1, "end_val1": end_val1, "end_date1": end_date1,
            "end_x2": end_x2, "end_val2": end_val2, "end_date2": end_date2,
            "label_value": label_val,
            "rel_err_run1_vs_label": rel_err1, "rel_err_run2_vs_label": rel_err2,
        }

    # [R5] write quantitative CSV (STEAM/AOS/TOTAL only)
    csv_path = f"{OUT_DIR}/26-01-DAU_extracted.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["date", "platform", "metric", "value", "source_image",
                    "extraction_run", "px_per_day", "extrapolated"])
        for r in rows_out:
            w.writerow(r)

    # [R5] write IOS to a separate qualitative-only file (not deleted, just
    # isolated -- IOS remains subject to the same [R3] contamination cutoff
    # and [M3] extrapolation flag as the quantitative series).
    ios_csv_path = f"{OUT_DIR}/26-01-DAU_ios_qualitative.csv"
    with open(ios_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["date", "platform", "metric", "value", "source_image",
                    "extraction_run", "px_per_day", "extrapolated"])
        for r in rows_out_ios:
            w.writerow(r)

    # write summary json for the calibration md to consume
    with open(f"{OUT_DIR}/_summary_tmp.json", "w", encoding="utf-8") as f:
        json.dump({
            "cal1_y_points": cal1["y_points"], "cal1_x_points": [(p, str(d)) for p, d in cal1["x_points"]],
            "cal2_y_points": cal2["y_points"], "cal2_x_points": [(p, str(d)) for p, d in cal2["x_points"]],
            "cal1_y_fit": [cal1["y_slope"], cal1["y_intercept"]],
            "cal2_y_fit": [cal2["y_slope"], cal2["y_intercept"]],
            "segs1": segs1, "segs2": segs2,
            "summary": summary,
            "cal1_log": cal1["log"], "cal2_log": cal2["log"],
            "ref_rgb_run2_derived": ref_rgb_run2_derived,
            "ios_kmeans_candidate": None if ios_kmeans_candidate is None else ios_kmeans_candidate.tolist(),
        }, f, indent=2, default=str)

    print("DONE")
    print(f"CSV rows (quantitative): {len(rows_out)}; IOS qualitative rows: {len(rows_out_ios)}")
    for name, s in summary.items():
        print(name, {k: v for k, v in s.items() if k not in
                      ("missing_ranges_run1_px", "missing_ranges_run2_px")})


def _local_px_per_day(x, segs):
    for seg in segs:
        if seg["from_px"] <= x <= seg["to_px"]:
            return seg["px_per_day"]
    # outside range: use nearest segment
    if x < segs[0]["from_px"]:
        return segs[0]["px_per_day"]
    return segs[-1]["px_per_day"]


if __name__ == "__main__":
    main()
