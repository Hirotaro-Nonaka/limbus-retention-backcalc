"""
v2 Phase 1 digitization: docs/26-01-MAU.png (3rd anniversary broadcast, 2026-02-27 stream,
MAU chart, X-axis labelled 2024-01..2025-07 but per protocol §2.5 shifted +6 months
to true 2024-07..2026-01, per Digitization_Protocol_frozen.md).

Lessons applied FROM THE START (per coordinator instruction), based on the 26-01-DAU
audit round-trip (see 26-01-DAU_calibration.md / 26-01-DAU_usage_conditions.md):
  - Reference colors for run1 and run2 are independently measured from the start
    (different region AND different statistic), not shared-then-fixed post-hoc.
  - Run-length tie-breaks in largest-run selection use mean color/hue distance
    (the DAU [bug] where Python's max() picked the wrong of two equal-length runs).
  - A mechanical, value-blind contamination cutoff is determined BEFORE looking at
    any series values (based on color-match tapering + known annotation bounding
    boxes), and applied uniformly to all series.
  - An `extrapolated` boolean column flags points outside the labelled axis range.
  - Vertices (monthly kink points) are detected and reported separately, since MAU
    is genuinely monthly data connected by straight (official) linear interpolation
    -- the vertices are the real data; the pixel-column trace between vertices is
    reading the CHART's own linear interpolation, not additional real information.

Run 1: RGB-space fixed-tolerance color matching (nearest-reference tie-break +
       run-length tie-break by mean distance); reference colors sampled from a
       LOCAL WINDOW (x=310-370, early-chart, well-separated series) via per-column
       argmin-distance-to-LEGEND-color, median over columns.
Run 2: HSV-space hue matching with per-series saturation gating; reference colors
       derived via UNSUPERVISED k-means (scipy.cluster.vq.kmeans2, k=6, seed=2) over
       chromatic pixels across the WHOLE plot interior, cluster-to-series assignment
       by nearest hue to the four canonical warm-palette anchors (0/28/37/47 deg),
       independent of run1's local-window colors.
"""
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.cluster.vq import kmeans2
from datetime import date, timedelta
import csv
import json

IMG_PATH = r"C:\Users\rhiro\MacData\Riga\Limbus-Paper\docs\26-01-MAU.png"
OUT_DIR = r"C:\Users\rhiro\MacData\Riga\Limbus-Paper\v2\digitized"
SOURCE_IMAGE = "26-01-MAU.png"

LEGEND_RGB = {
    "STEAM": (153, 1, 0),
    "IOS": (162, 108, 56),
    "AOS": (237, 202, 158),
    "TOTAL": (243, 195, 1),
}

# REF_RGB_RUN1: local-window (x=310-370) argmin-to-legend median, same method
# family as 26-01-DAU's run1 (but independently re-sampled on THIS image).
REF_RGB_RUN1 = {
    "STEAM": (130, 5, 6),
    "IOS": (154, 111, 69),
    "AOS": (229, 205, 173),
    "TOTAL": (233, 191, 24),
}

# REF_RGB_RUN2: NOT used directly (kept only for reference/discussion) --
# main() calls derive_ref_rgb_run2(arr) live and uses that as the single
# source of truth, to avoid the module constant silently drifting from what
# the documented derivation method actually produces (caught during dev:
# an earlier hand-copied version of this dict, from an exploratory run with
# slightly different region bounds and a greedy hue-only assignment, no
# longer matched the live derivation once the assignment was fixed to a
# global-optimal hue+saturation match).

LABEL_VALUES = {  # chart-printed end-value labels, for validation only
    "STEAM": 792250,
    "TOTAL": 1295599,
    "AOS": 626948,
    "IOS": 205263,
}

# Plot interior scan range. Y-axis vertical line at x=248-252 (excluded, start
# scan at 256). X-axis horizontal line at y=483-487 (excluded, Y_SCAN_BOTTOM=483).
Y_SCAN_TOP = 165
Y_SCAN_BOTTOM = 483
X_SCAN_LEFT = 256
X_SCAN_RIGHT = 792  # just past last real color match (~786), before label text (~800+)

# [Mechanical, value-blind] contamination cutoff: determined by finding, for
# each of the 4 series independently, the rightmost pixel column with ANY
# color match (loose tolerance 35) -- STEAM/TOTAL taper at x=785-786, IOS/AOS
# show later "matches" at x=847-848 that turned out (on pixel inspection) to
# be the printed end-value label glyphs (e.g. "626,948"), not chart lines.
# Cutoff set just past the STEAM/TOTAL taper point and confirmed against the
# green hand-annotation bounding box (x=575-799, y=150-222, covering "Un-",
# "130", and the TOTAL end-value underline) which starts overlapping the
# TOTAL line's high values before x=787. This was determined BEFORE tracing
# any series' full trajectory or comparing to labels.
X_CONTAMINATION_CUTOFF = 787

# [M3-style] Extrapolation flags: last real X-axis calibration anchor pixels.
X_FIRST_ANCHOR_PX = 280.0   # label 2024-01 -> corrected 2024-07-01
X_LAST_ANCHOR_PX = 788.0    # label 2025-07 -> corrected 2026-01-01

# [MAU-1 audit fix] IOS run-run RMSE (2.96%) is in the same order as the
# other 3 series (0.29-0.92%) -- qualitatively different from 26-01-DAU's
# IOS (20.7%) -- so IOS is included in the main quantitative CSV for this
# chart. calibration.md documents this per-chart decision explicitly so it
# is not confused with a general "IOS is quantitative" policy change.
QUANTITATIVE_SERIES_DEFAULT = ["STEAM", "AOS", "TOTAL", "IOS"]


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
    """Reproducible independent derivation of run2's reference colors.
    Unsupervised k-means (k=6, deterministic given fixed seed+input pixels)
    over chromatic, non-green pixels in the whole plot interior.

    [Fix, found during development, calibration-stage]: k=6 > 4 series, so
    some clusters are "extra" (antialiasing sub-populations, e.g. a darker
    STEAM fringe distinct from its core-color cluster). A GREEDY nearest-hue
    -only assignment (processing anchors in a fixed order, each claiming its
    nearest unclaimed cluster) is unstable: because IOS/AOS/an unused blend
    cluster all sit within ~5 degrees of hue of each other, greedy order
    could let AOS's true cluster get pre-empted, leaving IOS badly matched.
    We instead do a GLOBAL OPTIMAL assignment (brute force over all ways to
    choose+order 4-of-6 clusters, trivial at this size) minimizing total
    weighted (hue, saturation) distance to the LEGEND's own (hue,sat)
    anchors -- saturation strongly discriminates AOS (desaturated cream)
    from IOS (medium sat) even when their hues are only ~5 degrees apart.
    Anchors come from LEGEND_RGB (neutral ground truth), not from run1's
    derived colors, preserving run1/run2 independence.
    """
    sub = arr[Y_SCAN_TOP:Y_SCAN_BOTTOM, X_SCAN_LEFT:X_SCAN_RIGHT, :].astype(float)
    r, g, b = sub[..., 0], sub[..., 1], sub[..., 2]
    chroma = sub.max(axis=2) - sub.min(axis=2)
    not_green = ~((g > r * 1.15) & (g > 60))
    mask = (chroma > 30) & not_green
    pixels = sub[mask]
    centroids, labels = kmeans2(pixels, 6, minit="++", seed=2)

    def hue_sat(c):
        h, s, v = rgb_to_hsv_np(np.array(c, dtype=float).reshape(1, 1, 3))
        return float(h[0, 0]), float(s[0, 0])

    cluster_hs = [hue_sat(c) for c in centroids]
    anchors = {}
    for name, ref in LEGEND_RGB.items():
        h, s = hue_sat(np.array(ref, dtype=float))
        anchors[name] = (h, s)

    import itertools
    names = list(anchors.keys())
    best_assignment = None
    best_cost = None
    sat_weight = 60.0  # degrees-per-unit-saturation weighting
    for combo in itertools.permutations(range(len(centroids)), len(names)):
        cost = 0.0
        for name, ci in zip(names, combo):
            ah, asat = anchors[name]
            ch, csat = cluster_hs[ci]
            hd = min(abs(ch - ah), 360 - abs(ch - ah))
            cost += hd + sat_weight * abs(csat - asat)
        if best_cost is None or cost < best_cost:
            best_cost = cost
            best_assignment = combo

    result = {}
    for name, ci in zip(names, best_assignment):
        # [Fix, found during development]: rounding the centroid to integer
        # RGB before use shifted IOS's hue from ~29.6deg to ~35.5deg (close
        # to AOS's ~37deg anchor), because integer rounding of 3 correlated
        # channels is NOT hue-preserving at low chroma. This caused IOS
        # run2 to false-match AOS-adjacent antialiasing pixels (end value
        # off by ~196%). Keep full float precision for the reference color
        # actually used in hue/sat matching.
        result[name] = tuple(float(x) for x in centroids[ci])
    return result


# ===========================================================================
# CALIBRATION
# ===========================================================================
def _max_gap_threshold(profile):
    srt = np.sort(profile)
    diffs = np.diff(srt)
    gapidx = int(np.argmax(diffs))
    return (srt[gapidx] + srt[gapidx + 1]) / 2.0


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


def shift_6_months(d):
    """[Protocol §2.5, frozen] +6 month axis-mislabel correction for 26-01-MAU."""
    y, m = d.year, d.month + 6
    if m > 12:
        y += 1
        m -= 12
    return date(y, m, d.day)


def calibrate_run1(arr):
    """run1 calibration: per-row/col MAX-brightness profile, mean+std threshold,
    manual consecutive-index clustering (same family as 26-01-DAU run1, fixed
    from the start to use MAX not SUM profile per the DAU '0'-label lesson)."""
    log = []
    band = arr[150:500, 170:245, :].astype(float)
    rowmax = band.sum(axis=2).max(axis=1)
    thresh = rowmax.mean() + 0.3 * rowmax.std()
    rows = np.where(rowmax > thresh)[0]
    clusters = _cluster_1d(rows, gap=3)
    y_labels_expected = [1500000, 1000000, 500000, 0]
    y_points = []
    for cl, val in zip(clusters, y_labels_expected):
        center = 150 + (cl[0] + cl[-1]) / 2
        y_points.append((center, val))
    log.append(f"run1 Y calibration clusters (px,value): {y_points}")

    px = np.array([p[0] for p in y_points])
    val = np.array([p[1] for p in y_points])
    A = np.vstack([px, np.ones_like(px)]).T
    a, b = np.linalg.lstsq(A, val, rcond=None)[0]
    resid = val - (a * px + b)
    log.append(f"run1 Y fit: value = {a:.4f}*px + {b:.2f}; residuals={resid.tolist()}")

    bandx = arr[490:545, 240:850, :].astype(float)
    colmax = bandx.sum(axis=2).max(axis=0)
    threshx = _max_gap_threshold(colmax)
    idx = np.where(colmax > threshx)[0]
    clustersx = _cluster_1d(idx, gap=3)
    centers = [240 + (cl[0] + cl[-1]) / 2 for cl in clustersx]
    centers = [c for c in centers if 255 < c < 800]  # drop y-axis-line bleed (~240) and character/glow decoration (>800)
    label_dates_raw = [date(2024, 1, 1), date(2024, 7, 1), date(2025, 1, 1), date(2025, 7, 1)]
    label_dates_corrected = [shift_6_months(d) for d in label_dates_raw]
    assert len(centers) == 4, f"run1 expected 4 x-label clusters, got {centers}"
    x_points = list(zip(centers, label_dates_corrected))
    log.append(f"run1 X calibration clusters (px, corrected date): {x_points}")
    log.append(f"run1 X raw (uncorrected) labels: {list(zip(centers, label_dates_raw))}")

    return {"y_slope": a, "y_intercept": b, "y_points": y_points, "x_points": x_points, "log": log}


def calibrate_run2(arr):
    """run2 calibration: max-gap (Otsu-like) threshold + scipy.ndimage connected
    components -- independent detection algorithm from run1 (same family as
    26-01-DAU run2)."""
    log = []
    band = arr[150:500, 170:245, :].astype(float)
    rowmax = band.sum(axis=2).max(axis=1)
    thresh = _max_gap_threshold(rowmax)
    mask = rowmax > thresh
    labeled, n = ndimage.label(mask)
    comps = ndimage.find_objects(labeled)
    centers_raw = sorted(150 + (sl[0].start + sl[0].stop - 1) / 2 for sl in comps if sl is not None)
    y_labels_expected = [1500000, 1000000, 500000, 0]
    log.append(f"run2 Y raw component centers: {centers_raw}")
    assert len(centers_raw) == 4, f"run2 expected 4 y clusters got {centers_raw}"
    y_points = list(zip(centers_raw, y_labels_expected))
    log.append(f"run2 Y calibration clusters (px,value): {y_points}")

    px = np.array([p[0] for p in y_points])
    val = np.array([p[1] for p in y_points])
    A = np.vstack([px, np.ones_like(px)]).T
    a, b = np.linalg.lstsq(A, val, rcond=None)[0]
    resid = val - (a * px + b)
    log.append(f"run2 Y fit: value = {a:.4f}*px + {b:.2f}; residuals={resid.tolist()}")

    bandx = arr[490:545, 240:850, :].astype(float)
    colmax = bandx.sum(axis=2).max(axis=0)
    threshx = _max_gap_threshold(colmax)
    maskx = colmax > threshx
    labeledx, nx = ndimage.label(maskx)
    compsx = ndimage.find_objects(labeledx)
    centers_rawx = sorted(240 + (sl[0].start + sl[0].stop - 1) / 2 for sl in compsx if sl is not None)
    centers_rawx = [c for c in centers_rawx if 255 < c < 800]
    log.append(f"run2 X raw component centers: {centers_rawx}")
    label_dates_raw = [date(2024, 1, 1), date(2024, 7, 1), date(2025, 1, 1), date(2025, 7, 1)]
    label_dates_corrected = [shift_6_months(d) for d in label_dates_raw]
    assert len(centers_rawx) == 4, f"run2 expected 4 x clusters got {centers_rawx}"
    x_points = list(zip(centers_rawx, label_dates_corrected))
    log.append(f"run2 X calibration clusters (px, corrected date): {x_points}")

    return {"y_slope": a, "y_intercept": b, "y_points": y_points, "x_points": x_points, "log": log}


def make_px_to_date_affine(x_points):
    """Single global affine px->date fit (NOT piecewise) -- MAU's 4 calibration
    points, after the +6mo correction, have near-uniform day-per-pixel spacing
    (residuals <1.2 days across all 4 points when fit as one line), unlike
    26-01-DAU's X axis which was genuinely non-affine. This affine-vs-piecewise
    decision was made at the calibration stage (residual check below, before
    any series values were read) and is logged as the DAU-vs-MAU CONTRAST:
    same protocol step, different outcome because the underlying axis behaves
    differently -- this is NOT cherry-picking a nicer method, it's what each
    chart's own label spacing supports.
    """
    xs = np.array([p[0] for p in x_points])
    epoch = x_points[0][1]
    days = np.array([(d - epoch).days for d in [p[1] for p in x_points]], dtype=float)
    A = np.vstack([xs, np.ones_like(xs)]).T
    slope, intercept = np.linalg.lstsq(A, days, rcond=None)[0]  # days per px
    pred = slope * xs + intercept
    resid = days - pred

    def px_to_date(px):
        d = slope * px + intercept
        return epoch + timedelta(days=d)

    px_per_day = 1.0 / slope
    return px_to_date, {"slope_days_per_px": float(slope), "intercept_days": float(intercept),
                         "px_per_day": float(px_per_day), "residuals_days": resid.tolist(),
                         "epoch": str(epoch)}


# ===========================================================================
# Series tracing (same core logic as 26-01-DAU, with the tie-break fixes
# baked in from the start)
# ===========================================================================
def largest_run_median(ys, dist_col=None):
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


def trace_run2(arr, ref_rgb, hue_tol=4.0, sat_tol=0.15, sat_tol_by_series=None,
                val_tol_by_series=None):
    """[Fix, found during development]: MAU's actual rendered IOS/AOS pixel
    clusters differ in hue by only ~0.2deg (35.12 vs 35.29) -- MUCH closer
    than DAU's palette -- so hue provides essentially ZERO discrimination
    between them here; separation depends entirely on saturation (gap 0.29)
    and value/brightness (gap 0.26, AOS brighter). A stray antialiasing
    pixel elsewhere in the image landed inside IOS's default sat_tol=0.15
    window and produced a >190% end-value error (traced to pixel (357,770),
    color (88,70,45), which is not on the real IOS line at all). We add an
    optional per-series value(brightness) gate as a second discriminator
    alongside saturation, since IOS/AOS are well-separated in value even
    though barely separated in hue.
    """
    if sat_tol_by_series is None:
        sat_tol_by_series = {}
    if val_tol_by_series is None:
        val_tol_by_series = {}
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
        this_val_tol = val_tol_by_series.get(name, None)
        hue_dist = np.minimum(np.abs(h - rh), 360 - np.abs(h - rh))
        sv_gate = (np.abs(s - rs) < this_sat_tol) & (v > 0.15)
        if this_val_tol is not None:
            sv_gate = sv_gate & (np.abs(v - rv) < this_val_tol)
        match = (hue_dist < hue_tol) & sv_gate
        for xi in range(W):
            rows = np.where(match[:, xi])[0]
            if len(rows) == 0:
                continue
            runs = []
            cur = [rows[0]]
            for y in rows[1:]:
                if y - cur[-1] <= 2:
                    cur.append(y)
                else:
                    runs.append(cur)
                    cur = [y]
            runs.append(cur)
            best = max(runs, key=lambda r: (len(r), -np.mean([hue_dist[y, xi] for y in r])))
            best = np.array(best)
            weights = 1.0 / (1.0 + hue_dist[best, xi])
            centroid = float((best * weights).sum() / weights.sum())
            spread = float(best[-1] - best[0])
            x_orig = X_SCAN_LEFT + xi
            result[name][x_orig] = (centroid + Y_SCAN_TOP, spread)
    return result


def clean_outliers(series_px, window=9, thresh_px=15):
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
    return ranges


def _rdp(points, epsilon):
    """Ramer-Douglas-Peucker polyline simplification. points: list of (x,y)
    tuples with x,y on COMPARABLE scales (caller must normalize)."""
    if len(points) < 3:
        return points
    start, end = np.array(points[0]), np.array(points[-1])
    line = end - start
    line_len = np.hypot(*line)
    pts = np.array(points[1:-1])
    if line_len == 0:
        dists = np.hypot(pts[:, 0] - start[0], pts[:, 1] - start[1])
    else:
        dists = np.abs(line[0] * (start[1] - pts[:, 1]) - (start[0] - pts[:, 0]) * line[1]) / line_len
    if len(dists) == 0:
        return [tuple(start), tuple(end)]
    idx = int(np.argmax(dists))
    if dists[idx] > epsilon:
        left = _rdp(points[:idx + 2], epsilon)
        right = _rdp(points[idx + 1:], epsilon)
        return left[:-1] + right
    return [tuple(start), tuple(end)]


def detect_vertices(cleaned, epsilon_norm=4.0):
    """Detect monthly kink points via Ramer-Douglas-Peucker polyline
    simplification.

    [Fix, found during development]: a first attempt used local-extrema
    detection (sign change in smoothed first difference) directly on
    (pixel-x, people-value) pairs; this badly under/over-counted (e.g.
    STEAM run1: 2 vertices found vs. ~19-20 expected monthly points) --
    partly because flat multi-day plateaus in the RGB-tolerance trace
    (real: shallow-slope inter-vertex segments render as many identical
    rounded y pixels) don't have a clean single-pixel sign-change to key
    off of. Switched to RDP, which is designed for exactly this
    (approximate a polyline with fewer vertices within a distance budget).
    RDP also initially failed (found only 2 points down to epsilon=500)
    because pixel-index x (0-500ish) and people-value y (100,000s) are on
    wildly different scales, making the perpendicular-distance metric
    degenerate (effectively blind to y-deviations). Fixed by rescaling y
    to the same numeric range as x before running RDP. epsilon_norm=4 was
    chosen as a principled single value (not tuned per-series to hit the
    expected ~19-20 count): it is applied UNIFORMLY across all series/runs,
    and is of the same order as the run1-run2 RMSE converted to this
    normalized scale (~2.7 for STEAM). Resulting counts differ noticeably
    between run1 (close to 18-21, near the expected range) and run2
    (14-41, series-dependent) -- reported as-is, not adjusted to match
    the expected count.
    """
    xs = sorted(cleaned.keys())
    if len(xs) < 5:
        return []
    ys = np.array([cleaned[x] for x in xs], dtype=float)
    xs_arr = np.array(xs, dtype=float)
    x_span = xs_arr.max() - xs_arr.min()
    y_span = ys.max() - ys.min()
    if y_span == 0 or x_span == 0:
        return []
    y_scaled = (ys - ys.min()) / y_span * x_span
    pts = list(zip(xs_arr, y_scaled))
    simplified = _rdp(pts, epsilon_norm)
    simplified_x = set(round(p[0]) for p in simplified)
    # map back to original (x, y_pixel) pairs, excluding the two polyline
    # endpoints (not "kinks", just the trace boundary)
    xs_sorted = xs
    vertices = [(x, cleaned[x]) for x in xs_sorted if x in simplified_x]
    if len(vertices) >= 2:
        vertices = vertices[1:-1]
    return vertices


# ===========================================================================
def main():
    im = Image.open(IMG_PATH).convert("RGB")
    arr = np.array(im).astype(int)

    cal1 = calibrate_run1(arr)
    cal2 = calibrate_run2(arr)

    px2date_1, xfit1 = make_px_to_date_affine(cal1["x_points"])
    px2date_2, xfit2 = make_px_to_date_affine(cal2["x_points"])

    print("X affine fit run1:", xfit1)
    print("X affine fit run2:", xfit2)

    # [Fix] Use the LIVE-derived reference colors as the single source of
    # truth for trace_run2 (not a hardcoded module constant) -- avoids any
    # drift between the documented derivation method and what's actually
    # used, which we caught during development (module-level REF_RGB_RUN2
    # had been hand-copied from an earlier exploratory run with slightly
    # different pixel-region boundaries and a greedy hue-only assignment;
    # it silently no longer matched what derive_ref_rgb_run2() produces).
    ref_rgb_run2_derived = derive_ref_rgb_run2(arr)
    print("Derived REF_RGB_RUN2 (live, used for trace_run2):", ref_rgb_run2_derived)

    with open(f"{OUT_DIR}/_calib_stage_tmp.json", "w", encoding="utf-8") as f:
        json.dump({
            "cal1_y_points": cal1["y_points"], "cal1_x_points": [(p, str(d)) for p, d in cal1["x_points"]],
            "cal2_y_points": cal2["y_points"], "cal2_x_points": [(p, str(d)) for p, d in cal2["x_points"]],
            "cal1_y_fit": [cal1["y_slope"], cal1["y_intercept"]],
            "cal2_y_fit": [cal2["y_slope"], cal2["y_intercept"]],
            "xfit1": xfit1, "xfit2": xfit2,
            "cal1_log": cal1["log"], "cal2_log": cal2["log"],
            "ref_rgb_run2_derived": ref_rgb_run2_derived,
        }, f, indent=2, default=str)
    print("Calibration stage saved to _calib_stage_tmp.json")

    raw1 = trace_run1(arr, REF_RGB_RUN1, tol=40)
    raw2 = trace_run2(arr, ref_rgb_run2_derived, hue_tol=4.0, sat_tol=0.15,
                       sat_tol_by_series={"STEAM": 0.30, "TOTAL": 0.30, "IOS": 0.10, "AOS": 0.10},
                       val_tol_by_series={"IOS": 0.15, "AOS": 0.15})

    all_x_full = list(range(X_SCAN_LEFT, X_SCAN_RIGHT))

    rows_out = []
    rows_out_ios = []
    summary = {}
    vertices_all = {}

    for name in list(REF_RGB_RUN1.keys()):
        cleaned1_all, rej1 = clean_outliers(raw1[name])
        cleaned2_all, rej2 = clean_outliers(raw2[name])

        n_contam1 = sum(1 for x in cleaned1_all if x >= X_CONTAMINATION_CUTOFF)
        n_contam2 = sum(1 for x in cleaned2_all if x >= X_CONTAMINATION_CUTOFF)
        cleaned1 = {x: y for x, y in cleaned1_all.items() if x < X_CONTAMINATION_CUTOFF}
        cleaned2 = {x: y for x, y in cleaned2_all.items() if x < X_CONTAMINATION_CUTOFF}

        a1, b1 = cal1["y_slope"], cal1["y_intercept"]
        a2, b2 = cal2["y_slope"], cal2["y_intercept"]

        vtx1 = detect_vertices(cleaned1)
        vtx2 = detect_vertices(cleaned2)
        vertices_all[name] = {
            "run1": [(x, str(px2date_1(x)), round(a1 * y + b1, 1)) for x, y in vtx1],
            "run2": [(x, str(px2date_2(x)), round(a2 * y + b2, 1)) for x, y in vtx2],
        }

        # [MAU-1] IOS goes into BOTH the main quantitative CSV (rows_out)
        # AND the standalone archive (rows_out_ios, unchanged content/role
        # from before -- kept as a complete-series archive file per audit
        # instruction "ios.csvは残してよい").
        target_lists = [rows_out] if name != "IOS" else [rows_out, rows_out_ios]

        for x, y in sorted(cleaned1.items()):
            val = a1 * y + b1
            d = px2date_1(x)
            extrapolated = (x > X_LAST_ANCHOR_PX) or (x < X_FIRST_ANCHOR_PX)
            row = [str(d), name, "MAU", round(val, 1), SOURCE_IMAGE, 1,
                   round(xfit1["px_per_day"], 4), extrapolated]
            for tl in target_lists:
                tl.append(row)
        for x, y in sorted(cleaned2.items()):
            val = a2 * y + b2
            d = px2date_2(x)
            extrapolated = (x > X_LAST_ANCHOR_PX) or (x < X_FIRST_ANCHOR_PX)
            row = [str(d), name, "MAU", round(val, 1), SOURCE_IMAGE, 2,
                   round(xfit2["px_per_day"], 4), extrapolated]
            for tl in target_lists:
                tl.append(row)

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

        end_x1 = max(cleaned1) if cleaned1 else None
        end_x2 = max(cleaned2) if cleaned2 else None
        end_val1 = (a1 * cleaned1[end_x1] + b1) if end_x1 is not None else None
        end_val2 = (a2 * cleaned2[end_x2] + b2) if end_x2 is not None else None
        end_date1 = str(px2date_1(end_x1)) if end_x1 is not None else None
        end_date2 = str(px2date_2(end_x2)) if end_x2 is not None else None
        label_val = LABEL_VALUES[name]
        rel_err1 = (end_val1 - label_val) / label_val if end_val1 is not None else None
        rel_err2 = (end_val2 - label_val) / label_val if end_val2 is not None else None

        missing1 = find_missing_ranges(all_x_full, cleaned1.keys())
        missing2 = find_missing_ranges(all_x_full, cleaned2.keys())

        summary[name] = {
            "n_points_run1": len(cleaned1), "n_points_run2": len(cleaned2),
            "n_vertices_run1": len(vtx1), "n_vertices_run2": len(vtx2),
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

    csv_path = f"{OUT_DIR}/26-01-MAU_extracted.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["date", "platform", "metric", "value", "source_image",
                    "extraction_run", "px_per_day", "extrapolated"])
        for r in rows_out:
            w.writerow(r)

    ios_csv_path = f"{OUT_DIR}/26-01-MAU_ios.csv"
    with open(ios_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["date", "platform", "metric", "value", "source_image",
                    "extraction_run", "px_per_day", "extrapolated"])
        for r in rows_out_ios:
            w.writerow(r)

    with open(f"{OUT_DIR}/_summary_tmp_mau.json", "w", encoding="utf-8") as f:
        json.dump({
            "cal1_y_points": cal1["y_points"], "cal1_x_points": [(p, str(d)) for p, d in cal1["x_points"]],
            "cal2_y_points": cal2["y_points"], "cal2_x_points": [(p, str(d)) for p, d in cal2["x_points"]],
            "cal1_y_fit": [cal1["y_slope"], cal1["y_intercept"]],
            "cal2_y_fit": [cal2["y_slope"], cal2["y_intercept"]],
            "xfit1": xfit1, "xfit2": xfit2,
            "summary": summary,
            "cal1_log": cal1["log"], "cal2_log": cal2["log"],
            "ref_rgb_run2_derived": ref_rgb_run2_derived,
            "vertices": vertices_all,
        }, f, indent=2, default=str)

    print("DONE")
    print(f"CSV rows (main, all 4 series incl. IOS): {len(rows_out)}; "
          f"IOS archive rows (26-01-MAU_ios.csv, duplicate of IOS subset): {len(rows_out_ios)}")
    for name, s in summary.items():
        print(name, {k: v for k, v in s.items() if k not in
                      ("missing_ranges_run1_px", "missing_ranges_run2_px")})


if __name__ == "__main__":
    main()
