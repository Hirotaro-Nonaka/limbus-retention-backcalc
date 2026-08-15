"""
v2 Phase 2 Task (b): beta (playtime-boost) empirical identification.

r(t) = CCU(t)*24/(h*DAU_official(t)) - 1, regressed on core.py's frozen event-window
profile g(t): r(t) ~ a + beta*g(t)  (declared regression form, protocol frozen doc §3(b)).

Two h-handling variants (both mandatory per protocol):
  - h fixed at 2.2 (protocol-declared value): r(t) computed directly with h=2.2, then
    OLS r ~ a + b*g. Theoretical expectation: a=0 if h=2.2 is exactly correct; slope b = beta_hat.
  - h free: work with q(t) = CCU(t)*24/DAU_official(t) = h*(1+beta*g(t)). OLS q ~ a + b*g gives
    a = h_hat, b = h_hat*beta_hat  =>  beta_hat = b/a. Any h-mismatch is absorbed by the intercept.
  - h=2.128 (v1 fitted baseline value) also reported as a reference fixed-h case (not one of the
    two mandatory variants, but requested as a cross-check).

Frozen constraints (declared before running, not to be revisited post-hoc):
  - Event window / g(t): core.py's Fitter.g (max-clipped, exp3 unit-peak pulse, W=7 short/14 long).
  - Interval 1 (2023-08-01..2024-07-01) EXCLUDED from the event-window regression (date-uncertainty
    too large per digitization protocol changelog). Descriptive r(t) trend stats still computed for
    all 4 intervals, flagged accordingly.
  - DAU_official = STEAM run1 (primary). Reading-error floor 8.20% is propagated (see bootstrap note).
  - CCU = the same weekday-adjusted Steam CCU series core.py's Fitter is built on (f.y / adj), i.e.
    "v1と同じSteam CCU系列" per the frozen protocol. This is a forced-by-consistency choice (not
    post-hoc): using the frozen model's own CCU series keeps r(t) directly comparable to core.py's
    own beta*g(t) construction (Fitter uses adj as self.y, and g uses the same event calendar).

[TB-1 audit correction, declared BEFORE running per audit instruction]
The h=2.2-fixed regression residuals (sorted by date, intervals 2-4) show strong lag-1
autocorrelation (rho1=0.7095, Durbin-Watson=0.579, n_eff=387*(1-rho)/(1+rho)~=65.8), so the
original iid-row bootstrap CI understates uncertainty (audit estimate: ~2.4x too narrow). Two
block-bootstrap schemes are run and BOTH reported:
  (i)  EVENT-WINDOW BLOCKS (declared PRIMARY): each merged event window (g>0 span, extended
       backward to absorb the preceding non-event gap since the prior event's window end) forms
       one block. Rationale for primary status: the autocorrelation is structurally tied to the
       event-pulse process itself (g(t)); blocks that respect "one event cycle" as the natural
       correlation unit are the more faithful representation of the actual dependence structure,
       vs. an arbitrary fixed calendar length.
  (ii) FIXED 14-DAY CALENDAR BLOCKS (declared SECONDARY / robustness cross-check): block length
       chosen from the residual decorrelation scale (rho1=0.71 decaying with ~1.48-day average
       row spacing implies a decorrelation horizon of roughly 7-10 calendar days; 14 days is a
       conservative >=1x multiple of that horizon, and matches the length explicitly suggested in
       the audit instruction).
Both use block-level resampling with replacement (5000 reps), each draw also carries the same
independent 8.20%-sd reading-error perturbation used in the original (now superseded) bootstrap.
Point estimates (beta_hat=0.443 fixed-h / 0.541 free-h) are NOT changed by this correction --
only the propagated uncertainty (CI width) changes.
"""
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJ / 'scripts'))
import core  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent
DIGI_CSV = PROJ / 'v2' / 'digitized' / '26-01-DAU_extracted.csv'
RNG_SEED = 20260719

INTERVAL_BOUNDS = [
    ('2023-08-01', '2024-07-01', 1),
    ('2024-07-01', '2025-01-01', 2),
    ('2025-01-01', '2025-07-01', 3),
    ('2025-07-01', '2026-01-31', 4),
]
READ_ERR_STEAM_FLOOR = 0.0820


def assign_interval(d):
    for lo, hi, iid in INTERVAL_BOUNDS:
        if pd.Timestamp(lo) <= d < pd.Timestamp(hi):
            return iid
    return np.nan


def build_event_window_blocks(min_date, max_date):
    """Merged event-window blocks for scheme (i), PRIMARY per [TB-1] declaration.
    Each block = [end of previous block + 1 day, end of this (merged) event window],
    i.e. the event window plus all preceding non-event days back to the prior block's end.
    The final block absorbs any trailing non-event days up to max_date. Non-overlapping,
    contiguous, covers [min_date, max_date] exactly."""
    windows = []
    for name, (dstr, is_long) in core.LIMBUS_EVENTS.items():
        d0 = pd.Timestamp(dstr)
        if min_date <= d0 <= max_date:
            W = 14 if is_long else 7
            windows.append((d0, d0 + pd.Timedelta(days=W - 1)))
    for dstr in core.SPIKES:
        d0 = pd.Timestamp(dstr)
        if min_date <= d0 <= max_date:
            windows.append((d0, d0 + pd.Timedelta(days=6)))
    windows.sort()
    merged = []
    for s, e in windows:
        if merged and s <= merged[-1][1] + pd.Timedelta(days=1):
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    blocks = []
    prev_end = min_date - pd.Timedelta(days=1)
    for s, e in merged:
        blocks.append((prev_end + pd.Timedelta(days=1), e))
        prev_end = e
    if prev_end < max_date:
        blocks.append((prev_end + pd.Timedelta(days=1), max_date))
    return blocks


def build_fixed_length_blocks(min_date, max_date, length_days=14):
    """Fixed-length calendar blocks for scheme (ii), SECONDARY per [TB-1] declaration."""
    blocks = []
    d = min_date
    while d <= max_date:
        e = min(d + pd.Timedelta(days=length_days - 1), max_date)
        blocks.append((d, e))
        d = e + pd.Timedelta(days=1)
    return blocks


def assign_block_ids(dates, blocks):
    """dates: pd.Series/array of Timestamps. blocks: list of (start,end) contiguous, sorted.
    Returns integer block id per row (index into blocks)."""
    ids = np.full(len(dates), -1, dtype=int)
    dates = pd.DatetimeIndex(dates)
    for k, (s, e) in enumerate(blocks):
        mask = (dates >= s) & (dates <= e)
        ids[np.asarray(mask)] = k
    return ids


def block_bootstrap_beta(reg_data, block_ids, h_fixed, read_err_sd, n_boot, rng):
    """Resample blocks (by id) with replacement; within each draw also apply the same
    independent reading-error perturbation as the original iid bootstrap. Returns
    (boot_beta_fixed, boot_beta_free) arrays."""
    ccu_all = reg_data['ccu'].values.astype(float)
    dau_all = reg_data['official_dau'].values.astype(float)
    g_all = reg_data['g'].values.astype(float)
    n_blocks = block_ids.max() + 1
    row_idx_by_block = [np.where(block_ids == k)[0] for k in range(n_blocks)]

    boot_beta_fixed = np.empty(n_boot)
    boot_beta_free = np.empty(n_boot)
    for b in range(n_boot):
        chosen = rng.integers(0, n_blocks, size=n_blocks)
        idx = np.concatenate([row_idx_by_block[k] for k in chosen]) if n_blocks > 0 else np.array([], dtype=int)
        n_b = len(idx)
        delta = rng.normal(0.0, read_err_sd, size=n_b)
        dau_pert = dau_all[idx] * (1.0 + delta)
        ccu_b = ccu_all[idx]
        g_b = g_all[idx]
        r_b = ccu_b * 24.0 / (h_fixed * dau_pert) - 1.0
        _, slope_fixed, _, _, _ = ols_with_intercept(g_b, r_b)
        boot_beta_fixed[b] = slope_fixed
        q_b = ccu_b * 24.0 / dau_pert
        a_b, s_b, _, _, _ = ols_with_intercept(g_b, q_b)
        boot_beta_free[b] = s_b / a_b if a_b != 0 else np.nan
    return boot_beta_fixed, boot_beta_free


def ols_with_intercept(x, y):
    """returns (intercept, slope, se_intercept, se_slope, r2)"""
    X = np.column_stack([np.ones_like(x), x])
    beta_hat, resid, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta_hat
    n = len(y)
    dof = max(n - 2, 1)
    sigma2 = np.sum((y - yhat) ** 2) / dof
    XtX_inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(sigma2 * XtX_inv))
    ss_tot = np.sum((y - y.mean()) ** 2)
    ss_res = np.sum((y - yhat) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return beta_hat[0], beta_hat[1], se[0], se[1], r2


def main():
    # ---------- 1. Load frozen model g(t) and CCU series ----------
    f, adj = core.limbus_setup()  # f.g is the frozen event-window profile (core.py definition)
    g_series = pd.Series(f.g, index=adj.index, name='g')
    ccu_series = adj.copy()  # weekday-adjusted Steam CCU, same series Fitter.solve() uses as self.y
    ccu_series.name = 'ccu'

    print(f'[setup] g(t) built from {len(f.events)} frozen LIMBUS_EVENTS+SPIKES entries, '
          f'range [{f.g.min():.3f},{f.g.max():.3f}], CCU series {adj.index.min().date()}..{adj.index.max().date()}')

    # ---------- 2. Merge with digitized STEAM run1 (and run2 for reference) ----------
    digi = pd.read_csv(DIGI_CSV)
    digi['date'] = pd.to_datetime(digi['date'])

    all_rows = []
    for run in (1, 2):
        sub = digi[(digi['platform'] == 'STEAM') & (digi['extraction_run'] == run)].copy()
        sub = sub.sort_values('date')
        sub['ccu'] = sub['date'].map(ccu_series)
        sub['g'] = sub['date'].map(g_series)
        sub['extraction_run'] = run
        sub = sub.dropna(subset=['ccu', 'g'])
        sub['interval_id'] = sub['date'].apply(assign_interval)
        all_rows.append(sub[['date', 'extraction_run', 'value', 'ccu', 'g', 'interval_id']])
    merged = pd.concat(all_rows, ignore_index=True).rename(columns={'value': 'official_dau'})

    # r(t) for the two mandatory h variants + v1-fit reference h
    H_FIXED = 2.2
    H_V1FIT = 2.128  # reference only (v1 baseline fitted value, reproduced in task a: h=2.1281)
    merged['r_h2.2'] = merged['ccu'] * 24.0 / (H_FIXED * merged['official_dau']) - 1.0
    merged['r_h2.128'] = merged['ccu'] * 24.0 / (H_V1FIT * merged['official_dau']) - 1.0
    merged['q'] = merged['ccu'] * 24.0 / merged['official_dau']  # for free-h variant
    merged['event_flag'] = merged['g'] > 0

    merged.to_csv(OUT_DIR / 'r_series.csv', index=False)
    print(f'[merge] {len(merged)} matched rows (STEAM run1+run2, dates within CCU index range)')

    # ---------- 3. Descriptive r(t) trend by interval (ALL 4 intervals, run1 only, h=2.2) ----------
    r1 = merged[merged['extraction_run'] == 1].copy()
    print('[descriptive] r(t) [h=2.2 fixed] median by interval, STEAM run1 (interval 1 flagged low-confidence):')
    trend_rows = []
    for lo_d, hi_d, iid in INTERVAL_BOUNDS:
        sub = r1[r1['interval_id'] == iid]
        if len(sub) == 0:
            continue
        row = dict(interval=iid, date_range=f'{lo_d}..{hi_d}', n=len(sub),
                   median_r=sub['r_h2.2'].median(), mean_r=sub['r_h2.2'].mean(),
                   median_r_eventdays=sub.loc[sub['event_flag'], 'r_h2.2'].median() if sub['event_flag'].any() else np.nan,
                   median_r_noneventdays=sub.loc[~sub['event_flag'], 'r_h2.2'].median() if (~sub['event_flag']).any() else np.nan,
                   n_event_days=int(sub['event_flag'].sum()))
        trend_rows.append(row)
        print(f"  interval {iid} ({lo_d}..{hi_d}) n={len(sub)} median_r={row['median_r']:.4f} "
              f"mean_r={row['mean_r']:.4f} n_event_days={row['n_event_days']} "
              f"median_r|event={row['median_r_eventdays']} median_r|non-event={row['median_r_noneventdays']}")
    pd.DataFrame(trend_rows).to_csv(OUT_DIR / 'r_trend_by_interval.csv', index=False)

    # ---------- 4. Identifiability check: event-window coverage in intervals 2-4 ----------
    reg_data = r1[r1['interval_id'].isin([2, 3, 4])].copy().sort_values('date')
    reg_data = reg_data.reset_index(drop=True)

    # observed spacing of digitized STEAM run1 points within intervals 2-4
    gaps = reg_data['date'].diff().dt.days.dropna()
    print(f'[resolution] STEAM run1 intervals 2-4: n={len(reg_data)} points, '
          f'date-gap median={gaps.median():.2f}d mean={gaps.mean():.2f}d '
          f'(cf. protocol px/day~0.69 -> ~1.45 d/px)')

    # per-event coverage: how many digitized points fall inside each event's g>0 window
    idx_to_date = {i: d for d, i in {d: i for i, d in enumerate(adj.index)}.items()}
    date_to_idx = {d: i for i, d in enumerate(adj.index)}
    ev_cover = []
    for name, (dstr, is_long) in core.LIMBUS_EVENTS.items():
        d0 = pd.Timestamp(dstr)
        if not (pd.Timestamp('2024-07-01') <= d0 < pd.Timestamp('2026-01-31')):
            continue
        W = 14 if is_long else 7
        win_end = d0 + pd.Timedelta(days=W - 1)
        n_pts = int(((reg_data['date'] >= d0) & (reg_data['date'] <= win_end)).sum())
        ev_cover.append(dict(event=name, date=dstr, is_long=is_long, window_days=W, n_digitized_points=n_pts))
    for dstr in core.SPIKES:
        d0 = pd.Timestamp(dstr)
        if not (pd.Timestamp('2024-07-01') <= d0 < pd.Timestamp('2026-01-31')):
            continue
        W = 7
        win_end = d0 + pd.Timedelta(days=W - 1)
        n_pts = int(((reg_data['date'] >= d0) & (reg_data['date'] <= win_end)).sum())
        ev_cover.append(dict(event=f'SPIKE_{dstr}', date=dstr, is_long=False, window_days=W, n_digitized_points=n_pts))
    ev_cover_df = pd.DataFrame(ev_cover)
    ev_cover_df.to_csv(OUT_DIR / 'event_window_coverage.csv', index=False)
    cov_counts = ev_cover_df['n_digitized_points'].value_counts().sort_index()
    n_events_total = len(ev_cover_df)
    n_events_0pt = int((ev_cover_df['n_digitized_points'] == 0).sum())
    n_events_1pt = int((ev_cover_df['n_digitized_points'] == 1).sum())
    n_events_2plus = int((ev_cover_df['n_digitized_points'] >= 2).sum())
    print(f'[resolution] events in intervals 2-4: {n_events_total} total; '
          f'{n_events_0pt} with 0 digitized points in-window, {n_events_1pt} with exactly 1, '
          f'{n_events_2plus} with >=2 (needed to trace within-window decay shape)')
    print(f'[resolution] coverage distribution: {dict(cov_counts)}')

    n_event_rows = int(reg_data['event_flag'].sum())
    n_nonevent_rows = int((~reg_data['event_flag']).sum())
    print(f'[regression data] intervals 2-4: n={len(reg_data)} total, '
          f'{n_event_rows} on event-window days (g>0), {n_nonevent_rows} off-window (g=0)')

    # ---------- 5. Regression: h fixed at 2.2 ----------
    x = reg_data['g'].values.astype(float)
    y_r22 = reg_data['r_h2.2'].values.astype(float)
    a_r22, b_r22, se_a_r22, se_b_r22, r2_r22 = ols_with_intercept(x, y_r22)
    print(f'[h=2.2 fixed] intercept a={a_r22:.5f} (se={se_a_r22:.5f}), '
          f'beta_hat=slope b={b_r22:.5f} (se={se_b_r22:.5f}), R2={r2_r22:.4f}')

    # h=2.128 reference
    y_r_v1 = reg_data['r_h2.128'].values.astype(float)
    a_v1, b_v1, se_a_v1, se_b_v1, r2_v1 = ols_with_intercept(x, y_r_v1)
    print(f'[h=2.128 v1-fit reference] intercept a={a_v1:.5f} (se={se_a_v1:.5f}), '
          f'beta_hat=slope b={b_v1:.5f} (se={se_b_v1:.5f}), R2={r2_v1:.4f}')

    # ---------- 6. Regression: h free (q = h*(1+beta*g)) ----------
    y_q = reg_data['q'].values.astype(float)
    a_q, b_q, se_a_q, se_b_q, r2_q = ols_with_intercept(x, y_q)
    h_hat = a_q
    beta_hat_free = b_q / a_q if a_q != 0 else np.nan
    print(f'[h free] intercept a=h_hat={a_q:.5f} (se={se_a_q:.5f}), slope b={b_q:.5f} (se={se_b_q:.5f}), '
          f'R2={r2_q:.4f}  =>  beta_hat = b/a = {beta_hat_free:.5f}')

    # ---------- 7. [SUPERSEDED, not for publication] iid-row bootstrap ----------
    # Audit [TB-1]: residuals show strong lag-1 autocorrelation (rho1=0.7095, DW=0.579,
    # n_eff~=65.8 vs n=387), so this iid bootstrap understates uncertainty (~2.4x too narrow).
    # Kept only as a "pre-correction reference" value in the JSON/summary; NOT the reported CI.
    rng = np.random.default_rng(RNG_SEED)
    n_boot = 5000
    n = len(reg_data)
    ccu_arr = reg_data['ccu'].values.astype(float)
    dau_arr = reg_data['official_dau'].values.astype(float)
    g_arr = reg_data['g'].values.astype(float)

    boot_beta_fixed22 = np.empty(n_boot)
    boot_beta_free = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        # reading-error perturbation: independent multiplicative noise per resampled point,
        # sd = 8.20% (STEAM run1 floor). Treating as independent is a CONSERVATIVE choice: the
        # usage_conditions note the 8.20% floor includes a common-mode (axis-calibration) component
        # that would mostly shift the intercept and largely cancel out of the slope; modeling it as
        # fully independent per-point instead inflates the propagated slope (beta) uncertainty,
        # i.e. yields a wider, more conservative CI than the true common-mode-heavy error would.
        delta = rng.normal(0.0, READ_ERR_STEAM_FLOOR, size=n)
        dau_pert = dau_arr[idx] * (1.0 + delta)
        ccu_b = ccu_arr[idx]
        g_b = g_arr[idx]
        r_b = ccu_b * 24.0 / (H_FIXED * dau_pert) - 1.0
        _, slope_fixed, _, _, _ = ols_with_intercept(g_b, r_b)
        boot_beta_fixed22[b] = slope_fixed

        q_b = ccu_b * 24.0 / dau_pert
        a_b, s_b, _, _, _ = ols_with_intercept(g_b, q_b)
        boot_beta_free[b] = s_b / a_b if a_b != 0 else np.nan

    ci_fixed_iid_UNPUBLISHED = np.nanpercentile(boot_beta_fixed22, [2.5, 97.5])
    ci_free_iid_UNPUBLISHED = np.nanpercentile(boot_beta_free, [2.5, 97.5])
    print(f'[SUPERSEDED iid bootstrap n={n_boot}, NOT PUBLISHED, autocorrelation-naive] '
          f'h=2.2 fixed CI=[{ci_fixed_iid_UNPUBLISHED[0]:.5f}, {ci_fixed_iid_UNPUBLISHED[1]:.5f}]; '
          f'h free CI=[{ci_free_iid_UNPUBLISHED[0]:.5f}, {ci_free_iid_UNPUBLISHED[1]:.5f}]')

    # ---------- 7b. [TB-1] Block bootstrap (published CIs) ----------
    min_date, max_date = reg_data['date'].min(), reg_data['date'].max()
    dates_arr = reg_data['date'].values

    # (i) PRIMARY: event-window blocks
    ev_blocks = build_event_window_blocks(min_date, max_date)
    ev_block_ids = assign_block_ids(dates_arr, ev_blocks)
    rng_ev = np.random.default_rng(RNG_SEED + 1)
    boot_fixed_ev, boot_free_ev = block_bootstrap_beta(reg_data, ev_block_ids, H_FIXED,
                                                         READ_ERR_STEAM_FLOOR, n_boot, rng_ev)
    ci_fixed_ev = np.nanpercentile(boot_fixed_ev, [2.5, 97.5])
    ci_free_ev = np.nanpercentile(boot_free_ev, [2.5, 97.5])
    print(f'[block-bootstrap PRIMARY: event-window blocks, K={len(ev_blocks)}, n={n_boot}] '
          f'h=2.2 fixed: beta_hat={b_r22:.5f}, 95% CI=[{ci_fixed_ev[0]:.5f}, {ci_fixed_ev[1]:.5f}]')
    print(f'[block-bootstrap PRIMARY: event-window blocks] '
          f'h free: beta_hat={beta_hat_free:.5f}, 95% CI=[{ci_free_ev[0]:.5f}, {ci_free_ev[1]:.5f}]')

    # (ii) SECONDARY: fixed 14-day calendar blocks
    fx_blocks = build_fixed_length_blocks(min_date, max_date, length_days=14)
    fx_block_ids = assign_block_ids(dates_arr, fx_blocks)
    rng_fx = np.random.default_rng(RNG_SEED + 2)
    boot_fixed_fx, boot_free_fx = block_bootstrap_beta(reg_data, fx_block_ids, H_FIXED,
                                                         READ_ERR_STEAM_FLOOR, n_boot, rng_fx)
    ci_fixed_fx = np.nanpercentile(boot_fixed_fx, [2.5, 97.5])
    ci_free_fx = np.nanpercentile(boot_free_fx, [2.5, 97.5])
    print(f'[block-bootstrap SECONDARY: fixed 14-day blocks, K={len(fx_blocks)}, n={n_boot}] '
          f'h=2.2 fixed: beta_hat={b_r22:.5f}, 95% CI=[{ci_fixed_fx[0]:.5f}, {ci_fixed_fx[1]:.5f}]')
    print(f'[block-bootstrap SECONDARY: fixed 14-day blocks] '
          f'h free: beta_hat={beta_hat_free:.5f}, 95% CI=[{ci_free_fx[0]:.5f}, {ci_free_fx[1]:.5f}]')

    ci_fixed = ci_fixed_ev  # PUBLISHED primary CI
    ci_free = ci_free_ev    # PUBLISHED primary CI

    # ---------- 8. Identifiability verdict (based on PUBLISHED primary block-bootstrap CI) ----------
    ci_fixed_excludes_zero = not (ci_fixed[0] <= 0 <= ci_fixed[1])
    ci_free_excludes_zero = not (ci_free[0] <= 0 <= ci_free[1])
    ci_fixed_fx_excludes_zero = not (ci_fixed_fx[0] <= 0 <= ci_fixed_fx[1])
    ci_free_fx_excludes_zero = not (ci_free_fx[0] <= 0 <= ci_free_fx[1])
    identifiable = ci_fixed_excludes_zero and ci_free_excludes_zero
    print(f'[verdict] (primary, event-window blocks) h=2.2 CI excludes 0: {ci_fixed_excludes_zero}; '
          f'h free CI excludes 0: {ci_free_excludes_zero}')
    print(f'[verdict] (secondary, 14-day blocks) h=2.2 CI excludes 0: {ci_fixed_fx_excludes_zero}; '
          f'h free CI excludes 0: {ci_free_fx_excludes_zero}')
    print(f'[verdict] beta identifiable at this resolution (primary CI): {identifiable}')

    # ---------- 9. Compare to v1 D180-integrated band 3.0-6.6% ----------
    d180_band = (0.030, 0.066)
    print(f'[D180 band] v1 integrated band: {d180_band}')
    if identifiable:
        # naive translation note only -- beta here is an instantaneous CCU-boost fraction during
        # event windows, NOT directly the same units as v1's D180 playtime-contamination band
        # (which integrates a similar g(t) over the fit horizon). Report both point estimate & CI,
        # let the summary discuss the (lack of) direct unit correspondence. [TB-4]: this translation
        # is a rough guide only, not an independently-circulated result figure.
        print(f'  beta_hat (h=2.2) = {b_r22:.4f}, 95% CI (primary, event-window blocks) '
              f'[{ci_fixed[0]:.4f}, {ci_fixed[1]:.4f}]')
        print(f'  beta_hat (h free) = {beta_hat_free:.4f}, 95% CI (primary, event-window blocks) '
              f'[{ci_free[0]:.4f}, {ci_free[1]:.4f}]')

    # ---------- Summary JSON ----------
    summary = dict(
        n_matched_run1=int((merged['extraction_run'] == 1).sum()),
        n_matched_run2=int((merged['extraction_run'] == 2).sum()),
        r_trend_by_interval=trend_rows,
        resolution=dict(n_points_intervals2to4=len(reg_data),
                         gap_median_days=float(gaps.median()), gap_mean_days=float(gaps.mean()),
                         n_events_total=n_events_total, n_events_0pt=n_events_0pt,
                         n_events_1pt=n_events_1pt, n_events_2plus=n_events_2plus,
                         n_event_rows=n_event_rows, n_nonevent_rows=n_nonevent_rows),
        regression_h_fixed_2_2=dict(intercept=a_r22, se_intercept=se_a_r22,
                                     beta_hat=b_r22, se_beta=se_b_r22, r2=r2_r22,
                                     bootstrap_ci95_PUBLISHED_primary_eventblocks=list(ci_fixed),
                                     bootstrap_ci95_secondary_14day_blocks=list(ci_fixed_fx),
                                     bootstrap_ci95_iid_SUPERSEDED_unpublished=list(ci_fixed_iid_UNPUBLISHED)),
        regression_h_v1fit_2_128_reference=dict(intercept=a_v1, se_intercept=se_a_v1,
                                                  beta_hat=b_v1, se_beta=se_b_v1, r2=r2_v1),
        regression_h_free=dict(h_hat=a_q, se_h_hat=se_a_q, slope=b_q, se_slope=se_b_q, r2=r2_q,
                                beta_hat=beta_hat_free,
                                bootstrap_ci95_PUBLISHED_primary_eventblocks=list(ci_free),
                                bootstrap_ci95_secondary_14day_blocks=list(ci_free_fx),
                                bootstrap_ci95_iid_SUPERSEDED_unpublished=list(ci_free_iid_UNPUBLISHED)),
        residual_autocorrelation_diagnostic=dict(rho_lag1=0.709450403320763, durbin_watson=0.5794234033927281,
                                                   n=387, n_eff_approx=65.77710221743465),
        block_bootstrap_scheme=dict(primary='event_window_blocks', primary_n_blocks=len(ev_blocks),
                                     secondary='fixed_14day_blocks', secondary_n_blocks=len(fx_blocks)),
        identifiable=bool(identifiable),
        v1_d180_band=list(d180_band),
        reading_error_floor_steam=READ_ERR_STEAM_FLOOR,
        bootstrap_n=n_boot,
        bootstrap_seed=RNG_SEED,
    )
    with open(OUT_DIR / 'task_b_summary_raw.json', 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=2, default=str)
    print('[done] wrote r_series.csv, r_trend_by_interval.csv, event_window_coverage.csv, task_b_summary_raw.json')


if __name__ == '__main__':
    main()
