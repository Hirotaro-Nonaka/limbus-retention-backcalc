"""
v2 Phase 2 Task (a): 軌跡検証 (trajectory validation).

Compares the v1 frozen-baseline model's implied Steam DAU trajectory
(2023-08..2026-01) against the digitized official DAU chart series
(v2/digitized/26-01-DAU_extracted.csv).

Rules (frozen, see v2/Digitization_Protocol_frozen.md §3(a) and
v2/digitized/26-01-DAU_usage_conditions.md):
  - Primary series: STEAM run1. STEAM run2 reported as error-band reference (down-weighted).
  - TOTAL run1: shape-only comparison (level-matched), NOT used for component-sum checks.
  - iOS: excluded entirely (not present in this DAU csv anyway).
  - No BOUNDS changes, no re-tuning core.py. Baseline reproduction is verified first.
  - Relative RMSE definition matches core.py Fitter.report(): sqrt(mean(((model-official)/official)**2))
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJ / 'scripts'))
import core  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent
DIGI_CSV = PROJ / 'v2' / 'digitized' / '26-01-DAU_extracted.csv'

# Pre-declared interval boundaries (from the frozen protocol, matches the
# 5-point piecewise-linear breakpoints used in digitization calibration).
INTERVAL_BOUNDS = [
    ('2023-08-01', '2024-07-01', 1),  # low-resolution region, flagged separately
    ('2024-07-01', '2025-01-01', 2),
    ('2025-01-01', '2025-07-01', 3),
    ('2025-07-01', '2026-01-31', 4),  # includes trailing points up to 2026-01-25
]


def assign_interval(d: pd.Timestamp) -> int:
    for lo, hi, iid in INTERVAL_BOUNDS:
        if pd.Timestamp(lo) <= d < pd.Timestamp(hi):
            return iid
    return np.nan


def rel_rmse(model, official):
    model = np.asarray(model, dtype=float)
    official = np.asarray(official, dtype=float)
    return float(np.sqrt(np.mean(((model - official) / official) ** 2)))


def bias_pct(model, official):
    model = np.asarray(model, dtype=float)
    official = np.asarray(official, dtype=float)
    return float(np.mean((model - official) / official) * 100)


def main():
    # ---------- 1. Reproduce v1 frozen baseline fit ----------
    f, adj = core.limbus_setup()
    th, res = f.fit()
    out = f.report(th)
    baseline_rmse = out['rmse']
    print(f'[baseline] relative RMSE = {baseline_rmse:.6f} '
          f'(expected reproducibility range 0.1305-0.1314 per project memory)')
    print(f"[baseline] lam={out['lam']:.4f} k={out['k']:.4f} p={out['p']:.5f} "
          f"c={out['c']:.4f} tau={out['tau']:.2f} h={out['h']:.4f}")
    print(f"[baseline] D1={out['D1']:.4f} D7={out['D7']:.4f} D30={out['D30']:.4f} D180={out['D180']:.5f}")
    print(f"[baseline] base={out['base']:.6f} jan_dau={out.get('jan_dau'):.1f} jan_mau={out.get('jan_mau'):.1f}")

    # Boundary-sticking check (self-check requirement; report explicitly, do not alter BOUNDS)
    lo = np.array([b[0] for b in f.BOUNDS])
    hi = np.array([b[1] for b in f.BOUNDS])
    th_clip = np.clip(th, lo, hi)
    at_bound = np.isclose(th_clip, lo, rtol=1e-6, atol=1e-8) | np.isclose(th_clip, hi, rtol=1e-6, atol=1e-8)
    names = ['log(lam)', 'k', 'p', 'c', 'log(tau)', 'h']
    stuck = [n for n, s in zip(names, at_bound) if s]
    print(f'[baseline] parameters at BOUNDS edge: {stuck if stuck else "none"}')

    if not (0.129 <= baseline_rmse <= 0.133):
        print('[WARNING] baseline relative RMSE outside the expected 0.1305-0.1314 reproducibility '
              'band (looser check 0.129-0.133 also failed) -- STOP and investigate before proceeding.')
        sys.exit(1)

    # ---------- 2. Full model DAU trajectory (pre-CCU-conversion active-user count) ----------
    # report()['pred'] is CCU = modelDAU * h/24 * (1+beta*g). We need modelDAU = X @ a directly.
    lam, k, p, c, tau, h = out['lam'], out['k'], out['p'], out['c'], out['tau'], out['h']
    a_full = np.concatenate([[out['base']], out['amps']])
    R = core.R_kernel(f.T, lam, k, p, c, tau)
    X = core.build_X(f.T, f.events, R, f.shape, f.W[0], f.W[1])
    model_dau_full = X @ a_full  # length f.T, aligned to adj.index

    model_series = pd.Series(model_dau_full, index=adj.index, name='model_dau')
    model_series.to_csv(OUT_DIR / 'model_dau_full_trajectory.csv', header=True, index_label='date')

    # ---------- 3. Merge with digitized series (exact date match, no interpolation) ----------
    digi = pd.read_csv(DIGI_CSV)
    digi['date'] = pd.to_datetime(digi['date'])

    records = []
    for platform, run in [('STEAM', 1), ('STEAM', 2), ('TOTAL', 1)]:
        sub = digi[(digi['platform'] == platform) & (digi['extraction_run'] == run)].copy()
        sub = sub.sort_values('date')
        matched_model = []
        for d in sub['date']:
            if d in model_series.index:
                matched_model.append(model_series.loc[d])
            else:
                matched_model.append(np.nan)
        sub['model_dau'] = matched_model
        sub['interval_id'] = sub['date'].apply(assign_interval)
        sub = sub.dropna(subset=['model_dau'])
        records.append(sub[['date', 'platform', 'extraction_run', 'value', 'model_dau', 'interval_id']])

    comp = pd.concat(records, ignore_index=True)
    comp = comp.rename(columns={'value': 'official_dau'})
    comp.to_csv(OUT_DIR / 'trajectory_comparison.csv', index=False)
    print(f'[merge] {len(comp)} matched rows across STEAM run1/run2, TOTAL run1 '
          f'(unmatched dates outside model date-index dropped: '
          f"{sum(len(digi[(digi['platform']==p)&(digi['extraction_run']==r)]) for p,r in [('STEAM',1),('STEAM',2),('TOTAL',1)]) - len(comp)})")

    # ---------- 4. STEAM run1: primary metrics ----------
    s1 = comp[(comp['platform'] == 'STEAM') & (comp['extraction_run'] == 1)]
    overall_rmse = rel_rmse(s1['model_dau'], s1['official_dau'])
    overall_bias = bias_pct(s1['model_dau'], s1['official_dau'])
    print(f'[STEAM run1] n={len(s1)} overall relative RMSE={overall_rmse:.6f} bias%={overall_bias:.3f}')

    interval_rows = []
    for lo_d, hi_d, iid in INTERVAL_BOUNDS:
        sub = s1[s1['interval_id'] == iid]
        if len(sub) == 0:
            continue
        r = rel_rmse(sub['model_dau'], sub['official_dau'])
        b = bias_pct(sub['model_dau'], sub['official_dau'])
        interval_rows.append(dict(interval=iid, date_range=f'{lo_d}..{hi_d}', n=len(sub),
                                   rel_rmse=r, bias_pct=b))
        print(f'[STEAM run1] interval {iid} ({lo_d}..{hi_d}) n={len(sub)} '
              f'rel_rmse={r:.6f} bias%={b:.3f}')

    pd.DataFrame(interval_rows).to_csv(OUT_DIR / 'interval_breakdown_steam_run1.csv', index=False)

    # ---------- 5. STEAM run2 (reference / error-band, down-weighted) ----------
    s2 = comp[(comp['platform'] == 'STEAM') & (comp['extraction_run'] == 2)]
    run2_rows = []
    if len(s2) > 0:
        r2_overall = rel_rmse(s2['model_dau'], s2['official_dau'])
        b2_overall = bias_pct(s2['model_dau'], s2['official_dau'])
        print(f'[STEAM run2] n={len(s2)} overall relative RMSE={r2_overall:.6f} bias%={b2_overall:.3f}')
        for lo_d, hi_d, iid in INTERVAL_BOUNDS:
            sub = s2[s2['interval_id'] == iid]
            if len(sub) == 0:
                continue
            r = rel_rmse(sub['model_dau'], sub['official_dau'])
            b = bias_pct(sub['model_dau'], sub['official_dau'])
            run2_rows.append(dict(interval=iid, date_range=f'{lo_d}..{hi_d}', n=len(sub),
                                   rel_rmse=r, bias_pct=b))
            print(f'[STEAM run2] interval {iid} n={len(sub)} rel_rmse={r:.6f} bias%={b:.3f}')
    pd.DataFrame(run2_rows).to_csv(OUT_DIR / 'interval_breakdown_steam_run2.csv', index=False)

    # ---------- 6. TOTAL run1: shape-only comparison (level-matched) ----------
    t1 = comp[(comp['platform'] == 'TOTAL') & (comp['extraction_run'] == 1)].copy()
    total_summary = {}
    if len(t1) > 0:
        # Level-match: scale model trajectory by ratio of means over the matched dates (NOT a
        # component-sum check -- TOTAL is a distinct, non-additive official series; here we only
        # ask whether the *shape* (trend) of the model-Steam trajectory tracks TOTAL's shape).
        scale = t1['official_dau'].mean() / t1['model_dau'].mean()
        t1['model_dau_scaled'] = t1['model_dau'] * scale
        shape_rmse = rel_rmse(t1['model_dau_scaled'], t1['official_dau'])
        raw_level_ratio = (t1['model_dau'].mean() / t1['official_dau'].mean())
        total_summary = dict(n=len(t1), scale_factor=scale,
                              raw_model_to_official_level_ratio=raw_level_ratio,
                              shape_rmse_after_levelmatch=shape_rmse)
        print(f"[TOTAL run1 shape-only] n={len(t1)} scale_factor={scale:.4f} "
              f"raw_level_ratio(model/official)={raw_level_ratio:.4f} "
              f"shape_rmse_after_levelmatch={shape_rmse:.6f}")
        t1.to_csv(OUT_DIR / 'total_run1_shape_comparison.csv', index=False)

    # ---------- 7. Time-pattern: early vs late half within STEAM run1 (excluding interval 1) ----------
    s1_reliable = s1[s1['interval_id'] != 1].sort_values('date')
    if len(s1_reliable) > 4:
        mid = len(s1_reliable) // 2
        early = s1_reliable.iloc[:mid]
        late = s1_reliable.iloc[mid:]
        early_bias = bias_pct(early['model_dau'], early['official_dau'])
        late_bias = bias_pct(late['model_dau'], late['official_dau'])
        print(f'[time pattern, intervals 2-4] early-half bias%={early_bias:.3f} '
              f'(n={len(early)}, {early["date"].min().date()}..{early["date"].max().date()}); '
              f'late-half bias%={late_bias:.3f} (n={len(late)}, '
              f'{late["date"].min().date()}..{late["date"].max().date()})')
    else:
        early_bias = late_bias = np.nan

    # ---------- 8. Reading-error propagation ----------
    READ_ERR_STEAM_FLOOR = 0.0820  # from usage_conditions.md, floor, near-term biased
    exceeds = abs(overall_bias) / 100.0 > READ_ERR_STEAM_FLOOR
    print(f'[reading-error check] STEAM run1 |overall bias%|={abs(overall_bias):.3f}% vs '
          f'reading-error floor={READ_ERR_STEAM_FLOOR*100:.2f}% -> '
          f'{"EXCEEDS reading-error floor" if exceeds else "within reading-error floor"}')
    print(f'[reading-error check] STEAM run1 overall relative RMSE={overall_rmse:.6f} '
          f'({overall_rmse*100:.3f}%) vs reading-error floor {READ_ERR_STEAM_FLOOR*100:.2f}%')

    # ---------- Summary ----------
    summary = dict(
        baseline_rmse=baseline_rmse,
        baseline_params=dict(lam=lam, k=k, p=p, c=c, tau=tau, h=h),
        params_at_bounds=stuck,
        v1_holdout_dau_error_pct=-5.4,
        steam_run1_overall_rel_rmse=overall_rmse,
        steam_run1_overall_bias_pct=overall_bias,
        steam_run1_intervals=interval_rows,
        steam_run2_overall_rel_rmse=(r2_overall if len(s2) > 0 else None),
        steam_run2_overall_bias_pct=(b2_overall if len(s2) > 0 else None),
        steam_run2_intervals=run2_rows,
        total_run1_shape=total_summary,
        early_half_bias_pct_intervals2to4=early_bias,
        late_half_bias_pct_intervals2to4=late_bias,
        reading_error_floor_steam=READ_ERR_STEAM_FLOOR,
        exceeds_reading_error_floor=bool(exceeds),
    )
    import json
    with open(OUT_DIR / 'task_a_summary_raw.json', 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=2, default=str)
    print('[done] wrote trajectory_comparison.csv, interval breakdown CSVs, '
          'total_run1_shape_comparison.csv, model_dau_full_trajectory.csv, task_a_summary_raw.json')


if __name__ == '__main__':
    main()
