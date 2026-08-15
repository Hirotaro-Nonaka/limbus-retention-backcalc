"""
v2 Phase 3 Task (c): anchor marginal-value experiment.

Frozen protocol: v2/Digitization_Protocol_frozen.md sec.0(c) and sec.3(c).
Anchors are added in fixed chronological order (24-10 -> 25-01 -> 25-05 -> 25-10 -> 26-01),
stage k uses the first k anchors. CCU data, event calendar (core.LIMBUS_EVENTS + core.SPIKES)
and BOUNDS (core.Fitter.BOUNDS) are held FIXED across all stages (v1-frozen; not touched here).
Only the anchor constraint set changes between stages -- this isolates the marginal identification
value of each additional disclosed anchor point, per the declared experiment design.

---------------------------------------------------------------------------------------------
[FORCED EXTENSION, logged per the frozen protocol -- "data-constraint-forced", NOT a
post-hoc method change]: core.py's Fitter class supports exactly ONE anchor point (a single
jan_mask/dau/mau/ratio in its `anchors` dict). This task requires up to FIVE simultaneous anchor
points (stage 5). core.py is FROZEN and must not be edited (per the frozen protocol), so this
script defines MultiAnchorFitter(core.Fitter), a subclass that generalizes solve()/report() to
loop over a LIST of anchor dicts instead of a single one, while reusing core.py's kernels
(R_kernel, A_kernel), design-matrix builders (build_X, build_X_kernel), event calendar, BOUNDS,
and the two-stage Powell->Nelder-Mead optimizer (via the inherited .fit()/.obj_b()/.obj()).
The per-anchor weighting convention is an EXACT replication of core.py's single-anchor convention,
applied independently to each anchor:
  - DAU constraint row weighted at 5% relative tolerance (matches core.py's 0.05)
  - MAU constraint row weighted at 3% relative tolerance (matches core.py's 0.03)
  - ratio penalty term ((dau/mau - r0)/rs)^2 with rs=0.02 (matches core.py's ratio width 0.02);
    r0 is computed as that anchor's OWN dau/mau ratio (e.g. anchor 26-01: 360639/792250=0.4552,
    which is exactly how core.py's hardcoded r0=0.455 for the single 26-01 anchor was derived --
    so this is a literal generalization of the same rule, not a new rule).
This is the ONLY way to run the declared 1/2/3/4/5-anchor stage sequence at all; without it the
experiment (the entire point of task c) cannot be executed. Decided BEFORE any stage was fit.
---------------------------------------------------------------------------------------------

Anchor date convention [logged, decided before fitting]: for each broadcast, the anchor's
jan_mask uses the CALENDAR MONTH NAMED BY THE BROADCAST LABEL (e.g. "24-10" -> October 2024,
"26-01" -> January 2026), exactly mirroring core.py's own convention for the 26-01 anchor
(jan_mask = 2026-01-01..2026-01-31, matching the "26-01" label). This is applied uniformly to
all 5 broadcasts rather than using each chart's precise terminal x-axis date, because:
  - the exact terminal date is unambiguous only for 25-05 (axis explicitly labeled "2025.05.31")
    and 25-10 (axis explicitly labeled "2025.11.09", i.e. the "25-10" broadcast's data actually
    extends into early November -- see broadcast_anchors.md note);
  - for 24-10 and 25-01 no explicit terminal date is printed (only the tick "2024.10.01" /
    "2025.01.01"; the curve visibly extends a few pixels past the tick in both charts).
  Using the uniform "label month" rule keeps the anchor convention IDENTICAL in kind across all
  5 stages (reproducing exactly what core.py already does for 26-01) instead of mixing exact
  dates (where available) with inferred ones (where not) -- which would make the marginal-value
  comparison across stages confounded by a convention change rather than by anchor count alone.
  Consequence: the 25-10 anchor's mask is October 2025 (NOT extended to include the early-November
  tail visible on that chart); this is a deliberate, documented choice, not an oversight.
"""
import sys
import time
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize, nnls

PROJ = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJ / 'scripts'))
import core  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent
RNG_SEED = 20260720
N_BOOT = 100          # matches v1 boot.py: moving-block bootstrap, 100 replicates
BLOCK = 14            # matches v1 boot.py: block size B=14
PARAM_NAMES = ['log(lam)', 'k', 'p', 'c', 'log(tau)', 'h']

# ---------------------------------------------------------------------------------------------
# Anchor specs, chronological order fixed by protocol sec.3(c). Values from
# v2/task_c/broadcast_anchors.md, independently re-verified against docs/*.png (see
# v2/task_c/anchor_verification.md) -- all 8 charts (24-10/25-01/25-05/25-10 x DAU/MAU) matched
# exactly, no discrepancies.
# ---------------------------------------------------------------------------------------------
ANCHOR_SPECS = [
    dict(label='24-10', month_start='2024-10-01', month_end='2024-10-31', dau=175024.0, mau=327107.0),
    dict(label='25-01', month_start='2025-01-01', month_end='2025-01-31', dau=185661.0, mau=408939.0),
    dict(label='25-05', month_start='2025-05-01', month_end='2025-05-31', dau=279651.0, mau=520084.0),
    dict(label='25-10', month_start='2025-10-01', month_end='2025-10-31', dau=418242.0, mau=728973.0),
    dict(label='26-01', month_start='2026-01-01', month_end='2026-01-31', dau=360639.0, mau=792250.0),
]


# ---------------------------------------------------------------------------------------------
# Data / event calendar loading -- IDENTICAL reproduction of core.limbus_setup()'s internal
# logic, minus the single hardcoded 'anchors' dict (which we build per-stage below instead).
# CCU data and event calendar are therefore byte-identical to what core.py itself would use.
# ---------------------------------------------------------------------------------------------
def load_series_and_events():
    s = core.load_steamdb_daily(core.DATA_DIR / 'Limbus_steamdb_chart_1973530.csv')
    adj, _ = core.weekday_adjust(s)
    idx = {d: i for i, d in enumerate(adj.index)}
    events = []
    for name, (d, lng) in core.LIMBUS_EVENTS.items():
        events.append((idx[pd.Timestamp(d)], lng))
    for d in core.SPIKES:
        events.append((idx[pd.Timestamp(d)], False))
    return adj, events


def build_anchor_list(adj, k):
    """First k anchors (chronological), each with its boolean month-mask over adj.index."""
    out = []
    for spec in ANCHOR_SPECS[:k]:
        mask = np.asarray((adj.index >= spec['month_start']) & (adj.index <= spec['month_end']))
        out.append(dict(label=spec['label'], mask=mask, dau=spec['dau'], mau=spec['mau'],
                         n_days_in_mask=int(mask.sum())))
    return out


# ---------------------------------------------------------------------------------------------
# MultiAnchorFitter: generalizes core.Fitter to N simultaneous anchors (see module docstring).
# Inherits BOUNDS, unpack(), obj(), obj_b(), fit() unchanged from core.Fitter.
# ---------------------------------------------------------------------------------------------
class MultiAnchorFitter(core.Fitter):
    def __init__(self, y, events, anchors_list, **kwargs):
        kwargs.setdefault('anchors', None)  # disable base class's single-anchor path entirely
        super().__init__(y, events, **kwargs)
        self.anchors_list = anchors_list or []

    def solve(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        R = core.R_kernel(self.T, lam, k, p, c, tau)
        X = core.build_X(self.T, self.events, R, self.shape, self.W[0], self.W[1])
        mult = (h / 24.0) * (1.0 + self.beta * self.g)
        rows = [X * (mult / self.y)[:, None]]
        rhs = [np.ones(self.T)]
        A = Xm = None
        if self.anchors_list:
            A = core.A_kernel(self.T, lam, k, p, tau)
            Xm = core.build_X_kernel(self.T, self.events, A, self.shape, self.W[0], self.W[1])
            for anc in self.anchors_list:
                m = anc['mask']
                dau_row = X[m].mean(axis=0)
                mau_row = Xm[m].mean(axis=0)
                rows.append((dau_row / (0.05 * anc['dau']))[None, :])
                rhs.append(np.array([1 / 0.05]))
                rows.append((mau_row / (0.03 * anc['mau']))[None, :])
                rhs.append(np.array([1 / 0.03]))
        Xw = np.vstack(rows)
        yw = np.concatenate(rhs)
        a, rnorm = nnls(Xw, yw, maxiter=10 * Xw.shape[1])
        loss = rnorm ** 2
        if self.anchors_list:
            for anc in self.anchors_list:
                m = anc['mask']
                dau = X[m].mean(axis=0) @ a
                mau = Xm[m].mean(axis=0) @ a
                r0 = anc['dau'] / anc['mau']
                rs = 0.02
                loss += ((dau / mau - r0) / rs) ** 2
        if self.use_h:
            loss += ((h - self.h_prior[0]) / self.h_prior[1]) ** 2
        return loss, a

    def report(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        loss, a = self.solve(th)
        R = core.R_kernel(max(self.T, 400), lam, k, p, c, tau)
        X = core.build_X(self.T, self.events, core.R_kernel(self.T, lam, k, p, c, tau),
                          self.shape, self.W[0], self.W[1])
        mult = (h / 24.0) * (1.0 + self.beta * self.g)
        pred = (X @ a) * mult
        rmse = np.sqrt(np.mean(((pred - self.y) / self.y) ** 2))
        out = dict(D1=R[1], D7=R[7], D30=R[30], D180=R[180],
                   D365=R[365] if len(R) > 365 else np.nan,
                   lam=lam, k=k, p=p, c=c, tau=tau, h=h, rmse=rmse, loss=loss,
                   base=a[0], amps=a[1:], pred=pred)
        if self.anchors_list:
            A = core.A_kernel(self.T, lam, k, p, tau)
            Xm = core.build_X_kernel(self.T, self.events, A, self.shape, self.W[0], self.W[1])
            per_anchor = []
            for anc in self.anchors_list:
                m = anc['mask']
                pred_dau = float(X[m].mean(axis=0) @ a)
                pred_mau = float(Xm[m].mean(axis=0) @ a)
                per_anchor.append(dict(label=anc['label'], pred_dau=pred_dau, pred_mau=pred_mau,
                                        target_dau=anc['dau'], target_mau=anc['mau']))
            out['per_anchor'] = per_anchor
        return out


# ---------------------------------------------------------------------------------------------
# Multistart design: same 5 initial points used for EVERY stage (documented, matches v1's
# convention of a fixed th0 set spanning: core default (scripts/task_b.py DEFAULT_TH0), two
# beta-sweep near-optimum starts (task_b.py extra_02/extra_03), and two of task_c.py (7-title)
# th0_original_region() spread points). Two-stage Powell->Nelder-Mead per start, inherited
# unchanged from core.Fitter.fit().
# ---------------------------------------------------------------------------------------------
def th0_list():
    return [
        np.array([np.log(1.5), 0.30, 0.023, 0.90, np.log(6000.0), 2.2]),   # core default
        np.array([np.log(1.6), 0.35, 0.05, 0.85, np.log(2000.0), 2.14]),   # task_b extra_02
        np.array([np.log(1.3), 0.25, 0.10, 0.80, np.log(1500.0), 2.13]),   # task_b extra_03
        np.array([np.log(20.0), 1.0, 0.05, 0.9, np.log(1000.0), 2.2]),     # task_c th0_original_region[0]
        np.array([np.log(5.0), 0.5, 0.10, 0.5, np.log(3000.0), 2.2]),      # task_c th0_original_region[2]
    ]


def multistart_fit(f, th0s):
    lo = np.array([b[0] for b in f.BOUNDS])
    hi = np.array([b[1] for b in f.BOUNDS])
    results = []
    for th0 in th0s:
        th, res = f.fit(th0=np.asarray(th0, dtype=float))
        obj = f.obj_b(th)
        results.append((th, obj, np.asarray(th0, dtype=float)))
    results.sort(key=lambda t: t[1])
    return results


def bounds_stuck(th, f):
    lo = np.array([b[0] for b in f.BOUNDS])
    hi = np.array([b[1] for b in f.BOUNDS])
    thc = np.clip(th, lo, hi)
    at_lo = np.isclose(thc, lo, rtol=1e-6, atol=1e-8)
    at_hi = np.isclose(thc, hi, rtol=1e-6, atol=1e-8)
    stuck = []
    for i, n in enumerate(PARAM_NAMES):
        if at_lo[i]:
            stuck.append((n, 'lower', lo[i]))
        elif at_hi[i]:
            stuck.append((n, 'upper', hi[i]))
    return stuck


# ---------------------------------------------------------------------------------------------
# Bootstrap: moving-block bootstrap of RELATIVE RESIDUALS of the CCU fit (matches v1's
# docs/pastanalysis/from_other_chat/boot.py design exactly in kind: block=14, N=100, refit via
# a LIGHT Nelder-Mead-only local step from the best theta + small perturbation, NOT a full
# multistart -- v1's boot.py used exactly this "light refit" pattern for tractability).
# Anchor DAU/MAU targets are NOT resampled (v1's boot.py never resampled the anchor constraint
# either -- only the CCU series residuals are bootstrapped).
# ---------------------------------------------------------------------------------------------
def bootstrap_stage(f, th_best, rng, n_boot=N_BOOT, block=BLOCK):
    rep = f.report(th_best)
    pred = rep['pred']
    y = f.y
    resid = (y - pred) / pred  # relative residual: y = pred*(1+resid)
    n = len(y)
    lo = np.array([b[0] for b in f.BOUNDS])
    hi = np.array([b[1] for b in f.BOUNDS])

    rows = []
    for b in range(n_boot):
        starts = rng.integers(0, max(n - block, 1), size=n // block + 1)
        idx = np.concatenate([np.arange(s, s + block) for s in starts])[:n]
        idx = np.clip(idx, 0, n - 1)
        resid_b = resid[idx]
        y_b = np.maximum(pred * (1 + resid_b), 100.0)

        fb = MultiAnchorFitter(y_b, f.events, f.anchors_list, shape=f.shape, beta=f.beta,
                                h_prior=f.h_prior, use_h=f.use_h, W=f.W)
        th0_b = np.clip(th_best + rng.normal(0, 0.05, 6), lo, hi)
        res = minimize(fb.obj_b, th0_b, method='Nelder-Mead',
                        options=dict(maxiter=250, xatol=2e-3, fatol=1e-6))
        th_b = np.clip(res.x, lo, hi)
        rep_b = fb.report(th_b)
        rows.append(dict(D1=rep_b['D1'], D7=rep_b['D7'], D30=rep_b['D30'], D180=rep_b['D180'],
                          lam=rep_b['lam'], k=rep_b['k'], p=rep_b['p'], c=rep_b['c'],
                          tau=rep_b['tau'], h=rep_b['h'], base=rep_b['base'], rmse=rep_b['rmse']))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------
# Stage runner
# ---------------------------------------------------------------------------------------------
def run_stage(k, adj, events, rng):
    t0 = time.time()
    anchors_list = build_anchor_list(adj, k)
    y = adj.values.astype(float)
    f = MultiAnchorFitter(y, events, anchors_list, shape='exp3', beta=0.0,
                           h_prior=(2.2, 0.2), use_h=True, W=(7, 14))

    results = multistart_fit(f, th0_list())
    best_th, best_obj, best_th0 = results[0]
    rep = f.report(best_th)
    stuck = bounds_stuck(best_th, f)

    d180_restart = [f.report(th)['D180'] for th, _, _ in results]
    obj_restart = [o for _, o, _ in results]

    print(f'[stage {k}] anchors={[a["label"] for a in anchors_list]} '
          f'n_days_in_mask={[a["n_days_in_mask"] for a in anchors_list]} '
          f'best_obj={best_obj:.4f} restart_obj_range=[{min(obj_restart):.4f},{max(obj_restart):.4f}] '
          f'D180_restart_range=[{min(d180_restart):.5f},{max(d180_restart):.5f}] '
          f'stuck={stuck} fit_time={time.time()-t0:.1f}s', flush=True)

    t1 = time.time()
    boot_df = bootstrap_stage(f, best_th, rng)
    boot_df.to_csv(OUT_DIR / f'bootstrap_stage{k}.csv', index=False)
    print(f'[stage {k}] bootstrap n={len(boot_df)} time={time.time()-t1:.1f}s', flush=True)

    d180_ci_lo, d180_ci_hi = np.percentile(boot_df['D180'], [5, 95])
    tau_lo_boot = np.percentile(boot_df['tau'], 5)
    base_boot_lo, base_boot_hi = np.percentile(boot_df['base'], [5, 95])
    tau_lo_bound = np.exp(f.BOUNDS[4][0])
    tau_hi_bound = np.exp(f.BOUNDS[4][1])

    row = dict(
        stage=k, n_anchors=k, anchors=','.join(a['label'] for a in anchors_list),
        lam=rep['lam'], k_param=rep['k'], p=rep['p'], c=rep['c'], tau=rep['tau'], h=rep['h'],
        base=rep['base'], D1=rep['D1'], D7=rep['D7'], D30=rep['D30'], D180=rep['D180'],
        rmse=rep['rmse'], obj_best=best_obj,
        n_bounds_stuck=len(stuck), bounds_stuck_detail=';'.join(f'{n}@{side}' for n, side, _ in stuck),
        n_restarts=len(results),
        obj_restart_lo=min(obj_restart), obj_restart_hi=max(obj_restart),
        D180_restart_lo=min(d180_restart), D180_restart_hi=max(d180_restart),
        D180_boot_median=float(np.median(boot_df['D180'])),
        D180_ci90_lo=float(d180_ci_lo), D180_ci90_hi=float(d180_ci_hi),
        D180_ci90_width=float(d180_ci_hi - d180_ci_lo),
        tau_boot_p5=float(tau_lo_boot), tau_bound_lo=tau_lo_bound, tau_bound_hi=tau_hi_bound,
        tau_boot_at_lower_bound=bool(np.isclose(tau_lo_boot, tau_lo_bound, rtol=1e-3)),
        base_boot_p5=float(base_boot_lo), base_boot_p95=float(base_boot_hi),
        base_degenerate=bool(np.isclose(rep['base'], 0.0, atol=1e-6)),
        n_boot=len(boot_df), block=BLOCK,
        fit_time_s=time.time() - t0,
    )
    for anc, per in zip(anchors_list, rep['per_anchor']):
        pass  # per-anchor fit residuals already embedded in rep['per_anchor']; not part of the
              # 4 declared quantities, kept out of stage_results.csv row (recorded separately below)

    return row, results, rep, anchors_list


def _done_stages():
    """Resumability check (added after a 2026-07-20 crash mid-stage-5-bootstrap; the background
    process was externally KILLED by the sandbox -- no Python exception/traceback in run_log.txt,
    no memory-error signature -- most likely a background-process wall-clock cap in this
    environment, not a bug in the fit/bootstrap logic itself, since stages 1-4's bootstraps of
    identical size/design completed cleanly). A stage counts as 'done' only if BOTH its
    bootstrap_stage{k}.csv exists on disk AND stage_results.csv already has a row for it (defends
    against a bootstrap CSV left over from a stage whose stage_results row never got written,
    which is exactly the stage-5 failure mode we hit)."""
    csv_path = OUT_DIR / 'stage_results.csv'
    if not csv_path.exists():
        return set(), [], [], []
    existing = pd.read_csv(csv_path)
    rows = existing.to_dict('records')
    restart_rows = (pd.read_csv(OUT_DIR / 'restart_detail.csv').to_dict('records')
                     if (OUT_DIR / 'restart_detail.csv').exists() else [])
    per_anchor_rows = (pd.read_csv(OUT_DIR / 'per_anchor_fit_check.csv').to_dict('records')
                        if (OUT_DIR / 'per_anchor_fit_check.csv').exists() else [])
    done = {int(r['stage']) for r in rows if (OUT_DIR / f"bootstrap_stage{int(r['stage'])}.csv").exists()}
    return done, rows, restart_rows, per_anchor_rows


def main():
    adj, events = load_series_and_events()
    print(f'[setup] T={len(adj)} {adj.index.min().date()}..{adj.index.max().date()}, '
          f'n_events={len(events)}', flush=True)

    done_stages, rows, restart_rows, per_anchor_rows = _done_stages()
    if done_stages:
        print(f'[resume] stages already complete on disk (bootstrap CSV + stage_results row '
              f'both present): {sorted(done_stages)} -- reusing, not recomputing.', flush=True)

    # Seed-handling note [logged, per coordinator instruction after the 2026-07-20 crash]:
    # the ORIGINAL (pre-crash) run used ONE np.random.default_rng(RNG_SEED) instance threaded
    # sequentially through stages 1..5 (stage k's bootstrap draws are whatever the shared stream
    # was at when stage k started). That run completed stages 1-4's bootstraps intact (kept as-is,
    # NOT recomputed) before being killed mid-stage-5-bootstrap, so the shared stream's exact state
    # at the start of stage 5 cannot be recovered (it lived only in the killed process's memory).
    # Any stage re-run after a crash therefore uses its OWN independent seed
    # default_rng(RNG_SEED + 5000 + k) instead of trying to fake a continuation of the lost shared
    # stream. This is the SAME bootstrap design in kind (moving-block, block=14, N=100,
    # np.random.default_rng) as stages 1-4 -- only the seed-stream provenance differs, which does
    # not affect the validity of a given stage's own CI (each stage's bootstrap is a self-contained
    # resampling procedure). This is a data/infra-constraint-forced change (crash recovery), not a
    # post-hoc result-driven choice: decided before looking at any stage-5 bootstrap numbers.
    for k in range(1, 6):
        if k in done_stages:
            continue
        is_resume = len(done_stages) > 0
        rng = np.random.default_rng(RNG_SEED + 5000 + k) if is_resume else np.random.default_rng(RNG_SEED)
        if is_resume:
            print(f'[resume] stage {k}: independent rerun seed = {RNG_SEED + 5000 + k} '
                  f'(see seed-handling note above)', flush=True)
        row, results, rep, anchors_list = run_stage(k, adj, events, rng)
        rows = [r for r in rows if int(r['stage']) != k] + [row]
        rows.sort(key=lambda r: int(r['stage']))
        restart_rows = [r for r in restart_rows if int(r['stage']) != k]
        for i, (th, obj, th0) in enumerate(results):
            restart_rows.append(dict(stage=k, restart_idx=i, obj=obj,
                                      log_lam=th[0], k_param=th[1], p=th[2], c=th[3],
                                      log_tau=th[4], h=th[5],
                                      th0_log_lam=th0[0], th0_k=th0[1], th0_p=th0[2],
                                      th0_c=th0[3], th0_log_tau=th0[4], th0_h=th0[5]))
        per_anchor_rows = [r for r in per_anchor_rows if int(r['stage']) != k]
        for per in rep['per_anchor']:
            per_anchor_rows.append(dict(stage=k, **per))
        # incremental save after every stage (session-limit precaution)
        pd.DataFrame(rows).to_csv(OUT_DIR / 'stage_results.csv', index=False)
        pd.DataFrame(restart_rows).to_csv(OUT_DIR / 'restart_detail.csv', index=False)
        pd.DataFrame(per_anchor_rows).to_csv(OUT_DIR / 'per_anchor_fit_check.csv', index=False)

    print('[done] wrote stage_results.csv, restart_detail.csv, per_anchor_fit_check.csv, '
          'bootstrap_stage{1..5}.csv', flush=True)


# ===============================================================================================
# [TC-1] / [TC-2] post-audit diagnostics (added 2026-07-20, per coordinator's corrective
# instruction; changelog entry already recorded in v2/Digitization_Protocol_frozen.md). These are
# ADDITIONS ONLY -- run_stage()/main()/MultiAnchorFitter above (the already-audited stage 1-5
# pipeline) are UNCHANGED. Both diagnostics are POINT-ESTIMATE ONLY (5-restart multistart, no
# bootstrap) per the coordinator's spec, to keep runtime bounded.
# ===============================================================================================

# ---------- [TC-1] boundary-relaxation diagnostic (stages 3,4,5) ----------
# Same flavor as v1 Sec.3.5(c): DIAGNOSTIC-ONLY bounds relaxation, core.py's frozen BOUNDS itself
# is untouched -- this is a separate list used only inside this diagnostic function, exactly like
# v1's task_c.py BOUNDS_RELAXED/BOUNDS_RELAXED2 (which also never touched core.Fitter.BOUNDS).
#   - tau lower bound: 200 -> 30 days (probes whether stage 4/5's tau=200 lower-bound-stuck value
#     is a genuine short-tau regime or an artifact of the search floor)
#   - k upper bound: 1.5 -> 3.0 (probes stage 4/5's k=1.5 upper-bound-stuck value)
#   - c upper bound NOT relaxed (stays 1.0): physical ceiling (long-tail plateau share cannot
#     exceed 100%), same reasoning as v1's BOUNDS_RELAXED2 comment in scripts/task_c.py.
#   - lam, p, h bounds unchanged from core.Fitter.BOUNDS.
BOUNDS_TC1_RELAXED = [
    (np.log(0.1), np.log(30)),   # lam: unchanged
    (0.08, 3.0),                  # k: upper 1.5 -> 3.0
    (1e-3, 0.4),                   # p: unchanged
    (0.2, 1.0),                    # c: unchanged (physical ceiling, not relaxed)
    (np.log(30), np.log(6000)),   # tau: lower 200 -> 30 days
    (1.0, 4.0),                    # h: unchanged
]


def run_relaxation_diagnostic_stage(k, adj, events):
    t0 = time.time()
    anchors_list = build_anchor_list(adj, k)
    y = adj.values.astype(float)
    f = MultiAnchorFitter(y, events, anchors_list, shape='exp3', beta=0.0,
                           h_prior=(2.2, 0.2), use_h=True, W=(7, 14))
    f.BOUNDS = BOUNDS_TC1_RELAXED  # diagnostic-only override, instance attribute (class/core.py untouched)

    results = multistart_fit(f, th0_list())
    best_th, best_obj, best_th0 = results[0]
    rep = f.report(best_th)
    stuck = bounds_stuck(best_th, f)
    d180_restart = [f.report(th)['D180'] for th, _, _ in results]
    obj_restart = [o for _, o, _ in results]

    row = dict(
        stage=k, anchors=','.join(a['label'] for a in anchors_list),
        tau_bound_lo_relaxed=float(np.exp(f.BOUNDS[4][0])), tau_bound_hi_relaxed=float(np.exp(f.BOUNDS[4][1])),
        k_bound_hi_relaxed=float(f.BOUNDS[1][1]),
        lam=rep['lam'], k_param=rep['k'], p=rep['p'], c=rep['c'], tau=rep['tau'], h=rep['h'],
        D1=rep['D1'], D7=rep['D7'], D30=rep['D30'], D180=rep['D180'], rmse=rep['rmse'],
        obj_best=best_obj, n_bounds_stuck=len(stuck),
        bounds_stuck_detail=';'.join(f'{n}@{side}' for n, side, _ in stuck),
        obj_restart_lo=min(obj_restart), obj_restart_hi=max(obj_restart),
        D180_restart_lo=min(d180_restart), D180_restart_hi=max(d180_restart),
        fit_time_s=time.time() - t0,
    )
    print(f'[TC-1 relax] stage {k}: tau={rep["tau"]:.2f}d (relaxed bounds [30,6000]) '
          f'k={rep["k"]:.4f} (relaxed upper 3.0) D180={rep["D180"]:.5f} stuck={stuck} '
          f'time={row["fit_time_s"]:.1f}s', flush=True)
    return row


def run_tc1_relaxation_diagnostic(stages=(3, 4, 5)):
    adj, events = load_series_and_events()
    out_path = OUT_DIR / 'diagnostics_relaxation.csv'
    existing = pd.read_csv(out_path).to_dict('records') if out_path.exists() else []
    done = {int(r['stage']) for r in existing}
    rows = list(existing)
    for k in stages:
        if k in done:
            print(f'[TC-1 relax] stage {k} already in {out_path.name}, skipping.', flush=True)
            continue
        row = run_relaxation_diagnostic_stage(k, adj, events)
        rows.append(row)
        rows.sort(key=lambda r: int(r['stage']))
        pd.DataFrame(rows).to_csv(out_path, index=False)  # incremental save after EVERY fit
    print(f'[TC-1 relax] done, wrote {out_path}', flush=True)
    return pd.DataFrame(rows)


# ---------- [TC-2] regime-flip sensitivity diagnostic (stages 4,5) ----------
# [FORCED EXTENSION, logged]: to give the Collab event a 28-day window while every other
# is_long event keeps core.py's default 14 days, WindowOverrideFitter takes an EXPLICIT
# per-event (day_index, W) list instead of core.Fitter's (day_index, is_long) + fixed W=(7,14)
# convention. This is needed ONLY for this diagnostic (condition i/iii); core.py's LIMBUS_EVENTS
# calendar itself, and the main stage 1-5 pipeline above, are untouched. Duplicates
# MultiAnchorFitter.solve()/report() rather than modifying it, to keep zero risk to the
# already-audited stage pipeline (coordinator instruction: "既存のステージ実行部は変更しない").
def build_events_window_list(adj, collab_window=14):
    idx = {d: i for i, d in enumerate(adj.index)}
    out = []
    for name, (d, lng) in core.LIMBUS_EVENTS.items():
        W = collab_window if name == 'Collab' else (14 if lng else 7)
        out.append((idx[pd.Timestamp(d)], W))
    for d in core.SPIKES:
        out.append((idx[pd.Timestamp(d)], 7))
    return out


def build_anchor_list_tc2(adj, k, nov_override_2510=False):
    """Same as build_anchor_list(), except (condition ii/iii) the 25-10 anchor's mask month
    is moved from October 2025 (protocol default, see module docstring) to November 2025 --
    probing whether the STEAM DAU=418,242 label (chart terminal date 2025-11-09, see
    anchor_verification.md) is better matched to the broadcast's data end month than to the
    "25-10" label's nominal month. Diagnostic-only; the protocol's default October convention is
    NOT changed for the main stage 1-5 pipeline."""
    out = []
    for spec in ANCHOR_SPECS[:k]:
        s = dict(spec)
        label = s['label']
        if nov_override_2510 and label == '25-10':
            s['month_start'], s['month_end'] = '2025-11-01', '2025-11-30'
            label = '25-10(Nov-mask)'
        mask = np.asarray((adj.index >= s['month_start']) & (adj.index <= s['month_end']))
        out.append(dict(label=label, mask=mask, dau=s['dau'], mau=s['mau'], n_days_in_mask=int(mask.sum())))
    return out


class WindowOverrideFitter(MultiAnchorFitter):
    """[TC-2 diagnostic only] events_w: explicit list of (day_index, W) pulse widths, replacing
    core.build_X's is_long-driven W=(7,14) convention for specific events (here: Collab's window).
    Duplicates MultiAnchorFitter.solve()/report() with core.build_X/build_X_kernel swapped for an
    explicit-window equivalent; everything else (NNLS anchor rows, ratio penalty, h prior, BOUNDS,
    optimizer) is identical to MultiAnchorFitter."""

    def __init__(self, y, events, events_w, anchors_list, **kwargs):
        super().__init__(y, events, anchors_list, **kwargs)
        self.events_w = events_w

    def _design_X(self, kernel):
        cols = [np.cumsum(kernel)]
        for t0, W in self.events_w:
            s = core.pulse_shape(self.shape, W) if self.shape != 'delta' else core.pulse_shape('delta', 1)
            col = np.zeros(self.T)
            for j, sj in enumerate(s):
                if t0 + j >= self.T:
                    break
                col[t0 + j:] += sj * kernel[:self.T - t0 - j]
            cols.append(col)
        return np.column_stack(cols)

    def solve(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        R = core.R_kernel(self.T, lam, k, p, c, tau)
        X = self._design_X(R)
        mult = (h / 24.0) * (1.0 + self.beta * self.g)
        rows = [X * (mult / self.y)[:, None]]
        rhs = [np.ones(self.T)]
        A = Xm = None
        if self.anchors_list:
            A = core.A_kernel(self.T, lam, k, p, tau)
            Xm = self._design_X(A)
            for anc in self.anchors_list:
                m = anc['mask']
                dau_row = X[m].mean(axis=0)
                mau_row = Xm[m].mean(axis=0)
                rows.append((dau_row / (0.05 * anc['dau']))[None, :])
                rhs.append(np.array([1 / 0.05]))
                rows.append((mau_row / (0.03 * anc['mau']))[None, :])
                rhs.append(np.array([1 / 0.03]))
        Xw = np.vstack(rows)
        yw = np.concatenate(rhs)
        a, rnorm = nnls(Xw, yw, maxiter=10 * Xw.shape[1])
        loss = rnorm ** 2
        if self.anchors_list:
            for anc in self.anchors_list:
                m = anc['mask']
                dau = X[m].mean(axis=0) @ a
                mau = Xm[m].mean(axis=0) @ a
                r0 = anc['dau'] / anc['mau']
                rs = 0.02
                loss += ((dau / mau - r0) / rs) ** 2
        if self.use_h:
            loss += ((h - self.h_prior[0]) / self.h_prior[1]) ** 2
        return loss, a

    def report(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        loss, a = self.solve(th)
        R = core.R_kernel(max(self.T, 400), lam, k, p, c, tau)
        X = self._design_X(core.R_kernel(self.T, lam, k, p, c, tau))
        mult = (h / 24.0) * (1.0 + self.beta * self.g)
        pred = (X @ a) * mult
        rmse = np.sqrt(np.mean(((pred - self.y) / self.y) ** 2))
        out = dict(D1=R[1], D7=R[7], D30=R[30], D180=R[180],
                   D365=R[365] if len(R) > 365 else np.nan,
                   lam=lam, k=k, p=p, c=c, tau=tau, h=h, rmse=rmse, loss=loss,
                   base=a[0], amps=a[1:], pred=pred)
        return out


TC2_CONDITIONS = [
    ('i_collab28', dict(collab_window=28, nov_override_2510=False)),
    ('ii_novmask', dict(collab_window=14, nov_override_2510=True)),
    ('iii_both', dict(collab_window=28, nov_override_2510=True)),
]


def run_sensitivity_diagnostic_cell(stage, cond_name, collab_window, nov_override_2510, adj):
    t0 = time.time()
    events_w = build_events_window_list(adj, collab_window=collab_window)
    idx = {d: i for i, d in enumerate(adj.index)}
    events = []
    for name, (d, lng) in core.LIMBUS_EVENTS.items():
        events.append((idx[pd.Timestamp(d)], lng))
    for d in core.SPIKES:
        events.append((idx[pd.Timestamp(d)], False))
    anchors_list = build_anchor_list_tc2(adj, stage, nov_override_2510=nov_override_2510)
    y = adj.values.astype(float)
    f = WindowOverrideFitter(y, events, events_w, anchors_list, shape='exp3', beta=0.0,
                              h_prior=(2.2, 0.2), use_h=True, W=(7, 14))
    # BOUNDS unchanged (original core.Fitter.BOUNDS) -- TC-2 probes calendar/mask, not bounds.

    results = multistart_fit(f, th0_list())
    best_th, best_obj, best_th0 = results[0]
    rep = f.report(best_th)
    stuck = bounds_stuck(best_th, f)
    tau_regime = ('lower' if any(n == 'log(tau)' and side == 'lower' for n, side, _ in stuck)
                  else ('upper' if any(n == 'log(tau)' and side == 'upper' for n, side, _ in stuck) else 'interior'))
    d180_restart = [f.report(th)['D180'] for th, _, _ in results]
    obj_restart = [o for _, o, _ in results]

    row = dict(
        stage=stage, condition=cond_name, collab_window=collab_window,
        anchor_2510_month=('2025-11' if nov_override_2510 else '2025-10'),
        anchors=','.join(a['label'] for a in anchors_list),
        lam=rep['lam'], k_param=rep['k'], p=rep['p'], c=rep['c'], tau=rep['tau'], h=rep['h'],
        D1=rep['D1'], D7=rep['D7'], D30=rep['D30'], D180=rep['D180'], rmse=rep['rmse'],
        obj_best=best_obj, n_bounds_stuck=len(stuck),
        bounds_stuck_detail=';'.join(f'{n}@{side}' for n, side, _ in stuck),
        tau_regime=tau_regime,
        obj_restart_lo=min(obj_restart), obj_restart_hi=max(obj_restart),
        D180_restart_lo=min(d180_restart), D180_restart_hi=max(d180_restart),
        fit_time_s=time.time() - t0,
    )
    print(f'[TC-2 sens] stage {stage} cond={cond_name} (collab_W={collab_window}d, '
          f'25-10 mask={row["anchor_2510_month"]}): tau={rep["tau"]:.2f}d ({tau_regime}) '
          f'k={rep["k"]:.4f} D180={rep["D180"]:.5f} stuck={stuck} time={row["fit_time_s"]:.1f}s',
          flush=True)
    return row


def run_tc2_sensitivity_diagnostic(stages=(4, 5)):
    adj, events = load_series_and_events()
    out_path = OUT_DIR / 'diagnostics_sensitivity.csv'
    existing = pd.read_csv(out_path).to_dict('records') if out_path.exists() else []
    done = {(int(r['stage']), r['condition']) for r in existing}
    rows = list(existing)
    for k in stages:
        for cond_name, kwargs in TC2_CONDITIONS:
            if (k, cond_name) in done:
                print(f'[TC-2 sens] stage {k} cond={cond_name} already in {out_path.name}, skipping.',
                      flush=True)
                continue
            row = run_sensitivity_diagnostic_cell(k, cond_name, adj=adj, **kwargs)
            rows.append(row)
            rows.sort(key=lambda r: (int(r['stage']), r['condition']))
            pd.DataFrame(rows).to_csv(out_path, index=False)  # incremental save after EVERY fit
    print(f'[TC-2 sens] done, wrote {out_path}', flush=True)
    return pd.DataFrame(rows)


def run_diagnostics():
    print('[diagnostics] TC-1 (boundary relaxation, stages 3,4,5) starting...', flush=True)
    run_tc1_relaxation_diagnostic()
    print('[diagnostics] TC-2 (regime-flip sensitivity, stages 4,5) starting...', flush=True)
    run_tc2_sensitivity_diagnostic()
    print('[diagnostics] all done.', flush=True)


if __name__ == '__main__':
    if '--diagnostics' in sys.argv:
        run_diagnostics()
    elif '--tc1' in sys.argv:
        run_tc1_relaxation_diagnostic()
    elif '--tc2' in sys.argv:
        run_tc2_sensitivity_diagnostic()
    else:
        main()
