"""Task (b) driver — playtime-contamination beta sweep, reconstructed per Tasks_abc_Report.md sec.2.

Not part of core.py (frozen); core.py is unmodified. Uses core.limbus_setup(shape='exp3', beta=b)
directly (full Limbus anchors + 71-event calendar), refits all 6 Fitter params per beta.

Definitions (reconstructed to match outcome/taskB_beta.csv within numerical-optimization noise):
  - total_act_M: total modeled activations (arrivals) in millions = (base_per_day * T + sum(event
    amplitudes)) / 1e6. Each event-pulse column integrates to unit mass (pulse_shape normalizes to
    sum=1), so its NNLS amplitude IS the total activation count attributed to that event; the base
    column is a constant per-day arrival rate, so base_per_day * T is the cumulative base arrivals
    over the observed window.
  - f_share (reclassification share) -- RECONSTRUCTION CAVEAT (logged explicitly, not silently
    approximated): the exact original f_share formula behind outcome/taskB_beta.csv could NOT be
    uniquely reverse-engineered. The formula below is the best candidate found (monotonic in beta,
    correct order of magnitude, matches at beta=0 trivially), computed as: of the total "excess"
    CCU during event-active days (pred(t) minus the counterfactual base-only, no-boost CCU
    a[0]*cumsum(R)(t)*h/24), the fraction attributable to the playtime-boost multiplier
    (1+beta*g(t)) rather than to modeled extra arrivals:
        excess(t)          = pred(t) - a0 * cumsum(R)(t) * h/24            (t where g(t) > 0)
        playtime_boost(t)  = pred(t) * beta*g(t) / (1 + beta*g(t))         (t where g(t) > 0)
        f_share_proxy = sum_t playtime_boost(t) / sum_t excess(t)
    This UNDERESTIMATES outcome/taskB_beta.csv's reported f_share by a roughly-consistent ~2.2-2.5x
    factor at beta=0.05/0.2/0.3 (e.g. proxy=0.082 vs reported 0.18208 at beta=0.2). Several other
    candidate aggregations were tried (peak-day-only, y-weighted, unweighted/median over the event
    window, event-contribution-only denominators) and none matched either. Critically, the reported
    values (e.g. 0.18208 at beta=0.2) MATHEMATICALLY EXCEED the theoretical ceiling beta/(1+beta) =
    0.16667 that bounds ANY weighted average of beta*g(t)/(1+beta*g(t)) when g(t) is capped at 1
    (as core.Fitter's self.g literally is, via elementwise max across events) -- this is not a tuning
    problem, it proves the original f_share definition used a quantity NOT derivable from core.py's
    Fitter class alone. AUDIT NOTE (2026-07-10): the "uncapped/overlapping-event intensity" hypothesis
    was tested and ruled out numerically -- the uncapped event-sum g reaches at most 1.189 (only 6 days
    exceed 1) at the beta=0.2/0.3 optima, yielding f=0.0823/0.1165, nowhere near the original
    0.182/0.254. The original definition remains unknown; no mechanism is claimed. f_share_proxy is
    reported here as a directionally-consistent, reproducible proxy, explicitly NOT claimed to match
    taskB_beta.csv, and is NOT used as a headline sensitivity metric (beta itself is the axis).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from scipy.optimize import minimize
import core


def fit_multistart(beta, th0_list):
    f, adj = core.limbus_setup(shape='exp3', beta=beta)
    results = []
    for th0 in th0_list:
        th0 = np.asarray(th0, dtype=float)
        th, res = f.fit(th0=th0)
        obj = f.obj_b(th)
        results.append((th, obj, th0))
    results.sort(key=lambda t: t[1])
    return f, results


def report_extra(f, th):
    rep = f.report(th)
    a = np.concatenate([[rep['base']], rep['amps']])
    lam, k, p, c, tau, h = rep['lam'], rep['k'], rep['p'], rep['c'], rep['tau'], rep['h']
    R = core.R_kernel(f.T, lam, k, p, c, tau)
    X = core.build_X(f.T, f.events, R, f.shape, f.W[0], f.W[1])
    base_col = X[:, 0]
    pred = rep['pred']
    baseline = a[0] * base_col * (h / 24.0)
    g = f.g
    mask = g > 0
    beta_val = f.beta
    playtime_boost = pred * beta_val * g / (1 + beta_val * g)
    excess = pred - baseline
    denom = excess[mask].sum()
    f_share_proxy = playtime_boost[mask].sum() / denom if denom != 0 else 0.0
    total_act_M = (a[0] * f.T + a[1:].sum()) / 1e6
    ratio = None
    if f.anchors:
        ratio = rep.get('jan_dau', np.nan) / rep.get('jan_mau', np.nan) if rep.get('jan_mau', 0) else np.nan
    out = dict(D1=rep['D1'], D7=rep['D7'], D30=rep['D30'], D180=rep['D180'],
               rmse=rep['rmse'], h=h, lam=lam, k=k, p=p, c=c, tau=tau,
               f_share_proxy=f_share_proxy, total_act_M=total_act_M, ratio=ratio)
    return out


DEFAULT_TH0 = np.array([np.log(1.5), 0.30, 0.023, 0.90, np.log(6000.0), 2.2])


def th0_set(extra=None):
    starts = [DEFAULT_TH0.copy()]
    if extra:
        starts.extend(extra)
    return starts


def run(betas, th0_extra_by_beta, out_csv, restart_csv):
    rows, restart_rows = [], []
    for beta in betas:
        extras = th0_extra_by_beta.get(beta, [])
        th0_list = th0_set(extras)
        f, results = fit_multistart(beta, th0_list)
        best_th, best_obj, best_th0 = results[0]
        rep = report_extra(f, best_th)
        row = dict(beta=beta, **rep, n_restarts=len(results),
                   obj_restart_lo=min(o for _, o, _ in results),
                   obj_restart_hi=max(o for _, o, _ in results))
        row['obj_rel_spread'] = ((row['obj_restart_hi'] - row['obj_restart_lo']) / row['obj_restart_lo']
                                  if row['obj_restart_lo'] != 0 else float('nan'))
        d180_vals = []
        for th, obj, th0 in results:
            r = report_extra(f, th)
            d180_vals.append(r['D180'])
        row['D180_restart_lo'] = min(d180_vals)
        row['D180_restart_hi'] = max(d180_vals)
        rows.append(row)
        for i, (th, obj, th0) in enumerate(results):
            r = report_extra(f, th)
            restart_rows.append(dict(beta=beta, restart_idx=i, obj=obj,
                                      D1=r['D1'], D7=r['D7'], D30=r['D30'], D180=r['D180'],
                                      rmse=r['rmse'], h=r['h'], lam=r['lam'], k=r['k'], p=r['p'],
                                      c=r['c'], tau=r['tau'], f_share_proxy=r['f_share_proxy'],
                                      total_act_M=r['total_act_M'],
                                      th0_log_lam=th0[0], th0_k=th0[1], th0_p=th0[2],
                                      th0_c=th0[3], th0_log_tau=th0[4], th0_h=th0[5]))
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    pd.DataFrame(restart_rows).to_csv(restart_csv, index=False)
    return df


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--betas', required=True, help='comma-separated beta values')
    ap.add_argument('--out', required=True)
    ap.add_argument('--restart_out', required=True)
    args = ap.parse_args()
    betas = [float(x) for x in args.betas.split(',')]
    # near-optimum starts for beta=0.2 and beta=0.3 (md-reported region), offered as extra
    # multistart inits at every beta in this sweep, per the coordinator's Step-2 spec.
    extra_02 = np.array([np.log(1.6), 0.35, 0.05, 0.85, np.log(2000.0), 2.14])
    extra_03 = np.array([np.log(1.3), 0.25, 0.10, 0.80, np.log(1500.0), 2.13])
    th0_extra = {b: [extra_02, extra_03] for b in betas}
    run(betas, th0_extra, args.out, args.restart_out)
