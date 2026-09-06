"""
E3: isolate the anchor-WEIGHTING-magnitude effect alone (archive's 100x/31x/31x/10x
post-hoc penalties + anchor-agnostic nnls for 'a') while keeping CORE's own BOUNDS
(not archive's). This separates "objective weighting" from "search-range" as causes.
"""
import time, json, sys
import numpy as np
from pathlib import Path
from scipy.optimize import nnls, minimize

sys.path.insert(0, r"C:\Users\rhiro\MacData\Riga\Limbus-Paper\scripts")
import core

OUT = Path(__file__).resolve().parent

f, adj = core.limbus_setup()
T = f.T
events = f.events
y = f.y
N_NAMED = len(core.LIMBUS_EVENTS)
N_SPIKE = len(core.SPIKES)
n_cols = 1 + N_NAMED + N_SPIKE
spike_col_idx = np.arange(1 + N_NAMED, n_cols)
ridge = np.sqrt(1e-9)

anchors = f.anchors
jan_mask = anchors['jan_mask']
dau_target = anchors['dau']; mau_target = anchors['mau']
ratio0, ratio_s = anchors['ratio']

CORE_BOUNDS = f.BOUNDS  # unmodified core.py BOUNDS


def unpack(th):
    lam = np.exp(th[0]); k = th[1]; p = th[2]; c = th[3]; tau = np.exp(th[4]); h = th[5]
    return lam, k, p, c, tau, h


def solve_archiveweight_corebounds(th):
    lam, k, p, c, tau, h = unpack(th)
    R = core.R_kernel(T, lam, k, p, c, tau)
    X = core.build_X(T, events, R, 'exp3', 7, 14)
    mult = (h / 24.0)
    M = X * mult
    wgt = 1.0 / (0.06 * y)
    Mw = M * wgt[:, None]
    yw = y * wgt
    rr = np.zeros((len(spike_col_idx), n_cols))
    rr[np.arange(len(spike_col_idx)), spike_col_idx] = ridge
    Aeq = np.vstack([Mw, rr])
    beq = np.concatenate([yw, np.zeros(len(spike_col_idx))])
    a, _ = nnls(Aeq, beq, maxiter=10 * n_cols)
    pred = X @ a * mult
    rel = (pred - y) / y
    sse = np.sum((rel / 0.06) ** 2)
    A = core.A_kernel(T, lam, k, p, tau)
    Xm = core.build_X_kernel(T, events, A, 'exp3', 7, 14)
    dau_series = X @ a
    mau_series = Xm @ a
    dau_jan = dau_series[jan_mask].mean()
    mau_jan = mau_series[jan_mask].mean()
    dmr = dau_series[240:] / np.maximum(mau_series[240:], 1e-9)
    pen = 0.0
    pen += 100 * ((np.mean(dmr) - ratio0) / ratio_s) ** 2
    pen += 31 * ((dau_jan - dau_target) / (0.05 * dau_target)) ** 2
    pen += 31 * ((mau_jan - mau_target) / (0.03 * mau_target)) ** 2
    pen += 10 * ((h - 2.2) / 0.2) ** 2
    loss = sse + pen
    return loss, a, pred, dict(dau_jan=dau_jan, mau_jan=mau_jan, dm=np.mean(dmr), rel_rmse=float(np.sqrt(np.mean(rel**2))))


def obj(th):
    lo = np.array([b[0] for b in CORE_BOUNDS]); hi = np.array([b[1] for b in CORE_BOUNDS])
    thc = np.clip(th, lo, hi)
    pen = 1e6 * np.sum((th - thc) ** 2)
    try:
        loss, a, pred, info = solve_archiveweight_corebounds(thc)
        return loss + pen
    except Exception:
        return 1e12


starts_raw = [
    dict(lam=1.5, k=0.30, p=0.023, c=0.90, tau=6000.0, h=2.2),   # core default
    dict(lam=1.683, k=0.301, p=0.026, c=0.852, tau=2070.0, h=2.216),
    dict(lam=6.0, k=0.8, p=0.30, c=0.60, tau=700.0, h=2.2),
]

lo = np.array([b[0] for b in CORE_BOUNDS]); hi = np.array([b[1] for b in CORE_BOUNDS])
names = ['log(lam)', 'k', 'p', 'c', 'log(tau)', 'h']

results = []
best = None
for i, s in enumerate(starts_raw):
    th0 = np.clip(np.array([np.log(s['lam']), s['k'], s['p'], s['c'], np.log(s['tau']), s['h']]), lo, hi)
    t0 = time.time()
    res = minimize(obj, th0, method='Nelder-Mead', options=dict(maxiter=2000, xatol=1e-5, fatol=1e-6))
    thc = np.clip(res.x, lo, hi)
    loss, a, pred, info = solve_archiveweight_corebounds(thc)
    lam, k, p, c, tau, h = unpack(thc)
    dt = time.time() - t0
    print(f"[start {i}] loss={loss:.3f} rel_rmse={info['rel_rmse']:.6f} lam={lam:.3f} k={k:.3f} p={p:.4f} "
          f"c={c:.4f} tau={tau:.1f} h={h:.3f} dau_jan={info['dau_jan']:.0f} mau_jan={info['mau_jan']:.0f} "
          f"time={dt:.1f}s nfev={res.nfev}")
    results.append(dict(start=i, loss=float(loss), rel_rmse=info['rel_rmse'], lam=lam, k=k, p=p, c=c, tau=tau, h=h))
    if best is None or loss < best[0]:
        best = (loss, thc, info)

loss_b, th_b, info_b = best
lam, k, p, c, tau, h = unpack(th_b)
hits = []
for i, n in enumerate(names):
    if abs(th_b[i] - lo[i]) < 1e-4:
        hits.append(f'{n}=LOW')
    elif abs(th_b[i] - hi[i]) < 1e-4:
        hits.append(f'{n}=HIGH')
print('\n=== BEST (core bounds + archive weighting) ===')
print(f"loss={loss_b:.3f} rel_rmse={info_b['rel_rmse']:.6f} lam={lam:.3f} k={k:.3f} p={p:.4f} c={c:.4f} tau={tau:.1f} h={h:.3f}")
print('bound hits (core bounds):', hits)

with open(OUT / 'ablation3_coreboundsarchiveweight.json', 'w') as fh:
    json.dump(dict(all_starts=results, best=dict(loss=float(loss_b), rel_rmse=info_b['rel_rmse'], lam=lam, k=k, p=p,
                                                  c=c, tau=tau, h=h, bound_hits=hits)), fh, indent=2, default=str)
print('saved ablation3_coreboundsarchiveweight.json')
