"""
Ablation experiments to decompose the 0.133 (archive) vs 0.1305 (current core.py)
Limbus baseline relative-RMSE gap. Reads core.py, does NOT modify it. All bound/objective
changes are done via local monkey-patches / standalone functions in this script only.
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
spike_col_idx = np.arange(1 + N_NAMED, 1 + N_NAMED + N_SPIKE)  # column indices in X (0=Base)
print('spike cols range', spike_col_idx.min(), spike_col_idx.max(), 'ncols total expected', 1+N_NAMED+N_SPIKE)

results = {}

def rmse_from_pred(pred, yobs):
    return float(np.sqrt(np.mean(((pred - yobs) / yobs) ** 2)))

def bound_hits(th, bounds, names, tol=1e-6):
    lo = np.array([b[0] for b in bounds]); hi = np.array([b[1] for b in bounds])
    hits = []
    for i, n in enumerate(names):
        if abs(th[i] - lo[i]) < tol * max(1, abs(lo[i])) + 1e-9:
            hits.append(f'{n}=LOW')
        elif abs(th[i] - hi[i]) < tol * max(1, abs(hi[i])) + 1e-9:
            hits.append(f'{n}=HIGH')
    return hits

NAMES = ['log(lam)', 'k', 'p', 'c', 'log(tau)', 'h']

# ---------------- Baseline (reproduce, already known ~0.130457) ----------------
t0 = time.time()
th_base, res_base = f.fit()
out_base = f.report(th_base)
print(f"[E0 baseline] rmse={out_base['rmse']:.6f} params lam={out_base['lam']:.4f} k={out_base['k']:.4f} "
      f"p={out_base['p']:.5f} c={out_base['c']:.4f} tau={out_base['tau']:.2f} h={out_base['h']:.4f} "
      f"time={time.time()-t0:.1f}s")
hits0 = bound_hits(np.clip(th_base, [b[0] for b in f.BOUNDS], [b[1] for b in f.BOUNDS]), f.BOUNDS, NAMES)
print('[E0] bound hits:', hits0)
results['E0_baseline'] = dict(rmse=out_base['rmse'], lam=out_base['lam'], k=out_base['k'], p=out_base['p'],
                               c=out_base['c'], tau=out_base['tau'], h=out_base['h'], bound_hits=hits0)

# ---------------- E1: archive-like BOUNDS, core's normal fit() procedure, default start ----------------
ARCHIVE_BOUNDS = [(np.log(1.5), np.log(120)), (0.12, 3.0), (1e-4, 0.999), (1e-4, 0.999),
                   (np.log(60), np.log(6000)), (1.2, 3.5)]
f1 = core.Fitter(f.y, f.events, shape='exp3', anchors=f.anchors, beta=0.0, W=f.W)
f1.BOUNDS = ARCHIVE_BOUNDS
t0 = time.time()
th1, res1 = f1.fit()
out1 = f1.report(th1)
print(f"[E1 archive-bounds, core fit()] rmse={out1['rmse']:.6f} lam={out1['lam']:.4f} k={out1['k']:.4f} "
      f"p={out1['p']:.5f} c={out1['c']:.4f} tau={out1['tau']:.2f} h={out1['h']:.4f} time={time.time()-t0:.1f}s")
hits1 = bound_hits(th1, ARCHIVE_BOUNDS, NAMES)
print('[E1] bound hits:', hits1)
results['E1_archive_bounds'] = dict(rmse=out1['rmse'], lam=out1['lam'], k=out1['k'], p=out1['p'],
                                     c=out1['c'], tau=out1['tau'], h=out1['h'], bound_hits=hits1)

with open(OUT / 'ablation_partial1.json', 'w') as fh:
    json.dump(results, fh, indent=2, default=str)
print('=== partial checkpoint 1 saved ===')
