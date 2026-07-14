"""Task (c) driver — 7-title generalization, reconstructed per Tasks_abc_Report.md sec.3 frozen spec.

Not part of core.py's frozen kernel/BOUNDS definitions. This module:
  - detects events mechanically (21d rolling median ratio > 1.35, weekday-adjusted, min gap 5d,
    series head = launch-like 14d pulse)
  - adds a left-censoring "stock" NNLS column exp(-t/tau) sharing tau with the kernel's long-tail term
  - excludes the first 90 days from the loss (burn-in)
  - no anchors; h prior N(2.2,0.2) only
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import numpy as np, pandas as pd
from scipy.optimize import nnls, minimize
import core

DATA_DIR = core.DATA_DIR

TITLES = {
    'WarThunder': 'WT_steamdb_chart_236390.csv',
    'FF14':       'FF14_steamdb_chart_39210.csv',
    'RimWorld':   'RM_steamdb_chart_294100.csv',
    'Warframe':   'WF_steamdb_chart_230410.csv',
    'DbD':        'DbD_steamdb_chart_381210.csv',
    'PoE':        'PoE_steamdb_chart_238960.csv',
    'PoE2':       'PoE2_steamdb_chart_2694490.csv',
}

# Per-title left-censoring config (frozen per Tasks_abc_Report.md sec.3): titles whose SteamDB
# 'Average Players' window starts before their real launch (left-censored) get the stock
# regressor + 90d burn-in; PoE2 is observed from its actual launch (2024-12-06) so gets neither.
TITLE_CONFIG = {
    'WarThunder': dict(burn_in=90, use_stock=True),
    'FF14':       dict(burn_in=90, use_stock=True),
    'RimWorld':   dict(burn_in=90, use_stock=True),
    'Warframe':   dict(burn_in=90, use_stock=True),
    'DbD':        dict(burn_in=90, use_stock=True),
    'PoE':        dict(burn_in=90, use_stock=True),
    'PoE2':       dict(burn_in=0,  use_stock=False),
}

# ---------- frozen event detection ----------
def detect_events(adj, window=21, thresh=1.35, min_gap=5):
    med = adj.rolling(window, center=True, min_periods=(window//2+1)).median()
    ratio = adj / med
    cand = np.where(ratio.values > thresh)[0]
    cand = cand[cand > 0]  # day0 handled separately as launch pulse
    selected = []
    last = -10**9
    for idx in cand:
        if idx - last >= min_gap:
            selected.append(idx)
            last = idx
    events = [(0, True)]  # series head: launch-like pulse, 14d width
    for idx in selected:
        events.append((int(idx), False))
    return events


class TaskCFitter(core.Fitter):
    """Fitter variant: optionally adds a left-censoring stock column exp(-t/tau) to the NNLS
    design (use_stock=True, for titles whose observed window starts before their real launch),
    and excludes the first `burn_in` days from the fitted loss (rows only; kernel/report use full
    length). No anchors are supported/used here. use_stock=False (e.g. PoE2, observed from its
    actual launch) skips the stock column entirely -- not merely fits it to 0."""

    def __init__(self, y, events, burn_in=90, use_stock=True, **kwargs):
        kwargs.setdefault('anchors', None)
        kwargs.setdefault('beta', 0.0)
        super().__init__(y, events, **kwargs)
        self.burn_in = burn_in
        self.use_stock = use_stock

    def _design(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        R = core.R_kernel(self.T, lam, k, p, c, tau)
        X = core.build_X(self.T, self.events, R, self.shape, self.W[0], self.W[1])
        if self.use_stock:
            stock_col = np.exp(-np.arange(self.T, dtype=float) / tau)
            Xfull = np.column_stack([X, stock_col])
        else:
            Xfull = X
        return Xfull, R

    def solve(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        Xfull, R = self._design(th)
        mult = (h / 24.0) * np.ones(self.T)  # beta=0 fixed for task c
        bi = self.burn_in
        rows = Xfull[bi:] * (mult[bi:] / self.y[bi:])[:, None]
        rhs = np.ones(self.T - bi)
        a, rnorm = nnls(rows, rhs, maxiter=20 * rows.shape[1])
        loss = rnorm ** 2
        if self.use_h:
            loss += ((h - self.h_prior[0]) / self.h_prior[1]) ** 2
        return loss, a

    def report(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        loss, a = self.solve(th)
        Xfull, R = self._design(th)
        mult = (h / 24.0) * np.ones(self.T)
        pred = (Xfull @ a) * mult
        bi = self.burn_in
        rmse = np.sqrt(np.mean(((pred[bi:] - self.y[bi:]) / self.y[bi:]) ** 2))
        base_per_day = float(a[0])
        stock = float(a[-1]) if self.use_stock else 0.0
        amps = a[1:-1] if self.use_stock else a[1:]
        out = dict(D1=R[1], D7=R[7], D30=R[30], D180=R[180],
                    D365=R[365] if len(R) > 365 else np.nan,
                    lam=lam, k=k, p=p, c=c, tau=tau, h=h,
                    rmse=rmse, loss=loss, base_per_day=base_per_day, stock=stock,
                    amps=amps, pred=pred)
        return out


def load_title(fname):
    s = core.load_steamdb_daily(DATA_DIR / fname)
    adj, _ = core.weekday_adjust(s)
    adj = adj.dropna()
    return adj


def make_fitter(fname, bounds=None, burn_in=None, use_stock=None, title=None):
    """burn_in/use_stock default to TITLE_CONFIG[title] when title is given (preferred);
    explicit burn_in/use_stock args override for ad hoc use."""
    adj = load_title(fname)
    events = detect_events(adj)
    y = adj.values.astype(float)
    if title is not None and title in TITLE_CONFIG:
        cfg = TITLE_CONFIG[title]
        if burn_in is None:
            burn_in = cfg['burn_in']
        if use_stock is None:
            use_stock = cfg['use_stock']
    if burn_in is None:
        burn_in = 90
    if use_stock is None:
        use_stock = True
    f = TaskCFitter(y, events, burn_in=burn_in, use_stock=use_stock, shape='exp3', anchors=None,
                    beta=0.0, h_prior=(2.2, 0.2), use_h=True, W=(7, 14))
    if bounds is not None:
        f.BOUNDS = bounds
    return f, adj, events


def multistart_fit(f, th0_list):
    """Run f.fit() from multiple initial points (bypassing core.Fitter.fit's fixed default),
    return best (lowest obj_b) theta plus the list of all restart thetas/objs for stability check.
    Each result tuple is (th, obj, th0_used) so restart_idx <-> initial value is traceable."""
    results = []
    for th0 in th0_list:
        th0 = np.asarray(th0, dtype=float)
        lo = np.array([b[0] for b in f.BOUNDS]); hi = np.array([b[1] for b in f.BOUNDS])
        th0c = np.clip(th0, lo, hi)
        res = minimize(f.obj_b, th0c, method='Powell', bounds=f.BOUNDS,
                        options=dict(maxiter=3000, xtol=1e-5, ftol=1e-7))
        th = np.clip(res.x, lo, hi)
        res2 = minimize(f.obj_b, th, method='Nelder-Mead',
                         options=dict(maxiter=800, xatol=1e-5, fatol=1e-8))
        th2 = np.clip(res2.x, lo, hi)
        if f.obj_b(th2) < f.obj_b(th):
            th = th2
        results.append((th, f.obj_b(th), th0))
    results.sort(key=lambda t: t[1])
    return results


def th0_original_region():
    """Initial guesses near the taskC_7titles.csv (original-bounds) optimum region + a couple of
    generic starts, for step-1 reproduction (3 starts)."""
    return [
        np.array([np.log(20.0), 1.0, 0.05, 0.9, np.log(1000.0), 2.2]),
        np.array([np.log(1.5), 0.30, 0.023, 0.90, np.log(6000.0), 2.2]),
        np.array([np.log(5.0), 0.5, 0.10, 0.5, np.log(3000.0), 2.2]),
    ]


def th0_relaxed_region():
    """>=5 starts for the relaxed-bounds refit: includes near-old-optimum + newly opened region."""
    return [
        np.array([np.log(20.0), 1.0, 0.05, 0.9, np.log(1000.0), 2.2]),   # near old optimum
        np.array([np.log(1.5), 0.30, 0.023, 0.90, np.log(6000.0), 2.2]),  # core default
        np.array([np.log(60.0), 0.8, 0.2, 0.7, np.log(2000.0), 2.2]),    # new region: lam~60
        np.array([np.log(150.0), 0.6, 0.3, 0.6, np.log(3000.0), 2.2]),   # new region: lam~150
        np.array([np.log(10.0), 0.7, 0.5, 0.5, np.log(1500.0), 2.2]),    # new region: p~0.5
        np.array([np.log(100.0), 1.0, 0.4, 0.8, np.log(4000.0), 2.2]),   # new region: lam~100,p~0.4
    ]


def th0_relaxed2_region():
    """Diagnostic-relaxation (RELAXED2) starts: near-RELAXED-optimum per title is impossible to
    encode generically here (5 titles differ), so we use a title-agnostic spread of >=6 starts
    that explicitly includes the newly-opened c~0.1 region and k~2.0-2.5 region requested by the
    stats-auditor, in addition to points near the previous (RELAXED) optima and the core default."""
    return [
        np.array([np.log(20.0), 1.0, 0.05, 0.9, np.log(1000.0), 2.2]),    # near original-bounds optimum region
        np.array([np.log(1.5), 0.30, 0.023, 0.90, np.log(6000.0), 2.2]),  # core default
        np.array([np.log(60.0), 0.8, 0.2, 0.7, np.log(2000.0), 2.2]),     # RELAXED-era: lam~60
        np.array([np.log(100.0), 1.4, 0.4, 0.9, np.log(4000.0), 2.2]),    # RELAXED-era: lam~100,p~0.4
        np.array([np.log(30.0), 2.2, 0.1, 0.1, np.log(8000.0), 2.2]),     # NEW: c~0.1, k~2.2
        np.array([np.log(45.0), 2.5, 0.3, 0.1, np.log(15000.0), 2.2]),    # NEW: c~0.1, k~2.5, tau~15000
        np.array([np.log(25.0), 2.0, 0.5, 0.15, np.log(10000.0), 2.2]),   # NEW: k~2.0, p~0.5, c~0.15
    ]


BOUNDS_ORIGINAL = core.Fitter.BOUNDS
BOUNDS_RELAXED = [(np.log(0.1), np.log(200)), (0.08, 1.5), (1e-4, 0.6), (0.2, 1.0),
                   (np.log(200), np.log(6000)), (1.0, 4.0)]
# Diagnostic-only further relaxation (requested by stats-auditor after Step2 audit).
# Only c (lower bound), k (upper bound), tau (upper bound) are additionally loosened; lam/p keep
# the Step2 (RELAXED) range; c's physical upper bound of 1.0 (long-tail plateau share cannot exceed
# 100%) is NOT relaxed.
BOUNDS_RELAXED2 = [(np.log(0.1), np.log(200)), (0.08, 3.0), (1e-4, 0.6), (0.05, 1.0),
                    (np.log(200), np.log(20000)), (1.0, 4.0)]


def run(bounds, th0_fn, out_csv, restart_csv=None, titles=None):
    rows = []
    restart_rows = []
    title_items = TITLES.items() if titles is None else [(t, TITLES[t]) for t in titles]
    for title, fname in title_items:
        f, adj, events = make_fitter(fname, bounds=bounds, title=title)
        th0_list = th0_fn()
        results = multistart_fit(f, th0_list)
        best_th, best_obj, best_th0 = results[0]
        rep = f.report(best_th)
        n_events = len(events)
        row = dict(title=title, T=f.T, n_events=n_events,
                   D1=rep['D1'], D7=rep['D7'], D30=rep['D30'], D180=rep['D180'],
                   rmse=rep['rmse'], base_per_day=rep['base_per_day'], stock=rep['stock'],
                   lam=rep['lam'], k=rep['k'], p=rep['p'], c=rep['c'], tau=rep['tau'], h=rep['h'])
        # restart spread on D180 and on the objective (recompute for every restart theta)
        d180_vals = []
        obj_vals = []
        for th, obj, th0 in results:
            r = f.report(th)
            d180_vals.append(r['D180'])
            obj_vals.append(obj)
        row['D180_restart_lo'] = min(d180_vals)
        row['D180_restart_hi'] = max(d180_vals)
        row['obj_restart_lo'] = min(obj_vals)
        row['obj_restart_hi'] = max(obj_vals)
        row['obj_rel_spread'] = (max(obj_vals) - min(obj_vals)) / min(obj_vals) if min(obj_vals) != 0 else float('nan')
        row['n_restarts'] = len(results)
        rows.append(row)
        for i, (th, obj, th0) in enumerate(results):
            r = f.report(th)
            restart_rows.append(dict(title=title, restart_idx=i, obj=obj,
                                      D180=r['D180'], lam=r['lam'], k=r['k'], p=r['p'],
                                      c=r['c'], tau=r['tau'], stock=r['stock'],
                                      rmse=r['rmse'], h=r['h'],
                                      th0_log_lam=th0[0], th0_k=th0[1], th0_p=th0[2],
                                      th0_c=th0[3], th0_log_tau=th0[4], th0_h=th0[5]))
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    if restart_csv:
        pd.DataFrame(restart_rows).to_csv(restart_csv, index=False)
    return df


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--step', choices=['1', '2', '4'], required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--restart_out', default=None)
    ap.add_argument('--titles', default=None, help='comma-separated subset of TITLES keys')
    args = ap.parse_args()
    titles = args.titles.split(',') if args.titles else None
    if args.step == '1':
        run(BOUNDS_ORIGINAL, th0_original_region, args.out, args.restart_out, titles=titles)
    elif args.step == '2':
        run(BOUNDS_RELAXED, th0_relaxed_region, args.out, args.restart_out, titles=titles)
    else:
        run(BOUNDS_RELAXED2, th0_relaxed2_region, args.out, args.restart_out, titles=titles)
