"""Stock-column collinearity / non-negativity-binding diagnostics for task (c) TaskCFitter.

Investigates why the left-censoring "stock" NNLS column (exp(-t/tau)) collapses to exactly 0
for RimWorld/DbD but not for WarThunder/FF14/Warframe, at the RELAXED2-bounds best-fit theta
(outcome/taskC_boundscheck2.csv). Not part of core.py (frozen); core.py is unmodified.

IMPORTANT — basis/definition notes (per stats-auditor's reproducibility audit, round 2):
  - `corr_stock_base` and the two `cond_*` columns are computed on the RAW (unweighted),
    burn-in-excluded design matrix Xraw = Xfull[burn_in:] (columns: [base, event_1..event_N, stock]).
    This is NOT the same matrix NNLS actually solves (NNLS uses the row-reweighted system
    rows = Xfull[burn_in:] * (mult/y)[:,None]); it answers "is this a structural/geometric
    near-linear-dependence of the kernel-generated columns themselves", independent of which
    day-by-day CCU values happen to be observed.
  - `R2_stock_on_others_raw` and `R2_stock_on_others_weighted` are BOTH computed WITHOUT an
    intercept term (NNLS itself has no free intercept column — the "base" column already plays
    that structural role), and BOTH use the UNCENTERED R^2 definition
        R^2 = 1 - ||stock - proj(stock)||^2 / ||stock||^2
    i.e. the squared-cosine of the angle between the stock column and the subspace spanned by
    the remaining columns. Centered R^2 without an intercept is not reported here because it is
    not the quantity relevant to "is stock in the span of the other NNLS columns" and is sensitive
    to an arbitrary shift NNLS never performs. (A previous ad hoc script — since discarded before
    this rewrite — used raw-design-plus-an-explicit-intercept-column, i.e. affine centered R^2;
    that is a *different, less relevant* quantity and is why re-derivations using "raw centered",
    "raw uncentered", or "weighted centered" definitions could not reproduce the old CSV's numbers.
    This script fixes the definition and documents it so results are reproducible.)
  - `ols_stock_coef` (the item-4 non-negativity-binding test) is computed on the WEIGHTED system
    (rows = Xfull[burn_in:] * (mult/y)[:,None], rhs = ones), because that IS the literal linear
    system NNLS solves; the sign of its unconstrained (ordinary least squares, no non-negativity
    constraint) solution tells us directly whether the non-negativity constraint is binding.
  - `ols_stock_coef_tau_up5pct` / `_down5pct`: same OLS-sign test recomputed with tau perturbed by
    +5%/-5% (holding lam,k,p,c fixed), to check how close the sign is to flipping (marginal vs
    robust non-identifiability).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import task_c as tc


def _weighted_system(f, th):
    lam, k, p, c, tau, h = f.unpack(th)
    Xfull, R = f._design(th)
    bi = f.burn_in
    mult = (h / 24.0) * np.ones(f.T)
    rows_w = Xfull[bi:] * (mult[bi:] / f.y[bi:])[:, None]
    rhs_w = np.ones(f.T - bi)
    return rows_w, rhs_w


def _normed_cond(M):
    norms = np.linalg.norm(M, axis=0)
    norms = np.where(norms == 0, 1.0, norms)
    return np.linalg.cond(M / norms)


def _uncentered_r2_no_intercept(target, others):
    beta, *_ = np.linalg.lstsq(others, target, rcond=None)
    resid = target - others @ beta
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum(target ** 2))
    return 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')


def diagnose(title, lam, k, p, c, tau, h=2.2):
    f, adj, events = tc.make_fitter(tc.TITLES[title])
    th = np.array([np.log(lam), k, p, c, np.log(tau), h])
    Xfull, R = f._design(th)
    bi = f.burn_in
    Xraw = Xfull[bi:]                      # RAW, unweighted, burn-in-excluded
    stock_raw = Xraw[:, -1]
    other_raw = Xraw[:, :-1]
    base_raw = Xraw[:, 0]

    corr_stock_base = float(np.corrcoef(stock_raw, base_raw)[0, 1])
    R2_raw = _uncentered_r2_no_intercept(stock_raw, other_raw)

    rows_w, rhs_w = _weighted_system(f, th)
    stock_w = rows_w[:, -1]
    other_w = rows_w[:, :-1]
    R2_weighted = _uncentered_r2_no_intercept(stock_w, other_w)

    cond_full_raw = _normed_cond(Xraw)
    cond_without_stock_raw = _normed_cond(other_raw)

    beta_ols, *_ = np.linalg.lstsq(rows_w, rhs_w, rcond=None)
    ols_stock_coef = float(beta_ols[-1])

    def ols_coef_at_tau(tau_pert):
        th_p = np.array([np.log(lam), k, p, c, np.log(tau_pert), h])
        rows_p, rhs_p = _weighted_system(f, th_p)
        beta_p, *_ = np.linalg.lstsq(rows_p, rhs_p, rcond=None)
        return float(beta_p[-1])

    coef_up = ols_coef_at_tau(tau * 1.05)
    coef_down = ols_coef_at_tau(tau * 0.95)

    return dict(
        title=title, n_events=len(events),
        lam=lam, k=k, p=p, c=c, tau=tau,
        corr_stock_base_raw=corr_stock_base,
        R2_stock_on_others_raw_uncentered_nointercept=R2_raw,
        R2_stock_on_others_weighted_uncentered_nointercept=R2_weighted,
        cond_full_raw=cond_full_raw, cond_without_stock_raw=cond_without_stock_raw,
        ols_stock_coef_weighted=ols_stock_coef,
        ols_stock_negative=bool(ols_stock_coef < 0),
        ols_stock_coef_tau_up5pct=coef_up,
        ols_stock_coef_tau_down5pct=coef_down,
        sign_flip_tau_up5pct=bool((ols_stock_coef < 0) != (coef_up < 0)),
        sign_flip_tau_down5pct=bool((ols_stock_coef < 0) != (coef_down < 0)),
    )


if __name__ == '__main__':
    df_best = pd.read_csv(
        os.path.join(os.path.dirname(__file__), '..', 'outcome', 'taskC_boundscheck2.csv'))
    out_rows = []
    for _, r in df_best.iterrows():
        out_rows.append(diagnose(r['title'], r['lam'], r['k'], r['p'], r['c'], r['tau']))
    out_df = pd.DataFrame(out_rows)
    out_path = os.path.join(os.path.dirname(__file__), '..', 'outcome', 'taskC_stock_collinearity.csv')
    out_df.to_csv(out_path, index=False)
    pd.set_option('display.width', 220)
    pd.set_option('display.max_columns', None)
    print(out_df.to_string(index=False))
