"""Paper v1.1 figures — replaces outcome/fig_tasks_abc.png (kept, not deleted, for the record).

Generates:
  fig5_sensitivity.png  — (a) pulse-shape sensitivity (unchanged content vs old fig),
                          (b) playtime-contamination sweep, x-axis = beta (NOT f_share, which was
                          retracted as unreconstructible in the sec.2.1 audit).
  fig6_identifiability.png — D180 restart-range ("non-identifiability width") horizontal bars for
                          Limbus + 7 titles, log x-axis. This is explicitly a range/uncertainty
                          figure, not a point-estimate comparison (per the frozen policy of not
                          quoting point estimates for x/△-rated titles).

All numeric inputs are read from outcome/*.csv (no hardcoded fit results), EXCEPT the Limbus
point+90% CI (3.3%, [2.9%, 4.4%]), which is an external published reference value from
Limbus_Retention_Reanalysis.md / the paper draft, not a byproduct of this repo's scripts.

Run: py scripts/make_figs_paper.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'outcome')


# ---------------------------------------------------------------------------
# Figure 5: (a) pulse-shape sensitivity, (b) beta sweep (x-axis = beta, not f)
# ---------------------------------------------------------------------------
def fig5():
    taskA = pd.read_csv(os.path.join(OUT, 'taskA_shapes.csv'))
    baseline_row = taskA[taskA['variant'].str.startswith('exp3_W7/14')].iloc[0]
    d30_base, d180_base = baseline_row['D30'], baseline_row['D180']
    taskA = taskA[~taskA['variant'].str.startswith('exp3_W7/14')].copy()
    taskA['dD30_pct'] = (taskA['D30'] - d30_base) / d30_base * 100
    taskA['dD180_pct'] = (taskA['D180'] - d180_base) / d180_base * 100

    betaM = pd.read_csv(os.path.join(OUT, 'taskB_beta.csv'))
    betaF = pd.read_csv(os.path.join(OUT, 'taskB_beta_fine.csv'))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

    # --- panel (a) ---
    ax = axes[0]
    x = np.arange(len(taskA))
    w = 0.38
    ax.bar(x - w / 2, taskA['dD30_pct'], width=w, label='D30 % dev. vs exp3 baseline', color='#4C72B0')
    ax.bar(x + w / 2, taskA['dD180_pct'], width=w, label='D180 % dev. vs exp3 baseline', color='#DD8452')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(taskA['variant'], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('% deviation from exp3_W7/14 baseline')
    ax.set_title('(a) Pulse-shape sensitivity (task a)')
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.3)

    # --- panel (b): x-axis = beta ---
    ax = axes[1]
    lo, hi = 2.9, 4.4  # % — beta=0 published 90% CI (Limbus_Retention_Reanalysis.md)
    ax.axhspan(lo, hi, color='gray', alpha=0.18, label='beta=0 90% CI [2.9%, 4.4%]')
    ax.plot(betaM['beta'], betaM['D180'] * 100, 'o-', color='#4C72B0', ms=7,
            label='main grid (taskB_beta.csv)', zorder=3)
    ax.plot(betaF['beta'], betaF['D180'] * 100, 's', color='#C44E52', ms=5,
            label='fine grid (taskB_beta_fine.csv)', zorder=4)
    for _, r in betaM.iterrows():
        ax.annotate(f"{r['beta']:g}", (r['beta'], r['D180'] * 100), textcoords='offset points',
                    xytext=(0, 7), fontsize=7, ha='center', color='#4C72B0')
    for _, r in betaF.iterrows():
        ax.annotate(f"{r['beta']:g}", (r['beta'], r['D180'] * 100), textcoords='offset points',
                    xytext=(0, -11), fontsize=6.5, ha='center', color='#C44E52')
    ax.set_xlabel('beta (playtime-boost coefficient, exogenous)')
    ax.set_ylabel('D180 (%)')
    ax.set_title('(b) Playtime-contamination sweep (task b)')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3)

    fig.tight_layout()
    path = os.path.join(OUT, 'fig5_sensitivity.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Figure 6: D180 restart-range ("identifiability width") per title
# ---------------------------------------------------------------------------
def _rw_main_basin_range():
    df = pd.read_csv(os.path.join(OUT, 'taskC_boundscheck2_restarts.csv'))
    rw = df[df['title'] == 'RimWorld']
    main = rw[rw['restart_idx'] < 5]  # restarts 0-4: dominant basin (obj within ~0.001% of best)
    sub = rw[rw['restart_idx'] >= 5]  # restarts 5-6: sub-basin, rmse +6.3% inferior (see report)
    return main['D180'].min(), main['D180'].max(), sub['D180'].max()


def _dbd_envelope():
    relaxed = pd.read_csv(os.path.join(OUT, 'taskC_boundscheck_restarts.csv'))
    relaxed2 = pd.read_csv(os.path.join(OUT, 'taskC_boundscheck2.csv'))
    hi = relaxed[relaxed['title'] == 'DbD']['D180'].max()          # RELAXED-stage restart_hi
    lo = relaxed2[relaxed2['title'] == 'DbD']['D180_restart_lo'].iloc[0]  # RELAXED2-stage restart_lo
    return lo, hi


def fig6():
    bc2 = pd.read_csv(os.path.join(OUT, 'taskC_boundscheck2.csv')).set_index('title')
    poe = pd.read_csv(os.path.join(OUT, 'taskC_poe_boundscheck.csv')).set_index('title')

    rw_lo, rw_hi, rw_sub = _rw_main_basin_range()
    dbd_lo, dbd_hi = _dbd_envelope()

    rows = []
    rows.append(dict(title='Limbus', lo=0.029, hi=0.044, point=0.033, excluded=False,
                      note='anchored (DAU/MAU); external ref.'))
    rows.append(dict(title='PoE', lo=poe.loc['PoE', 'D180_restart_lo'], hi=poe.loc['PoE', 'D180_restart_hi'],
                      point=None, excluded=False, note='n_events=19(recon)/17(orig)'))
    rows.append(dict(title='FF14', lo=bc2.loc['FF14', 'D180_restart_lo'], hi=bc2.loc['FF14', 'D180_restart_hi'],
                      point=None, excluded=False, note='corner solution (p,c at floor)'))
    rows.append(dict(title='DbD', lo=dbd_lo, hi=dbd_hi, point=None, excluded=False,
                      note='c/p/tau non-identified; D180 relatively contained'))
    rows.append(dict(title='PoE2', lo=poe.loc['PoE2', 'D180_restart_lo'], hi=poe.loc['PoE2', 'D180_restart_hi'],
                      point=None, excluded=True, note='excluded: poor fit (rmse~0.34), tau at ceiling'))
    rows.append(dict(title='RimWorld', lo=rw_lo, hi=rw_hi, point=rw_sub, excluded=False,
                      note='sub-basin (open marker) rmse +6.3% inferior'))
    rows.append(dict(title='Warframe', lo=bc2.loc['Warframe', 'D180_restart_lo'],
                      hi=bc2.loc['Warframe', 'D180_restart_hi'], point=None, excluded=False,
                      note='p now at new upper bound 0.6'))
    rows.append(dict(title='WarThunder', lo=bc2.loc['WarThunder', 'D180_restart_lo'],
                      hi=bc2.loc['WarThunder', 'D180_restart_hi'], point=None, excluded=False,
                      note='multi-modal: sub-basins found within ~2% obj'))

    df = pd.DataFrame(rows)
    df['width_log'] = np.log(df['hi'] / df['lo'].clip(lower=1e-6))
    df = df.sort_values('width_log', ascending=True)
    # keep Limbus pinned at the top regardless of width ordering
    df = pd.concat([df[df.title == 'Limbus'], df[df.title != 'Limbus']])

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ypos = np.arange(len(df))[::-1]
    for y, (_, r) in zip(ypos, df.iterrows()):
        color = '0.6' if r['excluded'] else ('#55A868' if r['title'] == 'Limbus' else '#4C72B0')
        lo_pct, hi_pct = r['lo'] * 100, r['hi'] * 100
        ax.plot([lo_pct, hi_pct], [y, y], '-', color=color, lw=5, alpha=0.8, solid_capstyle='butt',
                zorder=2)
        # endpoint tick marks so near-zero-width ranges (e.g. FF14, PoE, Warframe) stay visible
        ax.plot([lo_pct, hi_pct], [y, y], '|', color=color, ms=11, mew=2.2, zorder=3)
        if r['title'] == 'Limbus':
            ax.plot(0.033 * 100, y, 'D', color=color, ms=8, zorder=5)
        if r['point'] is not None:
            ax.plot(r['point'] * 100, y, 'o', mfc='none', mec='#C44E52', mew=1.8, ms=9, zorder=5)
        far_x = max(hi_pct, (r['point'] or 0) * 100)
        ax.annotate(r['note'], (far_x, y), textcoords='offset points',
                    xytext=(8, 0), fontsize=7.5, va='center', color='0.3' if r['excluded'] else 'black')

    ax.set_yticks(ypos)
    ax.set_yticklabels(df['title'])
    ax.set_xscale('log')
    ax.set_xlim(3e-3, 2e3)
    ax.set_xlabel('D180 restart range, i.e. non-identifiability width (%, log scale)')
    ax.set_title('D180 restart-range by title (not point estimates)\n'
                  'filled bar = dominant basin; open circle = inferior sub-basin/point; gray = excluded')
    ax.grid(axis='x', which='both', alpha=0.3)
    fig.tight_layout()
    path = os.path.join(OUT, 'fig6_identifiability.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path, df


if __name__ == '__main__':
    p5 = fig5()
    p6, df6 = fig6()
    print('saved:', p5)
    print('saved:', p6)
    pd.set_option('display.width', 200)
    print(df6[['title', 'lo', 'hi', 'point', 'excluded', 'note']].to_string(index=False))
