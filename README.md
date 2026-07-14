# How Much Retention Can Be Estimated from Public CCU and a Single Official Anchor?

Reproduction materials for the preprint:

> Nonaka, H. (2026). *How Much Retention Can Be Estimated from Public CCU and a Single Official Anchor? Back-Calculation on Limbus Company and a Multi-Title Generalization Test.*
> - Japanese version (Jxiv): DOI to be added
> - English version (arXiv): ID to be added
> - Author ORCID: [0009-0009-6148-9974](https://orcid.org/0009-0009-6148-9974)

日本語版・英語版のプレプリント本文は `paper/` に同梱しています(公開サーバ上の版が正式版です)。

## What this repository contains

| Path | Contents |
|---|---|
| `paper/` | The preprint PDFs (JA/EN) and their Markdown sources |
| `scripts/` | All analysis code (see below) |
| `outcome/` | **Derived results** — fitted parameters, sensitivity sweeps, bound-relaxation diagnostics, bootstrap samples, event-inflow estimates, and the 6 paper figures. These are our own model outputs and are sufficient to verify every number in the paper without re-running anything |
| `data/steamdb/` | **Empty by design** — see "Raw data" below |

### Scripts

| File | Role |
|---|---|
| `core.py` | Shared model: retention kernel, pulse shapes, NNLS design, anchored Limbus fitter (paper §3.1–3.3) |
| `task_b.py` | Playtime-contamination β sweep (paper §3.5(b), §5.2) |
| `task_c.py` | 7-title generalization driver: frozen event detection, left-censoring stock column, burn-in, and the bound-relaxation diagnostics BOUNDS_RELAXED / BOUNDS_RELAXED2 (paper §3.4, §3.5(c), §6) |
| `stock_collinearity.py` | Stock/Base collinearity and binding-non-negativity diagnostics (paper §6.2) |
| `make_figs_paper.py` | Regenerates Figures 5–6 from the CSVs in `outcome/` |
| `make_pdf.py` | Renders the paper Markdown to PDF (Windows/Edge headless) |

## Raw data (not redistributed)

The raw daily CCU series come from third-party sites (SteamDB; raijin.gg archives) whose data we do not have the right to redistribute. To re-run the fits from scratch, obtain daily **average** CCU per title and place CSVs in `data/steamdb/` with the following names and columns:

- Required columns: `DateTime`, `Average Players` (UTF-8 with BOM tolerated; sub-daily rows are averaged per day by `core.load_steamdb_daily`)
- Expected filenames (see `task_c.TITLES` and `core.limbus_setup`):
  - `Limbus_steamdb_chart_1973530.csv`
  - `WT_steamdb_chart_236390.csv`, `FF14_steamdb_chart_39210.csv`, `RM_steamdb_chart_294100.csv`, `WF_steamdb_chart_230410.csv`, `DbD_steamdb_chart_381210.csv`, `PoE_steamdb_chart_238960.csv`, `PoE2_steamdb_chart_2694490.csv`

Note: `Average Players` histories are available only from 2022-09-24 onward for all titles (the left-censoring discussed in paper §3.4/§4.2). The paper's Limbus window is 2023-02-27 to 2026-07-02.

## Environment

- Python 3.14 (`py` launcher on Windows), numpy 2.5.1, pandas 3.0.3, scipy 1.18.0, matplotlib (figures), markdown + Microsoft Edge (PDF rendering only)
- Numerical note: the paper's reimplementation-baseline relative RMSE reproduces to ≈0.131 with small environment-dependent variation (±0.001 observed across pandas versions); see paper §3.5 note.

## Reproducing the main results

```bash
# 7-title fits: original frozen bounds, then the two diagnostic relaxations (paper §6)
py scripts/task_c.py --step 1 --out outcome/taskC_step1_verify.csv --restart_out outcome/taskC_step1_restarts.csv --titles WarThunder,FF14,RimWorld,Warframe,DbD
py scripts/task_c.py --step 2 --out outcome/taskC_boundscheck.csv  --restart_out outcome/taskC_boundscheck_restarts.csv --titles WarThunder,FF14,RimWorld,Warframe,DbD
py scripts/task_c.py --step 4 --out outcome/taskC_boundscheck2.csv --restart_out outcome/taskC_boundscheck2_restarts.csv --titles WarThunder,FF14,RimWorld,Warframe,DbD
py scripts/task_c.py --step 4 --out outcome/taskC_poe_boundscheck.csv --restart_out outcome/taskC_poe_boundscheck_restarts.csv --titles PoE,PoE2

# Playtime-contamination sweep (paper §5.2); fine grid used beta=0.2..0.35
py scripts/task_b.py --betas 0.2,0.225,0.25,0.275,0.3,0.325,0.35 --out outcome/taskB_beta_fine.csv --restart_out outcome/taskB_beta_fine_restarts.csv

# Stock/Base collinearity diagnostics (paper §6.2)
py scripts/stock_collinearity.py

# Figures 5-6
py scripts/make_figs_paper.py
```

## License

- Code (`scripts/`): MIT License (see `LICENSE`)
- Text and figures (`paper/`, `outcome/*.png`): CC BY 4.0
- Derived CSVs in `outcome/`: CC BY 4.0

## Citation

Until the DOIs are added above, please cite as:

> Nonaka, H. (2026). How much retention can be estimated from public CCU and a single official anchor? Back-calculation on Limbus Company and a multi-title generalization test. Preprint.
