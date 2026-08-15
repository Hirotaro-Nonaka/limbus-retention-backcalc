# v2 — What Does One More Disclosure Give an Outside Analyst?

Reproduction materials for the second preprint in this program:

> Nonaka, H. (2026). *What Does One More Disclosure Give an Outside Analyst? Empirical Identification of the Playtime Boost and Non-Identifiability of the Retention Tail from Limbus Company's Repeated DAU/MAU Disclosures.*
> - Japanese version: `paper/Limbus_Retention_v2_JA.pdf`
> - English version: `paper/Limbus_Retention_v2_EN.pdf`
> - Author ORCID: [0009-0009-6148-9974](https://orcid.org/0009-0009-6148-9974)

v2 is a single-title deep dive that builds on v1 (see the repository root). Where v1 back-calculated retention from external CCU plus a **single** terminal anchor, v2 exploits the discovery that Limbus Company disclosed platform-specific **DAU/MAU charts across five official broadcasts**, and asks experimentally: *how much does each additional disclosed data point add to an outside analyst's estimate?*

The preprint PDFs (JA/EN) and their Markdown sources are in `paper/` (the version on the public preprint server is the version of record).

## What this repository contains

| Path | Contents |
|---|---|
| `Digitization_Protocol_frozen.md` | The pre-extraction freeze: the three tasks (a)(b)(c) declared in advance, plus the change log documenting every deviation and its category (post-hoc method choice vs. data-constraint-forced). Methodological core against HARKing |
| `digitized/` | Digitized official DAU/MAU chart series (`*_extracted.csv`), the calibration notes that establish per-series reading-error floors, the usage-condition files that bind downstream use, and the extraction scripts |
| `task_a/` | **Trajectory validation** — the v1 frozen baseline's implied Steam DAU trajectory vs. the digitized official DAU chart; systematic underestimation and the constant-*h* confound (paper §3.1) |
| `task_b/` | **Empirical identification of β (playtime boost)** — r(t) regressed on the frozen event-window profile g(t), with HAC / block-bootstrap inference and the β-fixed refit (`task_b_beta_refit*.csv`); β is positively identified (paper §3.2) |
| `task_c/` | **Anchor marginal-value experiment** — anchors added in fixed chronological order, stage-by-stage; boundary-relaxation and regime-reversal diagnostics show the long-tail (D180) is **not** robustly identified (paper §3.3) |
| `task_d/` | **Theory** — closed-form stock/base collinearity that explains the non-identifiability, plus reconciliation-verification scripts (paper §4) |
| `paper/` | The preprint PDFs (JA/EN), their Markdown sources, and `sources/dart_2026_notes.md` (extraction notes for the financial figures in §5.3) |

Each task folder pairs its script(s) with its result CSVs and two Markdown files: a `*_summary.md` (what was computed) and a `*_usage_conditions.md` (the caveats and citation constraints that govern how each number may be reported). The appendix "Numerical source correspondence table" of the paper points into these files.

## Raw data (not redistributed)

Two kinds of raw input are **not** redistributed here:

- **The official chart images.** The DAU/MAU charts come from Project Moon's official broadcasts (screenshots). We do not redistribute the images; instead we provide our **digitized extractions** in `digitized/`, which are sufficient to reproduce every task-a/b/c number without the images.
- **The daily Steam CCU series** come from third-party sites (SteamDB; raijin.gg archives) whose data we do not have the right to redistribute — same as v1. To re-run from scratch, obtain daily **average** CCU and place the CSV as described in the repository root `README.md`.

The financial figures in §5.3 are sourced from a public DART filing: Project Moon Inc. (2026), 감사보고서 (Audit Report), 10th term, Financial Supervisory Service DART, filed 2026-04-01, rcpNo 20260401000553 — https://dart.fss.or.kr/dsaf001/main.do?rcpNo=20260401000553

## Environment

- Python 3.14 (`py` launcher on Windows), numpy, pandas, scipy, matplotlib (diagnostic plots), markdown + Microsoft Edge (PDF rendering only). Same stack as v1.
- Analysis code reuses the v1 model in the repository-root `scripts/` (`core.py` etc.). The v2 task scripts import `core` and hold its `BOUNDS` frozen; task-c defines `MultiAnchorFitter(core.Fitter)` to generalize the single-anchor fitter to N simultaneous anchors without editing the frozen core.

## Reproducing the main results

```bash
# from the repository root, with the v1 scripts/ on the path
py v2/task_a/run_task_a.py          # trajectory validation (§3.1)
py v2/task_b/run_task_b.py          # beta identification + beta-fixed refit (§3.2)
py v2/task_c/run_task_c.py          # anchor marginal-value stages + diagnostics (§3.3)
py v2/task_d/verify_reconciliation_E0E1_bounds.py   # theory-reconciliation checks (§4)
```

The derived CSVs already in each task folder are sufficient to verify every number in the paper without re-running anything.

## License

- Code: MIT License (see the repository-root `LICENSE`)
- Text, figures, and derived CSVs: CC BY 4.0
