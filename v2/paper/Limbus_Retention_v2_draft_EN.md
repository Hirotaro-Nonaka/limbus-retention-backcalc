# What Does One More Disclosure Give an Outside Analyst? Empirical Identification of the Playtime Boost and Non-Identifiability of the Retention Tail from Limbus Company's Repeated DAU/MAU Disclosures

**Author**: Hirotaro Nonaka (Riga Technical University) — ORCID: 0009-0009-6148-9974

**Keywords**: live-service games, retention estimation from external metrics, repeated DAU/MAU disclosure, empirical identification of the playtime boost (β), retention-tail non-identifiability, bound-relaxation diagnostic, sensitivity analysis

(Draft v2.0-draft-EN, 2026-07-20, translated from the JA draft. Ten sections plus appendix; finalized after independent audit (see the AI-use disclosure at the end).)

---

## Abstract

For an outside analyst who has only public CCU (concurrent users), what does one additional disclosure by the operator give to the estimation of a live-service game's retention? Using the DAU/MAU charts of Limbus Company — disclosed repeatedly across five broadcasts in official livestreams — as its material, this paper drills the external-estimation pipeline deep into a single title. By matching against these repeated disclosures, we positively identify the playtime boost β that prior work had assumed exogenously, in a form whose 95% confidence interval clearly excludes 0 (β=0.443, [0.327, 0.555]), and we derive the degeneracy of the space spanned by the two regressors as a closed form of the stock/base collinearity. Consistently with this, increasing the terminal anchors from 1 to 5 does not settle the long-term tail (D180) onto a single value; it swings among boundary-dependent regimes (6.72% under the frozen procedure, 4.53% after correcting a timepoint error — neither is cited as a definitive estimate). The informational value of a disclosure depends not on its count but on its kind.

---

## 1. Introduction

### 1.1 Background — from a single anchor to the discovery of repeated disclosure

Prior work (hereafter v1) examined whether an outside analyst without user-level telemetry can estimate the retention of a live-service game from a public concurrent-users (CCU) time series and a small number of official anchors alone. For Limbus Company (Project Moon, released February 2023), using SteamDB's daily CCU (about 3.3 years, 1,222 days) and a single official DAU/MAU measured anchor as of January 2026, it estimated — via a convolution decomposition of the same form as epidemiological back-calculation — an episode-basis D180 (the survival rate 180 days after activation) of 3.5%, and a combined band of 3.0–6.6% incorporating shape and playtime sensitivity. This estimate withstood the threefold validation of holdout, synthetic-data recovery, and independent back-calculation, but its calibration and evaluation were verified only in the vicinity of the single terminal anchor, and whether the model was consistent across the full period could not, in principle, be checked.

In the course of scrutinizing the provenance of this anchor, it was confirmed that Limbus Company had repeatedly disclosed per-platform DAU/MAU time-series charts across five broadcasts (2024-11-22, 2025-02-26, 2025-06-13, 2025-11-14, 2026-02-27) in official livestreams. What v1 used for fitting was only the single terminal value among these. This paper treats this body of repeated disclosures itself as primary material and drills the external-estimation pipeline that v1 built deep into Limbus alone. In place of generalizing across many titles, it takes as its axis how much disclosure has accumulated for a single title.

### 1.2 The question

The question of this paper is as follows. **What does one additional disclosure by the operator give to an outsider's estimation?** More concretely, it decomposes into three: (i) to what degree is a model calibrated on a single anchor consistent with the full period of the official DAU curve; (ii) can repeated disclosures move a sensitivity parameter such as the playtime boost β from an assumption to an empirical estimate; (iii) when anchor points are added one at a time, by how much does the identifiability of the long-term retention tail improve per disclosed point. Whereas v1 positioned the presence or absence of an anchor as the branch point of identifiability, this paper varies the quantity of anchors continuously and measures the marginal informational value of a disclosure.

### 1.3 Contents of this paper

This paper reports the following. First, we froze the procedure for digitizing the official DAU/MAU charts as full time series before numerical extraction, attached a lower bound on the reading error via dual extraction, and separated the stages of extraction, audit, and conditional adoption (§2.1, §2.2). Second, by matching against the official DAU curve, we show that the single-anchor-calibrated model systematically underestimates the further into the past it goes from the anchor, and that this deviation can be read as a constant-playtime confound (§3.1). Third, by regressing β from the event-window response of the ratio of official DAU to CCU, we positively identified β in a form whose 95% confidence interval clearly excludes 0 (§3.2). Fourth, we constructed a refit series that increases the terminal anchors from 1 to 5 in chronological order, and showed experimentally that D180 does not converge to a single value as anchors are added but swings among boundary-dependent regimes — that the long-term tail is not identified at this resolution of disclosure (§3.3). Fifth, we derive that this non-identification is intrinsic to the structure of the observed series, as a closed form of the stock/base collinearity (§4). The whole process was conducted under pre-registration, change-log operation, and separation of execution and audit (§2.2). In what follows, we integrate these results and discuss them from the standpoint that the informational value of a disclosure depends not on quantity but on kind (§5).

---

## 2.1 Numerical extraction from the official charts (digitization procedure)

In v1 only the single terminal point of the official DAU/MAU was used; in this paper we reconstruct it as a full time series. The details of the procedure are separated into the calibration records, and the constraints on downstream tasks into the usage-conditions files.

### 2.1.1 Material and target

The primary material is screenshots of the per-platform DAU/MAU charts repeatedly disclosed in official broadcasts. Charts of the same kind were disclosed across five broadcasts (2024-11-22, 2025-02-26, 2025-06-13, 2025-11-14, 2026-02-27), and we obtained 10 images in total for DAU and MAU. The target of this section is the two images from the 3rd Anniversary & Roadmap broadcast (2026-02-27): the DAU covers 2023-08 to 2026-01 and the MAU (after axis correction) 2024-07 to 2026-01, giving the longest observation window from a single source, which is why they were prioritized. The remaining 8 images are material for cross-broadcast checks and the anchor marginal-value experiment; in this section they are treated only as the object of pre-extraction cross-checking.

### 2.1.2 Pre-extraction freeze and change-log operation

This study froze the extraction and evaluation procedure **before** extracting numbers from the charts (the same discipline as the multi-title freezing protocol of the v1 paper). Changing the procedure or the evaluation criteria after seeing the extracted values (post-hoc method selection) is prohibited, and any subsequent deviation is recorded in the change log at the end of the freezing protocol, distinguishing "post-hoc method selection" from "forced by data constraints". The discovery of the x-axis property described below was made at the calibration stage (before series values were seen), the series-level exclusion of iOS was forced by data constraints, and the missing-value treatment of the terminal contamination zone was an exclusion of a known contamination zone based on an audit finding; all are recorded in the change log with their category made explicit, and the body of the procedure itself was not changed.

### 2.1.3 Extraction procedure

For each chart, the reference points of the x- and y-axes are taken from grid lines and axis labels, and the transformation from pixels to (date, headcount) is estimated. The pixel columns of the series color are scanned vertically to take the series center at each x-coordinate, and the effective resolution (px/day) is recorded. Where daily resolution is not obtained, we report at the resolution obtained and do not fake daily resolution by interpolation. Intervals where series separation breaks down (yellow–red crossings, anti-aliasing) are treated as missing (NaN) and not interpolated.

For estimating the reading error we use dual extraction. The same chart is extracted in two runs (run1/run2) made independent down to **the measurement method for the reference colors themselves**, in addition to the color space and the threshold algorithm, and the difference RMSE between the two is taken as the reading error. Because sharing the reference color would make the difference capture only the algorithm difference, run2's reference colors were independently re-derived from a different region and a different statistic (unsupervised k-means).

### 2.1.4 The property of the x-axis — a chart-specific discovery

At the calibration stage, the two images turned out to have different x-axis properties. In **26-01-DAU**, the pixel spacing of the five date labels is uniform (126 px), whereas the interval days are non-uniform (335/184/181/184 days), which cannot be represented by a single affine transformation. As a matter forced by data constraints, we adopted piecewise-linear interpolation. However, because the first interval (2023-08-01 to 2024-07-01, about 0.38 px/day) has a compression ratio about 1.8× that of the other intervals and the date-assignment error can reach several weeks at most, we prohibit its use for the time alignment of event windows (task b) and use only intervals 2–4 at weekly resolution.

In **26-01-MAU**, calibrating four x-axis points by the same procedure, the pixel spacing was nearly proportional to the day difference (the residual of a single straight-line fit to the four points is at most 1.12 days), and a single affine transformation held. In addition, because in this broadcast the MAU horizontal-axis labels were offset by six months due to a presenter error (the notation 25-07 is in fact 26-01), we interpret it with a +6-month shift, and internally verified it by the fact that, over the 18 months of overlap with the same broadcast's DAU, "MAU ≥ the maximum DAU within that month" had 0 violations across all 4 series. Whether piecewise-linear (DAU) or affine (MAU) is a response to chart-specific data constraints, and both were decided by residual diagnostics at the calibration stage before values were seen.

### 2.1.5 The lower-bound nature of the reading error

The run-to-run RMSE of the dual extraction was, for 26-01-DAU, STEAM 8.20%, AOS 5.33%, TOTAL 3.01%, and for 26-01-MAU, STEAM 0.92%, AOS 0.66%, TOTAL 0.29%, iOS 2.96%. These are treated as a **lower bound on the total reading error**. This is because the x- and y-axis calibration breakpoints of both runs rely on the same axis-label positions, so the axis-calibration error shared by both runs (common mode) does not appear in the run-to-run difference. In particular, the DAU STEAM 8.20% is a value whose run2 coverage was concentrated in the recent interval; we apply 8% to the full period and note that the older period may be even larger.

### 2.1.6 Quality control and asymmetric series judgments

**Missing-value treatment of the terminal contamination zone**: The terminal 2–3 pixel columns of 26-01-DAU (from 2026-01-27 onward) are a known contamination zone — an extrapolation region beyond the final axis label, overlaid with handwritten annotations, numeric label characters, and steep strokes — and were therefore set to missing across all series. Anchor cross-checking was performed mechanically at an interior point outside the contamination zone (x=785=2026-01-25), yielding STEAM −3.14%, AOS −2.88%, TOTAL −9.71% (run1). The TOTAL error is large because the official TOTAL is deduplicated and is not the simple sum of the platforms (label ratio 0.933), and because the terminal label value points to a timepoint inside the contamination zone. TOTAL is premised on not achieving exact agreement with its label value, and we attach ±10% to its absolute level.

**Asymmetric treatment of iOS**: The iOS series is excluded from quantitative use and downgraded to a qualitative mention in 26-01-DAU because its color separation from AOS is impossible (run-to-run RMSE 15.72–20.7%, about 2–5× that of the other series), whereas in 26-01-MAU it is adopted quantitatively because its run-to-run RMSE is 2.96%, on the same order as the other series. This asymmetry is an individual judgment based on the chart-specific feasibility of color separation, and when used across charts the asymmetry is stated explicitly. Because the anchor error of MAU's AOS and iOS remains in common mode (+4.23% and +7.35% respectively), we attach AOS ±5% and iOS ±8% to their absolute levels.

**Independent dual reading of anchor labels**: For the remaining 8 images (24-10/25-01/25-05/25-10 DAU/MAU), the terminal anchor values were cross-checked between two independent readings, and all 32 label values agreed (0 discrepancies). This is the material for the anchor series added chronologically in the anchor marginal-value experiment.

**Audit process**: The extraction result of each chart passes through the flow of extraction → audit → conditional adoption. 26-01-DAU was adopted after 2 audit rounds (conditional adoption, 5 items returned → adopted on re-audit), and 26-01-MAU after 3 returns (integration of iOS into the main CSV, verification of AOS halo contamination, visual confirmation of the AOS>Steam crossing interval). The downstream constraints are finalized in each usage-conditions file.

## 2.2 Verification design — pre-declaration and change-log operation

### 2.2.1 Pre-declaration of the three tasks

The uses of the digitized official time series were finalized into three tasks in the freezing protocol before beginning numerical extraction from the charts. (a) **Trajectory verification**: full-period comparison of the v1 model's DAU trajectory with the official DAU curve (v1 had only the single terminal point). (b) **Empirical identification of β**: estimate the playtime boost β from the event-window response of r(t) = CCU(t)·24/(h·DAU_official(t)) − 1. (c) **Anchor marginal-value experiment**: with a refit series that increases the anchor points from 1 to 5, quantify the increment in identifiability per disclosed point.

The evaluation methods were declared at the same time. (a) is a descriptive report of the relative RMSE (same definition as v1), and instead of setting a pass/fail threshold in advance, we note the v1 holdout DAU error of −5.4% as a reference point. (b) fixes the treatment of h to report both "fixed at 2.2" and "freely estimated", and states explicitly that if the event-day sample cannot separate them due to insufficient resolution, "β cannot be identified at this resolution" is accepted as the conclusion (the procedure is not changed until it can be identified). (c) fixes the order of anchor addition to chronological (24-10 → 25-01 → 25-05 → 25-10 → 26-01), pre-enumerates four quantities to report at each stage — the number of boundary-pinned parameters, the bootstrap CI width of D180, the lower bound of τ, and the state of base/stock — and then applies v1's bound-relaxation diagnostic at each stage. The known risks (the possibility that the effective resolution of ~0.8 px/day yields no daily resolution, the breakdown of series separation, drafting errors in the official chart itself) are also pre-disclosed in the same document. This configuration, which leaves no room to select the evaluation criteria after seeing the extracted values, is the countermeasure against post-hoc hypothesis construction (HARKing) in v2.

### 2.2.2 Track record of change-log operation

All post-freeze deviations are recorded in the change log at the end of the protocol, distinguishing "post-hoc method selection" from "forced by data constraints". The change log has 5 entries in total, including records from the calibration and extraction stages (the piecewise-linear treatment, the terminal missing-value treatment, the iOS exclusion, etc., described in §2.1); here we cite two of a different nature, bearing on the verification design. First, the case where the affine-transformation premise at freezing turned out not to hold for the x-axis of 26-01-DAU (see §2.1) and piecewise-linear interpolation was adopted. This was a discovery at the calibration stage before extraction of series values, was categorized as forced by data constraints, and its attendant usage restriction (prohibiting use of the first interval's dates for the event-window alignment of task (b)) was recorded at the same time. Second, the case where the pre-declared bound-relaxation diagnostic was not performed in the first execution of task (c). This was a deviation from the declaration due to an error in the delegation instructions; the audit detected it, and it was corrected by performing the relaxation diagnostic on stages 3–5 before adoption was finalized (see §3.3). What this operation reports is not zero deviations, but that the categorized recording and correction of deviations functioned.

### 2.2.3 Audit workflow

The numerical deliverable of each task is adopted only after passing a two-stage process that separates execution from audit. Execution and audit are carried out by multiple independent, separate-lineage models, enabling cross-model verification (see the AI-use disclosure at the end). In addition to numerical cross-checking by independent recomputation, the audit makes mandatory the detection of boundary pinning, degenerate estimates, and the apparent stability of the "small restart spread = stable" type, and adoption is finalized as conditional adoption with expression constraints (the usage-conditions file). All descriptions in the results sections of this paper (3.1–3.3) are written under the constraints of the usage conditions finalized in this audit.

---

## 3.1 Trajectory verification (task a) — from a single terminal point to full-period comparison

### 3.1.1 Setup

v1's model evaluation was limited to matching the single terminal point (the January 2026 anchor). Digitization makes possible the full-period matching (2023-08 to 2026-01) of the official DAU curve, which was impossible in v1. We first refit v1's frozen baseline (both BOUNDS and initial values unchanged) and obtained a relative RMSE of 0.130457. This nearly coincides with the lower end of the known environment-difference band (≈0.1305–0.1314) (undershooting by 0.00003, at the numerical-error level), confirming reproduction of the baseline. In this fit, both c=1.0 and τ=6000 are pinned at the BOUNDS upper limits, which is a known behavior already recorded in v1 (not newly arising in this task). The model DAU is the product of the design matrix and the NNLS coefficients X(t)@a (β=0, Steam alone), and the official side uses the main series STEAM run1 (n=509). The evaluation is descriptive as pre-declared (no pass/fail threshold).

### 3.1.2 Results — systematic underestimation and its shrinkage near the anchor

| Interval | Relative RMSE | Bias (model − official) |
|---|---|---|
| Full period | 28.34% | −16.68% |
| Interval 1 (23-08 to 24-07) | 42.40% (low confidence) | −39.05% |
| Interval 2 (to 25-01) | 24.14% | −16.00% |
| Interval 3 (to 25-07) | 22.57% | −7.77% |
| Interval 4 (to 26-01) | 19.59% | −5.41% |

The bias is negative (model < official) in every interval, and its magnitude is larger the further into the past it goes from the anchor, shrinking monotonically from −16.0% to −5.4% across intervals 2→4. Interval 1 (−39%) is a low-confidence value whose date-assignment error reaches up to several weeks, and the conclusion (systematic underestimation, shrinking near the anchor) holds only for intervals 2–4. The same direction and the same temporal pattern are reproduced in STEAM run2 (reference, −34.3% → −0.1% across intervals 1→4), and the shape relative RMSE of the TOTAL series after level matching (24.80%) is of the same magnitude as STEAM run1's full period (28.34%), supporting that the pattern is not STEAM-specific reading noise.

The magnitude of the deviation cannot be explained by reading error. Against the STEAM reading-error floor of 8.20% (treated as a lower bound) from dual extraction, the full-period relative RMSE is about 3.5× it. Only the bias of interval 4 (−5.41%) falls below the floor, but the relative RMSE of that interval (19.59%) itself exceeds the floor. Note that, because the digitized Steam series carries an extraction bias of −3.14% by label ratio, the reported underestimation magnitude is a lower bound relative to the true value. Also, at the terminal the model equals to slightly exceeds the digitization (model/official=1.101 at 2026-01-21, because the model is constrained to the label anchor), and the interval-4 mean of −5.4% is the result of the early-interval underestimation outweighing the terminal excess.

### 3.1.3 Interpretation — the deviation as a constant-h confound

This deviation pattern should be read not as a defect of the retention model but from the structure of the CCU→DAU conversion. In the β=0 baseline, model_dau = fitted CCU × 24/h (h=2.128, constant over the full period), so the comparison is effectively "CCU shape vs. official DAU shape". A bias that is negative over the full period and larger further in the past is consistent with a single constant h being unable to reconcile the CCU⇄DAU conversion across the full period — that is, with a time trend in per-person playtime (the domain that β, task b, §3.2, handles). That the retention parameters are already pinned at their upper limits at c=1.0, τ=6000 (retention set to the longest within the search range still fails to reach the initial DAU) reinforces this reading.

The implication for v1 is as follows. The interval-4 in-sample bias of −5.41% (digitization basis) and the v1 holdout error of −5.4% (official-label basis, out-of-sample) are quantities that differ doubly in definition, but they are of the same order (about −5%) and are consistent with the anchor being in interval 4. That is, v1's good terminal match can be explained not as evidence that the model was calibrated over the full period, but as a consequence of the evaluation point being near the anchor. A statement of the fit adequacy of a single-anchor-calibrated model becomes an overestimate unless the temporal distance from the anchor is made explicit.

## 3.2 Empirical identification of β (task b) — the playtime boost is positively identified

### 3.2.1 Setup and identifiability

v1 treated the playtime boost β (the degree to which per-person playtime increases at events) as an exogenous sensitivity parameter, and assumed β=0.3 as a "plausible cap" to set the upper limit of the D180 combined band. With the official DAU in hand, β can be estimated directly from data rather than assumed. We regress the ratio r(t) = CCU(t)·24/(h·DAU_official(t)) − 1 on the frozen event-window profile g(t) (the exp3 pulse of core.py, unchanged) with an intercept. Following the change-log constraint, the regression is performed on intervals 2–4, excluding the first interval (n=387), and the treatment of h reports both "fixed at 2.2" and "freely estimated" as pre-declared.

The pre-disclosed risk of "non-identifiability due to insufficient resolution" did not materialize. The measured resolution of intervals 2–4 is a median of 1.0 day/point; all 31 relevant event windows contain 2 or more digitized points (0 events with 0 or 1 points), which sufficed to follow the event response.

### 3.2.2 Estimation results — the CI clearly excludes 0

| Variant | β point estimate | 95% CI (event-window block bootstrap, main report) |
|---|---|---|
| h=2.2 fixed (declared) | 0.443 | [0.327, 0.555] |
| h freely estimated (h_hat=1.800) | 0.541 | [0.391, 0.682] |

The regression residuals have strong serial correlation (lag-1 autocorrelation 0.71), and the CI was built with a block bootstrap rather than the iid method that ignores this. The main method uses blocks with the event window as the unit (K=28, faithful to the dependence structure because the source of the serial correlation is tied to the event-pulse process itself), and the sub-method uses fixed 14-day blocks (K=41); in all four combinations of 2 block methods × 2 h variants, the 95% CI clearly excludes 0. With independent Newey–West HAC standard errors as well, se is flat over lags 14–60 days and the CI ≈ [0.33, 0.55] is reproduced, with no dependence on the block choice. The reading-error floor of 8.20% was added to each iteration as multiplicative noise independent across all points (a conservative choice that overestimates the common-mode component). **β is positively identified** — this is the definitive conclusion obtained from a single official-disclosure series at this resolution.

We also confirmed the robustness of the estimate. The individual per-interval regressions give β=0.32/0.55/0.50, and the fixed-effects model controlling for the secular trend in r(t) with interval dummies gives 0.458; all are of the same order as the main estimate 0.443, so the level trend in r(t) (median −0.41 → −0.09 across intervals 1→4, a descriptive statistic consistent with the constant-h confound finding of task a) is not a confound that artificially inflated β. Note that the h_hat=1.800 from free estimation is an aggregate value including the time trend, and we do not emphasize it as a point estimate of a constant h.

### 3.2.3 Bias direction and implications for the v1 band

The precise level of the point estimate carries residual bidirectional systematic error. In the upward direction there are the reading error of event correlation (run-to-run difference: 7.40% inside windows vs. 6.40% outside) and the resolution asymmetry between CCU (daily) and DAU (≈1.45 day/px), which the block bootstrap does not correct. In the downward direction, the date error of the g(t) boundary (±1.45 day/px) attenuates the coefficient toward 0 as measurement error in the explanatory variable. The net direction of the bias is indeterminate. Therefore, regarding the relation to v1's assumed cap, from the fact that the CI lower bound (0.327–0.391) exceeds or is comparable to 0.3, we state only that "external estimation suggests above 0.3, but an upward bias due to event-reading confounds is possible and it cannot be asserted definitively".

The implication for v1's D180 combined band [3.0%, 6.6%] is, contrary to the prior expectation (compression of the band), in the direction of an upward revision. We verified this not by linear interpolation of v1's sensitivity table but by a **formal refit with β fixed at the identified value** (all parameters refit at β=0.327/0.443/0.555; `v2/task_b/task_b_beta_refit.csv`). Because this task's β is the coefficient of the r(t) regression and uses the same quantity and the same event-window profile g as the β of the model's CCU=DAU·(h/24)·(1+β·g), the unit-correspondence premise that remained in the interpolation does not arise in this refit. The result resolves into two layers. **(i)** The ordering/direction — that D180 exceeds both β=0 (0.033) and v1's band (3.0–6.6%) — is robust. **(ii)** But the level is not pinned to a point: over β=0.33–0.56, D180 moves within a band of about 7.6–9.5%, and over this interval the short-term scale λ pins to the search upper bound (30 days) and becomes non-identified (relaxing the bound moves λ to 36–44 days and shifts D180 a further +3–4%, with relative RMSE essentially unchanged). That fixing β still leaves the tail level in a band via the non-identification of λ is a re-confirmation, from the β side, of the finding in §3.3 that the tail is not robustly identified. We therefore report this result as **"the direction of the upward revision is settled; the level is a band (no point estimate, λ non-identified)"** and do not circulate it as a single point estimate.

## 3.3 Anchor marginal-value experiment (task c) — the tail is not identified

### 3.3.1 Design

To measure the marginal informational value of one additional disclosure, we constructed a refit series adding the terminal anchors of the 5 broadcasts one at a time in chronological order (24-10 → 25-01 → 25-05 → 25-10 → 26-01). The four quantities to report at each stage (the number of boundary-pinned parameters, the bootstrap CI width of D180, the lower bound of τ, and the state of base/stock) were pre-enumerated (§2.2). The anchor values were cross-checked by independent dual reading of 32 labels across 8 charts, with 0 discrepancies. The anchor mask period was unified to "the calendar month the broadcast label points to" across all 5 stages (a replication of core.py's existing 26-01 anchor convention, to avoid a mixed treatment that would use precise terminal dates for only some). The model proper (BOUNDS, kernel, optimization) was inherited unchanged from v1's frozen core.py, with only the simultaneous injection of multiple anchors added as an extension.

### 3.3.2 Stage results — regime flip, not convergence

| Stage | Anchor | D180 | Boundary pinning | D180 90% CI width* |
|---|---|---|---|---|
| 1 | 24-10 | 4.86% | c@1.0, τ@upper 6000 | 1.52pp |
| 2 | +25-01 | 4.23% | same | 1.00pp |
| 3 | +25-05 | 3.67% | same | 1.19pp |
| 4 | +25-10 | 6.73% | k@1.5, c@1.0, τ@lower 200 | 1.37pp |
| 5 | +26-01 | 6.72% | same | 1.04pp |

\* The CI widths (n_boot=100, stage 5 with an independent seed) are for qualitative use only: the only thing to read is that "the CI width does not shrink monotonically as anchors are added", and no quantitative comparison of small differences is made.

The D180 point estimate falls across stages 1→3 (4.86% → 4.23% → 3.67%), then jumps discontinuously to 6.73% with the addition of the 4th point (25-10), and stays at 6.72% even at the 5th. At the same time the structure of the boundary pinning changes: in stages 1–3, τ is pinned at the upper limit of 6000 days (a regime in which no decay is visible within the observation period), whereas in stages 4–5, τ flips to the lower limit of 200 days and k is newly pinned at the upper limit of 1.5. In the multi-start diagnostic, the extremely small restart spread in stages 1–3 is the result of all initial values converging to the same boundary corner and is not adopted as evidence of identifiability (apparent stability due to boundary pinning). In stages 4–5 the restart spread clearly widens (relative objective difference 3.5% and 2.3%), changing to a structure in which multiple local solutions compete. base degenerated to 0 in all stages. Note that the structural-decomposition values (λ/k/p/c/τ) of stages 4–5 lie under simultaneous pinning of 3 parameters and are not trusted as a physical description of retention.

### 3.3.3 Diagnostics — disaggregating the cause of the relaxation and the flip

In the pre-declared bound-relaxation diagnostic (τ lower limit 200→30 days, k upper limit 1.5→3.0; the upper limit 1.0 of c is left in place as a physical upper limit), stage 3 reproduced both τ and k almost perfectly at their original values (D180 difference nearly 0%), confirming that the τ upper-limit pinning is not an appearance of insufficient search range. In stages 4–5, on the other hand, τ did not pin at the lower limit of 30 days but moved to interior values of 131–134 days (the original τ=200 was a product of insufficient search range), and k stayed pinned at the new upper limit of 3.0 (a runaway parameter that wants the upper limit by exactly as much as it is loosened). D180 falls from 6.7% to about 6.1%, but the jump itself between stages 3→4 remains after relaxation.

In the disaggregation of the cause of the flip (sensitivity analysis), making the collaboration-event window 28 days had no effect (D180 difference 0.02–0.03 pt). What mattered was the mask month of the 25-10 anchor: the 25-10 chart explicitly labels the terminal 2025-11-09, and because the anchor value 418,242 is the value of November 9, the October mask assigned by the unified convention is, for this anchor, an empirical timepoint error. Correcting the mask to November completely resolves the flip in stage 4 (τ returns toward the upper limit, D180=4.01%), and partially resolves it in stage 5 (the k pinning disappears but τ remains at the lower limit of 200 days, D180=4.53%). The correction of the timepoint error is not post-hoc method selection, but we also maintain transparent reporting of the value under the frozen procedure — that is, we always report both D180 = 6.72% (frozen procedure, contaminated by a known timepoint error) and 4.53% (after correction, more physically plausible), and neither is cited as a definitive estimate. The combination of "November mask × bound relaxation" has not been run, and whether the short τ (τ=200) remaining in stage 5 is a genuine survivor or a further boundary artifact has not been disaggregated — this cannot be asserted definitively here.

### 3.3.4 Conclusion — the marginal informational value for tail identification is approximately nil

The robust conclusion is neither 6.72% nor 4.53%, but "with CCU + anchors, D180 (the long-term retention tail) is not robustly identified". Adding anchors does not converge D180 to a single value; it merely selects one of the boundary-dependent regimes in which the choice of a single mask month swings D180 across 4.0–6.7%, which is an exposure of regime multiplicity. From the standpoint of the marginal informational value of an additional anchor, tail identification is unachieved even at 5 points, and the apparent stability of stages 1–3 is a fake produced by corner convergence. Under the frozen procedure, the addition of the 4th point induced a regime flip, but this is mainly an artifact of the mask timepoint error (§3.3.3), and even in the corrected series (3.67% → 4.01% → 4.53%) the tail remains unidentified. Hence the marginal informational value per additional disclosed point can be summarized, with respect to tail identification, as approximately nil. The mechanistic explanation of this experimental result (why anchors cannot constrain the tail) is left to the closed-form discussion of §4.

---

## 4. Theory — explaining tail non-identification via the stock/base collinearity

The anchor marginal-value experiment of §3.3 showed that even increasing anchors to 5 points does not converge D180 (the long-term tail) to a single value but makes it swing among boundary-dependent regimes. This section shows in closed form that this non-identification is not a malfunction of the fitting but an unidentifiability intrinsic to CCU-only observation. The material is the defining equations of the retention kernel and the design matrix (extracted from the implementation and audit-confirmed); below, using as an example the stock column (the left-censoring correction column) that does not appear in the Limbus-alone baseline but belongs to the same optimization system, we trace the degeneracy of the space spanned by the two regressors.

### 4.1 Basic identity — the base column asymptotes to an affine transform of the stock column

Write the retention kernel as R(d) = R_fast(d) + p·c·e^{−d/τ}. R_fast is the short-term component (time constant sufficiently smaller than τ), and the second term is the long-term component. The base column of the design matrix is the cumulative response of a constant inflow x_base(t) = Σ_{i=0}^{t} R(i), and the stock column is the decay of the population already present at the start of the observation window x_stock(t) = e^{−t/τ} (τ shared with the long-term component). The discrete sum of the long-term component closes as a geometric series,

  Σ_{i=0}^{t} p·c·e^{−i/τ} = G·(1 − e^{−(t+1)/τ}),  G ≡ p·c/(1 − e^{−1/τ}) ≈ p·c·τ (τ≫1)

The cumulative of the short-term component F(t) = Σ_{i=0}^{t} R_fast(i) saturates to a constant F_∞ once t exceeds several times the short-term time constant. Therefore, on the evaluation window W = [t₀, T] after burn-in (t ≥ t₀ = 90),

  x_base(t) = K − G′·x_stock(t) + O(ε),  K ≡ F_∞ + G,  G′ ≡ G·e^{−1/τ},  ε = |F(t)−F_∞|

holds. That is, on window W the space spanned by {x_base, x_stock} coincides with the space spanned by {1 (constant), e^{−t/τ}}, and the coefficient correspondence (a₀, a_s) → (a₀K, a_s − a₀G′) is always invertible as a triangular map. Therefore **the unidentifiability reduces entirely to "whether 1 and e^{−t/τ} can be distinguished on the window"**.

### 4.2 The two collapse modes

Evaluating u(t) = e^{−t/τ} on window W, identification breaks in both of the extreme τ bands.

**(i) When τ is large relative to the window length (τ ≫ T−t₀)**: u(t) ≈ e^{−t₀/τ}·(1 − (t−t₀)/τ) can be linearized, and the relative variation from the constant is CV(u) ≈ (T−t₀)/(τ·√12). Substituting, as an example, T−t₀ ≈ 1132 days and τ=6000 days (the upper-limit-pinned regime observed in v1 and v2 task c), CV ≈ 5.4%, so u becomes nearly constant and {1, u} numerically rank-drops. Then only the 1 degree of freedom of the combination a₀K + (a_s − a₀G′)·ū is identified, and a₀ and a_s cannot be separated. **(ii) When τ is small relative to the burn-in (τ ≲ t₀/ln(1/δ))**: on the window u(t) ≤ e^{−t₀/τ} ≈ 0 (e.g., ≈0.05 at τ=30 days), the stock column itself becomes a numerically zero column and a_s does not contribute to the objective. The band of identifiable τ is limited to the intermediate region sandwiched between these two modes (of order t₀ ≲ τ ≲ (T−t₀)/several times).

### 4.3 The non-negativity constraint of NNLS forces a corner solution

Under the near-rank-drop of mode (i), the objective is nearly flat along the degenerate direction of the a₀–a_s plane (the ridge where a₀K + (a_s−a₀G′)ū = constant). Under the non-negativity constraint a₀, a_s ≥ 0, when noise slightly tilts this ridge, the optimal solution falls at the intersection of the ridge and the constraint boundary, i.e., at a corner (a_s=0 or a₀=0). Which corner it falls into is determined by the sign of the tiny correlation between the residual and the u direction — effectively a coin toss of noise. Hence, that the coefficient on the degenerate side is 0 does not mean "that quantity is actually zero", and **the reproduction stability of the degeneracy (falling into the same corner every restart) is not evidence of identification** — because within the same data the sign of the ridge tilt is fixed.

### 4.4 Connection to §3.3 — why anchors cannot constrain the tail

This geometry explains the two phenomena observed in §3.3 as the same mechanism. First, the "apparent stability" in which multiple initial values fell into the same boundary corner (τ upper limit 6000, base=0) and the restart spread became extremely small in task c stages 1–3 is isomorphic to the corner forcing of §4.3 — the smallness of the spread is a consequence not of identifiability but of the flatness of the degenerate direction. Second, adding an anchor gives one constraint on the absolute level, but what is inseparable in mode (i) is the degrees of freedom that determine the shape of the long-term tail (a₀ and a_s, and thence the shape of τ). The constraint given by a level anchor lies along the already-identified combination (the level a₀K + (a_s−a₀G′)ū, nearly constant on the window), and has almost no component in the flat degenerate direction of the objective — the component in the degenerate direction is K·(ū − u(t_a)) = O(CV) (t_a the anchor timepoint), and is exactly zero when u(t_a)=ū. The difference in tail shape becomes manifest outside the observation window (t > T), and an in-window level anchor cannot separate it. Therefore, however many level anchors are added, as long as the τ band is on the mode-(i) side, the shape degrees of freedom of the tail are not constrained. That left-censoring correction can be identified from CCU alone only in the intermediate τ band is a necessary condition, not a sufficient condition of "identifiable if there is an anchor" — indeed task c could not identify the tail even with 5 anchors.

### 4.5 Limitations

This section's reduction has several premises, which limit the conclusion to "a closed-form explanation of the mechanism". First, the F-saturation assumption of §4.1 (t₀ ≫ short-term time constant) holds for Limbus-class parameters, but because an actual stock-target title may have λ pinned at its upper limit with R_fast unsaturated at t₀=90, the two-column affine reduction becomes a coarser approximation. Second, because each event column also has a τ tail through convolution with R, the constant and slow-decay directions are shared by the base/stock/event-column tails, so the actual quasi-null space can exceed two dimensions — the clean two-dimensional corner picture of §4.3 is a lower bound on the degeneracy. Third, because the implementation's NNLS is solved on a weighted design (dividing rows by the CCU that increases over time), the collinearity and the CV should strictly be evaluated in the weighted column space, and this section's unweighted CV is an approximation of that. The threshold of CV(u) (below what percentage it is practically rank-dropped) depends on the noise level, and this section gives no normative value.

---

## 5. Discussion

### 5.1 Integrating the three-piece set — the informational value of a disclosure depends not on quantity but on kind

The three results of this paper show a structure that cannot be captured by a single quantitative scale, regarding what an outsider gets per disclosed point. First, a single terminal anchor constrains its own vicinity: the good fit v1 obtained from the single terminal point is not evidence that the model was calibrated over the full period, but a consequence of the evaluation point being near the anchor; matched against the official DAU curve, systematic underestimation appears the further into the past it goes from the anchor, with the bias shrinking monotonically from −16.0% to −5.4% across intervals 2→4 (§3.1). Second, the full time series of repeated disclosures makes a sensitivity parameter that had been merely an assumption identifiable from data: the playtime boost β was positively identified in a form whose 95% confidence interval clearly excludes 0 (β=0.443 [0.327, 0.555]). Third, even so, the long-term tail is not identified: even increasing the terminal anchors to 5 points does not settle D180 onto a single value; it merely selects one of the boundary-dependent regimes in which the choice of a single mask month swings D180 across 4.0–6.7% (reporting both 6.72% under the frozen procedure and 4.53% after correcting the timepoint error), and the marginal informational value for tail identification is approximately nil (§3.3).

The axis running through these is that the informational value of a disclosure depends not on its quantity but on its kind. Even for the same "one disclosed point", the information that constrains the near-term level, the information that identifies the slope of the event response, and the information that constrains the shape of the long-term tail are different things, and however many terminal level anchors are stacked — because that constraint has almost no component in the degenerate direction, as the closed form of §4 shows — the shape degrees of freedom of the tail are not constrained. That β can be identified while the tail cannot is not a contradiction but a difference in the kind of information required: the former has information in the response on the short time scale of the event window, while the latter requires information in the shape of a time constant exceeding the observation window.

### 5.2 Re-evaluation of v1

This integration works not to negate v1's conclusions but to make its range of validity explicit. v1's D180 estimate (episode-basis 3.5%, combined band 3.0–6.6%) remains valid as an estimate that withstood the threefold validation near the terminal anchor. However, this paper's full-period comparison shows that the goodness of its fit cannot be read as evidence of full-period calibration: the interval-4 in-sample bias (−5.41%, digitization basis) and the v1 holdout error (−5.4%, official-label basis, out-of-sample) are quantities that differ doubly in definition but are of the same order, and both are consistent with the anchor being in interval 4. A statement of the fit adequacy of a single-anchor-calibrated model becomes an overestimate unless the temporal distance from the anchor is made explicit. For β as well, that the CI lower bound (0.327–0.391) exceeds or is comparable to v1's assumed cap of 0.3 is material suggesting that external estimation is above 0.3, but an upward bias due to event-reading confounds is possible and it cannot be asserted definitively. The implication for the v1 band is an upward revision, opposite to the prior expectation (compression). Under a formal refit with β fixed at the identified value (§3.2.3), the direction — that D180 exceeds both β=0 and v1's band — is robust, but the level stays in a band (about 7.6–9.5%) and cannot be pinned to a point, because fixing β still leaves the short-term scale λ non-identified at the search bound (isomorphic to the terminal-anchor finding of §3.3). The upward revision is settled as a direction, but sharpening the tail level is blocked by the non-identification of λ; it is not treated as an independent point estimate. On the whole, this paper is deflationary toward v1: it does not replace v1's figures but narrows its range of assertibility.

The downstream consequences of β extend beyond v1's D180 combined band to the three-layer person-basis figures of v1 §5.4. First, the permanent-attribution rate that v1 reported as "person-basis D180 ≈ 27% (20–39%)" is, when its derivation is traced, not a tail value of the survival function but effectively a snapshot ratio of an attributed population on the order of recent MAU divided by an estimate of lifetime unique installs (confirmable from the fact that, across the 3 cases in the reported table, the product of uniques and D180 — i.e., the numerator — is nearly constant). The reported band 20–39% propagates only the denominator (AppBrain's cumulative downloads ±30%) and does not include numerator-side uncertainty. Moreover, the basis for this quantity being "permanent" — the decay time constant of the monthly attribution kernel — is a boundary solution pinned to its lower limit, the same type of pattern this paper showed in task c: that "convergence agreement due to boundary pinning is not evidence of identification, and that lower bound can be an artifact of insufficient search range." Therefore the figure 27% itself does not move, but its designation should be read not as a survival-function value but as a snapshot ratio of "recent MAU to an estimate of lifetime uniques," and the "permanent" shape claim is unverified by this paper's diagnostic criteria.

Second, v1's "about 3.9 activations per person (3.0–5.5)" is the pure arithmetic of total activations (about 10.76 million) divided by lifetime uniques (≈2.78 million), and this total was obtained under β=0. Since this paper identified β=0.443 at a level exceeding v1's assumed cap of 0.3, the total activations are on the high side, and hence episodes per person is an upper-side value — the correction under β=0.443 is not computed in this paper, but its direction is a decrease below 3.9. However, the fact that total activations greatly exceed uniques (growth is driven by repeated return rather than new acquisition) is preserved after correction, and this is a decline in quantitative confidence rather than a retraction of the conclusion. The person-basis level (27%) requires only re-reading and involves no recomputation; what requires a formal re-derivation is episodes per person, which we defer to future work (§6).

### 5.3 Implications for disclosure practice (a qualitative discussion)

The above analysis depends on the fact that Limbus Company has repeatedly disclosed DAU/MAU across five broadcasts. As to why this level of disclosure is possible, we record one explanation as a qualitative, contextual discussion. **The discussion of this section is not a validated causal claim, and has no effect whatsoever on the results (§3), the method (§2), or the model.** According to the audit report filed with the Korean DART electronic-disclosure system (Project Moon Inc., 감사보고서 [audit report], 10th fiscal year, filed 2026-04-01, rcpNo 20260401000553), Project Moon's issued shares as of the end of the 10th fiscal year (2025-12-31) were composed solely of the representative, one individual, and treasury holdings, with all preferred shares converted into treasury shares, and there is effectively no external institutional-investor stake. Also, advertising expense is about 0.06% of operating revenue (about 0.10% in 2024), a structure in which customer acquisition does not depend on the advertising market. These two points work in the direction of reducing the cost of publicly disclosing user metrics rebounding on the company via the capital and advertising markets. That is, the explanation holds that a capital and cost structure in which external investors are nearly absent and dependence on advertising exposure is small minimizes the exposure cost of repeatedly disclosing metrics such as DAU/MAU. This is merely one qualitative interpretation, and has not been validated as a causal explanation of disclosure behavior. We do not enter into facts not stated in the report, such as the circumstances of the preferred-share buyback or the identification of former shareholders.

### 5.4 Summary

What the three-piece set shows is that the value of a disclosure to an outsider is determined not by the total number of disclosed points but by the alignment between the latent quantity one wants to constrain and the time scale that its information reaches. The near-term level, the short-scale response, and the long-term tail require different kinds of disclosure. This standpoint shows the necessity, when reporting estimation from external metrics, of making explicit not "how many anchors there are" but "which latent quantity is constrained by which disclosure".

---

## 6. Limitations and future work

**Undigested disclosure material.** What this paper digitized is the 2 DAU/MAU images of the 3rd Anniversary broadcast (2026-02-27); of the 10 images obtained across 5 broadcasts, the remaining 8 were used only as material for the anchor marginal-value experiment after independent dual reading of their terminal anchor values (all 32 labels agreeing), and full time-series extraction has not been reached. External verification of the reading error by cross-broadcast checks, and mutual verification of the calibration using the overlap intervals of multiple broadcasts, become possible in future extraction.

**Unit correspondence and refit of β (performed).** The β of §3.2 is the coefficient of the r(t) regression and uses the same quantity and the same g as the β of the model's (1+β·g). Hence the fixed-β formal refit (§3.2.3, `v2/task_b/task_b_beta_refit.csv`) does not incur the unit-correspondence premise (former TB-4) that remained in the interpolation. The refit finds that the direction of the upward revision beyond v1's band is robust, but the level of D180 stays in a band of about 7.6–9.5% over β=0.33–0.56, and over this interval the short-term scale λ pins to the search bound and becomes non-identified, so no point estimate is possible. Two limitations remain: (a) the multi-start coverage is small at each β, and sharpening the ends of the band requires more starting points; and (b) the level of β's point estimate itself carries residual bidirectional systematic error — an upward bias from event-correlation reading error and the CCU/DAU resolution asymmetry, and a downward bias from the g(t)-boundary date error (§3.2.3).

**Temporal alignment of event windows.** This paper's event windows are based on in-game run periods. However, the influence of events can extend outside the run period through prior announcements and previews, and through the regular recognition that periodic hosting creates in fans. This kind of timing mismatch cannot be fully captured by this paper's window design and remains as residual uncertainty of the event-response-based estimation (§3.2).

**Granularity of the uncertainty assessment.** The bootstrap is n_boot=100, and the CI width of D180 is limited to qualitative use (only the point that it does not shrink monotonically as anchors are added). Also, in the disaggregation of the cause of the flip in §3.3, the combination "November mask correction × bound relaxation" has not been run, and whether the short τ (τ=200) remaining in stage 5 is a genuine survivor or a further boundary artifact has not been disaggregated — this cannot be asserted definitively here.

**Unexplained residual of the original reconciliation.** Of the difference between the relative RMSE 0.133 of v1's original analysis and the 0.130457 of the current core.py, about 74% could be explained by the difference in the weighting structure of the objective, but the remaining about 26% (0.000662) has not been pinned down. The leading candidate cause is that the weekday-adjustment script no longer exists in the original archive and cannot be exactly reproduced, but it cannot be identified with the existing material.

**Person-basis metrics out of diagnostic scope.** This paper's diagnostics (trajectory verification, β identification, anchor marginal value, collinearity) target v1's episode-basis D180 and the inflow/retention structure, and do not touch the three-layer person-basis figures of v1 §5.4 (permanent-attribution rate ≈27%, episodes per person ≈3.9, monthly attribution kernel). As stated in §5.2, the person-basis level 27% requires re-reading as a snapshot ratio (the figure is unchanged), and episodes per person 3.9 becomes an upper-side value under β=0.443. The correction of the latter — the re-derivation of total activations and episodes per person under β=0.443 — is beyond the scope of this paper and is left as future work (v2.1).

**Restriction to a single title and single genre.** All results of this paper are based on a single biweekly-updated gacha-type live-service title. The very existence of repeated disclosure is exceptional (in v1's multi-title verification, no title had an official anchor), and how far this paper's findings on the marginal informational value of disclosure generalize to other operating rhythms and other genres cannot be verified until a title with similar repeated disclosure is found.

---

## 7. Conclusion

When an outsider estimates the retention of a live-service game from public CCU and official anchors, the marginal informational value of one additional disclosure differs greatly depending on the kind of latent quantity one wants to constrain. Using the repeated disclosures of the single title Limbus Company as material, this paper obtained three conclusions. First, a single terminal anchor constrains only its near-term level, and when matched against the official DAU curve over the full period, systematic underestimation appears the further into the past it goes from the anchor. Second, the time series of repeated disclosures makes the playtime boost β — which had been merely an assumption — positively identifiable from data (β=0.443, [0.327, 0.555]). Third, even so, the long-term tail (D180) is not identified even when the terminal anchors are increased to 5 points, but swings among boundary-dependent regimes (reporting both 6.72% and the corrected 4.53%). This non-identification is explained as a closed form of the stock/base collinearity intrinsic to CCU-only observation. Running through these is the standpoint that the informational value of a disclosure depends not on quantity but on kind. When reporting estimation from external metrics, one must make explicit not the number of anchors but which latent quantity is constrained by which disclosure. All results were obtained under a procedure frozen before extraction, with execution and audit separated.

---

## Disclosure of AI use

The conception, design, and judgments of this research are the author's. Large language models were used for numerical computation (model estimation and the like) and drafting assistance, and cross-verification was performed with multiple models. The author takes full responsibility for all contents.

---

## Appendix: Numerical source correspondence table (by section)

The provenance of the numbers and facts appearing in each section is shown, with section headings, in the order they appear in the section. Paths in the "File" column are relative to the public reproduction repository for this program (https://github.com/Hirotaro-Nonaka/limbus-retention-backcalc ; the v2 materials are under the `v2/` subtree). Items that pointed to section drafts of this paper itself have been replaced with the final section number after integration. Primary sources that cannot be redistributed (official broadcast chart images, financial filings, etc.) are cited via the extracted data or bibliographic record included in the repository.

### 1. Introduction

Provenance of the numbers and facts appearing in this section (number → file → location).

| Number / fact | File | Location |
|---|---|---|
| v1: external CCU + single terminal anchor; 1,222 days / about 3.3 years | `paper/Limbus_Retention_Paper_draft.md` | §1, §4.1 |
| v1: episode-basis D180=3.5%, combined band 3.0–6.6% | `paper/Limbus_Retention_Paper_draft.md` | D180 combined range |
| v1: threefold validation (holdout, synthetic data, independent back-calculation) | `paper/Limbus_Retention_Paper_draft.md` | §5.1 |
| Dates of the 5 broadcasts; repeated disclosure of DAU/MAU charts | Official broadcast schedule (public) | entire text |
| v1 fitting uses only the single terminal value | `paper/Limbus_Retention_Paper_draft.md` | §4.1 |
| Pre-extraction freeze of the 3 tasks (a)(b)(c) | `v2/Digitization_Protocol_frozen.md` | opening freeze declaration, §0 |

### 2.1 Digitization procedure

Provenance of the numbers appearing in this section (number → file → location). 24 items in total.

| Number | File | Location |
|---|---|---|
| Dates of the 5 broadcasts (2024-11-22 etc.); 10 images in total | Official broadcast schedule (public) | entire text |
| 3rd Anniversary broadcast 2026-02-27; DAU 2023-08–2026-01 / MAU 2024-07–2026-01 | `v2/Digitization_Protocol_frozen.md` | §1 Material |
| DAU x-axis: uniform 126 px, interval days 335/184/181/184 | `v2/digitized/26-01-DAU_calibration.md` | §1 table |
| DAU first interval about 0.38 px/day; prohibited for (b) | `v2/digitized/26-01-DAU_usage_conditions.md` | usage condition 4 |
| DAU first interval compression ratio about 1.8×; error up to several weeks | `v2/Digitization_Protocol_frozen.md` | change log 2026-07-19 |
| MAU x-axis affine holds; residual at most 1.12 days | `v2/digitized/26-01-MAU_calibration.md` | §0 x-axis |
| MAU +6-month axis-shift correction (label 25-07 = actual 26-01) | Official broadcast schedule (public) | broadcast schedule |
| MAU internal verification: 18-month overlap, 0 violations across all 4 series | `v2/digitized/26-01-MAU_usage_conditions.md` | usage condition 8 |
| DAU run-to-run RMSE STEAM 8.20 / AOS 5.33 / TOTAL 3.01% | `v2/digitized/26-01-DAU_calibration.md` | §5 [R4] table |
| DAU STEAM 8.20% recent-biased, treated as lower bound | `v2/digitized/26-01-DAU_usage_conditions.md` | usage condition 1 |
| MAU run-to-run RMSE STEAM 0.92 / AOS 0.66 / TOTAL 0.29 / iOS 2.96% | `v2/digitized/26-01-MAU_calibration.md` | §3 table |
| Run-to-run RMSE is a lower bound on total reading error (common mode) | `v2/digitized/26-01-DAU_calibration.md` | §5 [R4] note |
| DAU terminal contamination zone (from 2026-01-27) set to missing | `v2/Digitization_Protocol_frozen.md` | change log 2026-07-19 |
| DAU anchor x=785=2026-01-25: STEAM −3.14 / AOS −2.88 / TOTAL −9.71% (run1) | `v2/digitized/26-01-DAU_calibration.md` | §4 table |
| TOTAL non-additive; label ratio 0.933; ±10% | `v2/digitized/26-01-DAU_usage_conditions.md` | usage condition 3 |
| DAU iOS color separation impossible, RMSE 15.72–20.7%, excluded from quantitative use | `v2/digitized/26-01-DAU_calibration.md` | §5, §6 |
| MAU iOS RMSE 2.96%, quantitatively adopted (asymmetric) | `v2/digitized/26-01-MAU_usage_conditions.md` | usage condition 4 |
| MAU anchor error AOS +4.23 / iOS +7.35% (run1) | `v2/digitized/26-01-MAU_calibration.md` | §5 table |
| MAU AOS ±5% / iOS ±8% systematic uncertainty | `v2/digitized/26-01-MAU_usage_conditions.md` | usage condition 2 |
| MAU AOS>Steam crossing interval, visual confirmation = real data | `v2/digitized/26-01-MAU_calibration.md` | §8 |
| Independent re-derivation of run2 reference colors (k-means) | `v2/digitized/26-01-DAU_calibration.md` | §3 [R4] |
| Remaining 8 images' anchors, all 32 labels agree, 0 discrepancies | `v2/task_c/anchor_verification.md` | cross-check result |
| DAU audit 2 rounds, 5 returns | `v2/digitized/26-01-DAU_usage_conditions.md` | opening |
| MAU 3 returns (iOS integration / AOS halo / AOS>Steam) | `v2/digitized/26-01-MAU_calibration.md` | §7, §8, §9 |

### 2.2 Verification design

| Number / fact | File | Location |
|---|---|---|
| Definition of the 3 tasks (a)(b)(c); pre-extraction freeze | `v2/Digitization_Protocol_frozen.md` | opening freeze declaration, §0 |
| (a) descriptive report of relative RMSE; reference point −5.4% | `v2/Digitization_Protocol_frozen.md` | §3(a) |
| (b) both h=2.2 fixed / free reports; acceptance of non-identifiability | `v2/Digitization_Protocol_frozen.md` | §3(b) |
| (c) chronological fixing of addition order; pre-enumeration of 4 reported quantities; application of relaxation diagnostic | `v2/Digitization_Protocol_frozen.md` | §3(c) |
| Pre-disclosure of known risks (~0.8 px/day etc.) | `v2/Digitization_Protocol_frozen.md` | §4 |
| 2 categories of deviation (post-hoc method selection / forced by data constraints) | `v2/Digitization_Protocol_frozen.md` | opening freeze declaration |
| Forced piecewise-linear x-axis (calibration-stage discovery); first interval prohibited | `v2/Digitization_Protocol_frozen.md` | change log 2026-07-19 |
| Task c relaxation-diagnostic omission (delegation-instruction error) and correction | `v2/Digitization_Protocol_frozen.md` | change log 2026-07-20 |
| Making the audit mandatory: execution → audit → conditional adoption | `v2/Digitization_Protocol_frozen.md` | §3 audit |
| Audit detection of task c deviation (TC-1); correction before finalizing adoption | `v2/task_c/task_c_usage_conditions.md` | proper record as an experiment |

### 3.1 Trajectory verification (task a)

| Number | File | Location |
|---|---|---|
| Baseline relative RMSE 0.130457; environment-difference band 0.1305–0.1314; difference 0.00003 | `v2/task_a/task_a_summary.md` | §0 table |
| c=1.0, τ=6000 upper-limit pinning (v1 known behavior); h=2.128; β=0 | `v2/task_a/task_a_summary.md` | §0 table, self-check |
| Per-interval table (28.34/42.40/24.14/22.57/19.59%, −16.68/−39.05/−16.00/−7.77/−5.41%); n=509 | `v2/task_a/task_a_usage_conditions.md` | finalized figures table |
| Interval 1 low confidence (date error up to several weeks); conclusion holds for intervals 2–4 | `v2/task_a/task_a_usage_conditions.md` | recommendation 4 |
| run2 intervals 1→4: −34.263% → −0.126% (reference) | `v2/task_a/task_a_summary.md` | §3 table |
| TOTAL shape relative RMSE 24.80% | `v2/task_a/task_a_summary.md` | §4 |
| Reading-error floor 8.20%; full-period RMSE about 3.5×; only interval-4 bias below floor | `v2/task_a/task_a_summary.md` | §5 table |
| Extraction bias (label ratio −3.14%) → underestimation magnitude is a lower bound | `v2/task_a/task_a_usage_conditions.md` | recommendation 3 |
| 2026-01-21 model/official=1.101; interval 4 is a composite of early underestimation and terminal excess | `v2/task_a/task_a_usage_conditions.md` | recommendation 3 |
| model_dau = fitted CCU × 24/h; descriptive frame of the constant-h confound | `v2/task_a/task_a_usage_conditions.md` | [TA-1] |
| Definition difference from v1 holdout −5.4% (in/out-of-sample, digitization/label basis) | `v2/task_a/task_a_usage_conditions.md` | [TA-2] |
| Pre-declaration of descriptive evaluation (no pass/fail threshold) | `v2/Digitization_Protocol_frozen.md` | §3(a) |

### 3.2 Empirical identification of β (task b)

| Number | File | Location |
|---|---|---|
| r(t) definition; frozen reuse of g(t); regression with intercept; first interval excluded; n=387 | `v2/task_b/task_b_summary.md` | §0 |
| v1's β=0.3 cap assumption; D180 band [3.0%, 6.6%] | `v2/task_b/task_b_summary.md` | §5 |
| Resolution median 1.0 day/point; all 31 events with 2+ points in window | `v2/task_b/task_b_summary.md` | §2 |
| β=0.443 [0.327,0.555]; β=0.541 [0.391,0.682] (main report) | `v2/task_b/task_b_usage_conditions.md` | published-values table |
| lag-1 autocorrelation 0.71; event-window block K=28; fixed 14-day K=41; all 4 exclude 0 | `v2/task_b/task_b_summary.md` | §2.5 |
| HAC L=14–60 se flat; CI ≈ [0.33,0.55]; iid CI not published | `v2/task_b/task_b_usage_conditions.md` | published values, note |
| 8.20% reading error added independently at all points (conservative choice) | `v2/task_b/task_b_summary.md` | §4 |
| Per-interval β 0.32/0.55/0.50; fixed effects 0.458 (not a confound) | `v2/task_b/task_b_summary.md` | §3 robustness |
| r median trend intervals 1→4: −0.41 → −0.09 | `v2/task_b/task_b_usage_conditions.md` | items complied with |
| h_hat=1.800 is an aggregate value, not emphasized | `v2/task_b/task_b_usage_conditions.md` | items complied with |
| [TB-2] 7.40% inside / 6.40% outside windows; resolution asymmetry → possibility of upward bias | `v2/task_b/task_b_usage_conditions.md` | [TB-2] |
| [TB-3] ±1.45 day/px errors-in-variables → downward attenuation; net indeterminate | `v2/task_b/task_b_usage_conditions.md` | [TB-3] |
| The qualified expression "suggests above 0.3 but cannot be asserted definitively" | `v2/task_b/task_b_usage_conditions.md` | [TB-2] |
| Formal fixed-β refit: upward-revision direction robust; D180 band ≈7.6–9.5% with λ non-identified (no point estimate) | `v2/task_b/task_b_beta_refit.csv` | §3.2.3 |
| Suggests the direction of an upward revision (not compression) | `v2/task_b/task_b_summary.md` | §5 conclusion |

### 3.3 Anchor marginal-value experiment (task c)

| Number | File | Location |
|---|---|---|
| Chronological fixing of addition order; pre-enumeration of 4 reported quantities | `v2/Digitization_Protocol_frozen.md` | §3(c) |
| 8 charts, 32 labels dual reading, 0 discrepancies | `v2/task_c/anchor_verification.md` | cross-check result |
| Unified calendar-month mask convention (replication of 26-01 convention, avoiding mixed treatment) | `v2/task_c/task_c_summary.md` | §2 |
| MultiAnchorFitter extension (BOUNDS, kernel, optimization inherited) | `v2/task_c/task_c_summary.md` | §3 |
| 5-stage table (D180, boundary pinning, CI width) | `v2/task_c/task_c_usage_conditions.md` | finalized figures table |
| CI width qualitative-only (n_boot=100, stage 5 independent seed) | `v2/task_c/task_c_usage_conditions.md` | [TC-5] |
| τ flip (upper 6000 → lower 200); k@1.5 new pinning | `v2/task_c/task_c_summary.md` | §4 breakdown |
| Multi-start: 1–3 corner convergence (apparent stability); 4–5 obj relative difference 3.5%/2.3% | `v2/task_c/task_c_summary.md` | §5 |
| base degenerates to 0 in all stages | `v2/task_c/task_c_summary.md` | §7-4 |
| Distrust of structural-decomposition values (3-parameter simultaneous pinning) | `v2/task_c/task_c_usage_conditions.md` | [TC-1/d] |
| Relaxation diagnostic: stage 3 unchanged; stages 4–5 τ interior 131–134 days; k@3.0 re-pinning; D180 ≈ 6.1% | `v2/task_c/task_c_usage_conditions.md` | diagnostic (TC-1 relaxation) |
| Collaboration window 28 days = no effect (D180 difference 0.02–0.03 pt) | `v2/task_c/task_c_summary.md` | §8b [TC-2] |
| 25-10 terminal label 2025-11-09; anchor 418,242 is the Nov-09 value (timepoint error) | `v2/task_c/task_c_usage_conditions.md` | [TC-2/c] |
| November mask: stage 4 fully resolved (D180=4.01%); stage 5 partially resolved (D180=4.53%) | `v2/task_c/task_c_usage_conditions.md` | diagnostic (TC-2 sensitivity) |
| Report both 6.72%/4.53%; single citation prohibited | `v2/task_c/task_c_usage_conditions.md` | [TC-2/c] |
| "November mask × bound relaxation" not run; stage-5 residual short τ cannot be asserted | `v2/task_c/task_c_usage_conditions.md` | usage condition 6 |
| Regime multiplicity; 4.0–6.7% swing; "convergence" description prohibited | `v2/task_c/task_c_usage_conditions.md` | [TC-3] |
| Qualified conclusion on marginal value (unachieved, approximately nil to negative, false stability) | `v2/task_c/task_c_usage_conditions.md` | [TC-4] |

### 4. Theory

Provenance of the numbers and definitions appearing in this section (number → file → location).

| Number / definition | File | Location |
|---|---|---|
| R(d) decomposition; definition of base/stock columns (implementation extraction) | `v2/task_d/collinearity_note.md` | §1 |
| base=K−G′·x_stock+O(ε); reduction to the space {1, e^{−t/τ}}; invertible triangular map | `v2/task_d/collinearity_note.md` | §2 |
| Collapse mode (i) CV≈(T−t₀)/(τ√12); example CV≈5.4% (T−t₀≈1132, τ=6000) | `v2/task_d/collinearity_note.md` | §3(i) |
| Collapse mode (ii) τ≲t₀/ln(1/δ); zero-column (u≈0.05 at τ=30) | `v2/task_d/collinearity_note.md` | §3(ii) |
| NNLS corner forcing; coin toss of noise; reproduction stability is not evidence of identification | `v2/task_d/collinearity_note.md` | §4 |
| Necessary but not sufficient condition (unidentified even with 5 anchors); t₀/T design determines the identification band | `v2/task_d/collinearity_note.md` | §5 |
| Limitations (breakdown of F saturation, higher dimension from event-column tails, weighted geometry, non-normative CV threshold) | `v2/task_d/collinearity_note.md` | §6 |
| Confirmation of implementation consistency of the equation material | `v2/task_d/reconciliation_report.md` | §4 |
| Apparent stability from corner convergence in stages 1–3; τ upper 6000; base=0 degeneracy | `v2/task_c/task_c_summary.md` | §5, §7-4 |

### 5. Discussion

Provenance of the numbers and facts appearing in this section (number → file → location).

| Number / fact | File | Location |
|---|---|---|
| Intervals 2→4 bias −16.0% → −5.4%; shrinks near the anchor | `v2/task_a/task_a_usage_conditions.md` | finalized figures table |
| β=0.443 [0.327, 0.555]; CI clearly excludes 0 | `v2/task_b/task_b_usage_conditions.md` | published-values table |
| D180 4.0–6.7% swing; report both 6.72%/4.53%; marginal value approximately nil | `v2/task_c/task_c_usage_conditions.md` | [TC-3][TC-4][TC-2/c] |
| Level-anchor constraint has almost no component in the degenerate direction (connection from the closed form) | `v2/task_d/collinearity_note.md` | §4, §5/§4.4 |
| v1: D180=3.5%; combined band 3.0–6.6% | `paper/Limbus_Retention_Paper_draft.md` | D180 combined range |
| Definition difference between interval-4 bias −5.41% and v1 holdout −5.4% | `v2/task_a/task_a_usage_conditions.md` | [TA-2] |
| β "suggests above 0.3 but cannot be asserted definitively"; CI lower bound 0.327–0.391 | `v2/task_b/task_b_usage_conditions.md` | [TB-2] |
| Upward-revision direction to the v1 band robust under fixed-β refit; level a band (λ non-identified) | `v2/task_b/task_b_beta_refit.csv` | §5.2 |
| Zero external institutional-investor stake; all preferred shares converted to treasury shares (end of 10th fiscal year) | FSS DART rcpNo 20260401000553 (https://dart.fss.or.kr/dsaf001/main.do?rcpNo=20260401000553; extract: `v2/paper/sources/dart_2026_notes.md`) | Shareholder composition (Note 1, p.13) |
| Advertising expense / operating revenue ≈ 0.062% (2025); 0.098% (2024) | FSS DART rcpNo 20260401000553 (https://dart.fss.or.kr/dsaf001/main.do?rcpNo=20260401000553; extract: `v2/paper/sources/dart_2026_notes.md`) | Income statement (p.8) |
| Statement of qualitative, non-causal nature / prohibited-inference items | `v2/paper/sources/dart_2026_notes.md` | Usage / prohibitions |
| Person-basis D180 ≈27% (20–39%) is a snapshot ratio (numerator ≈ constant across 3 cases); band propagates denominator only; monthly kernel is a boundary solution | `paper/Limbus_Retention_Paper_draft.md` | §5.4 |
| Total activations ≈10.76M; uniques ≈2.78M; ≈3.9 assumes β=0 | `paper/Limbus_Retention_Paper_draft.md` | §5.4 |
| β=0.443 > assumed cap 0.3 → total activations overestimated → episodes per person on the upper side (inflow reduced by β) | `v2/task_b/task_b_usage_conditions.md`+`paper/Limbus_Retention_Paper_draft.md` | published-values table / §5.2, §8(1) |

### 6. Limitations and future work

Provenance of the numbers and facts appearing in this section (number → file → location).

| Number / fact | File | Location |
|---|---|---|
| 2 of 10 images digitized; remaining 8 anchors only; 32 labels agree | This paper (§2.1) | §2.1.1, §2.1.6 |
| β unit correspondence holds by construction (same form/g); fixed-β refit performed → band + λ non-identified | `v2/task_b/task_b_beta_refit.csv` | §3.2.3/§6 |
| n_boot=100; CI width for qualitative use | `v2/task_c/task_c_usage_conditions.md` | [TC-5] |
| "November mask × bound relaxation" not run; stage-5 residual short τ cannot be asserted | `v2/task_c/task_c_usage_conditions.md` | usage condition 6 |
| Original residual 26% (0.000662) unexplained; weekday-adjustment script missing; not identifiable | `v2/task_d/reconciliation_report.md` | §0, §3, §6 |
| No title with an official anchor in the multi-title verification | `paper/Limbus_Retention_Paper_draft.md` | §6.5 |
| Three-layer person-basis figures out of v2 diagnostic scope (0 mentions in v2); 27% is a snapshot ratio; 3.9 on the β upper side | `paper/Limbus_Retention_Paper_draft.md` | §5.4 (person-basis) |

### Abstract / 7. Conclusion

Provenance of the numbers and facts appearing in this section (number → file → location).

| Number / fact | File | Location |
|---|---|---|
| Repeated disclosure across 5 broadcasts; single-title deep dive | Official broadcast schedule (public) | entire text |
| β=0.443 [0.327, 0.555]; CI clearly excludes 0 | `v2/task_b/task_b_usage_conditions.md` | published-values table |
| Closed form of the stock/base collinearity | `v2/task_d/collinearity_note.md` | §2–4 |
| Report both D180 6.72%/4.53%; not cited as a definitive estimate | `v2/task_c/task_c_usage_conditions.md` | [TC-2/c][TC-4] |
| Near-anchor constraint; systematic underestimation over the full period | `v2/task_a/task_a_usage_conditions.md` | finalized figures table, [TA-2] |
| Pre-extraction freeze; separation of execution and audit | `v2/Digitization_Protocol_frozen.md` | opening freeze declaration, §3 audit |
