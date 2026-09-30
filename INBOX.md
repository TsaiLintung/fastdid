# INBOX

Drop-box for this repo: humans and agents working in other repos leave todo items or information here. Entries are processed only on explicit instruction — never act on this file unprompted. Newest entries first, each with a date and source, self-contained enough to act on without the conversation that produced it. Delete an entry once it is resolved (move anything durable into the repo's documentation first).

---

## 2026-09-15 — first-stage extensions that the effect-model second stage can inherit

Source: Claude Code session with the author, from the scope paragraph at the
end of Section 5 of `project/doubledid/latex/paper2/` (commit 961890f). The
second stage (`R/second_stage.R`, `R/effect_model.R`) reads only the cell
table: `att`, `inf_func`, `pg`, and the cohort label. Any first stage that
returns those, with each cell a mean difference between a dated history and
a baseline path, plugs in. The items below are ordered by the work they need.
None is started.

1. **Repeated cross-sections.** The CS estimators for repeated cross-sections
   (Callaway and Sant'Anna 2021, Section 4.2) as a first-stage path next to
   `estimate_did_rc()`. The cohort vector must be a group-level variable. The
   second stage is unchanged. Test: equality with the `did` package on a
   single event, then a known-truth double DiD on a repeated cross-section.
2. **Bounds on parallel-trends violations.** No new first stage. Export the
   second-stage weight matrix (`effect_diag$weights`) and the first-stage
   influence functions in the shape that `HonestDiD` takes: the event-study
   vector, its covariance, and the linear functional `l`. Rambachan and Roth
   (2023) then bound any second-stage estimate. Test: a single-event case
   reproduces `HonestDiD` on the CS event study.
3. **Nonlinear counterfactual mean** (Wooldridge 2023, Econometrics Journal
   26(3)). A first stage that fits the untreated mean on a transformed scale
   (Poisson, logit) and returns the cell as a level-scale mean difference with
   its influence function. The second stage is unchanged. Test: a known-truth
   count outcome.
4. **Triple differences** (Olden and Moen 2022, Econometrics Journal 25(3)).
   A first stage whose cell is the triple difference across a stratum
   variable, with its influence function. The second stage is unchanged.
   Test: a known-truth DGP with a stratum-specific trend.
5. **Treatment on in the first period.** The status-quo path as the baseline
   of every cell (de Chaisemartin and D'Haultfoeuille 2026, REStat 108(4)):
   condition the first stage on the period-one status, take the events as the
   switches after period one, and stop dropping always-treated units in
   `validate_dt()`. The least tested of the five. Test: a known-truth DGP with
   half the units on at period one.

Already covered and not on this list: covariates (the dr first stage), a known
anticipation horizon (`anticipation`, `anticipation2`), unbalanced panels
(`allow_unbalance_panel`). Not doable as a first-stage switch: no
not-yet-treated units (needs a joint fit of the baseline and the effect
model), continuous treatment, distributional parameters, fuzzy designs.



---

## 2026-09-17 — the over-identification test over-rejects under a true model with many cells

Source: Claude Code session with the author, from the ordered-events Monte
Carlo of `project/doubledid/scripts/5c-ordered-sim.R` (fastdid working tree,
branch `effect-model`, 1.1.0). `get_effect_overid()` in `R/effect_model.R`
computes the residual quadratic form in the pseudo-inverse of the
influence-function variance, with the degrees of freedom of the null space.
With N = 2,000 units, T = 10, and a stacked fit on about 150 cells (about
125 restrictions), the test rejects a true model far above the nominal 5%:
OVERID-NUMBERS. It has power (the no-carryover model under a carryover is
rejected in every draw), so the statistic is not broken, but the size is.
The likely cause is the many-restrictions problem: the n-by-n variance of
the fit residuals is estimated from N units and inverted. Candidates: a
finite-sample correction, a test on a low-dimensional set of contrasts (for
example the mean residual per state), or a multiplier bootstrap of the
statistic. Until then, document the test as a diagnostic, not a size-correct
test, in `R/fastdid.R` and in paper2's Proposition `prop:overid`. The
per-draw p-values are in `project/doubledid/est/osim_draws.csv`, column
`overid_p`.

## 2026-09-17 — the share influence of a stacked fit uses the fitted cells

Source: the same session. Under a stacked fit (`effect_fit = "state"` or
`"ordered"`), the aggregate of the reported contrasts over the fit cells with
the `pg` weights equals the aggregate of the raw cells whenever the columns
of the model span the indicator of those cells (the fit is a weighted
projection). The point estimate is then unbiased even when the cell-level
model is wrong. The influence of the cohort shares in the aggregation,
however, is computed from the fitted cell estimates, which vary less across
cohorts than the raw cells under a cohort-level misspecification, so the
standard error understates in that case (seen with a carryover level by
switch-off cohort under the state model: Monte Carlo s.d. about 0.32 against
an average s.e. of 0.19 at 10 replications). Decide whether the share
influence must use the raw cells, or whether the documentation states that
the standard error assumes the model. The design that showed it is not in
the paper; see `agent/ISSUES.md` at the mono root.
