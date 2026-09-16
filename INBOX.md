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


