# ISSUES

Decisions with their reasons, and traps to avoid.

### Key design decisions

1. **Sort once.** Coercion sorts the data by `(time, G, unit)`. Every later
   access is sequential and fast. Do not reorder the data after coercion.
2. **Two estimation paths.** `estimate_did_bp()` handles a balanced panel.
   `estimate_did_rc()` handles an unbalanced one.
3. **Influence functions.** The package computes a full semiparametric influence
   function at the 2x2 level, then aggregates it. The bootstrap runs through
   `BMisc::multiplier_bootstrap`.
4. **Locked lists.** A parameter object carries the `locked` S3 class. `$<-`,
   `[<-` and `[[<-` all raise an error, so nothing mutates it.
5. **Time normalization.** User time (for example 2000, 2005, 2010) becomes
   1, 2, 3 internally. `time_offset` and `time_step` recover the original scale.
6. **Two stages.** The first stage is Callaway-Sant'Anna on the cohort label.
   With `cohortvar2` the label is the G-string `"g1-g2-...-gM"`, and the first
   stage identifies the combined effect of all events in each cell. The second
   stage (`second_stage.R`) reads the cell table only. It returns one weight
   row per identified cell over the first-stage cells, so the estimate and its
   influence function are linear combinations of the first stage. The
   parallel-treatment-effects scheme uses three cases, chosen by the relative
   timing of g1 and g', and implements Theorem 3 of Tsai (2026).
7. **IPW caching.** When several outcomes share the same covariates, the
   propensity scores are computed once and reused.
8. **Effect models.** A formula in `effect_model` states a linear model for one
   component of a cell. `effect_model.R` fits it by weighted least squares on
   the clean cells, checks that each target row lies in the row space of the
   fit design, and writes the imputed value as a weight row over first-stage
   cells. The WLS weights are the cohort shares and are treated as fixed in
   the influence function. The `"parallel"` preset dispatches to the original
   scheme, so it is bit-identical to earlier versions.
9. **Switches and levels.** The ordered and the state fit model the k-th
   switch of one treatment. `add_switch_columns()` adds the date `S_k` and the
   level `L_k` of the k-th switch to the cell table. Without
   `exper$event_level` the cohort columns are in the order of the switches and
   the level is the parity. With it, each cohort column is one event "switch
   to this level", and the k-th switch is the k-th date of the cohort. The
   columns `G1..GM` stay in the order of the cohort label, because
   `get_weight_influence()` sorts the rows by them. Do not overwrite them with
   the sorted dates. The other result types read the target event from the
   first cohort column, so `event_level` needs `dynamic_event` or an
   `aggregate_scheme`.
10. **The strata path.** `control_type = "strata"` replaces the loop of 2x2
   estimates with `estimate_gtatt_outcome_strata()`. A stratum is one value of
   `covariatesvar`. Every cell is a function of the sums of the outcome by
   cohort row, stratum, and period, and the control cohorts of a cell are a
   range of rows, so prefix sums over the rows give every control mean. The
   influence function is a `strata_inf` object, not a matrix. Within a group
   of units with one cohort vector and stratum, it is linear in the outcome
   path of the unit, so the targets need only its coefficients
   (`strata_coef()`), never a matrix of the units times the cells. The
   terms are those of `estimate_did_bp()` with `control_type = "reg"` and
   one dummy for each stratum, with its
   small-group inflation, so the two paths agree to 1e-8
   (`inst/tinytest/test_9_strata.R`). Everything else is shared: validation,
   coercion, the base period, `get_control_range()`, the skip rules, the cell
   table, the second stage, and the standard errors. A result type with one
   target per cell (`group_time`) still builds a matrix of the targets times
   the cells, so a large panel needs a result type with few targets and
   custom weights.
11. **The first stage as an object.** `new_first_stage()` holds the cells,
   the cohort shares, and one uncertainty object. The second stage reads
   only this object, so `double_did()` runs on the cells of fastdid or of
   another estimator. The uncertainty has four classes:
   - `unit_inf`: the influence functions of the units (the 2x2 path, or
     `as_first_stage(inf_func =)`).
   - `strata_inf`: the strata object above.
   - `strata_vcov`: the factored covariance of the strata path, from
     `first_stage()`. It needs no unit data.
   - a dense covariance matrix of the cells (`as_first_stage(vcov =)`). The
     shares are fixed, and the bootstrap stops with an error.

   The map to the targets is lazy (`lin = list(W, Wm, cc)`): the weights
   on the cells, the weights on the means, and the influence of the shares,
   which is constant within a cohort. `lin_compose()` adds each step, and
   `lin_se()` evaluates the map once, for the targets.
12. **The factored covariance.** For a group k of the strata path,
   `strata_moments()` centers the outcome paths at the group mean, forms the
   moments of `(1, y)`, and keeps the factor `F_k` of their eigen
   decomposition. `F_k` has at most `min(n_k, T + 1)` rows, so the object is
   small. The covariance of a target is `sum_k |F_k theta_k|^2 / N^2`. Do
   not use the moments without the centering: at 150,000 units the standard
   error then differs by 2e-3 in relative terms.
