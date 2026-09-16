# TODO

Open questions and known limitations. Each entry names the code it affects.
Move an entry out of this file when it is resolved.

## Estimation

- **`R/estimate_did.R`, `estimate_did_bp()`, the outcome-regression branch.**
  The OR path refits the regression for every (g, t) cell. Optimize it with a
  better backend and a cache. `speedlm` is one candidate.

- **`R/estimate_did.R`, `estimate_did_rc()`, the `score_ps` term.** Confirm the
  `n / (n_pre + n_post)` factor on the propensity-score weights. A unit present
  in both the pre and the post period enters twice, so the factor is there to
  halve it. Verify that the same factor belongs in the influence function, where
  it is repeated.

- **`R/estimate_did.R`, `estimate_did_rc()`.** The function accepts a `cache`
  argument and never reads it. It still returns a cache. Decide whether the
  repeated-cross-section path can reuse a propensity score across outcomes. The
  weights depend on `inpre` and `inpost`, which depend on the current outcome's
  missingness pattern, so reuse is likely unsafe. If it is unsafe, drop the
  argument.

## Second stage

- **`R/effect_model.R`, `get_effect_scheme()`.** The weights of the WLS fit
  are the cohort shares `pg`, and the influence function treats them as fixed.
  The correction is second order under a correct model. If a Monte Carlo shows
  under-coverage, add the delta-method term through the pseudo-inverse, as
  `get_weight_influence_param()` does for the share weights of the parallel
  scheme.

- **`R/effect_model.R`, the filters of the fit set.** `double_control_option`
  and `control_option = "notyet"` restrict the control cohorts of the parallel
  scheme cell by cell. A global fit cannot apply them without per-target fits,
  so a formula ignores them with a warning. Decide whether a per-target local
  fit is worth its cost.

- **`R/estimate_gtatt.R`, `get_base_period()`.** The base period of a cohort
  uses `anticipation`, the horizon of event 1, even when a confounding event
  comes first. With one horizon per event the base period must be one before
  the earliest shifted date, `min_d (g^d - a_d) - 1`. The parallel scheme keeps
  the current rule. Decide before the first stage moves to that rule.

## Inference

- **`R/aggregate_gt.R`, `get_se()`, the analytic branch.** The standard error
  divides the sum of squared influence functions by `n`, not by `n - 1`. This
  matches the `did` package (Callaway and Sant'Anna), which fastdid is tested
  against. Decide whether to keep the match or to switch to `n - 1`. A switch
  breaks the equality tests in `inst/tinytest/test_2_compare_est.R`.

- **`R/aggregate_gt.R`, `get_se()`, the uniform band.** The critical value comes
  from `quantile(boot_tv, 1 - p$alpha)`. Confirm that `alpha` is the right knob
  for the band, and that the band and the pointwise interval must share it.

## Simulation

- **`R/sim_did.R`, `sim_did()`, the `vary_cov` branch.** Confirm that `xvar` is
  meant to confound the treatment. It is built from `pmin(G, time_period + 4)`,
  so it correlates with the cohort. The `pmin` also caps `G == Inf` for the
  never-treated units. Without the cap the covariate is infinite.

## Tests

- **`inst/tinytest/test_2_compare_est.R`.** fastdid and `did` do not agree on
  `control_type = "reg"` with a truly unbalanced panel. The comparison is
  disabled there. `inst/tinytest/test_99_coverage.R` covers the case instead,
  and that file is disabled with `if (FALSE)` because it is slow. Find the
  source of the disagreement.
