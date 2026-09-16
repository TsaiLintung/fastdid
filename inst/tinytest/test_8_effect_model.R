tol <- 1e-2 #allow 1% different between estimates
base <- list(timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y")
run <- function(data, ...) do.call(fastdid, c(list(data = data), base, list(...)))

#' A panel with unit and time fixed effects, known event effects, and small noise.
#' Event 1 follows `eff1` in event time. Event 2 grows by 0.5 a period. Event 3
#' grows by 3 a period. Each element of `cohorts` is the triple (g1, g2, g3).
make_double_dt <- function(cohorts, TT = 6, n = 500, seed = 7, eff1 = function(e) 1) {
  set.seed(seed)
  dt <- data.table::rbindlist(lapply(seq_along(cohorts), function(i) {
    g <- cohorts[[i]]
    data.table::CJ(unit = (i - 1) * n + seq_len(n), time = seq_len(TT))[
      , `:=`(G = g[1], G2 = g[2], G3 = g[3])]
  }))
  dt[, y := unit / 100 + time * 0.3 + stats::rnorm(.N, 0, 0.01) +
       data.table::fifelse(time >= G, eff1(time - G), 0) +
       data.table::fifelse(time >= G2, 0.5 * (time - G2 + 1), 0) +
       data.table::fifelse(time >= G3, 3 * (time - G3 + 1), 0)]
  dt[]
}

# backward compatibility ------------------------------------------------------

simdt <- sim_did(1e+03, 5, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1,
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt
inv <- unique(dt[, .(unit, G, G2)])
keep_units <- inv[, if(.N >= 2) .SD, by = .(G, G2)][, unit]
dt <- dt[unit %in% keep_units]

for (rt in c("group_time", "group_group_time", "dynamic_stagger", "dynamic")) {
  expect_equal(run(dt, result_type = rt, cohortvar2 = "G2"),
               run(dt, result_type = rt, cohortvar2 = "G2", effect_model = "parallel"),
               info = paste("the parallel preset is the default estimator,", rt))
}

# the parallel formula ---------------------------------------------------------

# with one clean cell per modeled cohort the WLS fit is the two-period comparison
parallel_formula <- ~ factor(gvec) + factor(gactive):factor(t)

dt_imp <- make_double_dt(list(c(2, 3, Inf), c(2, Inf, Inf), c(Inf, Inf, Inf)))
res_par <- run(dt_imp, result_type = "group_group_time", cohortvar2 = "G2")
res_frm <- run(dt_imp, result_type = "group_group_time", cohortvar2 = "G2", effect_model = parallel_formula)
expect_equal(res_frm[, att], res_par[, att], tolerance = 1e-8,
             info = "parallel formula equals the parallel preset, imputation case")
expect_equal(res_frm[, se], res_par[, se], tolerance = 0.2,
             info = "parallel formula se is close to the preset, imputation case")

dt_did <- make_double_dt(list(c(3, 2, Inf), c(Inf, 2, Inf), c(Inf, Inf, Inf)))
res_par <- run(dt_did, result_type = "group_group_time", cohortvar2 = "G2")
res_frm <- run(dt_did, result_type = "group_group_time", cohortvar2 = "G2", effect_model = parallel_formula)
expect_equal(res_frm[, att], res_par[, att], tolerance = 1e-8,
             info = "parallel formula equals the parallel preset, did case")

# with several clean cells the fit pools them, so the two agree within noise
res_par <- run(dt, result_type = "group_group_time", cohortvar2 = "G2")
res_frm <- run(dt, result_type = "group_group_time", cohortvar2 = "G2", effect_model = parallel_formula)
both <- merge(res_par, res_frm, by = c("cohort1", "cohort2", "time"))
expect_true(nrow(both) == nrow(res_par), info = "parallel formula reports the same cells as the preset")
expect_true(both[, max(abs(att.x - att.y) / pmax(se.x, se.y))] < 4, info = "parallel formula is within noise of the preset on simulated data")

# polynomial known truth -------------------------------------------------------

dt_poly <- make_double_dt(list(c(2, 5, Inf), c(2, Inf, Inf), c(3, Inf, Inf), c(Inf, Inf, Inf)), TT = 8,
                          eff1 = function(e) 1 + 0.5 * e + 0.1 * e^2)
res_poly <- run(dt_poly, result_type = "group_group_time", cohortvar2 = "G2",
                effect_model = ~ factor(gvec) + poly(e, 2), full = TRUE)
est <- res_poly$estimate[cohort1 == 2 & cohort2 == 5]
truth <- 1 + 0.5 * (est[, time] - 2) + 0.1 * (est[, time] - 2)^2
expect_equal(est[, att], truth, tolerance = tol, info = "quadratic effect model recovers the truth")
expect_equal(nrow(est), 7L, info = "quadratic effect model imputes every confounded period")

res_unr <- run(dt_poly, result_type = "group_group_time", cohortvar2 = "G2", effect_model = "unrestricted")
expect_equal(res_unr[cohort1 == 2 & cohort2 == 5, time], c(2, 3, 4),
             info = "unrestricted reports the clean cells only")

# diagnostics
diag <- res_poly$effect_diag$y
expect_equal(dim(diag$weights), c(nrow(res_poly$estimate), nrow(res_poly$gt_estimate$y$gt)),
             info = "diagnostics: one weight row per reported cell, one column per first-stage cell")
expect_true(all(c("G", "time", "event", "component", "estimable", "leverage", "e_gap") %in% names(diag$cells)),
            info = "diagnostics: cell table columns")
expect_equal(diag$cells[component == "target", e_gap], 1:4, info = "diagnostics: extrapolation gap in event time")
expect_equal(length(diag$fit$target$residuals$resid), diag$fit$target$n_fit, info = "diagnostics: one residual per fit cell")
expect_equal(diag$overid$target$df, diag$fit$target$n_fit - diag$fit$target$rank, info = "diagnostics: over-identification df")
expect_true(diag$overid$target$pvalue > 0.01, info = "diagnostics: the correct model is not rejected")

# rank deficiency and validation ----------------------------------------------

# a saturated model imputes nothing, so only the direct cells remain
expect_warning(res_sat <- run(dt_poly, result_type = "group_group_time", cohortvar2 = "G2",
                              effect_model = ~ 0 + factor(gvec):factor(t)),
               pattern = "not estimable", info = "a saturated model identifies nothing")
expect_equal(res_sat[cohort1 == 2 & cohort2 == 5, time], c(2, 3, 4), info = "a saturated model keeps the direct cells")
# a quadratic needs three clean cells of its own cohort
dt_two <- make_double_dt(list(c(2, 4, Inf), c(2, Inf, Inf), c(Inf, Inf, Inf)))
expect_warning(run(dt_two, result_type = "group_group_time", cohortvar2 = "G2",
                   effect_model = ~ factor(gvec) + factor(gvec):poly(e, 2)),
               pattern = "not estimable", info = "a cohort-specific polynomial with too few clean cells warns")
expect_error(run(dt_poly, result_type = "group_group_time", cohortvar2 = "G2", effect_model = ~ foo),
             pattern = "available variables", info = "unknown formula variable")
expect_error(run(dt_poly, result_type = "group_group_time", cohortvar2 = "G2", effect_model = y ~ e),
             info = "two-sided formula")
expect_error(run(dt_poly, result_type = "group_group_time", cohortvar2 = "G2", effect_model = "bad"),
             info = "unknown preset")
expect_error(run(dt_poly, result_type = "group_group_time", cohortvar2 = "G2", effect_fit = "joint"),
             pattern = "needs a formula", info = "a stacked fit needs a formula")
expect_error(run(dt_poly, result_type = "group_time", effect_model = ~ e),
             pattern = "cohortvar2", info = "a formula needs multiple events")
expect_error(run(dt_poly, result_type = "group_time", cohortvar2 = "G2", event_specific = FALSE, effect_model = ~ e),
             pattern = "event_specific", info = "a formula needs the event-specific effect")
expect_error(run(dt_poly, result_type = "dynamic_event", cohortvar2 = "G2", effect_model = ~ e),
             pattern = "ordered", info = "dynamic_event needs the ordered fit")
expect_error(run(dt_poly, result_type = "group_time", cohortvar2 = "G2", effect_model = ~ poly(e2, 2)),
             pattern = "effect model", info = "an infinite event time in the formula")

# ordered events ---------------------------------------------------------------

make_ordered <- function(cohorts, TT = 8, n = 400, seed = 3, f = function(e) 1 + 0.5 * e) {
  set.seed(seed)
  dt <- data.table::rbindlist(lapply(seq_along(cohorts), function(i) {
    g <- cohorts[[i]]
    data.table::CJ(unit = (i - 1) * n + seq_len(n), time = seq_len(TT))[, `:=`(G = g[1], G2 = g[2])]
  }))
  dt[, y := unit / 100 + time * 0.3 + stats::rnorm(.N, 0, 0.01) +
       data.table::fifelse(time >= G, f(time - G), 0) + data.table::fifelse(time >= G2, f(time - G2), 0)]
  dt[]
}
dt_ord <- make_ordered(list(c(2, 4, Inf), c(2, 5, Inf), c(3, Inf, Inf), c(Inf, Inf, Inf)))
res_ord <- run(dt_ord, result_type = "dynamic_event", cohortvar2 = "G2",
               effect_model = ~ 0 + factor(e), effect_fit = "ordered")
post <- res_ord[event_time >= 0]
expect_equal(post[, att], 1 + 0.5 * post[, event_time], tolerance = tol,
             info = "ordered fit recovers the shared event-time profile")
expect_silent(plot_did_dynamics(res_ord), info = "dynamic_event plots")

res_ord2 <- run(dt_ord, result_type = "group_group_time", cohortvar2 = "G2",
                effect_model = ~ 0 + factor(gvec):factor(event) + factor(e), effect_fit = "ordered")
est <- res_ord2[cohort1 == 2 & cohort2 == 4]
expect_equal(est[, att], 1 + 0.5 * (est[, time] - 2), tolerance = tol,
             info = "ordered fit with cohort levels recovers the first-event effect")

dt_bad <- data.table::copy(dt_ord)
dt_bad[G == 3, G2 := 2]
expect_error(run(dt_bad, result_type = "dynamic_event", cohortvar2 = "G2",
                 effect_model = ~ 0 + factor(e), effect_fit = "ordered"),
             pattern = "g1 < g2", info = "ordered fit rejects an unordered cohort")
expect_error(run(dt_ord, result_type = "dynamic_event", cohortvar2 = "G2", anticipation2 = 1,
                 effect_model = ~ 0 + factor(e), effect_fit = "ordered"),
             pattern = "anticipation", info = "ordered fit needs one anticipation horizon")

# joint fit --------------------------------------------------------------------

res_joint <- run(dt_imp, result_type = "group_group_time", cohortvar2 = "G2",
                 effect_model = ~ 0 + factor(gvec):factor(event) + factor(gown):factor(t):factor(event),
                 effect_fit = "joint")
est <- res_joint[cohort1 == 2 & cohort2 == 3 & time >= 3]
expect_equal(est[, att], rep(1, nrow(est)), tolerance = tol, info = "joint fit recovers the first-event effect")
