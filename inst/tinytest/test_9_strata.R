# control_type = "strata": the cells from the means by cohort, stratum, and period

simdt <- sim_did(3000, 6, cov = "int", seed = 3, second_cohort = TRUE)
dt <- simdt$dt
# the regression path needs one dummy for each stratum but the first
dums <- paste0("d", 2:5)
for (v in 2:5) dt[, (paste0("d", v)) := as.numeric(x == v)]

# the estimates of two calls, matched on every column but the estimates
match_est <- function(a, b) {
  keys <- setdiff(names(a), c("att", "se", "att_ciub", "att_cilb"))
  merge(a, b[, c(keys, "att", "se"), with = FALSE], by = keys, suffixes = c("", "_ref"))
}
run <- function(...) suppressWarnings(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", ...))

# without covariates, the strata path is the 2x2 path without covariates
for (rt in c("group_time", "dynamic")) {
  m <- match_est(run(control_type = "strata", result_type = rt), run(result_type = rt))
  expect_equal(nrow(m), nrow(run(result_type = rt)), info = paste("one estimate for each estimate,", rt))
  expect_equal(m$att, m$att_ref, tolerance = 1e-8, info = paste("no covariates, att,", rt))
  expect_equal(m$se, m$se_ref, tolerance = 1e-8, info = paste("no covariates, se,", rt))
}

# with strata, it is the regression path with one dummy for each stratum
for (rt in c("group_time", "dynamic", "simple")) {
  m <- match_est(run(control_type = "strata", covariatesvar = "x", result_type = rt),
                 run(control_type = "reg", covariatesvar = dums, result_type = rt))
  expect_equal(m$att, m$att_ref, tolerance = 1e-8, info = paste("strata, att,", rt))
  expect_equal(m$se, m$se_ref, tolerance = 1e-8, info = paste("strata, se,", rt))
}

# the same under the other control sets, an anticipation, and a varying base period
m <- match_est(run(control_type = "strata", covariatesvar = "x", control_option = "notyet", anticipation = 1),
               run(control_type = "reg", covariatesvar = dums, control_option = "notyet", anticipation = 1))
expect_equal(c(m$att, m$se), c(m$att_ref, m$se_ref), tolerance = 1e-8, info = "notyet with an anticipation")
m <- match_est(run(control_type = "strata", covariatesvar = "x", control_option = "never", base_period = "varying"),
               run(control_type = "reg", covariatesvar = dums, control_option = "never", base_period = "varying"))
expect_equal(c(m$att, m$se), c(m$att_ref, m$se_ref), tolerance = 1e-8, info = "never with a varying base period")

# two events: the first stage alone, and through the second stage
for (es in c(FALSE, TRUE)) {
  m <- match_est(run(control_type = "strata", covariatesvar = "x", cohortvar2 = "G2",
                     event_specific = es, result_type = "group_group_time"),
                 run(control_type = "reg", covariatesvar = dums, cohortvar2 = "G2",
                     event_specific = es, result_type = "group_group_time"))
  expect_true(nrow(m) > 0, info = paste("two events, some estimates, event_specific =", es))
  expect_equal(c(m$att, m$se), c(m$att_ref, m$se_ref), tolerance = 1e-8, info = paste("two events, event_specific =", es))
}

# combine_cells(): a one-cell target returns the cell, with the same standard error
res <- run(control_type = "strata", covariatesvar = "x", full = TRUE)
cells <- res$estimate[time >= cohort]
w <- cells[, .(target = paste(cohort, time), G = cohort, time, w = 1)]
est <- merge(combine_cells(res, w), cells[, .(target = paste(cohort, time), att, se_cell = se)], by = "target")
expect_equal(est$estimate, est$att, tolerance = 1e-10, info = "a one-cell target returns the cell")
expect_equal(est$se, est$se_cell, tolerance = 1e-8, info = "the se of a one-cell target is the se of the cell")

# a linear combination, and a weight on the means
w2 <- rbind(w[1:3][, target := "sum"], w[4][, `:=`(target = "sum", w = -2)])
r <- combine_cells(res, w2, return_if = TRUE)
expect_equal(r$estimates$estimate, sum(cells$att[1:3]) - 2 * cells$att[4], tolerance = 1e-10)
expect_equal(mean(r$inf_func[, 1]), 0, tolerance = 1e-10, info = "the influence function has mean zero")
wm <- data.table(target = "mean", G = Inf, time = c(1, 6), x = 2, w = c(-1, 1))
r3 <- combine_cells(res, w2[0], w_mean = wm, return_if = TRUE)
own <- dt[G == Inf & x == 2, .(y = mean(y)), by = time]
expect_equal(r3$estimates$estimate, own[time == 6, y] - own[time == 1, y], tolerance = 1e-10, info = "a weight on the means")
n2 <- dt[G == Inf & x == 2, uniqueN(unit)]
dy <- dt[G == Inf & x == 2][order(unit, time), .(d = y[time == 6] - y[time == 1]), by = unit]$d
expect_equal(r3$estimates$se, sd(dy) * sqrt((n2 - 1) / n2) / sqrt(n2), tolerance = 1e-8,
             info = "the se of a mean is the plug-in se")

# the cells before the first event are zero, and a cell that is not estimated stops
expect_equal(combine_cells(res, data.table(target = "pre", G = 5, time = 2, w = 1))$estimate, 0)
expect_error(combine_cells(res, data.table(target = "no", G = 99, time = 2, w = 1)), "not in the data")
expect_error(combine_cells(run(full = TRUE), w), "control_type")

# defenses
expect_error(run(control_type = "strata", weightvar = "x"), "weightvar")
expect_error(run(control_type = "strata", allow_unbalance_panel = TRUE), "balanced")
