# first_stage() and double_did(): the second stage on the estimates and their covariance

dt <- sim_did(2000, 6, seed = 2, second_cohort = TRUE)$dt
dt[, x := unit %% 3]

for (cv in list(NA, "x")) {
  res <- fastdid(dt, "time", "G", "unit", "y", cohortvar2 = "G2", event_specific = FALSE,
                 control_type = "strata", covariatesvar = cv, result_type = "simple", full = TRUE)
  fs <- first_stage(res)
  set.seed(1)
  w <- fs$estimates[kind == "cell"][sample(.N, 40)]
  w <- w[, .(target = sample(c("a", "b", "c"), .N, TRUE), G, G2, time, w = rnorm(.N))]
  wm <- fs$estimates[kind == "mean"][sample(.N, 30)]
  wm <- wm[, c("G", "G2", "time", if (!is.na(cv)) "x"), with = FALSE]
  wm[, `:=`(target = sample(c("a", "c", "d"), .N, TRUE), w = rnorm(.N))]

  # the factored covariance gives the estimates and the covariance of combine_cells()
  ref <- combine_cells(res, w, wm, return_if = TRUE)
  out <- double_did(fs, w_cell = w, w_mean = wm)
  m <- merge(ref$estimates, out, by = "target")
  info <- paste("covariates:", cv)
  expect_equal(m$estimate.y, m$estimate.x, tolerance = 1e-10, info = info)
  expect_equal(m$se.y, m$se.x, tolerance = 1e-8, info = info)
  V_if <- crossprod(ref$inf_func)[out$target, out$target] / nrow(ref$inf_func)^2
  expect_equal(unname(attr(out, "vcov")), unname(V_if), tolerance = 1e-8, info = info)

  # a dense covariance matrix of the cells gives the same targets
  if (is.na(cv)) {
    cells <- fs$estimates[kind == "cell", .(G, G2, time, estimate)]
    shares <- dt[time == 1, .(share = .N / uniqueN(dt$unit)), by = .(G, G2)]
    dense <- as_first_stage(cells, shares, c("G", "G2"), "time", vcov = vcov(fs))
    expect_equal(double_did(dense, event_specific = FALSE, w_cell = w), double_did(fs, w_cell = w),
                 tolerance = 1e-10, info = "the dense path")
  }
}

# the second stage and the aggregation of fastdid() run on its first stage
args <- list(dt, "time", "G", "unit", "y", cohortvar2 = "G2", result_type = "dynamic")
for (ct in c("ipw", "strata")) {
  ref <- do.call(fastdid, c(args, control_type = ct, full = TRUE))
  for (em in list("parallel", ~ factor(e))) {
    fit <- do.call(fastdid, c(args, control_type = ct, effect_model = list(em)))
    out <- double_did(first_stage(ref), effect_model = em)
    info <- paste(ct, deparse(em))
    expect_equal(out$att, fit$att, tolerance = 1e-10, info = info)
    expect_equal(out$se, fit$se, tolerance = 1e-8, info = info)
  }
}

# the cells of another estimator: the cells and the influence functions of the
# units of a 2x2 run, given as plain tables, give the result of fastdid()
fs <- first_stage(ref <- fastdid(dt, "time", "G", "unit", "y", cohortvar2 = "G2", full = TRUE))
cells <- fs$estimates[kind == "cell", .(G, G2, time, estimate)]
shares <- data_keys(fs$shares$G, 1L, c("G", "G2"), "time", fs$time_info)[, share := fs$shares$pg]
unit_cohort <- data_keys(fs$cohorts[fs$unc$cohort], 1L, c("G", "G2"), "time", fs$time_info)
ext <- as_first_stage(cells, shares, c("G", "G2"), "time", inf_func = fs$unc$inf,
                      unit_cohort = unit_cohort)
fit <- fastdid(dt, "time", "G", "unit", "y", cohortvar2 = "G2", result_type = "dynamic")
out <- double_did(ext, result_type = "dynamic")
expect_equal(out$att, fit$att, tolerance = 1e-10, info = "another estimator, influence functions")
expect_equal(out$se, fit$se, tolerance = 1e-8, info = "another estimator, influence functions")
ext_v <- as_first_stage(cells, shares, c("G", "G2"), "time", vcov = vcov(fs))
expect_equal(double_did(ext_v, result_type = "dynamic")$att, fit$att, tolerance = 1e-10,
             info = "another estimator, covariance matrix")
expect_error(double_did(ext_v, boot = TRUE), "bootstrap")

# defenses
bad <- data.table(target = "no", G = 99, G2 = 99, time = 100, w = 1)
expect_error(double_did(fs, w_cell = bad), "not in the data")
expect_error(first_stage(fastdid(dt, "time", "G", "unit", "y")), "full = TRUE")
