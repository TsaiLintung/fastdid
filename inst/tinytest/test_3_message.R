# setup --------------------------------------------------------

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+02, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

# wrong dt ---------------

expect_warning(fastdid(as.data.frame(dt), timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"),
             info = "data.frame")

dt2 <- data.table::copy(dt)
dt2 <- dt2[!is.infinite(G)]

expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"),
               info = "no never")

expect_error(fastdid(dt[time != 3], timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"),
             info = "missing time")

expect_warning(fastdid(dt[time != 3 | unit != 20], timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y", result_type = "group_time"),
               info = "non-balanced panel, missing")

extra_row <-  dt[unit == 20 & time == 3]
expect_error(fastdid(rbind(dt,extra_row), timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"),
             info = "non-balanced panel, extra")

dt2 <- data.table::copy(dt)
dt2[, x := 1]
expect_error(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                     covariatesvar = "x"),
             info = "covariates with no variation")

dt2 <- data.table::copy(dt)
dt2[G == Inf, x := x+10, by = "unit"]
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                     covariatesvar = "x", control_type = "ipw", control_option = "never"),
             info = "covariates with no overlap")

# warning for problematic dt ----------

dt2 <- data.table::copy(dt)
dt2[unit == 1 & time < 5, x := 3]
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                       covariatesvar = "x"),
               info = "time varying covariates is warned")

dt2 <- data.table::copy(dt)
dt2[unit == 1 & time == 4, time := NA]
dt2[time == 3 & unit > 30, y := NA]
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time"),
               info = "missing values")



# already_treated
dt_at <- data.table::copy(dt)
dt_at <- dt_at[time > min(G)]
expect_warning(fastdid(dt_at, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time"),
               info = "already treated")

dt2 <- data.table::copy(dt)
dt2[, x := as.factor(round(x))]
expect_error(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time", covariatesvar = "x"),
             info = "only numeric covariate")

# wrong arguments -------------------------------------------------

expect_error(fastdid::fastdid(dt, timevar = "time", cohortvar = "g", unitvar = "unit", outcomevar = "zz",  result_type = "group_time"),
             info = "wrong col name")

expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                     boot = FALSE, clustervar = "x"),
             info = "clustered but no boot")

expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                     boot = FALSE, clustervar = "x"),
             info = "clustered but no boot")

expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                     varycovariatesvar = "x2", allow_unbalanced_panel = TRUE),
             info = "unbalanced but vary")