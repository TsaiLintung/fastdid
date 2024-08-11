# setup --------------------------------------------------------

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+02, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

# message ---------------

expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                     boot = FALSE, clustervar = "x"),
             info = "clustered but no boot")

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
dt2[unit == 1 & time < 5, x := 3]
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                       covariatesvar = "x"),
               info = "time varying covariates is warned")

dt2 <- data.table::copy(dt)
dt2[unit == 1 & time == 4, time := NA]
dt2[time == 3 & unit > 30, y := NA]
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time"),
               info = "missing values")

expect_error(fastdid::fastdid(dt, timevar = "time", cohortvar = "g", unitvar = "unit", outcomevar = "zz",  result_type = "group_time"),
            info = "wrong col name")

# already_treated
dt_at <- data.table::copy(dt)
dt_at <- dt_at[time > min(G)]
expect_warning(fastdid(dt_at, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time"),
               info = "already treated")