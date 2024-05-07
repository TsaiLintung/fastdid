# setup --------------------------------------------------------

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+02, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

# basics ---------------------------------------------

# fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
#         control_type = "ipw",
#         covariatesvar = c("x", "x2"))

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"),
              info = "simple call")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      control_type = "ipw",
                      covariatesvar = c("x", "x2")),
              info = "covariates ipw")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      control_type = "reg",
                      covariatesvar = c("x", "x2")),
              info = "covariates reg")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      control_type = "dr",
                      covariatesvar = c("x", "x2")),
              info = "covariates dr")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      control_type = "dr",
                      covariatesvar = c("x", "x2"), varycovariatesvar = "xvar"),
              info = "with varying covariates dr")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      control_type = "dr",
                      varycovariatesvar = "xvar"),
              info = "varying covariates only dr")

expect_silent(fastdid(dt[G != 3], timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"),
              info = "missing cohort")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      boot = TRUE),
              info = "bootstrap")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      boot = TRUE, clustervar = "x"),
              info = "bootstrap clustered")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = c("y", "y2"),  result_type = "group_time",
                      boot = TRUE),
              info = "multiple outcome")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = c("y", "y2"),  result_type = "group_time",
                      boot = TRUE, covariatesvar = "x"),
              info = "multiple outcome with covariates")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", allow_unbalance_panel = TRUE),
              info = "balance panel false but dt is balance")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", allow_unbalance_panel = TRUE,
                      max_control_cohort_diff = 2),
              info = "max control cohort diff")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", allow_unbalance_panel = TRUE,
                      min_control_cohort_diff = 4),
              info = "min control cohort diff")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", allow_unbalance_panel = TRUE,
                      anticipation = 2),
              info = "anticipation")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", allow_unbalance_panel = TRUE,
                      base_period = "varying"),
              info = "baseperiod vary")


# dt that needs adjustment ---------------------------

base_result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time")

expect_equal(fastdid(dt[nrow(dt):1], timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"), base_result,
              info = "reversed")

dt2 <- copy(dt)
dt2[, time := time*2 + 3]
dt2[, G := G*2 + 3]
base_result2 <- copy(base_result)
base_result2[, cohort := cohort*2 + 3]
base_result2[, time := time*2 + 3]

expect_equal(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"), 
             base_result2,
             info = "time offset")

# unbalanced panel ----------------------------------------------------

dt2 <- copy(dt)
keep <- sample(c(rep(TRUE, 19),FALSE), dt2[,.N], TRUE)
dt2 <- dt2[keep]


expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", 
                      allow_unbalance_panel = TRUE),
              info = "balance panel false but dt is balance")

# throw error / warning at problematic dt -------------------------

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

dt2 <- copy(dt)
dt2[, x := 1]
expect_error(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                       covariatesvar = "x"),
               info = "covariates with no variation")

dt2 <- copy(dt)
dt2[unit == 1 & time < 5, x := 3]
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                     covariatesvar = "x"),
             info = "time varying covariates is warned")

dt2 <- copy(dt)
dt2[unit == 1 & time == 4, time := NA]
dt2[time == 3 & unit > 30, y := NA]
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time"),
               info = "missing values")
