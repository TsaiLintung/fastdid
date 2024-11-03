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

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", alpha = 0.01),
              info = "alternative alpha")

#get full result
full_res <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = c("y", "y2"),  result_type = "group_time", full = TRUE)
expect_equal(names(full_res), c("call", "estimate", "gt_estimate", "agg_inf_func", "agg_weight_matrix"),
             info = "full result")

units <- dt[, unique(unit)]
weights <- data.table::data.table(unit = units, w = rnorm(length(units), 1, 1))
dt2 <- dt |> merge(weights, by = "unit")
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      weightvar = "w"),
              info = "weighted")

# bootstrap part ------------------------

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      boot = TRUE),
              info = "bootstrap")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      boot = TRUE, cband = TRUE),
              info = "uniform confidence interval")

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
                      anticipation = 2),
              info = "anticipation")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", allow_unbalance_panel = TRUE,
                      base_period = "varying"),
              info = "baseperiod vary")

# dt that needs adjustment ---------------------------

base_result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time")

expect_equal(fastdid(dt[nrow(dt):1], timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"), base_result,
              info = "reversed")

dt2 <- data.table::copy(dt)
dt2 <- dt2[G < 4]
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", control_option = "notyet"),
              info = "limited treated group")


dt2 <- data.table::copy(dt)
dt2[, time := time*2 + 3]
dt2[, G := G*2 + 3]
base_result2 <- data.table::copy(base_result)
base_result2[, cohort := cohort*2 + 3]
base_result2[, time := time*2 + 3]

expect_equal(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time"), 
             base_result2,
             info = "time offset")

# unbalanced panel ----------------------------------------------------

dt2 <- data.table::copy(dt)
keep <- sample(c(rep(TRUE, 19),FALSE), dt2[,.N], TRUE)
dt2 <- dt2[keep]


expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", 
                      allow_unbalance_panel = TRUE),
              info = "balance panel false but dt is balance")

# internal behavor -------------------------

dt2 <- data.table::copy(dt)
data.table::setnames(dt2, c("time", "G", "unit"), c("t", "g", "u"))
expect_silent(fastdid(dt2, timevar = "t", cohortvar = "g", unitvar = "u",outcomevar = "y",  result_type = "group_time", alpha = 0.01),
              info = "other column names")

#make sure dt is not copied if copy == FALSE
#if(at_home()){
  dtc <- data.table::copy(dt)
  address <- tracemem(dtc) |> stringr::str_remove_all(">|<")
  out <- capture.output(fastdid(dtc, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time", copy = FALSE))
  #if a copy happen there will be a message at the start, so out[1] won't be the header. 
  expect_true(!stringr::str_detect(stringr::str_flatten_comma(out), address),
              info = "no unintentional copy")
#}

if(.Platform$OS.type == "unix" & at_home() & requireNamespace("parallel")){
  expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", parallel = TRUE),
                info = "parallel")
}
