tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+03, 5, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt

# basic tests ------------------------------------------------------------------

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2"),
              info = "double call")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", event_specific = FALSE),
              info = "double, combined effect")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE),
              info = "double, unbalanced panel")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", control_option = "never"),
              info = "double, only never")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", double_control_option = "never"),
              info = "double call, double never")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", double_control_option = "notyet"),
              info = "double cal, double notyet")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      covariatesvar = "x",
                      cohortvar2 = "G2"),
              info = "double, covariates")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_group_time",
                      cohortvar2 = "G2"),
              info = "double, ggt")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "dynamic_stagger",
                      cohortvar2 = "G2"),
              info = "double, dynamic_sq")

# non-standard data ---------------------------------------------------------

dt2 <- data.table::copy(dt)
dt2 <- dt2[G < 4]
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE, control_option = "notyet"),
              info = "double, limited treated group")


#all G2 > G
dt2 <- data.table::copy(dt)
dt2 <- dt2[G < G2]
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE, control_option = "notyet"),
              info = "double, G2 > G")

#all G > G2
dt2 <- data.table::copy(dt)
dt2 <- dt2[G2 < G]
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE, control_option = "notyet"),
              info = "double, G > G2")


dt2 <- data.table::copy(dt)
keep <- sample(c(rep(TRUE, 15),FALSE), dt2[,.N], TRUE)
dt2 <- dt2[keep]
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE),
              info = "double, unbalanced panel, unbalanced data")

dt2 <- data.table::copy(dt)
dt2[, time := time*2 + 3]
dt2[, G := G*2 + 3]
dt2[, G2 := G2*2 + 3]

expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2"),
             info = "time offset")

if(.Platform$OS.type == "unix" & at_home() & requireNamespace("parallel")){
  expect_silent(fastdid(dt, timevar = "time",cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                        cohortvar2 = "G2", parallel = TRUE),
                info = "parallel double")
}

