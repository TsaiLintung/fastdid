# setup --------------------------------------------------------

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+02, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

# did plot ---------------

#full result
fd_full <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "dynamic", full = TRUE)
expect_inherits(plot(fd_full), "ggplot", info = "from full plot")

fd<- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "dynamic")
expect_inherits(plot(fd), "ggplot", info = "dynamics plot")

fd<- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = c("y", "y2"),  result_type = "dynamic")
expect_inherits(plot(fd), "ggplot", info = "multi y plot")

fd<- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group")
expect_inherits(plot(fd), "ggplot", info = "group plot")

fd<- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "time")
expect_inherits(plot(fd), "ggplot", info = "time plot")

fd<- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y")
expect_error(plot(fd), info = "no plot for group-time")

#old syntax
fd<- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "dynamic")
expect_inherits(plot_did_dynamics(fd), "ggplot", info = "old function")

# diagnosis ---------------

simdt <- sim_did(1e+03, 10, cov = "no", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 stratify = FALSE, second_cohort = TRUE, confound_ratio = 1)
dt <- simdt$dt
diag <- diagnose_confound_event(dt, "time", "G", "G2")
expect_inherits(plot(diag), "ggplot", "plot diag")