# setup --------------------------------------------------------

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+02, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

# did plot ---------------

fd <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "dynamic")
expect_inherits(plot_did_dynamics(fd), "ggplot", info = "plot dynamics")

fd <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group")
expect_inherits(plot_did_dynamics(fd, margin = "group"), "ggplot", info = "plot group")

fd <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "time")
expect_inherits(plot_did_dynamics(fd, margin = "time"), "ggplot", info = "plot time")
