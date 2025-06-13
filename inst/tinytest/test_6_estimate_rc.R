# setup --------------------------------------------------------

simdt <- sim_did(100, 5, cov = "cont", hetero = "all", balanced = TRUE, seed = 1)
dt <- simdt$dt

# drop all controls in time 3 to trigger insufficient control observations
no_control <- dt[!(time == 3 & G >= 3)]
expect_warning(
  fastdid(no_control, timevar = "time", cohortvar = "G", unitvar = "unit",
          outcomevar = "y", result_type = "group_time", allow_unbalance_panel = TRUE),
  info = "rc did insufficient control"
)
