# setup --------------------------------------------------------

simdt <- sim_did(100, 5, cov = "cont", hetero = "all", balanced = TRUE, seed = 1)
dt <- simdt$dt

# drop all controls in time 3 to trigger insufficient control observations
no_control <- dt[!(time == 3 & D == 0)]
expect_warning(
  fastdid(no_control, timevar = "time", cohortvar = "G", unitvar = "unit",
          outcomevar = "y", result_type = "group_time", allow_unbalance_panel = TRUE),
  info = "rc did insufficient control"
)

# keep only a single control unit to force regression failure
control_unit <- no_control[D == 0, unique(unit)][1]
small_dt <- no_control[unit %in% c(control_unit, no_control[D == 1, unique(unit)])]
expect_warning(
  fastdid(small_dt, timevar = "time", cohortvar = "G", unitvar = "unit",
          outcomevar = "y", result_type = "group_time", allow_unbalance_panel = TRUE,
          control_type = "reg", covariatesvar = c("x", "x2")),
  info = "rc did regression fail"
)
