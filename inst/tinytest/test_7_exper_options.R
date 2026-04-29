# setup --------------------------------------------------------

simdt <- sim_did(1e+02, 10, cov = "cont", hetero = "all", balanced = TRUE, seed = 1)
dt <- simdt$dt

# baseline dynamic result for comparison
base <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                result_type = "dynamic")

# only_est_min / only_est_max: validation errors -------------------------

expect_error(
  fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          result_type = "group_time", exper = list(only_est_min = -3)),
  info = "only_est_min requires result_type == dynamic"
)

expect_error(
  fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          result_type = "group_time", exper = list(only_est_max = 3)),
  info = "only_est_max requires result_type == dynamic"
)

expect_error(
  fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          result_type = "dynamic", exper = list(only_est_min = 2, only_est_max = -2)),
  info = "only_est_min > only_est_max should error"
)

# only_est_min / only_est_max: double DiD error --------------------------

simdt2 <- sim_did(1e+02, 10, second_cohort = TRUE, seed = 1)
dt2 <- simdt2$dt

expect_error(
  fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          cohortvar2 = "G2", result_type = "dynamic",
          exper = list(only_est_min = -3)),
  info = "only_est_min cannot be used with double DiD"
)

# only_est_min / only_est_max: estimates within range are unchanged ------

# restrict to event times [-3, 3]
res_min_max <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                       result_type = "dynamic",
                       exper = list(only_est_min = -3, only_est_max = 3))

in_range <- base[event_time >= -3 & event_time <= 3]
setorder(in_range, event_time)
setorder(res_min_max, event_time)

expect_equal(res_min_max$att, in_range$att, tolerance = 1e-10,
             info = "att within range matches baseline when using only_est_min/only_est_max")

# result contains only event times within the requested range
expect_true(all(res_min_max$event_time >= -3 & res_min_max$event_time <= 3),
            info = "only event times in range are returned")

expect_true(!all(base$event_time >= -3 & base$event_time <= 3),
            info = "baseline contains event times outside range (sanity check)")

# only_est_min only
res_min <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                   result_type = "dynamic",
                   exper = list(only_est_min = 0))

in_range_min <- base[event_time >= 0]
setorder(in_range_min, event_time)
setorder(res_min, event_time)

expect_equal(res_min$att, in_range_min$att, tolerance = 1e-10,
             info = "att within range matches baseline when using only_est_min alone")

# only_est_max only
res_max <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                   result_type = "dynamic",
                   exper = list(only_est_max = 2))

in_range_max <- base[event_time <= 2]
setorder(in_range_max, event_time)
setorder(res_max, event_time)

expect_equal(res_max$att, in_range_max$att, tolerance = 1e-10,
             info = "att within range matches baseline when using only_est_max alone")
