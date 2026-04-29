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


dt2 <- data.table::copy(dt)
dt2 <- dt2[!is.infinite(G)]
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                       control_option = "both"),
               info = "no never but both")

dt2 <- data.table::copy(dt)
dt2 <- dt2[!is.infinite(G)]
expect_error(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                       control_option = "never"),
               info = "no never but never")

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

# irregular time steps --------------------------------------------

dt_irregular <- data.table::copy(dt)
dt_irregular[time == 5, time := 10]
dt_irregular[time == 6, time := 12]
dt_irregular[time == 7, time := 15]
expect_error(fastdid(dt_irregular, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
             info = "irregular time steps")

# Non-uniform time intervals (e.g., 1, 2, 4, 8) 
dt_irregular2 <- data.table::copy(dt)
dt_irregular2[, time := c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)[time]]
dt_irregular2[, G := ifelse(is.infinite(G), Inf, c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)[G])]
expect_error(fastdid(dt_irregular2, timevar = "time_new", cohortvar = "G_new", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
             info = "exponential time steps")

# balanced_event_time validation ----------------------------------

# balanced_event_time with wrong result_type
expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_time", balanced_event_time = 3),
             info = "balanced_event_time with non-dynamic result_type")

# balanced_event_time that's too large
expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "dynamic", balanced_event_time = 100),
             info = "balanced_event_time exceeds max event time")

# balanced_event_time negative
expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "dynamic", balanced_event_time = -1),
             info = "balanced_event_time negative")

# double DiD result_type validation -------------------------------

# dynamic_stagger without cohortvar2
expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "dynamic_stagger"),
             info = "dynamic_stagger without double DiD")

# group_group_time without cohortvar2
expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_group_time"),
             info = "group_group_time without double DiD")

# empty dataset issues --------------------------------------------

# Empty dataset after filtering
dt_empty <- dt[unit > 1000]
expect_error(fastdid(dt_empty, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
             info = "empty dataset")

# All units dropped due to missing periods
dt_all_missing <- data.table::copy(dt)
dt_all_missing <- dt_all_missing[time != 5]
expect_error(fastdid(dt_all_missing, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
             info = "all units missing periods - should suggest allow_unbalance_panel")

# IPW extreme weights warning -------------------------------------

# Create data with poor overlap to trigger IPW weight warning
dt_extreme <- data.table::copy(dt)
dt_extreme[is.infinite(G), x := x + 100]  # Extreme separation
expect_warning(fastdid(dt_extreme, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                       result_type = "group_time", covariatesvar = "x", control_type = "ipw"),
               info = "extreme IPW weights warning")

# cohort not found issues -----------------------------------------

# This should be caught gracefully with informative warning
dt_missing_cohort <- data.table::copy(dt)
dt_missing_cohort <- dt_missing_cohort[G != 5]  # Remove entire cohort
expect_silent(fastdid(dt_missing_cohort, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
              info = "missing cohort handled gracefully")

# invalid time structure ------------------------------------------

# No time periods
dt_no_time <- dt[time == 1]
dt_no_time[, time := 1]  # All same time
expect_error(fastdid(dt_no_time, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
             info = "no time variation")

# negative anticipation parameters --------------------------------

expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_time", anticipation = -2),
             info = "negative anticipation")

expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_time", anticipation2 = -3),
             info = "negative anticipation2")

# unbalanced panel with DR ----------------------------------------

expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_time", allow_unbalance_panel = TRUE, control_type = "dr"),
             info = "unbalanced panel with DR not supported")

# duplicate unit-time observations --------------------------------

dt_dup <- rbind(dt, dt[1:5])
expect_error(fastdid(dt_dup, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
             info = "duplicate unit-time observations")

# non-numeric variables -------------------------------------------

dt_char <- data.table::copy(dt)
dt_char[, time := as.character(time)]
expect_error(fastdid(dt_char, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
             info = "non-numeric time")

dt_char2 <- data.table::copy(dt)
dt_char2[, G := as.character(G)]
expect_error(fastdid(dt_char2, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
             info = "non-numeric cohort")

# overlapping covariate names -------------------------------------

expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_time", covariatesvar = c("x", "x2"), varycovariatesvar = c("x2", "xvar")),
             info = "overlapping covariate names")

# varname collision -----------------------------------------------

dt_collision <- data.table::copy(dt)
expect_error(fastdid(dt_collision, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_time", weightvar = "x", covariatesvar = "x"),
             info = "duplicate varnames not allowed")

# parallel on non-unix -------------------------------------------

if (.Platform$OS.type != "unix") {
  expect_error(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                       result_type = "group_time", parallel = TRUE),
               info = "parallel on non-unix system")
}

# edge case: single cohort ----------------------------------------

dt_single <- dt[G == 2 | is.infinite(G)]
expect_silent(fastdid(dt_single, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
              info = "single treated cohort works")

# edge case: very few observations per group ----------------------

dt_small <- dt[unit <= 5]
expect_silent(fastdid(dt_small, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
              info = "small sample works")

# test detailed diagnostics in balanced panel warning -------------

dt_detailed <- data.table::copy(dt)
dt_detailed <- dt_detailed[!(unit %in% c(1, 2, 3) & time == 5)]
expect_warning(fastdid(dt_detailed, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
               info = "detailed balanced panel diagnostics")

# zero weights validation -----------------------------------------

dt_zero_weights <- data.table::copy(dt)
dt_zero_weights[, w := 0]
expect_error(fastdid(dt_zero_weights, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_time", weightvar = "w"),
             info = "zero weights should error")


# time varying weight -------------

dt_vary_weights <- data.table::copy(dt)
dt_vary_weights[, w := runif(.N, 0.1, 0.9)]
expect_error(fastdid(dt_vary_weights, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                     result_type = "group_time", weightvar = "w"),
             info = "time-varying weights should error")

# test with time gaps that are uniform ----------------------------

dt_gaps <- data.table::copy(dt)
dt_gaps[, time := time * 2]  # Make time 2, 4, 6, 8, ...
dt_gaps[!is.infinite(G), G := G * 2]
expect_silent(fastdid(dt_gaps, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
              info = "uniform gaps in time periods work")

# test with time starting at non-1 value --------------------------

dt_offset <- data.table::copy(dt)
dt_offset[, time := time + 1990]  # Years like 1991, 1992, ...
dt_offset[!is.infinite(G), G := G + 1990]
expect_silent(fastdid(dt_offset, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
              info = "time offset works")

# Inf and NaN values in numeric columns ---------------------------

dt_inf <- data.table::copy(dt)
dt_inf[unit == 1 & time == 1, x := Inf]
expect_warning(fastdid(dt_inf, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                       result_type = "group_time", covariatesvar = "x"),
               info = "Inf in covariates removed as missing")

dt_nan <- data.table::copy(dt)
dt_nan[unit == 1 & time == 1, x := NaN]
expect_warning(fastdid(dt_nan, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                       result_type = "group_time", covariatesvar = "x"),
               info = "NaN in covariates removed as missing")

# large floating point time values that need coercion -------------

dt_float <- data.table::copy(dt)
dt_float[, time := time + 0.1]
dt_float[!is.infinite(G), G := G + 0.1]
expect_silent(fastdid(dt_float, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
              info = "floating point time coerced to integer")

# test with only two time periods ----------------------------------

dt_two_periods <- dt[time %in% c(1, 2)]
expect_silent(fastdid(dt_two_periods, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
              info = "two time periods work")

# test robust to different data.table key -------------------------

dt_keyed <- data.table::copy(dt)
data.table::setkey(dt_keyed, unit, time)
expect_silent(fastdid(dt_keyed, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", result_type = "group_time"),
              info = "data.table with key works")

# test copy=FALSE option ------------------------------------------

dt_no_copy <- data.table::copy(dt)
expect_silent(fastdid(dt_no_copy, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                      result_type = "group_time", copy = FALSE),
              info = "copy=FALSE works")
