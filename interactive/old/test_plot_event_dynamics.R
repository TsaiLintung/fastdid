

dt <- sim_did(100, 10)[["dt"]]

event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                 cohortvar = "G",
                                 control_group = "both", verbose = FALSE)

event_est <- get_event_result(event_panel, variable = "y", 
                              trends = FALSE, mem.clean = FALSE, result_type = "dynamic")

expect_silent(event_est |> plot_event_dynamics(), info = "raw plot")

rm(dt, event_panel, event_est)

dt <- sim_did(100, 10, stratify = TRUE)[["dt"]]

event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                 cohortvar = "G",
                                 covariate_base_stratify = "s",
                                 control_group = "both", verbose = FALSE)

event_est <- get_event_result(event_panel, variable = "y", 
                              trends = FALSE, mem.clean = FALSE, result_type = "dynamic")

expect_silent(event_est |> plot_event_dynamics(), info = "plot with stratify")

rm(dt, event_panel, event_est)

dt <- sim_did(100, 10, stratify = TRUE)[["dt"]]

event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                 cohortvar = "G",
                                 covariate_base_stratify = "s",
                                 control_group = "both", verbose = FALSE)

event_est <- get_event_result(event_panel, variable = "y", 
                              trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")

expect_silent(event_est |> plot_event_dynamics(), info = "plot with stratify and cohorts")

rm(dt, event_panel, event_est)