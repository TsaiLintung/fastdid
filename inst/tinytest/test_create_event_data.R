dt <- sim_did(100, 10)[["dt"]]
expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                               cohortvar = "G",
                                               control_group = "both", verbose = FALSE),
              info = "create_event_data base setting")
dt <-  sim_did(100, 10, cov = "int", stratify = TRUE)[["dt"]]
rm(dt, event_panel)


dt <-  sim_did(100, 10)[["dt"]]
expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                                cohortvar = "G",
                                                control_group = "never", verbose = FALSE),
               info = "create_event_data never only")
rm(dt, event_panel)


dt <-  sim_did(100, 10, untreated_prop = 0)[["dt"]]
expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                                cohortvar = "G",
                                                control_group = "later", verbose = FALSE),
               info = "create_event_data later only")
rm(dt, event_panel)


dt <-  sim_did(100, 10, untreated_prop = 0)[["dt"]]
expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                               cohortvar = "G",
                                               covariate_base_balance = "x",
                                               covariate_base_stratify = "s",
                                               control_group = "both", verbose = FALSE),
              info = "create_event_data cov and balance")
rm(dt, event_panel)


dt <-  sim_did(100, 10, balanced = FALSE)[["dt"]]
expect_warning(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                               cohortvar = "G",
                                               control_group = "both", verbose = FALSE),
              info = "warning when eventdata not balanced")
rm(dt, event_panel)


dt <-  sim_did(100, 10, cov = "int", stratify = TRUE)[["dt"]]
expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                                cohortvar = "G",
                                                covariate_base_balance_linear = "x",
                                                covariate_base_balance_linear_subset  = "s",
                                                control_group = "both", verbose = FALSE),
               info = "create_event_data covariate_base_balance_linear")
rm(dt, event_panel)

dt <-  sim_did(100, 10)[["dt"]]
expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                               cohortvar = "G",
                                               balanced_panel = FALSE,
                                               control_group = "both", verbose = FALSE),
              info = "create_event_data no balance panel so checking common supprt")
rm(dt, event_panel)

dt <-  sim_did(100, 10)[["dt"]]
expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                               cohortvar = "G",
                                               balanced_panel = FALSE,
                                               control_group = "both", verbose = FALSE, combine = FALSE),
              info = "create_event_data no combine")
rm(dt, event_panel)

dt <-  sim_did(100, 10)[["dt"]]
dt[, h := TRUE]
dt[, g := TRUE]
dt[, j := TRUE]
expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                               cohortvar = "G",
                                               check_not_treated = TRUE,
                                               stratify_by_cohort = TRUE,
                                               lower_event_time = -7,
                                               upper_event_time = 7,
                                               min_control_gap = 2,
                                               max_control_gap = 9,
                                               treat_criteria = "h",
                                               base_restrict = "g",
                                               base_restrict_treat = "j",
                                               control_group = "both", verbose = FALSE),
              info = "create_event_data no catchall for arguments, at least no error")
rm(dt, event_panel)
