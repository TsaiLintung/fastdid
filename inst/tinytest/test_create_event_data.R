dt <-  sim_did(5, 5, seed = 1, treatment_assign = "uniform", stratify = FALSE)[["dt"]]

expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                               cohortvar = "G",
                                               covariate_base_balance = "x",
                                               covariate_base_stratify = "s",
                                               control_group = "both", verbose = FALSE))

dt <-  sim_did(5, 5, seed = 1, treatment_assign = "uniform", stratify = TRUE)[["dt"]]

expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                               cohortvar = "G",
                                               covariate_base_balance = "x",
                                               covariate_base_stratify = "s",
                                               control_group = "both", verbose = FALSE))