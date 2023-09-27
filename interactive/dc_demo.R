
rm(list = ls())
gc()

source("interactive/eventcode.R")

dt <- sim_did(1000, 10, cov = "int", hetero = "all", balanced = FALSE, second_outcome = TRUE, stratify = TRUE)[["dt"]]

event_panel <- create_event_data(dt, timevar = "time", unitvar =  "unit", cohortvar = "G",
                                 covariate_base_stratify = "s",
                                 covariate_base_balance = "x",
                                 control_group = "both")

#this is a bit awkward ;(
cohort_obs <- get_cohort_size(dt, cohortvar = "G", unitvar = "unit")

event_est <- get_event_result(event_panel, variable = c("y", "y2"), 
                              result_type = "dynamic", cohort_obs = cohort_obs)

event_est |> plot_event_dynamics()
