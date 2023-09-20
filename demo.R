
rm(list = ls())
gc()

library(profvis)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("source/sim_did.R")
source("source/utils.R")
source("source/setup.R")
source("source/get_event_result.R")
source("source/create_event_data.R")
source("source/plot_event_dynamics.R")

# simple ---------------------------------------------------------------------

simdt <- sim_did(1000, 10, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

event_panel <- dt %>% create_event_data(timevar = "time", unitvar =  "unit", 
                                        cohortvar = "G",
                                        covariate_base_balance = "x",
                                        covariate_base_stratify = "s",
                                        balanced_panel = TRUE,
                                        control_group = "both", copy = FALSE)

event_est_ce <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")
event_est <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "dynamic")

event_est |> plot_event_dynamics()
