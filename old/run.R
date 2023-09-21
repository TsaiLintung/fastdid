
rm(list = ls())
gc()

library(profvis)
library(ggplot2)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("source/sim_did.R")
source("source/utils.R")
source("source/setup.R")
source("source/get_event_result.R")
source("source/create_event_data.R")
source("source/plot_event_dynamics.R")

# simple ---------------------------------------------------------------------

simdt <- sim_did(100000, 10, cov = "no", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = FALSE)
dt <- simdt$dt

profvis({
event_panel_list <- dt %>% create_event_data(timevar = "time", unitvar =  "unit", 
                                        cohortvar = "G",
                                        #covariate_base_balance = "x",
                                        #covariate_base_stratify = "s",
                                        balanced_panel = TRUE,
                                        control_group = "both", copy = FALSE, combine = FALSE)
#}) 
event_est <- get_event_result(event_panel_list, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")
})

event_est |> plot_event_dynamics()

event_est_ce <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")


