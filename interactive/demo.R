
rm(list = ls())
gc()

library(profvis)
library(devtools)
library(peakRAM)

#set the wd to the source folder
setwd("~/GitHub/EventStudyCode")

load_all()

# load event code ---------------------------------------------------------------------

# simple ---------------------------------------------------------------------

simdt <- sim_did(1000, 10, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

event_panel <- dt %>% create_event_data(timevar = "time", unitvar =  "unit", 
                                        cohortvar = "G",
                                        covariate_base_balance = "x",
                                        covariate_base_stratify = "s",
                                        balanced_panel = TRUE,
                                        control_group = "both", copy = FALSE)

event_est <- get_event_result(event_panel, variable = c("y"), trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")

event_est |> plot_event_dynamics()

# multiple outcome ----------------------------------------------------------------------

simdt <- sim_did(1000, 10, cov = "int", hetero = "all", balanced = FALSE, second_outcome = TRUE, seed = 1, stratify = TRUE)
dt <- simdt$dt

event_panel <- dt %>% create_event_data(timevar = "time", unitvar =  "unit", 
                                        cohortvar = "G",
                                        covariate_base_balance = "x",
                                        covariate_base_stratify = "s",
                                        balanced_panel = TRUE,
                                        control_group = "both", copy = FALSE)

event_est <- get_event_result(event_panel, variable = c("y", "y2"), trends = FALSE, mem.clean = FALSE, result_type = "dynamic")

event_est |> plot_event_dynamics()

# cohort-event time att -----------------------------------------------------------------

simdt <- sim_did(1000, 10, cov = "int", hetero = "all", balanced = FALSE, second_outcome = TRUE, seed = 1, stratify = TRUE)
dt <- simdt$dt

event_panel <- dt %>% create_event_data(timevar = "time", unitvar =  "unit", 
                                        cohortvar = "G",
                                        covariate_base_balance = "x",
                                        covariate_base_stratify = "s",
                                        balanced_panel = TRUE,
                                        control_group = "both", copy = FALSE)

event_est <- get_event_result(event_panel, variable = c("y"), trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")

event_est |> plot_event_dynamics()


# comparison with old event code ------------------------------------------------------

time_period <- 10
sample_size <- 100000

#new event code
simdt <- sim_did(sample_size, time_period, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

started.at <- proc.time()

new_ram_usage <- peakRAM({
  event_panel <- dt %>% create_event_data(timevar = "time", unitvar =  "unit", 
                                          cohortvar = "G",
                                          covariate_base_balance = "x",
                                          covariate_base_stratify = "s",
                                          balanced_panel = TRUE,
                                          control_group = "both", copy = FALSE)
  event_est <- get_event_result(event_panel, variable = c("y"), trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")
  
  ec_time <- timetaken(started.at)
})
gc()

#the old event code
source("interactive/source_raw/eventcode_IV.R")

simdt <- sim_did(sample_size, time_period, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

started.at <- proc.time()

old_ram_usage <- peakRAM({
  event_panel <- suppressWarnings(create_event_data_old(dt, timevar = "time", unitvar =  "unit", 
                                                        cohortvar = "G",
                                                        covariate_base_balance = "x",
                                                        covariate_base_stratify = "s",
                                                        balanced_panel = FALSE,
                                                        never_treat_action = "both"))
  event_panel <- construct_event_variables(event_panel)
  event_es <- get_result_dynamic(event_panel, variable = "y", trends = FALSE)
})
old_ec_time <- timetaken(started.at)

message("new time: ", ec_time) #33.4s
message("New peak ram: ", new_ram_usage$Peak_RAM_Used_MiB) #1534mb
message("old time: ", old_ec_time) #249s
message("Old peak ram: ", old_ram_usage$Peak_RAM_Used_MiB) #9783mb



