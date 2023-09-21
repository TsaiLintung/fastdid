
rm(list = ls())
gc()

library(profvis)

#set the wd to the source folder
setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("source/setup.R")

source("source/sim_did.R")
source("source/utils.R")
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

event_est <- get_event_result(event_panel, variable = c("y"), trends = FALSE, mem.clean = FALSE, result_type = "dynamic")

event_est |> plot_event_dynamics()

# time compared to old event code ------------------------------------------------------

period <- 10
sample_size <- 100000

#new event code
simdt <- sim_did(sample_size, period, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

started.at <- proc.time()

event_panel <- dt %>% create_event_data(timevar = "time", unitvar =  "unit", 
                                        cohortvar = "G",
                                        covariate_base_balance = "x",
                                        covariate_base_stratify = "s",
                                        balanced_panel = TRUE,
                                        control_group = "both", copy = FALSE, combine = FALSE)
event_est <- get_event_result(event_panel, variable = c("y"), trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")

ec_time <- timetaken(started.at)

gc()

#the old event code
source("source_raw/eventcode_IV.R")

simdt <- sim_did(sample_size, period, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

started.at <- proc.time()

event_panel <- suppressWarnings(create_event_data_old(dt, timevar = "time", unitvar =  "unit", 
                                                      cohortvar = "G",
                                                      covariate_base_balance = "x",
                                                      covariate_base_stratify = "s",
                                                      balanced_panel = FALSE,
                                                      never_treat_action = "both"))
event_panel <- construct_event_variables(event_panel)
event_es <- get_result_dynamic(event_panel, variable = "y", trends = FALSE)

old_ec_time <- timetaken(started.at)

message("old time: ", old_ec_time)
message("new time: ", ec_time)

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

# no combine estimation -----------------------------------------------------------------

#by deferring binding the cohorts, both time and memory can be saved, but it also means only cohort_event_time can be estimated easily. 

simdt <- sim_did(100000, 10, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

#combine as usual
default_ram_usage <- peakRAM({
  event_panel <- dt %>% create_event_data(timevar = "time", unitvar =  "unit", 
                                          cohortvar = "G",
                                          covariate_base_balance = "x",
                                          covariate_base_stratify = "s",
                                          balanced_panel = TRUE,
                                          control_group = "both", copy = FALSE)
  
  event_est <- get_event_result(event_panel, variable = c("y"), trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")
})

gc()

#old event code
simdt <- sim_did(sample_size, period, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

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
#no combine saves about 40% peak ram used
message("New peak ram: ", default_ram_usage$Peak_RAM_Used_MiB) #2206
message("Old peak ram: ", old_ram_usage$Peak_RAM_Used_MiB) #7998

