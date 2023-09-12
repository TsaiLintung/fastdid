
rm(list = ls())
gc()

library(tinytest)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("sim_did.R")

source("source/setup.R")
source("source/get_event_result.R")
source("source/create_event_data.R")


# test dynamic result ---------------------------------------------------------------------

simdt <- sim_did(100, 10, cov = "int", hetero = "dynamic", balanced = FALSE, second_outcome = FALSE, seed = 1)
dt <- simdt$dt

min_time <- -Inf
max_time <- Inf
y_name <- c("y")
t_name <- "time"
unit_name <- "unit"
cohort_name <- "G"
stratify_name <- "s"
balance_name <- "x"

#event_panel <- copy(dt)
event_panel <- dt

event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                        cohortvar = cohort_name,
                                        covariate_base_balance = balance_name,
                                        covariate_base_stratify = stratify_name,
                                        balanced_panel = TRUE,
                                        never_treat_action = "both")

dynamic_est <- get_result_dynamic(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)

ce <- get_result_cohort_event_time(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)