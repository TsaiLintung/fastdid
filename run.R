
rm(list = ls())
gc()

library(profvis)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("sim_did.R")
source("source/setup.R")
source("source/get_event_result.R")
source("source/create_event_data.R")

# test dynamic result ---------------------------------------------------------------------

simdt <- sim_did(1000, 10, cov = "int", hetero = "dynamic", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

min_time <- -Inf
max_time <- Inf
y_name <- c("y")
t_name <- "time"
unit_name <- "unit"
cohort_name <- "G"
balance_name <- "x"


event_panel <- dt %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                        cohortvar = cohort_name,
                                        covariate_base_balance = balance_name,
                                        balanced_panel = TRUE,
                                        control_group = "both", copy = FALSE)


event_est <- get_event_result(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE, result_type = "dynamic")

# stratified var ---------------------------------------------------------------------------------

simdt <- sim_did(100000, 10, cov = "int", hetero = "dynamic", balanced = FALSE, second_outcome = FALSE, seed = 1)
dt <- simdt$dt

min_time <- -Inf
max_time <- Inf
y_name <- c("y")
t_name <- "time"
unit_name <- "unit"
cohort_name <- "G"
balance_name <- "x"
stratify_name <- "s"

#event_panel <- copy(dt)
event_panel <- dt

profvis({
  event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                   cohortvar = cohort_name,
                                                   covariate_base_balance = balance_name,
                                                   covariate_base_stratify = stratify_name,
                                                   balanced_panel = TRUE,
                                                   control_group = "both")
})

profvis(
  event_est <- get_event_result(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE, result_type = "dynamic", separate_stratify = FALSE)
)