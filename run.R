
rm(list = ls())
gc()

library(profvis)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("sim_did.R")
source("source/utils.R")
source("source/setup.R")
source("source/get_event_result.R")
source("source/create_event_data.R")
source("source/plot_event_dynamics.R")


# simple ---------------------------------------------------------------------

simdt <- sim_did(1000, 10, cov = "int", hetero = "dynamic", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 na = "none")
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

cohort_pop <- dt[, .(pop = .N), by = "cohort"]

event_est_ce <- get_event_result(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")

event_est_ce <- event_est_ce |> merge(cohort_pop, by = "cohort")

event_est_ce_mean <- event_est_ce[, .(weighted_Estimate = sum(Estimate*pop)/sum(pop)), by = "event_time"]

event_est <- get_event_result(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE, result_type = "dynamic")

est_both <- merge(event_est, event_est_ce_mean, by = "event_time")
est_both[, est_diff := Estimate - weighted_Estimate]



# full  ---------------------------------------------------------------------------------

simdt <- sim_did(1000, 10, cov = "int", hetero = "dynamic", balanced = FALSE, second_outcome = TRUE, seed = 1, stratify = TRUE)
dt <- simdt$dt

min_time <- -Inf
max_time <- Inf
y_name <- c("y", "y2")
t_name <- "time"
unit_name <- "unit"
cohort_name <- "G"
balance_name <- "x"
stratify_name <- "s"

#event_panel <- copy(dt)
event_panel <- dt

event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                 cohortvar = cohort_name,
                                                 covariate_base_balance = balance_name,
                                                 covariate_base_stratify = stratify_name,
                                                 balanced_panel = TRUE,
                                                 control_group = "both")

event_est <-  get_event_result(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE, result_type = "dynamic")

event_est |> plot_event_dynamics()