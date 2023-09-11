
rm(list = ls())
gc()

library(profvis)
library(microbenchmark)
library(testthat)

setwd("~/GitHub/EventStudyCode")

# load event code ---------------------------------------------------------------------

source("sim_did.R")

source("source/setup.R")
source("source/preprocess.R")
source("source/estimation.R")
source("source/report.R")

# simulation ---------------------------------------------------------------------


time_period <- 10
sample_size <- 1000

simdt <- sim_did(sample_size, time_period, cov = "int", hetero = "dynamic", balanced = FALSE, second_outcome = FALSE, seed = 1)
dt <- simdt$dt
att <- simdt$att

  
event_panel <- copy(dt) #copying so that the original does not change

min_time <- -Inf
max_time <- Inf
y_name <- c("y")
t_name <- "time"
unit_name <- "unit"
cohort_name <- "G"
stratify_name <- "s"
balance_name <- "x"

event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                 cohortvar = cohort_name,
                                                 covariate_base_balance = balance_name,
                                                 covariate_base_stratify = stratify_name,
                                                 balanced_panel = TRUE,
                                                 never_treat_action = "both")

dynamic_est <- get_result_dynamic(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)

pooled_est <- get_result_pooled(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)
means_est <- get_result_means(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)
event_est <- get_result_covariates(event_panel, covariate = stratify_name, variable = y_name, trends = FALSE, mem.clean = FALSE)

#half way
att[, event_time := time-G]
att_dynamic <- att[!event_time %in% c(-1,time_period-1), .(attgt = mean(attgt)), by = "event_time"]
event_validate <- merge(dynamic_est, att_dynamic, by = "event_time")



att_comp <- validate_att_est(simdt$att, event_est)

plot_event_study(event_code_est)


# event code iV ----------------------------------------------------------------------------------
 
source("source_raw/eventcode_IV.R")

profvis({
  
  event_panel <- copy(dt) #copying so that the original does not change
  
  min_time <- -Inf
  max_time <- Inf
  y_name <- c("y")
  t_name <- "time"
  unit_name <- "unit"
  cohort_name <- "G"
  balance_covariate <- "x"
  
  event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                   cohortvar = cohort_name,
                                                   covariate_base_balance = balance_covariate,
                                                   never_treat_action = "both")
  
  event_panel <- event_panel %>% construct_event_variables(event_panel)
  
  event_est <- get_result_dynamic(event_panel, variable = y_name, trends = FALSE)
  
})


#did -------------------------------------------------------------------------------

est <- att_gt(yname = "y",
              tname = "time",
              idname = "unit",
              gname = "G",
              xformla = ~x,
              data = dt)

ratio <- validate_att_est(simdt$att, est$att, est$se)

