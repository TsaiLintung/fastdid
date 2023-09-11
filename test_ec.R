
rm(list = ls())
gc()

library(profvis)
library(microbenchmark)

setwd("~/GitHub/EventStudyCode")

source("sim_did.R")
source("source/setup.R")


# simulation ---------------------------------------------------------------------

#test with did
simdt <- sim_did(1000, 10, cov = "int", hetero = "dynamic", balanced = FALSE, second_outcome = FALSE)
dt <- simdt$dt

# new event code ---------------------------------------------------------------------


source("source/preprocess.R")
source("source/estimation.R")
source("source/report.R")

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
                                                   balanced_panel = TRUE,
                                                   never_treat_action = "both")

  event_est <- get_result_dynamic(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE)
  
})


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

