
rm(list = ls())
gc()

library(profvis)
library(microbenchmark)

setwd("~/GitHub/EventStudyCode")

source("sim_did.R")

# simulation ---------------------------------------------------------------------

#test with did

source("source/setup.R")
source("source/preprocess.R")
source("source/estimation.R")
source("source/report.R")

simdt <- sim_did(100000, 10, cov = "int", hetero = "dynamic", balanced = FALSE)
dt <- simdt$dt

# event code ---------------------------------------------------------------------

profvis({
  
  event_panel <- copy(dt) #copying so that the original does not change
  
  min_time <- -Inf
  max_time <- Inf
  y_name <- "y"
  t_name <- "time"
  unit_name <- "unit"
  cohort_name <- "G"
  balance_covariate <- "x"

  event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                   cohortvar = cohort_name,
                                                   covariate_base_balance = balance_covariate,
                                                   balanced_panel = TRUE,
                                                   never_treat_action = "both")
  
  event_panel <- construct_event_variables(event_panel)
  
  event_est <- get_result_dynamic(event_panel, variable = y_name, trends = FALSE)
  
})


profvis({

  event_est <- get_result_dynamic(event_panel, variable = y_name)
  
})



att_comp <- validate_att_est(simdt$att, event_code_est$att, event_code_est$se, type = "dynamic")

plot_event_study(event_code_est)

#did -------------------------------------------------------------------------------

est <- att_gt(yname = "y",
              tname = "time",
              idname = "unit",
              gname = "G",
              xformla = ~x,
              data = dt)

ratio <- validate_att_est(simdt$att, est$att, est$se)

