
rm(list = ls())
gc()

library(profvis)
library(data.table)
library(sandwich)
library(dplyr)
library(stringr)
library(fixest)
library(kit)
library(collapse)
library(doFuture)
library(purrr)
library(ast2ast)

setwd("~/GitHub/EventStudyCode")

source("sim_did.R")

#load functions

#library(did)
#library(DiDforBigData)


#source("~/GitHub/EventStudyCode/source_raw/EventCode/eventcode_helper.R")
#source("~/GitHub/EventStudyCode/source_raw/EventCode/eventcode_revised_MaxLouis_ver7.R")



# simulation ---------------------------------------------------------------------

#test with did


source("source/eventcode.R")
simdt <- sim_did(100000, 10, cov = "int", hetero = "dynamic")
dt <- simdt$dt

setDTthreads(0)
options(kit.nThread = getDTthreads())

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
                                                   outcomevar = y_name,
                                                   cohortvar = cohort_name,
                                                   covariate_base_balance = balance_covariate,
                                                   lower_event_time = min_time,
                                                   upper_event_time = max_time,
                                                   never_treat_action = "both")

  event_code_est <- get_result_dynamic(event_panel, y_name)
  
})

att_comp <- validate_att_est(simdt$att, event_code_est$att, event_code_est$se, type = "dynamic")

#did -------------------------------------------------------------------------------

est <- att_gt(yname = "y",
              tname = "time",
              idname = "unit",
              gname = "G",
              xformla = ~x,
              data = dt)

ratio <- validate_att_est(simdt$att, est$att, est$se)

