
rm(list = ls())
gc()

library(profvis)
library(bench)
library(data.table)
library(sandwich)
library(dplyr)
library(stringr)
library(fixest)
library(kit)

setwd("~/GitHub/EventStudyCode")

source("sim_did.R")
source("source/eventcode.R")

#load functions

#library(did)
#library(DiDforBigData)


#source("~/GitHub/EventStudyCode/source_raw/EventCode/eventcode_helper.R")
#source("~/GitHub/EventStudyCode/source_raw/EventCode/eventcode_revised_MaxLouis_ver7.R")



# simulation ---------------------------------------------------------------------

#test with did

source("source/eventcode_factor.R")

simdt <- sim_did(100000, 10, cov = "int", hetero = "dynamic")
dt <- simdt$dt

# event code ---------------------------------------------------------------------

profvis({
  event_panel <- copy(dt) #copying so that the original does not change
  
  min_time <- -8
  max_time <- 8
  y_name <- "y"
  t_name <- "time"
  unit_name <- "unit"
  cohort_name <- "G"
  balance_covariate <- "x"
  
  onset <- event_panel[, min(get(cohort_name)) - 1]
  event_panel[, onset_time := onset]
  
  event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                   cohortvar = cohort_name,
                                                   onset_agevar = "onset_time",
                                                   covariate_base_balance = balance_covariate,
                                                   covariate_base_stratify = c(),
                                                   never_treat_action = "both",
                                                   balanced_panel = FALSE)
  
  event_panel <- event_ATTs_head(event_panel, yname)
  
  event_code_est <- suppressMessages(get_result_dynamic(event_panel,min_time,max_time,y_name,trends = FALSE))
})

dynamic_att <- event_code_est[str_starts(variable, "treated"), ]
dynamic_att[, event_time := as.integer(str_remove_all(str_extract(variable, "y(.*?)\\."), "y|\\."))]
dynamic_att <- dynamic_att[, .(att = Estimate, att_se = `Std. Error`, event_time)]
setorder(dynamic_att, event_time)

ratio <- validate_att_est(simdt$att, dynamic_att$att, dynamic_att$att_se, type = "dynamic")

#did -------------------------------------------------------------------------------

est <- att_gt(yname = "y",
              tname = "time",
              idname = "unit",
              gname = "G",
              xformla = ~x,
              data = dt)

ratio <- validate_att_est(simdt$att, est$att, est$se)

