
rm(list = ls())
gc()

library(profvis)
library(bench)
library(data.table)
library(sandwich)
library(dplyr)
library(stringr)
library(fixest)

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

simdt <- sim_did(2000, 10, cov = "int", hetero = "dynamic")
dt <- simdt$dt

# event code ---------------------------------------------------------------------

started.at <- proc.time()
event_code_est <- estimate_event_dynamics(dt, start = -8, end = 8, 
                                          outcomes = "y", unitvar = "unit", timevar = "time", cohortvar = "G", control = "x",
                                          use_never_treat = TRUE)
timetaken(started.at)

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

