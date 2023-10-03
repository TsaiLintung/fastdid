
rm(list = ls())
gc()

library(profvis)
library(devtools)
library(peakRAM)
library(microbenchmark)
library(DRDID)

library(fastglm)
library(data.table)
library(collapse)
library(speedglm)
library(RcppArmadillo)
library(BMisc)


#set the wd to the source folder
setwd("~/GitHub/EventStudyCode")

load_all()
#load_all("~/GitHub/did")
library(did)


# load event code ---------------------------------------------------------------------

# simple ---------------------------------------------------------------------

simdt <- sim_did(1e+03, 5, cov = "no", hetero = "dynamic", balanced = TRUE, second_outcome = FALSE, seed = 333, stratify = FALSE)
dt <- simdt$dt

# started.at <- proc.time()
# results <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", result_type = "group_time")
# timetaken(started.at)



result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", result_type = "group_time", boot = TRUE,
                  clustervar = "G"
                  )


