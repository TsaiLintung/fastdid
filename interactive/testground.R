
rm(list = ls())
gc()

library(profvis)
library(devtools)
library(peakRAM)
library(microbenchmark)


#set the wd to the source folder
setwd("~/GitHub/EventStudyCode")

load_all()
load_all("~/GitHub/did")
#library(did)


# load event code ---------------------------------------------------------------------

# simple ---------------------------------------------------------------------

simdt <- sim_did(1e+03, 10, cov = "no", hetero = "dynamic", balanced = TRUE, second_outcome = FALSE, seed = 1, stratify = FALSE,
                 epsilon_size = 1)
dt <- simdt$dt

started.at <- proc.time()
results <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", result_type = "dynamic")
timetaken(started.at)


