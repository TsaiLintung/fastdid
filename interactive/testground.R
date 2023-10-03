
rm(list = ls())
gc()

library(profvis)
library(devtools)
library(peakRAM)
library(microbenchmark)
library(did)
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

source("R/sim_did.R")
source("interactive/newfuncs.R")


# load event code ---------------------------------------------------------------------

# simple ---------------------------------------------------------------------

simdt <- sim_did(1000000, 10, cov = "no", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, stratify = TRUE)
dt <- simdt$dt

started.at <- proc.time()
profvis(result <- fastdid(dt))


timetaken(started.at)



