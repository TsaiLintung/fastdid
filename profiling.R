
rm(list = ls())
gc()

library(profvis)
library(microbenchmark)
library(data.table)
library(sandwich)
library(dplyr)

setwd("~/GitHub/EventStudyCode")

source("source_raw/DiDforBigData/R/DiD_simulator.R")

sim <- SimDiD(sample_size = 10000)
dt <- sim$simdata
true_ATT <- dt$true_ATT

# DiDforBigData --------------------------------------------------------------

library(DiDforBigData)

source("source_raw/DiDforBigData/R/DiD_combine_cohorts.R")
source("source_raw/DiDforBigData/R/DiD_within_cohort.R")
source("source_raw/DiDforBigData/R/Utils.R")
source("source_raw/DiDforBigData/R/DiD_SEs.R")

varnames = list()
varnames$time_name = "year"
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"

profvis(DiD(dt, varnames))

# EventCode --------------------------------------------------------------------

source("source_raw/EventCode/eventcode_Max_ver7_raw.R")
source("source_raw/EventCode/eventcode_helper.R")

profvis(estimate <- estimate_event_dynamics(dt, start = -2, end = 2, outcomes = "Y", unitvar = "id", timevar = "year", cohortvar = "cohort", use_never_treat = FALSE)
)

microbenchmark(
  estimate <- DiD(dt, varnames),
  estimate <- estimate_event_dynamics(dt, start = -9, end = 6, outcomes = "Y", unitvar = "id", timevar = "year", cohortvar = "cohort", use_never_treat = TRUE),
  times = 3
)
