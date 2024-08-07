#setwd("~/Documents/GitHub/fastdid")

library(here)
library(devtools)
library(tinytest)
library(roxygen2)
library(profvis)

setwd(here())

load_all()

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+05, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt

result <-fastdid(data = dt, 
                 timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                 result_type = "group_time", exper = list(cohortvar2 = "G2", event_specific = TRUE))


tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+04, 20, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt
profvis(result <-fastdid(data = dt, 
                         timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                         result_type = "group_time", exper = list(cohortvar2 = "G2", event_specific = TRUE)), interval = 0.005)