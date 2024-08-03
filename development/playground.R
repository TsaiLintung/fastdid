#setwd("~/Documents/GitHub/fastdid")

library(devtools)
library(tinytest)
library(roxygen2)

load_all()


tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+02, 5, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt

result <-fastdid(data = dt, 
                 timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", copy = TRUE,
                 result_type = "group_time", exper = list(cohortvar2 = "G2", filtervar = NA), validate = FALSE)