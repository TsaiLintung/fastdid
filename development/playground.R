setwd("~/Documents/GitHub/fastdid")

library(devtools)
library(tinytest)
library(roxygen2)

load_all()


tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+05, 30, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

profvis(result <-fastdid(data = dt, 
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                  result_type = "group_time", copy = FALSE))