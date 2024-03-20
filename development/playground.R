setwd("~/GitHub/fastdid")

library(devtools)
library(tinytest)
library(roxygen2)

load_all()


tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+02, 4, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt


fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time")

fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_type = "ipw",
                  varycovariatesvar = c("xvar"))

fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
        control_type = "ipw",
        varycovariatesvar = c("x"))

fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
        control_type = "ipw",
        covariatesvar = c("x"))