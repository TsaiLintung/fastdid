setwd("~/GitHub/fastdid")

library(devtools)
library(tinytest)
library(roxygen2)

load_all()

#gen data
tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+02, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt



myfunc <- function(x, #covariates
                   y, w, g, t, control_bool){
  
 return(NULL)
}

fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
        control_type = "reg",
        covariatesvar = c("x", "x2"), or_func = myfunc)

fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time")