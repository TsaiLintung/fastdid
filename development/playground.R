#setwd("~/Documents/GitHub/fastdid")

library(here)
library(devtools)
library(tinytest)
library(roxygen2)
library(profvis)

setwd(here())

load_all()



tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+03, 4, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt


dt2 <- data.table::copy(dt)
dt2 <- dt2[G < 4]


res <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
        cohortvar2 = "G2", control_option = "notyet")