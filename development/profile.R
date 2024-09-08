library(profvis)

setDTthreads(0)
options(mc.cores = getDTthreads())

#simple
simdt <- sim_did(1e+06, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

profvis({
  result <- fastdid(dt, 
                    timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time")
})

gc()

#simple parallel
simdt <- sim_did(1e+06, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

profvis({
  result <- fastdid(dt, 
                    timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time", parallel = TRUE)
})

#double
simdt <- sim_did(1e+06, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt

profvis({
  result <- fastdid(dt, 
                    timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                    cohortvar2 = "G2")
})

#double parallel
simdt <- sim_did(1e+06, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE)
dt <- simdt$dt

profvis({
  result <- fastdid(dt, 
                    timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time")
})