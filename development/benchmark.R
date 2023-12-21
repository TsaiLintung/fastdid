rm(list = ls())
gc()

library(devtools)
library(microbenchmark)
library(peakRAM)
library(profvis)

setwd("~/GitHub/EventStudyCode")

data.table::setDTthreads(0)

load_all()
load_all("~/GitHub/did")

# setup ------------------------------------------------------------------------------------ 

benchmark_fastdid <- function(){
  
  for(order in c(5, 6)){
    for(aup in c(FALSE, TRUE)){
      message("testing dataset of ", 10^order,  " units x 10 period, allow unbalance panel: ", aup)
      raw_dt <- sim_did(10^order, 10, seed = 1, cov = "cont", second_outcome = TRUE, second_cov = TRUE)[["dt"]]
      
      dt <- copy(raw_dt)
      gc()
      bm_time <- microbenchmark(result <- fastdid(dt, 
                                                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                                                  result_type = "group_time", copy = FALSE, allow_unbalance_panel = aup),
                                times = 1)
      message("baseline: ", bm_time$time/10^9)
      
      dt <- copy(raw_dt)
      gc()
      bm_time <- microbenchmark(result <- fastdid(dt, 
                                                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                                                  covariatesvar = c("x", "x2"), control_type = "ipw",
                                                  result_type = "group_time", copy = FALSE, allow_unbalance_panel = aup),
                                times = 1)
      message("ipw: ", bm_time$time/10^9)
      
      dt <- copy(raw_dt)
      gc()
      bm_time <- microbenchmark(result <- fastdid(dt, 
                                                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = c("y", "y2"), 
                                                  covariatesvar = c("x", "x2"), control_type = "ipw",
                                                  result_type = "group_time", copy = FALSE, allow_unbalance_panel = aup),
                                times = 1)
      message("ipw multiple outcome: ", bm_time$time/10^9)
      
      dt <- copy(raw_dt)
      gc()
      bm_time <- microbenchmark(result <- fastdid(dt, 
                                                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                                                  result_type = "group_time", copy = FALSE, allow_unbalance_panel = aup),
                                times = 1)
      message("baseline dynamic: ", bm_time$time/10^9)
    }
  }
}

# run ---------------------------------------------------------------------------------------------------------------

benchmark_fastdid()


