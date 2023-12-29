rm(list = ls())
gc()

library(peakRAM)

library(magrittr)
library(ggplot2)
library(fixest)

library(did)
#install_github("TsaiLintung/fastdid")
library(fastdid)
library(DiDforBigData)


data.table::setDTthreads(0)

# setup ------------------------------------------------------------------------------------

min_order <- 2
max_order <- 6

# functions for running the estimation -----------------------------------------------------

run_fastdid <- function(dt){
  result <- fastdid(dt, 
                    timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", 
                    result_type = "group_time", copy = FALSE, allow_unbalance_panel = FALSE)
}

run_did <- function(dt){
  result <- att_gt(yname = "y", tname = "time", idname = "unit", gname = "G", data = dt,
                   control_group = "notyettreated", bstrap = FALSE, cband = FALSE)
}

run_dfbd <- function(dt){
  varnames = list()
  varnames$time_name = "time"
  varnames$outcome_name = "y"
  varnames$cohort_name = "G"
  varnames$id_name = "unit"
  
  # estimate the ATT for all cohorts at event time 1 only
  result <- DiD(dt, varnames)
}

run_benchmark <- function(min_order, max_order){
  all_bm <- data.table()
  for(order in min_order:max_order){
    raw_dt <- sim_did(10^order, 10, seed = 1, cov = "cont", second_outcome = TRUE, second_cov = TRUE)[["dt"]]
    
    gc()
    
    dt <- copy(raw_dt)
    fast_bm <- peakRAM(run_fastdid(dt))
    
    gc()
    
    dt <- copy(raw_dt)
    did_bm <- peakRAM(run_did(dt))
    
    gc()
    
    dt <- copy(raw_dt)
    dfbd_bm <- peakRAM(run_dfbd(dt))
    
    gc()
    
    bm <- rbind(fast_bm, did_bm, dfbd_bm) |> as.data.table()
    bm[, order := order] 
    all_bm <- rbind(all_bm, bm)
    message(order)
  }
  
  return(all_bm)
  
}

bm_result <- run_benchmark(min_order, max_order)

# visualize benchmark result ------------------------------------------------------------------------------------ 

bm_result[Function_Call == "run_fastdid(dt)", pkg := "fastdid"]
bm_result[Function_Call == "run_dfbd(dt)", pkg := "dfbd"]
bm_result[Function_Call == "run_did(dt)", pkg := "did"]

setnames(bm_result, c("Elapsed_Time_sec", "Peak_RAM_Used_MiB"), c("time", "RAM"))

bm_result %>% ggplot(aes(x = order, y = time, color = pkg)) + geom_line() + geom_point() + coord_trans(y = "log10") + 
  ggtitle("Comparison of computing time")

bm_result %>% ggplot(aes(x = order, y = RAM, color = pkg)) + geom_line() + geom_point() + 
  ggtitle("Comparison of peak RAM used")

