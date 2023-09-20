
rm(list = ls())
gc()

library(peakRAM)
library(microbenchmark)
library(profvis)
library(ggplot2)

setwd("~/GitHub/EventStudyCode")



# load event code ---------------------------------------------------------------------

source("source/sim_did.R")
source("source/utils.R")
source("source/setup.R")
source("source/get_event_result.R")
source("source/create_event_data.R")
source("source/plot_event_dynamics.R")

#old
#source("source_raw/eventcode_IV.r")

library(did)

library(DiDforBigData)

# ------------------------------------------------------------------------------------ 

get_dt <- function(sample_size, time_period){
  simdt <- sim_did(sample_size, time_period, cov = "no", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1, stratify = FALSE)
  dt <- simdt$dt
  return(dt)
}

# functions for individual package call ---------------------------------------------

run_event_code <- function(sample_size, time_period){
  dt <- get_dt(sample_size, time_period)
  event_panel <- suppressWarnings(create_event_data(dt, timevar = "time", unitvar =  "unit", 
                                          cohortvar = "G",
                                          #covariate_base_balance = "x",
                                          covariate_base_stratify = "s",
                                          control_group = "both", copy = FALSE, verbose = FALSE))
  event_est_ce <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "dynamic")
}



run_dfbd <- function(sample_size, time_period){
  dt <- get_dt(sample_size, time_period)
  
  varnames = list()
  varnames$time_name = "time"
  varnames$outcome_name = "y"
  varnames$cohort_name = "G"
  varnames$id_name = "unit"
  # estimate the ATT for all cohorts at event time 1 only
  result <- DiD(dt, varnames)
  return(result)
}

# start benchmarks ----------------------------------------------------------------------

profvis(run_event_code(100000,10))
profvis(run_dfbd(100000,10))

#

t <- 10

all_bm <- data.table()
for(order in seq(2,5)){
  s = 10^order
  bm_time <- microbenchmark(
    run_event_code(s,t),
    #run_did(s,t),
    #run_old_event_code(s,t),
    run_dfbd(s,t),
    times = 3, setup = gc())
  
  bm_time <- bm_time |> as.data.table()
  bm_time <- bm_time[, .(time_mean = mean(time)/10^9,
                         time_sd = sd(time)/10^9), by = "expr"]
  bm_time[, order := order]
  all_bm <- rbind(all_bm, bm_time)
}

all_bm |> ggplot(aes(x = order, y = time_mean, color = expr)) + geom_point() + geom_line()


# not used

run_did <-  function(sample_size, time_period){
  dt <- get_dt(sample_size, time_period)
  did_est <- suppressWarnings(att_gt(yname = "y",
                                     tname = "time",
                                     idname = "unit",
                                     gname = "G",
                                     #xformla = ~x,
                                     data = dt,
                                     bstrap = FALSE))
}


run_old_event_code <-  function(sample_size, time_period){
  dt <- get_dt(sample_size, time_period)
  event_panel <- suppressWarnings(create_event_data_old(dt, timevar = "time", unitvar =  "unit", 
                                                        cohortvar = "G",
                                                        #covariate_base_balance = "x",
                                                        covariate_base_stratify = "s",
                                                        balanced_panel = FALSE,
                                                        never_treat_action = "both"))
  event_panel <- construct_event_variables(event_panel)
  event_es <- get_result_dynamic(event_panel, variable = "y", trends = FALSE)
}




