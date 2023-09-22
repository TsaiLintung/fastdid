
rm(list = ls())
gc()

library(peakRAM)
library(microbenchmark)
library(profvis)
library(devtools)

load_all()

setwd("~/GitHub/EventStudyCode")
library(DiDforBigData)

# ------------------------------------------------------------------------------------ 

# functions for individual package call ---------------------------------------------

run_event_code <- function(sample_size, time_period){
  dt <- dt <- sim_did(sample_size, time_period)[["dt"]]
  event_panel <- suppressWarnings(create_event_data(dt, timevar = "time", unitvar =  "unit", cohortvar = "G",
                                          control_group = "both", copy = FALSE, verbose = FALSE, combine = FALSE))
  event_est_ce <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")
}

run_dfbd <- function(sample_size, time_period){
  dt <- dt <- sim_did(sample_size, time_period)[["dt"]]
  
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

t <- 10

all_bm <- data.table()
for(order in seq(2,6)){
  s = 10^order
  bm_time <- microbenchmark(
    run_event_code(s,t),
    #run_did(s,t),
    #run_old_event_code(s,t),
    run_dfbd(s,t),
    times = 1, setup = gc())
  
  bm_time <- bm_time |> as.data.table()
  bm_time <- bm_time[, .(time_mean = mean(time)/10^9), by = "expr"]
  bm_time[, order := order]
  message(s)
  all_bm <- rbind(all_bm, bm_time)
}

all_bm |> fwrite("old/profile_log/est_time.csv")
all_bm |> ggplot(aes(x = order, y = time_mean, color = expr)) + geom_point() + geom_line() +
  labs(title = "running time, event code v.s. DiDForBigData")
ggsave("old/profile_log/est_time_comp.png")

