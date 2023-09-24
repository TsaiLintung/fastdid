
rm(list = ls())
gc()

library(peakRAM)
library(microbenchmark)
library(profvis)
library(devtools)

setwd("~/GitHub/EventStudyCode")
load_all()


library(DiDforBigData)
library(did)

# ------------------------------------------------------------------------------------ 

# functions for individual package call ---------------------------------------------

run_event_code <- function(sample_size = 100000, time_period = 10){
  dt <- dt <- sim_did(sample_size, time_period, seed = 1)[["dt"]]
  event_panel <- suppressWarnings(create_event_data(dt, timevar = "time", unitvar =  "unit", cohortvar = "G",
                                          control_group = "both", copy = FALSE, verbose = FALSE))
  event_est <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time")
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

run_did <- function(sample_size, time_period){
  dt <- sim_did(sample_size, time_period)[["dt"]]
  
  results <- att_gt(yname = "y",
                          tname = "time",
                          idname = "unit",
                          gname = "G",
                          data = dt,
                          bstrap = FALSE
  )
  return(results)
}

# start benchmarks ----------------------------------------------------------------------

t <- 10

all_bm <- data.table()
for(order in seq(2,6)){
  s = 10^order
  if(order == 6){
    bm_time <- microbenchmark(
      run_event_code(s,t),
      #run_did(s,t),
      #run_old_event_code(s,t),
      run_dfbd(s,t),
      times = 1, setup = gc())
  } else {
    bm_time <- microbenchmark(
      run_event_code(s,t),
      run_did(s,t),
      #run_old_event_code(s,t),
      run_dfbd(s,t),
      times = 1, setup = gc())
  }

  bm_time <- bm_time |> as.data.table()
  bm_time <- bm_time[, .(time_mean = mean(time)/10^9), by = "expr"]
  bm_time[, order := order]
  message(s)
  all_bm <- rbind(all_bm, bm_time)
}

# results --------------------------------------------------------------------------------

all_bm |> fwrite("interactive/plots/est_time.csv")
all_bm |> ggplot(aes(x = order, y = time_mean, color = expr)) + geom_point() + geom_line() +
  labs(title = "time comparison by size of unique ID")
ggsave("interactive/plots/est_time_comp.png")

