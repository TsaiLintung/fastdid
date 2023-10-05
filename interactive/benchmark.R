
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

run_fastdid <- function(dt){
  result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time")
}

run_dfbd <- function(dt){
  varnames = list()
  varnames$time_name = "time"
  varnames$outcome_name = "y"
  varnames$cohort_name = "G"
  varnames$id_name = "unit"
  # estimate the ATT for all cohorts at event time 1 only
  result <- DiD(dt, varnames)
  return(result)
}

run_did <- function(dt){
  results <- att_gt(yname = "y",
                    tname = "time",
                    idname = "unit",
                    gname = "G",
                    data = dt,
                    bstrap = FALSE
  )
  return(results)
}


max_order <- 6
t <- 10

# time benchmarks ----------------------------------------------------------------------

all_bm <- data.table()
for(order in seq(2,max_order)){

  dt <- sim_did(10^order, t, seed = 1)[["dt"]]
  
  if(order == 6){
    bm_time <- microbenchmark(
      run_fastdid(copy(dt)),
      #run_did(s,t),
      #run_old_event_code(s,t),
      run_dfbd(copy(dt)),
      times = 1, setup = gc())
  } else {
    bm_time <- microbenchmark(
      run_fastdid(copy(dt)),
      run_did(copy(dt)),
      #run_old_event_code(s,t),
      run_dfbd(copy(dt)),
      times = 1, setup = gc())
  }

  bm_time <- bm_time |> as.data.table()
  bm_time <- bm_time[, .(time_mean = mean(time)/10^9), by = "expr"]
  bm_time[, order := order]
  message(order)
  all_bm <- rbind(all_bm, bm_time)
}

# results --------------------------------------------------------------------------------

all_bm[expr == "run_fastdid(copy(dt))", package := "fastdid"]
all_bm[expr == "run_did(copy(dt))", package := "did"]
all_bm[expr == "run_dfbd(copy(dt))", package := "dfbd"]

all_bm |> fwrite("interactive/plots/est_time.csv")

all_bm |> ggplot(aes(x = order, y = time_mean, color = package)) + geom_point() + geom_line() +
  labs(title = "Performance comparison - computing time") + 
  xlab("unique ID (log10)") + ylab("seconds") + theme_bw()

ggsave("interactive/plots/est_time_comp.png")

# RAM benchmark -------------------------------------------------

t <- 10
all_bm_ram <- data.table()
for(order in seq(2,max_order)){
  
  raw_dt <- sim_did(10^order, t, seed = 1)[["dt"]]
  gc()
  
  dt <- copy(raw_dt)
  ram_bench <- peakRAM(run_fastdid(dt))
  all_bm_ram <- all_bm_ram |> rbind(data.table(package = "fastdid", peak_ram = ram_bench[4], order = order))
  gc()
  
  dt <- copy(raw_dt)
  ram_bench <- peakRAM(run_dfbd(dt))
  all_bm_ram <- all_bm_ram |> rbind(data.table(package = "dfbd", peak_ram = ram_bench[4], order = order))
  gc()
  
  if(order != 6){
    dt <- copy(raw_dt)
    ram_bench <- peakRAM(run_did(dt))
    all_bm_ram <- all_bm_ram |> rbind(data.table(package = "did", peak_ram = ram_bench[4], order = order))
    gc()
  }
  
  message(order)
}

setnames(all_bm_ram, "peak_ram.Peak_RAM_Used_MiB", "peak_ram")
all_bm_ram |> fwrite("interactive/plots/est_ram.csv")

all_bm_ram |> ggplot(aes(x = order, y = peak_ram, color = package)) + geom_point() + geom_line() +
  labs(title = "Performance comparison - peak RAM used") + 
  xlab("unique ID (log10)") + ylab("seconds") + theme_bw()

ggsave("interactive/plots/est_ram_comp.png")


