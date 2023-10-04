
rm(list = ls())
gc()

library(profvis)
library(devtools)
library(peakRAM)
library(microbenchmark)
library(DRDID)

library(fastglm)
library(data.table)
library(collapse)
library(speedglm)
library(RcppArmadillo)
library(BMisc)


#set the wd to the source folder
setwd("~/GitHub/EventStudyCode")

load_all()
#load_all("~/GitHub/did")
library(did)


# load event code ---------------------------------------------------------------------

# simple ---------------------------------------------------------------------




boot <- TRUE
clust <- "strat"

ratio <- c()
for(i in 1:5){
  ep <- i*0.1
  simdt <- sim_did(1e+02, 10, cov = "int", hetero = "dynamic", balanced = TRUE, second_outcome = FALSE, seed = i, stratify = TRUE,
                   epsilon_size = ep)
  dt <- simdt$dt
  
  dt[, y := y*100]
  dt[, strat := s*G*x]
  
  result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", result_type = "group_time", boot = boot,
                   clustervar = clust
  )
  
  
  did_result <- att_gt(yname = "y",
                       gname = "G",
                       idname = "unit",
                       tname = "time",
                       data = dt,
                       #xformla = ~x,
                       base_period = "universal",
                       control_group = "notyettreated",
                       est_method = "ipw",
                       clustervars = clust,
                       bstrap = boot,
                       biters = 10000,
                       cband = FALSE)
  
  did_result_dt <- data.table(cohort = did_result$group, time = did_result$t, did_att = did_result$att, did_se = did_result$se)
  compare <- did_result_dt |> merge(result, by = c("cohort", "time"), all = TRUE) 
  ratio <- c(ratio, compare[,sum(did_se, na.rm = TRUE)/sum(se, na.rm = TRUE)])
}

mean(ratio)

