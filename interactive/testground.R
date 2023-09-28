
rm(list = ls())
gc()

library(profvis)
library(devtools)
library(peakRAM)
library(microbenchmark)
library(did)

library(fastglm)
library(data.table)
library(collapse)
library(speedglm)
library(RcppArmadillo)

#set the wd to the source folder
setwd("~/GitHub/EventStudyCode")

load_all()

source("R/sim_did.R")
source("interactive/source_raw/drdidipw.R")


# load event code ---------------------------------------------------------------------

# simple ---------------------------------------------------------------------

simdt <- sim_did(10000, 10, cov = "no", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, stratify = FALSE)
dt <- simdt$dt

started.at <- proc.time()

covariatesvar <- c()

# preprocess

setorder(dt, time, G, unit)

id_size <- dt[, uniqueN(unit)]
time_periods <- dt[, unique(time)]
time_size <- length(time_periods)
cohorts <- dt[, unique(G)]
treated_cohorts <- cohorts[!is.infinite(cohorts)]
dt_inv <- dt[1:id_size]
cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]

if(length(covariatesvar)>0){
  covariates <- dt_inv[,.SD, .SDcols = covariatesvar] 
} else {covariates <- c()}

get_cohort_pos <- function(cohort_sizes, start_cohort, end_cohort = start_cohort){
  if(!start_cohort %in% cohort_sizes[,unique(G)]|!end_cohort %in% cohort_sizes[,unique(G)]) {stop("cohort not in cohort_sizes")}
  start <- cohort_sizes[G < start_cohort, sum(cohort_size)]+1
  end <- cohort_sizes[G <= end_cohort, sum(cohort_size)]
  return(start:end)
}

outcomes <- list()
for(i in 1:time_size){
  outcomes[[i]] <- dt[seq((i-1)*id_size+1,i*id_size), .(y)]
}

# attgt  -------------------------------------------------------------------------------------------------------------- 

control_option <- "notyet"

if(!Inf %in% cohorts & control_option != "notyet"){
  warning("no never-treated availble, switching to not-yet-treated control")
  control_option <- "notyet"
}

max_control_cohort <- switch (control_option,
  "notyet" = max(treated_cohorts),
  "both" = Inf,
  "never" = Inf
)

last_coef <- NULL
gt_att <- data.table()
gt_inf_func <- data.table(unit = dt_inv[, unit])
for(g in cohorts){
  for(t in time_periods){
    message(g,t)
    base_period <- g-1
    min_control_cohort <- ifelse(control_option == "never", Inf, t+1) #not-yet treated / never treated
    if(t == base_period){next} #no treatmenteffect for the base period
    if(g == Inf){next} #no treatmenteffect for never treated
    if(t >= max_control_cohort){next} #no control available if no never treated at the end
    
    did_setup <- rep(NA, id_size)
    did_setup[get_cohort_pos(cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
    did_setup[get_cohort_pos(cohort_sizes, g)] <- 1 #treated cannot be controls, must assign treated after control to overwrite
    
    pre_outcome <- outcomes[[base_period]]
    post_outcome <- outcomes[[t]] 
    
    cohort_did <- data.table(D = did_setup, post = post_outcome, pre = pre_outcome, cov = covariates)
    cohort_did[, weights := 1]
    
    
    results <- estimate_did(cohort_did, last_coef)
    last_coef <- results$logit_coef
    
    gt_att <- rbind(gt_att, data.table(G = g, time = t, att = results$att))
    gt_inf_func[[paste0(g, ".", t)]] <- results$inf_func
    
  }
}

timetaken(started.at)

# aggregate  -------------------------------------------------------------------------------------------

names(gt_inf_func)

gt_att <- gt_att |> merge(cohort_sizes, by = "G")
gt_att[, event_time := time-G]
event_est <- gt_att[, .(att = sum(att*cohort_size)/sum(cohort_size)), by = "event_time"]


event_est <- data.table()
for(event_time in (min(time_periods)-max(treated_cohorts)):(max(time_periods)-min(treated_cohorts))){
  message(event_time)
 
}

