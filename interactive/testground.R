
rm(list = ls())
gc()

library(profvis)
library(devtools)
library(peakRAM)
library(microbenchmark)
library(did)
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

source("R/sim_did.R")
source("interactive/source_raw/drdidipw.R")


# load event code ---------------------------------------------------------------------

# simple ---------------------------------------------------------------------

simdt <- sim_did(1000, 10, cov = "no", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, stratify = FALSE)
dt <- simdt$dt

started.at <- proc.time()

covariatesvar <- c()
boot <- FALSE

# preprocess ----------------------------------------------------------

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
  #if(!start_cohort %in% cohort_sizes[,unique(G)]|!end_cohort %in% cohort_sizes[,unique(G)]) {stop("cohort not in cohort_sizes")}
  start <- cohort_sizes[G < start_cohort, sum(cohort_size)]+1
  end <- cohort_sizes[G <= end_cohort, sum(cohort_size)]
  return(start:end)
}

outcomes <- list()
for(i in 1:time_size){
  outcomes[[i]] <- dt[seq((i-1)*id_size+1,i*id_size), .(y)]
}

# attgt  -------------------------------------------------------------------------------------------------------------- 

control_option <- "both"
weights <- rep(1,id_size)

if(!Inf %in% cohorts & control_option != "notyet"){
  warning("no never-treated availble, switching to not-yet-treated control")
  control_option <- "notyet"
}
max_control_cohort <- ifelse(control_option == "notyet", max(treated_cohorts), Inf) 

last_coef <- NULL
gt_att <- data.table()
gt_inf_func <- data.table(unit = dt_inv[, unit])
for(g in cohorts){
  for(t in time_periods){
    
    #setup and checks
    base_period <- g-1
    min_control_cohort <- ifelse(control_option == "never", Inf, max(t+1, base_period+1)) #not-yet treated / never treated in both base and "treated" period
    if(t == base_period){next} #no treatmenteffect for the base period
    if(g >= max_control_cohort){next} #no treatmenteffect for never treated or the last treated cohort (notyet)
    if(t >= max_control_cohort){next} #no control available if the last cohort is treated too
    
    #select the right cohorts
    did_setup <- rep(NA, id_size)
    did_setup[get_cohort_pos(cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
    did_setup[get_cohort_pos(cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
    
    #select the post and pre outcome
    pre_outcome <- outcomes[[base_period]]
    post_outcome <- outcomes[[t]] 
    
    #construct the 2x2 dataset
    cohort_did <- data.table(D = did_setup, post = post_outcome, pre = pre_outcome, cov = covariates)
    cohort_did[, weights := weights] #the default weight, should allow overwrite next
    
    #estimate did
    results <- estimate_did(cohort_did, last_coef)
    last_coef <- results$logit_coef
    
    #collect the result
    gt_att <- rbind(gt_att, data.table(G = g, time = t, att = results$att))
    gt_inf_func[[paste0(g, ".", t)]] <- results$inf_func
    
  }
}
gt_inf_func[, unit := NULL]

# construct weight -------------------------------------------------------------------------------------

weights <- data.table()
group_time <- gt_att[,.(G, time)] |> merge(cohort_sizes, by = "G")

gt_count <- group_time[, .N]

result_type <- "group_time"
balanced_composition <- FALSE

bool_to_pn <- function(x){ifelse(x, 1, -1)}

if (result_type == "dynamic") {
  group_time[, target := time-G]
} else if (result_type == "group") {
  group_time[, target := G*(bool_to_pn(time>=G))]
} else if (result_type == "time") {
  group_time[, target := time*(bool_to_pn(time>=G))]
} else if (result_type == "group_time"){
  group_time[, target := G*max(time)+time]
}

min_target<- group_time[, min(target)]
max_target<- group_time[, max(target)]

targets <- c()
for(tar in min_target:max_target){
  group_time[, agg_weight := 0]
  total_size <- group_time[target == tar, sum(cohort_size)]
  group_time[, weight := ifelse(target == tar, cohort_size/total_size, 0)]
  target_weights <- group_time[, .(weight)] |> transpose()
  names(target_weights) <- names(gt_inf_func)
  
  targets <- c(targets, tar)
  weights <- rbind(weights, target_weights)
}

# aggregate  -------------------------------------------------------------------------------------------

agg_att <- as.matrix(weights) %*% as.vector(gt_att[, att]) |> as.vector()
inf_matrix <- as.matrix(gt_inf_func) %*% t(as.matrix(weights))

# get bootstrap se ------------------------------------------------------------------------------------------

clustervar <- "unit"

if(boot){
  top_quant <- 0.75
  bot_quant <- 0.25
  boot_results <- BMisc::multiplier_bootstrap(as.matrix(inf_matrix), biters = 1000) %>% as.data.table()
  boot_top <- boot_results[, lapply(.SD, function(x) quantile(x, top_quant, type=1, na.rm = TRUE)),]
  boot_bot <- boot_results[, lapply(.SD, function(x) quantile(x, bot_quant, type=1, na.rm = TRUE)),]
  
  dt_se <- rbind(boot_bot, boot_top) %>% transpose()
  names(dt_se) <- c("boot_bot", "boot_top")
  dt_se <- dt_se[,(boot_top-boot_bot)/(qnorm(top_quant) - qnorm(bot_quant))]
} else {
  inf_matrix <- inf_matrix %>% as.data.table()
  dt_se <- inf_matrix[, lapply(.SD, function(x) sd(x)/sqrt(length(x)))] %>% as.vector()
}

# gather results -------------------------------------------------------------------------------------------------

did_result <- data.table(target = targets, att = agg_att, se = dt_se)


#compare with did -------------------------------------------------------------------------------------------------

#construct the 2x2 dataset
# 
# did_result <- att_gt(yname = "y",
#                      gname = "G",
#                      idname = "unit",
#                      tname = "time",
#                      data = dt,
#                      #xformla = ~x,
#                      base_period = "universal",
#                      control_group = "notyettreated",
#                      est_method = "ipw",
#                      clustervars = clustervar)
# 
# did_result_dt <- data.table(G = did_result$group, time = did_result$t, did_att = did_result$att, did_se = did_result$se)
# 
# compare <- did_result_dt |> merge(dt_se, by = c("G", "time"), all = TRUE) |> merge(gt_att, by = c("G", "time"), all = TRUE)
# 
# compare[, se_diff := se-did_se]
# compare[, att_diff := att-did_att]
# compare


timetaken(started.at)
did_result