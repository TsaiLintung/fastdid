#' Fast Staggered DID Estimation
#'
#' Performs Difference-in-Differences (DID) estimation fast.
#'
#' @param data A data.table containing the panel data.
#' @param timevar The name of the time variable.
#' @param cohortvar The name of the cohort (group) variable.
#' @param unitvar The name of the unit (id) variable.
#' @param outcomevar The name of the outcome variable.
#' @param control_option The control units used for the DiD estimates. Default is "both".
#' @param result_type A character string indicating the type of result to be returned. Default is "group_time".
#' @param balanced_event_time A numeric scalar that indicates the max event time to balance the cohort composition, only meaningful when result_type == "dynamic". Default is NA
#' @param control_type The method for controlling for covariates. "ipw" for inverse probability weighting, "reg" for outcome regression, or "dr" for doubly-robust
#' @param allow_unbalance_panel Whether allow unbalance panel as input (if false will coerce the dataset to a balanced panel). Default is FALSE 
#' @param boot Logical, indicating whether bootstrapping should be performed. Default is FALSE
#' @param biters The number of bootstrap iterations. Only relevant if boot = TRUE. Default is 1000.
#' @param cband Logical, indicate whether to use uniform confidence band or point-wise, defulat is FALSE (use point-wise)
#' @param alpha The significance level, default is 0.05
#' @param weightvar The name of the weight variable, if not specified will cluster on unit level (optional).
#' @param clustervar The name of the cluster variable, can only be used when boot == TRUE (optional).
#' @param covariatesvar A character vector containing the names of time-invariant covariate variables (optional).
#' @param varycovariatesvar A character vector containing the names of time-varying covariate variables (optional).
#' @param copy whether to copy the dataset before processing, set to false to speed up the process, but the input data will be altered.
#' @param validate whether to validate the dataset before processing.
#' @param anticipation periods with aniticipation (delta in CS, default is 0, reference period is g - delta - 1).
#' @param exper the list of experimental features, for features that are not in CSDID originally. Generally less tested. 
#' @param base_period same as did
#' 
#' @import data.table parglm stringr dreamerr BMisc 
#' @importFrom stats quantile vcov sd binomial fitted qnorm rnorm as.formula
#' @importFrom collapse allNA fnrow whichNA fnunique fsum na_insert
#' @return A data.table containing the estimated treatment effects and standard errors.
#' @export
#'
#' @examples
#' # simulated data
#' simdt <- sim_did(1e+03, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE)
#' dt <- simdt$dt
#' 
#' #basic call
#' result <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                   unitvar = "unit", outcomevar = "y",  
#'                   result_type = "group_time")
#' 
#' #control for covariates
#' result2 <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                    unitvar = "unit", outcomevar = "y",  
#'                    result_type = "group_time",
#'                    covariatesvar = c("x", "x2"))
#'                   
#' #bootstrap and clustering
#' result3 <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                    unitvar = "unit", outcomevar = "y",  
#'                    result_type = "group_time",
#'                    boot = TRUE, clustervar = "x")
#'
#' #estimate for multiple outcomes
#' result4 <- fastdid(data = dt, #the dataset
#'                    timevar = "time", cohortvar = "G", unitvar = "unit", 
#'                    outcomevar = c("y", "y2"), #name of the outcome columns
#'                    result_type = "group_time") 
#'
#' @keywords difference-in-differences fast computation panel data estimation did
fastdid <- function(data,
                    timevar, cohortvar, unitvar, outcomevar, 
                    control_option="both",result_type="group_time", balanced_event_time = NA,
                    control_type = "ipw", allow_unbalance_panel = FALSE, boot=FALSE, biters = 1000, cband = FALSE, alpha = 0.05,
                    weightvar=NA,clustervar=NA, covariatesvar = NA, varycovariatesvar = NA, 
                    copy = TRUE, validate = TRUE,
                    anticipation = 0,  base_period = "universal",
                    exper = NULL){

  # validation --------------------------------------------------------
  
  if(!is.data.table(data)){
    warning("coercing input into a data.table.")
    data <- as.data.table(data)
  } 
  if(copy){dt <- copy(data)} else {dt <- data}
  
  # validate arguments
  p <- as.list(environment()) #collect everything besides data
  p$data <- NULL
  p$dt <- NULL
  p$exper <- get_exper_default(p$exper)
  validate_argument(dt, p)

  # validate and throw away not legal data 
  
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))
  dt <- validate_dt(dt, p)
  
  # preprocess -----------------------------------------------------------
  
  #make dt conform to the WLOG assumptions of fastdid
  coerce_result <- coerce_dt(dt, p) #also changed dt
  dt <- coerce_result$dt
  p <- coerce_result$p
  
  # get auxiliary data
  aux <- get_auxdata(dt, p)
  
  # main part  -------------------------------------------------

  gt_result_list <- estimate_gtatt(aux, p)

  all_result <- aggregate_gt(gt_result_list, aux, p)
  
  #convert "targets" back to meaningful parameter identifiers like cohort 1 post, time 2 post 
  all_result <- convert_targets(all_result, p) 
  
  return(all_result)
  
}

# small steps ----------------------------------------------------------------------

get_exper_default <- function(exper){
  exper_args <- c("filtervar", "min_dynamic", "max_dynamic", "min_control_cohort_diff", "max_control_cohort_diff")
  for(arg in exper_args){
    if(is.null(exper[[arg]])){
      exper[[arg]] <- NA
    }
  }
  return(exper)
}

coerce_dt <- function(dt, p){


  #chcek if there is availble never-treated group
  if(!is.infinite(dt[, max(G)]) & p$control_option != "notyet"){
    warning("no never-treated availble, switching to not-yet-treated control")
    p$control_option <- "notyet"
  }
  
  if(p$allow_unbalance_panel){
    dt_inv_raw <- dt[dt[, .I[1], by = unit]$V1]
    setorder(dt_inv_raw, G)
    dt_inv_raw[, new_unit := 1:.N] #let unit start from 1 .... N, useful for knowing which unit is missing
    dt <- dt |> merge(dt_inv_raw[,.(unit, new_unit)], by = "unit")
    dt[, unit := new_unit]
  }
  
  setorder(dt, time, G, unit) #sort the dataset essential for the sort-once-quick-access 
  #deal with time, coerice time to 1,2,3,4,5.......
  time_periods <- dt[, unique(time)]
  time_size <- length(time_periods)
  
  time_offset <- min(time_periods) - 1 #assume time starts at 1, first is min after sort :)
  if(time_offset != 0){
    dt[, G := G-time_offset]
    dt[, time := time-time_offset]
    time_periods <- time_periods - time_offset
  }
  
  time_step <- 1 #time may not jump at 1
  if(any(time_periods[2:length(time_periods)] - time_periods[1:length(time_periods)-1] != 1)){
    time_step <- time_periods[2]-time_periods[1]
    time_periods <- (time_periods-1)/time_step+1
    if(any(time_periods[2:length(time_periods)] - time_periods[1:length(time_periods)-1] != 1)){stop("time step is not uniform")}
    dt[G != 1, G := (G-1)/time_step+1]
    dt[time != 1, time := (time-1)/time_step+1]
  }
  
  #add the information to p
  p$time_step = time_step
  p$time_offset = time_offset
  
  if(nrow(dt) == 0){
    stop("no data after coercing the dataset")
  }
  
  return(list(dt = dt, p = p))
  
}

get_auxdata <- function(dt, p){

  time_periods <- dt[, unique(time)]
  id_size <- dt[, uniqueN(unit)]
  
  setorder(dt, time, G, unit) #sort the dataset essential for the sort-once-quick-access 
  
  #construct the outcomes list for fast access later
  #loop for multiple outcome
  outcomes_list <- list()
  for(outcol in p$outcomevar){
    outcomes <- list()
    
    if(!p$allow_unbalance_panel){
      for(i in time_periods){
        start <- (i-1)*id_size+1
        end <- i*id_size
        outcomes[[i]] <- dt[seq(start,end), get(outcol)]
      }
    } else {
      
      for(i in time_periods){
        #populate a outcome vector of length N with outcome data in the right place
        #NA is data gone or missing, will be addressed in estimate_did_rc
        outcome_period <- rep(NA, id_size)
        data_pos <- dt[time == i, unit] #units observed in i
        outcome_period[data_pos] <- dt[time == i, get(outcol)]
        outcomes[[i]] <- outcome_period
      }
      
    }
    
    outcomes_list[[outcol]] <- outcomes
  }
  
  #the time-invariant parts 
  if(!p$allow_unbalance_panel){
    dt_inv <- dt[1:id_size]
  } else {
    dt_inv <- dt[dt[, .I[1], by = unit]$V1] #the first observation
    setorder(dt_inv, unit) #can't move this outside
  }
  
  cohorts <- dt_inv[, unique(G)]
  cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]
  
  # the optional columns
  varycovariates <- list()
  if(!allNA(p$varycovariatesvar)){
    for(i in time_periods){
      start <- (i-1)*id_size+1
      end <- i*id_size
      varycovariates[[i]] <- dt[seq(start,end), .SD, .SDcols = p$varycovariatesvar]
    }
  } else {
    varycovariates <- NA
  }
  
  # filters
  filters <- list()
  if(!is.na(p$exper$filtervar)){
    for(i in time_periods){
      start <- (i-1)*id_size+1
      end <- i*id_size
      filters[[i]] <- dt[seq(start,end), .SD, .SDcols = p$exper$filtervar]
    }
  } else {
    filters <- NA
  }
  
  if(!allNA(p$covariatesvar)){
    covariates <- dt_inv[,.SD, .SDcols = p$covariatesvar]
  } else {
    covariates <- NA
  }
  
  if(!is.na(p$clustervar)){
    cluster <- dt_inv[, .SD, .SDcols = p$clustervar] |> unlist()
  } else {
    cluster <- NA
  }
  
  if(!is.na(p$weightvar)){
    weights <- dt_inv[, .SD, .SDcols = p$weightvar] |> unlist()
  } else {
    weights <- rep(1, id_size)
  }
  
  aux <- list(time_periods = time_periods,
              id_size = id_size,
              outcomes_list = outcomes_list,
              dt_inv= dt_inv,
              cohorts = cohorts,
              cohort_sizes = cohort_sizes, 
              varycovariates = varycovariates,
              covariates = covariates,
              cluster = cluster,
              weights = weights,
              filters = filters)
  
  return(aux)
  
}

convert_targets <- function(results, p){
  
  result_type <- p$result_type

  if(result_type == "dynamic"){
    setnames(results, "target", "event_time")
    
  } else if (result_type == "cohort"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), p$time_offset, p$time_step)]
    setnames(results, "target", "cohort")
    
  } else if (result_type == "calendar"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), p$time_offset, p$time_step)]
    setnames(results, "target", "time")
    
  } else if (result_type == "group_time"){
    
    results[, cohort := as.numeric(str_split_i(target, "\\.", 1))]
    results[, time :=  as.numeric(str_split_i(target, "\\.", 2))]
    
    #recover the time
    results[, cohort := recover_time(cohort, p$time_offset, p$time_step)]
    results[, time := recover_time(time, p$time_offset, p$time_step)]
    
    results[, target := NULL]
    
  } else if (result_type == "simple") {
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := NULL]
  } 
  return(results)
}

recover_time <- function(time, time_offset, time_step){
  return(((time-1)*time_step)+1+time_offset)
}