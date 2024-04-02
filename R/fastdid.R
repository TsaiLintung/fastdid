#' Fast Staggered DID Estimation
#'
#' Performs Difference-in-Differences (DID) estimation fast.
#'
#' @param data A data table containing the panel data.
#' @param timevar The name of the time variable.
#' @param cohortvar The name of the cohort (group) variable.
#' @param unitvar The name of the unit (id) variable.
#' @param outcomevar The name of the outcome variable.
#' @param control_option The control units used for the DiD estimates. Default is "both".
#' @param result_type A character string indicating the type of result to be returned. Default is "group_time".
#' @param balanced_event_time A numeric scalar that indicates the max event time to balance the cohort composition, only meaningful when result_type == "dynamic". Default is NA
#' @param control_type The method for controlling for covariates. "ipw" for inverse probability weighting, "reg" for outcome regression, or "dr" for doubly-robust
#' @param allow_unbalance_panel Whether allow unbalance panel as input (if false will coerce the dataset to a balanced panel). Default is FALSE 
#' @param boot Logical, indicating whether bootstrapping should be performed. Default is FALSE.
#' @param biters The number of bootstrap iterations. Only relevant if boot = TRUE. Default is 1000.
#' @param weightvar The name of the weight variable (optional).
#' @param clustervar The name of the cluster variable, can only be used when boot == TRUE (optional).
#' @param covariatesvar A character vector containing the names of time-invariant covariate variables (optional).
#' @param varycovariatesvar A character vector containing the names of time-varying covariate variables (optional).
#' @param copy whether to copy the dataset before processing, set to false to speed up the process, but the input data will be altered.
#' @param validate whether to validate the dataset before processing.
#' 
#' @import data.table parglm stringr collapse dreamerr BMisc 
#' @importFrom stats quantile vcov sd binomial fitted qnorm rnorm as.formula
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
                    control_type = "ipw", allow_unbalance_panel = FALSE, boot=FALSE, biters = 1000,
                    weightvar=NA,clustervar=NA,covariatesvar = NA, varycovariatesvar = NA,
                    copy = TRUE, validate = TRUE,
                    max_control_cohort_diff = NA
                    ){
  
  # validate arguments --------------------------------------------------------

  if(!is.data.table(data)){
    warning("coercing input into a data.table.")
    data <- as.data.table(data)
  } 
  if(copy){dt <- copy(data)} else {dt <- data}
  
  dt_names <- names(dt)
  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_set_arg(timevar, unitvar, cohortvar, "match", .choices = dt_names, .message = name_message)
  
  covariate_message <- "__ARG__ must be NA or a character vector which are all names of columns from the dataset."
  check_set_arg(varycovariatesvar, covariatesvar, outcomevar, 
            "NA | multi match", .choices = dt_names, .message = covariate_message)
  
  checkvar_message <- "__ARG__ must be NA or a character scalar if a name of columns from the dataset."
  check_set_arg(weightvar, clustervar,
            "NA | match", .choices = dt_names, .message = checkvar_message)
  
  check_set_arg(control_option, "match", .choices = c("both", "never", "notyet")) #kinda bad names since did's notyet include both notyet and never
  check_set_arg(control_type, "match", .choices = c("ipw", "reg", "dr")) 
  check_arg(copy, validate, boot, allow_unbalance_panel, "scalar logical")
  
  if(!is.na(balanced_event_time)){
    if(result_type != "dynamic"){stop("balanced_event_time is only meaningful with result_type == 'dynamic'")}
    check_arg(balanced_event_time, "numeric scalar")
  }
  if(allow_unbalance_panel == TRUE & control_type %in% c("dr", "reg")){
    stop("fastdid currently only supprts ipw when allowing for unbalanced panels.")
  }
  if(allow_unbalance_panel == TRUE & !allNA(varycovariatesvar)){
    stop("fastdid currently only supprts time varying covariates when allowing for unbalanced panels.")
  }
  if(any(covariatesvar %in% varycovariatesvar) & !allNA(varycovariatesvar) & !allNA(covariatesvar)){
    stop("time-varying var and invariant var have overlaps.")
  }
  if(!boot & !allNA(clustervar)){
    stop("clustering only available with bootstrap")
  }
  
  p <- list(timevar = timevar,
            cohortvar = cohortvar,
            unitvar = unitvar,
            outcomevar =  outcomevar,
            weightvar = weightvar,
            clustervar = clustervar,
            covariatesvar = covariatesvar,
            varycovariatesvar = varycovariatesvar,
            control_option = control_option,
            result_type = result_type,
            balanced_event_time = balanced_event_time,
            control_type = control_type,
            allow_unbalance_panel = allow_unbalance_panel,
            boot = boot, 
            biters = biters,
            max_control_cohort_diff = max_control_cohort_diff)
  
  
  # validate data -----------------------------------------------------
  
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))
  
  if(validate){
    varnames <- c("time", "G", "unit", outcomevar,weightvar,clustervar,covariatesvar,varycovariatesvar)
    dt <- validate_did(dt, varnames, p)
  }
  
  # preprocess -----------------------------------------------------------
  
  #make dt conform to the WLOG assumptions of fastdid
  coerce_result <- coerce_dt(dt, p) #also changed dt
  dt <- coerce_result$dt
  # get auxiliary data
  aux <- get_auxdata(dt, p)
  
  # main part  -------------------------------------------------

  # estimate attgt
  gt_result_list <- estimate_gtatt(aux, p)
  
  # aggregate the result by outcome
  all_result <- rbindlist(lapply(gt_result_list, function(x){
    result <- aggregate_gt(x, aux, p)
    #convert "targets" back to meaningful parameter identifiers like cohort 1 post, time 2 post 
    result <- result |> convert_targets(result_type, coerce_result$time_change)
    return(result)
  }))

  return(all_result)
  
}

# small steps ----------------------------------------------------------------------

coerce_dt <- function(dt, p){
  
  #change to int before sorting
  if(!is.numeric(dt[, G])){
    dt[, G := as.numeric(G)]
  }
  if(!is.numeric(dt[, time])){
    dt[, time := as.numeric(time)] 
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
  
  return(list(dt = dt,
              time_change = list(time_step = time_step,
                            max_time = max(time_periods),
                            time_offset = time_offset)))
  
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
  
  if(!allNA(p$covariatesvar)){
    covariates <- cbind(const = -1, dt_inv[,.SD, .SDcols = p$covariatesvar])
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
              weights = weights)
  
  return(aux)
  
}

convert_targets <- function(results, result_type, t){
  
  if(result_type == "dynamic"){
    setnames(results, "target", "event_time")
    
  } else if (result_type == "cohort"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), t$time_offset, t$time_step)]
    setnames(results, "target", "cohort")
    
  } else if (result_type == "calendar"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), t$time_offset, t$time_step)]
    setnames(results, "target", "time")
    
  } else if (result_type == "group_time"){
    
    results[, cohort := floor((target-1)/t$max_time)]
    results[, time := (target-cohort*t$max_time)]
    
    #recover the time
    results[, cohort := recover_time(cohort, t$time_offset, t$time_step)]
    results[, time := recover_time(time, t$time_offset, t$time_step)]
    
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