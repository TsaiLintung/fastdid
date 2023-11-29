#' Fast Staggered DID Estimation
#'
#' Performs Difference-in-Differences (DID) estimation fast.
#'
#' @param dt A data table containing the panel data.
#' @param timevar The name of the time variable.
#' @param cohortvar The name of the cohort (group) variable.
#' @param unitvar The name of the unit (id) variable.
#' @param outcomevar The name of the outcome variable.
#' @param control_option The control units used for the DiD estimates. Default is "both".
#' @param result_type A character string indicating the type of result to be returned. Default is "group_time".
#' @param balanced_event_time A numeric scalar that indicates the max event time to balance the cohort composition, only meaningful when result_type == "dynamic". Default is NULL
#' @param boot Logical, indicating whether bootstrapping should be performed. Default is FALSE.
#' @param biters The number of bootstrap iterations. Only relevant if boot = TRUE. Default is 1000.
#' @param weightvar The name of the weight variable (optional).
#' @param clustervar The name of the cluster variable, can only be used when boot == TRUE (optional).
#' @param covariatesvar A character vector containing the names of covariate variables (optional).
#' @param copy whether to copy the dataset before processing, set to true if the original dataset is to be re-used.
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
#' result <- fastdid(dt, timevar = "time", cohortvar = "G", 
#'                   unitvar = "unit", outcomevar = "y",  
#'                   result_type = "group_time")
#' 
#' #control for covariates
#' result2 <- fastdid(dt, timevar = "time", cohortvar = "G", 
#'                    unitvar = "unit", outcomevar = "y",  
#'                    result_type = "group_time",
#'                    covariatesvar = c("x", "x2"))
#'                   
#' #bootstrap and clustering
#' result3 <- fastdid(dt, timevar = "time", cohortvar = "G", 
#'                    unitvar = "unit", outcomevar = "y",  
#'                    result_type = "group_time",
#'                    boot = TRUE, clustervar = "x")
#'
#' #estimate for multiple outcomes
#' result4 <- fastdid(dt, #the dataset
#'                    timevar = "time", cohortvar = "G", unitvar = "unit", 
#'                    outcomevar = c("y", "y2"), #name of the outcome columns
#'                    result_type = "group_time") 
#'
#' @keywords difference-in-differences fast computation panel data estimation did
fastdid <- function(dt,
                    timevar, cohortvar, unitvar, outcomevar, 
                    control_option="both",result_type="group_time", balanced_event_time = NULL,
                    boot=FALSE, biters = 1000,
                    weightvar=NULL,clustervar=NULL,covariatesvar = NULL,
                    copy = FALSE, validate = TRUE
                    ){
  
  
  # validation arguments --------------------------------------------------------
  
  if(copy){dt <- copy(dt)}
  
  if(!is.data.table(dt)) stop("rawdata must be a data.table")
  
  dt_names <- names(dt)
  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_arg(timevar, unitvar, cohortvar, "scalar charin", .choices = dt_names, .message = name_message)
  
  covariate_message <- "__ARG__ must be NULL or a character vector which are all names of columns from the dataset."
  check_arg(covariatesvar, outcomevar, 
            "NULL | multi charin", .choices = dt_names, .message = covariate_message)
  
  checkvar_message <- "__ARG__ must be NULL or a character scalar if a name of columns from the dataset."
  check_arg(weightvar, clustervar,
            "NULL | scalar charin", .choices = dt_names, .message = checkvar_message)
  
  check_arg(control_option, "scalar charin", .choices = c("both", "never", "notyet")) #kinda bad names since did's notyet include both notyet and never
  check_arg(copy, validate, "scalar logical")
  
  if(!is.null(balanced_event_time)){
    if(result_type != "dynamic"){stop("balanced_event_time is only meaningful with result_type == 'dynamic'")}
    check_arg(balanced_event_time, "numeric scalar")
  }
  
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))

  # validate data -----------------------------------------------------
  
  
  if(validate){
    varnames <- c("time", "G", "unit", outcomevar,weightvar,clustervar,covariatesvar)
    dt <- validate_did(dt, covariatesvar, varnames, balanced_event_time)
  }
  
  # preprocess -----------------------------------------------------------
  
  #make dt conform to the WLOG assumptions of fastdid
  
  #change to int before sorting
  if(!is.numeric(dt[, G])){
    dt[, G := as.numeric(G)]
  }
  if(!is.numeric(dt[, time])){
    dt[, time := as.numeric(time)] 
  }
  
  setorder(dt, time, G, unit) #sort the dataset essential for the sort-once-quick-access 
  
  #deal with time
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
  
  #the outcomes list for fast access later
  id_size <- dt[, uniqueN(unit)]
  
  outcomes_list <- list()
  for(outcol in outcomevar){
    outcomes <- list()
    for(i in time_periods){
      start <- (i-1)*id_size+1
      end <- i*id_size
      outcomes[[i]] <- dt[seq(start,end), get(outcol)]
    }
    outcomes_list[[outcol]] <- outcomes
  }

  #the time-invariant parts 
  dt_inv <- dt[1:id_size]
  cohorts <- dt_inv[, unique(G)]
  cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]

  # the optional columns
  if(!is.null(covariatesvar)){
    covariates <- dt_inv[,.SD, .SDcols = covariatesvar]
    ipw_formula <- paste0(covariatesvar, collapse = "+")
  } else {
    covariates <- NULL
    ipw_formula <- NULL
  }
  
  if(!is.null(clustervar)){
    cluster <- dt_inv[, .SD, .SDcols = clustervar] |> unlist()
  } else {cluster <- NULL}
  
  if(!is.null(weightvar)){
    weights <- dt_inv[, .SD, .SDcols = weightvar] |> unlist()
  } else {weights <- rep(1,id_size)}
  
  # main part  -------------------------------------------------
  
  # attgt
  gt_result_list <- estimate_gtatt(outcomes_list, outcomevar, covariates, ipw_formula, weights,
                                   cohort_sizes,cohorts,id_size,time_periods, #info about the dt
                                   control_option)
  
  all_result <- data.table()
  for(outcol in outcomevar){
    gt_result <- gt_result_list[[outcol]]
    
    # aggregate att and inf function
    agg_result <- aggregate_gt(gt_result, cohort_sizes, 
                               weights, dt_inv[, G],
                               result_type, balanced_event_time)
    
    #get se from the influence function
    agg_se <- get_se(agg_result$inf_matrix, boot, biters, cluster)
    
    # post process -----------------------------------------------
    
    result <- data.table(agg_result$targets, agg_result$agg_att, agg_se)
    names(result) <- c("target", "att", "se")
    
    #convert "targets" back to meaningful parameter identifiers like cohort 1 post, time 2 post 
    result <- result |> convert_targets(result_type, time_offset, time_step, max(time_periods))
    result[, outcome := outcol]
    all_result <- rbind(all_result, result)
    
    rm(result)
    
  }

  return(all_result)
  
}

convert_targets <- function(results, result_type, 
                            time_offset, time_step, max_time){
  
  if(result_type == "dynamic"){
    setnames(results, "target", "event_time")
    
  } else if (result_type == "cohort"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), time_offset, time_step)]
    setnames(results, "target", "cohort")
    
  } else if (result_type == "calendar"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), time_offset, time_step)]
    setnames(results, "target", "time")
    
  } else if (result_type == "group_time"){
    
    results[, cohort := floor((target-1)/max_time)]
    results[, time := (target-cohort*max_time)]
    
    #recover the time
    results[, cohort := recover_time(cohort, time_offset, time_step)]
    results[, time := recover_time(time, time_offset, time_step)]
    
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