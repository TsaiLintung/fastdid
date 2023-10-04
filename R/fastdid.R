#' Fast DID Estimation
#'
#' This function performs Difference-in-Differences (DID) estimation using fast computation techniques.
#'
#' @param dt A data table containing the panel data.
#' @param timevar The name of the time variable.
#' @param cohortvar The name of the cohort variable.
#' @param unitvar The name of the unit (group) variable.
#' @param control_option A character string indicating the control option for the DID estimation. Default is "both".
#' @param result_type A character string indicating the type of result to be returned. Default is "group_time".
#' @param boot Logical, indicating whether bootstrapping should be performed. Default is FALSE.
#' @param biters The number of bootstrap iterations. Only relevant if boot = TRUE. Default is 1000.
#' @param weightvar The name of the weight variable (optional).
#' @param clustervar The name of the cluster variable (optional).
#' @param covariatesvar A character vector containing the names of covariate variables (optional).
#' 
#' @import data.table speedglm stringr magrittr collapse
#' 
#' @return A data table containing the estimated treatment effects and standard errors.
#' @export
#'
#' @examples
#' # Example usage of fastdid function
#' result <- fastdid(data, "time", "cohort", "unit")
#'
#' @seealso
#' \code{\link{estimate_gtatt}}, \code{\link{aggregate_gt_result}}, \code{\link{get_se}}, \code{\link{convert_targets}}
#'
#' @keywords difference-in-differences fast computation panel data estimation
fastdid <- function(dt,
                    timevar, cohortvar, unitvar,
                    control_option="both",result_type="group_time", 
                    boot=FALSE, biters = 1000,
                    weightvar=NULL,clustervar=NULL,covariatesvar = c()
                    ){
  
  
  # validation arguments --------------------------------------------------------
  
  
  # preprocess -------------------------------------------------------------------
  
  #make the dt conform to innocolus assumptions of fastdid
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))
  
  setorder(dt, time, G, unit) #sort the dataset essential for the sort-once-quick-access 
  
  #time offset
  time_offset <- dt[1,time] - 1 #assume time starts at 1, first is min after sort :)
  dt[, G := G-time_offset]
  dt[, time := time-time_offset]
  
  #get info about the dataset
  time_periods <- dt[, unique(time)]
  time_size <- length(time_periods)
  id_size <- dt[, uniqueN(unit)]
  
  #the time-invariant parts (cohort, covariates, weight, cluster)
  dt_inv <- dt[1:id_size]
  
  #more info about the dataset
  cohorts <- dt_inv[, unique(G)]
  cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]
  
  #the outcomes list for fast access later
  outcomes <- list()
  for(i in time_periods){
    start <- (i-1)*id_size+1
    end <- i*id_size
    outcomes[[i]] <- dt[seq(start,end), .(y)]
  }
  
  # the optional columns
  if(length(covariatesvar)>0){
    covariates <- dt_inv[,.SD, .SDcols = covariatesvar]
  } else {covariates <- c()}
  
  if(!is.null(clustervar)){
    cluster <- dt_inv[, .SD, .SDcols = clustervar] |> unlist()
  } else {cluster <- NULL}
  
  if(!is.null(weightvar)){
    weights <- dt_inv[, .SD, .SDcols = weightvar] |> unlist()
  } else {weights <- rep(1,id_size)}
  
  # validate data -----------------------------------------------------
  
  if(any(sapply(covariates, sd) == 0)){stop("some covariates have no variation")}
  
  # main part  -------------------------------------------------
  
  # attgt
  gt_result <- estimate_gtatt(outcomes, covariates, weights,
                              cohort_sizes,cohorts,id_size,time_periods, #info about the dt
                              control_option)
  
  # aggregate att and inf function
  agg_result <- aggregate_gt_result(gt_result, cohort_sizes,
                                    result_type)
  
  #influence function -> se
  agg_se <- get_se(agg_result$inf_matrix, boot, biters, cluster)
  
  # post process -----------------------------------------------
  
  result <- data.table(target = agg_result$targets, att = agg_result$agg_att, se = agg_se)
  
  #convert "targets" back to meaningful parameter identifiers like cohort 1 post, time 2 post 
  result <- result |> convert_targets(result_type, time_offset, max(time_periods))
  
  
  #recover the original dt
  dt[, G := G+time_offset]
  dt[, time := time+time_offset]
  setnames(dt, c("time", "G", "unit"), c(timevar, cohortvar, unitvar))
  
  return(result)
  
}

convert_targets <- function(results, result_type, 
                            time_offset, max_time){
  
  if(result_type == "dynamic"){
    setnames(results, "target", "event_time")
    
  } else if (result_type == "cohort"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := abs(target)+time_offset]
    setnames(results, "target", "cohort")
    
  } else if (result_type == "calendar"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := abs(target)+time_offset]
    setnames(results, "target", "time")
    
  } else if (result_type == "group_time"){
    
    results[, cohort := floor((target-1)/max_time)]
    results[, time := (target-cohort*max_time)]
    
    #recover the time
    results[, cohort := cohort + time_offset]
    results[, time := time + time_offset]
    
    results[, target := NULL]
    
  }
  return(results)
}