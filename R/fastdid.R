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
#' @param full return the full result, like the influence function, call, etc,. Default is false. 
#' @param parallel whether to use parallization (only available on unix systesm like Mac or Linux.)
#' @param cohortvar2 The name of the second cohort (group) variable.
#' @param event_specific Whether to recover event_specific treatment effect or report combined effect. 
#' 
#' @import data.table stringr dreamerr ggplot2
#' @importFrom stats quantile vcov sd binomial fitted qnorm rnorm as.formula weighted.mean
#' @importFrom collapse allNA fnrow whichNA fnunique fsum na_insert
#' @importFrom parallel mclapply
#' @importFrom BMisc multiplier_bootstrap
#' @importFrom parglm parglm.fit parglm.control
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
                    exper = NULL, full = FALSE, parallel = FALSE, 
                    cohortvar2 = NA, event_specific = TRUE){
  
  # preprocess --------------------------------------------------------
  
  if(!is.data.table(data)){
    warning("coercing input into a data.table.")
    data <- as.data.table(data)
  } 
  if(copy){dt <- copy(data)} else {dt <- data}
  
  # validate arguments
  p <- as.list(environment()) #collect everything besides data
  p$data <- NULL
  p$dt <- NULL
  
  exper_args <- c("filtervar", "filtervar_post", "only_balance_2by2",
                  "aggregate_scheme", "max_control_cohort_diff")
  p$exper <- get_exper_default(p$exper, exper_args)
  class(p) <- "locked" #no more changes!
  validate_argument(dt, p)
  
  #change name for main columns
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))
  if(!is.na(p$cohortvar2)){setnames(dt, p$cohortvar2, "G2")}
  
  # validate and throw away not legal data 
  dt <- validate_dt(dt, p)
  
  #make dt conform to the WLOG assumptions of fastdid
  coerce_result <- coerce_dt(dt, p) #also changed dt
  
  # get auxiliary data
  aux <- get_auxdata(coerce_result$dt, p)

  # main estimation  -------------------------------------------------

  gt_result_list <- estimate_gtatt(aux, p)
  agg_result <- aggregate_gt(gt_result_list, aux, p)
  
  # post process -------------------------------------------
  
  #convert "targets" back to meaningful parameters
  est_results <- convert_targets(agg_result$est, p, coerce_result$t) 
  
  if(!p$full){
    return(est_results)
  } else {
    full_result <- list(
      call = p,
      estimate = est_results,
      gt_estimate = gt_result_list,
      agg_inf_func = agg_result$inf_func,
      agg_weight_matrix = agg_result$agg_weight_matrix
    )
    class(full_result) <- c("fastdid_result", class(full_result))
    return(full_result)
  }
}