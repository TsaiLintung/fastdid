#' Fast Staggered DID Estimation
#'
#' Performs Difference-in-Differences (DID) estimation.
#'
#' @param data data.table, the dataset.
#' @param timevar character, name of the time variable.
#' @param cohortvar character, name of the cohort (group) variable.
#' @param unitvar character, name of the unit (id) variable.
#' @param outcomevar character vector, name(s) of the outcome variable(s).
#' @param control_option character, control units used for the DiD estimates, options are "both", "never", or "notyet". 
#' @param result_type character, type of result to return, options are "group_time", "time", "group", "simple", "dynamic" (time since event), "group_group_time", or "dynamic_stagger". 
#' @param balanced_event_time number, max event time to balance the cohort composition.
#' @param control_type character, estimator for controlling for covariates, options are "ipw" (inverse probability weighting), "reg" (outcome regression), or "dr" (doubly-robust).
#' @param allow_unbalance_panel logical, allow unbalance panel as input or coerce dataset into one. 
#' @param boot logical, whether to use bootstrap standard error. 
#' @param biters number, bootstrap iterations. Default is 1000.
#' @param cband logical, whether to use uniform confidence band or point-wise.
#' @param alpha number, the significance level. Default is 0.05.
#' @param weightvar character, name of the weight variable.
#' @param clustervar character, name of the cluster variable. 
#' @param covariatesvar character vector, names of time-invariant covariate variables. 
#' @param varycovariatesvar character vector, names of time-varying covariate variables.
#' @param copy logical, whether to copy the dataset. 
#' @param validate logical, whether to validate the dataset. 
#' @param anticipation number, periods with anticipation.
#' @param anticipation2 number, periods with anticipation for the second event.
#' @param exper list, arguments for experimental features.
#' @param base_period character, type of base period in pre-preiods, options are "universal", or "varying".
#' @param full logical, whether to return the full result (influence function, call, weighting scheme, etc,.). 
#' @param parallel logical, whether to use parallization on unix system. 
#' @param cohortvar2 character, name of the second cohort (group) variable.
#' @param event_specific logical, whether to recover target treatment effect or use combined effect.
#' @param double_control_option character, control units used for the double DiD, options are "both", "never", or "notyet". 
#' 
#' @import data.table stringr dreamerr ggplot2
#' @importFrom stats quantile vcov sd binomial fitted qnorm rnorm as.formula weighted.mean
#' @importFrom collapse allNA fnrow whichNA fnunique fsum na_insert
#' @importFrom parallel mclapply
#' @importFrom BMisc multiplier_bootstrap
#' @importFrom parglm parglm.fit parglm.control
#' @return A data.table containing the estimated treatment effects and standard errors or a list of all results when `full == TRUE`.
#' @export
#'
#' @details
#' `balanced_event_time` is only meaningful when `result_type == "dynamic`.
#' 
#' `result_type` as `group-group-time` and `dynamic staggered` is only meaningful when using double did.
#' 
#' `biter` and `clustervar` is only used when `boot == TRUE`.
#'
#' @examples
#' # simulated data
#' simdt <- sim_did(1e+02, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE, seed = 1)
#' dt <- simdt$dt
#' 
#' #basic call
#' result <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                   unitvar = "unit", outcomevar = "y",  
#'                   result_type = "group_time")
#'
#' @keywords difference-in-differences fast computation panel data estimation did
fastdid <- function(data,
                    timevar, cohortvar, unitvar, outcomevar, 
                    control_option="both",result_type="group_time", balanced_event_time = NA,
                    control_type = "ipw", allow_unbalance_panel = FALSE, boot=FALSE, biters = 1000, cband = FALSE, alpha = 0.05,
                    weightvar=NA,clustervar=NA, covariatesvar = NA, varycovariatesvar = NA, 
                    copy = TRUE, validate = TRUE,
                   anticipation = 0, anticipation2 = 0, base_period = "universal",
                   exper = NULL, full = FALSE, parallel = FALSE,
                    cohortvar2 = NA, event_specific = TRUE, double_control_option="both"){
  
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