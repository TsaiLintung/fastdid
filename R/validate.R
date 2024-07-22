validate_argument <- function(p, dt_names){

  #release p
  for(name in names(p)){
    assign(name, p[[name]])
  }
  
  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_set_arg(timevar, unitvar, cohortvar, "match", .choices = dt_names, .message = name_message)
  
  covariate_message <- "__ARG__ must be NA or a character vector which are all names of columns from the dataset."
  check_set_arg(varycovariatesvar, covariatesvar, outcomevar, 
                "NA | multi match", .choices = dt_names, .message = covariate_message)
  
  checkvar_message <- "__ARG__ must be NA or a character scalar if a name of columns from the dataset."
  check_set_arg(weightvar, clustervar, filtervar,
                "NA | match", .choices = dt_names, .message = checkvar_message)
  
  check_set_arg(control_option, "match", .choices = c("both", "never", "notyet")) #kinda bad names since did's notyet include both notyet and never
  check_set_arg(control_type, "match", .choices = c("ipw", "reg", "dr")) 
  check_set_arg(base_period, "match", .choices = c("varying", "universal"))
  check_arg(copy, validate, boot, allow_unbalance_panel, "scalar logical")
  check_arg(max_control_cohort_diff, min_control_cohort_diff, anticipation, "scalar numeric")
  
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
  
  # coerce non-sensible option
  if(!is.na(clustervar) && unitvar == clustervar){clustervar <- NA} #cluster on id anyway, would cause error otherwise
  if((!is.infinite(max_control_cohort_diff) | !is.infinite(min_control_cohort_diff)) & control_option == "never"){
    warning("control_cohort_diff can only be used with not yet")
    p$control_option <- "notyet"
  }
}

validate_dt <- function(dt,varnames,p){

  raw_unit_size <- dt[, uniqueN(unit)]
  raw_time_size <- dt[, uniqueN(time)]
  
  if(!is.na(p$balanced_event_time)){
    if(p$balanced_event_time > dt[, max(time-G)]){stop("balanced_event_time is larger than the max event time in the data")}
  }
  
  if(!is.na(p$filtervar) && !is.logical(dt[[p$filtervar]])){
    stop("filter var needs to be a logical column")
  }
  
  #doesn't allow missing value for now
  for(col in varnames){
    if(is.na(col)){next}
    na_obs <- whichNA(dt[[col]])
    if(length(na_obs) != 0){
      warning("missing values detected in ", col, ", removing ", length(na_obs), " observation.")
      dt <- dt[!na_obs]
    }
  }
  
  if(!allNA(p$covariatesvar) && uniqueN(dt, by = c("unit", p$covariatesvar)) > raw_unit_size){
    warning("some covariates is time-varying, fastdid only use the first observation for covariates.")
  }
  
  
  if(!allNA(p$covariatesvar)|!allNA(p$varycovariatesvar)){
    for(cov in c(p$covariatesvar, p$varycovariatesvar)){
      if(is.na(cov)){next}
      #check covaraites is not constant  
      if(fnunique(dt[, get(cov)[1], by = "unit"][, V1]) == 1)stop(cov, " have no variation")
    }
  }
  
  #check balanced panel
  #check if any is dup
  if(anyDuplicated(dt[, .(unit, time)])){
    dup_id <- dt[duplicated(dt[,.(unit, time)]), unique(unit)]
    stop(length(dup_id), " units is observed more than once in a period.")
  }
  
  #check if any is missing
  if(!p$allow_unbalance_panel){
    unit_count <- dt[, .(count = .N), by = unit]
    if(any(unit_count[, count < raw_time_size])){
      mis_unit <- unit_count[count < raw_time_size]
      warning(nrow(mis_unit), " units is missing in some periods, enforcing balanced panel by dropping them")
      dt <- dt[!unit %in% mis_unit[, unit]]
    }
  }
  return(dt)
}