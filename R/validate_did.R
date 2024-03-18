validate_did <- function(dt,covariatesvar, varycovariatesvar,varnames, balanced_event_time, allow_unbalance_panel){
  raw_unit_size <- dt[, uniqueN(unit)]
  raw_time_size <- dt[, uniqueN(time)]
  
  
  if(!is.na(balanced_event_time)){
    if(balanced_event_time > dt[, max(time-G)]){stop("balanced_event_time is larger than the max event time in the data")}
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
  
  
  if(!allNA(covariatesvar) && uniqueN(dt, by = c("unit", covariatesvar)) > raw_unit_size){
    warning("some covariates is time-varying, fastdid only use the first observation for covariates.")
  }
  
  
  if(!allNA(covariatesvar)|!allNA(varycovariatesvar)){
    for(cov in c(covariatesvar, varycovariatesvar)){
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
  if(!allow_unbalance_panel){
    unit_count <- dt[, .(count = .N), by = unit]
    if(any(unit_count[, count < raw_time_size])){
      mis_unit <- unit_count[count < raw_time_size]
      warning(nrow(mis_unit), " units is missing in some periods, enforcing balanced panel by dropping them")
      dt <- dt[!unit %in% mis_unit[, unit]]
    }
  }
  return(dt)
}