validate_did <- function(dt,covariatesvar,varnames){
  raw_unit_size <- dt[, uniqueN(unit)]
  raw_time_size <- dt[, uniqueN(time)]
  
  #doesn't allow missing value for now
  for(col in c(covariatesvar, varnames)){
    if(is.null(col)){next}
    na_obs <- whichNA(dt[[col]])
    if(length(na_obs) != 0){
      warning("missing values detected in ", col, ", removing ", length(na_obs), " observation.")
      dt <- dt[!na_obs]
    }
  }
  
  if(!is.null(covariatesvar)){
    
    for(cov in covariatesvar){
      #doesn't allow time varying covariates
      if(!all(dt[, uniqueN(get(cov)) == 1, by = "unit"])){stop(cov, " is time-varying")}
      #check covaraites is not constant  
      if(fnunique(dt[, get(cov)[1], by = "unit"][, V1]) == 1)stop(cov, " have no variation")
    }

  }
  
  #check balanced panel
  #check if any is dup
  if(anyDuplicated(dt[, .(unit, time)])){
    dup <- duplicated(dt[,.(unit, time)])
    warning(nrow(dup), " units is observed more than once in some periods, enforcing balanced panel by dropping them")
    dt <- dt[!unit %in% dup[, unit]]
  }
  
  #check if any is missing
  unit_count <- dt[, .(count = .N), by = unit]
  if(any(unit_count[, count < raw_time_size])){
    mis_unit <- unit_count[count < raw_time_size]
    warning(nrow(mis_unit), " units is missing in some periods, enforcing balanced panel by dropping them")
    dt <- dt[!unit %in% mis_unit[, unit]]
  }
  return(dt)
}