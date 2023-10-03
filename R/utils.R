get_att_in_ci_ratio <- function(merged_att){
  
  merged_att[, ci_ub := Estimate+`Std. Error`*1.96]
  merged_att[, ci_lb := Estimate-`Std. Error`*1.96]
  merged_att[, s :=  (function(x){str_sub(str_trim(x), str_length(str_trim(x)), str_length(str_trim(x)))})(variable)]
  merged_att[, attgt := attgt*as.numeric(s)]
  merged_att[, par_in_ci := (attgt <= ci_ub & attgt >= ci_lb)]
  ratio <- merged_att[, mean(par_in_ci)]
  
  return(ratio)
  
}

set_max_thread <- function(){
  setDTthreads(0)
  options(kit.nThread = getDTthreads())
  setFixest_nthreads(getDTthreads())
}

reverse_col <- function(x){
  return(x[,ncol(x):1])
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