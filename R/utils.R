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