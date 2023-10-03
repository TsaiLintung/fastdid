fastdid <- function(dt,
                    timevar, cohortvar, unitvar,
                    weightvar=NULL,clustervar=NULL,covariatesvar = c(),
                    control_option="both",result_type="group_time,
                    boot=FALSE){
  
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))
  
  setorder(dt, time, G, unit)
  
  time_periods <- dt[, unique(time)]
  time_size <- length(time_periods)
  id_size <- dt[, uniqueN(unit)]
  
  dt_inv <- dt[1:id_size]
  cohorts <- dt_inv[, unique(G)]
  cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]
  
  outcomes <- list()
  for(i in time_periods){
    start <- (i-1)*id_size+1
    end <- i*id_size
    outcomes[[i]] <- dt[seq(start,end), .(y)]
  }
  
  # fill 
  
  if(length(covariatesvar)>0){
    covariates <- dt_inv[,.SD, .SDcols = covariatesvar]
  } else {covariates <- c()}
  
  if(!is.null(clustervar)){
    cluster <- dt_inv[,.SD, .SDcols = clustervar]
  } else {cluster <- NULL}
  
  if(!is.null(weightvar)){
    weights <- dt_inv[, .SD, .SDcols = weightvar]
  } else {weights <- rep(1,id_size)}
  
  
  # validate data -----------------------------------------------------
  
  if(any(sapply(covariates, sd) == 0)){stop("some covariates have no variation")}
  
  # the main process -------------------------------------------------
  
  
  # attgt
  gt_result <- estimate_gtatt(outcomes, covariates, weights,
                               cohort_sizes,cohorts,id_size,time_periods, #info about the dt
                               control_option)
  
  # aggregate att and inf function
  agg_result <- aggregate_gt_result(gt_result, cohort_sizes,
                                    result_type)
  
  
  #influence function -> se
  agg_se <- get_se(agg_result$inf_matrix, boot, biters = 1000, cluster)
  
  #gather results
  result <- data.table(target = agg_result$targets, att = agg_result$agg_att, se = agg_se)
  
  setnames(dt, c("time", "G", "unit"), c(timevar, cohortvar, unitvar))
  return(result)
}