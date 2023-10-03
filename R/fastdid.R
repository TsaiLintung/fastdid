fastdid <- function(dt,
                    timevar, cohortvar, unitvar,
                    weightvar=NULL,clustervar=NULL,covariatesvar = c(),
                    control_option="both",result_type="group_time",
                    boot=FALSE){
  
  
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
  agg_se <- get_se(agg_result$inf_matrix, boot, biters = 1000, cluster)
  
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