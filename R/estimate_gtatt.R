estimate_gtatt <- function(aux, p){
  caches <- list()
  outcome_results <- list()
  for(outcol in p$outcomevar){
    y <- aux$outcomes_list[[outcol]]
    out_result <- estimate_gtatt_outcome(y, aux, p, caches)
    out_result$est$outname <- outcol
    outcome_results[[outcol]] <- out_result$est
    caches <- out_result$caches
  }
  return(outcome_results)
}

#gtatt for each outcome
estimate_gtatt_outcome <- function(y, aux, p, caches) {

    gt_all <- expand.grid(g = aux$cohorts, t = aux$time_periods, stringsAsFactors = FALSE) |> transpose() |> as.list() #first loop t then g
    
    #main estimation 
    if(!p$parallel){
      gt_results <- lapply(gt_all, estimate_gtatt_outcome_gt, y, aux, p, caches)
    } else {
      gt_results <- mclapply(gt_all, estimate_gtatt_outcome_gt, y, aux, p, caches, mc.cores = getDTthreads())
    }
    
    #post process
    gt_results <- gt_results[which(!sapply(gt_results, is.null))] #remove the ones with no valid didsetup
    if(length(gt_results) == 0){stop("no valid group-times att to compute")}
    
    gt <- lapply(gt_results, function(x) {x$gt}) |> as.data.table() |> transpose()
    names(gt) <- c("G", "time")
    gt[, time := as.integer(time)]
    
    
    
    gt_att <- lapply(gt_results, function(x) {x$result$att})
    gt_inf_func <- lapply(gt_results, function(x) {x$result$inf_func})
    caches <- lapply(gt_results, function(x) {x$result$cache})
    
    gt_names <- gt[,paste0(G,".",time)]
    names(gt_att) <- gt_names
    names(gt_inf_func) <- gt_names
    names(caches) <- gt_names
    
    gt_inf_func <- do.call(cbind, gt_inf_func)
    gt_att <- do.call(cbind, gt_att) |> t()

  return(list(est = list(gt = gt, att = gt_att, inf_func = gt_inf_func), caches = caches))    
}

#gtatt for each outcome, each gt
estimate_gtatt_outcome_gt <- function(gt, y, aux, p, caches){
  
  g <- gt[1]
  t <- as.numeric(gt[2])
  
  gt_name <- paste0(g,".",t)
  base_period <- get_base_period(g,t,p)
  did_setup <- get_did_setup(g,t, base_period, aux, p)
  
  if(is.null(did_setup)){return(NULL)} #no gtatt if no did setup
  
  #covariates matrix
  covvars <- get_covvars(base_period, t, aux, p)
  
  #the 2x2 dataset
  cohort_did <- data.table(did_setup, y[[t]], y[[base_period]], aux$weights)
  names(cohort_did) <- c("D", "post.y", "pre.y", "weights")
  
  # estimate --------------------
  result <- tryCatch(estimate_did(dt_did = cohort_did, covvars, p, caches[[gt_name]]),
                     error = function(e){stop("2-by-2 DiD failed for internal group-time ",g , 
                                              "-", t, ": ", e)})
  return(list(gt = gt, result = result))
  
}


get_base_period <- function(g,t,p){
  g <- ming(g) #for two period
  if(p$base_period == "universal"){
    base_period <- g-1-p$anticipation
  } else {
    base_period <- ifelse(t>=g, g-1-p$anticipation, t-1)
  }
  return(base_period)
}

get_did_setup <- function(g, t, base_period, aux, p){
  
  treated_cohorts <- aux$cohorts[!is.infinite(ming(aux$cohorts))]
  
  #get the range of cohorts
  if(p$control_option == "never"){
    min_control_cohort <- Inf
  } else {
    min_control_cohort <- max(t, base_period)+p$anticipation+1
  }
  max_control_cohort <- ifelse(p$control_option == "notyet", max(ming(treated_cohorts)), Inf) 
  
  #experimental
  if(!is.na(p$exper$max_control_cohort_diff)){
    max_control_cohort <- min(g+p$exper$max_control_cohort_diff, max(treated_cohorts))
  } 
  if(!is.na(p$exper$min_control_cohort_diff)){
    min_control_cohort <- max(g+p$exper$min_control_cohort_diff, min(treated_cohorts))
  } 
  if((!is.na(p$exper$max_dynamic))){
    if(t-ming(g) > p$exper$max_dynamic){return(NULL)}
  }
  if(!is.na(p$exper$min_dynamic)){
    if(t-ming(g) < p$exper$min_dynamic){return(NULL)}
  }
  
  # invalid gt
  if(t == base_period | #no treatment effect for the base period
     base_period < min(aux$time_periods) | #no treatment effect for the first period, since base period is not observed
     base_period > max(aux$time_periods) |
     ming(g) >= min(max_control_cohort) | #no treatment effect for never treated or the last treated cohort (for not yet notyet)
     t >= max_control_cohort | #no control available if the last cohort is treated too
     min_control_cohort > max_control_cohort){ #no control avalilble, most likely due to anticipation
    return(NULL)
  } 
  
  #select the control and treated cohorts
  did_setup <- rep(NA, aux$id_size)
  did_setup[get_control_pos(aux$cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
  did_setup[get_treat_pos(aux$cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
  
  if(!is.na(p$exper$filtervar)){
    did_setup[!aux$filters[[base_period]]] <- NA #only use units with filter == TRUE at base period
   
  }
  if(!is.na(p$exper$filtervar_post)){
    did_setup[!aux$filters[[t]]] <- NA #only use units with filter == TRUE at target period
  }
  
  if(all(is.na(did_setup))){return(NULL)}
  return(did_setup)
}

get_control_pos <- function(cohort_sizes, start_cohort, end_cohort = start_cohort){
  start <- cohort_sizes[ming(G) < start_cohort, sum(cohort_size)]+1 
  end <- cohort_sizes[ming(G) <= end_cohort, sum(cohort_size)]
  return(start:end)
}

get_treat_pos <- function(cohort_sizes, treat_cohort){
  index <- which(cohort_sizes[,G] == treat_cohort)
  start <- ifelse(index == 1, 1, cohort_sizes[1:(index-1), sum(cohort_size)]+1)
  end <- cohort_sizes[1:index, sum(cohort_size)]
  return(start:end)
}

get_covvars <- function(base_period, t, aux, p){
  
  if(all(is.na(p$covariatesvar)) & all(is.na(p$varycovariatesvar))){return(NA)}
  covvars <- data.table()
  
  #add time-varying covariates
  if(!all(is.na(p$varycovariatesvar))){
    
    precov <- aux$varycovariates[[base_period]]
    names(precov) <- paste0("pre_", names(precov))
    postcov <- aux$varycovariates[[t]]-aux$varycovariates[[base_period]]
    names(postcov) <- paste0("post_", names(aux$varycovariates[[t]]))
    covvars <- cbind(covvars, cbind(precov, postcov))
  }
  
  #add time-invariant covariates
  if(!all(is.na(p$covariatesvar))){
    covvars <- cbind(aux$covariates, covvars)
  }
  
  #add constant
  covvars <- as.matrix(cbind(const = -1, covvars))
  return(covvars)
}
