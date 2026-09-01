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

#' gtatt for each outcome.
#' @noRd
estimate_gtatt_outcome <- function(y, aux, p, caches) {
    
    treated_cohort <- aux$cohorts[!is.infinite(ming(aux$cohorts))] #otherwise would try to calculate the pre-period of nevertreated in varying base period lol
    gt_all <- expand.grid(g = treated_cohort, t = aux$time_periods, stringsAsFactors = FALSE) |> transpose() |> as.list() #first loop t then g
    
    #main estimation 
    if(!p$parallel){
      gt_results <- lapply(gt_all, estimate_gtatt_outcome_gt, y, aux, p, caches)
    } else {
      gt_results <- mclapply(gt_all, estimate_gtatt_outcome_gt, y, aux, p, caches)
    }
    
    #post process
    gt_results <- gt_results[which(!sapply(gt_results, is.null))] #remove the ones with no valid didsetup

    # cells skipped for a too-small group are reported once, not one warning per cell
    ess_skip <- vapply(gt_results, function(x) isTRUE(x$skip_ess), logical(1))
    if(any(ess_skip)){
      warning(sum(ess_skip), " group-time(s) skipped: a group has fewer than 2 effective units, so the variance is not estimable")
      gt_results <- gt_results[!ess_skip]
    }
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

#' gtatt for each outcome, each gt.
#' @noRd
estimate_gtatt_outcome_gt <- function(gt, y, aux, p, caches){
  
  g <- gt[1]
  t <- as.numeric(gt[2])
  
  # skip g-t outside requested event-time range
  event_time <- t - ming(g)
  if (!is.na(p$exper$only_est_min) && event_time < p$exper$only_est_min) return(NULL)
  if (!is.na(p$exper$only_est_max) && event_time > p$exper$only_est_max) return(NULL)

  #find base time
  gt_name <- paste0(g,".",t)
  base_period <- get_base_period(g,t,p)
  if(t == base_period || #no treatment effect for the base period
     !base_period %in% aux$time_periods){ #base period out of bounds
    return(NULL)
  } 
  #find treatment and control group
  did_setup <- get_did_setup(g,t, base_period, aux, p)
  valid_tc_groups <- any(did_setup == 1) && any(did_setup == 0) #if takes up too much time, consider use collapse anyv, but right now quite ok
  if(!isTRUE(valid_tc_groups)){return(NULL)} #no treatment group or control group #isTRUE for na as false
  
  #covariates matrix
  covvars <- get_covvars(base_period, t, aux, p)
  
  #the 2x2 dataset
  cohort_did <- data.table(did_setup, y[[t]], y[[base_period]], aux$weights)
  names(cohort_did) <- c("D", "post.y", "pre.y", "weights")

  # skip if no valid observation in either period when allowing for
  # unbalanced panels
  if(sum(!is.na(cohort_did$pre.y)) == 0 || sum(!is.na(cohort_did$post.y)) == 0){
    return(NULL)
  }

  # estimate --------------------
  result <- tryCatch(estimate_did(dt_did = cohort_did, covvars, p, caches[[gt_name]]),
                     error = function(e){
                       # small-group skips are counted and reported once by the caller
                       if(grepl("fewer than 2 effective units", e$message, fixed = TRUE)){
                         return("skip_ess")
                       }
                       warning("Skipping group-time ", g, "-", t,
                               ": ", e$message)
                       return(NULL)
                     })
  if(is.null(result)){return(NULL)}
  if(identical(result, "skip_ess")){return(list(gt = gt, skip_ess = TRUE))}
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

  min_control_cohort <- ifelse(p$control_option == "never", Inf, max(t, base_period)+p$anticipation+1)
  max_control_cohort <- ifelse(p$control_option == "notyet", max(ming(treated_cohorts)), Inf) 

  if(!is.na(p$exper$max_control_cohort_diff)){
    max_control_cohort <- min(g+p$exper$max_control_cohort_diff, max_control_cohort)
  } 
  
  #select the control and treated cohorts
  did_setup <- rep(NA, aux$id_size)
  if(allNA(p$cohortvar2) || p$anticipation == p$anticipation2){
    # one horizon for every event, so the cohorts in range are contiguous after the sort
    control_pos <- get_control_pos(aux$cohort_sizes, min_control_cohort, max_control_cohort)
  } else {
    # each event has its own horizon, so screen the cohorts one by one
    control_pos <- get_control_pos_event(aux$cohort_sizes, t, base_period, max_control_cohort, p)
  }
  did_setup[control_pos] <- 0
  did_setup[get_treat_pos(aux$cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
  
  if(!is.na(p$exper$filtervar)){
    did_setup[!aux$filters[[base_period]]] <- NA #only use units with filter == TRUE at base period
   
  }
  if(!is.na(p$exper$filtervar_post)){
    did_setup[!aux$filters[[t]]] <- NA #only use units with filter == TRUE at target period
  }
  
  if(length(did_setup) != aux$id_size){stop("internal bug: something wrong with did setup (again?)")}
  
  return(did_setup)
}

get_control_pos <- function(cohort_sizes, start_cohort, end_cohort = start_cohort){
  start <- cohort_sizes[ming(G) < start_cohort, sum(cohort_size)]+1 
  end <- cohort_sizes[ming(G) <= end_cohort, sum(cohort_size)]
  
  # Validate sequence parameters
  if(start > end || start <= 0 || end <= 0){
    return(c())  # Return empty vector when no valid control cohorts
  }
  return(seq(start, end, by = 1))
}

#' Control positions when the events have different anticipation horizons.
#'
#' A cohort is a valid not-yet-treated control only if every event is further
#' away than the horizon of that event. The first event uses `anticipation`, the
#' confounding events use `anticipation2`. The cohorts in range are not
#' contiguous after the sort, so screen them one by one.
#'
#' @param cohort_sizes the cohort table, in the order of the unit array.
#' @param t,base_period the two periods of the 2x2.
#' @param max_control_cohort the not-yet-treated upper bound.
#' @return the positions of the control units in the unit array.
#' @noRd
get_control_pos_event <- function(cohort_sizes, t, base_period, max_control_cohort, p){
  GG <- cohort_sizes[, G]
  last <- max(t, base_period)

  if(p$control_option == "never"){
    keep <- is.infinite(ming(GG))
  } else {
    keep <- g1(GG) > last + p$anticipation
    M <- 1L + length(p$cohortvar2)
    for(d in 2:M){
      keep <- keep & (gd(GG, d) > last + p$anticipation2)
    }
  }
  keep <- keep & (ming(GG) <= max_control_cohort)
  if(!any(keep)){return(c())}

  # each cohort is one contiguous block of the unit array
  end <- cumsum(cohort_sizes[, cohort_size])
  start <- end - cohort_sizes[, cohort_size] + 1
  return(unlist(lapply(which(keep), function(i) seq(start[i], end[i]))))
}

get_treat_pos <- function(cohort_sizes, treat_cohort){ #need to separate for double did to match exact g-g-t
  index <- which(cohort_sizes[,G] == treat_cohort)
  if(length(index) == 0){
    stop("Cohort ", treat_cohort, " not found in cohort_sizes")
  }
  start <- ifelse(index == 1, 1, cohort_sizes[1:(index-1), sum(cohort_size)]+1)
  end <- cohort_sizes[1:index, sum(cohort_size)]
  
  # Validate sequence parameters
  if(start > end || start <= 0 || end <= 0){
    return(c())  # Return empty vector when no valid treat positions
  }
  return(seq(start, end, by = 1))
}

get_covvars <- function(base_period, t, aux, p){
  
  if(allNA(p$covariatesvar) && allNA(p$varycovariatesvar)){return(NA)}
  covvars <- data.table()
  
  #add time-varying covariates
  if(!allNA(p$varycovariatesvar)){
    
    precov <- aux$varycovariates[[base_period]]
    names(precov) <- paste0("pre_", names(precov))
    postcov <- aux$varycovariates[[t]]-aux$varycovariates[[base_period]]
    names(postcov) <- paste0("post_", names(aux$varycovariates[[t]]))
    covvars <- cbind(covvars, cbind(precov, postcov))
  }
  
  #add time-invariant covariates
  if(!allNA(p$covariatesvar)){
    covvars <- cbind(aux$covariates, covvars)
  }
  
  #add constant
  covvars <- as.matrix(cbind(const = -1, covvars))
  return(covvars)
}
