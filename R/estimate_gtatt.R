estimate_gtatt <- function(aux, p) {

  if(!Inf %in% aux$cohorts & p$control_option != "notyet"){
    warning("no never-treated availble, switching to not-yet-treated control")
    p$control_option <- "notyet"
  }
  
  treated_cohorts <- aux$cohorts[!is.infinite(aux$cohorts)]
  max_control_cohort <- ifelse(p$control_option == "notyet", max(treated_cohorts), Inf) 
  
  cache_ps_fit_list <- list() #for the first outcome won't be able to use cache, empty list returns null for call like cache_ps_fit[["1.3"]]
  cache_hess_list <- list()
  outcome_result_list <- list()
  for(outcol in p$outcomevar){
    
    outcomes <- aux$outcomes_list[[outcol]]
    last_coef <- NULL
    gt <- data.table()
    gt_att <- c()
    gt_inf_func <- data.table(placeholder = rep(NA, aux$id_size)) #populate with NA
    
    for(t in aux$time_periods){
      for(g in aux$cohorts){
        
        gt_name <- paste0(g, ".", t)
        
        #setup and checks
        base_period <- g-1
        min_control_cohort <- ifelse(p$control_option  == "never", Inf, max(t+1, base_period+1)) #not-yet treated / never treated in both base and "treated" period
        if(t == base_period){next} #no treatment effect for the base period
        if(base_period < min(aux$time_periods)){next} #no treatment effect for the first period, since base period is not observed
        if(g >= max_control_cohort){next} #no treatment effect for never treated or the last treated cohort (for not yet notyet)
        if(t >= max_control_cohort){next} #no control available if the last cohort is treated too
        
        #select the control and treated cohorts
        did_setup <- rep(NA, aux$id_size)
        did_setup[get_cohort_pos(aux$cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
        did_setup[get_cohort_pos(aux$cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
        
        #construct the covariates matrix
        covvars <- get_covvars(aux$covariates, aux$varycovariates, base_period, t)
        
        #construct the 2x2 dataset
        cohort_did <- data.table(did_setup, outcomes[[t]], outcomes[[base_period]], aux$weights)
        names(cohort_did) <- c("D", "post.y", "pre.y", "weights")
        
        #estimate did
        if(!p$allow_unbalance_panel){
          result <- estimate_did(cohort_did, covvars, p$control_type, 
                                 last_coef, cache_ps_fit_list[[gt_name]], cache_hess_list[[gt_name]]) #cache
        } else {
          result <- estimate_did_rc(cohort_did, covvars, p$control_type, 
                                    last_coef, cache_ps_fit_list[[gt_name]], cache_hess_list[[gt_name]]) #cache
        }
       
        #collect the result
        last_coef <- result$logit_coef #for faster convergence in next iter
        gt <- rbind(gt, data.table(G = g, time = t)) #the sequence matters for the weights
        gt_att <- c(gt_att, att = result$att)
        gt_inf_func[[gt_name]] <- result$inf_func
        
        #assign cache for next outcome
        if(is.null(cache_ps_fit_list[[gt_name]])){cache_ps_fit_list[[gt_name]] <- result$cache_ps_fit}
        if(is.null(cache_hess_list[[gt_name]])){cache_hess_list[[gt_name]] <- result$cache_hess}
        rm(result)
        
      }
    }
    
    gt_inf_func[,placeholder := NULL]
    names(gt_att) <- names(gt_inf_func)
    outcome_result_list[[outcol]] <- list(att = gt_att, inf_func = gt_inf_func, gt = gt, outname = outcol)

  }
  
  return(outcome_result_list)
}

get_cohort_pos <- function(cohort_sizes, start_cohort, end_cohort = start_cohort){
  start <- cohort_sizes[G < start_cohort, sum(cohort_size)]+1
  end <- cohort_sizes[G <= end_cohort, sum(cohort_size)]
  return(start:end)
}

get_covvars <- function(covariates, varycovariates, base_period, t){
  
  if(is.list(varycovariates)){
    
    precov <- varycovariates[[base_period]]
    names(precov) <- paste0("pre_", names(precov))
    postcov <- varycovariates[[t]]-varycovariates[[base_period]]
    names(postcov) <- paste0("post_", names(varycovariates[[t]]))
    
    if(is.data.table(covariates)){
      covvars <- cbind(covariates, cbind(precov, postcov)) #if covariates is na
    } else {
      covvars <- cbind(const = -1, cbind(precov, postcov)) #if covariates is na
    }
    covvars <- as.matrix(covvars)
    
  } else {
    if(is.data.table(covariates)){
      covvars <- as.matrix(covariates) #will be NA if covariatesvar have nothing
    } else {
      covvars <- NA
    }
  }
  return(covvars)
}
