estimate_gtatt <- function(outcomes_list, outcomevar, covariates, ipw_formula, weights, 
                           cohort_sizes,cohorts,id_size,time_periods,
                           control_option) {

  treated_cohorts <- cohorts[!is.infinite(cohorts)]
  
  if(!Inf %in% cohorts & control_option != "notyet"){
    warning("no never-treated availble, switching to not-yet-treated control")
    control_option <- "notyet"
  }
  max_control_cohort <- ifelse(control_option == "notyet", max(treated_cohorts), Inf) 
  
  cache_ps_fit_list <- list() #for the first outcome won't be able to use cache, empty list returns null for call like cache_ps_fit[["1.3"]]
  cache_hess_list <- list()
  outcome_result_list <- list()
  
  for(outcol in outcomevar){
    
    outcomes <- outcomes_list[[outcol]]
    last_coef <- NULL
    gt <- data.table()
    gt_att <- c()
    gt_inf_func <- data.table(placeholder = rep(NA, id_size))
    
    for(t in time_periods){
      for(g in cohorts){
        
        gt_name <- paste0(g, ".", t)
        
        #setup and checks
        base_period <- g-1
        min_control_cohort <- ifelse(control_option == "never", Inf, max(t+1, base_period+1)) #not-yet treated / never treated in both base and "treated" period
        if(t == base_period){next} #no treatmenteffect for the base period
        if(base_period < min(time_periods)){next} #no treatmenteffect for the first period, since base period is not observed
        if(g >= max_control_cohort){next} #no treatmenteffect for never treated or the last treated cohort (notyet)
        if(t >= max_control_cohort){next} #no control available if the last cohort is treated too
        
        #select the control and treated cohorts
        did_setup <- rep(NA, id_size)
        did_setup[get_cohort_pos(cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
        did_setup[get_cohort_pos(cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
        
        #construct the 2x2 dataset
        cohort_did <- data.table(did_setup, outcomes[[t]], outcomes[[base_period]], weights, covariates)
        setnames(cohort_did, c("did_setup", "V2", "V3", "weights"), c("D", "post.y", "pre.y", "weights"))
        
        #estimate did
        result <- estimate_did(cohort_did, ipw_formula, last_coef, cache_ps_fit_list[[gt_name]], cache_hess_list[[gt_name]])
      
        #collect the result
        last_coef <- result$logit_coef
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
    outcome_result <- list(att = gt_att, inf_func = gt_inf_func, gt = gt)
    outcome_result_list[[outcol]] <- outcome_result
    
    rm(outcome_result)

  }
  
  return(outcome_result_list)
}

get_cohort_pos <- function(cohort_sizes, start_cohort, end_cohort = start_cohort){
  #if(!start_cohort %in% cohort_sizes[,unique(G)]|!end_cohort %in% cohort_sizes[,unique(G)]) {stop("cohort not in cohort_sizes")}
  start <- cohort_sizes[G < start_cohort, sum(cohort_size)]+1
  end <- cohort_sizes[G <= end_cohort, sum(cohort_size)]
  return(start:end)
}
