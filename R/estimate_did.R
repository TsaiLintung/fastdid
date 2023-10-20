estimate_did <- function(dt_did, ipw_formula, 
                         last_coef = NULL, cache_ps_fit, cache_hess){
  oldn <- dt_did[, .N]
  data_pos <- whichNA(dt_did[, D], invert = TRUE)
  dt_did <- dt_did[data_pos]
  n <- dt_did[, .N]

  ipw <- !is.null(ipw_formula)
  if(ipw){
    
    dt_did[, const := 1]

    
    if(is.null(cache_ps_fit)|is.null(cache_hess)){ #if no cache, calcuate ipw

      #estimate the logit
      logit_formula <- paste0("D ~", ipw_formula)
      # prop_score_est <- suppressWarnings(speedglm(as.formula(logit_formula), 
      #                                             data = dt_did,
      #                                             family = stats::binomial(), fitted = TRUE, start = last_coef,
      #                                             weights = dt_did[, weights]))
      
      prop_score_est <- suppressWarnings(parglm(as.formula(logit_formula),
                                                  data = dt_did,
                                                  family = stats::binomial(), start = last_coef,
                                                  weights = dt_did[, weights],
                                                  control = parglm.control(nthreads = getDTthreads())))
      
      covvars <-  c("const", names(prop_score_est$coefficients)[names(prop_score_est$coefficients) != "(Intercept)"])
      hess <- stats::vcov(prop_score_est) * n #for the influence function
      
      logit_coef <-  prop_score_est$coefficients 
      logit_coef[is.na(logit_coef)|abs(logit_coef) > 1e10] <- 0 #put extreme value and na to 0
      prop_score_fit <- fitted(prop_score_est)
      if(max(prop_score_fit) >= 1){warning(paste0("support overlap condition violated for some group_time"))}
      prop_score_fit <- pmin(1-1e-16, prop_score_fit) #for the ipw
      
    } else { #when using multiple outcome, ipw cache can be reused
      hess <- cache_hess
      prop_score_fit <- cache_ps_fit
      covvars <-  c("const", colnames(cache_hess)[colnames(cache_hess) != "(Intercept)"])
      logit_coef <- NULL #won't be needing the approximate cache
    }
    
    dt_did[, ps := prop_score_fit]
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*ps*(1-D)/(1-ps)]
    
  } else {
    
    prop_score_fit <- rep(1,n)
    logit_coef <- NULL
    hess <- NULL
    
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*(1-D)]
    
  }
  
  
  #did with ipw
  
  dt_did[, delta_y := post.y-pre.y]
  
  dt_did[, att_treat := treat_ipw_weight*delta_y]
  dt_did[, att_cont := cont_ipw_weight*delta_y]
  
  weighted_treat_delta <- dt_did[,sum(att_treat)/sum(treat_ipw_weight)]
  weighted_cont_delta <- dt_did[,sum(att_cont)/sum(cont_ipw_weight)]
  
  att <- weighted_treat_delta - weighted_cont_delta
  
  # influence from ipw
  if(ipw){
    score_ps <- as.matrix(dt_did[, weights*(D-ps)*.SD, .SDcols = covvars]) 
    asym_linear_ps <- score_ps %*% hess |> reverse_col()
    M2 <- dt_did[, .(cont_ipw_weight*(delta_y-weighted_cont_delta)*.SD), .SDcols = covvars][, lapply(.SD, mean), .SDcols = covvars] |> 
      as.matrix() |> reverse_col()
    inf_cont_ipw <- asym_linear_ps %*% as.matrix(M2)
  } else {inf_cont_ipw <- 0}
  
  # influence from did
  inf_cont_did <- dt_did[, att_cont - cont_ipw_weight*weighted_cont_delta]
  inf_cont <- (inf_cont_did+inf_cont_ipw)/dt_did[, mean(cont_ipw_weight)]
  inf_treat <- dt_did[, (att_treat - treat_ipw_weight*weighted_treat_delta)/mean(treat_ipw_weight)]
  
  inf_func_no_na <- inf_treat - inf_cont
  
  inf_func <- rep(0, oldn) #the default needs to be 0 for the matrix multiplication
  inf_func_no_na <- inf_func_no_na * oldn / n #adjust the value such that mean over the whole id size give the right result
  inf_func[data_pos] <- inf_func_no_na
  return(list(att = att, inf_func = inf_func, logit_coef = logit_coef, #for next gt
              cache_ps_fit = prop_score_fit, cache_hess = hess)) #for next outcome
}
