estimate_did <- function(dt_did, last_coef = NULL){
  oldn <- dt_did[, .N]
  data_pos <- whichNA(dt_did[, D], invert = TRUE)
  dt_did <- dt_did[data_pos]
  n <- dt_did[, .N]
  
  covvars <- names(dt_did)[str_starts(names(dt_did), "cov")]
  if(length(covvars) > 0){
    dt_did[, const := 1]
    covvars <- c(covvars, "const")
    ipw <- TRUE
  } else {
    ipw <- FALSE
  }
  
  
  #fastglm not possible
  #the intercept is essential
  if(ipw){
    prop_score_est <- suppressWarnings(speedglm(dt_did[,D] ~ as.matrix(dt_did[,.SD, .SDcols = covvars]), 
                                                family = binomial(), fitted = TRUE, start = last_coef,
                                                weights = dt_did[, weights]))
    logit_coef <-  prop_score_est$coefficients |> nafill(fill = 0)
    prop_score_fit <- fitted(prop_score_est)
    prop_score_fit <- pmin(1-1e-16, prop_score_fit)
    
    dt_did[, ps := prop_score_fit]
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*ps*(1-D)/(1-ps)]
    
  } else {
    prop_score_fit <- rep(1,n)
    logit_coef <- NULL
    
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
    score_ps <- as.matrix(dt_did[, weights*(D-ps)*.SD, .SDcols = covvars]) |> reverse_col()
    hess <- vcov(prop_score_est) * n
    asym_linear_ps <- score_ps %*% hess
    M2 <- dt_did[, .(cont_ipw_weight*(delta_y-weighted_cont_delta)*.SD), .SDcols = covvars][, lapply(.SD, mean), .SDcols = covvars] |> 
      as.matrix()|> reverse_col()
    inf_cont_ipw <- asym_linear_ps %*% as.matrix(M2)
  } else {inf_cont_ipw <- 0}

  # influence from did
  inf_cont_did <- dt_did[, att_cont - cont_ipw_weight*weighted_cont_delta]
  inf_cont <- (inf_cont_did+inf_cont_ipw)/dt_did[, mean(cont_ipw_weight)]
  inf_treat <- dt_did[, (att_treat - treat_ipw_weight*weighted_treat_delta)/mean(treat_ipw_weight)]
  
  inf_func_no_na <- inf_treat - inf_cont
  
  inf_func <- rep(0, oldn)
  inf_func[data_pos] <- inf_func_no_na
  return(list(att = att, inf_func = inf_func, logit_coef = logit_coef))
}
