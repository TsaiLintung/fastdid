estimate_did <- function(dt_did, covvars, p, cache){
  
  #estimate did
  param <- as.list(environment())
  if(!p$allow_unbalance_panel){
    result <- do.call(estimate_did_bp, param)
  } else {
    result <- do.call(estimate_did_rc, param)
  }
  return(result)
}

estimate_did_bp <- function(dt_did, covvars, p, cache){
  
  # preprocess --------
  oldn <- dt_did[, .N]
  data_pos <-  which(dt_did[, !is.na(D)])
  dt_did <- dt_did[data_pos]
  n <- dt_did[, .N]
  
  if(is.matrix(covvars)){
    ipw <- p$control_type %in% c("ipw", "dr") 
    or <- p$control_type %in% c("reg", "dr")
    covvars <- covvars[data_pos,] 
  } else {
    ipw <- FALSE
    or <- FALSE
  }
  
  # ipw --------

  if(ipw){
    if(is.null(cache)){ #if no cache, calcuate ipw
      #estimate the logit
      prop_score_est <- suppressWarnings(parglm::parglm.fit(covvars, dt_did[, D],
                                                    family = stats::binomial(), 
                                                    weights = dt_did[, weights],
                                                    control = parglm.control(nthreads = ifelse(p$parallel, 1, getDTthreads())), #no parallel if already parallel
                                                    intercept = FALSE))
      class(prop_score_est) <- "glm" #trick the vcov function to think that this is a glm object to dispatch the write method
      #const is implicitly put into the ipw formula, need to incorporate it manually

      logit_coef <-  prop_score_est$coefficients
      
      if(anyNA(logit_coef)){
        warning("some propensity score estimation resulted in NA coefficients, likely cause by perfect colinearity")
      }
      
      logit_coef[is.na(logit_coef)|abs(logit_coef) > 1e10] <- 0 #put extreme value and na to 0
      prop_score_fit <- fitted(prop_score_est)
      if(max(prop_score_fit) >= 1-1e-10){warning(paste0("extreme propensity score: ", max(prop_score_fit), ", support overlap is likely to be violated"))} #<=0 (only in control) is fine for ATT since it is just not used 
      prop_score_fit <- pmin(1-1e-10, prop_score_fit) #for the ipw

      hess <- stats::vcov(prop_score_est) * n #for the influence function
      hess[is.na(hess)|abs(hess) > 1e10] <- 0
      
      cache <- list(hess = hess, ps = prop_score_fit)

    } else { #when using multiple outcome, ipw cache can be reused
      hess <- cache$hess
      prop_score_fit <- cache$ps
      logit_coef <- NA #won't be needing the approximate cache
    }

    #get the results into the main did dt
    dt_did[, ps := prop_score_fit]
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*ps*(1-D)/(1-ps)]

  } else {

    prop_score_fit <- rep(1,n)
    logit_coef <- NA
    hess <- NA

    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*(1-D)]

  }

  #delta y is needed for PR
  dt_did[, delta_y := post.y-pre.y]

  # or --------

  if(or){

    #TODO: this should be optimized with better backend and some caching

    #should change to speedlm or something

    control_bool <- dt_did[, D==0]
    reg_coef <- stats::coef(stats::lm.wfit(x = covvars[control_bool,], y = dt_did[control_bool,delta_y],
                                           w = dt_did[control_bool,weights]))

    if(anyNA(reg_coef)){
      stop("some outcome regression resulted in NA coefficients, likely cause by perfect colinearity")
    }

    #the control function from outcome regression
    dt_did[, or_delta := as.vector(tcrossprod(reg_coef, covvars))]

  } else {
    dt_did[, or_delta := 0]
  }

  #did --------

  dt_did[, att_treat := treat_ipw_weight*(delta_y-or_delta)] #minus the OR adjust
  dt_did[, att_cont := cont_ipw_weight*(delta_y-or_delta)]

  weighted_treat_delta <- dt_did[,sum(att_treat)/sum(treat_ipw_weight)]
  weighted_cont_delta <- dt_did[,sum(att_cont)/sum(cont_ipw_weight)]

  att <- weighted_treat_delta - weighted_cont_delta
  
  # influence --------

  # influence from ipw
  if(ipw){

    M2 <- colMeans(dt_did[, cont_ipw_weight*(delta_y-weighted_cont_delta-or_delta)] * covvars)

    score_ps <- dt_did[, weights*(D-ps)] * covvars
    
    asym_linear_ps <- score_ps %*% hess

    #ipw for control
    inf_cont_ipw <- asym_linear_ps %*% as.matrix(M2)

  } else {inf_cont_ipw <- 0}

  if(or){

    M1 <- colMeans(dt_did[, treat_ipw_weight] * covvars)
    M3 <- colMeans(dt_did[, cont_ipw_weight] * covvars)

    or_x <- dt_did[, weights*(1-D)] * covvars
    or_ex <- dt_did[, weights*(1-D)*(delta_y - or_delta)] * covvars
    XpX <- crossprod(or_x, covvars)/n

    #calculate alrw = eX (XpX)^-1 by solve XpX*alrw = ex, much faster since avoided inv
    asym_linear_or <- t(solve(XpX, t(or_ex)))# |> reverse_col()

    #or for treat
    inf_treat_or <- -asym_linear_or %*% M1 #a negative sign here, since or_delta is subtracted from the att

    #or for control
    inf_cont_or <- -asym_linear_or %*% M3

  } else {
    inf_treat_or <- 0
    inf_cont_or <- 0
  }

  # influence from did
  inf_cont_did <- dt_did[, att_cont - cont_ipw_weight*weighted_cont_delta]
  inf_treat_did <-  dt_did[, (att_treat - treat_ipw_weight*weighted_treat_delta)]


  #get overall influence function
  #if(dt_did[, mean(cont_ipw_weight)] < 1e-10){warning("little/no overlap in covariates between control and treat group, estimates are unstable.")}
  inf_cont <- (inf_cont_did+inf_cont_ipw+inf_cont_or)/dt_did[, mean(cont_ipw_weight)]
  inf_treat <- (inf_treat_did+inf_treat_or)/dt_did[,mean(treat_ipw_weight)]
  inf_func_no_na <- inf_treat - inf_cont

  #post process (fill zeros for irrelevant ones)
  inf_func <- rep(0, oldn) #the default needs to be 0 for the matrix multiplication
  inf_func_no_na <- inf_func_no_na * oldn / n #adjust the value such that mean over the whole id size give the right result
  inf_func[data_pos] <- inf_func_no_na
  
  return(list(att = att, inf_func = inf_func, cache = list(ps = prop_score_fit, hess = hess))) #for next outcome
}

estimate_did_rc <- function(dt_did, covvars, p, cache){
  
  #TODO: skip if not enough valid data
  
  # preprocess --------
  
  
  oldn <- dt_did[, .N]
  data_pos <-  which(dt_did[, !is.na(D)])
  dt_did <- dt_did[data_pos]
  n <- dt_did[, .N]
  
  #separate the dataset into pre and post
  dt_did[, inpre := as.numeric(!is.na(pre.y))]
  dt_did[, inpost := as.numeric(!is.na(post.y))]
  n_pre <- dt_did[, sum(!is.na(pre.y))]
  n_post <- dt_did[, sum(!is.na(post.y))]
  
  sum_weight_pre <- dt_did[, sum(inpre*weights)]
  sum_weight_post <- dt_did[, sum(inpost*weights)]
  
  if(is.matrix(covvars)){
    ipw <- p$control_type %in% c("ipw", "dr") 
    or <- p$control_type %in% c("reg", "dr")
    covvars <- covvars[data_pos,] 
  } else {
    ipw <- FALSE
    or <- FALSE
  }
  
  # ipw --------
  
  if(ipw){
    
    #no caching since composition changes by period
    
    #estimate the logit
    prop_score_est <- suppressWarnings(parglm.fit(covvars, dt_did[, D],
                                                  family = stats::binomial(),
                                                  weights = dt_did[, weights*(inpre+inpost)*n/(n_pre+n_post)], #when seen in both pre and post have double weight
                                                  control = parglm.control(nthreads = ifelse(p$parallel, 1, getDTthreads())),
                                                  intercept = FALSE)) #*(inpre+inpost)
    class(prop_score_est) <- "glm" #trick the vcov function to think that this is a glm object to dispatch the write method
    #const is implicitly put into the ipw formula, need to incorporate it manually
    
    #for the influence, will be cached
    hess <- stats::vcov(prop_score_est) * n #for the influence function
    
    logit_coef <-  prop_score_est$coefficients 
    logit_coef[is.na(logit_coef)|abs(logit_coef) > 1e10] <- 0 #put extreme value and na to 0
    prop_score_fit <- fitted(prop_score_est)
    if(max(prop_score_fit) >= 1){warning(paste0("support overlap condition violated for some group_time"))}
    prop_score_fit <- pmin(1-1e-16, prop_score_fit) #for the ipw
    
    #get the results into the main did dt
    dt_did[, ps := prop_score_fit]
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*ps*(1-D)/(1-ps)]
    
  } else {
    
    prop_score_fit <- rep(1,n)
    logit_coef <- NA
    hess <- NA
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*(1-D)]
    
  }
  
  # or --------
  
  if(or){
    
    control_bool_post <- dt_did[, D==0 & inpost] #control group and have obs in post period
    control_bool_pre <- dt_did[, D==0 & inpre]
    reg_coef_post <- stats::coef(stats::lm.wfit(x = covvars[control_bool_post,], y = dt_did[control_bool_post,post.y],
                                                w = dt_did[control_bool_post,weights]))

    reg_coef_pre <- stats::coef(stats::lm.wfit(x = covvars[control_bool_pre,], y = dt_did[control_bool_pre,pre.y],
                                               w = dt_did[control_bool_pre,weights]))

    if(anyNA(reg_coef_post)|anyNA(reg_coef_pre)){
      stop("some outcome regression resulted in NA coefficients, likely cause by perfect colinearity")
    }

    #the control function from outcome regression
    dt_did[, or_delta_post := as.vector(tcrossprod(reg_coef_post, covvars))]
    dt_did[, or_delta_pre := as.vector(tcrossprod(reg_coef_pre, covvars))]
    
  } else {
    dt_did[, or_delta_post := 0]
    dt_did[, or_delta_pre := 0]
  }
  
  #did --------
  
  #mean weight
  mean_wcpo <-  dt_did[,sum(cont_ipw_weight*inpost)/n_post]
  mean_wtpo <-  dt_did[,sum(treat_ipw_weight*inpost)/n_post]
  mean_wcpr <-  dt_did[,sum(cont_ipw_weight*inpre)/n_pre]
  mean_wtpr <-  dt_did[,sum(treat_ipw_weight*inpre)/n_pre]
  
  #delta y is needed for PR
  dt_did[, att_treat_post := treat_ipw_weight*(post.y-or_delta_post)/mean_wtpo] #minus the OR adjust
  dt_did[, att_cont_post :=  cont_ipw_weight*(post.y-or_delta_post)/mean_wcpo]
  dt_did[, att_treat_pre := treat_ipw_weight*(pre.y-or_delta_pre)/mean_wtpr] #minus the OR adjust
  dt_did[, att_cont_pre := cont_ipw_weight*(pre.y-or_delta_pre)/mean_wcpr]
  
  weighted_treat_post <- dt_did[,mean(att_treat_post, na.rm = TRUE)]
  weighted_cont_post <- dt_did[,mean(att_cont_post, na.rm = TRUE)]
  weighted_treat_pre <- dt_did[,mean(att_treat_pre, na.rm = TRUE)]
  weighted_cont_pre <- dt_did[,mean(att_cont_pre, na.rm = TRUE)]
  
  att <- (weighted_treat_post - weighted_treat_pre) - (weighted_cont_post - weighted_cont_pre)
  
  # influence --------
  
  # influence from ipw
  if(ipw){
    
    # a bit unsure about this part
    M2_post <- colSums(dt_did[, inpost*cont_ipw_weight*(post.y-weighted_cont_post-or_delta_post)/n] * covvars, na.rm = TRUE) / mean_wcpo
    M2_pre <- colSums(dt_did[, inpre*cont_ipw_weight*(pre.y-weighted_cont_pre-or_delta_pre)/n] * covvars, na.rm = TRUE) / mean_wcpr
    
    #not sure about /2
    score_ps <- dt_did[, weights*(inpre+inpost)*n/(n_pre+n_post)*(D-ps)] * covvars#weight is doubled for observed in both post and pre
    asym_linear_ps <- score_ps %*% hess 
    
    #ipw for control
    inf_cont_ipw_post <- asym_linear_ps %*% M2_post 
    inf_cont_ipw_pre  <- asym_linear_ps %*% M2_pre
    
    
  } else {
    inf_cont_ipw_post <- 0
    inf_cont_ipw_pre <- 0
  }
  
  
  if(or){

    M1_post <- colSums(dt_did[, inpost*treat_ipw_weight/n] * covvars, na.rm = TRUE) / mean_wtpo
    M1_pre <- colSums(dt_did[, inpre*treat_ipw_weight/n] * covvars, na.rm = TRUE) / mean_wtpr
    M3_post <- colSums(dt_did[, inpost*cont_ipw_weight/n] * covvars, na.rm = TRUE) / mean_wcpo
    M3_pre <- colSums(dt_did[, inpre*cont_ipw_weight/n] * covvars, na.rm = TRUE) / mean_wcpr

    or_x_post <- dt_did[, inpost*weights*(1-D)] * covvars
    or_x_pre <- dt_did[, inpre*weights*(1-D)] * covvars
    or_ex_post <- dt_did[, inpost*weights*(1-D)*(post.y - or_delta_post)] * covvars
    or_ex_pre <- dt_did[, inpre*weights*(1-D)*(pre.y - or_delta_pre)] * covvars
    XpX_post <- crossprod(or_x_post, covvars)/n_post
    XpX_pre <- crossprod(or_x_pre, covvars)/n_pre

    #calculate alrw = eX (XpX)^-1 by solve XpX*alrw = ex, much faster since avoided inv
    asym_linear_or_post <- t(solve(XpX_post, t(or_ex_post)))
    asym_linear_or_pre <- t(solve(XpX_pre, t(or_ex_pre)))

    #or for treat
    inf_treat_or_post <- -asym_linear_or_post %*% M1_post #a negative sign here, since or_delta is subtracted from the att, THE PROBLEM
    inf_treat_or_pre <- -asym_linear_or_pre %*% M1_pre

    #or for control
    inf_cont_or_post <- -asym_linear_or_post %*% M3_post
    inf_cont_or_pre <- -asym_linear_or_pre %*% M3_pre
    
  } else {
    inf_treat_or_post <- 0
    inf_treat_or_pre <- 0
    inf_cont_or_post <- 0
    inf_cont_or_pre <- 0
  }
  
  # influence from did
  inf_cont_did_post <- dt_did[, att_cont_post - cont_ipw_weight*inpost*weighted_cont_post/mean_wcpo]
  inf_treat_did_post <-  dt_did[, att_treat_post - treat_ipw_weight*inpost*weighted_treat_post/mean_wtpo]
  inf_cont_did_pre <- dt_did[, att_cont_pre - cont_ipw_weight*inpre*weighted_cont_pre/mean_wcpr]
  inf_treat_did_pre <-  dt_did[, att_treat_pre -  treat_ipw_weight*inpre*weighted_treat_pre/mean_wtpr]
  
  #fill zero to avoid NA from addition
  inf_cont_did_post[is.na(inf_cont_did_post)] <- 0
  inf_treat_did_post[is.na(inf_treat_did_post)] <- 0
  inf_cont_did_pre[is.na(inf_cont_did_pre)] <- 0
  inf_treat_did_pre[is.na(inf_treat_did_pre)] <- 0
  
  #get overall influence function
  inf_cont_post <- inf_cont_did_post+inf_cont_ipw_post+inf_cont_or_post
  inf_treat_post <- inf_treat_did_post+inf_treat_or_post
  inf_cont_pre <- inf_cont_did_pre+inf_cont_ipw_pre+inf_cont_or_pre
  inf_treat_pre <- inf_treat_did_pre+inf_treat_or_pre
  
  #post process
  inf_func_no_na_post <- (inf_treat_post - inf_cont_post) * oldn / n_post #adjust the value such that mean over the whole id size give the right result
  inf_func_no_na_post[is.na(inf_func_no_na_post)] <- 0 #fill 0 for NA part (no influce if not in this gt)
  
  inf_func_no_na_pre <- (inf_treat_pre - inf_cont_pre) * oldn / n_pre #adjust the value such that mean over the whole id size give the right result
  inf_func_no_na_pre[is.na(inf_func_no_na_pre)] <- 0
  
  inf_func <- rep(0, oldn) #the default needs to be 0 for the matrix multiplication
  inf_func[data_pos] <- inf_func_no_na_post - inf_func_no_na_pre

  return(list(att = att, inf_func = inf_func, cache = list(ps = prop_score_fit, hess = hess))) #for next outcome
}

# utilities ------

reverse_col <- function(x){
  return(x[,ncol(x):1])
}

