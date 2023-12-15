estimate_did <- function(dt_did, covnames, control_type,
                         last_coef = NULL, cache_ps_fit, cache_hess){

  # preprocess --------

  oldn <- dt_did[, .N]
  data_pos <-  which(dt_did[, !is.na(D)])
  dt_did <- dt_did[data_pos]
  n <- dt_did[, .N]
  if(is.null(covnames)){
    covvars <- NULL
  } else {
    covvars <- as.matrix(dt_did[,.SD, .SDcols = covnames])
  }

  ipw <- control_type %in% c("ipw", "dr") & !is.null(covvars)
  or <- control_type %in% c("reg", "dr") & !is.null(covvars) #OR is REG

  # ipw --------

  if(ipw){
    if(is.null(cache_ps_fit)|is.null(cache_hess)){ #if no cache, calcuate ipw

      #estimate the logit
      prop_score_est <- suppressWarnings(parglm.fit(covvars, dt_did[, D],
                                                    family = stats::binomial(), start = last_coef,
                                                    weights = dt_did[, weights],
                                                    control = parglm.control(nthreads = getDTthreads()),
                                                    intercept = FALSE))
      class(prop_score_est) <- "glm" #trick the vcov function to think that this is a glm object to dispatch the write method
      #const is implicitly put into the ipw formula, need to incorporate it manually
      hess <- stats::vcov(prop_score_est) * n #for the influence function

      logit_coef <-  prop_score_est$coefficients
      logit_coef[is.na(logit_coef)|abs(logit_coef) > 1e10] <- 0 #put extreme value and na to 0
      prop_score_fit <- fitted(prop_score_est)
      if(max(prop_score_fit) >= 1){warning(paste0("support overlap condition violated for some group_time"))}
      prop_score_fit <- pmin(1-1e-16, prop_score_fit) #for the ipw


    } else { #when using multiple outcome, ipw cache can be reused
      hess <- cache_hess
      prop_score_fit <- cache_ps_fit
      logit_coef <- NULL #won't be needing the approximate cache
    }

    #get the results into the main did dt
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
  inf_cont <- (inf_cont_did+inf_cont_ipw+inf_cont_or)/dt_did[, mean(cont_ipw_weight)]
  inf_treat <- (inf_treat_did+inf_treat_or)/dt_did[,mean(treat_ipw_weight)]
  inf_func_no_na <- inf_treat - inf_cont

  #post process (fill zeros for irrelevant ones)
  inf_func <- rep(0, oldn) #the default needs to be 0 for the matrix multiplication
  inf_func_no_na <- inf_func_no_na * oldn / n #adjust the value such that mean over the whole id size give the right result
  inf_func[data_pos] <- inf_func_no_na

  return(list(att = att, inf_func = inf_func, logit_coef = logit_coef, #for next gt
              cache_ps_fit = prop_score_fit, cache_hess = hess)) #for next outcome
}
