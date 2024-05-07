#2024-05-07
message('loading fastdid source ver. ver: 0.9.3 (timing), date: 2024-05-07')
require(data.table);
 require(stringr);
 require(BMisc);
 require(collapse);
 require(dreamerr);
 require(parglm);
aggregate_gt <- function(gt_result, aux, p){
  

  #release the stuff
  id_cohorts <- aux$dt_inv[, G]
  
  id_dt <- data.table(weight = aux$weights/sum(aux$weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group_time <- gt_result$gt |> merge(pg_dt, by = "G")
  
  setorder(group_time, time, G) #change the order to match the order in gtatt
  
  gt_result$inf_func <- as.matrix(gt_result$inf_func)
  
  if(p$result_type == "group_time"){
    
    #don't need to do anything
    targets <- group_time[, unique(G*max(time)+time)]
    inf_matrix <- gt_result$inf_func
    agg_att <- as.vector(gt_result$att)
    
  } else {
    
    #get which gt(s) is a part of the aggregated param
    agg_sch <- get_aggregate_scheme(group_time, p$result_type, aux$weights, id_cohorts, p$balanced_event_time)
    targets <- agg_sch$targets
    
    #aggregated att
    agg_att <- agg_sch$agg_weights %*% gt_result$att
    
    #get the influence from weight estimation
    #this needs to be optimized!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    inf_weights <- sapply(asplit(agg_sch$agg_weights, 1), function (x){
      get_weight_influence(x, gt_result$att, aux$weights, id_cohorts, group_time[, .(G, time)])
    })
    
    #aggregated influence function
    inf_matrix <- (gt_result$inf_func %*% t(agg_sch$agg_weights)) + inf_weights 
    
  }
  
  #get se from influence function
  agg_se <- get_se(inf_matrix, p$boot, p$biters, aux$cluster, p$clustervar)
  
  # post process
  result <- data.table(targets, agg_att, agg_se)
  names(result) <- c("target", "att", "se")
  result[,outcome := gt_result$outname]
  
  return(result)
}

get_aggregate_scheme <- function(group_time, result_type, weights, id_cohorts, balanced_event_time){
  
  #browser()
  
  agg_weights <- data.table()
  gt_count <- group_time[, .N]
  
  bool_to_pn <- function(x){ifelse(x, 1, -1)}
  
  #choose the target based on aggregation type
  if (result_type == "dynamic") {
    group_time[, target := time-G]
  } else if (result_type == "group") {
    group_time[, target := G*(bool_to_pn(time>=G))] # group * treated
  } else if (result_type == "time") {
    group_time[, target := time*(bool_to_pn(time>=G))] #calendar time * treated
  } else if (result_type == "simple") {
    group_time[, target := bool_to_pn(time>=G)] #treated / not treated
  }  
  
  targets <- group_time[, unique(target)]
  targets <- sort(as.integer(targets))
  
  #for balanced cohort composition in dynamic setting
  #a cohort us only used if it is seen for all dynamic time
  if(result_type == "dynamic" & !is.na(balanced_event_time)){

    cohorts <- group_time[, .(max_et = max(time-G),
                              min_et = min(time-G)), by = "G"]
    cohorts[, used := max_et >= balanced_event_time] #the max
    if(!cohorts[, any(used)]){stop("balanced_comp_range outside avalible range")}
    group_time[, used := G %in% cohorts[used == TRUE, G]]
    
    targets <- targets[targets <= balanced_event_time & targets >= cohorts[used == TRUE, min(min_et)]]
    
  } else{group_time[, used := TRUE]}
  
  for(tar in targets){ #the order matters
    
    group_time[, targeted := target == tar & used]
    
    total_pg <- group_time[targeted == TRUE, sum(pg)] #all gt that fits in the target
    group_time[, weight := ifelse(targeted, pg/total_pg, 0)] #weight is 0 if not a target
    target_weights <- group_time[, .(weight)] |> transpose()
    
    group_time[, targeted := NULL]
    
    agg_weights <- rbind(agg_weights, target_weights)
  }
  
  return(list(agg_weights = as.matrix(agg_weights), #a matrix of each target and gt's weight in it 
              targets = targets)) 
}

get_weight_influence <- function(agg_weights, gt_att, weights, id_cohorts, group) {

  keepers <- which(agg_weights > 0)
  
  id_dt <- data.table(weight = weights/sum(weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group <- group |> merge(pg_dt, by = "G")
  
  group[, time := as.integer(time)]
  group[, G := as.integer(G)]
  setorder(group, time, G)

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (weights*BMisc::TorF(id_cohorts == group[k,G]) - group[k,pg]) /
      sum(group[keepers,pg])
  })
  # effect of estimating weights in the denominator
  if2 <- base::rowSums(sapply(keepers, function(k) {
    weights*BMisc::TorF(id_cohorts == group[k,G]) - group[k,pg]
  })) %*%
    t(group[keepers,pg]/(sum(group[keepers,pg])^2))
  # return the influence function for the weights
  inf_weight <- (if1 - if2) %*% as.vector(gt_att[keepers])
  inf_weight[abs(inf_weight) < sqrt(.Machine$double.eps)*10] <- 0 #fill zero 
  return(inf_weight)
}

get_se <- function(inf_matrix, boot, biters, cluster, clustervar) {
  
  if(boot){
    
    top_quant <- 0.75
    bot_quant <- 0.25
    if(!allNA(clustervar)){
      #take average within the cluster
      cluster_n <- stats::aggregate(cluster, by=list(cluster), length)[,2]
      inf_matrix <- fsum(inf_matrix, cluster) / cluster_n #the mean without 0 for each cluster of each setting
    }
    
    boot_results <- BMisc::multiplier_bootstrap(inf_matrix, biters = biters) %>% as.data.table()
    
    boot_top <- boot_results[, lapply(.SD, function(x) stats::quantile(x, top_quant, type=1, na.rm = TRUE))]
    boot_bot <- boot_results[, lapply(.SD, function(x) stats::quantile(x, bot_quant, type=1, na.rm = TRUE))]
    
    dt_se <- rbind(boot_bot, boot_top) %>% transpose()
    names(dt_se) <- c("boot_bot", "boot_top")
    
    se <- dt_se[,(boot_top-boot_bot)/(qnorm(top_quant) - qnorm(bot_quant))]
    se[se < sqrt(.Machine$double.eps)*10] <- NA
    
  } else {

    inf_matrix <- inf_matrix  |> as.data.table()
    se <- inf_matrix[, lapply(.SD, function(x) sqrt(sum(x^2, na.rm = TRUE)/length(x)^2))] %>% as.vector() #should maybe use n-1 but did use n
    
  }
  return(unlist(se))
}


estimate_did <- function(dt_did, covvars, control_type,
                         last_coef = NULL, cache_ps_fit, cache_hess){
  
  # preprocess --------
  oldn <- dt_did[, .N]
  data_pos <-  which(dt_did[, !is.na(D)])
  dt_did <- dt_did[data_pos]
  n <- dt_did[, .N]
  
  if(is.matrix(covvars)){
    ipw <- control_type %in% c("ipw", "dr") 
    or <- control_type %in% c("reg", "dr")
    covvars <- covvars[data_pos,] 
  } else {
    ipw <- FALSE
    or <- FALSE
  }
  
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

      logit_coef <-  prop_score_est$coefficients
      
      if(anyNA(logit_coef)){
        warning("some propensity score estimation resulted in NA coefficients, likely cause by perfect colinearity")
      }
      
      logit_coef[is.na(logit_coef)|abs(logit_coef) > 1e10] <- 0 #put extreme value and na to 0
      prop_score_fit <- fitted(prop_score_est)
      if(max(prop_score_fit) >= 1){warning(paste0("support overlap condition violated for some group_time"))}
      prop_score_fit <- pmin(1-1e-16, prop_score_fit) #for the ipw

      hess <- stats::vcov(prop_score_est) * n #for the influence function
      hess[is.na(hess)|abs(hess) > 1e10] <- 0

    } else { #when using multiple outcome, ipw cache can be reused
      hess <- cache_hess
      prop_score_fit <- cache_ps_fit
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

  return(list(att = att, inf_func = inf_func, logit_coef = logit_coef, #for next gt
              cache_ps_fit = prop_score_fit, cache_hess = hess)) #for next outcome
}

estimate_did_rc <- function(dt_did, covvars, control_type,
                         last_coef = NULL, cache_ps_fit, cache_hess){
  
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
    ipw <- control_type %in% c("ipw", "dr") 
    or <- control_type %in% c("reg", "dr")
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
                                                  family = stats::binomial(), start = last_coef,
                                                  weights = dt_did[, weights*(inpre+inpost)*n/(n_pre+n_post)], #when seen in both pre and post have double weight
                                                  control = parglm.control(nthreads = getDTthreads()),
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
    
    stop("or in RC should not be called right now")

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
    stop("or in RC should not be called right now")
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
    inf_treat_or_post <- -asym_linear_or_post %*% M1_post #a negative sign here, since or_delta is subtracted from the att
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

  return(list(att = att, inf_func = inf_func, logit_coef = logit_coef, #for next gt
              cache_ps_fit = prop_score_fit, cache_hess = hess)) #for next outcome
}

estimate_gtatt <- function(aux, p) {

  #chcek if there is availble never-treated group
  if(!Inf %in% aux$cohorts & p$control_option != "notyet"){
    warning("no never-treated availble, switching to not-yet-treated control")
    p$control_option <- "notyet"
  }
  
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
        if(p$base_period == "universal"){
          base_period <- g-1-p$anticipation
        } else {
          base_period <- ifelse(t>=g, g-1-p$anticipation, t-1)
        }
        
        did_setup <- get_did_setup(g, t, base_period, aux, p)
        if(is.null(did_setup)){next} #no gtatt if no did setup
        
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

get_did_setup <- function(g, t, base_period, aux, p){
  
  treated_cohorts <- aux$cohorts[!is.infinite(aux$cohorts)]
  
  #setup and checks
  
  if(p$control_option == "never"){
    min_control_cohort <- Inf
  } else {
    if(!is.infinite(p$min_control_cohort_diff)){
      min_control_cohort <- max(g+p$min_control_cohort_diff, min(treated_cohorts))
    } else {
      min_control_cohort <- max(t, base_period)+p$anticipation+1
    }
  }

  if(!is.infinite(p$max_control_cohort_diff)){
    max_control_cohort <- min(g+p$max_control_cohort_diff, max(treated_cohorts))
  } else {
    max_control_cohort <- ifelse(p$control_option == "notyet", max(treated_cohorts), Inf) 
  }
  
  if(t == base_period | #no treatment effect for the base period
     base_period < min(aux$time_periods) | #no treatment effect for the first period, since base period is not observed
     g >= max_control_cohort | #no treatment effect for never treated or the last treated cohort (for not yet notyet)
     t >= max_control_cohort){ #no control available if the last cohort is treated too
    return(NULL)
  } else {
    #select the control and treated cohorts
    did_setup <- rep(NA, aux$id_size)
    did_setup[get_cohort_pos(aux$cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
    did_setup[get_cohort_pos(aux$cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
    
    if(!is.na(p$filtervar)){
      did_setup[!aux$filters[[base_period]]] <- NA #only use units with filter == TRUE at base period
    }
    
    return(did_setup)
  }
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

#' Fast Staggered DID Estimation
#'
#' Performs Difference-in-Differences (DID) estimation fast.
#'
#' @param data A data table containing the panel data.
#' @param timevar The name of the time variable.
#' @param cohortvar The name of the cohort (group) variable.
#' @param unitvar The name of the unit (id) variable.
#' @param outcomevar The name of the outcome variable.
#' @param control_option The control units used for the DiD estimates. Default is "both".
#' @param result_type A character string indicating the type of result to be returned. Default is "group_time".
#' @param balanced_event_time A numeric scalar that indicates the max event time to balance the cohort composition, only meaningful when result_type == "dynamic". Default is NA
#' @param control_type The method for controlling for covariates. "ipw" for inverse probability weighting, "reg" for outcome regression, or "dr" for doubly-robust
#' @param allow_unbalance_panel Whether allow unbalance panel as input (if false will coerce the dataset to a balanced panel). Default is FALSE 
#' @param boot Logical, indicating whether bootstrapping should be performed. Default is FALSE.
#' @param biters The number of bootstrap iterations. Only relevant if boot = TRUE. Default is 1000.
#' @param weightvar The name of the weight variable (optional).
#' @param clustervar The name of the cluster variable, can only be used when boot == TRUE (optional).
#' @param covariatesvar A character vector containing the names of time-invariant covariate variables (optional).
#' @param varycovariatesvar A character vector containing the names of time-varying covariate variables (optional).
#' @param filtervar A logical vector that specifies which units to use (can be time varying, units will only be included in 2x2 when filtervar is TRUE is the base period)
#' @param copy whether to copy the dataset before processing, set to false to speed up the process, but the input data will be altered.
#' @param validate whether to validate the dataset before processing.
#' @param anticipation periods with aniticipation (delta in CS, default is 0, reference period is g - delta - 1).
#' @param min_control_cohort_diff the min cohort difference between treated and control group 
#' @param max_control_cohort_diff the max cohort difference between treated and control group 
#' @param base_period same as did
#' 
#' @import data.table parglm stringr dreamerr BMisc 
#' @importFrom stats quantile vcov sd binomial fitted qnorm rnorm as.formula
#' @importFrom collapse allNA fnrow whichNA fnunique fsum na_insert
#' @return A data.table containing the estimated treatment effects and standard errors.
#' @export
#'
#' @examples
#' # simulated data
#' simdt <- sim_did(1e+03, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE)
#' dt <- simdt$dt
#' 
#' #basic call
#' result <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                   unitvar = "unit", outcomevar = "y",  
#'                   result_type = "group_time")
#' 
#' #control for covariates
#' result2 <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                    unitvar = "unit", outcomevar = "y",  
#'                    result_type = "group_time",
#'                    covariatesvar = c("x", "x2"))
#'                   
#' #bootstrap and clustering
#' result3 <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                    unitvar = "unit", outcomevar = "y",  
#'                    result_type = "group_time",
#'                    boot = TRUE, clustervar = "x")
#'
#' #estimate for multiple outcomes
#' result4 <- fastdid(data = dt, #the dataset
#'                    timevar = "time", cohortvar = "G", unitvar = "unit", 
#'                    outcomevar = c("y", "y2"), #name of the outcome columns
#'                    result_type = "group_time") 
#'
#' @keywords difference-in-differences fast computation panel data estimation did
fastdid <- function(data,
                    timevar, cohortvar, unitvar, outcomevar, 
                    control_option="both",result_type="group_time", balanced_event_time = NA,
                    control_type = "ipw", allow_unbalance_panel = FALSE, boot=FALSE, biters = 1000,
                    weightvar=NA,clustervar=NA, covariatesvar = NA, varycovariatesvar = NA, filtervar = NA,
                    copy = TRUE, validate = TRUE,
                    max_control_cohort_diff = Inf, anticipation = 0, min_control_cohort_diff = -Inf, base_period = "universal"
                    ){
  
  # validate arguments --------------------------------------------------------

  if(!is.data.table(data)){
    warning("coercing input into a data.table.")
    data <- as.data.table(data)
  } 
  if(copy){dt <- copy(data)} else {dt <- data}
  
  dt_names <- names(dt)
  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_set_arg(timevar, unitvar, cohortvar, "match", .choices = dt_names, .message = name_message)
  
  covariate_message <- "__ARG__ must be NA or a character vector which are all names of columns from the dataset."
  check_set_arg(varycovariatesvar, covariatesvar, outcomevar, 
            "NA | multi match", .choices = dt_names, .message = covariate_message)
  
  checkvar_message <- "__ARG__ must be NA or a character scalar if a name of columns from the dataset."
  check_set_arg(weightvar, clustervar, filtervar,
            "NA | match", .choices = dt_names, .message = checkvar_message)
  
  check_set_arg(control_option, "match", .choices = c("both", "never", "notyet")) #kinda bad names since did's notyet include both notyet and never
  check_set_arg(control_type, "match", .choices = c("ipw", "reg", "dr")) 
  check_set_arg(base_period, "match", .choices = c("varying", "universal"))
  check_arg(copy, validate, boot, allow_unbalance_panel, "scalar logical")
  check_arg(max_control_cohort_diff, min_control_cohort_diff, anticipation, "scalar numeric")
  
  if(!is.na(balanced_event_time)){
    if(result_type != "dynamic"){stop("balanced_event_time is only meaningful with result_type == 'dynamic'")}
    check_arg(balanced_event_time, "numeric scalar")
  }
  if(allow_unbalance_panel == TRUE & control_type %in% c("dr", "reg")){
    stop("fastdid currently only supprts ipw when allowing for unbalanced panels.")
  }
  if(allow_unbalance_panel == TRUE & !allNA(varycovariatesvar)){
    stop("fastdid currently only supprts time varying covariates when allowing for unbalanced panels.")
  }
  if(any(covariatesvar %in% varycovariatesvar) & !allNA(varycovariatesvar) & !allNA(covariatesvar)){
    stop("time-varying var and invariant var have overlaps.")
  }
  if(!boot & !allNA(clustervar)){
    stop("clustering only available with bootstrap")
  }
  
  # coerce non-sensible option
  
  if((!is.infinite(max_control_cohort_diff) | !is.infinite(min_control_cohort_diff)) & control_option == "never"){
    warning("control_cohort_diff can only be used with not yet")
    p$control_option <- "notyet"
  }
  
  p <- list(timevar = timevar,
            cohortvar = cohortvar,
            unitvar = unitvar,
            outcomevar =  outcomevar,
            weightvar = weightvar,
            clustervar = clustervar,
            filtervar = filtervar,
            covariatesvar = covariatesvar,
            varycovariatesvar = varycovariatesvar,
            control_option = control_option,
            result_type = result_type,
            balanced_event_time = balanced_event_time,
            control_type = control_type,
            allow_unbalance_panel = allow_unbalance_panel,
            boot = boot, 
            biters = biters,
            max_control_cohort_diff = max_control_cohort_diff,
            min_control_cohort_diff = min_control_cohort_diff,
            anticipation = anticipation, 
            base_period = base_period)
  
  
  # validate data -----------------------------------------------------
  
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))
  
  if(validate){
    varnames <- c("time", "G", "unit",outcomevar,weightvar,clustervar,covariatesvar,varycovariatesvar,filtervar)
    dt <- validate_did(dt, varnames, p)
  }
  
  # preprocess -----------------------------------------------------------
  
  #make dt conform to the WLOG assumptions of fastdid
  coerce_result <- coerce_dt(dt, p) #also changed dt
  dt <- coerce_result$dt
  # get auxiliary data
  aux <- get_auxdata(dt, p)
  
  # main part  -------------------------------------------------

  # estimate attgt
  gt_result_list <- estimate_gtatt(aux, p)
  
  # aggregate the result by outcome
  all_result <- rbindlist(lapply(gt_result_list, function(x){
    result <- aggregate_gt(x, aux, p)
    #convert "targets" back to meaningful parameter identifiers like cohort 1 post, time 2 post 
    result <- result |> convert_targets(result_type, coerce_result$time_change)
    return(result)
  }))

  return(all_result)
  
}

# small steps ----------------------------------------------------------------------

coerce_dt <- function(dt, p){
  
  #change to int before sorting
  if(!is.numeric(dt[, G])){
    dt[, G := as.numeric(G)]
  }
  if(!is.numeric(dt[, time])){
    dt[, time := as.numeric(time)] 
  }
  
  if(p$allow_unbalance_panel){
    dt_inv_raw <- dt[dt[, .I[1], by = unit]$V1]
    setorder(dt_inv_raw, G)
    dt_inv_raw[, new_unit := 1:.N] #let unit start from 1 .... N, useful for knowing which unit is missing
    dt <- dt |> merge(dt_inv_raw[,.(unit, new_unit)], by = "unit")
    dt[, unit := new_unit]
  }
  
  setorder(dt, time, G, unit) #sort the dataset essential for the sort-once-quick-access 
  #deal with time, coerice time to 1,2,3,4,5.......
  time_periods <- dt[, unique(time)]
  time_size <- length(time_periods)
  
  time_offset <- min(time_periods) - 1 #assume time starts at 1, first is min after sort :)
  if(time_offset != 0){
    dt[, G := G-time_offset]
    dt[, time := time-time_offset]
    time_periods <- time_periods - time_offset
  }
  
  time_step <- 1 #time may not jump at 1
  if(any(time_periods[2:length(time_periods)] - time_periods[1:length(time_periods)-1] != 1)){
    time_step <- time_periods[2]-time_periods[1]
    time_periods <- (time_periods-1)/time_step+1
    if(any(time_periods[2:length(time_periods)] - time_periods[1:length(time_periods)-1] != 1)){stop("time step is not uniform")}
    dt[G != 1, G := (G-1)/time_step+1]
    dt[time != 1, time := (time-1)/time_step+1]
  }
  
  return(list(dt = dt,
              time_change = list(time_step = time_step,
                            max_time = max(time_periods),
                            time_offset = time_offset)))
  
}

get_auxdata <- function(dt, p){

  time_periods <- dt[, unique(time)]
  id_size <- dt[, uniqueN(unit)]
  
  setorder(dt, time, G, unit) #sort the dataset essential for the sort-once-quick-access 
  
  #construct the outcomes list for fast access later
  #loop for multiple outcome
  outcomes_list <- list()
  for(outcol in p$outcomevar){
    outcomes <- list()
    
    if(!p$allow_unbalance_panel){
      for(i in time_periods){
        start <- (i-1)*id_size+1
        end <- i*id_size
        outcomes[[i]] <- dt[seq(start,end), get(outcol)]
      }
    } else {
      
      for(i in time_periods){
        #populate a outcome vector of length N with outcome data in the right place
        #NA is data gone or missing, will be addressed in estimate_did_rc
        outcome_period <- rep(NA, id_size)
        data_pos <- dt[time == i, unit] #units observed in i
        outcome_period[data_pos] <- dt[time == i, get(outcol)]
        outcomes[[i]] <- outcome_period
      }
      
    }
    
    outcomes_list[[outcol]] <- outcomes
  }
  
  #the time-invariant parts 
  if(!p$allow_unbalance_panel){
    dt_inv <- dt[1:id_size]
  } else {
    dt_inv <- dt[dt[, .I[1], by = unit]$V1] #the first observation
    setorder(dt_inv, unit) #can't move this outside
  }
  
  cohorts <- dt_inv[, unique(G)]
  cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]
  
  # the optional columns
  varycovariates <- list()
  if(!allNA(p$varycovariatesvar)){
    for(i in time_periods){
      start <- (i-1)*id_size+1
      end <- i*id_size
      varycovariates[[i]] <- dt[seq(start,end), .SD, .SDcols = p$varycovariatesvar]
    }
  } else {
    varycovariates <- NA
  }
  
  # filters
  filters <- list()
  if(!is.na(p$filtervar)){
    for(i in time_periods){
      start <- (i-1)*id_size+1
      end <- i*id_size
      filters[[i]] <- dt[seq(start,end), .SD, .SDcols = p$filtervar]
    }
  } else {
    filters <- NA
  }
  
  if(!allNA(p$covariatesvar)){
    covariates <- cbind(const = -1, dt_inv[,.SD, .SDcols = p$covariatesvar])
  } else {
    covariates <- NA
  }
  
  if(!is.na(p$clustervar)){
    cluster <- dt_inv[, .SD, .SDcols = p$clustervar] |> unlist()
  } else {
    cluster <- NA
  }
  
  if(!is.na(p$weightvar)){
    weights <- dt_inv[, .SD, .SDcols = p$weightvar] |> unlist()
  } else {
    weights <- rep(1, id_size)
  }
  
  aux <- list(time_periods = time_periods,
              id_size = id_size,
              outcomes_list = outcomes_list,
              dt_inv= dt_inv,
              cohorts = cohorts,
              cohort_sizes = cohort_sizes, 
              varycovariates = varycovariates,
              covariates = covariates,
              cluster = cluster,
              weights = weights,
              filters = filters)
  
  return(aux)
  
}

convert_targets <- function(results, result_type, t){
  
  if(result_type == "dynamic"){
    setnames(results, "target", "event_time")
    
  } else if (result_type == "cohort"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), t$time_offset, t$time_step)]
    setnames(results, "target", "cohort")
    
  } else if (result_type == "calendar"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), t$time_offset, t$time_step)]
    setnames(results, "target", "time")
    
  } else if (result_type == "group_time"){
    
    results[, cohort := floor((target-1)/t$max_time)]
    results[, time := (target-cohort*t$max_time)]
    
    #recover the time
    results[, cohort := recover_time(cohort, t$time_offset, t$time_step)]
    results[, time := recover_time(time, t$time_offset, t$time_step)]
    
    results[, target := NULL]
    
  } else if (result_type == "simple") {
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := NULL]
  } 
  return(results)
}

recover_time <- function(time, time_offset, time_step){
  return(((time-1)*time_step)+1+time_offset)
}

NULL

if(FALSE){
  text <- ". agg_weight att att_cont att_treat attgt cohort cohort_size
    conf_lwb conf_upb const cont_ipw_weight count delta_y element_rect
    element_text event_time pg placeholder post.y pre.y ps s se target
    tau time_fe treat_ipw_weight treat_latent type unit unit_fe weight x
    x2 x_trend y y0 y1 y2 time outcome V1 att_cont_post att_cont_pre att_treat_post att_treat_pre inpost
     inpre max_et min_et new_unit or_delta or_delta_post or_delta_pre
     targeted used"
  text <- text |> str_remove_all("\\\n") |>str_split(" ") |> unlist()
  text <- text[text!=""]
  text <- text |> str_flatten(collapse = "','")
  text <- paste0("c('", text, "')")
}

# quiets concerns of R CMD check re: the .'s that appear in pipelines and data.table variables
utils::globalVariables(c('.','agg_weight','att','att_cont','att_treat','attgt','cohort','cohort_size','conf_lwb','conf_upb',
                         'const','cont_ipw_weight','count','delta_y','element_rect','element_text','event_time','pg','placeholder',
                         'post.y','pre.y','ps','s','se','target','tau','time_fe',
                         'treat_ipw_weight','treat_latent','type','unit','unit_fe','weight','x','x2',
                         'x_trend','y','y0','y1','y2', 'time', 'weights', 'outcome', "G", "D",
                         'V1','att_cont_post','att_cont_pre','att_treat_post','att_treat_pre','inpost','inpre','max_et','min_et','new_unit','or_delta','or_delta_post','or_delta_pre','targeted','used'))


#' Create an event study plot for Difference-in-Differences (DiD) analysis.
#'
#' This function generates an event study plot based on the results of a DiD analysis.
#'
#' @param dt A data table containing the results of the DiD analysis. It should include columns for 'att' (average treatment effect), 'se' (standard error), and 'event_time' (time points).
#' @param graphname A character string specifying the title of the plot (default is "event study plot").
#' @param note A character string for adding additional notes or comments to the plot (default is empty).
#' @param base_time The time point representing the base period (default is -1).
#' @param significance_level The significance level for confidence intervals (default is 0.05).
#'
#' @return A ggplot2 object representing the event study plot.
#' @export

plot_did_dynamics <-function(dt, 
                             graphname = "event study plot", note = "", base_time = -1, significance_level = 0.05#, 
                             #stratify_offset =0.1
){
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("The ggplot2 package must be installed to use plotting functions")
    #Either exit or do something without rgl
    return(NULL)
  }
  

  #add the base period
  if("outcome" %in% names(dt)){
    base_row <- data.table(att = 0, se = 0, event_time = base_time, outcome = dt[, unique(outcome)])
  } else {
    base_row <- data.table(att = 0, se = 0, event_time = base_time)
  }
  
  #dt <- dt |> add_base_row(base_row)
  dt <- dt |> rbind(base_row, fill = TRUE)
  
  dt[, conf_upb := att + qnorm(1-significance_level/2) * se]
  dt[, conf_lwb := att - qnorm(1-significance_level/2) * se]
  

  figure <- dt |> 
    ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dashed', col = 'red')
  
  figure <- figure + ggplot2::geom_line(ggplot2::aes(x = event_time, y = att), color = "black") + 
    ggplot2::geom_point(ggplot2::aes(x = event_time, y = att), color = "black") +
    ggplot2::geom_errorbar(ggplot2::aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                width = 0.1, linetype = "dashed")

  if("outcome" %in% names(dt)){
    figure <- figure + ggplot2::facet_wrap(~outcome, scales = "free")
  }

  figure <- figure + ggplot2::theme_classic()

  return(figure)
  
}

#' Simulate a Difference-in-Differences (DiD) dataset
#'
#' Simulates a dataset for a Difference-in-Differences analysis with various customizable options.
#'
#' @param sample_size The number of units in the dataset.
#' @param time_period The number of time periods in the dataset.
#' @param untreated_prop The proportion of untreated units.
#' @param epsilon_size The standard deviation for the error term in potential outcomes.
#' @param cov The type of covariate to include ("no", "int", or "cont").
#' @param hetero The type of heterogeneity in treatment effects ("all" or "dynamic").
#' @param second_outcome Whether to include a second outcome variable.
#' @param second_cov Whether to include a second covariate.
#' @param na Whether to generate missing data ("none", "y", "x", or "both").
#' @param balanced Whether to balance the dataset by random sampling.
#' @param seed Seed for random number generation.
#' @param stratify Whether to stratify the dataset based on a binary covariate.
#' @param treatment_assign The method for treatment assignment ("latent" or "uniform").
#'
#' @return A list containing the simulated dataset (dt) and the treatment effect values (att).
#'
#' @examples
#' # Simulate a DiD dataset with default settings
#' data <- sim_did(sample_size = 100, time_period = 5)
#'
#' # Simulate a DiD dataset with customized settings
#' data <- sim_did(sample_size = 200, time_period = 8, cov = "int", hetero = "dynamic")
#'
#' @export
sim_did <- function(sample_size, time_period, untreated_prop = 0.3, epsilon_size = 0.001,
                    cov = "no", hetero = "all", second_outcome = FALSE, second_cov = FALSE, vary_cov = FALSE, na = "none", 
                    balanced = TRUE, seed = NA, stratify = FALSE, treatment_assign = "latent"){
  
  if(!is.na(seed)){set.seed(seed)}
  
  #unit  -------------
  dt_i <- data.table(unit = 1:sample_size)
  if(cov == "int"){
    dt_i[, x := sample.int(5, sample_size, replace = TRUE)] 
  } else if (cov == "no"){
    dt_i[, x := 1] 
  } else if (cov == "cont"){
    dt_i[, x := rnorm(sample_size)]
  }
  
  if(second_cov){
    dt_i[, x2 := rnorm(sample_size)]
  } else {dt_i[, x2 := 1]}
  
  if(stratify){
    dt_i[, s := fifelse(rnorm(sample_size) > 0, 1, 2)]
  } else {
    dt_i[, s := 1]
  }

  #treatment assignment ---------------------
  
  #assign treated group based on a latent related to X
  
  if(treatment_assign == "latent"){
    dt_i[, treat_latent := x*0.2 + x2*0.2 + rnorm(sample_size)] #unit with larger X tend to be treated and treated earlier
    untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
    dt_i[treat_latent <= untreated_thres, G := Inf] #unit with low latent is never treated
    
    cohort_prop <- (1-untreated_prop)/(time_period-1)
    last_treat_thres <- untreated_thres
    for(t in time_period:2){ #unit start getting treated in t = 2
      treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop*(time_period - t + 1))
      dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
      last_treat_thres <- treat_thres
    }
    rm(t)
  } else if (treatment_assign == "uniform"){
    #when treatment is set to 'uniform', untreated propensity is fixed
    dt_i[,G := floor((unit-1)/(sample_size/time_period))]
    dt_i[G < 2, G := Inf]
  }

  #assign unit FE
  dt_i[, unit_fe := rnorm(sample_size)]
  
  #time ------------------
  
  dt_t <- data.table(time = 1:time_period)
  dt_t[, time_fe := rnorm(time_period)]
  dt_t[, x_trend := rnorm(time_period)]
  
  #panel --------------------------
  dt <- CJ(unit = 1:sample_size, time = 1:time_period)
  dt <- dt |> merge(dt_i, by = "unit")
  dt <- dt |> merge(dt_t, by = "time")
  
  dt[, D := as.integer(time >= G)]
  
  #add time_varying covariates
  if(vary_cov){
    dt[, xvar := pmin(G, time_period+4)*time+rnorm(sample_size*time_period, 0,30)] #should be confounding....?
  } else {
    dt[, xvar := 1]
  }
  
  #untreated potential outcomes
  dt[, y0 := unit_fe + time_fe + (x+x2)*x_trend + xvar + rnorm(sample_size*time_period, sd = epsilon_size)]
  
  #generate gtatt
  att <- CJ(G = 1:time_period, time = 1:time_period)
  if(hetero == "all"){
    att[, attgt := rnorm(time_period*time_period, mean = 2, sd = 1)]
  } else if (hetero == "dynamic"){
    for(event_t in 0:max(att[,time-G])){
      att[time - G == event_t, attgt := rnorm(1, mean = 2, sd = 1)]
    }
  }
  att[time < G, attgt := 0] #no anticipation
  
  dt <- dt |> merge(att, by = c("time", "G"), all.x = TRUE, all.y = FALSE)
  dt[is.na(attgt), attgt := 0]
  dt[, tau := attgt*s]
  
  #potential outcome
  dt[, y1 := y0 + tau]
  dt[, y := y1*D + y0*(1-D)]
  dt <- dt[, .SD, .SDcols = c("time", "G", "unit", "x", "x2", "y", "s", "xvar")]
  
  #additional -----------------
  
  if(na == "y"){
    dt[, y := na_insert(y)]
  } else if (na == "x") {
    dt[, x := na_insert(x)]
  } else if( na == "both"){
    dt[, y := na_insert(y)]
    dt[, x := na_insert(x)]
  }
  
  if(balanced == FALSE){
    size <- fnrow(dt)
    dt <- dt[sample(1:size, size*0.99)]
  }
  if(second_outcome == TRUE){  dt[, y2 := y + 1 + rnorm(fnrow(dt), 0, 0.1)]}
  
  return(list(dt = dt, att = att))
  
}





set_max_thread <- function(){
  data.table::setDTthreads(0)
  options(kit.nThread = getDTthreads())
}

reverse_col <- function(x){
  return(x[,ncol(x):1])
}

validate_did <- function(dt,varnames,p){

  raw_unit_size <- dt[, uniqueN(unit)]
  raw_time_size <- dt[, uniqueN(time)]
  
  if(!is.na(p$balanced_event_time)){
    if(p$balanced_event_time > dt[, max(time-G)]){stop("balanced_event_time is larger than the max event time in the data")}
  }
  
  if(!is.na(p$filtervar) && !is.logical(dt[[p$filtervar]])){
    stop("filter var needs to be a logical column")
  }
  
  #doesn't allow missing value for now
  for(col in varnames){
    if(is.na(col)){next}
    na_obs <- whichNA(dt[[col]])
    if(length(na_obs) != 0){
      warning("missing values detected in ", col, ", removing ", length(na_obs), " observation.")
      dt <- dt[!na_obs]
    }
  }
  
  if(!allNA(p$covariatesvar) && uniqueN(dt, by = c("unit", p$covariatesvar)) > raw_unit_size){
    warning("some covariates is time-varying, fastdid only use the first observation for covariates.")
  }
  
  
  if(!allNA(p$covariatesvar)|!allNA(p$varycovariatesvar)){
    for(cov in c(p$covariatesvar, p$varycovariatesvar)){
      if(is.na(cov)){next}
      #check covaraites is not constant  
      if(fnunique(dt[, get(cov)[1], by = "unit"][, V1]) == 1)stop(cov, " have no variation")
    }
  }
  
  #check balanced panel
  #check if any is dup
  if(anyDuplicated(dt[, .(unit, time)])){
    dup_id <- dt[duplicated(dt[,.(unit, time)]), unique(unit)]
    stop(length(dup_id), " units is observed more than once in a period.")
  }
  
  #check if any is missing
  if(!p$allow_unbalance_panel){
    unit_count <- dt[, .(count = .N), by = unit]
    if(any(unit_count[, count < raw_time_size])){
      mis_unit <- unit_count[count < raw_time_size]
      warning(nrow(mis_unit), " units is missing in some periods, enforcing balanced panel by dropping them")
      dt <- dt[!unit %in% mis_unit[, unit]]
    }
  }
  return(dt)
}

