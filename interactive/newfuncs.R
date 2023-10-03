
fastdid <- function(dt,
                    weightvar=NULL,clustervar=NULL,covariatesvar = c(),
                    control_option="both",result_type="group_time",balanced_composition=FALSE,
                    boot=FALSE){
  setorder(dt, time, G, unit)
  
  time_periods <- dt[, unique(time)]
  time_size <- length(time_periods)
  id_size <- dt[, uniqueN(unit)]
  
  dt_inv <- dt[1:id_size]
  cohorts <- dt_inv[, unique(G)]
  cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]
  
  outcomes <- list()
  for(i in time_periods){
    start <- (i-1)*id_size+1
    end <- i*id_size
    outcomes[[i]] <- dt[seq(start,end), .(y)]
  }
  
  # fill 
  
  if(length(covariatesvar)>0){
    covariates <- dt_inv[,.SD, .SDcols = covariatesvar]
  } else {covariates <- c()}
  
  if(!is.null(clustervar)){
    cluster <- dt_inv[,.SD, .SDcols = clustervar]
  } else {cluster <- NULL}
  
  if(!is.null(weightvar)){
    weights <- dt_inv[, .SD, .SDcols = weightvar]
  } else {weights <- rep(1,id_size)}
  
  
  # validate data -----------------------------------------------------
  
  if(any(sapply(covariates, sd) == 0)){stop("some covariates have no variation")}
  
  # the main process -------------------------------------------------
  
  
  # attgt
  gt_result <- estimate_gt_att(outcomes, covariates, weights,
                               cohort_sizes,cohorts,id_size,time_periods, #info about the dt
                               control_option)
  
  # aggregate att and inf function
  agg_result <- aggregate_gt_result(gt_result, cohort_sizes,
                                    result_type)
  
  
  #influence function -> se
  agg_se <- get_se(agg_result$inf_matrix, boot, biters = 1000, cluster)
  
  #gather results
  result <- data.table(target = agg_result$targets, att = agg_result$agg_att, se = agg_se)
  return(result)
}


get_cohort_pos <- function(cohort_sizes, start_cohort, end_cohort = start_cohort){
  #if(!start_cohort %in% cohort_sizes[,unique(G)]|!end_cohort %in% cohort_sizes[,unique(G)]) {stop("cohort not in cohort_sizes")}
  start <- cohort_sizes[G < start_cohort, sum(cohort_size)]+1
  end <- cohort_sizes[G <= end_cohort, sum(cohort_size)]
  return(start:end)
}

estimate_gt_att <- function(outcomes, covariates, weights, 
                            cohort_sizes,cohorts,id_size,time_periods,
                            control_option) {
  
  treated_cohorts <- cohorts[!is.infinite(cohorts)]
  
  if(!Inf %in% cohorts & control_option != "notyet"){
    warning("no never-treated availble, switching to not-yet-treated control")
    control_option <- "notyet"
  }
  max_control_cohort <- ifelse(control_option == "notyet", max(treated_cohorts), Inf) 
  
  last_coef <- NULL
  gt_att <- data.table()
  gt_inf_func <- data.table(placeholder = rep(NA, id_size))
  for(t in time_periods){
    for(g in cohorts){
      
      #setup and checks
      base_period <- g-1
      min_control_cohort <- ifelse(control_option == "never", Inf, max(t+1, base_period+1)) #not-yet treated / never treated in both base and "treated" period
      if(t == base_period){next} #no treatmenteffect for the base period
      if(g >= max_control_cohort){next} #no treatmenteffect for never treated or the last treated cohort (notyet)
      if(t >= max_control_cohort){next} #no control available if the last cohort is treated too
      
      #select the right cohorts
      did_setup <- rep(NA, id_size)
      did_setup[get_cohort_pos(cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
      did_setup[get_cohort_pos(cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
      
      #select the post and pre outcome
      pre_outcome <- outcomes[[base_period]]
      post_outcome <- outcomes[[t]] 
      
      #construct the 2x2 dataset
      cohort_did <- data.table(D = did_setup, post = post_outcome, pre = pre_outcome, cov = covariates)
      cohort_did[, weights := weights] #the default weight, should allow overwrite next
      
      #estimate did
      results <- estimate_did(cohort_did, last_coef)
      last_coef <- results$logit_coef
      
      #collect the result
      gt_att <- rbind(gt_att, data.table(G = g, time = t, att = results$att))
      gt_inf_func[[paste0(g, ".", t)]] <- results$inf_func
      
    }
  }
  gt_inf_func[, placeholder := NULL]
  return(list(att = gt_att, inf_func = gt_inf_func))
}

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

get_aggregate_scheme <- function(group_time, result_type, balanced_composition){
  weights <- data.table()
  gt_count <- group_time[, .N]
  
  if(balanced_composition){stop("balanced_composition not implemented yet")}
  
  bool_to_pn <- function(x){ifelse(x, 1, -1)}
  
  if (result_type == "dynamic") {
    group_time[, target := time-G]
  } else if (result_type == "group") {
    group_time[, target := G*(bool_to_pn(time>=G))]
  } else if (result_type == "time") {
    group_time[, target := time*(bool_to_pn(time>=G))]
  } 
  
  min_target<- group_time[, min(target)]
  max_target<- group_time[, max(target)]
  
  targets <- group_time[, unique(target)]
  for(tar in targets){
    group_time[, agg_weight := 0]
    total_size <- group_time[target == tar, sum(cohort_size)]
    group_time[, weight := ifelse(target == tar, cohort_size/total_size, 0)]
    target_weights <- group_time[, .(weight)] |> transpose()
    names(target_weights) <- names(gt_inf_func)
    
    weights <- rbind(weights, target_weights)
  }
  return(list(weights = weights, targets = targets))
}

aggregate_gt_result <- function(gt_result, cohort_sizes,
                                result_type){
  
  gt_att <- gt_result$att
  gt_inf_func <- gt_result$inf_func
  
  group_time <- gt_att[,.(G, time)] |> merge(cohort_sizes, by = "G")
  
  if(result_type == "group_time"){
    targets <- group_time[, unique(G*max(time)+time)]
    inf_matrix <- as.matrix(gt_inf_func)
    agg_att <- as.vector(gt_att[, att])
    
  } else {
    
    
    agg_sch <- get_aggregate_scheme(group_time, result_type, balanced_composition)
    targets <- agg_sch$targets
    weights <- agg_sch$weights
    
    inf_matrix <- as.matrix(gt_inf_func) %*% t(as.matrix(weights))
    agg_att <- as.matrix(weights) %*% as.vector(gt_att[, att]) |> as.vector()
  }
  return(list(inf_matrix = inf_matrix, agg_att = agg_att, targets = targets))
}

get_se <- function(inf_matrix, boot, biters, cluster) {
  
  if(boot){
    
    top_quant <- 0.75
    bot_quant <- 0.25
    if(!is.null(cluster)){
      
      #aggregate the influence function by cluster
      inf_matrix <- inf_matrix |> as.data.table()
      inf_matrix[, cluster := cluster]
      inf_matrix <- inf_matrix[, lapply(.SD, mean), by = "cluster", .SDcols = names(inf_matrix)[names(inf_matrix) != "cluster"]] 
      inf_matrix[, cluster := NULL]
      inf_matrix <- inf_matrix |> as.matrix()
    }
    
    boot_results <- BMisc::multiplier_bootstrap(inf_matrix, biters = biters) %>% as.data.table()
    boot_top <- boot_results[, lapply(.SD, function(x) quantile(x, top_quant, type=1, na.rm = TRUE)),]
    boot_bot <- boot_results[, lapply(.SD, function(x) quantile(x, bot_quant, type=1, na.rm = TRUE)),]
    
    dt_se <- rbind(boot_bot, boot_top) %>% transpose()
    names(dt_se) <- c("boot_bot", "boot_top")
    dt_se[, n_adjust := nrow(inf_matrix)/colSums(inf_matrix != 0)]
    se <- dt_se[,(boot_top-boot_bot)/(qnorm(top_quant) - qnorm(bot_quant))*n_adjust]
    
  } else {
    
    if(!is.null(cluster)){stop("clustering only available with bootstrap")}
    
    inf_matrix <- inf_matrix %>% as.data.table()
    se <- inf_matrix[, lapply(.SD, function(x) sd(x, na.rm = TRUE)*sqrt(length(x)-1)/length(x[x!=0]))] %>% as.vector() #should maybe use n-1 but did use n
    
  }
  return(unlist(se))
}

