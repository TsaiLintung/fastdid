#2023-10-06
aggregate_gt <- function(gt_result, cohort_sizes, 
                         id_weights, id_cohorts,
                         result_type){
  
  gt_att <- gt_result$att
  gt_inf_func <- gt_result$inf_func
  gt <- gt_result$gt
  
  id_dt <- data.table(weight = id_weights/sum(id_weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group_time <- gt |> merge(pg_dt, by = "G")
  
  setorder(group_time, time, G) #change the order to match the order in gtatt
  
  gt_inf_func <- as.matrix(gt_inf_func)
  
  if(result_type == "group_time"){
    
    targets <- group_time[, unique(G*max(time)+time)]
    inf_matrix <- gt_inf_func
    agg_att <- as.vector(gt_att)
    
  } else {
    
    agg_sch <- get_aggregate_scheme(group_time, result_type, id_weights, id_cohorts)
    targets <- agg_sch$targets
    weights <- as.matrix(agg_sch$weights)
    
    #aggregated att
    agg_att <- weights %*% gt_att
    
    
    
    #get the influence from weight estimation
    #this needs to be optimized!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    inf_weights <- sapply(asplit(weights, 1), function (x){
      get_weight_influence(x, gt_att, id_weights, id_cohorts, group_time[, .(G, time)])
    })
    
    #aggregated influence function
    inf_matrix <- (gt_inf_func %*% t(weights)) + inf_weights 

  }
  return(list(inf_matrix = inf_matrix, agg_att = agg_att, targets = targets))
}

get_aggregate_scheme <- function(group_time, result_type, id_weights, id_cohorts){
  
  weights <- data.table()
  gt_count <- group_time[, .N]
  
  bool_to_pn <- function(x){ifelse(x, 1, -1)}
  
  if (result_type == "dynamic") {
    group_time[, target := time-G]
  } else if (result_type == "group") {
    group_time[, target := G*(bool_to_pn(time>=G))]
  } else if (result_type == "time") {
    group_time[, target := time*(bool_to_pn(time>=G))]
  } else if (result_type == "simple") {
    group_time[, target := bool_to_pn(time>=G)]
  } 
  
  min_target<- group_time[, min(target)]
  max_target<- group_time[, max(target)]
  
  targets <- group_time[, unique(target)]
  targets <- sort(as.integer(targets))
  for(tar in targets){ #the order matters
    
    group_time[, agg_weight := 0]
    total_pg <- group_time[target == tar, sum(pg)]
    group_time[, weight := ifelse(target == tar, pg/total_pg, 0)]
    target_weights <- group_time[, .(weight)] |> transpose()
    
    weights <- rbind(weights, target_weights)
  }
  
  return(list(weights = weights, targets = targets))
}

get_weight_influence <- function(agg_weights, gt_att, id_weights, id_cohorts, group) {

  keepers <- which(agg_weights > 0)
  
  id_dt <- data.table(weight = id_weights/sum(id_weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group <- group |> merge(pg_dt, by = "G")
  
  group[, time := as.integer(time)]
  group[, G := as.integer(G)]
  setorder(group, time, G)

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (id_weights*BMisc::TorF(id_cohorts == group[k,G]) - group[k,pg]) /
      sum(group[keepers,pg])
  })
  # effect of estimating weights in the denominator
  if2 <- base::rowSums(sapply(keepers, function(k) {
    id_weights*BMisc::TorF(id_cohorts == group[k,G]) - group[k,pg]
  })) %*%
    t(group[keepers,pg]/(sum(group[keepers,pg])^2))
  # return the influence function for the weights
  inf_weight <- (if1 - if2) %*% as.vector(gt_att[keepers])
  inf_weight[abs(inf_weight) < sqrt(.Machine$double.eps)*10] <- 0 #fill zero 
  return(inf_weight)
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
    logit_coef <-  prop_score_est$coefficients 
    logit_coef[is.na(logit_coef)|abs(logit_coef) > 1e10] <- 0 #put extreme value and na to 0
    
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
  
  inf_func <- rep(0, oldn) #the default needs to be 0 for the matrix multiplication
  inf_func_no_na <- inf_func_no_na * oldn / n #adjust the value such that mean over the whole id size give the right result
  inf_func[data_pos] <- inf_func_no_na
  return(list(att = att, inf_func = inf_func, logit_coef = logit_coef))
}

estimate_gtatt <- function(outcomes, covariates, weights, 
                           cohort_sizes,cohorts,id_size,time_periods,
                           control_option) {
  
  treated_cohorts <- cohorts[!is.infinite(cohorts)]
  
  if(!Inf %in% cohorts & control_option != "notyet"){
    warning("no never-treated availble, switching to not-yet-treated control")
    control_option <- "notyet"
  }
  max_control_cohort <- ifelse(control_option == "notyet", max(treated_cohorts), Inf) 
  
  last_coef <- NULL
  gt <- data.table()
  gt_att <- c()
  gt_inf_func <- data.table(placeholder = rep(NA, id_size))
  for(t in time_periods){
    for(g in cohorts){
      
      #setup and checks
      base_period <- g-1
      min_control_cohort <- ifelse(control_option == "never", Inf, max(t+1, base_period+1)) #not-yet treated / never treated in both base and "treated" period
      if(t == base_period){next} #no treatmenteffect for the base period
      if(g >= max_control_cohort){next} #no treatmenteffect for never treated or the last treated cohort (notyet)
      if(t >= max_control_cohort){next} #no control available if the last cohort is treated too
      
      #select the control and treated cohorts
      did_setup <- rep(NA, id_size)
      did_setup[get_cohort_pos(cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
      did_setup[get_cohort_pos(cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
      
      #construct the 2x2 dataset
      cohort_did <- data.table(D = did_setup, post = outcomes[[t]], pre = outcomes[[base_period]], 
                               cov = covariates, weights = weights)
      
      #estimate did
      results <- estimate_did(cohort_did, last_coef)
      last_coef <- results$logit_coef
      
      #collect the result
      gt <- rbind(gt, data.table(G = g, time = t)) #the sequence matters for the weights
      gt_att <- c(gt_att, att = results$att)
      gt_inf_func[[paste0(g, ".", t)]] <- results$inf_func
      
    }
  }
  gt_inf_func[,placeholder := NULL]
  names(gt_att) <- names(gt_inf_func)
  return(list(att = gt_att, inf_func = gt_inf_func, gt = gt))
}

get_cohort_pos <- function(cohort_sizes, start_cohort, end_cohort = start_cohort){
  #if(!start_cohort %in% cohort_sizes[,unique(G)]|!end_cohort %in% cohort_sizes[,unique(G)]) {stop("cohort not in cohort_sizes")}
  start <- cohort_sizes[G < start_cohort, sum(cohort_size)]+1
  end <- cohort_sizes[G <= end_cohort, sum(cohort_size)]
  return(start:end)
}

#' Fast Staggered DID Estimation
#'
#' Performs Difference-in-Differences (DID) estimation fast.
#'
#' @param dt A data table containing the panel data.
#' @param timevar The name of the time variable.
#' @param cohortvar The name of the cohort (group) variable.
#' @param unitvar The name of the unit (id) variable.
#' @param outcomevar The name of the outcome variable.
#' @param control_option The control units used for the DiD estimates. Default is "both".
#' @param result_type A character string indicating the type of result to be returned. Default is "group_time".
#' @param boot Logical, indicating whether bootstrapping should be performed. Default is FALSE.
#' @param biters The number of bootstrap iterations. Only relevant if boot = TRUE. Default is 1000.
#' @param weightvar The name of the weight variable (optional).
#' @param clustervar The name of the cluster variable, can only be used when boot == TRUE (optional).
#' @param covariatesvar A character vector containing the names of covariate variables (optional).
#' @param copy whether to copy the dataset before processing, set to true if the original dataset is to be re-used.
#' @param validate whether to validate the dataset before processing.
#' 
#' @import data.table speedglm stringr collapse dreamerr BMisc
#' 
#' @return A data.table containing the estimated treatment effects and standard errors.
#' @export
#'
#' @examples
#' # simulated data
#' simdt <- sim_did(1e+03, 10, cov = "cont", second_cov = TRUE)
#' dt <- simdt$dt
#' 
#' #basic call
#' result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time")
#' 
#' #control for covariates
#' result2 <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
#'                   covariatesvar = c("x", "x2"))
#'                   
#' #bootstrap and clustering
#' result3 <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
#'                   boot = TRUE, clustervar = "x")
#'
#' @keywords difference-in-differences fast computation panel data estimation did
fastdid <- function(dt,
                    timevar, cohortvar, unitvar, outcomevar, 
                    control_option="both",result_type="group_time", 
                    boot=FALSE, biters = 1000,
                    weightvar=NULL,clustervar=NULL,covariatesvar = NULL,
                    copy = FALSE, validate = TRUE
                    ){
  
  
  # validation arguments --------------------------------------------------------
  
  if(copy){dt <- copy(dt)}
  
  if(!is.data.table(dt)) stop("rawdata must be a data.table")
  
  dt_names <- names(dt)
  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_arg(timevar, unitvar, cohortvar, "scalar charin", .choices = dt_names, .message = name_message)
  
  covariate_message <- "__ARG__ must be NULL or a character vector which are all names of columns from the dataset."
  check_arg(covariatesvar, 
            "NULL | multi charin", .choices = dt_names, .message = covariate_message)
  
  checkvar_message <- "__ARG__ must be NULL or a character scalar if a name of columns from the dataset."
  check_arg(weightvar, clustervar,
            "NULL | scalar charin", .choices = dt_names, .message = checkvar_message)
  
  check_arg(control_option, "scalar charin", .choices = c("both", "never", "notyet")) #kinda bad since did's notyet include both notyet and never
  check_arg(copy, validate, "scalar logical")
  
  
  setnames(dt, c(timevar, cohortvar, unitvar, outcomevar), c("time", "G", "unit", "y"))

  # validate data -----------------------------------------------------
  
  if(validate){
    dt <- validate_did(dt, covariatesvar)
  }
  
  # preprocess -----------------------------------------------------------
  
  #make dt conform to the WLOG assumptions of fastdid
  
  #change to int before sorting
  if(!is.numeric(dt[, G])){
    dt[, G := as.numeric(G)]
  }
  if(!is.numeric(dt[, time])){
    dt[, time := as.numeric(time)] 
  }
  
  setorder(dt, time, G, unit) #sort the dataset essential for the sort-once-quick-access 
  
  #deal with time
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
  
  #the outcomes list for fast access later
  id_size <- dt[, uniqueN(unit)]
  outcomes <- list()
  for(i in time_periods){
    start <- (i-1)*id_size+1
    end <- i*id_size
    outcomes[[i]] <- dt[seq(start,end), .(y)]
  }

  #the time-invariant parts 
  dt_inv <- dt[1:id_size]
  cohorts <- dt_inv[, unique(G)]
  cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]

  # the optional columns
  if(!is.null(covariatesvar)){
    covariates <- dt_inv[,.SD, .SDcols = covariatesvar]
  } else {covariates <- NULL}
  
  if(!is.null(clustervar)){
    cluster <- dt_inv[, .SD, .SDcols = clustervar] |> unlist()
  } else {cluster <- NULL}
  
  if(!is.null(weightvar)){
    weights <- dt_inv[, .SD, .SDcols = weightvar] |> unlist()
  } else {weights <- rep(1,id_size)}
  
  # main part  -------------------------------------------------
  
  # attgt
  gt_result <- estimate_gtatt(outcomes, covariates, weights,
                              cohort_sizes,cohorts,id_size,time_periods, #info about the dt
                              control_option)
  
  # aggregate att and inf function
  agg_result <- aggregate_gt(gt_result, cohort_sizes, 
                             weights, dt_inv[, G],
                             result_type)
  
  #get se from the influence function
  agg_se <- get_se(agg_result$inf_matrix, boot, biters, cluster)
  
  # post process -----------------------------------------------

  result <- data.table(agg_result$targets, agg_result$agg_att, agg_se)
  names(result) <- c("target", "att", "se")
  
  #convert "targets" back to meaningful parameter identifiers like cohort 1 post, time 2 post 
  result <- result |> convert_targets(result_type, time_offset, time_step, max(time_periods))
  
  return(result)
  
}

convert_targets <- function(results, result_type, 
                            time_offset, time_step, max_time){
  
  if(result_type == "dynamic"){
    setnames(results, "target", "event_time")
    
  } else if (result_type == "cohort"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), time_offset, time_step)]
    setnames(results, "target", "cohort")
    
  } else if (result_type == "calendar"){
    
    results[, type := ifelse(target >= 0, "post", "pre")]
    results[, target := recover_time(abs(target), time_offset, time_step)]
    setnames(results, "target", "time")
    
  } else if (result_type == "group_time"){
    
    results[, cohort := floor((target-1)/max_time)]
    results[, time := (target-cohort*max_time)]
    
    #recover the time
    results[, cohort := recover_time(cohort, time_offset, time_step)]
    results[, time := recover_time(time, time_offset, time_step)]
    
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

get_se <- function(inf_matrix, boot, biters, cluster) {
  
  if(boot){

    top_quant <- 0.75
    bot_quant <- 0.25
    if(!is.null(cluster)){
      cluster_n <- aggregate(cluster, by=list(cluster), length)[,2]
      inf_matrix <- fsum(inf_matrix, cluster) / cluster_n #the mean without 0 for each cluster of each setting
    }

    boot_results <- BMisc::multiplier_bootstrap(inf_matrix, biters = biters) %>% as.data.table()
    
    boot_top <- boot_results[, lapply(.SD, function(x) quantile(x, top_quant, type=1, na.rm = TRUE)),]
    boot_bot <- boot_results[, lapply(.SD, function(x) quantile(x, bot_quant, type=1, na.rm = TRUE)),]
    
    dt_se <- rbind(boot_bot, boot_top) %>% transpose()
    names(dt_se) <- c("boot_bot", "boot_top")

    se <- dt_se[,(boot_top-boot_bot)/(qnorm(top_quant) - qnorm(bot_quant))]
    se[se < sqrt(.Machine$double.eps)*10] <- NA
    
  } else {
    if(!is.null(cluster)){stop("clustering only available with bootstrap")}
    inf_matrix <- inf_matrix  |> as.data.table()
    se <- inf_matrix[, lapply(.SD, function(x) sqrt(sum(x^2, na.rm = TRUE)/length(x)^2))] %>% as.vector() #should maybe use n-1 but did use n
    
  }
  return(unlist(se))
}


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
#'
#' @import ggplot2

plot_did_dynamics <-function(dt, 
                             graphname = "event study plot", note = "", base_time = -1, significance_level = 0.05, 
                             stratify_offset =0.1
){


  #add the base period
  base_row <- data.table(att = 0, se = 0, event_time = base_time)
  #dt <- dt |> add_base_row(base_row)
  dt <- dt |> rbind(base_row)
  
  dt[, conf_upb := att + qnorm(1-significance_level/2) * se]
  dt[, conf_lwb := att - qnorm(1-significance_level/2) * se]
  
  #add some offset
  # if("stratify" %in% names(dt)){
  #   dt[, event_time := event_time + stratify_offset*(as.integer(stratify)-1)]
  # }
  
  figure <- dt |> 
    ggplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'red')
  
  # if("stratify" %in% names(dt)){
  #   figure <- figure + geom_line(aes(x = event_time, y = att, color = stratify)) + 
  #     geom_point(aes(x = event_time, y = att, color = stratify)) +
  #     geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
  #                   width = 0.1, linetype = "dashed")
  # } else {
    figure <- figure + geom_line(aes(x = event_time, y = att), color = "black") + 
    geom_point(aes(x = event_time, y = att), color = "black") +
    geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                  width = 0.1, linetype = "dashed")
  # }
  
  # if("outcome" %in% names(dt)){
  #   figure <- figure + facet_wrap(~ outcome, scales = "free")
  # }
  # if("cohort" %in% names(dt)){
  #   figure <- figure + facet_wrap(~ cohort, scales = "free")
  # }

  figure <- figure +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.background = element_rect(linetype = "dashed", color = "black"),
          legend.box = "horizontal",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0)) +
    labs(title = graphname, subtitle = note)
  
  return(figure)
  
}

# add_base_row <- function(dt, base_row){
#   
#   col_names <- names(dt)
#   opt_cols <- c("outcome", "stratify", "cohort")
#   
#   for(opt_col in opt_cols){
#     if(opt_col %in% col_names){
#       
#       row_list <- list()
#       
#       for (col_name in dt[, unique(get(opt_col))]) {
#         base_row_sub <- copy(base_row)
#         base_row_sub[, c(opt_col) := col_name]
#         row_list <- c(row_list, list(base_row_sub))
#       }
#       
#     } 
#     base_row <- rbindlist(row_list)
#   }
# 
#   dt <- dt |> rbind(base_row)
#   
#   return(dt)
# }

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
                    cov = "no", hetero = "all", second_outcome = FALSE, second_cov = FALSE, na = "none", 
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
  
  #untreated potential outcomes
  dt[, y0 := unit_fe + time_fe + x*x_trend + x2*x_trend + rnorm(sample_size*time_period, sd = epsilon_size)]
  
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
  dt <- dt[, .SD, .SDcols = c("time", "G", "unit", "x", "x2", "y", "s")]
  
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
  setDTthreads(0)
  options(kit.nThread = getDTthreads())
  setFixest_nthreads(getDTthreads())
}

reverse_col <- function(x){
  return(x[,ncol(x):1])
}

validate_did <- function(dt, covariatesvar){
  raw_unit_size <- dt[, uniqueN(unit)]
  raw_time_size <- dt[, uniqueN(time)]
  
  #doesn't allow missing value for now
  cols <- names(dt)
  for(col in cols){
    na_obs <- whichNA(dt[[col]])
    if(length(na_obs) != 0){
      warning("missing values detected in ", col, ", removing ", length(na_obs), " observation.")
      dt <- dt[!na_obs]
    }
  }
  
  if(!is.null(covariatesvar)){
    
    #doesn't allow time varying covariates
    for(cov in covariatesvar){
      if(!all(dt[, all(get(cov) == first(get(cov))), by = "unit"])){stop(cov, " is time-varying")}
    }
    
    #check covaraites is not constant  
    if(any(sapply(dt[1:raw_unit_size, .SD, .SDcols = covariatesvar], sd) == 0)){stop("some covariates have no variation")}
  }
  
  #check balanced panel
  #check if any is dup
  if(anyDuplicated(dt[, .(unit, time)])){
    dup <- duplicated(dt[,.(unit, time)])
    warning(nrow(dup), " units is observed more than once in some periods, enforcing balanced panel by dropping them")
    dt <- dt[!unit %in% dup[, unit]]
  }
  
  #check if any is missing
  unit_count <- dt[, .(count = .N), by = unit]
  if(any(unit_count[, count < raw_time_size])){
    mis_unit <- unit_count[count < raw_time_size]
    warning(nrow(mis_unit), " units is missing in some periods, enforcing balanced panel by dropping them")
    dt <- dt[!unit %in% mis_unit[, unit]]
  }
  return(dt)
}

