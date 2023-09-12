
test_create_event_data <- function(){
  
  dt <- generate_sim_dt()[["dt"]]
  
  event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                          cohortvar = "G",
                                          covariate_base_balance = "x",
                                          covariate_base_stratify = "s",
                                          balanced_panel = TRUE,
                                          never_treat_action = "both")
  
  expect_equal(nrow(event_panel), 5472,
               info = "nrow after create_event_panel")
  
}

test_dynamic <- function(p){
  
  #generate estimation and att
  sim_dt <- generate_sim_dt(p, "dynamic")
  att <- sim_dt$att
  dt <- sim_dt$dt
  dynamic_est <- generate_est(dt, p, "dynamic")
  
  #process att 
  att[, event_time := time-G]
  att_dynamic <- att[!event_time %in% c(-1,p$time_period-1), .(attgt = mean(attgt)), by = "event_time"]
  
  #merge both
  merged_att <- merge(dynamic_est, att_dynamic, by = "event_time")
  
  ratio <- get_att_in_ci_ratio(merged_att)

  return(expect_equal(ratio, 1, info = "att in dynamic's CI"))
  
}

test_cohort_event_time <- function(p){
  
  #generate estimation and att
  sim_dt <- generate_sim_dt(p, "cohort_event_time")
  att <- sim_dt$att
  dt <- sim_dt$dt
  cohort_time_est <- generate_est(dt, p, "cohort_event_time")
  
  #process att
  att[, event_time := time-G]
  cohort_time_att <- att[!event_time %in% c(-1,p$time_period-1),]
  
  #process est
  get_cohort <- function(x){
    start <- str_locate_all(x, "\\.")[[1]][1,1]
    end <- str_locate_all(x, "\\.")[[1]][2,1]
    return(as.numeric(str_sub(x, start + 1, end - 1)))
  }
  cohort_time_est[,G:= lapply(variable, get_cohort)]
  cohort_time_est[, G:= unlist(G)]
  
  #merge both
  merged_att <- merge(cohort_time_est, cohort_time_att, by = c("event_time", "G"))
  
  ratio <- get_att_in_ci_ratio(merged_att)
  
  return(expect_equal(ratio, 1, info = "att in cohort_event_time's CI"))
  
}

# helper function --------------------

generate_sim_dt <- function(p = list(sample_size = 100, time_period = 10),type = "dynamic"){
  hetero_type <- ifelse(type == "cohort_event_time", "all", "dynamic")
  simdt <- sim_did(p$sample_size, p$time_period, cov = "int", hetero = hetero_type, balanced = FALSE, second_outcome = FALSE, seed = 1)
  dt <- simdt$dt
  att <- simdt$att
  return(list(dt = dt, att = att))
}

generate_est <- function(dt, p,type){
  
  event_panel <- suppressWarnings(create_event_data(dt, timevar = p$t_name, unitvar = p$unit_name, 
                                          cohortvar = p$cohort_name,
                                          covariate_base_balance = p$balance_name,
                                          covariate_base_stratify = p$stratify_name,
                                          balanced_panel = TRUE,
                                          never_treat_action = "both"))
  
  est <- suppressMessages(suppressWarnings(
    get_event_result(event_panel, variable = p$y_name, trends = FALSE, mem.clean = FALSE, result_type = type)
    ))
  
  return(est)
  
}

get_att_in_ci_ratio <- function(merged_att){
  
  merged_att[, ci_ub := Estimate+`Std. Error`*1.96]
  merged_att[, ci_lb := Estimate-`Std. Error`*1.96]
  merged_att[, s :=  (function(x){str_sub(str_trim(x), str_length(str_trim(x)), str_length(str_trim(x)))})(variable)]
  merged_att[, attgt := attgt*as.numeric(s)]
  merged_att[, par_in_ci := (attgt <= ci_ub & attgt >= ci_lb)]
  ratio <- merged_att[, mean(par_in_ci)]
  
  return(ratio)
  
}

# test ctatt result -----------------------------------------------------------------------------------

temp <- function(){
  
  simdt <- sim_did(sample_size, time_period, cov = "int", hetero = "all", balanced = FALSE, second_outcome = FALSE, seed = 1)
  dt <- simdt$dt
  att_ct <- simdt$att
  
  event_panel <- copy(dt) #copying so that the original does not change
  
  event_panel <- event_panel %>% create_event_data(timevar = t_name, unitvar = unit_name, 
                                                   cohortvar = cohort_name,
                                                   covariate_base_balance = balance_name,
                                                   covariate_base_stratify = stratify_name,
                                                   balanced_panel = TRUE,
                                                   never_treat_action = "both")

  means <- get_event_result(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE, result_type = "means")
  pooled <- get_event_result(event_panel, variable = y_name, trends = FALSE, mem.clean = FALSE, result_type = "pooled")
  
}


