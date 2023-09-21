# test create event data -----------------------------------

test_create_event_data <- function(){
  
  dt <-  sim_did(5, 5, seed = 1, treatment_assign = "uniform", stratify = FALSE)[["dt"]]
  
  expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                       cohortvar = "G",
                                       covariate_base_balance = "x",
                                       covariate_base_stratify = "s",
                                       control_group = "both", verbose = FALSE))
  
  dt <-  sim_did(5, 5, seed = 1, treatment_assign = "uniform", stratify = TRUE)[["dt"]]
  
  expect_silent(event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                                 cohortvar = "G",
                                                 covariate_base_balance = "x",
                                                 covariate_base_stratify = "s",
                                                 control_group = "both", verbose = FALSE))

  
}

test_get_event_result <- function(){
  
  dt <-  sim_did(100, 10, seed = 1, treatment_assign = "uniform", stratify = TRUE, second_outcome = TRUE)[["dt"]]
  
  event_panel <- create_event_data(dt, timevar = "time", unitvar = "unit", 
                                   cohortvar = "G",
                                   covariate_base_balance = "x",
                                   covariate_base_stratify = "s",
                                   control_group = "both", verbose = FALSE)
  expect_silent(event_est <- get_event_result(event_panel, variable = c("y", "y2"), 
                                              trends = FALSE, mem.clean = FALSE, result_type = "dynamic", 
                                              separate_cohort_time = FALSE, separate_stratify = FALSE),
                info = "get_event_result dynamic, no split")
  expect_silent(event_est <- get_event_result(event_panel, variable = c("y", "y2"), 
                                              trends = FALSE, mem.clean = FALSE, result_type = "dynamic"),
                info = "get_event_result dynamic, split by stratify and cohort event time")
  expect_silent(event_est <- get_event_result(event_panel, variable = c("y", "y2"), 
                                              trends = FALSE, mem.clean = FALSE, result_type = "cohort_event_time", separate_cohort_time = TRUE),
                info = "get_event_result, cohort event time, split by stratify and cohort event time")
  
  
}


test_plot_event_dynamics <- function(){
  #generate estimation and att
  dt <- generate_sim_dt()[["dt"]]
  dynamic_est <- generate_est(dt, "dynamic")
  expect_silent(dynamic_est |> plot_event_dynamics(), info = "no error in single outcome 2 stratify")
}

test_no_combine <- function(){
  dt <- generate_sim_dt()[["dt"]]
  event_panel <- suppressWarnings(create_event_data(dt, timevar = "time", unitvar =  "unit", 
                                               cohortvar = "G",
                                               balanced_panel = TRUE,
                                               control_group = "both", copy = FALSE, combine = FALSE, verbose = FALSE))
  expect_silent(event_est <- get_event_result(event_panel, variable = "y", trends = FALSE, mem.clean = FALSE, 
                                              result_type = "cohort_event_time"),
                info = "not combining works.")
}

# test get result -------------------------------------

test_dynamic_est <- function(p){
  
  #generate estimation and att
  sim_dt <- generate_sim_dt(p, "dynamic")
  att <- sim_dt$att
  dt <- sim_dt$dt
  dynamic_est <- generate_est(dt, "dynamic")
  
  #process att 
  att[, event_time := time-G]
  att_dynamic <- att[!event_time %in% c(-1,p$time_period-1), .(attgt = mean(attgt)), by = "event_time"]
  
  #merge both
  merged_att <- merge(dynamic_est, att_dynamic, by = "event_time")
  
  ratio <- get_att_in_ci_ratio(merged_att)

  return(expect_equal(ratio, 1, info = "att in dynamic's CI"))
  
}

test_cohort_event_time_est <- function(p){
  
  #generate estimation and att
  sim_dt <- generate_sim_dt(p, "cohort_event_time")
  att <- sim_dt$att
  dt <- sim_dt$dt
  cohort_time_est <- generate_est(dt, "cohort_event_time")
  
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

test_get_result_means <- function(){
  
  #generate estimation and att
  sim_dt <- generate_sim_dt()
  att <- sim_dt$att
  dt <- sim_dt$dt
  means_est <- generate_est(dt, "means")
  
  return(expect_equal(nrow(means_est), dt[, uniqueN(s)], info = "means nrow is stratify count"))
  
}

test_dynamic_cohort_event_time_consistent <- function(){
  
  #generate estimation and att
  sim_dt <- generate_sim_dt(stratify = FALSE, balanced = TRUE)
  att <- sim_dt$att
  dt <- sim_dt$dt
  
  cohort_pop <- dt[, .(pop = .N), by = "G"]
  
  event_est_ce <- generate_est(dt, "cohort_event_time")
  event_est_ce <- event_est_ce |> merge(cohort_pop, by.x = "cohort", by.y = "G")
  event_est_ce_mean <- event_est_ce[, .(weighted_Estimate = sum(Estimate*pop)/sum(pop)), by = "event_time"]
  
  event_est <-  generate_est(dt, "dynamic")
  
  est_both <- merge(event_est, event_est_ce_mean, by = "event_time")
  est_both[, est_diff := Estimate - weighted_Estimate]
  
  return(expect_equal(est_both[, Estimate], est_both[, weighted_Estimate], tolerance = 1e-13,
                      info = "Dynamic estimate is equal to pop-weighted average of cohort event_time estimate."))
  
}



# helper function --------------------

generate_sim_dt <- function(p = list(sample_size = 100, time_period = 10),type = "dynamic", seed = 1, stratify = TRUE, balanced = FALSE){
  hetero_type <- ifelse(type == "cohort_event_time", "all", "dynamic")
  simdt <- sim_did(p$sample_size, p$time_period, cov = "int", hetero = hetero_type, balanced = balanced, second_outcome = FALSE, seed = seed,
                   stratify = stratify)
  dt <- simdt$dt
  att <- simdt$att
  return(list(dt = dt, att = att))
}

generate_est <- function(dt, type, p = list(t_name = "time", unit_name = "unit",
                                            cohort_name = "G", balance_name = "x", stratify_name = "s",
                                            y_name = "y")){
  
  event_panel <- suppressWarnings(create_event_data(dt, timevar = p$t_name, unitvar = p$unit_name, 
                                          cohortvar = p$cohort_name,
                                          covariate_base_balance = p$balance_name,
                                          covariate_base_stratify = p$stratify_name,
                                          control_group = "both", verbose = FALSE))
  
  est <- get_event_result(event_panel, variable = p$y_name, trends = FALSE, mem.clean = FALSE, result_type = type)

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

